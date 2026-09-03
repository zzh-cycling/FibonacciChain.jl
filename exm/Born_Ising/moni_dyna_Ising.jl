using FibonacciChain
using JLD
using LinearAlgebra
using Random
using Statistics

include(joinpath(@__DIR__, "config.jl"))

const ISING_DATA_ROOT = joinpath("exm", "data", "Bulk_measure", "Ising")
const ISING_TOPO_SECTOR_DATA_ROOT =
    joinpath("exm", "data", "Bulk_measure", "Ising_topo_sector")

"""
    initial_topo_sector_state(L) -> Vector{Float64}

Construct the normalized initial state

    |ψ₀⟩ ∝ |+⟩^⊗L + |GHZ⟩,
    |GHZ⟩ = (|0⟩^⊗L + |1⟩^⊗L)/√2.

Since `⟨+|^⊗L GHZ⟩ = 2^((1-L)/2)`, its normalization is
`√(2 + 2^((3-L)/2))`.
"""
function initial_topo_sector_state(L::Integer)
    L >= 1 || throw(ArgumentError("L must be positive"))
    model = ising_model(Int(L))
    dim = length(anyon_basis(model))
    plus_state = fill(inv(sqrt(dim)), dim)
    ghz_state = zeros(Float64, dim)
    ghz_state[1] = inv(sqrt(2))
    ghz_state[end] = inv(sqrt(2))

    state = plus_state + ghz_state
    expected_norm_squared = 2 + 2.0^((3 - L) / 2)
    isapprox(dot(state, state), expected_norm_squared; atol = 1e-12, rtol = 1e-12) ||
        error("Unexpected norm for the topological-sector initial state")
    normalize!(state)
    return state
end

function _validate_ising_run(L::Integer, τ_idx::Integer, t::Integer)
    L >= 1 || throw(ArgumentError("L must be positive"))
    τ_idx in eachindex(τlis) || throw(BoundsError(τlis, τ_idx))
    t >= 1 || throw(ArgumentError("t must be a positive number of periods"))
    return Int(L), Int(τ_idx), Int(t)
end

function _trajectory_path(
    data_root::AbstractString,
    L::Integer,
    τ_idx::Integer,
    t::Integer,
    seed::Integer,
)
    return joinpath(
        data_root,
        "L$(L)",
        "gammaind$(τ_idx)",
        "t$(t)_samples$(seed).jld",
    )
end

function _samples_generate_ising(
    L::Integer,
    τ_idx::Integer,
    seed::Integer,
    t::Integer,
    initial_state::AbstractVector,
    initial_state_label::AbstractString,
    data_root::AbstractString,
)
    L, τ_idx, t = _validate_ising_run(L, τ_idx, t)
    model = ising_model(L)
    length(initial_state) == length(anyon_basis(model)) ||
        throw(DimensionMismatch("initial state has the wrong Hilbert-space dimension"))
    isapprox(norm(initial_state), 1; atol = 1e-12, rtol = 1e-12) ||
        throw(ArgumentError("initial state must be normalized"))

    τ = τlis[τ_idx]
    config = MeasureConfig(
        τ = τ,
        mode = :Born,
        t₂ = t,
        rng = MersenneTwister(seed),
    )
    @time outcome = bulk_evolution(model, collect(initial_state), config)

    n_layers = FibonacciChain.layers_per_period(model)
    size(outcome.samples, 1) == n_layers * t ||
        error("Unexpected trajectory depth")
    length(outcome.free_energys) == n_layers * t ||
        error("Unexpected free-energy time-series length")
    length(outcome.entanglement_entropys) == t ||
        error("Unexpected entanglement time-series length")

    output_path = _trajectory_path(data_root, L, τ_idx, t, seed)
    mkpath(dirname(output_path))
    save(
        output_path,
        "initial_state",
        String(initial_state_label),
        "L",
        L,
        "τ_idx",
        τ_idx,
        "τ",
        τ,
        "gamma",
        tanh(τ),
        "t",
        t,
        "layers_per_period",
        n_layers,
        "sample",
        outcome.samples,
        "sample_free_energy",
        outcome.free_energys,
        "seed",
        Int(seed),
        "halfchain_EE_tlis",
        outcome.entanglement_entropys,
        "final_EElis",
        anyon_eelis(model, outcome.state),
    )
    return output_path
end

"""
    samples_generate(L, τ_idx, seed, t=get_cfg_params_Born(τ_idx, L)[1])

Generate a Born trajectory from `|+⟩^⊗L`. Here and in the output filename,
`t` is the number of complete two-layer measurement periods.
"""
function samples_generate(
    L::Integer,
    τ_idx::Integer,
    seed::Integer,
    t::Integer = get_cfg_params_Born(τ_idx, L)[1];
    data_root::AbstractString = ISING_DATA_ROOT,
)
    model = ising_model(Int(L))
    initial_state = ones(Float64, length(anyon_basis(model)))
    normalize!(initial_state)
    return _samples_generate_ising(
        L,
        τ_idx,
        seed,
        t,
        initial_state,
        "all_plus",
        data_root,
    )
end

"""
    samples_generate_topo_sector(
        L, τ_idx, seed, t=get_cfg_params_Born(τ_idx, L)[1]
    )

Generate a Born trajectory from [`initial_topo_sector_state`](@ref). Here and
in the output filename, `t` is the number of complete two-layer measurement
periods.
"""
function samples_generate_topo_sector(
    L::Integer,
    τ_idx::Integer,
    seed::Integer,
    t::Integer = get_cfg_params_Born(τ_idx, L)[1];
    data_root::AbstractString = ISING_TOPO_SECTOR_DATA_ROOT,
)
    return _samples_generate_ising(
        L,
        τ_idx,
        seed,
        t,
        initial_topo_sector_state(L),
        "topo_sector",
        data_root,
    )
end

function _mean_and_stderr(values::AbstractMatrix)
    sample_count = size(values, 1)
    average = vec(mean(values; dims = 1))
    stderr = sample_count == 1 ?
             zeros(size(values, 2)) : vec(std(values; dims = 1)) / sqrt(sample_count)
    return average, stderr
end

"""
    late_time_drift(values; fraction=0.25)

Compare the penultimate and final time blocks of an ensemble array whose rows
are trajectories and columns are periods. The returned drift is the ensemble
mean of `mean(final block) - mean(previous block)`, with its standard error.
"""
function late_time_drift(values::AbstractMatrix; fraction::Real = 0.25)
    0 < fraction <= 0.5 || throw(ArgumentError("fraction must lie in (0, 0.5]"))
    samples_num, t = size(values)
    block = max(1, floor(Int, fraction * t))
    previous = (t - 2block + 1):(t - block)
    final = (t - block + 1):t
    per_sample = vec(mean(values[:, final]; dims = 2) - mean(values[:, previous]; dims = 2))
    drift = mean(per_sample)
    stderr = samples_num == 1 ? 0.0 : std(per_sample) / sqrt(samples_num)
    return (; drift, stderr, previous, final)
end


function _samples_collect_process_data_ising(
    L::Integer,
    τ_idx::Integer,
    t::Integer,
    initial_state_label::AbstractString,
    data_root::AbstractString,
)
    L, τ_idx, t = _validate_ising_run(L, τ_idx, t)
    dir_path = joinpath(data_root, "L$(L)", "gammaind$(τ_idx)")
    isdir(dir_path) || error("Data directory does not exist: $dir_path")
    existing_files = sort(filter(
        file -> startswith(file, "t$(t)_samples") && endswith(file, ".jld"),
        readdir(dir_path),
    ))
    samples_num = length(existing_files)
    samples_num > 0 || error(
        "No trajectory files with L=$L, τ_idx=$τ_idx, t=$t in $dir_path",
    )
    println("collecting $samples_num $(initial_state_label) sample files")

    n_layers = FibonacciChain.layers_per_period(ising_model(L))
    n_layers == 2 || error("This collector assumes the two-layer Ising circuit")
    ensemble_seed = zeros(Int, samples_num)
    ensemble_EE_dynamics = zeros(Float64, samples_num, t)
    ensemble_final_EElis = zeros(Float64, samples_num, L - 1)
    ensemble_FE_dynamics = zeros(Float64, samples_num, t)

    for (i, fname) in enumerate(existing_files)
        data = load(joinpath(dir_path, fname))
        get(data, "initial_state", initial_state_label) == initial_state_label ||
            error("Unexpected initial state in $fname")
        Int(get(data, "L", L)) == L || error("Inconsistent L in $fname")
        Int(get(data, "τ_idx", τ_idx)) == τ_idx ||
            error("Inconsistent τ_idx in $fname")
        Int(get(data, "t", t)) == t || error("Inconsistent period count in $fname")

        sample = data["sample"]
        sample_free_energy = Float64.(data["sample_free_energy"])
        halfchain_EE_tlis = Float64.(data["halfchain_EE_tlis"])
        final_EElis = Float64.(data["final_EElis"])
        size(sample, 1) == n_layers * t ||
            error("Inconsistent sample depth in $fname")
        length(sample_free_energy) == n_layers * t ||
            error("Inconsistent free-energy length in $fname")
        length(halfchain_EE_tlis) == t ||
            error("Inconsistent entanglement time-series length in $fname")
        length(final_EElis) == L - 1 ||
            error("Inconsistent final entropy profile in $fname")

        ensemble_seed[i] = Int(data["seed"])
        ensemble_EE_dynamics[i, :] = halfchain_EE_tlis
        ensemble_final_EElis[i, :] = final_EElis
        # Keep the free-energy normalization per layer while using period as
        # the time coordinate.
        ensemble_FE_dynamics[i, :] =
            (sample_free_energy[1:2:end] + sample_free_energy[2:2:end]) / 2
    end
    check_duplicates(ensemble_seed)

    average_EE_tlis, stderr_EE_tlis = _mean_and_stderr(ensemble_EE_dynamics)
    bulk_meanEElis, ensemble_stderr_EElis = _mean_and_stderr(ensemble_final_EElis)
    time_FElis, time_FEstderr = _mean_and_stderr(ensemble_FE_dynamics)

    _, _, configured_periods = get_cfg_params_Born(τ_idx, L)
    averaging_periods = collect(configured_periods)
    isempty(averaging_periods) && error("The configured averaging window is empty")
    first(averaging_periods) >= 1 && last(averaging_periods) <= t ||
        error("The configured averaging window lies outside 1:t")
    time_average_free_energy = vec(mean(
        ensemble_FE_dynamics[:, averaging_periods];
        dims = 2,
    ))
    bulk_FE = mean(time_average_free_energy)
    bulk_FE_stderr = samples_num == 1 ?
                     0.0 : std(time_average_free_energy) / sqrt(samples_num)

    ee_drift = late_time_drift(ensemble_EE_dynamics)
    free_energy_drift = late_time_drift(ensemble_FE_dynamics)
    output_path = joinpath(
        data_root,
        "L$(L)",
        "EE_FEdynamics_L$(L)_gamma$(τ_idx)_t$(t).jld2",
    )
    save(
        output_path,
        "initial_state",
        String(initial_state_label),
        "L",
        L,
        "τ_idx",
        τ_idx,
        "τ",
        τlis[τ_idx],
        "gamma",
        tanh(τlis[τ_idx]),
        "t",
        t,
        "layers_per_period",
        n_layers,
        "samples_num",
        samples_num,
        "average_EE_tlis",
        average_EE_tlis,
        "stderr_EE_tlis",
        stderr_EE_tlis,
        "bulk_meanEElis",
        bulk_meanEElis,
        "ensemble_stderr_EElis",
        ensemble_stderr_EElis,
        "averaging_periods",
        averaging_periods,
        "time_average_free_energy",
        time_average_free_energy,
        "bulk_FE",
        bulk_FE,
        "bulk_FE_stderr",
        bulk_FE_stderr,
        "time_FEstderr",
        time_FEstderr,
        "time_FElis",
        time_FElis,
        "ensemble_seed",
        ensemble_seed,
        "late_time_EE_drift",
        ee_drift.drift,
        "late_time_EE_drift_stderr",
        ee_drift.stderr,
        "late_time_FE_drift",
        free_energy_drift.drift,
        "late_time_FE_drift_stderr",
        free_energy_drift.stderr,
        "late_time_previous_periods",
        collect(ee_drift.previous),
        "late_time_final_periods",
        collect(ee_drift.final),
    )
    return output_path
end

"""Collect `|+⟩^⊗L` trajectories; `t` is the number of periods."""
function samples_collect_process_data(
    L::Integer,
    τ_idx::Integer,
    t::Integer = get_cfg_params_Born(τ_idx, L)[1];
    data_root::AbstractString = ISING_DATA_ROOT,
)
    return _samples_collect_process_data_ising(
        L,
        τ_idx,
        t,
        "all_plus",
        data_root,
    )
end

"""Collect topological-sector trajectories; `t` is the number of periods."""
function samples_collect_process_data_topo_sector(
    L::Integer,
    τ_idx::Integer,
    t::Integer = get_cfg_params_Born(τ_idx, L)[1];
    data_root::AbstractString = ISING_TOPO_SECTOR_DATA_ROOT,
)
    return _samples_collect_process_data_ising(
        L,
        τ_idx,
        t,
        "topo_sector",
        data_root,
    )
end


function print_usage()
    println("Here t is the number of complete measurement periods.")
    println("Usage:")
    println("  julia --project=. exm/Born_Ising/moni_dyna_Ising.jl L τ_idx seed [t]")
    println("  julia --project=. exm/Born_Ising/moni_dyna_Ising.jl collect L τ_idx [t]")
    println("  julia --project=. exm/Born_Ising/moni_dyna_Ising.jl topo_sector L τ_idx seed [t]")
    println("  julia --project=. exm/Born_Ising/moni_dyna_Ising.jl collect_topo_sector L τ_idx [t]")
end

function _configured_or_explicit_time(args, index::Int, τ_idx::Int, L::Int)
    return length(args) >= index ? parse(Int, args[index]) : get_cfg_params_Born(τ_idx, L)[1]
end

if isempty(ARGS)
    print_usage()
elseif lowercase(ARGS[1]) == "collect"
    length(ARGS) in (3, 4) || error("Usage: collect L τ_idx [t]")
    L = parse(Int, ARGS[2])
    τ_idx = parse(Int, ARGS[3])
    t = _configured_or_explicit_time(ARGS, 4, τ_idx, L)
    println("saved: $(samples_collect_process_data(L, τ_idx, t))")
elseif lowercase(ARGS[1]) == "topo_sector"
    length(ARGS) in (4, 5) || error("Usage: topo_sector L τ_idx seed [t]")
    L = parse(Int, ARGS[2])
    τ_idx = parse(Int, ARGS[3])
    seed = parse(Int, ARGS[4])
    t = _configured_or_explicit_time(ARGS, 5, τ_idx, L)
    println("saved: $(samples_generate_topo_sector(L, τ_idx, seed, t))")
elseif lowercase(ARGS[1]) == "collect_topo_sector"
    length(ARGS) in (3, 4) || error("Usage: collect_topo_sector L τ_idx [t]")
    L = parse(Int, ARGS[2])
    τ_idx = parse(Int, ARGS[3])
    t = _configured_or_explicit_time(ARGS, 4, τ_idx, L)
    println("saved: $(samples_collect_process_data_topo_sector(L, τ_idx, t))")
else
    length(ARGS) in (3, 4) || error("Usage: L τ_idx seed [t]")
    L = parse(Int, ARGS[1])
    τ_idx = parse(Int, ARGS[2])
    seed = parse(Int, ARGS[3])
    t = _configured_or_explicit_time(ARGS, 4, τ_idx, L)
    println("saved: $(samples_generate(L, τ_idx, seed, t))")
end
