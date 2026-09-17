using FibonacciChain
using ITensorMPS
using ITensors
using JLD2
using LinearAlgebra
using Random
using Statistics

# Sector-resolved Lyapunov spectra along Born trajectories (y=1 or y=τ).
# Merged from lyapunov_spectrum_y1.jl and lyapunov_spectrum_ytau.jl;
# existing data roots, file names and numerical JLD2 fields are preserved.
# New runs mark initial_frame_orthonormalized=true. Older MPS files lack this
# flag and include a sector-dependent initial-frame contribution to λ(t).
#
# Shared γ/τ grid (τlis), fib_model, and the evolution-time tables
# (get_cfg_params_Born / get_mps_params_Born); measurement strengths are
# selected by integer index τ_idx, as in transfer_matrix.jl
include(joinpath(@__DIR__, "config.jl"))

# In the repository normalization the y=1 sector has topological charge
# eigenvalue ϕ, while the y=τ sector has eigenvalue -1/ϕ.
const SECTOR_SPECS = Dict{Symbol,NamedTuple}(
    :y1 => (;
        title = "y=1",
        sector = :trivial,
        dir_name = "lyapunov_spectrum_y1",
        file_tag = "lyapunov_y1",
        sector_label = "y=1",
        y_eigenvalue = (1 + √5) / 2,
        other_eigenvalue = -2 / (1 + √5),
        min_χ = 1,
    ),
    :ytau => (;
        title = "y=τ",
        sector = :tau,
        dir_name = "lyapunov_spectrum_ytau",
        file_tag = "lyapunov_ytau",
        sector_label = "y=tau",
        y_eigenvalue = -2 / (1 + √5),
        other_eigenvalue = (1 + √5) / 2,
        # The y=τ projected all-τ MPS needs bond dimension at least 5.
        min_χ = 5,
    ),
)

function sector_spec(sector::Symbol)
    haskey(SECTOR_SPECS, sector) ||
        error("Unknown sector: $sector (expected :y1 or :ytau)")
    return SECTOR_SPECS[sector]
end

"""
    initial_sector_state(model, spec) -> Vector

Project the all-`τ` fusion path into the sector selected by `spec` and
normalize it. This is the state used to sample the Born trajectory. With `y`
the sector eigenvalue and `ȳ` the other one, `P = (Y - ȳ I)/(y - ȳ)`.
"""
function initial_sector_state(model::AnyonModel{FibonacciAnyon}, spec)
    model.pbc || error("The topological $(spec.title) sector requires periodic boundaries")
    Y = topological_charge_operator(model)
    P = (Y - spec.other_eigenvalue * I(size(Y, 1))) /
        (spec.y_eigenvalue - spec.other_eigenvalue)
    state = zeros(Float64, length(anyon_basis(model)))
    state[1] = 1.0
    projected = P * state
    weight = real(dot(projected, projected))
    weight > 1e-14 || error("The all-τ path has zero $(spec.title) weight")
    return projected / sqrt(weight)
end

"""
    initial_sector_mps(model, spec; cutoff=1e-14, maxdim=32)

MPS version of [`initial_sector_state`](@ref): apply the topological charge
MPO to the all-`τ` product state and add the identity part,

    P|00⋯0⟩ ∝ (Y - ȳ I)|00⋯0⟩,

so no dense Y matrix is ever constructed and large L is accessible.
"""
function initial_sector_mps(
    model::AnyonModel{FibonacciAnyon},
    spec;
    cutoff::Float64 = 1e-14,
    maxdim::Int = 32,
)
    model.pbc || error("The topological $(spec.title) sector requires periodic boundaries")
    N = model.N
    N >= 2 || error("initial_sector_mps requires at least two sites")
    sites = siteinds("Qubit", N)

    vacuum = productMPS(sites, fill("0", N))
    Y_mpo = topological_charge_mpo(sites; pbc = true)
    y_on_vacuum = apply(Y_mpo, vacuum; cutoff = cutoff, maxdim = maxdim)

    vacuum[1] = -spec.other_eigenvalue * vacuum[1]
    projected = add(y_on_vacuum, vacuum; cutoff = cutoff, maxdim = maxdim)
    normalize!(projected)
    return projected, sites
end

"""
    lyapunov_sector_time(backend, L, τ_idx) -> Int

Evolution time `t` (number of measurement periods, two layers each) taken from
the shared tables in `config.jl`, as in `transfer_matrix.jl`.  Both config
helpers return this period count directly.
"""
function lyapunov_sector_time(backend::Symbol, L::Integer, τ_idx::Integer)
    if backend == :exact
        t, _, _ = get_cfg_params_Born(τ_idx, L)
        return t
    elseif backend == :mps
        time_in_L, _, _ = get_mps_params_Born(τ_idx, L)
        return time_in_L * L
    end
    error("backend must be :exact or :mps")
end

"""
    simulate_sector_lyapunov(spec, L, τ, t; n_states=10, trajectory_seed=1)

Generate a Born trajectory of `t` periods from the sector-projected all-`τ`
state (exact state-vector evolution), then compute the Lyapunov spectrum of
the corresponding transfer-matrix product restricted to the full sector via
`lyapunov_spectrum_topological_sector`.
"""
function simulate_sector_lyapunov(
    spec,
    L::Int,
    τ::Real,
    t::Int;
    n_states::Int = 10,
    trajectory_seed::Int = 1,
)
    t >= 1 || throw(ArgumentError("t must be positive"))
    model = fib_model(L)
    initial_state = initial_sector_state(model, spec)
    born_config = MeasureConfig(
        τ = Float64(τ),
        t₂ = t,
        mode = :Born,
        rng = MersenneTwister(trajectory_seed),
        enable_τ_eff = false,
        track_y_expectation = true,
    )
    trajectory = bulk_evolution(model, initial_state, born_config)
    spectrum = lyapunov_spectrum_topological_sector(
        model,
        Float64(τ),
        trajectory.samples;
        sector = spec.sector,
        n_states = n_states,
    )

    maximum(abs.(trajectory.y_expectation_values .- spec.y_eigenvalue)) < 1e-5 ||
        error("the Born trajectory left the $(spec.title) sector")
    maximum(spectrum.sector_leakage) < 1e-9 ||
        error("the Lyapunov frame left the $(spec.title) sector")
    return trajectory, (
        local_log_stretches = spectrum.local_log_stretches,
        lyapunov_exponents = spectrum.lyapunov_exponents,
        free_energy_spectrum = spectrum.free_energy_spectrum,
        sector_dimension = spectrum.sector_dimension,
        sector_leakage = spectrum.sector_leakage,
    )
end

"""
    simulate_sector_lyapunov_mps(spec, L, τ, t, χ; n_states=10, trajectory_seed=1)

MPS version of [`simulate_sector_lyapunov`](@ref) for large L. The Born
trajectory of `t` periods is generated from the Y-MPO-projected initial state
(`initial_sector_mps`), with the sector monitored on-the-fly by
`topological_charge_mpo`. The Lyapunov spectrum is computed by
`lyapunov_spectrum_mps` with `sector = spec.sector`, which projects the
propagated frame back into the sector after every period via the Y MPO; since
the transfer matrix commutes with `Y`, this only removes MPS truncation
leakage into the other sector and the unphysical fusion-path subspace.
"""
function simulate_sector_lyapunov_mps(
    spec,
    L::Int,
    τ::Real,
    t::Int,
    χ::Int;
    n_states::Int = 10,
    trajectory_seed::Int = 1,
)
    t >= 1 || throw(ArgumentError("t must be positive"))
    χ >= spec.min_χ ||
        error("χ must be at least $(spec.min_χ) to represent the $(spec.title) initial MPS")
    model = fib_model(L)
    initial_state, sites = initial_sector_mps(model, spec; cutoff = 1e-12, maxdim = χ)
    born_config = MeasureConfig(
        τ = Float64(τ),
        t₂ = t,
        mode = :Born,
        rng = MersenneTwister(trajectory_seed),
        enable_τ_eff = false,
        cutoff = 1e-12,
        maxdim = χ,
        track_y_expectation = true,
    )
    trajectory = bulk_evolution(model, sites, initial_state, born_config)

    spectrum_mps = lyapunov_spectrum_mps(
        model,
        sites,
        Float64(τ),
        trajectory.samples;
        n_states = n_states,
        sector = spec.sector,
        cutoff = 1e-12,
        maxdim = χ,
    )
    local_log_stretches = -spectrum_mps
    t_out = size(local_log_stretches, 2)
    lyapunov_exponents =
        cumsum(local_log_stretches; dims = 2) ./ reshape(collect(1:t_out), 1, :)

    maximum(abs.(trajectory.y_expectation_values .- spec.y_eigenvalue)) < 1e-4 ||
        error("the MPS Born trajectory left the $(spec.title) sector")
    return trajectory, (
        local_log_stretches = local_log_stretches,
        lyapunov_exponents = lyapunov_exponents,
        free_energy_spectrum = -lyapunov_exponents,
        sector_dimension = 0,  # not available without the dense Y
        sector_leakage = Float64[],  # tracked instead via y_expectation_values
    )
end


"""
    lyapunov_sector_dir(spec, backend, L, τ_idx; χ=nothing) -> String

Directory holding the per-seed data files of this script, following the
`transfer_matrix.jl` layout: `exm/data/Bulk_measure/<spec.dir_name>/
L\$(L)/gammaind\$(τ_idx)` for the exact backend, with an additional
`chi\$(χ)` level for the MPS backend. Paths are relative to the repository
root.
"""
function lyapunov_sector_dir(
    spec,
    backend::Symbol,
    L::Integer,
    τ_idx::Integer;
    χ::Union{Nothing,Integer} = nothing,
)
    backend in (:exact, :mps) || throw(ArgumentError("backend must be :exact or :mps"))
    backend == :mps && χ === nothing &&
        throw(ArgumentError("χ is required for the MPS backend"))
    path = joinpath(
        "exm",
        "data",
        "Bulk_measure",
        spec.dir_name,
        "L$(L)",
        "gammaind$(τ_idx)",
    )
    return backend == :mps ? joinpath(path, "chi$(χ)") : path
end

"""
    collect_lyapunov_sector(spec, backend, L, τ_idx; χ=nothing)

Aggregate the per-seed sector Lyapunov spectra produced by the `exact`/`mps`
modes of this script into a single ensemble file: element-wise mean and
standard error (over seeds) of `local_log_stretches`, `lyapunov_exponents`,
`free_energy_spectrum` and `y_expectation_values`, plus the per-seed final
exponents. The evolution time `t` is taken from [`lyapunov_sector_time`](@ref).
Following `transfer_matrix.jl`, the ensemble file is saved one level above the
per-seed directory, in `.../<spec.dir_name>/L\$(L)`. Returns the output path.
"""
function collect_lyapunov_sector(
    spec,
    backend::Symbol,
    L::Integer,
    τ_idx::Integer;
    χ::Union{Nothing,Integer} = nothing,
)
    backend in (:exact, :mps) || throw(ArgumentError("backend must be :exact or :mps"))
    backend == :mps && χ === nothing &&
        throw(ArgumentError("collecting the mps backend requires χ"))
    t = lyapunov_sector_time(backend, L, τ_idx)
    data_dir = lyapunov_sector_dir(spec, backend, L, τ_idx; χ = χ)
    files = sort(filter(
        file -> startswith(file, "$(spec.file_tag)_L$(L)_t$(t)_seed") &&
                endswith(file, ".jld2"),
        readdir(data_dir),
    ))
    isempty(files) &&
        error("No $backend files with L=$L, τ_idx=$τ_idx, t=$t found in $data_dir")

    first_data = JLD2.load(joinpath(data_dir, first(files)))
    Int(first_data["t"]) == t || error("Inconsistent evolution length")
    n_states = Int(first_data["n_states"])
    initial_frame_orthonormalized = get(
        first_data, "initial_frame_orthonormalized", backend == :exact,
    )

    samples_num = length(files)
    seeds = zeros(Int, samples_num)
    log_stretches = Vector{Matrix{Float64}}(undef, samples_num)
    exponents = Vector{Matrix{Float64}}(undef, samples_num)
    free_energies = Vector{Matrix{Float64}}(undef, samples_num)
    y_expectations = Vector{Vector{Float64}}(undef, samples_num)
    sector_leakages = Vector{Vector{Float64}}(undef, samples_num)

    for (i, file) in enumerate(files)
        data = JLD2.load(joinpath(data_dir, file))
        String(data["backend"]) == String(backend) ||
            error("Unexpected backend in $file")
        data["topological_sector"] == spec.sector_label ||
            error("Unexpected sector in $file")
        Int(data["L"]) == L && Int(data["t"]) == t ||
            error("Inconsistent system size or evolution length in $file")
        Int(data["τ_idx"]) == τ_idx || error("Inconsistent τ_idx in $file")
        Int(data["n_states"]) == n_states || error("Inconsistent n_states in $file")
        get(data, "initial_frame_orthonormalized", backend == :exact) ==
            initial_frame_orthonormalized || error(
            "Cannot mix legacy and corrected initial-frame normalization in $file",
        )
        seeds[i] = Int(data["trajectory_seed"])
        log_stretches[i] = Float64.(data["local_log_stretches"])
        exponents[i] = Float64.(data["lyapunov_exponents"])
        free_energies[i] = Float64.(data["free_energy_spectrum"])
        y_expectations[i] = Float64.(data["y_expectation_values"])
        sector_leakages[i] = Float64.(get(data, "sector_leakage", Float64[]))
    end
    length(unique(seeds)) == samples_num || error("Duplicate trajectory seeds found")
    all(size(m) == (n_states, t) for m in exponents) ||
        error("Inconsistent spectrum shapes")
    all(length(y) == t for y in y_expectations) ||
        error("Inconsistent y_expectation lengths")

    stderr_of(mats) = std(mats; corrected = false) ./ √samples_num
    final_exponents = hcat([m[:, end] for m in exponents]...)
    y_matrix = reduce(hcat, y_expectations)  # t × samples_num

    out_dir = joinpath("exm", "data", "Bulk_measure", spec.dir_name, "L$(L)")
    mkpath(out_dir)
    output_path = joinpath(
        out_dir,
        backend == :mps ?
        "ensemble_$(spec.file_tag)_L$(L)_gammaind$(τ_idx)_t$(t)_chi$(χ).jld2" :
        "ensemble_$(spec.file_tag)_L$(L)_gammaind$(τ_idx)_t$(t).jld2",
    )
    jldsave(
        output_path;
        backend = String(backend),
        initial_frame_orthonormalized = initial_frame_orthonormalized,
        topological_sector = spec.sector_label,
        y_eigenvalue = spec.y_eigenvalue,
        L = Int(L),
        τ_idx = Int(τ_idx),
        τ = τlis[τ_idx],
        gamma = tanh(τlis[τ_idx]),
        χ = isnothing(χ) ? 0 : Int(χ),
        t = t,
        n_states = n_states,
        samples_num = samples_num,
        ensemble_seeds = seeds,
        mean_local_log_stretches = mean(log_stretches),
        stderr_local_log_stretches = stderr_of(log_stretches),
        mean_lyapunov_exponents = mean(exponents),
        stderr_lyapunov_exponents = stderr_of(exponents),
        mean_free_energy_spectrum = mean(free_energies),
        stderr_free_energy_spectrum = stderr_of(free_energies),
        mean_y_expectation = vec(mean(y_matrix; dims = 2)),
        stderr_y_expectation =
            vec(std(y_matrix; dims = 2, corrected = false)) ./ √samples_num,
        final_exponents_per_seed = final_exponents,
        mean_final_exponents = vec(mean(final_exponents; dims = 2)),
        stderr_final_exponents = vec(std(final_exponents; dims = 2, corrected = false)) ./
                                 √samples_num,
        sector_leakages = sector_leakages,
    )
    println("collected $samples_num seed files from $data_dir")
    return output_path
end


# =============================================================================
# Shell interface
# =============================================================================
if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia --project=. exm/Bulk_measure/lyapunov_spectrum_sector.jl sector backend [args...]")
    println("")
    println("sector is y1 or ytau. Backends (measurement strength by index τ_idx")
    println("into config.jl's τlis; evolution time t in periods from")
    println("get_cfg_params_Born/get_mps_params_Born):")
    println("")
    println("  exact L τ_idx n_states seed [output.jld2]")
    println("      Born trajectory (exact state vector) from the sector-projected")
    println("      all-τ state, plus the sector Lyapunov spectrum of the")
    println("      transfer-matrix product along it.")
    println("")
    println("  mps L τ_idx chi n_states seed [output.jld2]")
    println("      MPS version for large L; the sector is enforced and monitored")
    println("      through the topological-charge MPO.")
    println("")
    println("  collect exact L τ_idx")
    println("  collect mps L τ_idx chi")
    println("      Aggregate the per-seed spectra into an ensemble file in the L directory.")
    println("")
    println("Examples:")
    println("  julia --project=. exm/Bulk_measure/lyapunov_spectrum_sector.jl y1 exact 8 10 3 1")
    println("  julia --project=. exm/Bulk_measure/lyapunov_spectrum_sector.jl ytau collect exact 8 10")
else
    length(ARGS) >= 2 || error("Usage: lyapunov_spectrum_sector.jl sector backend [args...]")
    spec = sector_spec(Symbol(lowercase(ARGS[1])))
    backend = Symbol(lowercase(ARGS[2]))

    if backend == :collect
        # -------------------------------------------------------------------
        # collect: ensemble-average the per-seed spectra
        # -------------------------------------------------------------------
        length(ARGS) >= 3 || error("collect requires a sub-backend: exact or mps")
        sub = Symbol(lowercase(ARGS[3]))
        if sub == :exact
            length(ARGS) == 5 || error("Usage: SECTOR collect exact L τ_idx")
            output_path = collect_lyapunov_sector(
                spec,
                :exact,
                parse(Int, ARGS[4]),
                parse(Int, ARGS[5]),
            )
        elseif sub == :mps
            length(ARGS) == 6 || error("Usage: SECTOR collect mps L τ_idx chi")
            output_path = collect_lyapunov_sector(
                spec,
                :mps,
                parse(Int, ARGS[4]),
                parse(Int, ARGS[5]);
                χ = parse(Int, ARGS[6]),
            )
        else
            error("Unknown collect backend: $sub")
        end
        println("saved: $output_path")

    elseif backend == :exact || backend == :mps
        # -------------------------------------------------------------------
        # exact/mps: single Born trajectory + sector Lyapunov spectrum
        # -------------------------------------------------------------------
        if backend == :exact
            length(ARGS) in (6, 7) ||
                error("Usage: SECTOR exact L τ_idx n_states seed [output.jld2]")
            L = parse(Int, ARGS[3])
            τ_idx = parse(Int, ARGS[4])
            n_states = parse(Int, ARGS[5])
            seed = parse(Int, ARGS[6])
            χ = 0
        else
            length(ARGS) in (7, 8) ||
                error("Usage: SECTOR mps L τ_idx chi n_states seed [output.jld2]")
            L = parse(Int, ARGS[3])
            τ_idx = parse(Int, ARGS[4])
            χ = parse(Int, ARGS[5])
            n_states = parse(Int, ARGS[6])
            seed = parse(Int, ARGS[7])
        end
        τ = τlis[τ_idx]
        t = lyapunov_sector_time(backend, L, τ_idx)

        println("=== $(spec.title) Lyapunov spectrum ($backend backend) ===")
        println("L = $L, τ_idx = $τ_idx, τ = $τ, γ = $(tanh(τ))" *
                (backend == :mps ? ", χ = $χ" : ""))
        println("t = $t, n_states = $n_states, seed = $seed")

        trajectory, spectrum = if backend == :exact
            simulate_sector_lyapunov(spec, L, τ, t; n_states = n_states, trajectory_seed = seed)
        else
            simulate_sector_lyapunov_mps(
                spec, L, τ, t, χ; n_states = n_states, trajectory_seed = seed,
            )
        end

        has_custom_output = (backend == :exact && length(ARGS) == 7) ||
                            (backend == :mps && length(ARGS) == 8)
        output_path = has_custom_output ? ARGS[end] : joinpath(
            lyapunov_sector_dir(spec, backend, L, τ_idx; χ = backend == :mps ? χ : nothing),
            "$(spec.file_tag)_L$(L)_t$(t)_seed$(seed).jld2",
        )
        mkpath(dirname(output_path))
        jldsave(
            output_path;
            backend = String(backend),
            initial_frame_orthonormalized = true,
            L = L,
            τ_idx = τ_idx,
            τ = τ,
            gamma = tanh(τ),
            χ = χ,
            t = t,
            n_states = n_states,
            trajectory_seed = seed,
            topological_sector = spec.sector_label,
            y_eigenvalue = spec.y_eigenvalue,
            sector_dimension = spectrum.sector_dimension,
            samples = trajectory.samples,
            local_log_stretches = spectrum.local_log_stretches,
            lyapunov_exponents = spectrum.lyapunov_exponents,
            free_energy_spectrum = spectrum.free_energy_spectrum,
            sector_leakage = spectrum.sector_leakage,
            y_expectation_values = trajectory.y_expectation_values,
        )

        println("saved: $output_path")
        spectrum.sector_dimension > 0 &&
            println("$(spec.title) sector dimension: $(spectrum.sector_dimension)")
        !isempty(spectrum.sector_leakage) &&
            println("max sector leakage: $(maximum(spectrum.sector_leakage))")
        println("final Lyapunov exponents: $(spectrum.lyapunov_exponents[:, end])")
        println("final free-energy spectrum: $(spectrum.free_energy_spectrum[:, end])")
    else
        error("Unknown backend: $backend")
    end
end
