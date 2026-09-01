using FibonacciChain
using ITensorMPS
using ITensors
using JLD2
using LinearAlgebra
using Random
using Statistics

# Shared γ/τ grid (τlis), fib_model, and the evolution-time tables
# (get_cfg_params_Born / get_mps_params_Born); measurement strengths are
# selected by integer index τ_idx, as in transfer_matrix.jl
include(joinpath(@__DIR__, "config.jl"))

"""
    initial_ytau_state(model) -> Vector

Project the all-`τ` fusion path into the `y=τ` sector and normalize it. This is
the state used to sample the Born trajectory. In the repository normalization
the `y=1` sector has topological charge eigenvalue `ϕ`, while the `y=τ` sector
has eigenvalue `-1/ϕ`, so `Pτ = (ϕ I - Y) / (ϕ + ϕ⁻¹)`.
"""
function initial_ytau_state(model::AnyonModel{FibonacciAnyon})
    model.pbc || error("The topological y=τ sector requires periodic boundaries")
    ϕ = (1 + √5) / 2
    Y = topological_charge_operator(model)
    Pτ = (ϕ * I(size(Y, 1)) - Y) / (ϕ + inv(ϕ))
    state = zeros(Float64, length(anyon_basis(model)))
    state[1] = 1.0
    projected = Pτ * state
    weight = real(dot(projected, projected))
    weight > 1e-14 || error("The all-τ path has zero y=τ weight")
    return projected / sqrt(weight)
end

"""
    initial_ytau_mps(model; cutoff=1e-14, maxdim=32)

MPS version of [`initial_ytau_state`](@ref): apply the topological charge MPO to
the all-`τ` product state and add the identity part,

    Pτ|00⋯0⟩ ∝ (ϕ I - Y)|00⋯0⟩,

so no dense Y matrix is ever constructed and large L is accessible.
"""
function initial_ytau_mps(model::AnyonModel{FibonacciAnyon}; cutoff::Float64 = 1e-14, maxdim::Int = 32)
    model.pbc || error("The topological y=τ sector requires periodic boundaries")
    N = model.N
    N >= 2 || error("initial_ytau_mps requires at least two sites")
    sites = siteinds("Qubit", N)

    vacuum = productMPS(sites, fill("0", N))
    Y_mpo = topological_charge_mpo(sites; pbc = true)
    y_on_vacuum = apply(Y_mpo, vacuum; cutoff = cutoff, maxdim = maxdim)

    ϕ = (1 + √5) / 2
    y_on_vacuum[1] = -y_on_vacuum[1]
    vacuum[1] = ϕ * vacuum[1]
    projected = add(y_on_vacuum, vacuum; cutoff = cutoff, maxdim = maxdim)
    normalize!(projected)
    return projected, sites
end

"""
    lyapunov_ytau_time(backend, L, τ_idx) -> Int

Evolution time `t` (number of measurement periods, two layers each) taken from
the shared tables in `config.jl`, as in `transfer_matrix.jl`: the exact backend
uses `div(D, 2)` from `get_cfg_params_Born`, the MPS backend uses
`time_in_L * L` from `get_mps_params_Born`.
"""
function lyapunov_ytau_time(backend::Symbol, L::Integer, τ_idx::Integer)
    if backend == :exact
        D, _, _ = get_cfg_params_Born(τ_idx, L)
        return div(D, 2)
    elseif backend == :mps
        time_in_L, _, _ = get_mps_params_Born(τ_idx, L)
        return time_in_L * L
    end
    error("backend must be :exact or :mps")
end

"""
    simulate_ytau_lyapunov(L, τ, t; n_states=10, trajectory_seed=1)

Generate a Born trajectory of `t` periods from `Pτ|τ⋯τ⟩` (exact state-vector
evolution), then compute the Lyapunov spectrum of the corresponding
transfer-matrix product restricted to the full `y=τ` sector via
`lyapunov_spectrum_topological_sector`.
"""
function simulate_ytau_lyapunov(
    L::Int,
    τ::Real,
    t::Int;
    n_states::Int = 10,
    trajectory_seed::Int = 1,
)
    t >= 1 || throw(ArgumentError("t must be positive"))
    model = fib_model(L)
    initial_state = initial_ytau_state(model)
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
        sector = :tau,
        n_states = n_states,
    )

    ϕ = (1 + √5) / 2
    maximum(abs.(trajectory.y_expectation_values .+ inv(ϕ))) < 1e-5 ||
        error("the Born trajectory left the y=τ sector")
    maximum(spectrum.sector_leakage) < 1e-9 ||
        error("the Lyapunov frame left the y=τ sector")
    return trajectory, (
        local_log_stretches = spectrum.local_log_stretches,
        lyapunov_exponents = spectrum.lyapunov_exponents,
        free_energy_spectrum = spectrum.free_energy_spectrum,
        sector_dimension = spectrum.sector_dimension,
        sector_leakage = spectrum.sector_leakage,
    )
end

"""
    simulate_ytau_lyapunov_mps(L, τ, t, χ; n_states=10, trajectory_seed=1)

MPS version of [`simulate_ytau_lyapunov`](@ref) for large L. The Born trajectory
of `t` periods is generated from the Y-MPO-projected initial state
(`initial_ytau_mps`), with the `y=τ` sector monitored on-the-fly by
`topological_charge_mpo`. The Lyapunov spectrum is computed by
`lyapunov_spectrum_mps` with `sector = :tau`, which projects the propagated
frame back into the `y=τ` sector after every period via the Y MPO; since the
transfer matrix commutes with `Y`, this only removes MPS truncation leakage
into the `y=1` sector and the unphysical fusion-path subspace.
"""
function simulate_ytau_lyapunov_mps(
    L::Int,
    τ::Real,
    t::Int,
    χ::Int;
    n_states::Int = 10,
    trajectory_seed::Int = 1,
)
    t >= 1 || throw(ArgumentError("t must be positive"))
    χ >= 5 || error("χ must be at least 5 to represent the y=τ initial MPS")
    model = fib_model(L)
    initial_state, sites = initial_ytau_mps(model; cutoff = 1e-12, maxdim = χ)
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
        sector = :tau,
        cutoff = 1e-12,
        maxdim = χ,
    )
    local_log_stretches = -spectrum_mps
    t_out = size(local_log_stretches, 2)
    lyapunov_exponents =
        cumsum(local_log_stretches; dims = 2) ./ reshape(collect(1:t_out), 1, :)

    ϕ = (1 + √5) / 2
    maximum(abs.(trajectory.y_expectation_values .+ inv(ϕ))) < 1e-4 ||
        error("the MPS Born trajectory left the y=τ sector")
    return trajectory, (
        local_log_stretches = local_log_stretches,
        lyapunov_exponents = lyapunov_exponents,
        free_energy_spectrum = -lyapunov_exponents,
        sector_dimension = 0,  # not available without the dense Y
        sector_leakage = Float64[],  # tracked instead via y_expectation_values
    )
end


"""
    lyapunov_ytau_dir(backend, L, τ_idx; χ=nothing) -> String

Directory holding the per-seed data files of this script, following the
`transfer_matrix.jl` layout: `exm/data/Bulk_measure/lyapunov_spectrum_ytau/
L\$(L)/gammaind\$(τ_idx)` for the exact backend, with an additional
`chi\$(χ)` level for the MPS backend. Paths are relative to the repository
root.
"""
function lyapunov_ytau_dir(
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
        "lyapunov_spectrum_ytau",
        "L$(L)",
        "gammaind$(τ_idx)",
    )
    return backend == :mps ? joinpath(path, "chi$(χ)") : path
end

"""
    collect_lyapunov_ytau(backend, L, τ_idx; χ=nothing)

Aggregate the per-seed `y=τ` Lyapunov spectra produced by the `exact`/`mps`
modes of this script into a single ensemble file: element-wise mean and
standard error (over seeds) of `local_log_stretches`, `lyapunov_exponents`,
`free_energy_spectrum` and `y_expectation_values`, plus the per-seed final
exponents. The evolution time `t` is taken from [`lyapunov_ytau_time`](@ref).
Following `transfer_matrix.jl`, the ensemble file is saved one level above the
per-seed directory, in `.../lyapunov_spectrum_ytau/L\$(L)`. Returns the output
path.
"""
function collect_lyapunov_ytau(
    backend::Symbol,
    L::Integer,
    τ_idx::Integer;
    χ::Union{Nothing,Integer} = nothing,
)
    backend in (:exact, :mps) || throw(ArgumentError("backend must be :exact or :mps"))
    backend == :mps && χ === nothing &&
        throw(ArgumentError("collecting the mps backend requires χ"))
    t = lyapunov_ytau_time(backend, L, τ_idx)
    data_dir = lyapunov_ytau_dir(backend, L, τ_idx; χ = χ)
    files = sort(filter(
        file -> startswith(file, "lyapunov_ytau_L$(L)_t$(t)_seed") &&
                endswith(file, ".jld2"),
        readdir(data_dir),
    ))
    isempty(files) &&
        error("No $backend files with L=$L, τ_idx=$τ_idx, t=$t found in $data_dir")

    first_data = JLD2.load(joinpath(data_dir, first(files)))
    Int(first_data["t"]) == t || error("Inconsistent evolution length")
    n_states = Int(first_data["n_states"])

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
        data["topological_sector"] == "y=tau" || error("Unexpected sector in $file")
        Int(data["L"]) == L && Int(data["t"]) == t ||
            error("Inconsistent system size or evolution length in $file")
        Int(data["τ_idx"]) == τ_idx || error("Inconsistent τ_idx in $file")
        Int(data["n_states"]) == n_states || error("Inconsistent n_states in $file")
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

    out_dir = joinpath("exm", "data", "Bulk_measure", "lyapunov_spectrum_ytau", "L$(L)")
    mkpath(out_dir)
    output_path = joinpath(
        out_dir,
        backend == :mps ?
        "ensemble_lyapunov_ytau_L$(L)_gammaind$(τ_idx)_t$(t)_chi$(χ).jld2" :
        "ensemble_lyapunov_ytau_L$(L)_gammaind$(τ_idx)_t$(t).jld2",
    )
    jldsave(
        output_path;
        backend = String(backend),
        topological_sector = "y=tau",
        y_eigenvalue = -2 / (1 + √5),
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
    println("Usage: julia --project=. exm/Bulk_measure/lyapunov_spectrum_ytau.jl backend [args...]")
    println("")
    println("Backends (measurement strength by index τ_idx into config.jl's τlis;")
    println("evolution time t in periods from get_cfg_params_Born/get_mps_params_Born):")
    println("")
    println("  exact L τ_idx n_states seed [output.jld2]")
    println("      Born trajectory (exact state vector) from Pτ|τ⋯τ⟩, plus the y=τ-sector")
    println("      Lyapunov spectrum of the transfer-matrix product along it.")
    println("")
    println("  mps L τ_idx chi n_states seed [output.jld2]")
    println("      MPS version for large L; the y=τ sector is enforced and monitored")
    println("      through the topological-charge MPO.")
    println("")
    println("  collect exact L τ_idx")
    println("  collect mps L τ_idx chi")
    println("      Aggregate the per-seed spectra into an ensemble file in the L directory.")
    println("")
    println("Examples:")
    println("  julia --project=. exm/Bulk_measure/lyapunov_spectrum_ytau.jl exact 8 10 3 1")
    println("  julia --project=. exm/Bulk_measure/lyapunov_spectrum_ytau.jl collect exact 8 10")
else
    backend = Symbol(lowercase(ARGS[1]))

    if backend == :collect
        # -------------------------------------------------------------------
        # collect: ensemble-average the per-seed spectra
        # -------------------------------------------------------------------
        length(ARGS) >= 2 || error("collect requires a sub-backend: exact or mps")
        sub = Symbol(lowercase(ARGS[2]))
        if sub == :exact
            length(ARGS) == 4 || error("Usage: collect exact L τ_idx")
            output_path = collect_lyapunov_ytau(
                :exact,
                parse(Int, ARGS[3]),
                parse(Int, ARGS[4]),
            )
        elseif sub == :mps
            length(ARGS) == 5 || error("Usage: collect mps L τ_idx chi")
            output_path = collect_lyapunov_ytau(
                :mps,
                parse(Int, ARGS[3]),
                parse(Int, ARGS[4]);
                χ = parse(Int, ARGS[5]),
            )
        else
            error("Unknown collect backend: $sub")
        end
        println("saved: $output_path")

    elseif backend == :exact || backend == :mps
        # -------------------------------------------------------------------
        # exact/mps: single Born trajectory + y=τ-sector Lyapunov spectrum
        # -------------------------------------------------------------------
        if backend == :exact
            length(ARGS) in (5, 6) ||
                error("Usage: exact L τ_idx n_states seed [output.jld2]")
            L = parse(Int, ARGS[2])
            τ_idx = parse(Int, ARGS[3])
            n_states = parse(Int, ARGS[4])
            seed = parse(Int, ARGS[5])
            χ = 0
        else
            length(ARGS) in (6, 7) ||
                error("Usage: mps L τ_idx chi n_states seed [output.jld2]")
            L = parse(Int, ARGS[2])
            τ_idx = parse(Int, ARGS[3])
            χ = parse(Int, ARGS[4])
            n_states = parse(Int, ARGS[5])
            seed = parse(Int, ARGS[6])
        end
        τ = τlis[τ_idx]
        t = lyapunov_ytau_time(backend, L, τ_idx)

        println("=== y=τ Lyapunov spectrum ($backend backend) ===")
        println("L = $L, τ_idx = $τ_idx, τ = $τ, γ = $(tanh(τ))" *
                (backend == :mps ? ", χ = $χ" : ""))
        println("t = $t, n_states = $n_states, seed = $seed")

        trajectory, spectrum = if backend == :exact
            simulate_ytau_lyapunov(L, τ, t; n_states = n_states, trajectory_seed = seed)
        else
            simulate_ytau_lyapunov_mps(L, τ, t, χ; n_states = n_states, trajectory_seed = seed)
        end

        has_custom_output = (backend == :exact && length(ARGS) == 6) ||
                            (backend == :mps && length(ARGS) == 7)
        output_path = has_custom_output ? ARGS[end] : joinpath(
            lyapunov_ytau_dir(backend, L, τ_idx; χ = backend == :mps ? χ : nothing),
            "lyapunov_ytau_L$(L)_t$(t)_seed$(seed).jld2",
        )
        mkpath(dirname(output_path))
        jldsave(
            output_path;
            backend = String(backend),
            L = L,
            τ_idx = τ_idx,
            τ = τ,
            gamma = tanh(τ),
            χ = χ,
            t = t,
            n_states = n_states,
            trajectory_seed = seed,
            topological_sector = "y=tau",
            y_eigenvalue = -2 / (1 + √5),
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
            println("y=τ sector dimension: $(spectrum.sector_dimension)")
        !isempty(spectrum.sector_leakage) &&
            println("max sector leakage: $(maximum(spectrum.sector_leakage))")
        println("final Lyapunov exponents: $(spectrum.lyapunov_exponents[:, end])")
        println("final free-energy spectrum: $(spectrum.free_energy_spectrum[:, end])")
    else
        error("Unknown backend: $backend")
    end
end
