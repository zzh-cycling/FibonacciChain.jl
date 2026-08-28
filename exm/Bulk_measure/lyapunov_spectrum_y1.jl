using FibonacciChain
using ITensorMPS
using ITensors
using JLD2
using LinearAlgebra
using Random
using Statistics

"""
    initial_y1_state(model) -> Vector

Project the all-`τ` fusion path into the `y=1` sector and normalize it. This is
the state used to sample the Born trajectory. In the repository normalization
the `y=1` sector has topological charge eigenvalue `ϕ`, while the `y=τ` sector
has eigenvalue `-1/ϕ`, so `P1 = (Y + ϕ⁻¹ I) / (ϕ + ϕ⁻¹)`.
"""
function initial_y1_state(model::AnyonModel{FibonacciAnyon})
    model.pbc || error("The topological y=1 sector requires periodic boundaries")
    ϕ = (1 + √5) / 2
    Y = topological_charge_operator(model)
    P1 = (Y + inv(ϕ) * I(size(Y, 1))) / (ϕ + inv(ϕ))
    state = zeros(Float64, length(anyon_basis(model)))
    state[1] = 1.0
    projected = P1 * state
    weight = real(dot(projected, projected))
    weight > 1e-14 || error("The all-τ path has zero y=1 weight")
    return projected / sqrt(weight)
end

"""
    initial_y1_mps(model; cutoff=1e-14, maxdim=32)

MPS version of [`initial_y1_state`](@ref): apply the topological charge MPO to
the all-τ product state and add the identity part,

    P1|00⋯0⟩ ∝ (Y + ϕ⁻¹ I)|00⋯0⟩,

so no dense Y matrix is ever constructed and large L is accessible.
"""
function initial_y1_mps(model::AnyonModel{FibonacciAnyon}; cutoff::Float64 = 1e-14, maxdim::Int = 32)
    model.pbc || error("The topological y=1 sector requires periodic boundaries")
    N = model.N
    N >= 2 || error("initial_y1_mps requires at least two sites")
    sites = siteinds("Qubit", N)

    vacuum = productMPS(sites, fill("0", N))
    Y_mpo = topological_charge_mpo(sites; pbc = true)
    y_on_vacuum = apply(Y_mpo, vacuum; cutoff = cutoff, maxdim = maxdim)

    ϕ = (1 + √5) / 2
    vacuum[1] = inv(ϕ) * vacuum[1]
    projected = add(y_on_vacuum, vacuum; cutoff = cutoff, maxdim = maxdim)
    normalize!(projected)
    return projected, sites
end

"""
    simulate_y1_lyapunov(L, τ, periods; n_states=10, trajectory_seed=1)

Generate a Born trajectory from `P1|τ⋯τ⟩` (exact state-vector evolution), then
compute the Lyapunov spectrum of the corresponding transfer-matrix product
restricted to the full `y=1` sector via `lyapunov_spectrum_topological_sector`.
"""
function simulate_y1_lyapunov(
    L::Int,
    τ::Real,
    periods::Int;
    n_states::Int = 10,
    trajectory_seed::Int = 1,
)
    periods >= 1 || throw(ArgumentError("periods must be positive"))
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    initial_state = initial_y1_state(model)
    born_config = MeasureConfig(
        τ = Float64(τ),
        t₂ = periods,
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
        sector = :trivial,
        n_states = n_states,
    )

    ϕ = (1 + √5) / 2
    maximum(abs.(trajectory.y_expectation_values .- ϕ)) < 1e-5 ||
        error("the Born trajectory left the y=1 sector")
    maximum(spectrum.sector_leakage) < 1e-9 ||
        error("the Lyapunov frame left the y=1 sector")
    return trajectory, (
        local_log_stretches = spectrum.local_log_stretches,
        lyapunov_exponents = spectrum.lyapunov_exponents,
        free_energy_spectrum = spectrum.free_energy_spectrum,
        sector_dimension = spectrum.sector_dimension,
        sector_leakage = spectrum.sector_leakage,
    )
end

"""
    simulate_y1_lyapunov_mps(L, τ, periods, χ; n_states=10, trajectory_seed=1)

MPS version of [`simulate_y1_lyapunov`](@ref) for large L. The Born trajectory
is generated from the Y-MPO-projected initial state (`initial_y1_mps`), with
the `y=1` sector monitored on-the-fly by `topological_charge_mpo`. The
Lyapunov spectrum is computed by `lyapunov_spectrum_mps` with
`sector = :trivial`, which projects the propagated frame back into the `y=1`
sector after every period via the Y MPO; since the transfer matrix commutes
with `Y`, this only removes MPS truncation leakage into the `y=τ` sector.
"""
function simulate_y1_lyapunov_mps(
    L::Int,
    τ::Real,
    periods::Int,
    χ::Int;
    n_states::Int = 10,
    trajectory_seed::Int = 1,
)
    periods >= 1 || throw(ArgumentError("periods must be positive"))
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    initial_state, sites = initial_y1_mps(model; cutoff = 1e-12, maxdim = χ)
    born_config = MeasureConfig(
        τ = Float64(τ),
        t₂ = periods,
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
        sector = :trivial,
        cutoff = 1e-12,
        maxdim = χ,
    )
    local_log_stretches = -spectrum_mps
    periods_out = size(local_log_stretches, 2)
    lyapunov_exponents =
        cumsum(local_log_stretches; dims = 2) ./ reshape(collect(1:periods_out), 1, :)

    ϕ = (1 + √5) / 2
    maximum(abs.(trajectory.y_expectation_values .- ϕ)) < 1e-4 ||
        error("the MPS Born trajectory left the y=1 sector")
    return trajectory, (
        local_log_stretches = local_log_stretches,
        lyapunov_exponents = lyapunov_exponents,
        free_energy_spectrum = -lyapunov_exponents,
        sector_dimension = 0,  # not available without the dense Y
        sector_leakage = Float64[],  # tracked instead via y_expectation_values
    )
end


"""
    lyapunov_y1_dir(backend) -> String

Directory holding the per-seed data files of this script for `backend`
(`:exact` or `:mps`); MPS files live in the `mps` subdirectory.
"""
function lyapunov_y1_dir(backend::Symbol)
    dir = joinpath(@__DIR__, "..", "data", "Bulk_measure", "lyapunov_spectrum_y1")
    return backend == :mps ? joinpath(dir, "mps") : dir
end

"""
    collect_lyapunov_y1(backend, L, τ, periods; χ=nothing)

Aggregate the per-seed `y=1` Lyapunov spectra produced by the `exact`/`mps`
modes of this script into a single ensemble file: element-wise mean and
standard error (over seeds) of `local_log_stretches`, `lyapunov_exponents`,
`free_energy_spectrum` and `y_expectation_values`, plus the per-seed final
exponents. Returns the output path.
"""
function collect_lyapunov_y1(
    backend::Symbol,
    L::Integer,
    τ::Real,
    periods::Integer;
    χ::Union{Nothing,Integer} = nothing,
)
    backend in (:exact, :mps) || throw(ArgumentError("backend must be :exact or :mps"))
    backend == :mps && χ === nothing &&
        throw(ArgumentError("collecting the mps backend requires χ"))
    data_dir = lyapunov_y1_dir(backend)
    prefix = backend == :mps ? "L$(L)_tau$(τ)_chi$(χ)_" : "L$(L)_tau$(τ)_"
    files = sort(filter(
        file -> startswith(file, prefix) &&
                contains(file, "_periods$(periods)_seed") &&
                endswith(file, ".jld2"),
        readdir(data_dir),
    ))
    isempty(files) &&
        error("No $backend files with L=$L, τ=$τ, periods=$periods found in $data_dir")

    first_data = JLD2.load(joinpath(data_dir, first(files)))
    Int(first_data["periods"]) == periods || error("Inconsistent evolution length")
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
        data["topological_sector"] == "y=1" || error("Unexpected sector in $file")
        Int(data["L"]) == L && Int(data["periods"]) == periods ||
            error("Inconsistent system size or evolution length in $file")
        Int(data["n_states"]) == n_states || error("Inconsistent n_states in $file")
        seeds[i] = Int(data["trajectory_seed"])
        log_stretches[i] = Float64.(data["local_log_stretches"])
        exponents[i] = Float64.(data["lyapunov_exponents"])
        free_energies[i] = Float64.(data["free_energy_spectrum"])
        y_expectations[i] = Float64.(data["y_expectation_values"])
        sector_leakages[i] = Float64.(get(data, "sector_leakage", Float64[]))
    end
    length(unique(seeds)) == samples_num || error("Duplicate trajectory seeds found")
    all(size(m) == (n_states, periods) for m in exponents) ||
        error("Inconsistent spectrum shapes")
    all(length(y) == periods for y in y_expectations) ||
        error("Inconsistent y_expectation lengths")

    stderr_of(mats) = std(mats; corrected = false) ./ √samples_num
    final_exponents = hcat([m[:, end] for m in exponents]...)
    y_matrix = reduce(hcat, y_expectations)  # periods × samples_num

    output_path = joinpath(
        data_dir,
        backend == :mps ?
        "collect_L$(L)_tau$(τ)_chi$(χ)_periods$(periods)_nstates$(n_states).jld2" :
        "collect_L$(L)_tau$(τ)_periods$(periods)_nstates$(n_states).jld2",
    )
    jldsave(
        output_path;
        backend = String(backend),
        topological_sector = "y=1",
        y_eigenvalue = (1 + √5) / 2,
        L = Int(L),
        τ = Float64(τ),
        χ = isnothing(χ) ? 0 : Int(χ),
        periods = Int(periods),
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
    return output_path
end


function main(args = ARGS)
    length(args) >= 1 || error("Usage: see below")
    backend = Symbol(lowercase(args[1]))
    usage =
        "Usage:\n" *
        "  julia --project=. exm/Bulk_measure/lyapunov_spectrum_y1.jl exact L tau periods n_states seed [output.jld2]\n" *
        "  julia --project=. exm/Bulk_measure/lyapunov_spectrum_y1.jl mps L tau chi periods n_states seed [output.jld2]\n" *
        "  julia --project=. exm/Bulk_measure/lyapunov_spectrum_y1.jl collect exact L tau periods\n" *
        "  julia --project=. exm/Bulk_measure/lyapunov_spectrum_y1.jl collect mps L tau chi periods"

    if backend == :collect
        length(args) >= 2 || error(usage)
        sub = Symbol(lowercase(args[2]))
        if sub == :exact
            length(args) == 5 || error(usage)
            output_path = collect_lyapunov_y1(
                :exact,
                parse(Int, args[3]),
                parse(Float64, args[4]),
                parse(Int, args[5]),
            )
        elseif sub == :mps
            length(args) == 6 || error(usage)
            output_path = collect_lyapunov_y1(
                :mps,
                parse(Int, args[3]),
                parse(Float64, args[4]),
                parse(Int, args[6]);
                χ = parse(Int, args[5]),
            )
        else
            error(usage)
        end
        println("saved: $output_path")
        return output_path
    end

    if backend == :exact
        length(args) in (6, 7) || error(usage)
        L = parse(Int, args[2])
        τ = parse(Float64, args[3])
        periods = parse(Int, args[4])
        n_states = parse(Int, args[5])
        seed = parse(Int, args[6])
        trajectory, spectrum = simulate_y1_lyapunov(
            L,
            τ,
            periods;
            n_states = n_states,
            trajectory_seed = seed,
        )
        χ = 0
        output_path = length(args) == 7 ? args[7] : joinpath(
            lyapunov_y1_dir(:exact),
            "L$(L)_tau$(τ)_periods$(periods)_seed$(seed).jld2",
        )
    elseif backend == :mps
        length(args) in (7, 8) || error(usage)
        L = parse(Int, args[2])
        τ = parse(Float64, args[3])
        χ = parse(Int, args[4])
        periods = parse(Int, args[5])
        n_states = parse(Int, args[6])
        seed = parse(Int, args[7])
        trajectory, spectrum = simulate_y1_lyapunov_mps(
            L,
            τ,
            periods,
            χ;
            n_states = n_states,
            trajectory_seed = seed,
        )
        output_path = length(args) == 8 ? args[8] : joinpath(
            lyapunov_y1_dir(:mps),
            "L$(L)_tau$(τ)_chi$(χ)_periods$(periods)_seed$(seed).jld2",
        )
    else
        error(usage)
    end

    mkpath(dirname(output_path))
    jldsave(
        output_path;
        backend = String(backend),
        L = L,
        τ = τ,
        χ = χ,
        periods = periods,
        n_states = n_states,
        trajectory_seed = seed,
        topological_sector = "y=1",
        y_eigenvalue = (1 + √5) / 2,
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
        println("y=1 sector dimension: $(spectrum.sector_dimension)")
    !isempty(spectrum.sector_leakage) &&
        println("max sector leakage: $(maximum(spectrum.sector_leakage))")
    println("final Lyapunov exponents: $(spectrum.lyapunov_exponents[:, end])")
    println("final free-energy spectrum: $(spectrum.free_energy_spectrum[:, end])")
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
