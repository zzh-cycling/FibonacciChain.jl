using FibonacciChain
using JLD2
using LinearAlgebra
using Random
using Statistics

include(joinpath(@__DIR__, "config.jl"))

const MOMENTUM_DATA_ROOT =
    joinpath("exm", "data", "Bulk_measure", "momentum_spectrum")

function print_spectrum(label, spectrum)
    println("\n=== $label ===")
    println("rank    -log|μ|       k-index       k/π            Y")
    for index in eachindex(spectrum.free_energies)
        println(
            lpad(index, 4),
            "  ",
            lpad(round(spectrum.free_energies[index]; digits = 8), 12),
            "  ",
            lpad(spectrum.momentum_indices[index], 7),
            "  ",
            lpad(round(spectrum.momenta[index] / π; digits = 6), 10),
            "  ",
            round(spectrum.topological_charges[index]; digits = 8),
        )
    end
    println("[G,U] residual: $(spectrum.commutator_residual)")
    println("G²-T²U residual: $(spectrum.glide_relation_residual)")
end

"""
Benchmark the physical one-site momentum convention on uniform exact
post-selection trajectories. The TCI and Potts spectra are resolved with the
brickwork glide symmetry rather than `T²`, which would fold `k=0` and `k=π`.
For a Potts benchmark choose `L` divisible by three so that the low pair at
`k=±2π/3` is represented exactly (for example `L=12`).
"""
function benchmark_postselection_momentum(
    L::Int,
    τ_idx::Int,
    n_states::Int;
    output_path::Union{Nothing,String} = nothing,
)
    τ = τlis[τ_idx]
    model = fib_model(L)
    tci = postselection_glide_spectrum(model, τ, true; n_states = n_states)
    potts = postselection_glide_spectrum(model, τ, false; n_states = n_states)

    maximum(tci.eigenvector_residuals) < 1e-10 ||
        error("TCI transfer/glide eigenvectors failed the residual check")
    maximum(potts.eigenvector_residuals) < 1e-10 ||
        error("Potts transfer/glide eigenvectors failed the residual check")
    maximum(tci.momentum_quantization_residuals) < 1e-10 ||
        error("TCI glide phases are not quantized physical momenta")
    maximum(potts.momentum_quantization_residuals) < 1e-10 ||
        error("Potts glide phases are not quantized physical momenta")
    if n_states >= 4
        expected_tci = [0, div(L, 2), 0, div(L, 2)]
        tci.momentum_indices[1:4] == expected_tci || error(
            "TCI benchmark failed: expected first momenta $expected_tci, got " *
            "$(tci.momentum_indices[1:4])",
        )
    end
    if n_states >= 3 && L % 3 == 0
        expected_potts_pair = Set((div(L, 3), 2 * div(L, 3)))
        Set(potts.momentum_indices[2:3]) == expected_potts_pair || error(
            "Potts benchmark failed: expected the first pair at ±2π/3, got " *
            "$(potts.momentum_indices[2:3])",
        )
    end

    print_spectrum("TCI / antiferromagnetic post-selection", tci)
    print_spectrum("three-state Potts / ferromagnetic post-selection", potts)

    path = isnothing(output_path) ? joinpath(
        MOMENTUM_DATA_ROOT,
        "postselection_glide_L$(L)_gammaind$(τ_idx).jld2",
    ) : output_path
    mkpath(dirname(path))
    jldsave(
        path;
        L = L,
        τ_idx = τ_idx,
        τ = τ,
        gamma = tanh(τ),
        n_states = n_states,
        tci = tci,
        potts = potts,
    )
    return path, tci, potts
end

function sector_data(sector::Symbol)
    ϕ = (1 + √5) / 2
    if sector == :y1
        return (label = "y=1", api = :trivial, eigenvalue = ϕ, other = -inv(ϕ))
    elseif sector == :ytau
        return (label = "y=tau", api = :tau, eigenvalue = -inv(ϕ), other = ϕ)
    end
    error("sector must be y1 or ytau")
end

function initial_sector_state(model, spec)
    Y = topological_charge_operator(model)
    projector = (Y - spec.other * I(size(Y, 1))) /
                (spec.eigenvalue - spec.other)
    state = zeros(Float64, length(anyon_basis(model)))
    state[1] = 1.0
    state = projector * state
    normalize!(state)
    return state
end

"""
Generate one exact Born trajectory in a fixed topological sector and record
the physical one-site momentum distribution of every sector-resolved
Lyapunov/QR vector after every period. `average_start` is the first period
included in the saved time average. For a near-degenerate group `A`, average
`momentum_weights[:, A, :]` over the state index as well as trajectories/time;
this equals `Tr(P_A P_k)/dim(A)` and is invariant under QR rotations in `A`.
"""
function born_momentum(
    sector::Symbol,
    L::Int,
    τ_idx::Int,
    periods::Int,
    n_states::Int,
    seed::Int;
    average_start::Int = max(1, periods ÷ 2),
    output_path::Union{Nothing,String} = nothing,
)
    1 <= average_start <= periods ||
        throw(ArgumentError("average_start must lie in 1:periods"))
    spec = sector_data(sector)
    τ = τlis[τ_idx]
    model = fib_model(L)
    initial_state = initial_sector_state(model, spec)
    config = MeasureConfig(
        τ = τ,
        t₂ = periods,
        mode = :Born,
        rng = MersenneTwister(seed),
        enable_τ_eff = false,
        track_y_expectation = true,
    )
    trajectory = bulk_evolution(model, initial_state, config)
    spectrum = lyapunov_spectrum_topological_sector(
        model,
        τ,
        trajectory.samples;
        sector = spec.api,
        n_states = n_states,
        track_momentum = true,
    )

    maximum(abs.(trajectory.y_expectation_values .- spec.eigenvalue)) < 1e-5 ||
        error("Born trajectory left the $(spec.label) sector")
    maximum(spectrum.sector_leakage) < 1e-9 ||
        error("Lyapunov frame left the $(spec.label) sector")

    time_average = dropdims(
        mean(spectrum.momentum_weights[:, :, average_start:end]; dims = 3);
        dims = 3,
    )
    path = isnothing(output_path) ? joinpath(
        MOMENTUM_DATA_ROOT,
        "$(sector)",
        "L$(L)",
        "gammaind$(τ_idx)",
        "momentum_L$(L)_t$(periods)_seed$(seed).jld2",
    ) : output_path
    mkpath(dirname(path))
    jldsave(
        path;
        backend = "exact",
        L = L,
        τ_idx = τ_idx,
        τ = τ,
        gamma = tanh(τ),
        periods = periods,
        n_states = n_states,
        trajectory_seed = seed,
        topological_sector = spec.label,
        y_eigenvalue = spec.eigenvalue,
        average_start = average_start,
        samples = trajectory.samples,
        y_expectation_values = trajectory.y_expectation_values,
        local_log_stretches = spectrum.local_log_stretches,
        lyapunov_exponents = spectrum.lyapunov_exponents,
        free_energy_spectrum = spectrum.free_energy_spectrum,
        sector_leakage = spectrum.sector_leakage,
        momentum_weights = spectrum.momentum_weights,
        time_averaged_momentum_weights = time_average,
    )

    println("\n=== Exact Born momentum: $(spec.label) ===")
    println("L=$L, τ_idx=$τ_idx, γ=$(tanh(τ)), periods=$periods, seed=$seed")
    for state_index = 1:n_states
        dominant_k = argmax(@view time_average[:, state_index]) - 1
        dominant_weight = time_average[dominant_k + 1, state_index]
        println(
            "state $state_index: k=$dominant_k (k/π=$(round(2 * dominant_k / L; digits=6))), " *
            "weight=$(round(dominant_weight; digits=6))",
        )
    end
    return path, trajectory, spectrum
end

"""
Collect every per-seed exact Born momentum file for fixed
`(sector, L, τ_idx, periods)`. The ensemble output contains the mean and
standard error of the already time-averaged weights, plus the final finite-time
Lyapunov exponents for every seed.
"""
function collect_born_momentum(
    sector::Symbol,
    L::Int,
    τ_idx::Int,
    periods::Int;
    output_path::Union{Nothing,String} = nothing,
)
    spec = sector_data(sector)
    data_dir = joinpath(
        MOMENTUM_DATA_ROOT,
        "$(sector)",
        "L$(L)",
        "gammaind$(τ_idx)",
    )
    isdir(data_dir) || error("data directory does not exist: $data_dir")
    prefix = "momentum_L$(L)_t$(periods)_seed"
    files = sort(filter(
        file -> startswith(file, prefix) && endswith(file, ".jld2"),
        readdir(data_dir),
    ))
    isempty(files) && error("no per-seed files matching $prefix in $data_dir")

    weights = Matrix{Float64}[]
    final_exponents = Vector{Float64}[]
    seeds = Int[]
    average_starts = Int[]
    n_states = 0
    for file in files
        data = JLD2.load(joinpath(data_dir, file))
        String(data["backend"]) == "exact" || error("non-exact backend in $file")
        Int(data["L"]) == L || error("inconsistent L in $file")
        Int(data["τ_idx"]) == τ_idx || error("inconsistent τ_idx in $file")
        Int(data["periods"]) == periods || error("inconsistent periods in $file")
        String(data["topological_sector"]) == spec.label ||
            error("inconsistent sector in $file")
        this_n_states = Int(data["n_states"])
        n_states == 0 && (n_states = this_n_states)
        this_n_states == n_states || error("inconsistent n_states in $file")
        push!(weights, Float64.(data["time_averaged_momentum_weights"]))
        push!(final_exponents, Float64.(data["lyapunov_exponents"][:, end]))
        push!(seeds, Int(data["trajectory_seed"]))
        push!(average_starts, Int(data["average_start"]))
    end
    length(unique(seeds)) == length(seeds) || error("duplicate trajectory seeds")
    length(unique(average_starts)) == 1 || error(
        "cannot collect files with different time-averaging windows: $average_starts",
    )

    weight_stack = cat(weights...; dims = 3)
    mean_weights = dropdims(mean(weight_stack; dims = 3); dims = 3)
    stderr_weights = dropdims(
        std(weight_stack; dims = 3, corrected = false);
        dims = 3,
    ) ./ √length(files)
    exponent_matrix = hcat(final_exponents...)

    path = isnothing(output_path) ? joinpath(
        data_dir,
        "ensemble_momentum_L$(L)_t$(periods).jld2",
    ) : output_path
    mkpath(dirname(path))
    jldsave(
        path;
        backend = "exact",
        L = L,
        τ_idx = τ_idx,
        τ = τlis[τ_idx],
        gamma = tanh(τlis[τ_idx]),
        periods = periods,
        n_states = n_states,
        topological_sector = spec.label,
        y_eigenvalue = spec.eigenvalue,
        average_start = only(unique(average_starts)),
        samples_num = length(files),
        ensemble_seeds = seeds,
        mean_momentum_weights = mean_weights,
        stderr_momentum_weights = stderr_weights,
        final_lyapunov_exponents_per_seed = exponent_matrix,
        mean_final_lyapunov_exponents = vec(mean(exponent_matrix; dims = 2)),
        stderr_final_lyapunov_exponents =
            vec(std(exponent_matrix; dims = 2, corrected = false)) ./ √length(files),
    )
    println("collected $(length(files)) exact Born trajectories from $data_dir")
    return path, mean_weights, stderr_weights
end

function usage()
    println("Exact momentum workflows:")
    println("  julia --project=. exm/Bulk_measure/momentum_spectrum.jl benchmark L τ_idx n_states [output.jld2]")
    println("  julia --project=. exm/Bulk_measure/momentum_spectrum.jl born SECTOR L τ_idx periods n_states seed [average_start] [output.jld2]")
    println("  julia --project=. exm/Bulk_measure/momentum_spectrum.jl collect SECTOR L τ_idx periods [output.jld2]")
    println("SECTOR is y1 or ytau. Use L divisible by 3 for the Potts benchmark.")
end

function main(args)
    isempty(args) && return usage()
    mode = Symbol(lowercase(args[1]))
    if mode == :benchmark
        length(args) in (4, 5) || error("benchmark requires L τ_idx n_states [output]")
        path, _, _ = benchmark_postselection_momentum(
            parse(Int, args[2]),
            parse(Int, args[3]),
            parse(Int, args[4]);
            output_path = length(args) == 5 ? args[5] : nothing,
        )
        println("saved: $path")
    elseif mode == :born
        length(args) in (7, 8, 9) || error(
            "born requires SECTOR L τ_idx periods n_states seed [average_start] [output]",
        )
        periods = parse(Int, args[5])
        path, _, _ = born_momentum(
            Symbol(lowercase(args[2])),
            parse(Int, args[3]),
            parse(Int, args[4]),
            periods,
            parse(Int, args[6]),
            parse(Int, args[7]);
            average_start = length(args) >= 8 ? parse(Int, args[8]) : max(1, periods ÷ 2),
            output_path = length(args) == 9 ? args[9] : nothing,
        )
        println("saved: $path")
    elseif mode == :collect
        length(args) in (5, 6) || error(
            "collect requires SECTOR L τ_idx periods [output]",
        )
        path, _, _ = collect_born_momentum(
            Symbol(lowercase(args[2])),
            parse(Int, args[3]),
            parse(Int, args[4]),
            parse(Int, args[5]);
            output_path = length(args) == 6 ? args[6] : nothing,
        )
        println("saved: $path")
    else
        error("unknown mode: $mode")
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
