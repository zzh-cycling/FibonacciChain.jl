using Distributed
using FibonacciChain
using ITensorMPS
using ITensors
using JLD2
using LinearAlgebra
using Random
using Statistics

# Monitored dynamics from the sector-projected all-τ state (y=1 or y=τ).
# Merged from monitored_dynamics_y1.jl and monitored_dynamics_ytau.jl;
# data roots, file names and JLD2 schemas are unchanged.
#
# Launch trajectory generation with `julia -p N` to use N worker processes.
const BULK_MEASURE_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BULK_MEASURE_CONFIG)

@everywhere begin
    using FibonacciChain
    using ITensorMPS
    using ITensors
    using JLD2
    using LinearAlgebra
    using Random
    using Statistics

    # In the repository normalization the y=1 sector has topological charge
    # eigenvalue ϕ, while the y=τ sector has eigenvalue -1/ϕ.  With `y` the
    # sector eigenvalue and `ȳ` the other one, the projector is
    # P = (Y - ȳ I)/(y - ȳ); the y=1 weight of the all-τ path is
    # w₁(N) = (1 + ϕ (-ϕ⁻¹)^N)/(1 + ϕ²), and the y=τ weight is 1 - w₁(N).
    const SECTOR_DYN_SPECS = Dict{Symbol,NamedTuple}(
        :y1 => (;
            title = "y=1",
            data_root = joinpath("exm", "data", "Bulk_measure", "monitored_dynamics_y1"),
            file_suffix = "y1",
            sector_label = "y=1",
            y_eigenvalue = (1 + √5) / 2,
            other_eigenvalue = -2 / (1 + √5),
            vacuum_sector_weight =
                N -> (1 + ((1 + √5) / 2) * (-2 / (1 + √5))^N) / (1 + ((1 + √5) / 2)^2),
            # Historical asymmetry kept: only the y=τ MPS trajectory enforces
            # the Fibonacci constraint during evolution.
            enforce_constraint = false,
        ),
        :ytau => (;
            title = "y=τ",
            data_root = joinpath("exm", "data", "Bulk_measure", "monitored_dynamics_ytau"),
            file_suffix = "ytau",
            sector_label = "y=tau",
            y_eigenvalue = -2 / (1 + √5),
            other_eigenvalue = (1 + √5) / 2,
            vacuum_sector_weight =
                N -> 1 - (1 + ((1 + √5) / 2) * (-2 / (1 + √5))^N) / (1 + ((1 + √5) / 2)^2),
            enforce_constraint = true,
        ),
    )

    function sector_dyn_spec(sector::Symbol)
        haskey(SECTOR_DYN_SPECS, sector) ||
            error("Unknown sector: $sector (expected :y1 or :ytau)")
        return SECTOR_DYN_SPECS[sector]
    end

    """Return the normalized projection of the all-τ path into the sector."""
    function initial_sector_exact(model, spec)
        model.pbc ||
            error("The topological $(spec.title) sector requires periodic boundaries")
        basis = anyon_basis(model)
        state = zeros(Float64, length(basis))
        state[1] = 1.0

        # Y is symmetric, so the row returned by topological_symmetry_basismap
        # is also the single column Y|00⋯0⟩. Avoid constructing the full dense
        # Y matrix, whose size grows exponentially with L.
        y_on_state = topological_symmetry_basismap(model, first(basis))
        projected =
            (y_on_state .- spec.other_eigenvalue .* state) ./
            (spec.y_eigenvalue - spec.other_eigenvalue)
        projected_norm = norm(projected)
        projected_norm > 1e-14 ||
            error("The initial state has zero $(spec.title) weight")

        expected_weight = spec.vacuum_sector_weight(model.N)
        isapprox(projected_norm^2, expected_weight; atol = 1e-10, rtol = 1e-10) ||
            error(
                "Unexpected $(spec.title) weight $(projected_norm^2); " *
                "expected $expected_weight",
            )
        projected ./= projected_norm
        return projected
    end

    """
        initial_sector_mps(model, spec; cutoff=1e-14, maxdim=32)

    Construct the same projected all-τ state as [`initial_sector_exact`](@ref),
    directly as an MPS, by applying the topological charge MPO to the all-τ
    product state:

        P|00⋯0⟩ ∝ (Y - ȳ I)|00⋯0⟩.

    This avoids constructing the exponentially large dense topological-charge
    operator and remains suitable for the large-L MPS simulations.
    """
    function initial_sector_mps(model, spec; cutoff::Float64 = 1e-14, maxdim::Int = 32)
        model.pbc ||
            error("The topological $(spec.title) sector requires periodic boundaries")
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

    function sector_dyn_trajectory_dir(
        spec,
        backend::Symbol,
        L::Integer,
        τ_idx::Integer;
        χ::Union{Nothing,Integer} = nothing,
    )
        backend ∈ (:exact, :mps) || error("backend must be :exact or :mps")
        backend == :mps && isnothing(χ) && error("χ is required for the MPS backend")
        path = joinpath(
            spec.data_root,
            String(backend),
            "L$(L)",
            "gammaind$(τ_idx)",
        )
        return backend == :mps ? joinpath(path, "chi$(χ)") : path
    end

    function sector_dyn_default_periods(backend::Symbol, L::Integer, τ_idx::Integer)
        if backend == :exact
            D, _, _ = get_cfg_params_Born(τ_idx, L)
            return div(D, 2)
        elseif backend == :mps
            time_in_L, _, _ = get_mps_params_Born(τ_idx, L)
            return time_in_L * L
        end
        error("backend must be :exact or :mps")
    end

    function simulate_sector_exact(
        spec,
        L::Integer,
        τ_idx::Integer,
        seed::Integer;
        periods::Integer = sector_dyn_default_periods(:exact, L, τ_idx),
    )
        iseven(L) || error("L must be even, got $L")
        τ = τlis_ext[τ_idx]
        model = fib_model(L)
        initial_state = initial_sector_exact(model, spec)
        config = MeasureConfig(
            τ = τ,
            mode = :Born,
            t₂ = periods,
            rng = MersenneTwister(seed),
        )
        outcome = bulk_evolution(model, initial_state, config)
        return (
            state = outcome.state,
            samples = outcome.samples,
            free_energy = outcome.free_energys,
            halfchain_entropy = outcome.entanglement_entropys,
            final_entropy_profile = anyon_eelis(model, outcome.state),
            initial_bond_dimension = 0
        )
    end

    function simulate_sector_mps(
        spec,
        L::Integer,
        τ_idx::Integer,
        χ::Integer,
        seed::Integer;
        periods::Integer = sector_dyn_default_periods(:mps, L, τ_idx),
    )
        iseven(L) || error("L must be even, got $L")
        χ >= 5 ||
            error("χ must be at least 5 to represent the exact $(spec.title) initial MPS")
        τ = τlis_ext[τ_idx]
        model = fib_model(L)
        initial_state, sites = initial_sector_mps(model, spec; maxdim = χ)
        initial_bond_dimension = maxlinkdim(initial_state)
        config = MeasureConfig(
            τ = τ,
            mode = :Born,
            t₂ = periods,
            rng = MersenneTwister(seed),
            cutoff = 1e-12,
            maxdim = χ,
            enforce_fibonacci_constraint = spec.enforce_constraint,
            track_y_expectation = true,
        )
        outcome = bulk_evolution(model, sites, initial_state, config)
        maximum(abs.(outcome.y_expectation_values .- spec.y_eigenvalue)) < 1e-4 ||
            error("the MPS trajectory left the $(spec.title) sector")
        return (
            state = outcome.state,
            samples = outcome.samples,
            free_energy = outcome.free_energys,
            halfchain_entropy = outcome.entanglement_entropys,
            y_expectation = outcome.y_expectation_values,
            final_entropy_profile = anyon_eelis(model, outcome.state),
            initial_bond_dimension = initial_bond_dimension,
        )
    end

    function save_sector_trajectory(
        spec,
        backend::Symbol,
        L::Integer,
        τ_idx::Integer,
        seed::Integer,
        result;
        χ::Union{Nothing,Integer} = nothing,
    )
        out_dir = sector_dyn_trajectory_dir(spec, backend, L, τ_idx; χ = χ)
        mkpath(out_dir)
        periods = length(result.halfchain_entropy)
        output_path =
            joinpath(out_dir, "periods$(periods)_trajectory_seed$(seed).jld2")
        JLD2.jldsave(
            output_path;
            backend = String(backend),
            topological_sector = spec.sector_label,
            L = Int(L),
            τ_idx = Int(τ_idx),
            τ = τlis_ext[τ_idx],
            gamma = tanh(τlis_ext[τ_idx]),
            χ = isnothing(χ) ? 0 : Int(χ),
            seed = Int(seed),
            periods = periods,
            initial_bond_dimension = result.initial_bond_dimension,
            sample = result.samples,
            sample_free_energy = result.free_energy,
            halfchain_EE_tlis = result.halfchain_entropy,
            y_expectation_values = get(result, :y_expectation, Float32[]),
            final_EElis = result.final_entropy_profile,
        )
        return output_path
    end

    function run_sector_task(spec, task)
        backend, L, τ_idx, χ, seed, periods = task
        try
            result = if backend == :exact
                simulate_sector_exact(spec, L, τ_idx, seed; periods = periods)
            else
                simulate_sector_mps(spec, L, τ_idx, χ, seed; periods = periods)
            end
            output_path = save_sector_trajectory(
                spec,
                backend,
                L,
                τ_idx,
                seed,
                result;
                χ = backend == :mps ? χ : nothing,
            )
            return (backend, L, τ_idx, seed, :success, output_path)
        catch error
            return (backend, L, τ_idx, seed, :failed, error)
        end
    end

    function collect_sector_trajectories(
        spec,
        backend::Symbol,
        L::Integer,
        τ_idx::Integer;
        χ::Union{Nothing,Integer} = nothing,
        periods::Integer = sector_dyn_default_periods(backend, L, τ_idx),
    )
        data_dir = sector_dyn_trajectory_dir(spec, backend, L, τ_idx; χ = χ)
        files = sort(filter(
            file -> startswith(file, "periods$(periods)_trajectory_seed") &&
                    endswith(file, ".jld2"),
            readdir(data_dir),
        ))
        isempty(files) && error(
            "No $(spec.title) trajectory files with periods=$periods found in $data_dir",
        )

        first_data = JLD2.load(joinpath(data_dir, first(files)))
        Int(first_data["periods"]) == periods || error("Inconsistent evolution length")
        samples_num = length(files)
        entropy_dynamics = zeros(Float64, samples_num, periods)
        final_profiles = zeros(Float64, samples_num, L - 1)
        free_energies = Vector{Vector{Float64}}(undef, samples_num)
        seeds = zeros(Int, samples_num)

        for (i, file) in enumerate(files)
            data = JLD2.load(joinpath(data_dir, file))
            data["topological_sector"] == spec.sector_label ||
                error("Unexpected sector in $file")
            Int(data["periods"]) == periods ||
                error("Inconsistent evolution length in $file")
            entropy_dynamics[i, :] = Float64.(data["halfchain_EE_tlis"])
            final_profiles[i, :] = Float64.(data["final_EElis"])
            free_energies[i] = Float64.(data["sample_free_energy"])
            seeds[i] = Int(data["seed"])
        end
        length(unique(seeds)) == samples_num || error("Duplicate trajectory seeds found")
        all(length(fe) == length(first(free_energies)) for fe in free_energies) ||
            error("Inconsistent free-energy trajectory lengths")

        free_energy_matrix = hcat(free_energies...)
        average_EE_tlis = vec(mean(entropy_dynamics; dims = 1))
        stderr_EE_tlis =
            vec(std(entropy_dynamics; dims = 1, corrected = false)) ./ √samples_num
        bulk_meanEElis = vec(mean(final_profiles; dims = 1))
        ensemble_stderr_EElis =
            vec(std(final_profiles; dims = 1, corrected = false)) ./ √samples_num
        time_FElis = vec(mean(free_energy_matrix; dims = 2))
        time_FEstderr =
            vec(std(free_energy_matrix; dims = 2, corrected = false)) ./ √samples_num

        # Use the shared backend-specific average window from config.jl.  It
        # already fixes the Fibonacci-layer parity and stops before the
        # terminal half-strength layer generated by default
        # `enable_τ_eff=true`.
        total_layers = size(free_energy_matrix, 1)
        average_window = if backend == :exact
            get_cfg_params_Born(τ_idx, L)[3]
        elseif backend == :mps
            get_mps_params_Born(τ_idx, L)[3]
        else
            error("backend must be :exact or :mps")
        end
        isempty(average_window) &&
            error("The configured free-energy averaging window is empty")
        first(average_window) >= 1 && last(average_window) <= total_layers ||
            error(
                "The configured free-energy averaging window $(average_window) " *
                "lies outside 1:$total_layers for L=$L, periods=$periods",
            )
        trajectory_FE = vec(mean(free_energy_matrix[average_window, :]; dims = 1))
        bulk_FE = mean(trajectory_FE)
        bulk_FE_stderr = std(trajectory_FE; corrected = false) / √samples_num

        output_path = joinpath(
            data_dir,
            "EE_FEdynamics_L$(L)_gamma$(τ_idx)_periods$(periods)_$(spec.file_suffix).jld2",
        )
        JLD2.jldsave(
            output_path;
            backend = String(backend),
            topological_sector = spec.sector_label,
            L = Int(L),
            τ_idx = Int(τ_idx),
            τ = τlis_ext[τ_idx],
            gamma = tanh(τlis_ext[τ_idx]),
            χ = isnothing(χ) ? 0 : Int(χ),
            periods = periods,
            samples_num = samples_num,
            average_EE_tlis = average_EE_tlis,
            stderr_EE_tlis = stderr_EE_tlis,
            bulk_meanEElis = bulk_meanEElis,
            ensemble_stderr_EElis = ensemble_stderr_EElis,
            time_average_free_energy = trajectory_FE,
            bulk_FE = bulk_FE,
            bulk_FE_stderr = bulk_FE_stderr,
            time_FEstderr = time_FEstderr,
            time_FElis = time_FElis,
            ensemble_seed = seeds,
        )
        return output_path
    end
end

function print_sector_dyn_usage()
    script = "exm/Bulk_measure/monitored_dynamics_sector.jl"
    println("Usage:")
    println("  julia --project=. -p N $script SECTOR exact L τ_idx seed_start seed_end [periods]")
    println("  julia --project=. -p N $script SECTOR mps L τ_idx χ seed_start seed_end [periods]")
    println("  julia --project=. $script SECTOR collect_exact L τ_idx [periods]")
    println("  julia --project=. $script SECTOR collect_mps L τ_idx χ [periods]")
    println("  SECTOR is either y1 or ytau")
end

if isempty(ARGS)
    print_sector_dyn_usage()
else
    length(ARGS) >= 2 || (print_sector_dyn_usage(); error("Not enough arguments"))
    spec = sector_dyn_spec(Symbol(lowercase(ARGS[1])))
    action = Symbol(lowercase(ARGS[2]))
    if action ∈ (:exact, :mps)
        minimum_args = action == :exact ? 6 : 7
        length(ARGS) >= minimum_args || error("Not enough arguments for $action")
        L = parse(Int, ARGS[3])
        τ_idx = parse(Int, ARGS[4])
        χ = action == :mps ? parse(Int, ARGS[5]) : 0
        seed_offset = action == :mps ? 5 : 4
        seed_start = parse(Int, ARGS[seed_offset + 1])
        seed_end = parse(Int, ARGS[seed_offset + 2])
        default_periods = sector_dyn_default_periods(action, L, τ_idx)
        periods = length(ARGS) >= seed_offset + 3 ?
                  parse(Int, ARGS[seed_offset + 3]) : default_periods
        periods >= 1 || error("periods must be positive")

        tasks = [(action, L, τ_idx, χ, seed, periods) for seed in seed_start:seed_end]
        println("=== $(spec.title) monitored dynamics ($(action)) ===")
        println("L=$L, τ_idx=$τ_idx, periods=$periods, trajectories=$(length(tasks))")
        action == :mps && println("χ=$χ")
        println("workers=$(nworkers())")
        results = pmap(task -> run_sector_task(spec, task), tasks; batch_size = 1)
        failed = filter(result -> result[5] != :success, results)
        println("Successes: $(length(results) - length(failed))")
        println("Failures: $(length(failed))")
        for result in failed
            println("Failed seed $(result[4]): $(result[6])")
        end
    elseif action == :collect_exact
        length(ARGS) ∈ (4, 5) || error("collect_exact requires L τ_idx [periods]")
        L = parse(Int, ARGS[3])
        τ_idx = parse(Int, ARGS[4])
        periods = length(ARGS) == 5 ?
                  parse(Int, ARGS[5]) : sector_dyn_default_periods(:exact, L, τ_idx)
        output_path = collect_sector_trajectories(
            spec,
            :exact,
            L,
            τ_idx;
            periods = periods,
        )
        println("Saved $output_path")
    elseif action == :collect_mps
        length(ARGS) ∈ (5, 6) || error("collect_mps requires L τ_idx χ [periods]")
        L = parse(Int, ARGS[3])
        τ_idx = parse(Int, ARGS[4])
        χ = parse(Int, ARGS[5])
        periods = length(ARGS) == 6 ?
                  parse(Int, ARGS[6]) : sector_dyn_default_periods(:mps, L, τ_idx)
        output_path = collect_sector_trajectories(
            spec,
            :mps,
            L,
            τ_idx;
            χ = χ,
            periods = periods,
        )
        println("Saved $output_path")
    else
        print_sector_dyn_usage()
        error("Unknown action: $(ARGS[2])")
    end
end
