using Distributed
using FibonacciChain
using ITensorMPS
using ITensors
using JLD2
using LinearAlgebra
using Random
using Statistics

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

    const YTAU_DYNAMICS_DATA_ROOT =
        joinpath("exm", "data", "Bulk_measure", "monitored_dynamics_ytau")

    """Return the normalized projection of the all-τ path into the y=τ sector."""
    function initial_ytau_exact(model)
        model.pbc || error("The topological y=τ sector requires periodic boundaries")
        basis = anyon_basis(model)
        state = zeros(Float64, length(basis))
        state[1] = 1.0

        ϕ = (1 + √5) / 2
        yτ = -inv(ϕ)
        # Y is symmetric, so the row returned by topological_symmetry_basismap
        # is also the single column Y|00⋯0⟩. Avoid constructing the full dense
        # Y matrix, whose size grows exponentially with L.
        y_on_state = topological_symmetry_basismap(model, first(basis))
        # Y has eigenvalues y₁=ϕ and yτ=-1/ϕ, hence
        # Pτ = (y₁ I - Y)/(y₁-yτ).
        projected = (ϕ .* state .- y_on_state) ./ (ϕ - yτ)
        projected_norm = norm(projected)
        projected_norm > 1e-14 || error("The initial state has zero y=τ weight")

        y1_weight = (1 + ϕ * (-inv(ϕ))^model.N) / (1 + ϕ^2)
        expected_weight = 1 - y1_weight
        isapprox(projected_norm^2, expected_weight; atol = 1e-10, rtol = 1e-10) ||
            error(
                "Unexpected y=τ weight $(projected_norm^2); expected $expected_weight",
            )
        projected ./= projected_norm
        return projected
    end

    # Matrix element of Y|00⋯0⟩ for two adjacent output fusion-path bits.
    # In the repository convention 0=τ and 1=vacuum.
    function _y_vacuum_edge_weight(x::Int, xnext::Int)
        ϕ = (1 + √5) / 2
        x == 0 && xnext == 0 && return -inv(ϕ)
        x == 0 && xnext == 1 && return inv(sqrt(ϕ))
        x == 1 && xnext == 0 && return 1.0
        return 0.0
    end

    """
        initial_ytau_mps(model; cutoff=1e-14, maxdim=32)

    Construct the same projected all-τ state as [`initial_ytau_exact`](@ref),
    directly as an MPS. The amplitudes of `Y|00⋯0⟩` form a periodic
    nearest-neighbor weight with bond dimension two. Resolving its two boundary
    values and combining it with `ϕ|00⋯0⟩` implements

        Pτ|00⋯0⟩ ∝ (ϕ I - Y)|00⋯0⟩.

    This avoids constructing the exponentially large dense topological-charge
    operator and remains suitable for the large-L MPS simulations.
    """
    function initial_ytau_mps(model; cutoff::Float64 = 1e-14, maxdim::Int = 32)
        model.pbc || error("The topological y=τ sector requires periodic boundaries")
        N = model.N
        N >= 2 || error("initial_ytau_mps requires at least two sites")
        sites = siteinds("Qubit", N)

        function fixed_boundary_mps(boundary_bit::Int)
            links = [Index(2, "Link,ytau_init,l=$i") for i in 1:(N - 1)]
            tensors = Vector{ITensor}(undef, N)

            first_tensor = ITensor(sites[1], links[1])
            for next_bit in 0:1
                first_tensor[
                    sites[1] => boundary_bit + 1,
                    links[1] => next_bit + 1,
                ] = _y_vacuum_edge_weight(boundary_bit, next_bit)
            end
            tensors[1] = first_tensor

            for site in 2:(N - 1)
                tensor = ITensor(links[site - 1], sites[site], links[site])
                for bit in 0:1, next_bit in 0:1
                    tensor[
                        links[site - 1] => bit + 1,
                        sites[site] => bit + 1,
                        links[site] => next_bit + 1,
                    ] = _y_vacuum_edge_weight(bit, next_bit)
                end
                tensors[site] = tensor
            end

            last_tensor = ITensor(links[end], sites[end])
            for bit in 0:1
                last_tensor[
                    links[end] => bit + 1,
                    sites[end] => bit + 1,
                ] = _y_vacuum_edge_weight(bit, boundary_bit)
            end
            tensors[end] = last_tensor
            return MPS(tensors)
        end

        y_on_vacuum = add(
            fixed_boundary_mps(0),
            fixed_boundary_mps(1);
            cutoff = cutoff,
            maxdim = maxdim,
        )
        vacuum = productMPS(sites, fill("0", N))
        ϕ = (1 + √5) / 2
        y_on_vacuum[1] = -y_on_vacuum[1]
        vacuum[1] = ϕ * vacuum[1]
        projected = add(y_on_vacuum, vacuum; cutoff = cutoff, maxdim = maxdim)
        normalize!(projected)
        return projected, sites
    end

    function ytau_trajectory_dir(
        backend::Symbol,
        L::Integer,
        τ_idx::Integer;
        χ::Union{Nothing,Integer} = nothing,
    )
        backend ∈ (:exact, :mps) || error("backend must be :exact or :mps")
        backend == :mps && isnothing(χ) && error("χ is required for the MPS backend")
        path = joinpath(
            YTAU_DYNAMICS_DATA_ROOT,
            String(backend),
            "L$(L)",
            "gammaind$(τ_idx)",
        )
        return backend == :mps ? joinpath(path, "chi$(χ)") : path
    end

    function ytau_default_periods(backend::Symbol, L::Integer, τ_idx::Integer)
        if backend == :exact
            D, _, _ = get_cfg_params_Born(τ_idx, L)
            return div(D, 2)
        elseif backend == :mps
            time_in_L, _, _ = get_mps_params_Born(τ_idx, L)
            return time_in_L * L
        end
        error("backend must be :exact or :mps")
    end

    function simulate_ytau_exact(
        L::Integer,
        τ_idx::Integer,
        seed::Integer;
        periods::Integer = ytau_default_periods(:exact, L, τ_idx),
    )
        iseven(L) || error("L must be even, got $L")
        τ = τlis_ext[τ_idx]
        model = fib_model(L)
        initial_state = initial_ytau_exact(model)
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

    function simulate_ytau_mps(
        L::Integer,
        τ_idx::Integer,
        χ::Integer,
        seed::Integer;
        periods::Integer = ytau_default_periods(:mps, L, τ_idx),
    )
        iseven(L) || error("L must be even, got $L")
        χ >= 5 || error("χ must be at least 5 to represent the exact y=τ initial MPS")
        τ = τlis_ext[τ_idx]
        model = fib_model(L)
        initial_state, sites = initial_ytau_mps(model; maxdim = χ)
        initial_bond_dimension = maxlinkdim(initial_state)
        config = MeasureConfig(
            τ = τ,
            mode = :Born,
            t₂ = periods,
            rng = MersenneTwister(seed),
            cutoff = 1e-12,
            maxdim = χ,
            enforce_fibonacci_constraint = true,
        )
        outcome = bulk_evolution(model, sites, initial_state, config)
        return (
            state = outcome.state,
            samples = outcome.samples,
            free_energy = outcome.free_energys,
            halfchain_entropy = outcome.entanglement_entropys,
            final_entropy_profile = anyon_eelis(model, outcome.state),
            initial_bond_dimension = initial_bond_dimension,
        )
    end

    function save_ytau_trajectory(
        backend::Symbol,
        L::Integer,
        τ_idx::Integer,
        seed::Integer,
        result;
        χ::Union{Nothing,Integer} = nothing,
    )
        out_dir = ytau_trajectory_dir(backend, L, τ_idx; χ = χ)
        mkpath(out_dir)
        periods = length(result.halfchain_entropy)
        output_path =
            joinpath(out_dir, "periods$(periods)_trajectory_seed$(seed).jld2")
        JLD2.jldsave(
            output_path;
            backend = String(backend),
            topological_sector = "y=tau",
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
            final_EElis = result.final_entropy_profile,
        )
        return output_path
    end

    function run_ytau_task(task)
        backend, L, τ_idx, χ, seed, periods = task
        try
            result = if backend == :exact
                simulate_ytau_exact(L, τ_idx, seed; periods = periods)
            else
                simulate_ytau_mps(L, τ_idx, χ, seed; periods = periods)
            end
            output_path = save_ytau_trajectory(
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

    function collect_ytau_trajectories(
        backend::Symbol,
        L::Integer,
        τ_idx::Integer;
        χ::Union{Nothing,Integer} = nothing,
        periods::Integer = ytau_default_periods(backend, L, τ_idx),
    )
        data_dir = ytau_trajectory_dir(backend, L, τ_idx; χ = χ)
        files = sort(filter(
            file -> startswith(file, "periods$(periods)_trajectory_seed") &&
                    endswith(file, ".jld2"),
            readdir(data_dir),
        ))
        isempty(files) &&
            error("No y=τ trajectory files with periods=$periods found in $data_dir")

        first_data = JLD2.load(joinpath(data_dir, first(files)))
        Int(first_data["periods"]) == periods || error("Inconsistent evolution length")
        samples_num = length(files)
        entropy_dynamics = zeros(Float64, samples_num, periods)
        final_profiles = zeros(Float64, samples_num, L - 1)
        free_energies = Vector{Vector{Float64}}(undef, samples_num)
        seeds = zeros(Int, samples_num)

        for (i, file) in enumerate(files)
            data = JLD2.load(joinpath(data_dir, file))
            data["topological_sector"] == "y=tau" ||
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

        # Use the latter half of the record for a backend-independent steady-state
        # estimate. The full time series is saved for alternative fitting windows.
        steady_window = (div(size(free_energy_matrix, 1), 2) + 1):size(free_energy_matrix, 1)
        trajectory_FE = vec(mean(free_energy_matrix[steady_window, :]; dims = 1))
        bulk_FE = mean(trajectory_FE)
        bulk_FE_stderr = std(trajectory_FE; corrected = false) / √samples_num

        output_path = joinpath(
            data_dir,
            "EE_FEdynamics_L$(L)_gamma$(τ_idx)_periods$(periods)_ytau.jld2",
        )
        JLD2.jldsave(
            output_path;
            backend = String(backend),
            topological_sector = "y=tau",
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

function print_ytau_usage()
    println("Usage:")
    println("  julia --project=. -p N exm/Bulk_measure/monitored_dynamics_ytau.jl exact L τ_idx seed_start seed_end [periods]")
    println("  julia --project=. -p N exm/Bulk_measure/monitored_dynamics_ytau.jl mps L τ_idx χ seed_start seed_end [periods]")
    println("  julia --project=. exm/Bulk_measure/monitored_dynamics_ytau.jl collect_exact L τ_idx [periods]")
    println("  julia --project=. exm/Bulk_measure/monitored_dynamics_ytau.jl collect_mps L τ_idx χ [periods]")
end

if isempty(ARGS)
    print_ytau_usage()
else
    action = Symbol(lowercase(ARGS[1]))
    if action ∈ (:exact, :mps)
        minimum_args = action == :exact ? 5 : 6
        length(ARGS) >= minimum_args || error("Not enough arguments for $action")
        L = parse(Int, ARGS[2])
        τ_idx = parse(Int, ARGS[3])
        χ = action == :mps ? parse(Int, ARGS[4]) : 0
        seed_offset = action == :mps ? 4 : 3
        seed_start = parse(Int, ARGS[seed_offset + 1])
        seed_end = parse(Int, ARGS[seed_offset + 2])
        default_periods = ytau_default_periods(action, L, τ_idx)
        periods = length(ARGS) >= seed_offset + 3 ?
                  parse(Int, ARGS[seed_offset + 3]) : default_periods
        periods >= 1 || error("periods must be positive")

        tasks = [(action, L, τ_idx, χ, seed, periods) for seed in seed_start:seed_end]
        println("=== y=τ monitored dynamics ($(action)) ===")
        println("L=$L, τ_idx=$τ_idx, periods=$periods, trajectories=$(length(tasks))")
        action == :mps && println("χ=$χ")
        println("workers=$(nworkers())")
        results = pmap(run_ytau_task, tasks; batch_size = 1)
        failed = filter(result -> result[5] != :success, results)
        println("Successes: $(length(results) - length(failed))")
        println("Failures: $(length(failed))")
        for result in failed
            println("Failed seed $(result[4]): $(result[6])")
        end
    elseif action == :collect_exact
        length(ARGS) ∈ (3, 4) || error("collect_exact requires L τ_idx [periods]")
        L = parse(Int, ARGS[2])
        τ_idx = parse(Int, ARGS[3])
        periods = length(ARGS) == 4 ?
                  parse(Int, ARGS[4]) : ytau_default_periods(:exact, L, τ_idx)
        output_path = collect_ytau_trajectories(
            :exact,
            L,
            τ_idx;
            periods = periods,
        )
        println("Saved $output_path")
    elseif action == :collect_mps
        length(ARGS) ∈ (4, 5) || error("collect_mps requires L τ_idx χ [periods]")
        L = parse(Int, ARGS[2])
        τ_idx = parse(Int, ARGS[3])
        χ = parse(Int, ARGS[4])
        periods = length(ARGS) == 5 ?
                  parse(Int, ARGS[5]) : ytau_default_periods(:mps, L, τ_idx)
        output_path = collect_ytau_trajectories(
            :mps,
            L,
            τ_idx;
            χ = χ,
            periods = periods,
        )
        println("Saved $output_path")
    else
        print_ytau_usage()
        error("Unknown action: $(ARGS[1])")
    end
end
