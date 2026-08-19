using Distributed
using FibonacciChain
using ITensorMPS
using ITensors
using JLD2
using LinearAlgebra
using Random
using Statistics
using Arpack

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
    using Arpack

    const TCI_DYNAMICS_DATA_ROOT =
        joinpath("exm", "data", "Bulk_measure", "monitored_dynamics_tci")
    const TCI_DMRG_SWEEPS = 20
    const TCI_DMRG_CUTOFF = 1e-12
    const TCI_CONSTRAINT_PENALTY = 20.0
    const TCI_DMRG_SEED = 271828

    # Ground states are identical for all trajectories with the same parameters.
    # Cache one copy per worker instead of repeating diagonalization/DMRG per seed.
    const TCI_EXACT_GS_CACHE = Dict{Int,Any}()
    const TCI_MPS_GS_CACHE = Dict{Tuple{Int,Int,Int,Float64,Float64,Int},Any}()

    """Return the TCI ground state obtained from the dense anyon Hamiltonian."""
    function initial_tci_exact(model)
        model.pbc || error("The TCI initial state requires periodic boundaries")
        model.measure_operator == :Antiferro ||
            error("The TCI initial state requires measure_operator=:Antiferro")
        cached = get!(TCI_EXACT_GS_CACHE, model.N) do
            # `anyon_ham(model)` is the repository's dense-Hamiltonian API.
            energies, states = Arpack.eigs(anyon_ham_sparse(model), which = :SR, nev = 1)
            (; state = vec(states[:, 1]), energy = real(energies[1]))
        end
        return copy(cached.state), cached.energy
    end

    """Return the cached constrained-DMRG TCI ground state."""
    function initial_tci_mps(
        model;
        maxdim::Int,
        sweep_times::Int = TCI_DMRG_SWEEPS,
        cutoff::Float64 = TCI_DMRG_CUTOFF,
        constraint_penalty::Float64 = TCI_CONSTRAINT_PENALTY,
        dmrg_seed::Int = TCI_DMRG_SEED,
    )
        key = (
            model.N,
            maxdim,
            sweep_times,
            cutoff,
            constraint_penalty,
            dmrg_seed,
        )
        cached = get!(TCI_MPS_GS_CACHE, key) do
            state, energy = anyon_mps_gst(
                model;
                maxdim = maxdim,
                sweep_times = sweep_times,
                cutoff = cutoff,
                constraint_penalty = constraint_penalty,
                seed = dmrg_seed,
                outputlevel = 0,
            )
            violation = fibonacci_constraint_violation(state; pbc = model.pbc)
            (; state, sites = siteinds(state), energy, violation)
        end
        return (
            copy(cached.state),
            copy(cached.sites),
            cached.energy,
            cached.violation,
        )
    end

    function tci_trajectory_dir(
        backend::Symbol,
        L::Integer,
        tau_idx::Integer;
        chi::Union{Nothing,Integer} = nothing,
    )
        backend in (:exact, :mps) || error("backend must be :exact or :mps")
        backend == :mps && isnothing(chi) && error("chi is required for the MPS backend")
        source = backend == :exact ? "dense_eigensolver" : "dmrg"
        path = joinpath(
            TCI_DYNAMICS_DATA_ROOT,
            String(backend),
            source,
            "L$(L)",
            "gammaind$(tau_idx)",
        )
        return backend == :mps ? joinpath(path, "chi$(chi)") : path
    end

    function tci_default_periods(backend::Symbol, L::Integer, tau_idx::Integer)
        if backend == :exact
            D, _, _ = get_cfg_params_Born(tau_idx, L)
            return div(D, 2)
        elseif backend == :mps
            time_in_L, _, _ = get_mps_params_Born(tau_idx, L)
            return time_in_L * L
        end
        error("backend must be :exact or :mps")
    end

    function simulate_tci_exact(
        L::Integer,
        tau_idx::Integer,
        seed::Integer;
        periods::Integer = tci_default_periods(:exact, L, tau_idx),
    )
        iseven(L) || error("L must be even, got $L")
        τ = τlis_ext[tau_idx]
        model = fib_model(L)
        initial_state, initial_energy = initial_tci_exact(model)
        config = MeasureConfig(
            τ = τ,
            mode = :Born,
            t₂ = periods,
            rng = MersenneTwister(seed),
        )
        outcome = bulk_evolution(model, initial_state, config)
        return (
            samples = outcome.samples,
            free_energy = outcome.free_energys,
            halfchain_entropy = outcome.entanglement_entropys,
            final_entropy_profile = anyon_eelis(model, outcome.state),
            initial_energy = initial_energy,
            initial_mpo_energy = initial_energy,
            initial_constraint_violation = 0.0,
            initial_bond_dimension = 0,
        )
    end

    function simulate_tci_mps(
        L::Integer,
        tau_idx::Integer,
        chi::Integer,
        seed::Integer;
        periods::Integer = tci_default_periods(:mps, L, tau_idx),
    )
        iseven(L) || error("L must be even, got $L")
        τ = τlis_ext[tau_idx]
        model = fib_model(L)
        initial_state, sites, initial_energy, violation =
            initial_tci_mps(model; maxdim = chi)
        initial_bond_dimension = maxlinkdim(initial_state)
        config = MeasureConfig(
            τ = τ,
            mode = :Born,
            t₂ = periods,
            rng = MersenneTwister(seed),
            cutoff = 1e-12,
            maxdim = chi,
            enforce_fibonacci_constraint = true,
        )
        outcome = bulk_evolution(model, sites, initial_state, config)
        return (
            samples = outcome.samples,
            free_energy = outcome.free_energys,
            halfchain_entropy = outcome.entanglement_entropys,
            final_entropy_profile = anyon_eelis(model, outcome.state),
            initial_energy = initial_energy,
            initial_mpo_energy = initial_energy,
            initial_constraint_violation = violation,
            initial_bond_dimension = initial_bond_dimension,
        )
    end

    function save_tci_trajectory(
        backend::Symbol,
        L::Integer,
        tau_idx::Integer,
        seed::Integer,
        result;
        chi::Union{Nothing,Integer} = nothing,
    )
        out_dir = tci_trajectory_dir(backend, L, tau_idx; chi = chi)
        mkpath(out_dir)
        periods = length(result.halfchain_entropy)
        output_path =
            joinpath(out_dir, "periods$(periods)_trajectory_seed$(seed).jld2")
        JLD2.jldsave(
            output_path;
            backend = String(backend),
            initial_state = "TCI_GS",
            initial_state_method = backend == :exact ? "dense_eigensolver" : "dmrg",
            L = Int(L),
            tau_idx = Int(tau_idx),
            tau = τlis_ext[tau_idx],
            gamma = tanh(τlis_ext[tau_idx]),
            chi = isnothing(chi) ? 0 : Int(chi),
            seed = Int(seed),
            periods = periods,
            dmrg_sweeps = backend == :mps ? TCI_DMRG_SWEEPS : 0,
            dmrg_cutoff = backend == :mps ? TCI_DMRG_CUTOFF : 0.0,
            dmrg_seed = backend == :mps ? TCI_DMRG_SEED : 0,
            constraint_penalty = backend == :mps ? TCI_CONSTRAINT_PENALTY : 0.0,
            initial_energy = result.initial_energy,
            initial_mpo_energy = result.initial_mpo_energy,
            initial_constraint_violation = result.initial_constraint_violation,
            initial_bond_dimension = result.initial_bond_dimension,
            sample = result.samples,
            sample_free_energy = result.free_energy,
            halfchain_EE_tlis = result.halfchain_entropy,
            final_EElis = result.final_entropy_profile,
        )
        return output_path
    end

    function run_tci_task(task)
        backend, L, tau_idx, chi, seed, periods = task
        try
            result = if backend == :exact
                simulate_tci_exact(L, tau_idx, seed; periods = periods)
            else
                simulate_tci_mps(L, tau_idx, chi, seed; periods = periods)
            end
            output_path = save_tci_trajectory(
                backend,
                L,
                tau_idx,
                seed,
                result;
                chi = backend == :mps ? chi : nothing,
            )
            return (backend, L, tau_idx, seed, :success, output_path)
        catch err
            return (backend, L, tau_idx, seed, :failed, err)
        end
    end

    function collect_tci_trajectories(
        backend::Symbol,
        L::Integer,
        tau_idx::Integer;
        chi::Union{Nothing,Integer} = nothing,
        periods::Integer = tci_default_periods(backend, L, tau_idx),
    )
        data_dir = tci_trajectory_dir(backend, L, tau_idx; chi = chi)
        isdir(data_dir) || error("TCI trajectory directory does not exist: $data_dir")
        files = sort(filter(
            file -> startswith(file, "periods$(periods)_trajectory_seed") &&
                    endswith(file, ".jld2"),
            readdir(data_dir),
        ))
        isempty(files) &&
            error("No TCI trajectory files with periods=$periods found in $data_dir")

        samples_num = length(files)
        entropy_dynamics = zeros(Float64, samples_num, periods)
        final_profiles = zeros(Float64, samples_num, L - 1)
        free_energies = Vector{Vector{Float64}}(undef, samples_num)
        seeds = zeros(Int, samples_num)
        initial_energies = zeros(Float64, samples_num)
        initial_mpo_energies = zeros(Float64, samples_num)
        initial_violations = zeros(Float64, samples_num)
        initial_bond_dimensions = zeros(Int, samples_num)
        expected_method = backend == :exact ? "dense_eigensolver" : "dmrg"

        for (i, file) in enumerate(files)
            data = JLD2.load(joinpath(data_dir, file))
            data["initial_state"] == "TCI_GS" || error("Unexpected initial state in $file")
            data["initial_state_method"] == expected_method ||
                error("Unexpected initial-state method in $file")
            data["backend"] == String(backend) || error("Unexpected backend in $file")
            Int(data["L"]) == L || error("Inconsistent L in $file")
            Int(data["tau_idx"]) == tau_idx || error("Inconsistent tau_idx in $file")
            Int(data["periods"]) == periods ||
                error("Inconsistent evolution length in $file")
            backend == :mps && Int(data["chi"]) != chi &&
                error("Inconsistent chi in $file")
            if backend == :mps
                Int(data["dmrg_sweeps"]) == TCI_DMRG_SWEEPS ||
                    error("Inconsistent DMRG sweep count in $file")
                Float64(data["dmrg_cutoff"]) == TCI_DMRG_CUTOFF ||
                    error("Inconsistent DMRG cutoff in $file")
                Int(data["dmrg_seed"]) == TCI_DMRG_SEED ||
                    error("Inconsistent DMRG seed in $file")
                Float64(data["constraint_penalty"]) == TCI_CONSTRAINT_PENALTY ||
                    error("Inconsistent constraint penalty in $file")
            end
            entropy_dynamics[i, :] = Float64.(data["halfchain_EE_tlis"])
            final_profiles[i, :] = Float64.(data["final_EElis"])
            free_energies[i] = Float64.(data["sample_free_energy"])
            seeds[i] = Int(data["seed"])
            initial_energies[i] = Float64(data["initial_energy"])
            initial_mpo_energies[i] = Float64(data["initial_mpo_energy"])
            initial_violations[i] = Float64(data["initial_constraint_violation"])
            initial_bond_dimensions[i] = Int(data["initial_bond_dimension"])
        end
        length(unique(seeds)) == samples_num || error("Duplicate trajectory seeds found")
        all(length(fe) == length(first(free_energies)) for fe in free_energies) ||
            error("Inconsistent free-energy trajectory lengths")
        maximum(initial_energies) - minimum(initial_energies) <= 1e-6 ||
            error("Trajectories used inconsistent TCI ground states")
        length(first(free_energies)) == 2 * periods ||
            error("A Fibonacci trajectory must contain two free-energy entries per period")

        free_energy_matrix = hcat(free_energies...)
        average_EE_tlis = vec(mean(entropy_dynamics; dims = 1))
        stderr_EE_tlis =
            vec(std(entropy_dynamics; dims = 1, corrected = false)) ./ sqrt(samples_num)
        bulk_meanEElis = vec(mean(final_profiles; dims = 1))
        ensemble_stderr_EElis =
            vec(std(final_profiles; dims = 1, corrected = false)) ./ sqrt(samples_num)
        time_FElis = vec(mean(free_energy_matrix; dims = 2))
        time_FEstderr =
            vec(std(free_energy_matrix; dims = 2, corrected = false)) ./ sqrt(samples_num)

        steady_window =
            (div(size(free_energy_matrix, 1), 2) + 1):size(free_energy_matrix, 1)
        trajectory_FE = vec(mean(free_energy_matrix[steady_window, :]; dims = 1))
        bulk_FE = mean(trajectory_FE)
        bulk_FE_stderr = std(trajectory_FE; corrected = false) / sqrt(samples_num)

        output_path = joinpath(
            data_dir,
            "EE_FEdynamics_L$(L)_gamma$(tau_idx)_periods$(periods)_tci.jld2",
        )
        JLD2.jldsave(
            output_path;
            backend = String(backend),
            initial_state = "TCI_GS",
            initial_state_method = expected_method,
            L = Int(L),
            tau_idx = Int(tau_idx),
            tau = τlis_ext[tau_idx],
            gamma = tanh(τlis_ext[tau_idx]),
            chi = isnothing(chi) ? 0 : Int(chi),
            periods = periods,
            samples_num = samples_num,
            initial_energy = mean(initial_energies),
            initial_mpo_energy = mean(initial_mpo_energies),
            maximum_initial_constraint_violation = maximum(initial_violations),
            maximum_initial_bond_dimension = maximum(initial_bond_dimensions),
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

function print_tci_usage()
    println("Usage:")
    println("  julia --project=. -p N exm/Bulk_measure/monitored_dynamics_tci.jl exact L tau_idx seed_start seed_end [periods]")
    println("  julia --project=. -p N exm/Bulk_measure/monitored_dynamics_tci.jl mps L tau_idx chi seed_start seed_end [periods]")
    println("  julia --project=. exm/Bulk_measure/monitored_dynamics_tci.jl collect_exact L tau_idx [periods]")
    println("  julia --project=. exm/Bulk_measure/monitored_dynamics_tci.jl collect_mps L tau_idx chi [periods]")
end

if isempty(ARGS)
    print_tci_usage()
else
    action = Symbol(lowercase(ARGS[1]))
    if action in (:exact, :mps)
        minimum_args = action == :exact ? 5 : 6
        length(ARGS) >= minimum_args || error("Not enough arguments for $action")
        L = parse(Int, ARGS[2])
        tau_idx = parse(Int, ARGS[3])
        chi = action == :mps ? parse(Int, ARGS[4]) : 0
        seed_offset = action == :mps ? 4 : 3
        seed_start = parse(Int, ARGS[seed_offset + 1])
        seed_end = parse(Int, ARGS[seed_offset + 2])
        seed_start <= seed_end || error("seed_start must not exceed seed_end")
        default_periods = tci_default_periods(action, L, tau_idx)
        periods = length(ARGS) >= seed_offset + 3 ?
                  parse(Int, ARGS[seed_offset + 3]) : default_periods
        periods >= 1 || error("periods must be positive")

        tasks = [(action, L, tau_idx, chi, seed, periods) for seed in seed_start:seed_end]
        println("=== TCI-ground-state monitored dynamics ($(action)) ===")
        println("L=$L, tau_idx=$tau_idx, periods=$periods, trajectories=$(length(tasks))")
        action == :mps && println(
            "chi=$chi, DMRG sweeps=$TCI_DMRG_SWEEPS, constraint penalty=$TCI_CONSTRAINT_PENALTY",
        )
        println("workers=$(nworkers())")
        results = pmap(run_tci_task, tasks; batch_size = 1)
        failed = filter(result -> result[5] != :success, results)
        println("Successes: $(length(results) - length(failed))")
        println("Failures: $(length(failed))")
        for result in failed
            println("Failed seed $(result[4]): $(result[6])")
        end
    elseif action == :collect_exact
        length(ARGS) in (3, 4) || error("collect_exact requires L tau_idx [periods]")
        L = parse(Int, ARGS[2])
        tau_idx = parse(Int, ARGS[3])
        periods = length(ARGS) == 4 ?
                  parse(Int, ARGS[4]) : tci_default_periods(:exact, L, tau_idx)
        println("Saved $(collect_tci_trajectories(:exact, L, tau_idx; periods = periods))")
    elseif action == :collect_mps
        length(ARGS) in (4, 5) || error("collect_mps requires L tau_idx chi [periods]")
        L = parse(Int, ARGS[2])
        tau_idx = parse(Int, ARGS[3])
        chi = parse(Int, ARGS[4])
        periods = length(ARGS) == 5 ?
                  parse(Int, ARGS[5]) : tci_default_periods(:mps, L, tau_idx)
        println("Saved $(collect_tci_trajectories(:mps, L, tau_idx; chi = chi, periods = periods))")
    else
        print_tci_usage()
        error("Unknown action: $(ARGS[1])")
    end
end
