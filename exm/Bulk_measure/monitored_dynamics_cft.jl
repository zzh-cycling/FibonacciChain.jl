using Distributed
using FibonacciChain
using ITensorMPS
using ITensors
using JLD2
using LinearAlgebra
using Random
using Statistics
using Arpack

# Monitored dynamics starting from a CFT ground state (TCI or 3-state Potts).
# Merged from monitored_dynamics_tci.jl and monitored_dynamics_potts.jl;
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
    using Arpack

    const CFT_DMRG_CUTOFF = 1e-12
    const CFT_CONSTRAINT_PENALTY = 20.0
    const CFT_DMRG_SEED = 271828

    # Per-CFT settings.  Everything not listed here is shared between the two
    # CFTs: the dynamics always uses the Antiferro measurement operator, the
    # DMRG sweep schedule (maxdim ramp, cutoff, noise) is identical, and the
    # file formats are the pre-merge ones.
    const CFT_SPECS = Dict{Symbol,NamedTuple}(
        # TCI: Antiferro ground state, c = 7/10, checked at every size.
        # Fixed 3L DMRG sweeps (energy_tol = 0 disables early stopping).
        :tci => (;
            title = "TCI-ground-state",
            data_root =
                joinpath("exm", "data", "Bulk_measure", "monitored_dynamics_tci"),
            file_suffix = "tci",
            initial_state_tag = "TCI_GS",
            initial_measure_operator = :Antiferro,
            expected_central_charge = 0.7,
            central_charge_atol = 0.02,
            min_central_charge_check_size = 0,
            entanglement_mincut = L -> div(L, 4),
            dmrg_min_sweeps = nothing,
            dmrg_energy_tol = 0.0,
        ),
        # Potts: Ferro ground state, c = 4/5, checked only for L >= 14.
        # Convergence scans at chi=150 gave 72, 80, and 72 sweeps for
        # L=28, 30, and 32 at energy_tol=1e-11: use 2L as the minimum before
        # early stopping and retain 3L as a conservative cap.
        :potts => (;
            title = "Potts-ground-state",
            data_root =
                joinpath("exm", "data", "Bulk_measure", "monitored_dynamics_potts"),
            file_suffix = "potts",
            initial_state_tag = "POTTS_GS",
            initial_measure_operator = :Ferro,
            expected_central_charge = 0.8,
            central_charge_atol = 0.05,
            min_central_charge_check_size = 14,
            # The entropy profiles form two distinct finite-size sequences.
            # For L mod 3 = 0, firstcut=L/6 lies on a broad bulk-fit plateau.
            # For the L mod 3 = 1,2 sequence there is no deep-cut plateau;
            # firstcut=2 is the stable estimator approaching c=4/5.
            entanglement_mincut = L -> (
                L >= 4 || error("L must be at least 4 for a central-charge fit");
                L % 3 == 0 ? max(2, div(L, 6)) : 2
            ),
            dmrg_min_sweeps = L -> 2L,
            dmrg_energy_tol = 1e-11,
        ),
    )

    function cft_spec(cft::Symbol)
        haskey(CFT_SPECS, cft) || error("Unknown CFT: $cft (expected :tci or :potts)")
        return CFT_SPECS[cft]
    end

    # Ground states are identical for all trajectories with the same
    # parameters.  Cache one copy per worker instead of repeating
    # diagonalization/DMRG per seed.
    const CFT_EXACT_GS_CACHE = Dict{Tuple{Symbol,Int},Any}()
    const CFT_MPS_GS_CACHE = Dict{Tuple{Symbol,Int,Int,Int,Int,Float64,Float64,Float64,Int},Any}()

    cft_dmrg_sweeps(L::Integer) = 3L
    cft_dmrg_min_sweeps(spec, L::Integer) =
        isnothing(spec.dmrg_min_sweeps) ? cft_dmrg_sweeps(L) : spec.dmrg_min_sweeps(L)

    """Run constrained DMRG with an energy-based stopping criterion
    (`energy_tol = 0` disables early stopping, giving fixed `sweep_times`)."""
    function cft_mps_gst(
        model;
        sweep_times::Int,
        minimum_sweeps::Int,
        maxdim::Int,
        cutoff::Float64,
        constraint_penalty::Float64,
        energy_tol::Float64,
        seed::Int,
    )
        1 <= minimum_sweeps <= sweep_times ||
            error("minimum_sweeps must lie between 1 and sweep_times")
        sites = siteinds("Qubit", model.N)
        H = anyon_ham(model, sites)
        violation_mpo =
            FibonacciChain.fibonacci_constraint_violation_mpo(sites; pbc = model.pbc)
        penalty_mpo = copy(violation_mpo)
        penalty_mpo[1] *= constraint_penalty
        constrained_H = add(H, penalty_mpo; cutoff = min(cutoff, 1e-14))

        rng = MersenneTwister(seed)
        state = random_mps(rng, sites; linkdims = min(8, maxdim))
        projector =
            FibonacciChain.fibonacci_constraint_projector_mpo(sites; pbc = model.pbc)
        state = apply(
            projector,
            state;
            cutoff = min(cutoff, 1e-14),
            maxdim = maxdim,
        )
        normalize!(state)

        sweeps = Sweeps(sweep_times)
        setmaxdim!(sweeps, min(16, maxdim), min(32, maxdim), maxdim)
        setcutoff!(sweeps, cutoff)
        setnoise!(sweeps, 1e-5, 1e-6, 1e-8, 0.0)
        observer = DMRGObserver(; energy_tol = energy_tol, minsweeps = minimum_sweeps)
        _, state = dmrg(
            constrained_H,
            state,
            sweeps;
            observer = observer,
            outputlevel = 0,
        )
        normalize!(state)

        violation = fibonacci_constraint_violation(state; pbc = model.pbc)
        violation <= 1e-6 || error(
            "DMRG state has Fibonacci-constraint violation $violation; " *
            "increase maxdim/sweeps or constraint_penalty",
        )
        energy = real(inner(prime(state), H, state))
        energy_history = real.(energies(observer))
        completed_sweeps = length(energy_history)
        final_energy_change = completed_sweeps >= 2 ?
                              abs(energy_history[end] - energy_history[end-1]) : Inf
        converged = final_energy_change < energy_tol
        converged || completed_sweeps == sweep_times || error(
            "DMRG stopped before reaching either the tolerance or sweep cap",
        )
        return state, energy, completed_sweeps, final_energy_change, converged
    end

    """Return the CFT ground state obtained from the sparse anyon Hamiltonian."""
    function initial_cft_exact(model, spec)
        model.pbc || error("The $(spec.title) initial state requires periodic boundaries")
        model.measure_operator == spec.initial_measure_operator || error(
            "The $(spec.title) initial state requires " *
            "measure_operator=:$(spec.initial_measure_operator)",
        )
        cached = get!(CFT_EXACT_GS_CACHE, (Symbol(spec.file_suffix), model.N)) do
            # `anyon_ham_sparse(model)` is the repository's sparse-Hamiltonian API.
            energies, states = Arpack.eigs(anyon_ham_sparse(model), which = :SR, nev = 1)
            (; state = vec(states[:, 1]), energy = real(energies[1]))
        end
        return copy(cached.state), cached.energy
    end

    """Return the cached constrained-DMRG CFT ground state."""
    function initial_cft_mps(
        model,
        spec;
        maxdim::Int,
        sweep_times::Int = cft_dmrg_sweeps(model.N),
        minimum_sweeps::Int = cft_dmrg_min_sweeps(spec, model.N),
        cutoff::Float64 = CFT_DMRG_CUTOFF,
        constraint_penalty::Float64 = CFT_CONSTRAINT_PENALTY,
        energy_tol::Float64 = spec.dmrg_energy_tol,
        dmrg_seed::Int = CFT_DMRG_SEED,
    )
        key = (
            Symbol(spec.file_suffix),
            model.N,
            maxdim,
            sweep_times,
            minimum_sweeps,
            cutoff,
            constraint_penalty,
            energy_tol,
            dmrg_seed,
        )
        cached = get!(CFT_MPS_GS_CACHE, key) do
            state, energy, completed_sweeps, final_energy_change, converged = cft_mps_gst(
                model;
                maxdim = maxdim,
                sweep_times = sweep_times,
                minimum_sweeps = minimum_sweeps,
                cutoff = cutoff,
                constraint_penalty = constraint_penalty,
                energy_tol = energy_tol,
                seed = dmrg_seed,
            )
            violation = fibonacci_constraint_violation(state; pbc = model.pbc)
            eelis = Float64.(anyon_eelis(model, state))
            mincut = spec.entanglement_mincut(model.N)
            c = fit_central_charge_pbc(eelis; mincut = mincut)
            if model.N >= spec.min_central_charge_check_size
                expected = spec.expected_central_charge
                atol = spec.central_charge_atol
                abs(c - expected) <= atol || error(
                    "$(spec.title) DMRG ground state has fitted central charge c=$c, " *
                    "outside |c - $expected| <= $atol; " *
                    "increase maxdim/sweeps or constraint_penalty",
                )
            end
            (;
                state,
                sites = siteinds(state),
                energy,
                violation,
                eelis,
                c,
                mincut,
                completed_sweeps,
                final_energy_change,
                converged,
            )
        end
        return (
            copy(cached.state),
            copy(cached.sites),
            cached.energy,
            cached.violation,
            cached.eelis,
            cached.c,
            cached.mincut,
            cached.completed_sweeps,
            cached.final_energy_change,
            cached.converged,
        )
    end

    """Fit the PBC CFT form `S(l) = c*x(l) + const`, where
    `x(l) = log[(L/pi)sin(pi*l/L)]/3`, retaining cuts
    `mincut:(L-mincut)`.  Here `mincut` is the first included cut."""
    function fit_central_charge_pbc(eelis::Vector{Float64}; mincut::Int = 1)
        L = length(eelis) + 1
        idx = collect(mincut:L-mincut)
        length(idx) >= 3 || error("At least three entropy cuts are required")
        xs = log.((L / pi) .* sin.(pi .* idx ./ L)) ./ 3
        ys = eelis[idx]
        xm, ym = sum(xs) / length(xs), sum(ys) / length(ys)
        slope = sum((xs .- xm) .* (ys .- ym)) / sum((xs .- xm) .^ 2)
        return slope
    end

    function cft_trajectory_dir(
        spec,
        backend::Symbol,
        L::Integer,
        tau_idx::Integer;
        chi::Union{Nothing,Integer} = nothing,
    )
        backend in (:exact, :mps) || error("backend must be :exact or :mps")
        backend == :mps && isnothing(chi) && error("chi is required for the MPS backend")
        source = backend == :exact ? "dense_eigensolver" : "dmrg"
        path = joinpath(
            spec.data_root,
            String(backend),
            source,
            "L$(L)",
            "gammaind$(tau_idx)",
        )
        return backend == :mps ? joinpath(path, "chi$(chi)") : path
    end

    function cft_default_periods(backend::Symbol, L::Integer, tau_idx::Integer)
        if backend == :exact
            D, _, _ = get_cfg_params_Born(tau_idx, L)
            return div(D, 2)
        elseif backend == :mps
            time_in_L, _, _ = get_mps_params_Born(tau_idx, L)
            return time_in_L * L
        end
        error("backend must be :exact or :mps")
    end

    """The initial-state Hamiltonian uses the CFT's own measurement operator
    (Antiferro for TCI, Ferro for Potts); the monitored dynamics is always
    Antiferro (`fib_model`)."""
    initial_cft_model(spec, L::Integer) = AnyonModel(
        FibonacciAnyon(),
        L;
        pbc = true,
        measure_operator = spec.initial_measure_operator,
    )

    function simulate_cft_exact(
        spec,
        L::Integer,
        tau_idx::Integer,
        seed::Integer;
        periods::Integer = cft_default_periods(:exact, L, tau_idx),
    )
        iseven(L) || error("L must be even, got $L")
        τ = τlis_ext[tau_idx]
        dynamics_model = fib_model(L)
        initial_model = initial_cft_model(spec, L)
        initial_state, initial_energy = initial_cft_exact(initial_model, spec)
        Y = topological_charge_operator(initial_model)
        initial_y_expectation =
            real(dot(initial_state, Y * initial_state)) /
            real(dot(initial_state, initial_state))
        config = MeasureConfig(
            τ = τ,
            mode = :Born,
            t₂ = periods,
            rng = MersenneTwister(seed),
            track_y_expectation = true,
        )
        outcome = bulk_evolution(dynamics_model, initial_state, config)
        return (
            samples = outcome.samples,
            free_energy = outcome.free_energys,
            halfchain_entropy = outcome.entanglement_entropys,
            final_entropy_profile = anyon_eelis(dynamics_model, outcome.state),
            initial_energy = initial_energy,
            initial_mpo_energy = initial_energy,
            initial_constraint_violation = 0.0,
            initial_bond_dimension = 0,
            initial_y_expectation = initial_y_expectation,
            final_y_expectation = Float64(outcome.y_expectation_values[end]),
        )
    end

    function simulate_cft_mps(
        spec,
        L::Integer,
        tau_idx::Integer,
        chi::Integer,
        seed::Integer;
        periods::Integer = cft_default_periods(:mps, L, tau_idx),
    )
        iseven(L) || error("L must be even, got $L")
        τ = τlis_ext[tau_idx]
        dynamics_model = fib_model(L)
        (
            initial_state,
            sites,
            initial_energy,
            violation,
            initial_eelis,
            initial_c,
            initial_c_mincut,
            dmrg_completed_sweeps,
            dmrg_final_energy_change,
            dmrg_converged,
        ) = initial_cft_mps(initial_cft_model(spec, L), spec; maxdim = chi)
        initial_bond_dimension = maxlinkdim(initial_state)
        Y_mpo = topological_charge_mpo(sites; pbc = true)
        initial_y_expectation =
            real(inner(prime(initial_state), Y_mpo, initial_state)) /
            real(inner(initial_state, initial_state))
        config = MeasureConfig(
            τ = τ,
            mode = :Born,
            t₂ = periods,
            rng = MersenneTwister(seed),
            cutoff = 1e-12,
            maxdim = chi,
            enforce_fibonacci_constraint = true,
            track_y_expectation = true,
        )
        outcome = bulk_evolution(dynamics_model, sites, initial_state, config)
        return (
            samples = outcome.samples,
            free_energy = outcome.free_energys,
            halfchain_entropy = outcome.entanglement_entropys,
            final_entropy_profile = anyon_eelis(dynamics_model, outcome.state),
            initial_energy = initial_energy,
            initial_mpo_energy = initial_energy,
            initial_constraint_violation = violation,
            initial_bond_dimension = initial_bond_dimension,
            initial_eelis = initial_eelis,
            initial_central_charge = initial_c,
            initial_central_charge_mincut = initial_c_mincut,
            dmrg_completed_sweeps = dmrg_completed_sweeps,
            dmrg_final_energy_change = dmrg_final_energy_change,
            dmrg_converged = dmrg_converged,
            initial_y_expectation = initial_y_expectation,
            final_y_expectation = Float64(outcome.y_expectation_values[end]),
        )
    end

    function save_cft_trajectory(
        spec,
        backend::Symbol,
        L::Integer,
        tau_idx::Integer,
        seed::Integer,
        result;
        chi::Union{Nothing,Integer} = nothing,
    )
        out_dir = cft_trajectory_dir(spec, backend, L, tau_idx; chi = chi)
        mkpath(out_dir)
        periods = length(result.halfchain_entropy)
        output_path =
            joinpath(out_dir, "periods$(periods)_trajectory_seed$(seed).jld2")
        data = Dict{String,Any}(
            "backend" => String(backend),
            "initial_state" => spec.initial_state_tag,
            "initial_state_method" =>
                backend == :exact ? "dense_eigensolver" : "dmrg",
            "L" => Int(L),
            "tau_idx" => Int(tau_idx),
            "tau" => τlis_ext[tau_idx],
            "gamma" => tanh(τlis_ext[tau_idx]),
            "chi" => isnothing(chi) ? 0 : Int(chi),
            "seed" => Int(seed),
            "periods" => periods,
            "dmrg_sweeps" => backend == :mps ? result.dmrg_completed_sweeps : 0,
            "dmrg_cutoff" => backend == :mps ? CFT_DMRG_CUTOFF : 0.0,
            "dmrg_seed" => backend == :mps ? CFT_DMRG_SEED : 0,
            "constraint_penalty" =>
                backend == :mps ? CFT_CONSTRAINT_PENALTY : 0.0,
            "initial_energy" => result.initial_energy,
            "initial_mpo_energy" => result.initial_mpo_energy,
            "initial_constraint_violation" => result.initial_constraint_violation,
            "initial_bond_dimension" => result.initial_bond_dimension,
            "initial_eelis" => get(result, :initial_eelis, Float64[]),
            "initial_central_charge" => get(result, :initial_central_charge, NaN),
            "initial_y_expectation" => result.initial_y_expectation,
            "final_y_expectation" => result.final_y_expectation,
            "sample" => result.samples,
            "sample_free_energy" => result.free_energy,
            "halfchain_EE_tlis" => result.halfchain_entropy,
            "final_EElis" => result.final_entropy_profile,
        )
        if spec.file_suffix == "potts"
            data["initial_state_cft"] = "potts"
            data["initial_hamiltonian_measure_operator"] = "Ferro"
            data["dynamics_measure_operator"] = "Antiferro"
            data["dmrg_max_sweeps"] = backend == :mps ? cft_dmrg_sweeps(L) : 0
            data["dmrg_min_sweeps"] =
                backend == :mps ? cft_dmrg_min_sweeps(spec, L) : 0
            data["dmrg_energy_tolerance"] =
                backend == :mps ? spec.dmrg_energy_tol : 0.0
            data["dmrg_final_energy_change"] =
                get(result, :dmrg_final_energy_change, NaN)
            data["dmrg_converged"] = get(result, :dmrg_converged, false)
            data["initial_central_charge_mincut"] =
                get(result, :initial_central_charge_mincut, 0)
        end
        save(output_path, data)
        return output_path
    end

    function run_cft_task(spec, task)
        backend, L, tau_idx, chi, seed, periods = task
        try
            result = if backend == :exact
                simulate_cft_exact(spec, L, tau_idx, seed; periods = periods)
            else
                simulate_cft_mps(spec, L, tau_idx, chi, seed; periods = periods)
            end
            output_path = save_cft_trajectory(
                spec,
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

    function collect_cft_trajectories(
        spec,
        backend::Symbol,
        L::Integer,
        tau_idx::Integer;
        chi::Union{Nothing,Integer} = nothing,
        periods::Integer = cft_default_periods(backend, L, tau_idx),
    )
        data_dir = cft_trajectory_dir(spec, backend, L, tau_idx; chi = chi)
        isdir(data_dir) ||
            error("$(spec.title) trajectory directory does not exist: $data_dir")
        files = sort(filter(
            file -> startswith(file, "periods$(periods)_trajectory_seed") &&
                    endswith(file, ".jld2"),
            readdir(data_dir),
        ))
        isempty(files) && error(
            "No $(spec.title) trajectory files with periods=$periods found in $data_dir",
        )

        samples_num = length(files)
        entropy_dynamics = zeros(Float64, samples_num, periods)
        final_profiles = zeros(Float64, samples_num, L - 1)
        free_energies = Vector{Vector{Float64}}(undef, samples_num)
        seeds = zeros(Int, samples_num)
        initial_energies = zeros(Float64, samples_num)
        initial_mpo_energies = zeros(Float64, samples_num)
        initial_violations = zeros(Float64, samples_num)
        initial_bond_dimensions = zeros(Int, samples_num)
        initial_central_charges = fill(NaN, samples_num)
        initial_central_charge_mincuts = zeros(Int, samples_num)
        initial_y_expectations = zeros(Float64, samples_num)
        final_y_expectations = zeros(Float64, samples_num)
        dmrg_completed_sweeps = zeros(Int, samples_num)
        dmrg_final_energy_changes = fill(NaN, samples_num)
        expected_method = backend == :exact ? "dense_eigensolver" : "dmrg"
        is_potts = spec.file_suffix == "potts"

        for (i, file) in enumerate(files)
            data = JLD2.load(joinpath(data_dir, file))
            data["initial_state"] == spec.initial_state_tag ||
                error("Unexpected initial state in $file")
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
                Float64(data["dmrg_cutoff"]) == CFT_DMRG_CUTOFF ||
                    error("Inconsistent DMRG cutoff in $file")
                Int(data["dmrg_seed"]) == CFT_DMRG_SEED ||
                    error("Inconsistent DMRG seed in $file")
                Float64(data["constraint_penalty"]) == CFT_CONSTRAINT_PENALTY ||
                    error("Inconsistent constraint penalty in $file")
            end
            if is_potts
                data["initial_state_cft"] == "potts" ||
                    error("Unexpected initial-state CFT in $file")
                data["initial_hamiltonian_measure_operator"] == "Ferro" ||
                    error("Unexpected initial Hamiltonian in $file")
                data["dynamics_measure_operator"] == "Antiferro" ||
                    error("Unexpected dynamics operator in $file")
            elseif backend == :mps
                Int(data["dmrg_sweeps"]) == cft_dmrg_sweeps(L) ||
                    error("Inconsistent DMRG sweep count in $file")
            end
            if is_potts && backend == :mps
                completed_sweeps = Int(data["dmrg_sweeps"])
                max_sweeps = Int(data["dmrg_max_sweeps"])
                min_sweeps = Int(data["dmrg_min_sweeps"])
                max_sweeps == cft_dmrg_sweeps(L) ||
                    error("Inconsistent DMRG maximum sweep count in $file")
                min_sweeps == cft_dmrg_min_sweeps(spec, L) ||
                    error("Inconsistent DMRG minimum sweep count in $file")
                min_sweeps <= completed_sweeps <= max_sweeps ||
                    error("Invalid completed DMRG sweep count in $file")
                Float64(data["dmrg_energy_tolerance"]) == spec.dmrg_energy_tol ||
                    error("Inconsistent DMRG energy tolerance in $file")
                Bool(data["dmrg_converged"]) || completed_sweeps == max_sweeps ||
                    error("DMRG neither converged nor reached its sweep cap in $file")
                Int(data["initial_central_charge_mincut"]) ==
                    spec.entanglement_mincut(L) ||
                    error("Inconsistent central-charge mincut in $file")
            end
            entropy_dynamics[i, :] = Float64.(data["halfchain_EE_tlis"])
            final_profiles[i, :] = Float64.(data["final_EElis"])
            free_energies[i] = Float64.(data["sample_free_energy"])
            seeds[i] = Int(data["seed"])
            initial_energies[i] = Float64(data["initial_energy"])
            initial_mpo_energies[i] = Float64(data["initial_mpo_energy"])
            initial_violations[i] = Float64(data["initial_constraint_violation"])
            initial_bond_dimensions[i] = Int(data["initial_bond_dimension"])
            initial_central_charges[i] = Float64(data["initial_central_charge"])
            initial_y_expectations[i] = Float64(data["initial_y_expectation"])
            final_y_expectations[i] = Float64(data["final_y_expectation"])
            if is_potts
                initial_central_charge_mincuts[i] =
                    Int(data["initial_central_charge_mincut"])
                dmrg_completed_sweeps[i] = Int(data["dmrg_sweeps"])
                dmrg_final_energy_changes[i] = Float64(data["dmrg_final_energy_change"])
            end
        end
        length(unique(seeds)) == samples_num || error("Duplicate trajectory seeds found")
        all(length(fe) == length(first(free_energies)) for fe in free_energies) ||
            error("Inconsistent free-energy trajectory lengths")
        maximum(initial_energies) - minimum(initial_energies) <= 1e-6 ||
            error("Trajectories used inconsistent $(spec.title) ground states")
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

        # Use the shared backend-specific average window from config.jl.  It
        # already fixes the parity of the Fibonacci layers and stops before
        # the terminal half-strength layer generated by the default
        # `enable_τ_eff=true`.
        total_layers = size(free_energy_matrix, 1)
        average_window = if backend == :exact
            get_cfg_params_Born(tau_idx, L)[3]
        elseif backend == :mps
            get_mps_params_Born(tau_idx, L)[3]
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
        bulk_FE_stderr = std(trajectory_FE; corrected = false) / sqrt(samples_num)

        output_path = joinpath(
            data_dir,
            "EE_FEdynamics_L$(L)_gamma$(tau_idx)_periods$(periods)_$(spec.file_suffix).jld2",
        )
        data = Dict{String,Any}(
            "backend" => String(backend),
            "initial_state" => spec.initial_state_tag,
            "initial_state_method" => expected_method,
            "L" => Int(L),
            "tau_idx" => Int(tau_idx),
            "tau" => τlis_ext[tau_idx],
            "gamma" => tanh(τlis_ext[tau_idx]),
            "chi" => isnothing(chi) ? 0 : Int(chi),
            "periods" => periods,
            "samples_num" => samples_num,
            "initial_energy" => mean(initial_energies),
            "initial_mpo_energy" => mean(initial_mpo_energies),
            "maximum_initial_constraint_violation" => maximum(initial_violations),
            "maximum_initial_bond_dimension" => maximum(initial_bond_dimensions),
            "initial_central_charge" => first(initial_central_charges),
            "initial_y_expectation" => mean(initial_y_expectations),
            "final_y_expectation" => mean(final_y_expectations),
            "final_y_expectation_stderr" =>
                std(final_y_expectations; corrected = false) / sqrt(samples_num),
            "average_EE_tlis" => average_EE_tlis,
            "stderr_EE_tlis" => stderr_EE_tlis,
            "bulk_meanEElis" => bulk_meanEElis,
            "ensemble_stderr_EElis" => ensemble_stderr_EElis,
            "time_average_free_energy" => trajectory_FE,
            "bulk_FE" => bulk_FE,
            "bulk_FE_stderr" => bulk_FE_stderr,
            "time_FEstderr" => time_FEstderr,
            "time_FElis" => time_FElis,
            "ensemble_seed" => seeds,
        )
        if is_potts
            data["initial_state_cft"] = "potts"
            data["initial_hamiltonian_measure_operator"] = "Ferro"
            data["dynamics_measure_operator"] = "Antiferro"
            data["initial_central_charge_mincut"] =
                first(initial_central_charge_mincuts)
            data["dmrg_completed_sweeps"] = first(dmrg_completed_sweeps)
            data["dmrg_final_energy_change"] = first(dmrg_final_energy_changes)
        end
        save(output_path, data)
        return output_path
    end
end

function print_cft_usage()
    script = "exm/Bulk_measure/monitored_dynamics_cft.jl"
    println("Usage:")
    println("  julia --project=. -p N $script CFT exact L tau_idx seed_start seed_end [periods]")
    println("  julia --project=. -p N $script CFT mps L tau_idx chi seed_start seed_end [periods]")
    println("  julia --project=. $script CFT collect_exact L tau_idx [periods]")
    println("  julia --project=. $script CFT collect_mps L tau_idx chi [periods]")
    println("  CFT is either tci or potts")
end

if isempty(ARGS)
    print_cft_usage()
else
    length(ARGS) >= 2 || (print_cft_usage(); error("Not enough arguments"))
    spec = cft_spec(Symbol(lowercase(ARGS[1])))
    action = Symbol(lowercase(ARGS[2]))
    if action in (:exact, :mps)
        minimum_args = action == :exact ? 6 : 7
        length(ARGS) >= minimum_args || error("Not enough arguments for $action")
        L = parse(Int, ARGS[3])
        tau_idx = parse(Int, ARGS[4])
        chi = action == :mps ? parse(Int, ARGS[5]) : 0
        seed_offset = action == :mps ? 5 : 4
        seed_start = parse(Int, ARGS[seed_offset + 1])
        seed_end = parse(Int, ARGS[seed_offset + 2])
        seed_start <= seed_end || error("seed_start must not exceed seed_end")
        default_periods = cft_default_periods(action, L, tau_idx)
        periods = length(ARGS) >= seed_offset + 3 ?
                  parse(Int, ARGS[seed_offset + 3]) : default_periods
        periods >= 1 || error("periods must be positive")

        tasks = [(action, L, tau_idx, chi, seed, periods) for seed in seed_start:seed_end]
        println("=== $(spec.title) monitored dynamics ($(action)) ===")
        println("L=$L, tau_idx=$tau_idx, periods=$periods, trajectories=$(length(tasks))")
        action == :mps && println(
            "chi=$chi, DMRG sweeps=$(cft_dmrg_min_sweeps(spec, L)):$(cft_dmrg_sweeps(L)), " *
            "energy tolerance=$(spec.dmrg_energy_tol), " *
            "central-charge mincut=$(spec.entanglement_mincut(L)), " *
            "constraint penalty=$CFT_CONSTRAINT_PENALTY",
        )
        println("workers=$(nworkers())")
        results = pmap(task -> run_cft_task(spec, task), tasks; batch_size = 1)
        failed = filter(result -> result[5] != :success, results)
        println("Successes: $(length(results) - length(failed))")
        println("Failures: $(length(failed))")
        for result in failed
            println("Failed seed $(result[4]): $(result[6])")
        end
    elseif action == :collect_exact
        length(ARGS) in (4, 5) || error("collect_exact requires L tau_idx [periods]")
        L = parse(Int, ARGS[3])
        tau_idx = parse(Int, ARGS[4])
        periods = length(ARGS) == 5 ?
                  parse(Int, ARGS[5]) : cft_default_periods(:exact, L, tau_idx)
        println("Saved $(collect_cft_trajectories(spec, :exact, L, tau_idx; periods = periods))")
    elseif action == :collect_mps
        length(ARGS) in (5, 6) || error("collect_mps requires L tau_idx chi [periods]")
        L = parse(Int, ARGS[3])
        tau_idx = parse(Int, ARGS[4])
        chi = parse(Int, ARGS[5])
        periods = length(ARGS) == 6 ?
                  parse(Int, ARGS[6]) : cft_default_periods(:mps, L, tau_idx)
        println(
            "Saved $(collect_cft_trajectories(spec, :mps, L, tau_idx; chi = chi, periods = periods))",
        )
    else
        print_cft_usage()
        error("Unknown action: $(ARGS[2])")
    end
end
