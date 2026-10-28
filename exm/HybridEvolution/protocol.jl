using FibonacciChain
using LinearAlgebra
using Random
using Statistics
using JLD2

include("protocol_mps.jl")

"""Generate a coherent-state trajectory using the package hybrid evolution."""
function samples_generate(L::Int, p::Float64, periods::Int, seed::Int;
    stride::Int = 1, save_schedule::Bool = true,
)
    L >= 4 && iseven(L) || throw(ArgumentError("L must be even and >= 4"))
    periods >= 1 && stride >= 1 && seed >= 0 ||
        throw(ArgumentError("periods/stride must be positive; seed nonnegative"))
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    basis = anyon_basis(model)
    state = zeros(ComplexF64, length(basis))
    state[1] = 1
    iszero(basis[1].buf) || error("Expected the all-zero initial fusion path")
    initial_y = Fsymmetry_coef(model, basis[1], basis[1])
    phi = (1 + sqrt(5.0)) / 2
    initial_weight = (initial_y + inv(phi)) / (phi + inv(phi))
    0 < initial_weight < 1 || error("Initial state must have both sector components")

    config = HybridConfig(
        τ = Inf,
        t₂ = periods,
        mode = :Born,
        rng = MersenneTwister(seed),
        enable_τ_eff = false,
        track_y_expectation = true,
        p = p,
        random_angles = true,
    )
    outcome = bulk_evolution(model, state, config)
    times = sort!(unique!([0; collect(stride:stride:periods); periods]))
    entropy = [ee(anyon_rdm(model, collect(1:(L ÷ 2)), state));
        Float64.(outcome.entanglement_entropys)][times .+ 1]
    y_expectation = [initial_y; Float64.(outcome.y_expectation_values)][times .+ 1]
    free_energys, free_energy = free_energy_dynamics(outcome.free_energys, times, periods)
    reference_entropy_final = reference_entropy_from_y(last(y_expectation))
    reference_entropy = reference_entropy_from_y.(y_expectation)

    measurement_count = count(outcome.schedule.measurement_mask)
    schedule = save_schedule ? outcome.schedule : nothing
    final_observables = final_entanglement(model, outcome.state)
    return (; L, p, periods, seed, times, entropy, y_expectation,
        initial_weight, measurement_count, reference_entropy, reference_entropy_final,
        schedule, free_energys, free_energy, final_observables...)
end

"""Orthonormal trivial-sector purification branches `|v0⟩, |v1⟩`, each with
norm `1/√2` so the joint ancilla-system state is normalized."""
function purification_initial_states(model, sector = :trivial)
    P, _, dimension, _ = FibonacciChain._hybrid_sector_projector(model, sector)
    dimension >= 2 || error("Purification requires two states in the selected sector")
    factor = qr(P[:, 1:2])
    minimum(abs.(diag(factor.R))) > 1e-12 ||
        error("Initial reference branches are linearly dependent in sector $sector")
    return ComplexF64.(Matrix(factor.Q)[:, 1:2]) / sqrt(2)
end

"""Generate same-sector purification entropy dynamics with a record Born-sampled
from the joint ancilla-system state.

The purification ancilla is equally entangled with two orthogonal states in
the trivial topological sector. First evolve both branches through a common
unitary-only circuit for L^2 complete periods. Monitored dynamics then start
at t=0, with the measurement record drawn from the joint state's Born
distribution. The same RNG is continued across scrambling and monitoring.
The decay of the ensemble-averaged `purification_entropy` gives the
purification time. When `save_schedule=true`, return both `scramble_schedule`
and the monitored `schedule` for full replay.
"""
function purification_samples_generate(L::Int, p::Float64, periods::Int, seed::Int;
    stride::Int = 1, save_schedule::Bool = true,
)
    L >= 4 && iseven(L) || throw(ArgumentError("L must be even and >= 4"))
    periods >= 1 && stride >= 1 && seed >= 0 ||
        throw(ArgumentError("periods/stride must be positive; seed nonnegative"))
    0 <= p <= 1 || throw(ArgumentError("p must satisfy 0 ≤ p ≤ 1, got $p"))
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    purification_sector = :trivial
    states = purification_initial_states(model, purification_sector)
    rng = MersenneTwister(seed)
    scramble_periods = L^2
    scramble_config = HybridConfig(
        τ = Inf,
        t₂ = scramble_periods,
        mode = :Born,
        rng = rng,
        enable_τ_eff = false,
        p = 0.0,
        random_angles = true,
    )
    scrambled = hybrid_purification_evolution(
        model, states[:, 1], states[:, 2], scramble_config,
    )
    config = HybridConfig(
        τ = Inf,
        t₂ = periods,
        mode = :Born,
        rng = rng,
        enable_τ_eff = false,
        p = p,
        random_angles = true,
    )
    purified = hybrid_purification_evolution(
        model, scrambled.state0, scrambled.state1, config,
    )
    times = sort!(unique!([0; collect(stride:stride:periods); periods]))
    purification_entropy = [FibonacciChain.purification_entropy(scrambled.state0, scrambled.state1);
        Float64.(purified.purification_entropys)][times .+ 1]
    free_energys, free_energy = free_energy_dynamics(purified.free_energys, times, periods)
    measurement_count = count(purified.schedule.measurement_mask)
    schedule = save_schedule ? purified.schedule : nothing
    scramble_schedule = save_schedule ? scrambled.schedule : nothing
    return (; L, p, periods, seed, times, purification_entropy, purification_sector,
        scramble_periods, measurement_count, schedule, scramble_schedule,
        free_energys, free_energy)
end

function process_task(task)
    L, p, periods, seed, stride, save_schedule = task
    return samples_generate(L, p, periods, seed; stride, save_schedule)
end

function process_task_and_save(job)
    task, directory = job
    result = process_task(task)
    save_trajectory(directory, result)
    GC.gc()
    return result
end

function sharp(y, epsilon)
    phi = (1 + sqrt(5.0)) / 2
    return min(abs(y - phi), abs(y + inv(phi))) < epsilon
end

"""Reference-qubit entropy for orthogonal topological-charge sectors."""
function reference_entropy_from_y(y::Real)
    phi = (1 + sqrt(5.0)) / 2
    weight = (Float64(y) + inv(phi)) / (phi + inv(phi))
    isfinite(weight) && -1e-5 <= weight <= 1 + 1e-5 ||
        throw(ArgumentError("Unphysical Y expectation: $y"))
    weight = clamp(weight, 0.0, 1.0)
    return (weight == 0 || weight == 1) ? 0.0 :
        -weight * log(weight) - (1 - weight) * log1p(-weight)
end

standard_error(x) = length(x) > 1 ? std(x) / sqrt(length(x)) : NaN

"""Layer record surprisal and cumulative F(t)=-log P(record through t).

Times count complete even-then-odd periods; gate-choice and angle sampling
probabilities are excluded, exactly as in `bulk_evolution.free_energys`.
"""
function free_energy_dynamics(layer_values, times, periods)
    length(layer_values) == 2periods || error("Expected two free-energy layers per period")
    values = Float64.(layer_values)
    all(isfinite, values) || error("Nonfinite record free energy")
    cumulative = [0.0; cumsum(values)]
    return values, cumulative[2 .* times .+ 1]
end

const FREE_ENERGY_DEFINITION = "free_energys stores per-layer -log Born outcome probabilities; free_energy(t) sums the first 2t layers, with F(0)=0; gate-choice and angle sampling probabilities are excluded"

"""Mean and SEM time series of an observable matrix restricted to `mask` rows."""
function sector_stats(matrix, mask)
    subset = matrix[mask, :]
    return vec(mean(subset; dims = 1)),
        [standard_error(column) for column in eachcol(subset)]
end

"""Atomically save one completed trajectory, including its replay schedule when requested."""
function save_trajectory(directory, result)
    mkpath(directory)
    path = joinpath(directory, "trajectory_seed$(result.seed).jld2")
    ispath(path) && error("Trajectory already exists: $path")
    temporary = tempname(directory)
    try
        jldopen(temporary, "w") do file
            file["L"] = result.L
            file["p"] = result.p
            file["periods"] = result.periods
            file["seed"] = result.seed
            file["time"] = result.times
            file["S_half"] = result.entropy
            file["free_energys"] = result.free_energys
            file["free_energy"] = result.free_energy
            file["free_energy_definition"] = FREE_ENERGY_DEFINITION
            file["Y_expectation"] = result.y_expectation
            file["initial_weight"] = result.initial_weight
            file["measurement_count"] = result.measurement_count
            file["S_R"] = reference_entropy_from_y.(result.y_expectation)
            file["reference_entropy_final"] = result.reference_entropy_final
            if hasproperty(result, :purification_entropy)
                file["S_P"] = result.purification_entropy
                file["purification_sector"] = result.purification_sector
            end
            file["subsystem_sizes"] = result.subsystem_sizes
            file["S_subsystem_final"] = result.S_subsystem_final
            file["I3_partition"] = result.I3_partition
            file["I3_entropy_labels"] = ["A", "B", "C", "AB", "AC", "BC", "ABC"]
            file["I3_entropies"] = result.I3_entropies
            file["I3_final"] = result.I3_final
            if result.schedule !== nothing
                file["measurement_mask"] = result.schedule.measurement_mask
                file["outcomes"] = result.schedule.outcomes
                file["unitary_angles"] = result.schedule.unitary_angles
            end
            if hasproperty(result, :final_bond_dimension)
                file["final_bond_dimension"] = result.final_bond_dimension
                file["entropy_swap_cutoff"] = result.entropy_swap_cutoff
            end
        end
        mv(temporary, path)
    finally
        isfile(temporary) && rm(temporary)
    end
    return path
end

"""Save JLD2 datasets; observable matrices have axes (trajectory, time)."""
function save_ensemble(directory, results; epsilon = 0.05, fraction = 0.9)
    first_result = first(results)
    n = length(results)
    L, p, time = first_result.L, first_result.p, first_result.times
    all(r -> r.L == L && r.p == p && r.times == time, results) ||
        error("Ensemble parameters and observation times must agree")
    subsystem_sizes, I3_partition = first_result.subsystem_sizes, first_result.I3_partition
    all(r -> r.subsystem_sizes == subsystem_sizes && r.I3_partition == I3_partition,
        results) || error("Ensemble final-state partitions must agree")
    seeds = [r.seed for r in results]
    length(unique(seeds)) == n || error("Duplicate trajectory seeds")
    S_half = reduce(vcat, [permutedims(r.entropy) for r in results])
    Y_expectation = reduce(vcat, [permutedims(r.y_expectation) for r in results])
    S_R = reference_entropy_from_y.(Y_expectation)
    free_energys = reduce(vcat, [permutedims(r.free_energys) for r in results])
    free_energy = reduce(vcat, [permutedims(r.free_energy) for r in results])
    size(free_energys) == (n, 2first_result.periods) || error("Free-energy layer grid mismatch")
    size(free_energy) == size(S_half) || error("Free-energy observation grid mismatch")
    has_purification = hasproperty(first_result, :purification_entropy)
    all(r -> hasproperty(r, :purification_entropy) == has_purification, results) ||
        error("Ensemble purification observables must agree")
    if has_purification
        all(r -> r.purification_sector == first_result.purification_sector, results) ||
            error("Ensemble purification sectors must agree")
    end
    purification_data = has_purification ?
        (; S_P = reduce(vcat, [permutedims(r.purification_entropy) for r in results]),
            purification_sector = first_result.purification_sector) : (;)
    purification_stats = has_purification ?
        (; S_P_mean = vec(mean(purification_data.S_P; dims = 1)),
            S_P_sem = [standard_error(column) for column in eachcol(purification_data.S_P)],
            purification_sector = first_result.purification_sector) : (;)
    reference_entropy_final = [r.reference_entropy_final for r in results]
    S_subsystem_final = reduce(vcat, [permutedims(r.S_subsystem_final) for r in results])
    I3_entropies = reduce(vcat, [permutedims(r.I3_entropies) for r in results])
    I3_final = [r.I3_final for r in results]
    I3_entropy_labels = ["A", "B", "C", "AB", "AC", "BC", "ABC"]
    is_sharp = sharp.(Y_expectation, epsilon)
    jldsave(joinpath(directory, "trajectories.jld2");
        L, p, time, trajectory_seed = seeds, S_half, Y_expectation, is_sharp,
        epsilon_Y = epsilon, initial_weight = [r.initial_weight for r in results],
        measurement_count = [r.measurement_count for r in results],
        S_R, reference_entropy_final, subsystem_sizes, S_subsystem_final,
        I3_partition, I3_entropy_labels, I3_entropies, I3_final,
        free_energys, free_energy, free_energy_definition = FREE_ENERGY_DEFINITION,
        purification_data...)

    S_mean = vec(mean(S_half; dims = 1))
    fractions = vec(mean(is_sharp; dims = 1))
    phi = (1 + sqrt(5.0)) / 2
    sector1 = Y_expectation[:, end] .> (phi - inv(phi)) / 2
    S_mean_sector1, S_sem_sector1 = sector_stats(S_half, sector1)
    S_mean_sectortau, S_sem_sectortau = sector_stats(S_half, .!sector1)
    Y_mean_sector1, Y_sem_sector1 = sector_stats(Y_expectation, sector1)
    Y_mean_sectortau, Y_sem_sectortau = sector_stats(Y_expectation, .!sector1)
    jldsave(joinpath(directory, "summary.jld2");
        L, p, time, n, S_mean,
        free_energy_mean = vec(mean(free_energy; dims = 1)),
        free_energy_sem = [standard_error(column) for column in eachcol(free_energy)],
        free_energys_mean = vec(mean(free_energys; dims = 1)),
        free_energys_sem = [standard_error(column) for column in eachcol(free_energys)],
        free_energy_definition = FREE_ENERGY_DEFINITION,
        S_sem = [standard_error(column) for column in eachcol(S_half)],
        S_density = S_mean / L,
        Y_mean = vec(mean(Y_expectation; dims = 1)),
        Y_sem = [standard_error(column) for column in eachcol(Y_expectation)],
        S_R_mean = vec(mean(S_R; dims = 1)),
        S_R_sem = [standard_error(column) for column in eachcol(S_R)],
        reference_entropy_mean = mean(reference_entropy_final),
        reference_entropy_sem = standard_error(reference_entropy_final),
        subsystem_sizes,
        S_subsystem_final_mean = vec(mean(S_subsystem_final; dims = 1)),
        S_subsystem_final_sem = [standard_error(column) for column in eachcol(S_subsystem_final)],
        I3_partition, I3_entropy_labels,
        I3_entropies_mean = vec(mean(I3_entropies; dims = 1)),
        I3_entropies_sem = [standard_error(column) for column in eachcol(I3_entropies)],
        I3_final_mean = mean(I3_final), I3_final_sem = standard_error(I3_final),
        sharp_fraction = fractions,
        sharp_sem = [standard_error(column) for column in eachcol(is_sharp)],
        epsilon_Y = epsilon,
        sector_rule = "sector 1 if final-time m_Y > (phi - 1/phi)/2 = 0.5, else sector tau",
        n_sector1 = count(sector1), n_sectortau = count(.!sector1),
        S_mean_sector1, S_sem_sector1, S_mean_sectortau, S_sem_sectortau,
        Y_mean_sector1, Y_sem_sector1, Y_mean_sectortau, Y_sem_sectortau,
        purification_stats...)
    index = findfirst(>=(fraction), fractions)
    jldsave(joinpath(directory, "sharpening.jld2");
        L, p, epsilon_Y = epsilon, target_fraction = fraction,
        t_sharp = isnothing(index) ? NaN : time[index],
        censored = isnothing(index), last_time = last(time))
    if hasproperty(first_result, :final_bond_dimension)
        jldsave(joinpath(directory, "mps_diagnostics.jld2");
            L, p, trajectory_seed = seeds,
            final_bond_dimension = [r.final_bond_dimension for r in results],
            entropy_swap_cutoff = [r.entropy_swap_cutoff for r in results])
    end
    return directory
end

"""Collect only per-point averages into one portable JLD2 file."""
function save_averaged(output, sizes, rates, trajectories_per_point;
    filename = "averaged.jld2",
)
    path = joinpath(output, filename)
    ispath(path) && error("Averaged file already exists: $path")
    temporary = tempname(output)
    try
        jldopen(temporary, "w") do file
            file["sizes"] = collect(sizes)
            file["rates"] = collect(rates)
            file["trajectories_per_point"] = trajectories_per_point
            config_path = joinpath(output, "config.toml")
            if isfile(config_path)
                file["config_toml"] = read(config_path, String)
            end
            for L in sizes, p in rates
                directory = joinpath(output, "L$(L)_p$(p)")
                summary = load(joinpath(directory, "summary.jld2"))
                sharpening = load(joinpath(directory, "sharpening.jld2"))
                summary["L"] == L && summary["p"] == p &&
                    summary["n"] == trajectories_per_point ||
                    error("Incomplete or mismatched ensemble in $directory")
                haskey(summary, "reference_entropy_mean") ||
                    error("Missing reference entropy in $directory")
                prefix = "L$(L)/p$(p)"
                for (key, value) in summary
                    file["$prefix/summary/$key"] = value
                end
                for (key, value) in sharpening
                    file["$prefix/sharpening/$key"] = value
                end
            end
        end
        mv(temporary, path)
    finally
        isfile(temporary) && rm(temporary)
    end
    return path
end
