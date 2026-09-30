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
    reference_entropy_final = reference_entropy_from_y(last(y_expectation))
    measurement_count = count(outcome.schedule.measurement_mask)
    schedule = save_schedule ? outcome.schedule : nothing
    return (; L, p, periods, seed, times, entropy, y_expectation,
        initial_weight, measurement_count, reference_entropy_final, schedule)
end

function process_task(task)
    L, p, periods, seed, stride, save_schedule = task
    return samples_generate(L, p, periods, seed; stride, save_schedule)
end

function process_task_and_save(job)
    task, directory = job
    result = process_task(task)
    save_trajectory(directory, result)
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
            file["Y_expectation"] = result.y_expectation
            file["initial_weight"] = result.initial_weight
            file["measurement_count"] = result.measurement_count
            file["reference_entropy_final"] = result.reference_entropy_final
            if result.schedule !== nothing
                file["measurement_mask"] = result.schedule.measurement_mask
                file["outcomes"] = result.schedule.outcomes
                file["unitary_angles"] = result.schedule.unitary_angles
            end
            if hasproperty(result, :final_bond_dimension)
                file["final_bond_dimension"] = result.final_bond_dimension
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
    seeds = [r.seed for r in results]
    length(unique(seeds)) == n || error("Duplicate trajectory seeds")
    S_half = reduce(vcat, [permutedims(r.entropy) for r in results])
    Y_expectation = reduce(vcat, [permutedims(r.y_expectation) for r in results])
    reference_entropy_final = [r.reference_entropy_final for r in results]
    is_sharp = sharp.(Y_expectation, epsilon)
    jldsave(joinpath(directory, "trajectories.jld2");
        L, p, time, trajectory_seed = seeds, S_half, Y_expectation, is_sharp,
        epsilon_Y = epsilon, initial_weight = [r.initial_weight for r in results],
        measurement_count = [r.measurement_count for r in results],
        reference_entropy_final)

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
        S_sem = [standard_error(column) for column in eachcol(S_half)],
        S_density = S_mean / L,
        Y_mean = vec(mean(Y_expectation; dims = 1)),
        Y_sem = [standard_error(column) for column in eachcol(Y_expectation)],
        reference_entropy_mean = mean(reference_entropy_final),
        reference_entropy_sem = standard_error(reference_entropy_final),
        sharp_fraction = fractions,
        sharp_sem = [standard_error(column) for column in eachcol(is_sharp)],
        epsilon_Y = epsilon,
        sector_rule = "sector 1 if final-time m_Y > (phi - 1/phi)/2 = 0.5, else sector tau",
        n_sector1 = count(sector1), n_sectortau = count(.!sector1),
        S_mean_sector1, S_sem_sector1, S_mean_sectortau, S_sem_sectortau,
        Y_mean_sector1, Y_sem_sector1, Y_mean_sectortau, Y_sem_sectortau)
    index = findfirst(>=(fraction), fractions)
    jldsave(joinpath(directory, "sharpening.jld2");
        L, p, epsilon_Y = epsilon, target_fraction = fraction,
        t_sharp = isnothing(index) ? NaN : time[index],
        censored = isnothing(index), last_time = last(time))
    for r in results
        if r.schedule !== nothing
            jldsave(joinpath(directory, "schedule_seed$(r.seed).jld2");
                L = r.L, p = r.p, seed = r.seed,
                measurement_mask = r.schedule.measurement_mask,
                outcomes = r.schedule.outcomes,
                unitary_angles = r.schedule.unitary_angles)
        end
    end
    if hasproperty(first_result, :final_bond_dimension)
        jldsave(joinpath(directory, "mps_diagnostics.jld2");
            L, p, trajectory_seed = seeds,
            final_bond_dimension = [r.final_bond_dimension for r in results])
    end
    return directory
end

"""Collect only per-point averages into one portable JLD2 file."""
function save_averaged(output, sizes, rates, trajectories_per_point)
    path = joinpath(output, "averaged.jld2")
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
