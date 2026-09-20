using FibonacciChain
using LinearAlgebra
using Random
using Statistics
using JLD2

include("protocol_mps.jl")

"""Generate a coherent-state trajectory using the package hybrid evolution."""
function samples_generate(L::Int, p::Float64, periods::Int, seed::Int;
    stride::Int = 1, save_schedule::Bool = false,
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
    measurement_count = count(outcome.schedule.measurement_mask)
    schedule = save_schedule ? outcome.schedule : nothing
    return (; L, p, periods, seed, times, entropy, y_expectation,
        initial_weight, measurement_count, schedule)
end

function process_task(task)
    L, p, periods, seed, stride, save_schedule = task
    return samples_generate(L, p, periods, seed; stride, save_schedule)
end

function sharp(y, epsilon)
    phi = (1 + sqrt(5.0)) / 2
    return min(abs(y - phi), abs(y + inv(phi))) < epsilon
end

standard_error(x) = length(x) > 1 ? std(x) / sqrt(length(x)) : NaN

"""Save JLD2 datasets; observable matrices have axes (trajectory, time)."""
function save_ensemble(directory, results; epsilon = 0.05, fraction = 0.9)
    first_result = first(results)
    n = length(results)
    L, p, time = first_result.L, first_result.p, first_result.times
    all(r -> r.L == L && r.p == p && r.times == time, results) ||
        error("Ensemble parameters and observation times must agree")
    seeds = [r.seed for r in results]
    S_half = reduce(vcat, [permutedims(r.entropy) for r in results])
    Y_expectation = reduce(vcat, [permutedims(r.y_expectation) for r in results])
    is_sharp = sharp.(Y_expectation, epsilon)
    jldsave(joinpath(directory, "trajectories.jld2");
        L, p, time, trajectory_seed = seeds, S_half, Y_expectation, is_sharp,
        epsilon_Y = epsilon, initial_weight = [r.initial_weight for r in results],
        measurement_count = [r.measurement_count for r in results])

    S_mean = vec(mean(S_half; dims = 1))
    fractions = vec(mean(is_sharp; dims = 1))
    jldsave(joinpath(directory, "summary.jld2");
        L, p, time, n, S_mean,
        S_sem = [standard_error(column) for column in eachcol(S_half)],
        S_density = S_mean / L,
        Y_mean = vec(mean(Y_expectation; dims = 1)),
        Y_sem = [standard_error(column) for column in eachcol(Y_expectation)],
        sharp_fraction = fractions,
        sharp_sem = [standard_error(column) for column in eachcol(is_sharp)],
        epsilon_Y = epsilon)
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
