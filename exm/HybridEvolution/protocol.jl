using FibonacciChain
using LinearAlgebra
using Random
using Statistics
using Serialization

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

"""Write an ensemble at one (L,p); NaN sharpening time means right-censored."""
function save_ensemble(directory, results; epsilon = 0.05, fraction = 0.9)
    first_result = first(results)
    n = length(results)
    open(joinpath(directory, "trajectories.csv"), "w") do io
        println(io, "L,p,trajectory_seed,time,S_half,Y_expectation,is_sharp")
        for r in results, i in eachindex(r.times)
            println(io, join((r.L, r.p, r.seed, r.times[i], r.entropy[i],
                r.y_expectation[i], sharp(r.y_expectation[i], epsilon)), ','))
        end
    end
    fractions = Float64[]
    open(joinpath(directory, "summary.csv"), "w") do io
        println(io, "L,p,time,n,S_mean,S_sem,S_density,Y_mean,Y_sem,sharp_fraction,sharp_sem")
        for i in eachindex(first_result.times)
            s = [r.entropy[i] for r in results]
            y = [r.y_expectation[i] for r in results]
            indicators = [sharp(v, epsilon) for v in y]
            push!(fractions, mean(indicators))
            println(io, join((first_result.L, first_result.p, first_result.times[i],
                n, mean(s), standard_error(s), mean(s) / first_result.L,
                mean(y), standard_error(y), last(fractions),
                standard_error(indicators)), ','))
        end
    end
    index = findfirst(>=(fraction), fractions)
    open(joinpath(directory, "sharpening.csv"), "w") do io
        println(io, "L,p,epsilon_Y,target_fraction,t_sharp,censored,last_time")
        println(io, join((first_result.L, first_result.p, epsilon, fraction,
            isnothing(index) ? NaN : first_result.times[index], isnothing(index),
            last(first_result.times)), ','))
    end
    if first_result.schedule !== nothing
        for r in results
            serialize(joinpath(directory, "schedule_seed$(r.seed).jls"),
                (; L = r.L, p = r.p, seed = r.seed, schedule = r.schedule))
        end
    end
    if hasproperty(first_result, :final_bond_dimension)
        open(joinpath(directory, "mps_diagnostics.csv"), "w") do io
            println(io, "L,p,trajectory_seed,final_bond_dimension")
            for r in results
                println(io, join((r.L, r.p, r.seed, r.final_bond_dimension), ','))
            end
        end
    end
    return directory
end
