using FibonacciChain
using JLD
using Random
using Statistics
# include("../FitEntEntScal.jl")

binary_distribution(p, rng) = rand(rng) < p ? 1 : 0
function ps_prob_evolution(L::Int64, τ::Float64,  D::Int, prob::Float64; seed::Int=100, temp::Bool=false)
    # Generate a state based on the given sample and initial state
    # τ is the temperature parameter, initial_state is the initial state vector
    # sample is a matrix of binary values representing the state configuration
    gate_num = round(Int, D*L / 2)
    rng = MersenneTwister(seed)
    sample = reshape([binary_distribution(prob, rng) for _ in 1:gate_num], D, div(L, 2))
    initial_state = zeros(length(anyon_basis(L)))
    initial_state[1] = 1.0  # Set the first state as the initial state
    stlis = generate_state(τ, initial_state, sample, temp=temp)

    return stlis
end

function average_Born_sample_p(L::Int64, τ::Float64)
    indexlis = collect(1:2000)
    plis = zeros(length(indexlis))
    seedlis = zeros(length(indexlis))
    for (idx, index) in enumerate(indexlis)
        @show idx
        D, _, _ = get_system_params(τ, L)
        sample,seed = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$L/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "sample", "seed")
        p =  mean(reshape(sample, size(sample,1)*size(sample,2)))
        plis[idx] = p
        seedlis[idx] = seed
    end
    average_p = mean(plis)
    stderr = std(plis) / sqrt(length(indexlis))
    save("exm/data/Bulk_measure/ps_prob/L$(L)_τ$(τ)_Born_average.jld", "average_p", average_p, "stderr", stderr, "seedlis", seedlis)
    return average_p, stderr
end

function get_system_params(τ, L)
    cfg = Dict(
        atanh(0.1)  => (2500L, 1000, 1500L),
        atanh(0.2)  => (500L,  100, 250L),
        atanh(0.3)  => (120L,  48, 100L),
        atanh(0.4)  => (100L,  40, 80L),
        atanh(0.5)  => (80L,   32, 40L),
        atanh(0.6)  => (45L,   20, 30L),
        log(1 + √2) => (35L,   14, 20L),
        atanh(0.8)  => (25L,   10, 10L),
        atanh(0.9)  => (8L,    4, 4L),
        atanh(0.95) => (8L,    4, 4L),
        atanh(0.999)=> (5L,    2, 2L),
    )
    D, step, start = get(cfg, τ, (5L, 2, 2L))
    inds = collect(1:step:D)
    avg_range = start:D-5
    return D, inds, avg_range
end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 
seed_interval_lis = collect(1:100:2000)

if length(ARGS) == 0
    println("No arguments provided.")
else
    # seed=parse(Int64, ARGS[1])
    # stlis = [ps_prob_evolution(L, τ, D, p, seed=seed) for p in problis]
    # eelis = [anyon_eelis(L, st) for st in stlis]
    # save("exm/data/Bulk_measure/ps_prob_evolution/L$(L)_D$(div(D,L))_τ$(τ)_sample_$(seed).jld", "seed", seed, "eelis", eelis, "problis", problis)
    L = parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    τ = τlis[inds]
    println("Computed Born Sample average for L=$L, τ=$τ")
    average_Born_sample_p(L, τ)
end


# prob_eelis = zeros(L-1, length(problis))
# for i in 1:num_samples
#     if i % 100 == 0
#         println("Sample $i of $num_samples")
#     end
#     prob_eelis += hcat(eelis...)    
# end
# save("exm/Bulk_measure/data/ps_prob_sample_$(i).jld", "seed", collect(1:1000))
# prob_eelis ./= num_samples
