using FibonacciChain
include("../FitEntEntScal.jl")

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

τ = log(1 + √2) 
L=20;D=20L 
problis = collect(0.1:0.1:0.99)

# num_samples = 1000

if length(ARGS) == 0
    println("No arguments provided.")
else
    seed=parse(Int64, ARGS[1])
    stlis = [ps_prob_evolution(L, τ, D, p, seed=seed) for p in problis]
    eelis = [anyon_eelis(L, st) for st in stlis]
    save("exm/Bulk_measure/data/ps_prob_evolution/L$(L)_D$(div(D,L))_τ$(τ)_sample_$(seed).jld", "seed", seed, "eelis", eelis, "problis", problis)
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

# centlis = zeros(length(problis))
# centerrlis = zeros(length(problis))

# for (i, ee) in enumerate(eelis)
#     cent, fig = fitCCEntEntScal(ee, mincut=2, pbc=true)
#     centlis[i] = cent[1]
#     centerrlis[i] = cent[2]
#     savefig(fig, "./exm/fig/ps_prob_evocc_L$(L)_D$(D)_p$(problis[i]).pdf")    
# end