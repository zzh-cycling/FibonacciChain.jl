using FibonacciChain
using JLD
using Random
using Statistics


binary_distribution(p, rng) = rand(rng) < p ? 1 : 0
function ps_prob_evolution(L::Int64, τ::Float64, D::Int; seed::Int=100)
    # D is the number of layers, not period
    # Generate a state based on the given sample and initial state
    # τ is the temperature parameter, initial_state is the initial state vector
    # sample is a matrix of binary values representing the state configuration
    problis = collect(0.1:0.1:0.9)
    ee_plis = Vector{Vector{Float64}}(undef, length(problis))
    initial_state = zeros(length(anyon_basis(L)))
    initial_state[1] = 1.0  # Set the first state as the initial state
    gate_num = div(D*L, 2)

    for (idx, prob) in enumerate(problis)
        @show idx

        rng = MersenneTwister(seed)
        sample = reshape([binary_distribution(prob, rng) for _ in 1:gate_num], D, div(L, 2))
        stlis, sample_return = generate_state(τ, initial_state, sample)
        ee_plis[idx] = anyon_eelis(L, stlis[end])
    end
    save("exm/data/Bulk_measure/ps_prob_evolution/L$(L)/τ$(τ)/L$(L)_D$(div(D,L))_τ$(τ)_sample$(seed).jld", "seed", seed, "ee_plis", ee_plis, "problis", problis)
    return ee_plis
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

function process_ps_prob_evolution(L, τ, D)
    # Load the data
    samplelis = collect(1:1000)
    problis = collect(0.1:0.1:0.9)
    centlis = []
    seedlis = zeros(Int64, length(samplelis))
    # Process the data
    ee_problis = zeros(L-1, length(problis))
    for (i, sample) in enumerate(samplelis)
        ee_plis, seed = load("exm/data/Bulk_measure/ps_prob_evolution/L$(L)/τ$(τ)/L$(L)_D$(div(D,L))_τ$(τ)_sample$(sample).jld", "ee_plis", "seed")
        ee_problis += hcat(ee_plis...)
        # cent_problis[i, j] = 
        seedlis[i] = seed
    end

    ee_problis ./= length(samplelis)
    

    for i in eachindex(problis)
        push!(centlis, fitCCEntEntScal(vec(ee_problis[:, i]), mincut=2, pbc=true)[1])
    end

    save("exm/data/Bulk_measure/ps_prob_evolution/L$(L)/τ$(τ)/L$(L)_D$(div(D,L))_τ$(τ)_cent.jld", "centlis", centlis, "seedlis", seedlis)
    return centlis, seedlis
end

function process_data()
    Llis = collect(8:2:16)
    problis = collect(0.1:0.1:0.9)
    ixs = [1, 3, 4, 7, 9, 10, 12]
    for (inds, τ) in enumerate(τlis[ixs])
        cent_Lplis = zeros(length(Llis), length(problis))
        cent_stderrlis = zeros(length(Llis), length(problis))
        for (id, L) in enumerate(Llis)
            D, _, _ = get_system_params(τ, L)
            @show (L, τ)
            centlis= load("exm/data/Bulk_measure/ps_prob_evolution/L$(L)/τ$(τ)/L$(L)_D$(div(D,L))_τ$(τ)_cent.jld", "centlis")
            cent_Lplis[id, :] = [i[1] for i in centlis]
            cent_stderrlis[id, :] = [i[2] for i in centlis]
        end
        save("exm/data/Bulk_measure/ps_prob_evolution/centlis_L$(Llis[1])$(Llis[end])_τ$(τ).jld", "cent_Lplis", cent_Lplis, "cent_stderrlis", cent_stderrlis)
    end
end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 
seed_interval_lis = collect(1:100:2000)

if length(ARGS) == 0
    println("No arguments provided.")
else
    L = parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    seed=parse(Int64, ARGS[3])
    τ = τlis[inds]
    D, _, _ = get_system_params(τ, L)
    # println("Computed Born Sample average for L=$L, τ=$τ")
    println("Computed Post-selection average for seed=$seed, τ=$τ")
    # average_Born_sample_p(L, τ)
    ps_prob_evolution(L, τ, D; seed=seed)
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
