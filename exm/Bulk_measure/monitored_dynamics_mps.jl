using FibonacciChain
using ITensorMPS, ITensors
using JLD
using Statistics
using Random

function samples_generate(L::Int64, τ::Float64, index::Int64, seed::Int64, D::Int64=5L)
    rng = MersenneTwister(seed)
    
    ψ, sites = initial_mps(L)
    
    @time sample_measured_states, sample, sample_free_energy = mps_bulk_measurement(ψ, sites, L, τ, D; rng=rng, pbc=pbc, cutoff=1e-12, maxdim=500)
    
    halfchain_EE_tlis = [ee_mps(j, div(L,2)) for j in sample_measured_states]
    final_state = sample_measured_states[end]
    final_EElis = eelis_Fibo_mps(L, final_state)

    save("exm/data/Bulk_measure/monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld", 
    "sample", sample, "sample_free_energy", sample_free_energy, "seed", seed, 
    "halfchain_EE_tlis", halfchain_EE_tlis, "final_EElis ", final_EElis)
end

function sample_continue_calculate(L::Int64, τ::Float64, index::Int64, seed::Int64, D::Int64=5L, additional_layers::Int64=15L)
    rng = MersenneTwister(seed)
    
    sample, sample_free_energy, seed, halfchain_EE_tlis, final_EElis= load("exm/data/Bulk_measure/monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld", 
    "sample", "sample_free_energy", "seed", "halfchain_EE_tlis", "final_EElis")
    st = generate_state_mps(τ, st, sample, true, true) 
    sample_measured_states, sample, sample_free_energy = mps_bulk_measurement(st[end-1], sites, L, τ, additional_layers; rng=rng, pbc=pbc, cutoff=1e-12, maxdim=500) 
    final_state = sample_measured_states[end]

    new_halfchain_EE_tlis = vcat(halfchain_EE_tlis, [ee_mps(j, div(L,2)) for j in sample_measured_states])
    new_final_EElis = vcat(final_EElis, eelis_Fibo_mps(L, final_state))

    
    save("exm/data/Bulk_measure/monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld", 
    "sample", sample, "sample_free_energy", sample_free_energy, "seed", seed, 
    "halfchain_EE_tlis", new_halfchain_EE_tlis, "final_EElis ", new_final_EElis)
end

function samples_collect(L::Int64, τ::Float64, D::Int64=5L)
    samples_num = 2000
    ensemble = Vector{Matrix{Int}}(undef, samples_num)
    ensemble_free_energy = Vector{Vector{Float64}}(undef, samples_num)
    ensemble_seed = Vector{Int64}(undef, samples_num)
    ensemble_EE_dynamics= zeros(samples_num, D) 
    ensemble_final_EElis = zeros(samples_num, L-1)

     for i in 1:samples_num
        @show i
        @time sample, sample_free_energy, seed, halfchain_EE_tlis, final_EElis = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(i).jld", "sample", "sample_free_energy", "seed", "halfchain_EE_tlis", "final_EElis")
        ensemble[i] = sample
        ensemble_free_energy[i] = sample_free_energy
        ensemble_seed[i] = seed
        ensemble_EE_dynamics[i, :] = halfchain_EE_tlis
        ensemble_final_EElis[i, :] = final_EElis
    end

    bulk_meanEElis = mean(ensemble_final_EElis, dims=1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims=1)[:]
    ensemble_stderr_EElis = (std(ensemble_final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims=1) ./ sqrt(samples_num))[:]

    save("exm/data/Bulk_measure/monitored_dynamics_ensemble_L$(L)_τ$(τ)_D$(div(D,L)).jld", 
    "ensemble", ensemble, "ensemble_free_energy", ensemble_free_energy, "ensemble_seed", ensemble_seed,  
    "average_EE_tlis", average_EE_tlis, "stderr_EE_tlis", stderr_EE_tlis, 
    "bulk_meanEElis", bulk_meanEElis, "ensemble_stderr_EElis",ensemble_stderr_EElis)
end

function monitored_dynamics(L::Int64, τ::Float64, D::Int64=5L)
    ψ, sites = initial_mps(L)
    
    samples_num = 2000

    ensemble_EE_dynamics= zeros(samples_num, D) 
    ensemble_stderr_EElis = zeros(L-1)
    final_EElis = zeros(samples_num, L-1)

    all_FE_tlis = zeros(samples_num, D)
    final_FElis = zeros(samples_num)
    for i in 1:samples_num
        @show i
        @time sample_measured_states, sample, sample_free_energy = mps_bulk_measurement(ψ, sites, L, τ, D; rng=rng, pbc=pbc, cutoff=1e-12, maxdim=500)
    
        ensemble_EE_dynamics[i, :] = [ee_mps(j, div(L,2)) for j in sample_measured_states]
        final_state = sample_measured_states[end]
        final_EElis[i, :] = eelis_Fibo_mps(L, final_state)

        all_FE_tlis[i, :] = sample_free_energy
        final_FElis[i] = sample_free_energy[end]
    end

    bulk_meanEElis = mean(final_EElis, dims=1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims=1)[:]
    ensemble_stderr_EElis = (std(final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims=1) ./ sqrt(samples_num))[:]
    
    average_FE = mean(final_FElis)
    stderr_FE = std(final_FElis) / sqrt(samples_num)
    average_FE_tlis = mean(all_FE_tlis, dims=1)[:]
    stderr_FE_tlis = (std(all_FE_tlis, dims=1) ./ sqrt(samples_num))[:]

    
    return average_EE_tlis, stderr_EE_tlis, bulk_meanEElis, ensemble_stderr_EElis, 
           average_FE_tlis, stderr_FE_tlis, average_FE, stderr_FE
end

τ = 1000.0

if length(ARGS) == 0
    println("No arguments provided.")
else
    N=parse(Int64, ARGS[1])
    index=parse(Int64, ARGS[2])
    seed=parse(Int64, ARGS[3])
    println("Received argument: $N, $index")
    samples_generate(N, τ, index, seed)
    # samples_collect(N, τ)
end