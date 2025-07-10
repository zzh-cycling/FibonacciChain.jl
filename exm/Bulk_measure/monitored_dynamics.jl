using FibonacciChain
using JLD
using Statistics
using Random

function samples_generate(L::Int64, τ::Float64, index::Int64, seed::Int64, D::Int64=35L)
    rng = MersenneTwister(seed)
    
    st = zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    
    @time sample_measured_states, sample, sample_free_energy = Bulkmeasure(L, τ, st, D, rng) 
    
    save("exm/data/Bulk_measure/Samples_monitored_dynamics/L$(L)/τ$(τ)_D$(div(D,L))_Samples$(index).jld", "sample", sample, "sample_free_energy", sample_free_energy, "seed", seed)
    # return sample_measured_states, samples, sample_free_energy
end

function samples_collect(L::Int64, τ::Float64, D::Int64=35L)
    samples_num = 2000
    ensemble = Vector{Matrix{Int}}(undef, samples_num)
    ensemble_free_energy = Vector{Vector{Float64}}(undef, samples_num)
    ensemble_seed = Vector{Int64}(undef, samples_num)
     for i in 1:samples_num
        @show i
        @time sample, sample_free_energy, seed = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$(L)/τ$(τ)_D$(div(D,L))_Samples$(index).jld", "sample", "sample_free_energy", "seed")
        ensemble[i] = sample
        ensemble_free_energy[i] = sample_free_energy
        ensemble_seed[i] = sample_free_energy["seed"]
    end

    save("exm/data/Bulk_measure/monitored_dynamics_ensemble_L$(L)_τ$(τ)_D$(div(D,L)).jld", "ensemble", ensemble, "ensemble_free_energy", sensemble_free_energy, "ensemble_seed", ensemble_seed)
end

function samples_calculate(L::Int64, τ::Float64, index::Int64, D::Int64=35L)
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    
    sample, sample_free_energy, seed= load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$(L)/τ$(τ)_D$(div(D,L))_Samples$(index).jld", "sample", "sample_free_energy","seed")

    # evolution by a measurement gate, Born's rule
    sample_measured_states = Generate_state(τ, st, sample, true, true) 
    # sample_free_energy = D-dimensional vector of (-ln conditional_P_of_each_layer) = free energy contributed by each layer of L/2 gates (PBC)
    halfchain_EE_tlis = [eelis_Fibo_state(L, j)[div(L,2)] for j in sample_measured_states]
    final_state = sample_measured_states[end]
    final_EElis = eelis_Fibo_state(L, final_state)

    
    save("./exm/data/Bulk_measure/Observable_monitored_dynamics/L$(L)/τ$(τ)_D$(div(D,L))_Samples$(index).jld", "halfchain_EE_tlis", halfchain_EE_tlis, "final_EElis ", final_EElis, "seed", seed, "sample_free_energy", sample_free_energy)
end

function Observable_collect(L::Int64, τ::Float64, D::Int64=35L)
    samples_num = 2000

    bulk_meanEElis=zeros(L-1)
    ensemble_seed = Vector{Int64}(undef, samples_num)
    ensemble_EE_dynamics= zeros(samples_num, D) 
    final_EElis = zeros(samples_num, L-1)

    for i in 1:samples_num
        @show i
        halfchain_EE_tlis, final_EElis, seed, sample_free_energy = load("./exm/data/Bulk_measure/Observable_monitored_dynamics/L$(L)/τ$(τ)_D$(div(D,L))_Samples$(index).jld", "halfchain_EE_tlis", "final_EElis ", "seed",  "sample_free_energy")

        all_EE_tlis[i, :] = halfchain_EE_tlis
        final_EElis[i, :] = final_EElis
        ensemble_seed[i] = seed
    end

    bulk_meanEElis = mean(final_EElis, dims=1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims=1)[:]
    ensemble_stderr_EElis = (std(final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims=1) ./ sqrt(samples_num))[:]

    
    save("exm/data/Bulk_measure/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld", "average_EE_tlis", average_EE_tlis, "stderr_EE_tlis", stderr_EE_tlis, "bulk_meanEElis", bulk_meanEElis, "ensemble_stderr_EElis",ensemble_stderr_EElis, "ensemble_free_energy", ensemble_free_energy, "ensemble_seed", ensemble_seed)
end

function ensemble_calculate(L::Int64, τ::Float64, D::Int64=35L)
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    bulk_meanEElis=zeros(L-1)
    
    samples_num = 2000

    ensemble_EE_dynamics= zeros(samples_num, D) 
    final_EElis = zeros(samples_num, L-1)


    ensemble, ensemble_free_energy = load("./exm/data/Bulk_measure/monitored_dynamics_ensemble_L$(L)_τ$(τ)_D$(div(D,L)).jld", "ensemble", "ensemble_free_energy")

    for i in 1:samples_num
        @show i
        # evolution by a measurement gate, Born's rule
        sample_measured_states = Generate_state(τ, st, sample, true, true)  
        # sample_free_energy = D-dimensional vector of (-ln conditional_P_of_each_layer) = free energy contributed by each layer of L/2 gates (PBC)

        all_EE_tlis[i, :] = [eelis_Fibo_state(L, j)[div(L,2)] for j in sample_measured_states]
        final_state = sample_measured_states[end]
        final_EElis[i, :] = eelis_Fibo_state(L, final_state)
    end

    bulk_meanEElis = mean(final_EElis, dims=1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims=1)[:]
    ensemble_stderr_EElis = (std(final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims=1) ./ sqrt(samples_num))[:]

    
    save("exm/data/Bulk_measure/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld", "average_EE_tlis", average_EE_tlis, "stderr_EE_tlis", stderr_EE_tlis, "bulk_meanEElis", bulk_meanEElis, "ensemble_stderr_EElis",ensemble_stderr_EElis, "ensemble_free_energy", ensemble_free_energy)
end

function monitored_dynamics(L::Int64, τ::Float64, D::Int64=35L)
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    bulk_meanEElis=zeros(L-1)
    
    samples_num = 2000

    ensemble_EE_dynamics= zeros(samples_num, D) 
    ensemble_stderr_EElis = zeros(L-1)
    final_EElis = zeros(samples_num, L-1)

    all_FE_tlis = zeros(samples_num, D)
    final_FElis = zeros(samples_num)
    for i in 1:samples_num
        @show i
        sample_measured_states, sample, sample_free_energy = Bulkmeasure(L, τ, st, D) 
        ensemble_EE_dynamics[i, :] = [eelis_Fibo_state(L, j)[div(L,2)] for j in sample_measured_states]
        final_state = sample_measured_states[end]
        final_EElis[i, :] = eelis_Fibo_state(L, final_state)

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



## == monitored_dynamics == ##
# L = 10
τ = log(1+√2)
# D = 100L

# ensemble_EE_dynamics, average_EE_tlis, stderr_EE_tlis, all_EElis, bulk_meanEElis, ensemble_stderr_EElis = monitored_dynamics(L, τ, D)


if length(ARGS) == 0
    println("No arguments provided.")
else
    N=parse(Int64, ARGS[1])
    index=parse(Int64, ARGS[2])
    # seed=parse(Int64, ARGS[3])
    println("Received argument: $N, $index")
    # samples_generate(N, τ, index, seed)
    samples_calculate(N, τ, index)
end