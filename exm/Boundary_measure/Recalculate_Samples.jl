using FibonacciChain
using LinearAlgebra
using JLD
using Arpack
using Statistics

function myprint(io::IO, xs...)
    println(io, xs..., '\n')
    flush(io)
end

function Born_Sampling_EE_FE_tau_lis(N, num_samples::Int=1000)
    γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999], 1.0)
    τlis = vcat(atanh.(vcat(collect(0.0:0.05:0.95), [0.99, 0.999])), 1e3)

    @time energy, states = eigs(anyon_ham_sparse(N), nev=1, which=:SR)
    measurement_sites = collect(2:2:N)
    antiGS= states[:, 1]
    
    for (idx, τ) in enumerate(τlis[end])
        myprint(stdout, "N = $N, τ = $τ")
        sample_measured_states, samples, sample_weights = boundary_measure(N, τ, antiGS, measurement_sites, num_samples)
        save("./exm/data/Born_Samples_N$(N)_τ$(τ).jld", "sample_measured_states", sample_measured_states, "samples", samples, "sample_weights", sample_weights)
        

      
    end

    # save("./exm/data/Born_Sampling_EE_FE_tau_lis_N$(N).jld", "average_EE_tau_lis", average_EE_tau_lis, "variance_EE_tau_lis", variance_EE_tau_lis, "average_FE_tau_lis", average_FE_tau_lis, "variance_FE_tau_lis", variance_FE_tau_lis)

end 

function Born_Sampling_rewrite(N::Int64, num_samples::Int=1000)
    γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999], 1.0)
    τlis = vcat(atanh.(vcat(collect(0.0:0.05:0.95), [0.99, 0.999])), 1e3)

    total_mean_EE_tau_lis = Vector{Vector{Float64}}(undef, length(τlis))
    total_stdvar_EE_tau_lis = Vector{Vector{Float64}}(undef, length(τlis))
    total_mean_FE_tau_lis = Vector{Float64}(undef, length(τlis))
    total_stdvar_FE_tau_lis = Vector{Float64}(undef, length(τlis))

    mean_EE_tau_lis, stdvar_EE_tau_lis, mean_FE_tau_lis, stdvar_FE_tau_lis = load("exm/data/Born_Sampling_EE_FE_tau_lis_N$(N).jld", "average_EE_tau_lis", "variance_EE_tau_lis", "average_FE_tau_lis", "variance_FE_tau_lis")
    total_mean_EE_tau_lis[1:end-1] = mean_EE_tau_lis
    total_stdvar_EE_tau_lis[1:end-1] = stdvar_EE_tau_lis
    total_mean_FE_tau_lis[1:end-1] = mean_FE_tau_lis
    total_stdvar_FE_tau_lis[1:end-1] = stdvar_FE_tau_lis
    
    sample_measured_states, samples, sample_weights = load("exm/data/Born_Samples_N$(N)/Born_Samples_N$(N)_τ1000.0.jld", "sample_measured_states", "samples", "sample_weights")

    
    all_EE_values = []
    average_EE_lis=zeros(N-1)
    for i in sample_measured_states
        @time EE_lis = eelis_Fibo_state(N, i)
        myprint(stdout, "N=$(N), Entanglement entropy: $(EE_lis[12])")
        average_EE_lis+= EE_lis
        push!(all_EE_values, EE_lis)
    end

    myprint(stdout, "Number of samples: $(length(sample_measured_states))")
    average_EE_lis = average_EE_lis ./num_samples
    squared_diffs = [(x .- average_EE_lis).^2 for x in all_EE_values]
    stdvar_EE_lis = sqrt.(sum(squared_diffs) ./ num_samples)
    
    total_mean_EE_tau_lis[end] = average_EE_lis
    total_stdvar_EE_tau_lis[end] = stdvar_EE_lis
    # cent, fig = fitCCEntEntScal(average_EE_lis, mincut=2, pbc=true)
    # display(fig)

    random_variable = -log.(sample_weights)
    Sc=mean(random_variable)
    stdvariance_FE = sqrt(var(random_variable))

    total_mean_FE_tau_lis[end] = Sc
    total_stdvar_FE_tau_lis[end] = stdvariance_FE

    save("./exm/data/Boundary_measure/Born_Sampling_EE_FE_tau_lis_N$(N).jld", "mean_EE_tau_lis", total_mean_EE_tau_lis, 
         "stdvar_EE_tau_lis", total_stdvar_EE_tau_lis, 
         "mean_FE_tau_lis", total_mean_FE_tau_lis, 
         "stdvar_FE_tau_lis", total_stdvar_FE_tau_lis)
end

function Born_eelis(N::Int64, τ::Float64, pbc::Bool=true)
    γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999], 1.0)
    τlis = vcat(atanh.(vcat(collect(0.0:0.05:0.95), [0.99, 0.999])), 1e3)

    total_entropieslis = Vector{Vector{Vector{Float64}}}(undef, length(τlis))
    total_probabilitieslis = Vector{Vector{Float64}}(undef, length(τlis))
    total_trajectoriqeslis = Vector{Vector{Vector{Symbol}}}(undef, length(τlis))

    entropieslis, probabilitieslis, trajectorieslis = load("exm/data/Born_enumed_eelis_N$(N).jld", "entropieslis", "probabilitieslis", "trajectorieslis")
    total_entropieslis[1:end-1] = entropieslis
    total_probabilitieslis[1:end-1] = probabilitieslis
    total_trajectorieslis[1:end-1] = trajectorieslis

    energy, states = eigs(anyon_ham_sparse(N), nev=1, which=:SR)
    measurement_sites = collect(2:2:N)
    initial_state= states[:, 1]
    @time final_states, trajectories, probabilities = measurement_enumeration(
        N, τ, initial_state, measurement_sites, pbc)
    
    num_final_states = length(final_states)
    myprint(stdout, "N is $(N), Total number of final states: $(num_final_states)")
    myprint(stdout, "Expected number: $(2^length(measurement_sites))")
    
    total_prob = sum(probabilities)
    myprint(stdout, "Total probability: $(total_prob)")
    
    
    entropies = Vector{Vector{Float64}}(undef, num_final_states)
    for i in 1:num_final_states
        entropies[i] = eelis_Fibo_state(N, final_states[i], pbc)
    end
    

    # avg_eelis = sum(probabilities .* entropies)
    # println("Born-averaged entanglement entropy: $(avg_entropy)")
    
    
    total_entropieslis[end] = entropies
    total_probabilitieslis[end] = probabilities
    total_trajectorieslis[end] = trajectories
    save("./exm/data/Boundary_measure/Born_enumed_eelis_N$(N).jld", "entropieslis", total_entropieslis, "probabilitieslis", total_probabilitieslis, "trajectorieslis", total_trajectorieslis)
end

if length(ARGS) == 0
    println("No arguments provided.")
else
    for arg in ARGS
        println("Received argument: $arg")
        N=parse(Int64, arg)
        # Born_Sampling_EE_FE_tau_lis(N)
        # Born_eelis(N, 1e3)
        Born_Sampling_rewrite(N, 2000)
    end
end