using FibonacciChain
using LinearAlgebra
using JLD
using Arpack
# include("../FitEntEntScal.jl")

γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999], 1.0)
τlis = vcat(atanh.(vcat(collect(0.0:0.05:0.95), [0.99, 0.999])), 1e3)

# N= 22

function myprint(io::IO, xs...)
    println(io, xs..., '\n')
    flush(io)
end

function Born_Sampling_EE_FE_tau_lis(N, num_samples::Int=1000)
    mean_EE_tau_lis=Vector{Vector{Float64}}(undef, length(τlis))
    stderr_EE_tau_lis = Vector{Vector{Float64}}(undef, length(τlis))
    mean_FE_tau_lis = Vector{Float64}(undef, length(τlis))
    stderr_FE_tau_lis = Vector{Float64}(undef, length(τlis))

    model = AnyonModel(FibonacciAnyon(), N; pbc=true)
    @time energy, states = eigs(anyon_ham_sparse(model), nev=1, which=:SR)
    measurement_sites = collect(2:2:N)

    for (idx, τ) in enumerate(τlis)
        antiGS = states[:, 1]
        myprint(stdout, "N = $N, τ = $τ")
        
        # Initialize arrays to store results from all samples
        all_EE_values = zeros(num_samples, N-1)
        all_FE_values = zeros(num_samples)
        
        @time for i in 1:num_samples
            config = MeasureConfig(τ=τ, mode=:Born)
            outcome = boundary_evolution(model, antiGS, config)
            sample_measured_states = [outcome.state]
            samples = outcome.sample
            sample_weights = [exp(-outcome.free_energy)]
            save("./exm/data/Born_Samples_N$(N)_τ$(τ)_sample$(i).jld", "sample_measured_states", sample_measured_states, "samples", samples, "sample_weights", sample_weights)
            
            # Store EE values for this sample
            all_EE_values[i, :] = anyon_eelis(model, sample_measured_states[1])
            
            # Store FE value for this sample
            all_FE_values[i] = -log(sample_weights[1])
        end
        
        # Calculate statistics for EE
        average_EE_lis = mean(all_EE_values, dims=1)[:]
        stderr_EE_lis = std(all_EE_values, dims=1)[:] ./ sqrt(num_samples)
        
        mean_EE_tau_lis[idx] = average_EE_lis
        stderr_EE_tau_lis[idx] = stderr_EE_lis
        
        # Calculate statistics for FE
        average_FE = mean(all_FE_values)
        stderr_FE = std(all_FE_values) / sqrt(num_samples)
        
        mean_FE_tau_lis[idx] = average_FE
        stderr_FE_tau_lis[idx] = stderr_FE
    end

    save("./exm/data/Born_Sampling_EE_FE_tau_lis_N$(N).jld", "mean_EE_tau_lis", mean_EE_tau_lis, "stderr_EE_tau_lis", stderr_EE_tau_lis, "mean_FE_tau_lis", mean_FE_tau_lis, "stderr_FE_tau_lis", stderr_FE_tau_lis)
end

if length(ARGS) == 0
    println("No arguments provided.")
else
    for arg in ARGS
        println("Received argument: $arg")
        N=parse(Int64, arg)
        Born_Sampling_EE_FE_tau_lis(N)
    end
end
