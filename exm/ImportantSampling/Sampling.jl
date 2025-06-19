using FibonacciChain
using LinearAlgebra
using JLD
using Arpack
# include("../FitEntEntScal.jl")

# N= 22

function myprint(io::IO, xs...)
    println(io, xs..., '\n')
    flush(io)
end

function Born_Sampling_ee_tau_lis(N, num_samples::Int=1000)
    γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999])
    τlis = atanh.(γlis)
    average_EE_tau_lis=Vector{Vector{Float64}}(undef, length(τlis))
    variance_EE_tau_lis = Vector{Vector{Float64}}(undef, length(τlis))

    @time energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
    measurement_sites = collect(2:2:N)

    for (idx, τ) in enumerate(τlis)
        myprint(stdout, "N = $N, τ = $τ")
        antiGS= states[:, 1]
        sample_measured_states, samples, sample_weights = Sampling(N, τ, antiGS, measurement_sites, num_samples)

        all_EE_values = []
        average_EE_lis=zeros(N-1)
        for i in sample_measured_states
            EE_lis = eelis_Fibo_state(N, i)
            average_EE_lis+= EE_lis
            push!(all_EE_values, EE_lis)
        end
        average_EE_lis = average_EE_lis ./num_samples
        squared_diffs = [(x .- average_EE_lis).^2 for x in all_EE_values]
        variance_EE_lis = sqrt.(sum(squared_diffs) ./ num_samples)
        
        average_EE_tau_lis[idx] = average_EE_lis
        variance_EE_tau_lis[idx] = variance_EE_lis
        # cent, fig = fitCCEntEntScal(average_EE_lis, mincut=2, pbc=true)
        # display(fig)
    end

    save("./exm/data/Born_Sampling_ee_tau_lis_N$(N).jld", "average_EE_tau_lis", average_EE_tau_lis, "variance_EE_tau_lis", variance_EE_tau_lis)

end 

if length(ARGS) == 0
    println("No arguments provided.")
else
    for arg in ARGS
        println("Received argument: $arg")
        N=parse(Int64, arg)
        Born_Sampling_ee_tau_lis(N)
    end
end
