using FibonacciChain
using LinearAlgebra
using JLD
using Arpack
include("../FitEntEntScal.jl")

N= 18
γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999])
τlis = atanh.(γlis)

average_EE_tau_lis=Vector{Vector{Float64}}(undef, length(τlis))
for (idx, τ) in enumerate(τlis)
    @show τ
    energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
    antiGS= states[:, 1]

    measurement_sites = collect(2:2:N)
    sample_measured_states, samples, sample_weights = Sampling(N, τ, antiGS, measurement_sites)

    average_EE_lis=zeros(N-1)
    for i in sample_measured_states
        average_EE_lis+= eelis_Fibo_state(N, i)
    end
    average_EE_lis = average_EE_lis ./1000
    average_EE_tau_lis[idx] = average_EE_lis
    # cent, fig = fitCCEntEntScal(average_EE_lis, mincut=2, pbc=true)
    # display(fig)
end

save("./exm/data/Born_Sampling_ee_tau_lis_N$(N).jld", "average_EE_tau_lis", average_EE_tau_lis)

