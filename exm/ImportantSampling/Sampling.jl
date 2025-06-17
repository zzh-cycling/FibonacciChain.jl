using FibonacciChain
using LinearAlgebra
using JLD
using Arpack
include("../FitEntEntScal.jl")


function compute_sample_average(samples::Vector{Vector{Symbol}}, sample_weights::Vector{Float64}, observable_func::Function)
    num_samples = length(samples)
    total = 0.0
    
    for i in 1:num_samples
        observable_value = observable_func(samples[i])
        total += observable_value
    end
    
    return total / num_samples
end


N=16
energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
antiGS= states[:, 1]
