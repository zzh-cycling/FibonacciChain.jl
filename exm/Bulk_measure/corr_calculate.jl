using FibonacciChain
using JLD
using Statistics
using Random

sample = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L12/τ0.8813735870195429/D35_Samples1.jld", "sample")

N = 12
pbc = true
τ= log(1+√2)
measure_class = :Fibo
initial_state = zeros(length(Fibonacci_basis(BitStr{N, Int}, pbc, measure_class=measure_class)))
initial_state[1] = 1.0 # initial state is all zero state

statelis = generate_state(0.5, initial_state, sample, true, temp= true)

final_st= statelis[end]
spatial_corr = spatial_correlation(N, final_st, 1, 5, pbc)

temporal_corr = temporal_correlation(τ, initial_state, sample, div(N,2), 30N, 31N)