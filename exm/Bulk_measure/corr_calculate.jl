using FibonacciChain
using JLD
using Statistics
using Random
using BitBasis
using BenchmarkTools

sample = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L12/τ0.8813735870195429/D35_Samples1.jld", "sample")
halfchain_EE_tlis = load("exm/data/Bulk_measure/Observable_monitored_dynamics/L12/τ0.8813735870195429/D35_Samples1.jld", "halfchain_EE_tlis")

N = 12
pbc = true
τ= log(1+√2)
measure_class = :Fibo
initial_state = zeros(length(Fibonacci_basis(BitStr{N, Int}, pbc, measure_class=measure_class)))
initial_state[1] = 1.0 # initial state is all zero state

statelis = generate_state(0.5, initial_state, sample, true, temp= true)

final_st= statelis[end-20]
spatial_corr = spatial_correlation(N, final_st, 1, div(N,2), pbc)

timeslice1 = 35*8
temporal_corr_lis = [temporal_correlation(τ, initial_state, sample, div(N,2), timeslice1, j) for j in timeslice1+1:35*12-5]

# 0.8s per run

plot(collect(1:length(temporal_corr_lis))./N, temporal_corr_lis ./spatial_corr, xlabel=L"\Delta t /L", ylabel=L"g(0, \Delta t)/g_{space}", label=false)