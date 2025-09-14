using JLD
using FibonacciChain

N=10
energy, states = eigen(anyon_ham(N))
antiGS= states[:, 1]
τ = 1000.0
measurement_sites = collect(2:2:N)



sample_measured_states,  samples1000, sample_weights = boundary_measure(N, τ, antiGS, measurement_sites)
sample_measured_states, samples10000, sample_weights = load("exm/data/Born_Samples_N10_τ1000.0.jld", "sample_measured_states",  "samples","sample_weights")
sample_measured_states,  samples100000, sample_weights = boundary_measure(N, τ, antiGS, measurement_sites, 100000)

final_states, trajectories, probabilities = measurement_enumeration(N, τ, antiGS, measurement_sites)

binary_digits_enum = map.(symbol -> symbol == :p ? 0 : 1, trajectories)
decimal_value_enum = [sum(d * 2^(length(j) - i) for (i, d) in enumerate(j)) for j in binary_digits_enum]

binary_digits_samples10000 = map.(symbol -> symbol == :p ? 0 : 1, samples10000)
decimal_value_samples10000 = [sum(d * 2^(length(j) - i) for(i, d) in enumerate(j)) for j in binary_digits_samples10000]
count_lis10000 = [count(x->x==i, decimal_value_samples10000) for i in 0:31] ./10000

binary_digits_samples1000 = map.(symbol -> symbol == :p ? 0 : 1, samples1000)
decimal_value_samples1000 = [sum(d * 2^(length(j) - i) for(i, d) in enumerate(j)) for j in binary_digits_samples1000]
count_lis1000 = [count(x->x==i, decimal_value_samples1000) for i in 0:31] ./1000


binary_digits_samples100000 = map.(symbol -> symbol == :p ? 0 : 1, samples100000)
decimal_value_samples100000 = [sum(d * 2^(length(j) - i) for(i, d) in enumerate(j)) for j in binary_digits_samples100000]
count_lis100000 = [count(x->x==i, decimal_value_samples100000) for i in 0:31] ./100000

