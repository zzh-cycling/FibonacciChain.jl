using JLD
using FibonacciChain

N=10
model = AnyonModel(FibonacciAnyon(), N; pbc = true)
energy, states = eigen(anyon_ham(model))
antiGS = states[:, 1]
τ = 1000.0
measurement_sites = collect(2:2:N)

# Note: boundary_measure is now boundary_evolution, need to run multiple times for sampling
config = MeasureConfig(τ = τ, mode = :Born)

# Generate 1000 samples
samples1000 = Vector{BitVector}(undef, 1000)
sample_weights1000 = Vector{Float64}(undef, 1000)
sample_measured_states1000 = Vector{Vector{Float64}}(undef, 1000)
for i = 1:1000
    outcome = boundary_evolution(model, antiGS, config)
    sample_measured_states1000[i] = outcome.state
    samples1000[i] = outcome.sample
    sample_weights1000[i] = exp(-outcome.free_energy)
end

sample_measured_states, samples10000, sample_weights = load(
    "exm/data/Born_Samples_N10_τ1000.0.jld",
    "sample_measured_states",
    "samples",
    "sample_weights",
)

# Generate 100000 samples
samples100000 = Vector{BitVector}(undef, 100000)
sample_weights100000 = Vector{Float64}(undef, 100000)
for i = 1:100000
    outcome = boundary_evolution(model, antiGS, config)
    samples100000[i] = outcome.sample
    sample_weights100000[i] = exp(-outcome.free_energy)
end

final_states, trajectories, probabilities =
    measurement_enumeration(model, τ, antiGS, measurement_sites)

binary_digits_enum = map.(symbol -> symbol == :p ? 0 : 1, trajectories)
decimal_value_enum =
    [sum(d * 2^(length(j) - i) for (i, d) in enumerate(j)) for j in binary_digits_enum]

binary_digits_samples10000 = map.(symbol -> symbol == :p ? 0 : 1, samples10000)
decimal_value_samples10000 = [
    sum(d * 2^(length(j) - i) for (i, d) in enumerate(j)) for
    j in binary_digits_samples10000
]
count_lis10000 = [count(x->x==i, decimal_value_samples10000) for i = 0:31] ./ 10000

binary_digits_samples1000 = map.(symbol -> symbol == :p ? 0 : 1, samples1000)
decimal_value_samples1000 = [
    sum(d * 2^(length(j) - i) for (i, d) in enumerate(j)) for j in binary_digits_samples1000
]
count_lis1000 = [count(x->x==i, decimal_value_samples1000) for i = 0:31] ./ 1000


binary_digits_samples100000 = map.(symbol -> symbol == :p ? 0 : 1, samples100000)
decimal_value_samples100000 = [
    sum(d * 2^(length(j) - i) for (i, d) in enumerate(j)) for
    j in binary_digits_samples100000
]
count_lis100000 = [count(x->x==i, decimal_value_samples100000) for i = 0:31] ./ 100000
