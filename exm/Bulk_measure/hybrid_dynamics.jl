using FibonacciChain
using Random

# Every site of the usual staggered Fibonacci spacetime lattice is chosen
# independently: measurement with probability p, U(θ) with probability 1-p.
N = 12
p = 0.35
θ = π / 2

model = AnyonModel(FibonacciAnyon(), N; pbc = true)
state = zeros(length(anyon_basis(model)))
state[1] = 1.0

born_config = HybridConfig(
    τ = 1.0,
    t₂ = 20,
    mode = :Born,
    rng = MersenneTwister(42),
    enable_τ_eff = false,
    p = p,
    θ = θ,
)
trajectory = bulk_evolution(model, state, born_config)

# The gate mask and measurement outcomes together define the full spacetime
# realization and can be replayed exactly (or passed to the MPS backend).
replay_config = HybridConfig(
    τ = 1.0,
    t₂ = 20,
    mode = :sample,
    enable_τ_eff = false,
    p = p,
    θ = θ,
)
replayed = bulk_evolution(model, state, replay_config, trajectory.schedule)

@assert replayed.state ≈ trajectory.state
@assert replayed.free_energys ≈ trajectory.free_energys

println("measurement fraction = ", sum(trajectory.schedule.measurement_mask) / length(trajectory.schedule.measurement_mask))
println("final half-chain entropy = ", last(trajectory.entanglement_entropys))
