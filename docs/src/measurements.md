# Measurement and Dynamics

This section covers quantum measurement protocols and dynamical evolution of anyon chains.

## Single-Site Measurements

### Basic Measurement Operations


### Ladder System Measurements


## Measurement Protocols

### Boundary Measurements

Boundary measurements are performed at the edges of the chain, allowing for:
- Probabilistic measurement outcomes
- Post-selection of specific measurement results
- Free energy calculation for each trajectory

### Bulk Measurements

### State Evolution


## Measurement Trees and Enumeration


These functions provide systematic ways to explore all possible measurement outcomes and their probabilities.

## Noise Channels

### Probabilistic Braiding


The Choi map implements a noise channel that applies braiding operations with probability p:

```math
\mathcal{E}[\rho] = (1-p)\rho + p B[\rho]
```

where B is the braiding operation and ρ is the density matrix.

## Usage Examples

### Single Measurement

```julia
using FibonacciChain

N = 6
τ = 1.0  # evolution time parameter
model = AnyonModel(FibonacciAnyon(), N; pbc=true)
initial_state = normalize!(ones(ComplexF64, length(anyon_basis(model))))

# Perform measurement at site 2 with outcome 0 (+)
final_state = measuremap(model, τ, initial_state, 2, false)
```

### Measurement Layer Protocol

```julia
# Apply a single boundary measurement layer (Born mode)
out = boundary_evolution(
	model,
	initial_state,
	MeasureConfig(τ=τ, t₁=1, t₂=1, mode=:Born, rng=MersenneTwister()),
)
evolved_state = out.state
layer_sample = out.sample
layer_free_energy = out.free_energy
```

### Boundary Measurement Sampling

```julia
# Sample multiple layers with bulk evolution (Born rule)
D = 10  # number of time steps
config = MeasureConfig(τ=τ, t₁=1, t₂=D, mode=:Born, rng=MersenneTwister())
mo = bulk_evolution(model, initial_state, config)

# Analyze free energy across layers
using Statistics
mean_free_energy = mean(mo.free_energys)
```
