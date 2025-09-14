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
initial_state = normalize!(ones(ComplexF64, length(anyon_basis(N))))

# Perform measurement at site 2 with outcome 0 (+)
final_state = measuremap(N, τ, initial_state, 2, 0)
```

### Measurement Layer Protocol

```julia
# Define measurement sites (even sites for Fibonacci anyons)
measurement_sites = 2:2:N

# Generate random measurement layer
layer_sample = rand([0, 1], length(measurement_sites))

# Apply measurement layer
evolved_state = apply_measurement_layer!(N, initial_state, τ, layer_sample, 1)
```

### Boundary Measurement Sampling

```julia
# Generate 1000 measurement samples
measurement_sites = [2, 4, 6]
states, sequences, free_energies = boundary_measure(N, τ, initial_state, measurement_sites, 1000)

# Analyze free energy distribution
using Statistics
mean_free_energy = mean(free_energies)
```
