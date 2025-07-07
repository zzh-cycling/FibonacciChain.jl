# MPS-based Fibonacci Chain Implementation

This document describes the MPS (Matrix Product State) implementation for Fibonacci chain measurements using ITensor.

## Overview

The MPS implementation provides an efficient way to handle larger system sizes compared to exact diagonalization, while maintaining high accuracy for ground state properties and measurement dynamics.

## Key Features

### 1. Ground State Generation
- `fibonacci_mps_ground_state(N; pbc=true)`: Generates the ground state of the Fibonacci chain using DMRG
- Returns both the MPS state and ground state energy
- Supports periodic boundary conditions (PBC) and open boundary conditions (OBC)

### 2. Measurement Operators
- `measurement_operator_mps(sites, i, τ, sign; pbc=true)`: Creates measurement operators as MPO
- Supports both `:p` (plus) and `:m` (minus) measurement outcomes
- Handles boundary conditions and local measurement rules

### 3. State Evolution
- `apply_measurement_mps(ψ, sites, i, τ, sign; pbc=true)`: Applies measurement to MPS state
- Returns the evolved state and measurement probability
- Automatically normalizes the resulting state

### 4. Sampling and Enumeration
- `mps_sampling(ψ, sites, measurement_sites, τ; num_samples=1000, pbc=true)`: Monte Carlo sampling
- `mps_measurement_enumeration(ψ, sites, measurement_sites, τ; pbc=true)`: Exact enumeration
- `mps_bulk_measurement(ψ, sites, N, τ, D; pbc=true)`: Bulk measurement with alternating pattern

### 5. Entanglement Analysis
- `calculate_entanglement_entropy_mps(ψ, b)`: Calculates von Neumann entanglement entropy at bond b

## Usage Examples

### Basic Ground State Calculation

```julia
using FibonacciChain

# Generate ground state for 10-site chain
N = 10
ψ, energy = fibonacci_mps_ground_state(N; pbc=true)
println("Ground state energy: $energy")
println("Bond dimensions: $(linkdims(ψ))")
```

### Single Measurement

```julia
# Get site indices
sites = siteinds(ψ)

# Apply measurement at site 2
τ = 1.0
ψ_measured, prob = apply_measurement_mps(ψ, sites, 2, τ, :p; pbc=true)
println("Measurement probability: $prob")
```

### Sampling Measurements

```julia
# Define measurement sites
measurement_sites = [2, 4, 6, 8]
num_samples = 1000

# Perform sampling
samples, weights = mps_sampling(ψ, sites, measurement_sites, τ; 
                               num_samples=num_samples, pbc=true)

# Analyze results
println("Average weight: $(sum(weights)/length(weights))")
```

### Entanglement Entropy

```julia
# Calculate entanglement entropy at different cuts
for b in 2:(N-1)
    ee = calculate_entanglement_entropy_mps(ψ, b)
    println("EE at cut $b: $ee")
end
```

## Comparison with Exact Diagonalization

| Aspect | Exact Diagonalization | MPS |
|--------|----------------------|-----|
| System Size | Limited (~20 sites) | Large (>100 sites) |
| Memory Usage | Exponential | Polynomial |
| Ground State | Exact | Highly accurate |
| Measurements | Exact probabilities | Accurate probabilities |
| Entanglement | Exact | Accurate |

## Performance Considerations

### Bond Dimension
- Controls the accuracy vs. computational cost trade-off
- Default maximum bond dimension: 200
- Can be adjusted in DMRG settings

### Cutoff Parameters
- SVD cutoff: 1e-10 (default)
- Controls numerical precision
- Smaller cutoffs = higher accuracy but more memory

### System Size Scaling
- Memory: O(χ²d) where χ is bond dimension, d is physical dimension
- Time: O(Nχ³d²) for DMRG sweeps
- Measurement application: O(Nχ³)

## Advanced Usage

### Custom DMRG Parameters

```julia
function custom_ground_state(N; max_bond_dim=500, cutoff=1e-12)
    sites = siteinds("Fermion", N; conserve_qns=true)
    H = fibonacci_hamiltonian_mps(sites)
    
    ψ₀ = randomMPS(sites, ["0" for _ in 1:N])
    
    sweeps = Sweeps(15)
    setmaxdim!(sweeps, 50, 100, 200, max_bond_dim)
    setcutoff!(sweeps, cutoff)
    
    energy, ψ = dmrg(H, ψ₀, sweeps)
    return ψ, energy
end
```

### Bulk Measurement with Custom Pattern

```julia
function custom_bulk_measurement(ψ, sites, N, τ, pattern)
    current_state = copy(ψ)
    results = []
    
    for layer_sites in pattern
        layer_result = []
        for site in layer_sites
            ψ_p, prob_p = apply_measurement_mps(current_state, sites, site, τ, :p)
            # Apply measurement based on probability or deterministic choice
            # ... custom logic here
        end
        push!(results, layer_result)
    end
    
    return results
end
```

## Error Handling and Debugging

### Common Issues
1. **Bond dimension too small**: Increase `maxdim` in DMRG
2. **Numerical instability**: Decrease cutoff parameters
3. **Memory issues**: Reduce bond dimension or system size
4. **Slow convergence**: Increase number of DMRG sweeps

### Validation
- Check state normalization: `abs(inner(ψ, ψ) - 1.0) < 1e-10`
- Verify probability conservation: `prob_p + prob_m ≈ 1.0`
- Monitor bond dimensions: `linkdims(ψ)`

## Integration with Existing Code

The MPS implementation is designed to complement the existing exact diagonalization code:

```julia
# Compare exact vs MPS for small systems
N = 8
τ = 1.0

# Exact method (using existing functions)
exact_basis = Fibonacci_basis(N)
exact_gs = # ... your existing ground state calculation

# MPS method
ψ_mps, energy_mps = fibonacci_mps_ground_state(N)

# Compare results
println("Energy difference: $(energy_exact - energy_mps)")
```

## Future Extensions

1. **Time evolution**: Implement TEBD or other time evolution methods
2. **Finite temperature**: Add finite temperature MPS algorithms
3. **Disorder**: Include random disorder in the Hamiltonian
4. **Advanced measurements**: Implement more complex measurement protocols
5. **Optimization**: Further optimize for specific Fibonacci chain properties
