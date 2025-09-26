# Matrix Product State Methods

This section covers the efficient Matrix Product State (MPS) implementations for simulating larger anyon systems using the ITensors.jl library.

## Ground State Calculations


MPS methods allow simulation of much larger systems (N ~ 50-100) compared to exact diagonalization:

```julia
using FibonacciChain

# Find ground state for N=50 Fibonacci anyon chain
N = 50
ψ, E0 = anyon_mps_gst(N, pbc=true, maxdim=100, cutoff=1e-10)
```

## Hamiltonian Construction


The MPS Hamiltonian is constructed as a Matrix Product Operator (MPO) with three-body interactions:

```julia
# Create sites and Hamiltonian MPO
sites = siteinds("Qubit", N)
H = anyon_ham(sites, pbc=true)
```

## MPS Measurements

### Measurement Operators


### MPS-based Protocols



## Entanglement Entropy



MPS naturally provides access to entanglement structure:

```julia
# Calculate entanglement entropy between first L and last N-L sites
L = N ÷ 2
S = ee_mps(ψ, L)

# Get full entanglement profile
ee_profile = anyon_eelis(N, ψ)
```

## State Generation and Evolution



## Advantages of MPS Methods

### Computational Efficiency
- **Memory**: Scales as O(χ²N) instead of O(2ᴺ)
- **Time**: DMRG scales as O(χ³N) for ground states
- **Bond dimension**: χ ~ 50-200 typically sufficient

### Physical Insight
- Direct access to entanglement structure
- Natural handling of gapped phases
- Efficient for 1D systems with area law

## Usage Example: Large System Simulation

```julia
using FibonacciChain, ITensors

N = 80  # Much larger than exact diagonalization allows
maxdim = 200
cutoff = 1e-12

# Find ground state
ψ_gs, E0 = anyon_mps_gst(N, 
                                      pbc=true,
                                      sweep_times=30, 
                                      maxdim=maxdim,
                                      cutoff=cutoff)

# Calculate observables
ee_profile = anyon_eelis(N, ψ_gs)
central_ee = ee_profile[N÷2]

println("Ground state energy: $E0")
println("Central entanglement entropy: $central_ee")

# Apply measurement and evolve
τ = 0.5
measurement_site = N÷2
ψ_measured, prob = apply_measurement_mps(ψ_gs, sites, measurement_site, τ, 0)

println("Measurement probability: $prob")
```

## Performance Tips

1. **Bond dimension**: Start with small χ and increase gradually
2. **Convergence**: Monitor energy convergence during DMRG sweeps  
3. **Memory**: Use appropriate number precision (Float64 vs ComplexF64)
4. **Parallelization**: ITensors supports threading for large calculations
