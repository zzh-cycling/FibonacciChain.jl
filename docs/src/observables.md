# Observable Functions

This section covers the calculation of physical observables and correlation functions.

## Entanglement Measures

### Entanglement Entropy


The von Neumann entanglement entropy is a fundamental measure of quantum entanglement:

```math
S(\rho) = -\text{Tr}(\rho \log \rho)
```

### Usage Examples

```julia
using FibonacciChain

N = 8
H = anyon_ham(N, true)
eigenvals, eigenvecs = eigen(H)
ground_state = eigenvecs[:, 1]

# Calculate entanglement entropy profile
ee_profile = anyon_eelis(N, ground_state)

# Plot the profile
using Plots
plot(1:N-1, ee_profile, xlabel="Bipartition Size", ylabel="Entanglement Entropy")
```

## Symmetry Operations

### Translation and Inversion


These functions implement discrete symmetries of the anyon chain:

- **Translation**: Shifts all anyons by one lattice site
- **Inversion**: Reflects the chain configuration

## Braiding Operations


Braiding operations are fundamental to anyonic systems, implementing the exchange statistics that define non-Abelian anyons.

### Fibonacci Braiding

For Fibonacci anyons, braiding two τ particles gives:
```math
R_{\tau\tau} = e^{4\pi i/5} \left(\begin{matrix} \phi^{-1} & \phi^{-1/2} \\ \phi^{-1/2} & -\phi^{-1} \end{matrix}\right)
```
where φ = (1+√5)/2 is the golden ratio.

## Correlation Functions

### Spatial Correlations


### Temporal Correlations


The correlation functions measure quantum correlations between different regions or times:

```julia
# Spatial correlation between sites 1 and 4
corr_spatial = spatial_correlation(N, ground_state, 1, 4)

# Temporal correlation using reference qubits
model = AnyonModel(FibonacciAnyon(), N; pbc=true)
state_with_ref = add_reference_qubits(model, ground_state, 3)
corr_temporal = temporal_correlation(N, state_with_ref)
```
