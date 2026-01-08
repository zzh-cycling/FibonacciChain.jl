# Examples

This section provides comprehensive examples demonstrating the key capabilities of FibonacciChain.jl.

## Example 1: Ground State Properties

Let's calculate the ground state properties of a Fibonacci anyon chain and analyze its entanglement structure.

```julia
using FibonacciChain, LinearAlgebra, Plots

# System parameters
N = 8  # Chain length

# Create model and generate basis/Hamiltonian
model = AnyonModel(FibonacciAnyon(), N; pbc=true)
basis = anyon_basis(model)
H = anyon_ham(model)

println("Hilbert space dimension: $(length(basis))")

# Find ground state
eigenvals, eigenvecs = eigen(H)
E0 = eigenvals[1]
ψ_gs = eigenvecs[:, 1]

println("Ground state energy: $E0")

# Calculate entanglement entropy profile
ee_profile = anyon_eelis(model, ψ_gs)

# Plot entanglement entropy
plot(1:N-1, ee_profile, 
     xlabel="Subsystem Size", 
     ylabel="Entanglement Entropy",
     title="EE Profile for N=$N Fibonacci Chain",
     marker=:circle)
```

## Example 2: Measurement Protocol Simulation

This example demonstrates a quantum measurement protocol with post-selection.

```julia
using FibonacciChain, Random, LinearAlgebra

N = 6
τ = 1.0

# Create model and get ground state as initial state
model = AnyonModel(FibonacciAnyon(), N; pbc=true)
basis = anyon_basis(model)
H = anyon_ham(model)
eigenvals, eigenvecs = eigen(H)
initial_state = eigenvecs[:, 1]

# Configure measurement
Random.seed!(42)
measure_config = MeasureConfig(τ=τ, t₁=1, t₂=3, mode=:Born, rng=MersenneTwister(42))

# Perform bulk evolution with Born-rule sampling
println("Performing bulk evolution...")
outcome = bulk_evolution(model, initial_state, measure_config)

# Analyze results
println("Evolution layers: $(size(outcome.samples, 1))")
println("Total free energy: $(sum(outcome.free_energys))")

# Calculate entanglement of final state
final_state = outcome.states[end]
ee_final = anyon_eelis(model, final_state)
println("Final state EE at center: $(round(ee_final[N÷2], digits=3))")
```

## Example 3: Braiding Operations and Topological Properties

Explore the effects of braiding operations on anyon states.

```julia
using FibonacciChain, LinearAlgebra

N = 6

# Create model and get ground state
model = AnyonModel(FibonacciAnyon(), N; pbc=true)
basis = anyon_basis(model)
H = anyon_ham(model)
eigenvals, eigenvecs = eigen(H)
initial_state = ComplexF64.(eigenvecs[:, 1])

# Apply braiding operation at site 3
braided_state = braidingsqmap(model, initial_state, 3)

# Calculate overlap between original and braided states
overlap = abs(dot(initial_state, braided_state))
println("Braiding overlap: $(round(overlap, digits=4))")

# Study braiding under multiple operations
state = copy(initial_state)
overlaps = Float64[]

for i in 1:10
    state = braidingsqmap(model, state, 3)
    push!(overlaps, abs(dot(initial_state, state))^2)
end

# Plot braiding evolution
plot(1:10, overlaps,
     xlabel="Braiding Steps", 
     ylabel="Fidelity |⟨ψ₀|ψ(t)⟩|²",
     title="Evolution under Repeated Braiding",
     marker=:circle)
```

## Example 4: Large System with MPS

Simulate a larger system using Matrix Product State methods.

```julia
using FibonacciChain, ITensors, ITensorMPS

N = 20
maxdim = 100
cutoff = 1e-10

println("Finding ground state for N=$N system...")

# Create model
model = AnyonModel(FibonacciAnyon(), N; pbc=true)

# Find MPS ground state
ψ_gs, E0 = anyon_mps_gst(model; 
                         sweep_times=20,
                         maxdim=maxdim, 
                         cutoff=cutoff)

println("Ground state energy density: $(E0/N)")

# Calculate entanglement entropy profile
ee_profile = anyon_eelis(model, ψ_gs)

# Fit to logarithmic scaling (CFT prediction)
using LsqFit
logChord(l, L) = @. log(sin(π * l /L))/6
    
L = length(ee_profile) + 1

mincut = 2
lm(x,p) = @. p[1] * x + p[2]
xdata = logChord([1:L-1;],L); #log.(sin.(π .* [1:L-1;] ./L))./6
fit = curve_fit(lm, xdata[mincut:L-mincut], ee_profile[mincut:L-mincut], [0.5, 0.0])
fitparam = fit.param
cent = fitparam[1]

println("Effective central charge: $(round(cent, digits=3))")

# Plot with fit
plot(L_values, ee_values, label="MPS Data", marker=:circle)
plot!(L_values, log_model(L_values, fit_result.param), 
      label="Log Fit (c=$(round(c_eff, digits=2)))", 
      linestyle=:dash)
xlabel!("Subsystem Size L")
ylabel!("Entanglement Entropy")
title!("Entanglement Scaling for N=$N")
```

## Example 5: Reference Qubit Protocol

Demonstrate the reference qubit method for measuring correlations.

```julia
using FibonacciChain, LinearAlgebra, Random

N = 8
site1 = 1  # First site for reference qubit
site2 = N ÷ 2  # Second site for reference qubit

# Create model and get ground state
model = AnyonModel(FibonacciAnyon(), N; pbc=true)
H = anyon_ham(model)
eigenvals, eigenvecs = eigen(H)
ground_state = eigenvecs[:, 1]

# Add first reference qubit at site1
println("Adding reference qubit at site $site1...")
state_with_ref1 = add_reference_qubits(model, ground_state, site1)

# Add second reference qubit at site2
println("Adding reference qubit at site $site2...")
state_with_ref2 = add_reference_qubits(model, state_with_ref1, site2)

# Configure measurement protocol
τ = 0.5
Random.seed!(42)
measure_config = MeasureConfig(τ=τ, t₁=1, t₂=4, x₁=site1, x₂=site2, mode=:Born, rng=MersenneTwister(42))

# Perform bulk evolution without reference tracking first
outcome = bulk_evolution(model, ground_state, measure_config)
println("Bulk evolution completed with $(size(outcome.samples, 1)) layers")

# Perform reference evolution for correlation measurement
ref_outcome = reference_evolution(model, outcome.states, measure_config, outcome.samples)
println("Reference evolution layers: $(size(ref_outcome.samples, 1))")

# Calculate temporal correlation
corr_temporal = temporal_correlation(model, ref_outcome.states[end])
println("Temporal correlation: $(round(corr_temporal, digits=4))")

# Calculate spatial correlation in ground state for comparison
spatial_corr = spatial_correlation(model, ground_state, site1, site2)
println("Spatial correlation (sites $site1-$site2): $(round(spatial_corr, digits=4))")
```

# Calculate spatial correlation in ground state for comparison
```julia
spatial_corr = spatial_correlation(model, ground_state, site1, site2)
println("Spatial correlation (sites $site1-$site2): $(round(spatial_corr, digits=4))")
```

## Example 6: Ising Anyon Chain

Study an Ising anyon chain with transverse field.

```julia
using FibonacciChain, LinearAlgebra

N = 10

# Create Ising model
model = AnyonModel(IsingAnyon(), N; pbc=true)
basis = anyon_basis(model)
println("Ising chain Hilbert space dimension: $(length(basis))")

# Build Hamiltonian with different field strengths
h_values = 0.0:0.2:2.0
gap_values = Float64[]

for h in h_values
    H = anyon_ham(model; J=1.0, h=h)
    eigenvals = eigvals(H)
    gap = eigenvals[2] - eigenvals[1]  # Energy gap
    push!(gap_values, gap)
end

# Find critical point (minimum gap)
min_gap, min_idx = findmin(gap_values)
h_critical = h_values[min_idx]
println("Minimum gap at h ≈ $(round(h_critical, digits=2))")

# Plot phase diagram
plot(h_values, gap_values,
     xlabel="Transverse Field h", 
     ylabel="Energy Gap Δ",
     title="Phase Diagram for N=$N Ising Chain",
     marker=:circle)
```

These examples demonstrate the versatility of FibonacciChain.jl for studying quantum many-body physics with anyons, from basic calculations to advanced protocols involving measurements and correlations.
