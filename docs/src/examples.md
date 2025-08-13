# Examples

This section provides comprehensive examples demonstrating the key capabilities of FibonacciChain.jl.

## Example 1: Ground State Properties

Let's calculate the ground state properties of a Fibonacci anyon chain and analyze its entanglement structure.

```julia
using FibonacciChain, LinearAlgebra, Plots

# System parameters
N = 8  # Chain length
pbc = true  # Periodic boundary conditions

# Generate basis and Hamiltonian
basis = anyon_basis(N, pbc)
H = anyon_ham(N, pbc)

println("Hilbert space dimension: $(length(basis))")

# Find ground state
eigenvals, eigenvecs = eigen(H)
E0 = eigenvals[1]
ψ_gs = eigenvecs[:, 1]

println("Ground state energy: $E0")

# Calculate entanglement entropy profile
ee_profile = anyon_eelis(N, ψ_gs, pbc)

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
using FibonacciChain, Random

N = 6
τ = 1.0
num_samples = 1000

# Initialize random state  
Random.seed!(42)
initial_state = normalize!(rand(ComplexF64, length(anyon_basis(N))))

# Define measurement sites (even sites for Fibonacci anyons)
measurement_sites = [2, 4, 6]

# Perform boundary measurements
println("Performing $num_samples measurement samples...")
measured_states, sequences, free_energies = Boundary_measure(
    N, τ, initial_state, measurement_sites, num_samples
)

# Analyze results
using Statistics
mean_energy = mean(free_energies)
std_energy = std(free_energies)

println("Average free energy: $(round(mean_energy, digits=4)) ± $(round(std_energy, digits=4))")

# Post-select specific measurement outcome (all +1)
target_sequence = [0, 0, 0]  # All +1 measurements
selected_states = [state for (state, seq) in zip(measured_states, sequences) if seq == target_sequence]

if !isempty(selected_states)
    println("Found $(length(selected_states)) trajectories with target sequence")
    avg_selected_state = mean(selected_states)
    
    # Calculate entanglement of post-selected state
    ee_selected = anyon_eelis(N, avg_selected_state)
    println("Post-selected EE at center: $(round(ee_selected[N÷2], digits=3))")
end
```

## Example 3: Braiding Operations and Topological Properties

Explore the effects of braiding operations on anyon states.

```julia
using FibonacciChain

N = 6
initial_state = normalize!(ones(ComplexF64, length(anyon_basis(N))))

# Apply braiding operation at site 3
braided_state = braidingsqmap(N, initial_state, 3)

# Calculate overlap between original and braided states
overlap = abs(dot(initial_state, braided_state))
println("Braiding overlap: $(round(overlap, digits=4))")

# Study braiding under multiple operations
state = copy(initial_state)
overlaps = Float64[]

for i in 1:10
    state = braidingsqmap(N, state, 3)
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
using FibonacciChain, ITensors

N = 40
maxdim = 100
cutoff = 1e-10

println("Finding ground state for N=$N system...")

# Find MPS ground state
ψ_gs, E0 = fibonacci_mps_ground_state(N, 
                                      pbc=true,
                                      sweep_times=20,
                                      maxdim=maxdim, 
                                      cutoff=cutoff)

println("Ground state energy density: $(E0/N)")

# Calculate entanglement entropy profile
ee_profile = anyon_eelis_mps(ψ_gs)

# Analyze scaling behavior
L_values = 1:N÷2
ee_values = ee_profile[L_values]

# Fit to logarithmic scaling (CFT prediction)
using LsqFit
@. log_model(x, p) = p[1] * log(x) + p[2]
fit_result = curve_fit(log_model, Float64.(L_values), ee_values, [1.0, 0.0])
c_eff = 6 * fit_result.param[1]  # Effective central charge

println("Effective central charge: $(round(c_eff, digits=3))")

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
using FibonacciChain

N = 8
site = 4  # Site to probe

# Initialize ground state
H = anyon_ham(N, true)
eigenvals, eigenvecs = eigen(H)
ground_state = eigenvecs[:, 1]

# Add reference qubit at site 4
println("Adding reference qubit at site $site...")
state_with_ref = add_reference_qubits!(N, ground_state, site)

# Perform measurement protocol to create temporal correlation
τ = 0.5
sample_sequence = rand([0, 1], 5, N÷2)  # 5 time steps, N/2 measurement sites per step

# Generate forward evolution
forward_states = reference_generate_state(N, τ, state_with_ref, sample_sequence, temp=true)

# Calculate temporal correlation between time slices 1 and 4
time_corr = reference_evolution(τ, forward_states, sample_sequence, site, 1, 4)

println("Temporal correlation: $(round(real(time_corr), digits=4))")

# Compare with spatial correlation in ground state
spatial_corr = spatial_correlation(N, ground_state, 1, 5)
println("Spatial correlation (sites 1-5): $(round(spatial_corr, digits=4))")
```

## Example 6: Phase Transition Study

Investigate phase transitions by varying Hamiltonian parameters.

```julia
using FibonacciChain

N = 10
h_values = 0.0:0.1:2.0  # Magnetic field strength
gap_values = Float64[]

for h in h_values
    # Construct Hamiltonian with external field
    H = anyon_ham(N, true, h_field=h)  # Assuming h_field parameter exists
    
    eigenvals = eigvals(H)
    gap = eigenvals[2] - eigenvals[1]  # Energy gap
    push!(gap_values, gap)
end

# Plot phase diagram
plot(h_values, gap_values,
     xlabel="Field Strength h", 
     ylabel="Energy Gap Δ",
     title="Phase Diagram for N=$N",
     marker=:circle)
```

These examples demonstrate the versatility of FibonacciChain.jl for studying quantum many-body physics with anyons, from basic calculations to advanced protocols involving measurements and correlations.
