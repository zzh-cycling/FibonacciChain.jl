using FibonacciChain
using ITensorMPS

"""
Example usage of MPS-based Fibonacci chain measurements
"""

function mps_example()
    println("=== MPS-based Fibonacci Chain Example ===")
    
    # System parameters
    N = 10  # Number of sites
    τ = 1.0  # Measurement parameter
    pbc = true  # Periodic boundary conditions
    
    println("1. Generating MPS ground state for N=$N sites...")
    
    # Generate ground state using MPS
    ψ, energy = fibonacci_mps_ground_state(N; pbc=pbc)
    println("Ground state energy: $energy")
    println("MPS bond dimensions: $(linkdims(ψ))")
    
    # Get site indices for measurements
    sites = siteinds(ψ)
    
    println("\n2. Calculating entanglement entropy...")
    
    # Calculate entanglement entropy at different cuts
    for b in 2:(N-1)
        ee = calculate_entanglement_entropy_mps(ψ, b)
        println("Entanglement entropy at cut $b: $ee")
    end
    
    println("\n3. Performing single measurement...")
    
    # Apply a single measurement at site 2
    measurement_site = 2
    ψ_measured_p, prob_p = apply_measurement_mps(ψ, sites, measurement_site, τ, :p; pbc=pbc)
    ψ_measured_m, prob_m = apply_measurement_mps(ψ, sites, measurement_site, τ, :m; pbc=pbc)
    
    println("Measurement at site $measurement_site:")
    println("  P(+) = $prob_p")
    println("  P(-) = $prob_m")
    println("  Total probability = $(prob_p + prob_m)")
    
    println("\n4. Sampling measurements...")
    
    # Define measurement sites (boundary measurements)
    measurement_sites = [2, 4, 6, 8]
    num_samples = 100
    
    samples, weights = mps_sampling(ψ, sites, measurement_sites, τ; 
                                   num_samples=num_samples, pbc=pbc)
    
    println("Generated $num_samples measurement samples")
    println("Sample measurement outcomes (first 10):")
    for i in 1:min(10, length(samples))
        println("  Sample $i: $(samples[i]) (weight: $(weights[i]))")
    end
    
    # Calculate average weight
    avg_weight = sum(weights) / length(weights)
    println("Average sample weight: $avg_weight")
    
    println("\n5. Enumerating all measurement trajectories...")
    
    # Use fewer sites for enumeration to keep it manageable
    enum_sites = [2, 4]
    
    final_states, trajectories, probabilities = mps_measurement_enumeration(
        ψ, sites, enum_sites, τ; pbc=pbc)
    
    println("Total number of trajectories: $(length(trajectories))")
    println("All measurement trajectories:")
    for i in 1:length(trajectories)
        println("  $(trajectories[i]) → probability: $(probabilities[i])")
    end
    
    # Verify total probability
    total_prob = sum(probabilities)
    println("Total probability: $total_prob")
    
    println("\n6. Bulk measurements...")
    
    # Perform bulk measurements with alternating pattern
    D = 3  # Number of layers
    
    bulk_states, bulk_samples, bulk_weights = mps_bulk_measurement(
        ψ, sites, N, τ, D; pbc=pbc)
    
    println("Performed $D layers of bulk measurements")
    for layer in 1:D
        println("Layer $layer:")
        println("  Measurement pattern: $(bulk_samples[layer])")
        println("  Log probability: $(bulk_weights[layer])")
        println("  Final state bond dims: $(linkdims(bulk_states[layer]))")
    end
    
    println("\n=== Example completed ===")
end

function compare_exact_vs_mps()
    println("\n=== Comparing Exact vs MPS Methods ===")
    
    N = 8  # Small system for comparison
    τ = 1.0
    pbc = true
    
    println("System size: N=$N")
    
    # Generate MPS ground state
    ψ_mps, energy_mps = fibonacci_mps_ground_state(N; pbc=pbc)
    sites = siteinds(ψ_mps)
    
    # Compare with exact diagonalization (if feasible)
    try
        # This would use your existing exact methods
        println("MPS ground state energy: $energy_mps")
        println("MPS bond dimensions: $(linkdims(ψ_mps))")
        
        # Single measurement comparison
        measurement_site = 2
        ψ_p, prob_p_mps = apply_measurement_mps(ψ_mps, sites, measurement_site, τ, :p; pbc=pbc)
        ψ_m, prob_m_mps = apply_measurement_mps(ψ_mps, sites, measurement_site, τ, :m; pbc=pbc)
        
        println("\nMPS measurement probabilities at site $measurement_site:")
        println("  P(+) = $prob_p_mps")
        println("  P(-) = $prob_m_mps")
        
        # Calculate entanglement entropy
        println("\nMPS entanglement entropy:")
        for b in 2:(N-1)
            ee = calculate_entanglement_entropy_mps(ψ_mps, b)
            println("  Cut $b: $ee")
        end
        
    catch e
        println("Note: Exact comparison not available: $e")
    end
end

# Run examples if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    mps_example()
    compare_exact_vs_mps()
end
