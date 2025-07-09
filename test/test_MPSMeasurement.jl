using Test
using FibonacciChain
using ITensorMPS
using LinearAlgebra


@testset "MPS Ground State Generation" begin
    N = 6
    pbc = true
    
    # Test ground state generation
    ψ, energy = fibonacci_mps_ground_state(N; pbc=pbc)
    
    @test ψ isa MPS
    @test length(ψ) == N
    @test energy isa Real
    @test energy < 0  # Ground state should have negative energy
    
    # Test normalization
    @test abs(inner(ψ, ψ) - 1.0) < 1e-10
    
    # Test that bond dimensions are reasonable
    max_bond_dim = maximum(linkdims(ψ))
    @test max_bond_dim >= 1
    @test max_bond_dim <= 200  # Should not exceed our DMRG settings
end

@testset "Measurement Operators" begin
    N = 4
    sites = siteinds("Fermion", N)
    τ = 1.0
    
    # Test measurement operator creation
    M_p = measurement_operator_mps(sites, 2, τ, :p)
    M_m = measurement_operator_mps(sites, 2, τ, :m)
    
    @test M_p isa MPO
    @test M_m isa MPO
    @test length(M_p) == N
    @test length(M_m) == N
    
    # Test invalid inputs
    @test_throws AssertionError measurement_operator_mps(sites, 0, τ, :p)
    @test_throws AssertionError measurement_operator_mps(sites, N+1, τ, :p)
    @test_throws AssertionError measurement_operator_mps(sites, 2, τ, :invalid)
end

@testset "Single Measurement Application" begin
    N = 4
    pbc = true
    τ = 1.0
    
    # Generate ground state
    ψ, _ = fibonacci_mps_ground_state(N; pbc=pbc)
    sites = siteinds(ψ)
    
    # Apply measurements
    ψ_p, prob_p = apply_measurement_mps(ψ, sites, 2, τ, :p; pbc=pbc)
    ψ_m, prob_m = apply_measurement_mps(ψ, sites, 2, τ, :m; pbc=pbc)
    
    # Test that results are valid
    @test ψ_p isa MPS
    @test ψ_m isa MPS
    @test 0 <= prob_p <= 1
    @test 0 <= prob_m <= 1
    @test abs(prob_p + prob_m - 1.0) < 1e-10  # Probabilities should sum to 1
    
    # Test normalization of resulting states
    if prob_p > 1e-12
        @test abs(inner(ψ_p, ψ_p) - 1.0) < 1e-10
    end
    if prob_m > 1e-12
        @test abs(inner(ψ_m, ψ_m) - 1.0) < 1e-10
    end
end

@testset "Sampling Measurements" begin
    N = 6
    pbc = true
    τ = 1.0
    num_samples = 10
    
    # Generate ground state
    ψ, _ = fibonacci_mps_ground_state(N; pbc=pbc)
    sites = siteinds(ψ)
    
    # Define measurement sites
    measurement_sites = [2, 4]
    
    # Perform sampling
    samples, weights = mps_sampling(ψ, sites, measurement_sites, τ; 
                                   num_samples=num_samples, pbc=pbc)
    
    # Test outputs
    @test length(samples) == num_samples
    @test length(weights) == num_samples
    
    # Check each sample
    for i in 1:num_samples
        @test length(samples[i]) == length(measurement_sites)
        @test all(s -> s in [:p, :m], samples[i])
        @test weights[i] > 0
        @test weights[i] <= 1
    end
end

@testset "Measurement Enumeration" begin
    N = 4
    pbc = true
    τ = 1.0
    
    # Generate ground state
    ψ, _ = fibonacci_mps_ground_state(N; pbc=pbc)
    sites = siteinds(ψ)
    
    # Use minimal measurement sites for enumeration
    measurement_sites = [2]
    
    # Enumerate trajectories
    final_states, trajectories, probabilities = mps_measurement_enumeration(
        ψ, sites, measurement_sites, τ; pbc=pbc)
    
    # Should have exactly 2 trajectories for 1 measurement site
    @test length(final_states) == 2
    @test length(trajectories) == 2
    @test length(probabilities) == 2
    
    # Check trajectories
    trajectory_outcomes = [traj[1] for traj in trajectories]
    @test :p in trajectory_outcomes
    @test :m in trajectory_outcomes
    
    # Check probability normalization
    @test abs(sum(probabilities) - 1.0) < 1e-10
    
    # Check that all probabilities are positive
    @test all(p -> p > 0, probabilities)
end

@testset "Entanglement Entropy" begin
    N = 6
    pbc = true
    
    # Generate ground state
    ψ, _ = fibonacci_mps_ground_state(N; pbc=pbc)
    
    # Calculate entanglement entropy at different cuts
    for b in 2:(N-1)
        ee = calculate_entanglement_entropy_mps(ψ, b)
        @test ee >= 0  # Entanglement entropy should be non-negative
        @test isfinite(ee)
    end
end

@testset "Bulk Measurements" begin
    N = 6
    pbc = true
    τ = 1.0
    D = 2  # Small number of layers for testing
    
    # Generate ground state
    ψ, _ = fibonacci_mps_ground_state(N; pbc=pbc)
    sites = siteinds(ψ)
    
    # Perform bulk measurements
    bulk_states, bulk_samples, bulk_weights = mps_bulk_measurement(
        ψ, sites, N, τ, D; pbc=pbc)
    
    # Test outputs
    @test length(bulk_states) == D
    @test length(bulk_samples) == D
    @test length(bulk_weights) == D
    
    # Check each layer
    for layer in 1:D
        @test bulk_states[layer] isa MPS
        @test length(bulk_samples[layer]) == div(N, 2)
        @test all(s -> s in [:p, :m], bulk_samples[layer])
        @test bulk_weights[layer] >= 0
        @test isfinite(bulk_weights[layer])
    end
end

@testset "Parameter Validation" begin
    N = 4
    
    # Test invalid system sizes
    @test_throws BoundsError fibonacci_mps_ground_state(0)
    @test_throws BoundsError fibonacci_mps_ground_state(-1)
    
    # Generate valid state for further tests
    ψ, _ = fibonacci_mps_ground_state(N)
    sites = siteinds(ψ)
    
    # Test invalid measurement parameters
    @test_throws AssertionError apply_measurement_mps(ψ, sites, 0, 1.0, :p)
    @test_throws AssertionError apply_measurement_mps(ψ, sites, N+1, 1.0, :p)
    @test_throws AssertionError apply_measurement_mps(ψ, sites, 2, 1.0, :invalid)
    
    # Test edge cases for τ
    ψ_large_τ, prob_large_τ = apply_measurement_mps(ψ, sites, 2, 1e3, :p)
    @test isfinite(prob_large_τ)
    @test 0 <= prob_large_τ <= 1
    
    ψ_small_τ, prob_small_τ = apply_measurement_mps(ψ, sites, 2, 1e-3, :p)
    @test isfinite(prob_small_τ)
    @test 0 <= prob_small_τ <= 1
end
