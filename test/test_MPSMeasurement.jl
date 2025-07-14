using Test
using FibonacciChain
using ITensorMPS, ITensors
using LinearAlgebra

@testset "initial_mps" begin
    N = 6
    # Test initial MPS generation
    ψ = initial_mps(N)
    
    @test ψ isa MPS
    @test length(ψ) == N
    @test all(linkdims(ψ) .== 1)  # Initial state should have bond dimension 1
    
    # Test normalization, set to be 1
    @test abs(inner(ψ, ψ) - 1.0) < 1e-10
end

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
    sites = siteinds("Qubit", N)
    τ = 1.0
    
    idx=2
    # Test measurement operator creation
    M_p = measurement_operator_mps(sites, idx, τ, 0)
    M_m = measurement_operator_mps(sites, idx, τ, 1)

    s_im1 = sites[idx-1] # site 1
    s_i   = sites[idx]   # site 2
    s_ip1 = sites[idx+1] # site 3
    row_inds = (prime(s_im1), prime(s_i), prime(s_ip1))
    col_inds = (s_im1, s_i, s_ip1)
    permuted_M_p = permute(M_p, row_inds..., col_inds...)
    permuted_M_m = permute(M_m, row_inds..., col_inds...)
    M_pmatrix = reshape(permuted_M_p.tensor.storage, 8,8)
    M_mmatrix = reshape(permuted_M_m.tensor.storage, 8,8)

    @test M_pmatrix^2 + M_mmatrix^2 ≈ I(8)
    # @test  M_pmatrix == 
    #       [1.0 0.0 0.0 0.0; 
    #        0.0 1.0 0.0 0.0; 
    #        0.0 0.0 exp(-τ) 0.0; 
    #        0.0 0.0 0.0 exp(-τ)]
    # @test M_mmatrix == 
    #       [exp(-τ) 0.0 0.0 0.0; 
    #        0.0 exp(-τ) 0.0 0.0; 
    #        0.0 0.0 1.0 0.0; 
    #        0.0 0.0 0.0 1.0]


    # Test invalid inputs
    @test_throws AssertionError measurement_operator_mps(sites, 0, τ, 0)
    @test_throws AssertionError measurement_operator_mps(sites, N+1, τ, 0)
    @test_throws AssertionError measurement_operator_mps(sites, 2, τ, 2)
end

@testset "Single Measurement Application" begin
    N = 4
    pbc = true
    τ = 1.0
    
    sites = siteinds("Qubit", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]
    
    # Create MPS from product state
    ψ0 = randomMPS(sites, state)
    
    # Apply measurements
    ψ_p, prob_p = apply_measurement_mps(ψ0, sites, 1, τ, 0; pbc=pbc)
    ψ_m, prob_m = apply_measurement_mps(ψ0, sites, 1, τ, 1; pbc=pbc)

    st = zeros(length(Fibonacci_basis(N))); st[1] = 1.0
    state_after_p = measuremap(N, τ, st, 1, 0)
    p = state_after_p'*state_after_p

    @test prob_p ≈ p
    @test prob_m ≈ 1 - p  # Should be orthogonal to ψ_p

    if prob_p > 1e-12
        @test abs(inner(ψ_p, ψ_p) - 1.0) < 1e-10
    end
    if prob_m > 1e-12
        @test abs(inner(ψ_m, ψ_m) - 1.0) < 1e-10
    end
end


@testset "Measurement Enumeration" begin
    N = 4
    pbc = true
    τ = 1.0
    
    sites = siteinds("Qubit", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]
    state_exact = zeros(length(Fibonacci_basis(N)))
    state_exact[1] = 1.0  # Vacuum state
    # Use minimal measurement sites for enumeration
    measurement_sites = collect(2:2:N)
    
    ψ = randomMPS(sites, state)
    # Enumerate trajectories
    final_states, trajectories, probabilities = mps_measurement_enumeration(
        ψ, sites, measurement_sites, τ; pbc=pbc)

    # final_states_vector = map(x-> reduce(*, x).tensor.storage, final_states) # Noting its elements arranging is different from the final_states_exact, they choose different definition?
    
    final_states_exact, trajectories_exact, probabilities_exact = measurement_enumeration(N, τ, state_exact, measurement_sites)

    @test trajectories == trajectories_exact
    @test probabilities ≈ probabilities_exact
 
    # Check probability normalization
    @test abs(sum(probabilities) - 1.0) < 1e-10
    
    # Check that all probabilities are positive
    @test all(p -> p > 0, probabilities)
end

@testset "Boundary Measurements" begin
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
        @test all(s -> s in [0, 1], samples[i])
        @test weights[i] > 0
        @test weights[i] <= 1
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
        @test all(s -> s in [0, 1], bulk_samples[layer])
        @test bulk_weights[layer] >= 0
        @test isfinite(bulk_weights[layer])
    end
end

@testset "Generate_state_mps " begin
    
end

@testset "Entanglement Entropy" begin
    N = 6
    pbc = true
    
    sites = siteinds("Qubit", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]

    ψ = randomMPS(sites, state)
    # Calculate entanglement entropy at different cuts
    EElis = eelis_Fibo_mps(N, ψ)
    @test all(EElis .>= 0)  # Entanglement entropy should be non-negative
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
    @test_throws AssertionError apply_measurement_mps(ψ, sites, 0, 1.0, 0)
    @test_throws AssertionError apply_measurement_mps(ψ, sites, N+1, 1.0, 0)
    @test_throws AssertionError apply_measurement_mps(ψ, sites, 2, 1.0, 2)
    
    # Test edge cases for τ
    ψ_large_τ, prob_large_τ = apply_measurement_mps(ψ, sites, 2, 1e3, 0)
    @test isfinite(prob_large_τ)
    @test 0 <= prob_large_τ <= 1
    
    ψ_small_τ, prob_small_τ = apply_measurement_mps(ψ, sites, 2, 1e-3, 0)
    @test isfinite(prob_small_τ)
    @test 0 <= prob_small_τ <= 1
end
