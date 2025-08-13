using Test
using FibonacciChain
using ITensorMPS, ITensors
using LinearAlgebra
using Random

@testset "initial_mps" begin
    N = 6
    # Test initial MPS generation
    ψ, sites = initial_mps(N)
    
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
    permuted_M_p = ITensors.permute(M_p, row_inds..., col_inds...)
    permuted_M_m = ITensors.permute(M_m, row_inds..., col_inds...)
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

    st = zeros(length(anyon_basis(N))); st[1] = 1.0
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
    state_exact = zeros(length(anyon_basis(N)))
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
    
    ψ, sites = initial_mps(N)
    
    # Define measurement sites
    measurement_sites = collect(2:2:N)
    seed=10
    # Perform sampling
    st = zeros(length(anyon_basis(N))); st[1] = 1.0
    samples_mps, samples_free_energy_mps = mps_boundary_measure(ψ, sites, measurement_sites, τ; num_samples=num_samples, rng=MersenneTwister(seed), pbc=pbc)
    sample_measured_states, samples, samples_free_energy = Boundary_measure(N, τ, st, measurement_sites,num_samples, MersenneTwister(seed))

    @test samples_mps == samples
    @test samples_free_energy_mps ≈ samples_free_energy
end

@testset "Bulk Measurements" begin
    N = 6
    pbc = true
    τ = 1.0
    D = 2  # Small number of layers for testing
    
    ψ, sites = initial_mps(N)
    
    # Define measurement sites
    measurement_sites = collect(2:2:N)
    seed=10
    # Perform sampling
    st = zeros(length(anyon_basis(N))); st[1] = 1.0
    
    # Perform bulk measurements
    bulk_states, bulk_samples, bulk_free_energy = mps_bulk_measurement(
        ψ, sites, N, τ, D; rng=MersenneTwister(seed), pbc=pbc)
    
    bulk_states_exact, bulk_samples_exact, bulk_free_energy_exact = Bulkmeasure(N, τ, st, D, MersenneTwister(seed))

    @test bulk_samples == bulk_samples_exact
    @test bulk_free_energy ≈ bulk_free_energy_exact
end

@testset "apply_measurement_layer_mps" begin
    N = 6
    pbc = true
    τ = 1.0
    
    ψ, sites = initial_mps(N)

    seed=10
    # Perform sampling
    st = zeros(length(anyon_basis(N))); st[1] = 1.0
    
    # Apply measurement to a specific layer
    measurement_layer = 2
    bulk_samples = [1, 1, 1]
    
    ψ_layer= FibonacciChain.apply_measurement_layer_mps!(N, sites, ψ, τ, bulk_samples, measurement_layer; pbc=pbc)
    
    st_exact= FibonacciChain.apply_measurement_layer!(N, st, τ, bulk_samples, measurement_layer, pbc)

    inds = [i.buf for i in anyon_basis(N)] .+1

    ψ_dense = reduce(*, ψ_layer).tensor.storage[inds]
    @test ψ_dense[ψ_dense .>0] ≈ st_exact[st_exact .>0]
end 

@testset "Generate_state_mps " begin
    N = 6
    pbc = true
    
    ψ, sites = initial_mps(N)

    seed=10
    # Perform sampling
    st = zeros(length(anyon_basis(N))); st[1] = 1.0
    
    # Generate a specific state
    measurement_sites = collect(2:2:N)  # Example measurement sites
    τ = 1.0  # Example τ value
    bulk_samples = [1 1 1; 0 0 0]
    generated_state = generate_state_mps(τ, sites, ψ, bulk_samples, true; pbc= pbc)
    generated_state_exact = generate_state(τ, st, bulk_samples, pbc, temp = true)

    inds = [i.buf for i in anyon_basis(N)] .+1

    ψ_dense = [reduce(*, ψ_layer).tensor.storage[inds] for ψ_layer in generated_state]
    ψ_dense = [sort(ψ[ψ.>0]) for ψ in ψ_dense] 
    generated_state_exact = [sort(st_exact[st_exact .>0]) for st_exact in generated_state_exact]

    @test ψ_dense ≈ generated_state_exact
end

@testset "Entanglement Entropy" begin
    N = 6
    pbc = true
    
    sites = siteinds("Qubit", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]

    ψ = randomMPS(sites, state)
    # Calculate entanglement entropy at different cuts
    EElis = anyon_eelis_mps(N, ψ)
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

function samples_generate_mps(L::Int64, τ::Float64, seed::Int64, D::Int64=5L)
    rng = MersenneTwister(seed)
    
    ψ, sites = initial_mps(L)
    
    sample_measured_states, sample, sample_free_energy = mps_bulk_measurement(ψ, sites, L, τ, D;rng=rng, pbc=true) 
    
    halfchain_EE_tlis = [ee_mps(j, div(L,2)) for j in sample_measured_states]
    final_state = sample_measured_states[end]
    final_EElis = anyon_eelis_mps(L, final_state)

    return sample, sample_free_energy, final_EElis, halfchain_EE_tlis
end

function samples_generate(L::Int64, τ::Float64, seed::Int64, D::Int64=5L)
    rng = MersenneTwister(seed)
    
    st = zeros(length(anyon_basis(L)))
    st[1] = 1.0
    
    sample_measured_states, sample, sample_free_energy = Bulkmeasure(L, τ, st, D, rng, true) 
    
    halfchain_EE_tlis = [ee(anyon_rdm(L, collect(1:div(L,2)), j)) for j in sample_measured_states]
    final_state = sample_measured_states[end]
    final_EElis = anyon_eelis(L, final_state)

    return sample, sample_free_energy, final_EElis, halfchain_EE_tlis
end

@testset "Observable" begin
    L = 6
    τ = 1.0
    D = 5
    
    # Generate samples
    seed = 42
    sample, sample_free_energy, final_EElis, halfchain_EE_tlis = samples_generate(L, τ, seed, D)
    
    sample_mps, sample_free_energy_mps, final_EElis_mps, halfchain_EE_tlis_mps = samples_generate_mps(L, τ, seed, D)

    @test sample == sample_mps
    @test sample_free_energy ≈ sample_free_energy_mps
    @test final_EElis ≈ final_EElis_mps
    @test halfchain_EE_tlis ≈ halfchain_EE_tlis_mps 
end