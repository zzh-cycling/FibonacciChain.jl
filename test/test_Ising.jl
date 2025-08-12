using FibonacciChain
using Test
using LinearAlgebra
using BitBasis
using Arpack
using Random

@testset "Isingmap" begin
    N = 6
    T = BitStr{N, Int}
    state = T(0)
    idx = 2
    pbc = false
    output = FibonacciChain.Isingmap(T, state, idx)
    @test output == (state, T(bit"010000"), -1.0, -1.0)

    state = T(bit"001000")
    output = FibonacciChain.Isingmap(T, state, idx)
    @test output == (state, T(bit"011000"), 1.0, -1.0)

    state = T(bit"010100")
    output = FibonacciChain.Isingmap(T, state, idx)
    @test output == (state, T(bit"000100"), 1.0, -1.0)

    state = T(bit"111111")
    output = FibonacciChain.Isingmap(T, state, idx)
    @test output == (state, T(bit"101111"), -1.0, -1.0)

    state = T(bit"111111")
    output = FibonacciChain.Isingmap(T, state, 6, false)
    @test output == (T(bit"111110"), -1.0)

    # Test with periodic boundary conditions
    pbc = true
    output_pbc = FibonacciChain.Isingmap(T, state, N, pbc)
    @test output_pbc == (state, T(bit"111110"), -1.0, -1.0)
end

@testset "actingHamobc" begin
    N = 3
    T = BitStr{N, Int}
    output1 = FibonacciChain.actingHam(T, bit"000", false, measure_class=:IsingX, J=2.0, h=1.0) 
    states, weights = keys(output1), values(output1)
    @test [states...]== T.([bit"000", bit"100", bit"010", bit"001"])
    @test [weights...] ≈ [-4.0, -1.0, -1.0, -1.0]

    output2 = FibonacciChain.actingHam(T, bit"010",false, measure_class=:IsingX) 
    states, weights = keys(output2), values(output2)
    @test [states...]== T.([bit"000", bit"110", bit"010", bit"011"])
    @test [weights...] ≈ [-1.0, -1.0, 0.0, -1.0]

    output3 = FibonacciChain.actingHam(T, bit"001",false, measure_class=:IsingX) 
    states, weights = keys(output3), values(output3)
    @test [states...]== T.([bit"101", bit"000", bit"011", bit"001"])
    @test [weights...] ≈ [-1.0, -1.0, -1.0, 2.0]

    output4 = FibonacciChain.actingHam(T, bit"100",false, measure_class=:IsingX) 
    states, weights = keys(output4), values(output4)
    @test [states...]== T.([bit"000", bit"100", bit"110", bit"101"])
    @test [weights...] ≈ [-1.0, -2.0, -1.0, -1.0]

    output = FibonacciChain.actingHam(T, bit"101",false, measure_class=:IsingX)
    states, weights = keys(output), values(output)
    @test [states...]== T.([bit"101", bit"100", bit"111", bit"001"])
    @test [weights...] ≈ [2.0, -1.0, -1.0, -1.0]
end

@testset "actingHampbc" begin
    N = 3
    T = BitStr{N, Int}
    output1 = FibonacciChain.actingHam(T, bit"000", measure_class=:IsingX, J=2.0, h=1.0) 
    states, weights = keys(output1), values(output1)
    @test [states...]== T.([bit"000",bit"100", bit"010", bit"001"])
    @test [weights...] ≈ [-6.0, -1.0, -1.0, -1.0]

    output2 = FibonacciChain.actingHam(T, bit"010", measure_class=:IsingX) 
    states, weights = keys(output2), values(output2)
    @test [states...]== T.([bit"000", bit"110", bit"010", bit"011"])
    @test [weights...] ≈ [-1.0, -1.0, -1.0, -1.0]
    
    output3 = FibonacciChain.actingHam(T, bit"001", measure_class=:IsingX) 
    states, weights = keys(output3), values(output3)
    @test [states...]== T.([bit"101", bit"000", bit"011", bit"001"])
    @test [weights...] ≈ [-1.0, -1.0, -1.0, 3.0]

    output4 = FibonacciChain.actingHam(T, bit"100", measure_class=:IsingX) 
    states, weights = keys(output4), values(output4)
    @test [states...]== T.([bit"000",bit"100", bit"110", bit"101"])
    @test [weights...] ≈ [-1.0, -1.0, -1.0, -1.0]

    output = FibonacciChain.actingHam(BitStr{10}, bit"1000010000", measure_class=:IsingX)
    states, weights = keys(output), values(output)
    @test [states...] == BitStr{10}.([16, 560, 656, 532, 529, 784, 512, 536, 528, 530, 592])
    @test [weights...] ≈ vcat(-ones(8),[-4.0, -1.0, -1.0])
end

@testset "basis.jl" begin
    # Test the Fibonacci basis creation
    fib_basis = Fibonacci_basis(5, false, measure_class=:IsingX)
    @test length(fib_basis) == 32
    fib_basis = Fibonacci_basis(5, measure_class=:IsingX)
    @test length(fib_basis) == 32
    # Test the Fibonacci Hamiltonian
    fib_ham = Fibonacci_Ham(5, measure_class=:IsingX)
    @test size(fib_ham) == (32, 32)
    @test ishermitian(fib_ham)

    @test Fibonacci_Ham(3,false, measure_class=:IsingX) == [-2.0 -1.0 -1.0 0.0 -1.0 0.0 0.0 0.0; -1.0 2.0 0.0 -1.0 0.0 -1.0 0.0 0.0; -1.0 0.0 0.0 -1.0 0.0 0.0 -1.0 0.0; 0.0 -1.0 -1.0 0.0 0.0 0.0 0.0 -1.0; -1.0 0.0 0.0 0.0 -2.0 -1.0 -1.0 0.0; 0.0 -1.0 0.0 0.0 -1.0 2.0 0.0 -1.0; 0.0 0.0 -1.0 0.0 -1.0 0.0 0.0 -1.0; 0.0 0.0 0.0 -1.0 0.0 -1.0 -1.0 0.0]

    H = Fibonacci_Ham(3, measure_class=:IsingX)
    @test H == [-3.0 -1.0 -1.0 0.0 -1.0 0.0 0.0 0.0; -1.0 3.0 0.0 -1.0 0.0 -1.0 0.0 0.0; -1.0 0.0 -1.0 -1.0 0.0 0.0 -1.0 0.0; 0.0 -1.0 -1.0 1.0 0.0 0.0 0.0 -1.0; -1.0 0.0 0.0 0.0 -1.0 -1.0 -1.0 0.0; 0.0 -1.0 0.0 0.0 -1.0 1.0 0.0 -1.0; 0.0 0.0 -1.0 0.0 -1.0 0.0 1.0 -1.0; 0.0 0.0 0.0 -1.0 0.0 -1.0 -1.0 -1.0]

    gs = eigvecs(fib_ham)[:,1]
    # Test the reduced density matrix function
    rdm = rdm_Fibo(5, collect(1:1), gs, measure_class=:IsingX)
    @test size(rdm) == (2, 2)
    @test rdm ≈ [0.8191079332998358 0.37755708971817986; 0.37755708971817986 0.18089206670016397]
end

@testset "measure_basismap_IsingX" begin
    N = 3
    T = BitStr{N, Int}
    idx = 2
    pbc = false
    τ =0.0
    cstτ = cosh(τ/2) / √(2cosh(τ))
    sign = 0
    basis0 = [T(0b000), T(0b001), T(0b010), T(0b100), T(0b101)]

    output = measure_basismap.(T, τ, basis0, idx, sign, pbc, measure_class=:IsingX)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, 0.0)
    @test output[2] == (T(bit"001"), T(bit"011"), cstτ, 0.0)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, 0.0)
    @test output[4] == (T(bit"100"), T(bit"110"), cstτ, 0.0)
    @test output[5] == (T(bit"101"), T(bit"111"), cstτ, 0.0)

    sign = 1
    output2 = measure_basismap.(T, τ, basis0, idx, sign, pbc, measure_class=:IsingX)
    @test output2 == output

    τ = 1.0
    sign = 0
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc, measure_class=:IsingX)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, coef)
    @test output[2] == (T(bit"001"), T(bit"011"), cstτ, coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, coef)
    @test output[4] == (T(bit"100"), T(bit"110"), cstτ, coef)
    @test output[5] == (T(bit"101"), T(bit"111"), cstτ, coef)

    sign = 1
    coef = -sinh(τ/2) / √(2cosh(τ))
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc, measure_class=:IsingX)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, coef)
    @test output[2] == (T(bit"001"), T(bit"011"), cstτ, coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, coef)
    @test output[4] == (T(bit"100"), T(bit"110"), cstτ, coef)
    @test output[5] == (T(bit"101"), T(bit"111"), cstτ, coef)

    τ = 1e3
    sign = 0
    cstτ = 1/2
    coef = 1/2
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc, measure_class=:IsingX)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, coef)
    @test output[2] == (T(bit"001"), T(bit"011"), cstτ, coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, coef)
    @test output[4] == (T(bit"100"), T(bit"110"), cstτ, coef)
    @test output[5] == (T(bit"101"), T(bit"111"), cstτ, coef)


    idx = 3
    output = measure_basismap.(T, τ, basis0, idx, sign, true, measure_class=:IsingX)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"001"), cstτ, coef)
    @test output[2] == (T(bit"001"), T(bit"000"), cstτ, coef)
    @test output[3] == (T(bit"010"), T(bit"011"), cstτ, coef)
    @test output[4] == (T(bit"100"), T(bit"101"), cstτ, coef)
    @test output[5] == (T(bit"101"), T(bit"100"), cstτ, coef)

end

@testset "measure_basismap_IsingZZ" begin
    N = 3
    T = BitStr{N, Int}
    idx = 2
    pbc = false
    τ =0.0
    cstτ = cosh(τ/2) / √(2cosh(τ))
    sign = 0
    basis0 = [T(0b000), T(0b001), T(0b010), T(0b100), T(0b101)]


    output = measure_basismap.(T, τ, basis0, idx, sign, pbc, measure_class=:IsingZZ)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), cstτ)
    @test output[2] == (T(bit"001"), cstτ)
    @test output[3] == (T(bit"010"), cstτ)
    @test output[4] == (T(bit"100"), cstτ)
    @test output[5] == (T(bit"101"), cstτ)

    sign = 1
    output2 = measure_basismap.(T, τ, basis0, idx, sign, pbc, measure_class=:IsingZZ)
    @test output2 == output

    τ = 1.0
    sign = 0
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc, measure_class=:IsingZZ)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), cstτ+coef)
    @test output[2] == (T(bit"001"), cstτ-coef)
    @test output[3] == (T(bit"010"), cstτ-coef)
    @test output[4] == (T(bit"100"), cstτ+coef)
    @test output[5] == (T(bit"101"), cstτ-coef)

    sign = 1
    coef = -sinh(τ/2) / √(2cosh(τ))
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc, measure_class=:IsingZZ)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), cstτ+coef)
    @test output[2] == (T(bit"001"), cstτ-coef)
    @test output[3] == (T(bit"010"), cstτ-coef)
    @test output[4] == (T(bit"100"), cstτ+coef)
    @test output[5] == (T(bit"101"), cstτ-coef)

    idx=3
    τ = 1e3
    sign = 0
    cstτ = 1/2
    coef = 1/2
    output = measure_basismap.(T, τ, basis0, idx, sign, measure_class=:IsingZZ)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), 1.0)
    @test output[2] == (T(bit"001"), 0.0)
    @test output[3] == (T(bit"010"), 1.0)
    @test output[4] == (T(bit"100"), 0.0)
    @test output[5] == (T(bit"101"), 1.0)


    sign = 1
    output = measure_basismap.(T, τ, basis0, idx, sign, measure_class=:IsingZZ) # pbc is true by default
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), 0.0)
    @test output[2] == (T(bit"001"), 1.0)
    @test output[3] == (T(bit"010"), 0.0)
    @test output[4] == (T(bit"100"), 1.0)
    @test output[5] == (T(bit"101"), 0.0)

end

@testset "measure_matrix" begin
    N = 3
    T = BitStr{N, Int}
    
    ⊗(a,b) = kron(a, b)
    τ = 1.0
    idx = 2
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))

    σx = [0.0 1.0; 1.0 0.0]
    σz = [1.0 0.0; 0.0 -1.0]
    # measuring X
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σx ⊗ I(2)
    Mpobc = FibonacciChain.measure_matrix(T, τ, idx, 0, false, measure_class=:IsingX)
    @test Mpobc == expected_matrix 

    coef = -sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σx ⊗ I(2)
    Mmobc = FibonacciChain.measure_matrix(T, τ, idx, 1, false, measure_class=:IsingX)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(8) 

    # measuring ZZ
    coef = sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σz ⊗ σz
    Mpobc = FibonacciChain.measure_matrix(T, τ, idx, 0, false, measure_class=:IsingZZ)
    @test Mpobc == expected_matrix 

    coef = -sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σz ⊗ σz
    Mmobc = FibonacciChain.measure_matrix(T, τ, idx, 1, false, measure_class=:IsingZZ)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(8) 

    # Test with a different τ, idx
    idx = 3
    τ = 1000.0   
    cstτ = 0.5
    coef = 0.5      
    expected_matrix = cstτ* I(8) + coef * σz ⊗ I(2) ⊗ σz 
    Mpobc = FibonacciChain.measure_matrix(T, τ, idx, 0, measure_class=:IsingZZ)
    @test Mpobc == expected_matrix 

    coef = -0.5    
    expected_matrix = cstτ* I(8) + coef * σz ⊗ I(2) ⊗ σz 
    Mmobc = FibonacciChain.measure_matrix(T, τ, idx, 1, measure_class=:IsingZZ)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(8) 


end

@testset "measuremap_IsingX" begin
    N = 3
    T = BitStr{N, Int}
    τ = 1.0
    idx = 2
    sign = 0
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))

    state = fill(1.0,2^N)
    output = measuremap(T, τ, state, idx, sign, measure_class=:IsingX)        
    @test output == (cstτ+coef) .* ones(2^N)
    
    sign = 1
    coef = -sinh(τ/2) / √(2cosh(τ))
    output = measuremap(T, τ, state, idx, sign, measure_class=:IsingX)  
    @test output == (cstτ+coef) .* ones(2^N)
end

@testset "measuremap_IsingZZ" begin
    N = 3
    T = BitStr{N, Int}
    τ = 1.0
    idx = 2
    sign = 0
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))

    state = fill(1.0,2^N)
    output = measuremap(T, τ, state, idx, sign, measure_class=:IsingZZ)        
    @test output == [cstτ+coef, cstτ-coef,cstτ-coef, cstτ+coef, cstτ+coef, cstτ-coef,cstτ-coef, cstτ+coef]
    
    sign = 1
    coef = -sinh(τ/2) / √(2cosh(τ))
    output = measuremap(T, τ, state, idx, sign, measure_class=:IsingZZ)  
    @test output == [cstτ+coef, cstτ-coef,cstτ-coef, cstτ+coef, cstτ+coef, cstτ-coef,cstτ-coef, cstτ+coef]


    # Try with idx=3 pbc
    state = ones(2^N)
    output = measuremap(T, 1000.0, state, idx, sign, measure_class=:IsingZZ)
    @test output == [0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0]
end


@testset "measurement_enumeration" begin
    N=6
    st = zeros(length(Fibonacci_basis(N, measure_class=:IsingX)))
    st[1] = 1.0 
    τ = 0.0
    measurement_sites = collect(2:2:N)

    final_states, trajectories, probabilities = measurement_enumeration(N, τ, st, measurement_sites, measure_class=:IsingX)

    num_final_states = length(final_states)
    @test num_final_states == 2^length(measurement_sites)

    # all final states should be equally probable
    total_prob = sum(probabilities)
    @test isapprox(total_prob, 1.0, atol=1e-6)
    @test probabilities ≈ 1/8 .* ones(2^length(measurement_sites))
    @test sum(map(x->-x*log(x)/length(measurement_sites), probabilities)) ≈ log(2) # Shannon entropy non-measurement state

    final_states, trajectories, probabilities = measurement_enumeration(N, τ, st, measurement_sites, measure_class=:IsingZZ)

    num_final_states = length(final_states)
    @test num_final_states == 2^length(measurement_sites)

    total_prob = sum(probabilities)
    @test isapprox(total_prob, 1.0, atol=1e-6)
    @test probabilities ≈ 1/8 .* ones(2^length(measurement_sites))
    @test sum(map(x->-x*log(x)/length(measurement_sites), probabilities)) ≈ log(2) # Shannon entropy non-measurement state
end

@testset "Boundary_measure" begin
    N=6
    st = zeros(length(Fibonacci_basis(N, measure_class=:IsingX)))
    st[1] = 1.0 
    τ = 3.802
    measurement_sites = collect(2:2:N)

    sample_measured_states, samples, sample_free_energy = Boundary_measure(N, τ, st, measurement_sites, 500, measure_class=:IsingX)

    num_final_states = length(sample_measured_states)
    @test num_final_states == 500

    @test [0 for i in 2:2:N] in samples
    
end

@testset "Boundarypost_selection" begin
    N = 6
    τ = 1e3
    st = zeros(length(Fibonacci_basis(N, measure_class=:IsingX)))
    st[1] = 1.0 
 
    measurement_sites = collect(2:2:N)
    final_state_p, final_sequence_p, total_free_energy_p = Boundarypost_selection(N, τ, st, measurement_sites, 0, measure_class=:IsingX)

    # all final states should be equally probable, will give Nlog(2) free energy
    @test total_free_energy_p /length(measurement_sites) ≈ log(2) atol = 1e-6

    final_state_p, final_sequence_p, total_free_energy_p = Boundarypost_selection(N, τ, final_state_p, measurement_sites, 0, measure_class=:IsingZZ)

    @test total_free_energy_p /length(measurement_sites) ≈ log(2) atol = 1e-6
    @test final_state_p[1] == 1.0
end

@testset "Bulkmeasure" begin
    L = 6
    D = 2L
    τ = 1e3
    st = zeros(length(Fibonacci_basis(L, measure_class=:IsingX)))
    st[1] = 1.0

    sample_measured_states, samples, sample_free_energy = Bulkmeasure(L, 1000.0, st, D, MersenneTwister(100), measure_class=:IsingX) 
    EElis = [eelis_Fibo_state(L, state_t, measure_class=:IsingX)[div(L,2)] for state_t in sample_measured_states]
    @test size(samples) == (D, L)
    # Each layer will erase previous info.
    @test EElis ≈ [i % 2 == 1 ? 0.0 : log(2) for i in 1:D] atol = 1e-6
end

@testset "Bulkpost_selection" begin
    L = 6
    τ = 1000.0
    D = 10L
    pbc = true
    st=zeros(length(Fibonacci_basis(L, measure_class=:IsingX)))
    st[1] = 1.0
    average_EElis=zeros(L-1)

    EE_tlis = zeros(D)
    sample_measured_states, samples, sample_free_energy = Bulkpost_selection(L, τ, st, D, 0, pbc, measure_class=:IsingX)
    state_t = sample_measured_states[end]
    EE = eelis_Fibo_state(L, state_t, measure_class=:IsingX)
    @test samples[end] == fill(0, L)
    @test EE ≈ log(2)*ones(L-1) atol = 1e-4
end

@testset "apply_measurement_layer" begin
    N = 6
    τ = 1e3
    st = zeros(length(Fibonacci_basis(N, measure_class=:IsingX)))
    st[1] = 1.0

    sample_measured_states, samples, sample_free_energy = Bulkmeasure(N, τ, st, N, MersenneTwister(100), measure_class=:IsingX)
    state_t = sample_measured_states[end]

    new_state = FibonacciChain.apply_measurement_layer!(N, st, τ, samples[1,:], 1, measure_class=:IsingX)
    @test new_state ≈ sample_measured_states[1]
end

# @testset "generate_state" begin
#     N = 6
#     τ = 1e3
#     st=zeros(length(Fibonacci_basis(N, measure_class=:IsingX)))
#     st[1] = 1.0

#     sample_measured_states, samples, sample_free_energy = Boundary_measure(N, τ, st, collect(1:N), 10, measure_class=:IsingX)
#     state = generate_state(τ, st, samples[1], measure_class=:IsingX)
#     @test state ≈ sample_measured_states[1]

#     sample_measured_states, samples, sample_free_energy = Bulkmeasure(N, τ, st, N, MersenneTwister(100), measure_class=:IsingX)
#     state_t = generate_state(τ, st, samples, measure_class=:IsingX)
#     statelis = generate_state(τ, st, samples, true, temp = true, measure_class=:IsingX)
#     @test statelis ≈ sample_measured_states
#     @test state_t ≈ sample_measured_states[end]

# end