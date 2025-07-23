using FibonacciChain
using Test
using LinearAlgebra
using BitBasis
using Arpack

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

    # Test with periodic boundary conditions
    pbc = true
    output_pbc = FibonacciChain.Isingmap(T, state, N, pbc)
    @test output_pbc == (state, T(bit"111110"), -1.0, -1.0)
end

@testset "actingHamobc" begin
    output1 = FibonacciChain.actingHam(BitStr{3}, bit"000", false, measure_class=:Ising) 
    states, weights = keys(output1), values(output1)
    @test [states...]== BitStr{3}.([bit"000", bit"100", bit"010"])
    @test [weights...] ≈ [-2.0, -1.0, -1.0]

    output2 = FibonacciChain.actingHam(BitStr{3}, bit"010",false, measure_class=:Ising) 
    states, weights = keys(output2), values(output2)
    @test [states...]== BitStr{3}.([bit"000", bit"110", bit"010"])
    @test [weights...] ≈ [-1.0, -1.0, 0.0]

    output3 = FibonacciChain.actingHam(BitStr{3}, bit"001",false, measure_class=:Ising) 
    states, weights = keys(output3), values(output3)
    @test [states...]== BitStr{3}.([bit"101", bit"011", bit"001"])
    @test [weights...] ≈ [-1.0, -1.0, 2.0]

    output4 = FibonacciChain.actingHam(BitStr{3}, bit"100",false, measure_class=:Ising) 
    states, weights = keys(output4), values(output4)
    @test [states...]== BitStr{3}.([bit"100"])
    @test [weights...] ≈ [0.0]

    output = FibonacciChain.actingHam(BitStr{3}, bit"101",false, measure_class=:Ising)
    states, weights = keys(output), values(output)
    @test [states...]== BitStr{3}.([bit"101"])
    @test [weights...] ≈ [-1.0]
end

@testset "actingHampbc" begin
    ϕ = (1+√5)/2
    output1 = FibonacciChain.actingHam(BitStr{3}, bit"000", measure_class=:Ising) 
    states, weights = keys(output1), values(output1)
    @test [states...]== BitStr{3}.([bit"000",bit"100", bit"010", bit"001"])
    @test [weights...] ≈ [-3ϕ^(-1), -ϕ^(-3/2), -ϕ^(-3/2), -ϕ^(-3/2)]
    output2 = FibonacciChain.actingHam(BitStr{3}, bit"010", measure_class=:Ising) 
    states, weights = keys(output2), values(output2)
    @test [states...]== BitStr{3}.([bit"000", bit"010"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output3 = FibonacciChain.actingHam(BitStr{3}, bit"001", measure_class=:Ising) 
    states, weights = keys(output3), values(output3)
    @test [states...]== BitStr{3}.([bit"000", bit"001"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output4 = FibonacciChain.actingHam(BitStr{3}, bit"100", measure_class=:Ising) 
    states, weights = keys(output4), values(output4)
    @test [states...]== BitStr{3}.([bit"000",bit"100"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output = FibonacciChain.actingHam(BitStr{10}, bit"1000010000", measure_class=:Ising)
    states, weights = keys(output), values(output)
    @test [states...] == BitStr{10}.([bit"1000010000", bit"0000010000",bit"1010010000", bit"1000010010", bit"1000010100", bit"1000000000", bit"1001010000"])
    @test [weights...] ≈ vcat([-(4ϕ^(-1)+2ϕ^(-2))],fill(-ϕ^(-3/2),6))
end

@testset "basis.jl" begin
    # Test the Fibonacci basis creation
    fib_basis = Fibonacci_basis(5, false, measure_class=:Ising)
    @test length(fib_basis) == 32
    fib_basis = Fibonacci_basis(5, measure_class=:Ising)
    @test length(fib_basis) == 32
    # Test the Fibonacci Hamiltonian
    fib_ham = Fibonacci_Ham(5)
    @test size(fib_ham) == (32, 32)
    @test ishermitian(fib_ham)

    @test Fibonacci_Ham(3,false) == [-0.6180339887498948 0.0 -0.48586827175664565 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; -0.48586827175664565 0.0 -0.3819660112501051 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0]
    @test Fibonacci_Ham(3) == [-1.8541019662496843 -0.48586827175664565 -0.48586827175664565 -0.48586827175664565; -0.48586827175664565 -0.3819660112501051 0.0 0.0; -0.48586827175664565 0.0 -0.3819660112501051 0.0; -0.48586827175664565 0.0 0.0 -0.3819660112501051]
    # Test the reduced density matrix function
    # rdm = FibonacciChain.rdm_Fibo(fib_basis, 2)
    # @test size(rdm) == (2, 2)

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
    @test output[3] == (T(bit"010"), cstτ+coef)
    @test output[4] == (T(bit"100"), cstτ+coef)
    @test output[5] == (T(bit"101"), cstτ-coef)

    sign = 1
    coef = -sinh(τ/2) / √(2cosh(τ))
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc, measure_class=:IsingZZ)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), cstτ+coef)
    @test output[2] == (T(bit"001"), cstτ-coef)
    @test output[3] == (T(bit"010"), cstτ+coef)
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
    τ = 0.0   
    cstτ = 1/√2
    coef = 0.0      
    expected_matrix = cstτ* I(8) + coef * kron(kron(I(2), [0.0 1.0; 1.0 0.0]) , I(2))
    Mpobc = FibonacciChain.measure_matrix(T, τ, idx, 0,  measure_class=:IsingZZ)
    @test Mpobc == expected_matrix 

    coef = -sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ* I(8) + coef * kron(kron(I(2), [0.0 1.0; 1.0 0.0]) , I(2))
    Mmobc = FibonacciChain.measure_matrix(T, τ, idx, 1,  measure_class=:IsingZZ)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(8) 


end

@testset "measuremap" begin
    N = 3
    T = BitStr{N, Int}
    τ = 1.0
    idx = 2
    sign = 0
    cstτ = (exp(1)+1)/2√(exp(2)+1)
    coef = (exp(1)-1)/2√(exp(2)+1)
    ϕ = (1 + √5) / 2

    state = fill(1.0,4)
    output = measuremap(T, τ, state, idx, sign)        
    @test output == [cstτ+coef*(1-2ϕ^(-1))-2*coef*ϕ^(-3/2), cstτ+coef, cstτ+coef*(2ϕ^(-1)-1)-2*coef*ϕ^(-3/2), cstτ+coef]
    
    sign = 1
    coef = (1-exp(1))/2√(exp(2)+1)
    output = measuremap(T, τ, state, idx, sign)  
    @test output == [cstτ+coef*(1-2ϕ^(-1))-2*coef*ϕ^(-3/2), cstτ+coef, cstτ+coef*(2ϕ^(-1)-1)-2*coef*ϕ^(-3/2), cstτ+coef]

    # Test with a different state
    state = collect(1.0:4)
    output = measuremap(T, τ, state, idx, sign) 
    @test output == [cstτ+coef*(1-2ϕ^(-1))-6*coef*ϕ^(-3/2), 2(cstτ+coef), 3(cstτ+coef*(2ϕ^(-1)-1))-2*coef*ϕ^(-3/2), 4(cstτ+coef)]

    # Try with obc
    pbc = false
    state = collect(1.0:5)
    output = measuremap(T, τ, state, idx, sign, pbc)
    @test output == [cstτ+coef*(1-2ϕ^(-1))-6*coef*ϕ^(-3/2), 2(cstτ+coef), 3(cstτ+coef*(2ϕ^(-1)-1))-2*coef*ϕ^(-3/2), 4(cstτ+coef), 5(cstτ-coef)]
end

@testset "laddermeasuremap" begin
    N = 3
    T = BitStr{N, Int}
    τ = 1.0
    idx = 2


    sign = 1
    state = fill(1.0,16)
    output = laddermeasuremap(T, τ, state, idx, sign)  
    onechain_st = measuremap(T, τ, fill(1.0, 4), idx, sign)      
    @test output ≈ kron(onechain_st, onechain_st)
    
    sign = 1
    output = laddermeasuremap(T, τ, state, idx, sign)  
    onechain_st = measuremap(T, τ, fill(1.0, 4), idx, sign)
    @test output ≈ kron(onechain_st, onechain_st)
end

@testset "measurement_enumeration" begin
    N=6
    energy, states = eigen(Fibonacci_Ham(N))
    antiGS= states[:, 1]

    τ = 0.0
    measurement_sites = collect(2:2:N)
    
    final_states, trajectories, probabilities = measurement_enumeration(N, τ, antiGS, measurement_sites)

    num_final_states = length(final_states)
    @test num_final_states == 2^length(measurement_sites)

    total_prob = sum(probabilities)
    @test isapprox(total_prob, 1.0, atol=1e-6)

    @test sum(map(x->-x*log(x)/3, probabilities)) ≈ log(2) # Shannon entropy non-measurement state
end

@testset "Boundary_measure" begin
    N=6
    energy, states = Arpack.eigs(Fibonacci_Ham(N), nev=1, which=:SR)
    antiGS= states[:, 1]
    τ = 3.802
    measurement_sites = collect(2:2:N)
    
    sample_measured_states, samples, sample_free_energy = Boundary_measure(N, τ, antiGS, measurement_sites, 1000)

    num_final_states = length(sample_measured_states)
    @test num_final_states == 1000

    @test [0 for i in 2:2:N] in samples
    
end

@testset "Boundarypost_selection" begin
    N = 10
    τ = 1e3
    energy, states = Arpack.eigs(Fibonacci_Ham(N), nev=1, which=:SR)
    antiGS = states[:, 1]
    measurement_sites = collect(2:2:N)
    final_state_p, final_sequence_p, total_free_energy_p = Boundarypost_selection(N, τ, antiGS, measurement_sites, 0)
    
    @test total_free_energy_p /5 ≈ 1.1136495433981064 
end

@testset "Bulkmeasure" begin
    L = 10
    D = 2L
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0

    sample_measured_states, samples, sample_free_energy = Bulkmeasure(L, 1000.0, st, D) 
    EElis = [eelis_Fibo_state(L, state_t)[5] for state_t in sample_measured_states]
    @test size(samples) == (20, 5)
    @test EElis[1] ≈ 0.0 atol = 1e-4
    @test EElis[end] > 0.5098675501545762 
end

@testset "Bulkpost_selection" begin
    L = 10
    τ = 0.1
    D = 15L
    pbc = true
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    average_EElis=zeros(L-1)

    EE_tlis = zeros(D)
    sample_measured_states, samples, sample_free_energy = Bulkpost_selection(L, τ, st, D, 0, pbc)
    state_t = sample_measured_states[end]
    EE = eelis_Fibo_state(L, state_t)[5]
    @test samples[end] == fill(0, div(L,2))
    @test EE ≈ 0.8098675501545762 atol = 1e-4
end

@testset "generate_state" begin
    N = 10
    τ = 1e3
    energy, states = Arpack.eigs(Fibonacci_Ham(N), nev=1, which=:SR)
    antiGS = states[:, 1]
    measurement_sites = collect(2:2:N)
    
    sample_measured_states, samples, sample_free_energy = Boundary_measure(N, τ, antiGS, measurement_sites, 10)
    state = generate_state(τ, antiGS, samples[1])
    @test state ≈ sample_measured_states[1]

    st = zeros(length(Fibonacci_basis(N)))
    st[1] = 1.0

    sample_measured_states, samples, sample_free_energy = Bulkmeasure(N, τ, st, N)
    state_t = generate_state(τ, st, samples)
    statelis = generate_state(τ, st, samples, true, true)
    @test statelis ≈ sample_measured_states
    @test state_t ≈ sample_measured_states[end]

end

# Helper function to verify the distortion is working correctly
function verify_distortion(γ::Float64, original_traj::Vector{Int64}, distorted_traj::Vector{Int64})
    """
    Calculate the conditional probability P(s̃|s) for verification.
    """
    conditional_prob = 1.0
    for j in 1:length(original_traj)
        s_j = original_traj[j]
        s_tilde_j = distorted_traj[j]
        conditional_prob *= (1 + γ * s_tilde_j * s_j) / 2
    end
    return conditional_prob
end

@testset "bayes_distort" begin
    γ = 0.0
    original_traj = [1,0, 1, 0]
    probabilities = [0.25, 0.25, 0.25, 0.25]

    distorted_trajectories, distorted_probabilities = FibonacciChain.bayes_distort(γ, original_traj, probabilities)

    @test length(distorted_trajectories) == 2^4
    @test length(distorted_probabilities) == 2^4

    @test distorted_probabilities == 1/16 .* ones(16) 

    γ = 1.0
    original_traj = [1,0, 1, 0]
    
    probabilities = [0.25, 0.25, 0.25, 0.25]

    distorted_trajectories, distorted_probabilities = FibonacciChain.bayes_distort(γ, original_traj, probabilities)

    @test length(distorted_trajectories) == 2^4
    @test length(distorted_probabilities) == 2^4

    inds = findfirst(x -> x == [1, 1, 1, 1], distorted_trajectories)
    @test distorted_probabilities[inds] == (3/4)^2*(1/4)^2
end