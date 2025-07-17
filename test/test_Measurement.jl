using FibonacciChain
using Test
using LinearAlgebra
using BitBasis
using Arpack

@testset "measure_basismap" begin
    N = 3
    ϕ = (1 + √5) / 2
    T = BitStr{N, Int}
    state = T(0b000)
    idx = 2
    pbc = false
    τ =0.0
    cstτ = (exp(τ)+1)/2√(exp(2τ)+1)
    sign = 0
    basis0 = [T(0b000), T(0b001), T(0b010), T(0b100), T(0b101)]

    output = measure_basismap.(T, τ, basis0, idx, sign, pbc)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, 0.0)
    @test output[2] == (T(bit"001"), cstτ)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, 0.0)
    @test output[4] == (T(bit"100"), cstτ)
    @test output[5] == (T(bit"101"), cstτ)

    τ = 0.0
    sign = 1
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, 0.0)
    @test output[2] == (T(bit"001"), cstτ)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, 0.0)
    @test output[4] == (T(bit"100"), cstτ)
    @test output[5] == (T(bit"101"), cstτ)

    τ = 1.0
    sign = 0
    cstτ = (exp(τ)+1)/2√(exp(2τ)+1)
    coef = (exp(τ)-1)/2√(exp(2τ)+1)
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), cstτ+coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2))
    @test output[4] == (T(bit"100"), cstτ+coef)
    @test output[5] == (T(bit"101"), cstτ-coef)

    τ = 1.0
    sign = 1
    coef = (1-exp(τ))/2√(exp(2τ)+1)
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), cstτ+coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2))
    @test output[4] == (T(bit"100"), cstτ+coef)
    @test output[5] == (T(bit"101"), cstτ-coef)

    τ = 1e3
    sign = 0
    cstτ = 1/2
    coef = 1/2
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), cstτ+coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2))
    @test output[4] == (T(bit"100"), cstτ+coef)
    @test output[5] == (T(bit"101"), cstτ-coef)

    idx = 1
    output = measure_basismap.(T, τ, basis0, idx, sign) # pbc is true by default
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"100"), cstτ+coef*(1-2ϕ^(-1)),-2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), cstτ+coef)
    @test output[3] == (T(bit"010"), cstτ+coef)
    @test output[4] == (T(bit"100"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1),-2*coef*ϕ^(-3/2))
    @test output[5] === nothing

    idx = 3
    output = measure_basismap.(T, τ, basis0, idx, sign) # pbc is true by default
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"001"), cstτ+coef*(1-2ϕ^(-1)),-2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1),-2*coef*ϕ^(-3/2))
    @test output[3] == (T(bit"010"), cstτ+coef)
    @test output[4] == (T(bit"100"), cstτ+coef)
    @test output[5] === nothing
end

@testset "measure_matrix" begin
    N = 3
    T = BitStr{N, Int}
    ϕ = (1 + √5) / 2


    τ = 1.0
    idx = 2
    cstτ = (exp(1)+1)/2√(exp(2)+1)
    coef = (exp(1)-1)/2√(exp(2)+1)
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0 0.0;
        0.0 cstτ+coef 0.0 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0 0.0;
        0.0 0.0 0.0 cstτ+coef 0.0;
        0.0 0.0 0.0 0.0 cstτ-coef]
    Mpobc = FibonacciChain.measure_matrix(T, τ, idx, 0, false)
    @test Mpobc == expected_matrix 

    coef = (1-exp(1))/2√(exp(2)+1)
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0 0.0;
        0.0 cstτ+coef 0.0 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0 0.0;
        0.0 0.0 0.0 cstτ+coef 0.0;
        0.0 0.0 0.0 0.0 cstτ-coef
    ]
    Mmobc = FibonacciChain.measure_matrix(T, τ, idx, 1, false)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(5) 

    coef = (exp(1)-1)/2√(exp(2)+1)
    expected_matrix = [        
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0;
        0.0 cstτ+coef 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0;
        0.0 0.0 0.0 cstτ+coef
]
    Mppbc = FibonacciChain.measure_matrix(T, τ, idx, 0) # pbc
    @test Mppbc == expected_matrix 
    coef = (1-exp(1))/2√(exp(2)+1)
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0;
        0.0 cstτ+coef 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0;
        0.0 0.0 0.0 cstτ+coef
        ]
    Mmpbc = FibonacciChain.measure_matrix(T, τ, idx, 1) # pbc
    @test Mmpbc == expected_matrix    
    @test Mppbc^2+Mmpbc^2 ≈ I(4)

    # Test with a different τ
    τ = 0.0   
    cstτ = 1/√2
    coef = 0.0      
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0 0.0;
        0.0 cstτ+coef 0.0 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0 0.0;
        0.0 0.0 0.0 cstτ+coef 0.0;
        0.0 0.0 0.0 0.0 cstτ-coef
    ]
    Mpobc = FibonacciChain.measure_matrix(T, τ, idx, 0, false)
    @test Mpobc == expected_matrix 
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0 0.0;
        0.0 cstτ+coef 0.0 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0 0.0;
        0.0 0.0 0.0 cstτ+coef 0.0;
        0.0 0.0 0.0 0.0 cstτ-coef
    ]
    Mmobc = FibonacciChain.measure_matrix(T, τ, idx, 1, false)
    @test Mmobc == expected_matrix 
    @test Mpobc^2+Mmobc^2 ≈ I(5) 

    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0;
        0.0 cstτ+coef 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0;
        0.0 0.0 0.0 cstτ+coef
    ]
    Mppbc = FibonacciChain.measure_matrix(T, τ, idx, 0) # pbc
    @test Mppbc == expected_matrix 
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0;
        0.0 cstτ+coef 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0;
        0.0 0.0 0.0 cstτ+coef
    ]
    Mmpbc = FibonacciChain.measure_matrix(T, τ, idx, 1) # pbc
    @test Mmpbc == expected_matrix 
    @test Mppbc^2+Mmpbc^2 ≈ I(4) 


end

@testset "measure_matrix" begin
    # Test with a different idx， must be pbc.
    N = 3
    T = BitStr{N, Int}
    ϕ = (1 + √5) / 2
    τ = 1.0
    cstτ = (exp(1)+1)/2√(exp(2)+1)
    coef = (exp(1)-1)/2√(exp(2)+1)
    idx = 1

    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 0.0 -2*coef*ϕ^(-3/2);
        0.0 cstτ+coef 0.0 0.0;
        0.0 0.0 cstτ+coef 0.0 ;
        -2*coef*ϕ^(-3/2) 0.0 0.0 cstτ+coef*(2ϕ^(-1)-1)
    ]
    Mppbc = FibonacciChain.measure_matrix(T, τ, idx, 0) # pbc
    @test Mppbc == expected_matrix 
    coef = (1-exp(1))/2√(exp(2)+1)
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 0.0 -2*coef*ϕ^(-3/2);
        0.0 cstτ+coef 0.0 0.0;
        0.0 0.0 cstτ+coef 0.0 ;
        -2*coef*ϕ^(-3/2) 0.0 0.0 cstτ+coef*(2ϕ^(-1)-1)
        ]
    Mmpbc = FibonacciChain.measure_matrix(T, τ, idx, 1) # pbc
    @test Mmpbc == expected_matrix    
    @test Mppbc^2+Mmpbc^2 ≈ I(4)
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