using FibonacciChain
using Test
using LinearAlgebra
using BitBasis
using Arpack
using Random 
using LsqFit

@testset "measure_basismap" begin
    N = 3
    ϕ = (1 + √5) / 2
    T = BitStr{N, Int}
    state = T(0b000)
    idx = 2
    pbc = false
    τ =0.0
    cstτ = (exp(τ)+1)/2√(exp(2τ)+1)
    sign = false
    basis0 = [T(0b000), T(0b001), T(0b010), T(0b100), T(0b101)]
    model = AnyonModel(FibonacciAnyon(), N, pbc = pbc)

    output = measure_basismap.(Ref(model), τ, basis0, idx, sign) # s=0
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, 0.0)
    @test output[2] == (T(bit"001"), T(bit"001"), cstτ, 0.0)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, 0.0)
    @test output[4] == (T(bit"100"), T(bit"100"), cstτ, 0.0)
    @test output[5] == (T(bit"101"), T(bit"101"), cstτ, 0.0)

    sign = false
    output2 = measure_basismap.(Ref(model), τ, basis0, idx, sign) # s=0
    @test output2 == output
    
    τ = 1.0
    sign = false
    cstτ = (exp(τ)+1)/2√(exp(2τ)+1)
    coef = (exp(τ)-1)/2√(exp(2τ)+1)
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign) # s=0
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), T(bit"001"), cstτ+coef, 0.0)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2))
    @test output[4] == (T(bit"100"), T(bit"100"), cstτ+coef, 0.0)
    @test output[5] == (T(bit"101"), T(bit"101"), cstτ-coef, 0.0)

    sign = true
    coef = (1-exp(τ))/2√(exp(2τ)+1)
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign) # s=1
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), T(bit"001"), cstτ+coef, 0.0)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2))
    @test output[4] == (T(bit"100"), T(bit"100"), cstτ+coef, 0.0)
    @test output[5] == (T(bit"101"), T(bit"101"), cstτ-coef, 0.0)

    τ = 1e3
    sign = false
    cstτ = 1/2
    coef = 1/2
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign) # s=0
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), T(bit"001"), cstτ+coef, 0.0)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2))
    @test output[4] == (T(bit"100"), T(bit"100"), cstτ+coef, 0.0)
    @test output[5] == (T(bit"101"), T(bit"101"), cstτ-coef, 0.0)

    idx = 1
    model = AnyonModel(FibonacciAnyon(), N) # pbc is true by default
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign) # pbc is true by default, s=0
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"100"), cstτ+coef*(1-2ϕ^(-1)),-2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), T(bit"001"), cstτ+coef, 0.0)
    @test output[3] == (T(bit"010"), T(bit"010"), cstτ+coef, 0.0)
    @test output[4] == (T(bit"100"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1),-2*coef*ϕ^(-3/2))
    @test output[5] === nothing

    idx = 3
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign) # pbc is true by default, s=0
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"001"), cstτ+coef*(1-2ϕ^(-1)),-2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1),-2*coef*ϕ^(-3/2))
    @test output[3] == (T(bit"010"), T(bit"010"), cstτ+coef, 0.0)
    @test output[4] == (T(bit"100"), T(bit"100"), cstτ+coef, 0.0)
    @test output[5] === nothing

    model = AnyonModel(FibonacciAnyon(), N, pbc = true, measure_operator=:reset)
    output = measure_basismap.(Ref(model), 1000.0, basis0, idx, sign)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"000"), 1.0, 0.0)
    @test output[2] == (T(bit"001"), T(bit"001"), 0.0, 0.0)
    @test output[3] == (T(bit"010"), T(bit"010"), 1.0, 0.0)
    @test output[4] == (T(bit"100"), T(bit"100"), 1.0, 0.0)
    @test output[5] == (T(bit"101"), T(bit"101"), 0.0, 0.0) # Noting such basis didn't show in Fibonacci basis
end

@testset "measure_matrix" begin
    N = 3
    ϕ = (1 + √5) / 2
    model = AnyonModel(FibonacciAnyon(), N, pbc = false)

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
    Mpobc = FibonacciChain.measure_matrix(model, τ, idx, false) # s=0
    @test Mpobc == expected_matrix

    coef = (1-exp(1))/2√(exp(2)+1)
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0 0.0;
        0.0 cstτ+coef 0.0 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0 0.0;
        0.0 0.0 0.0 cstτ+coef 0.0;
        0.0 0.0 0.0 0.0 cstτ-coef
    ]
    Mmobc = FibonacciChain.measure_matrix(model, τ, idx, true) # s=1
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(5) 


    model = AnyonModel(FibonacciAnyon(), N) # pbc is true by default
    coef = (exp(1)-1)/2√(exp(2)+1)
    expected_matrix = [        
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0;
        0.0 cstτ+coef 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0;
        0.0 0.0 0.0 cstτ+coef
]
    Mppbc = FibonacciChain.measure_matrix(model, τ, idx, false) # pbc, s=0
    @test Mppbc == expected_matrix 
    coef = (1-exp(1))/2√(exp(2)+1)
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0;
        0.0 cstτ+coef 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0;
        0.0 0.0 0.0 cstτ+coef
    ]
    Mmpbc = FibonacciChain.measure_matrix(model, τ, idx, true) # pbc, s=1
    @test Mmpbc == expected_matrix    
    @test Mppbc^2+Mmpbc^2 ≈ I(4)

    # Test with a different τ
    τ = 0.0   
    cstτ = 1/√2
    coef = 0.0      
    model = AnyonModel(FibonacciAnyon(), N, pbc = false)
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0 0.0;
        0.0 cstτ+coef 0.0 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0 0.0;
        0.0 0.0 0.0 cstτ+coef 0.0;
        0.0 0.0 0.0 0.0 cstτ-coef
    ]
    Mpobc = FibonacciChain.measure_matrix(model, τ, idx, false) # s=0
    @test Mpobc == expected_matrix 
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0 0.0;
        0.0 cstτ+coef 0.0 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0 0.0;
        0.0 0.0 0.0 cstτ+coef 0.0;
        0.0 0.0 0.0 0.0 cstτ-coef
    ]
    Mmobc = FibonacciChain.measure_matrix(model, τ, idx, true) # s=1
    @test Mmobc == expected_matrix 
    @test Mpobc^2+Mmobc^2 ≈ I(5) 

    model = AnyonModel(FibonacciAnyon(), N) # pbc is true by default
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0;
        0.0 cstτ+coef 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0;
        0.0 0.0 0.0 cstτ+coef
    ]
    Mppbc = FibonacciChain.measure_matrix(model, τ, idx, false) # pbc, s=0
    @test Mppbc == expected_matrix 
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 -2*coef*ϕ^(-3/2) 0.0;
        0.0 cstτ+coef 0.0 0.0;
        -2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0;
        0.0 0.0 0.0 cstτ+coef
    ]
    Mmpbc = FibonacciChain.measure_matrix(model, τ, idx, true) # pbc, s=1
    @test Mmpbc == expected_matrix
    @test Mppbc^2+Mmpbc^2 ≈ I(4)

    model = AnyonModel(IsingAnyon(), N, pbc = true, measure_operator=:reset)
    Mppbc = FibonacciChain.measure_matrix(model, 1000.0, idx, false) 
    @test diag(Mppbc) ≈ [1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
    Mmpbc = FibonacciChain.measure_matrix(model, 1000.0, idx, true) # pbc
    @test diag(Mmpbc) ≈ [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
    # false = 0 → +1 eigenvalue, 1+ Id ⊕ Z ⊕ Id; true = 1 → -1 eigenvalue, 1 - Id ⊕ Z ⊕ Id
end

@testset "measure_matrix" begin
    # Test with a different idx， must be pbc.
    N = 3
    ϕ = (1 + √5) / 2
    τ = 1.0
    cstτ = (exp(1)+1)/2√(exp(2)+1)
    coef = (exp(1)-1)/2√(exp(2)+1)
    idx = 1


    model = AnyonModel(FibonacciAnyon(), N) # pbc is true by default
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 0.0 -2*coef*ϕ^(-3/2);
        0.0 cstτ+coef 0.0 0.0;
        0.0 0.0 cstτ+coef 0.0 ;
        -2*coef*ϕ^(-3/2) 0.0 0.0 cstτ+coef*(2ϕ^(-1)-1)
    ]
    Mppbc = FibonacciChain.measure_matrix(model, τ, idx, false) # pbc, s=0
    @test Mppbc == expected_matrix 
    coef = (1-exp(1))/2√(exp(2)+1)
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 0.0 -2*coef*ϕ^(-3/2);
        0.0 cstτ+coef 0.0 0.0;
        0.0 0.0 cstτ+coef 0.0 ;
        -2*coef*ϕ^(-3/2) 0.0 0.0 cstτ+coef*(2ϕ^(-1)-1)
    ]
    Mmpbc = FibonacciChain.measure_matrix(model, τ, idx, true) # pbc, s=1
    @test Mmpbc == expected_matrix    
    @test Mppbc^2+Mmpbc^2 ≈ I(4)
end

@testset "Temperley Lieb algebra" begin
    N = 8
    τ = 1000.0
    ϕ = (1 + √5) / 2
    model = AnyonModel(FibonacciAnyon(), N)
    Xlis = ϕ .* [FibonacciChain.measure_matrix(model, τ, idx, true) for idx in 1:N] # s=1

    # X_i ^2 = d X_i
    @test all(Xlis[i] * Xlis[i] ≈ ϕ .* Xlis[i] for i in 1:N)
    # X_i * X_{i+1} * X_i = X_i
    @test all(Xlis[i] * Xlis[i+1] * Xlis[i] ≈ Xlis[i] for i in 1:N-1)
    # X_i * X_{i-1} * X_i = X_i
    @test all(Xlis[i] * Xlis[i-1] * Xlis[i] ≈ Xlis[i] for i in 2:N)
    # [X_i, X_{j}] = 0, |i-j|>=2
    @test all(Xlis[i] * Xlis[j] ≈ Xlis[j] * Xlis[i] for i in 1:N for j in i+2:N-1)
end

@testset "measuremap" begin
    N = 3
    model = AnyonModel(FibonacciAnyon(), N)
    τ = 1.0
    idx = 2
    sign = false # s=0
    cstτ = (exp(1)+1)/2√(exp(2)+1)
    coef = (exp(1)-1)/2√(exp(2)+1)
    ϕ = (1 + √5) / 2

    state = fill(1.0,4)
    output = measuremap(model, τ, state, idx, sign)        
    @test output == [cstτ+coef*(1-2ϕ^(-1))-2*coef*ϕ^(-3/2), cstτ+coef, cstτ+coef*(2ϕ^(-1)-1)-2*coef*ϕ^(-3/2), cstτ+coef]
    
    sign = true # s=1
    coef = (1-exp(1))/2√(exp(2)+1)
    output = measuremap(model, τ, state, idx, sign)  
    @test output == [cstτ+coef*(1-2ϕ^(-1))-2*coef*ϕ^(-3/2), cstτ+coef, cstτ+coef*(2ϕ^(-1)-1)-2*coef*ϕ^(-3/2), cstτ+coef]

    # Test with a different state
    state = collect(1.0:4)
    output = measuremap(model, τ, state, idx, sign) 
    @test output == [cstτ+coef*(1-2ϕ^(-1))-6*coef*ϕ^(-3/2), 2(cstτ+coef), 3(cstτ+coef*(2ϕ^(-1)-1))-2*coef*ϕ^(-3/2), 4(cstτ+coef)]

    # Try with obc
    model = AnyonModel(FibonacciAnyon(), N, pbc = false)
    state = collect(1.0:5)
    output = measuremap(model, τ, state, idx, sign)
    @test output == [cstτ+coef*(1-2ϕ^(-1))-6*coef*ϕ^(-3/2), 2(cstτ+coef), 3(cstτ+coef*(2ϕ^(-1)-1))-2*coef*ϕ^(-3/2), 4(cstτ+coef), 5(cstτ-coef)]
end

@testset "laddermeasuremap" begin
    N = 3
    model = AnyonModel(FibonacciAnyon(), N)
    τ = 1.0
    idx = 2


    sign = true
    state = fill(1.0,16)
    output = laddermeasuremap(model, τ, state, idx, sign)  
    onechain_st = measuremap(model, τ, fill(1.0, 4), idx, sign)      
    @test output ≈ kron(onechain_st, onechain_st)

    sign = true # s=1
    output = laddermeasuremap(model, τ, state, idx, sign)
    onechain_st = measuremap(model, τ, fill(1.0, 4), idx, sign)
    @test output ≈ kron(onechain_st, onechain_st)
end

@testset "measurement_enumeration" begin
    N=6
    model = AnyonModel(FibonacciAnyon(), N)
    energy, states = eigen(anyon_ham(model))
    antiGS= states[:, 1]

    τ = 0.0
    measurement_sites, measure_operator = FibonacciChain._obtain_measurement_config(model, 2)

    final_states, trajectories, probabilities = measurement_enumeration(model, τ, antiGS, measurement_sites)

    num_final_states = length(final_states)
    @test num_final_states == 2^length(measurement_sites)

    total_prob = sum(probabilities)
    @test isapprox(total_prob, 1.0, atol=1e-6)

    @test sum(map(x->-x*log(x)/3, probabilities)) ≈ log(2) # Shannon entropy non-measurement state
end

@testset "_sample_layer" begin
    model = AnyonModel(FibonacciAnyon(), 6)
    τ = 1000.0
    state = zeros(length(anyon_basis(model))); state[1] = 1.0
    rng = MersenneTwister(42)
    measure_outcome = FibonacciChain._sample_layer(model, τ, state, rng = rng)
    @test measure_outcome.sample == [1, 0, 1]
    @test measure_outcome.free_energy ≈ 1.9248473002384137 atol=1e-6
end

@testset "_apply_measurement_layer" begin
    model = AnyonModel(FibonacciAnyon(), 6)
    τ = 1000.0
    state = zeros(length(anyon_basis(model))); state[1] = 1.0
    sample = BitVector(zeros(3))
    measure_outcome = FibonacciChain._apply_measurement_layer(model, τ, state, sample, layer_idx = 1)
    @test measure_outcome.free_energy ≈ 2.887270950357621 atol=1e-6
end

@testset "_born_measure" begin
    model = AnyonModel(FibonacciAnyon(), 6)
    t = 10
    measure_config = MeasureConfig(τ=1000.0, t₂=t, rng=MersenneTwister(42), mode=:Born)
    state = zeros(length(anyon_basis(model))); state[1] = 1.0

    measure_outcome = FibonacciChain._born_measure(model, state, measure_config)
    @test size(measure_outcome.samples) == (20, 3)
    @test measure_outcome.free_energys[end] ≈ 1.5009765892377303 atol=1e-6
end

@testset "_sample_measure" begin
    N = 6
    model = AnyonModel(IsingAnyon(), N, measure_operator=:X)
    t = 10
    measure_config = MeasureConfig(τ=1000.0, t₂=t, rng=MersenneTwister(42), mode=:sample)
    state = zeros(length(anyon_basis(model))); state[1] = 1.0
    samples = BitMatrix(zeros(2t, N))

    measure_outcome = FibonacciChain._sample_measure(model, state , samples, measure_config)
    @test measure_outcome.states[end][[1,64]] ≈ 1/√2 .* ones(2)
    @test measure_outcome.free_energys[end] ≈ 5log(2) atol=1e-6
end

@testset "boundary_evolution, bulk_evolution" begin
    N = 10
    τ = 1e3
    model = AnyonModel(FibonacciAnyon(), N)
    energy, states = Arpack.eigs(anyon_ham(model), nev=1, which=:SR)
    antiGS = states[:, 1]
    measure_config = MeasureConfig(τ=τ, t₂=1, rng=MersenneTwister(42), mode=:Born, enable_τ_eff = false)

    measure_outcome = boundary_evolution(model, antiGS, measure_config)
    sample_measured_state, sample, sample_free_energy = measure_outcome.state, measure_outcome.sample, measure_outcome.free_energy

    measure_config_boundary_generate = MeasureConfig(τ=τ, t₂=1, mode=:sample, enable_τ_eff = false)
    measure_outcome_boundary =  boundary_evolution(model, antiGS, measure_config_boundary_generate, sample)

    @test measure_outcome_boundary.state ≈ sample_measured_state
    @test measure_outcome_boundary.free_energy ≈ sample_free_energy atol=1e-6

    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0

    measure_config_bulk_evolution = MeasureConfig(τ=τ, t₂=N, rng=MersenneTwister(42), mode=:Born)
    measure_outcome3 = bulk_evolution(model, st, measure_config_bulk_evolution)
    sample_measured_states, sample_bulk, sample_free_energy = measure_outcome3.states, measure_outcome3.samples, measure_outcome3.free_energys

    measure_config_bulk_generate = MeasureConfig(τ=τ, t₂=N, mode=:sample)
    measure_outcome4 = bulk_evolution(model, st, measure_config_bulk_generate, sample_bulk)
    statelis, Fs = measure_outcome4.states, measure_outcome4.free_energys
    state = statelis[end]
    @test state ≈ sample_measured_states[end]
end

@testset "boundary_Born" begin
    N=6
    model = AnyonModel(FibonacciAnyon(), N)
    energy, states = Arpack.eigs(anyon_ham(model), nev=1, which=:SR)
    antiGS= states[:, 1]
    τ = 3.802
    
    measure_config = MeasureConfig(t₂=1, τ=τ, rng=MersenneTwister(100), mode=:Born, enable_τ_eff = false)
    # all samples are Matrix now
    measure_outcome = boundary_evolution(model, antiGS, measure_config)
    measured_state, sample, sample_free_energy = measure_outcome.state, measure_outcome.sample, measure_outcome.free_energy

    @test BitVector(ones(3)) == sample
    @test sample_free_energy ≈ 0.38031571763999733 atol=1e-6
end

@testset "boundary_post_selection" begin
    N = 10
    τ = 1e3
    measure_config = MeasureConfig(τ=τ, t₂=1, mode=:sample, enable_τ_eff = false)
    model = AnyonModel(FibonacciAnyon(), N)

    energy, states = Arpack.eigs(anyon_ham(model), nev=1, which=:SR)
    antiGS = states[:, 1]
    sample = BitVector(zeros(div(N,2)))
    measure_outcome = boundary_evolution(model, antiGS, measure_config, sample) # layer_idx default to 1

    @test measure_outcome.free_energy /5 ≈ 1.1136495433981064 atol = 1e-7
end

@testset "bulk_Born" begin
    L = 10
    t = 2L
    measure_config = MeasureConfig(τ=1000.0, t₂=t, rng=MersenneTwister(2), mode=:Born)
    model = AnyonModel(FibonacciAnyon(), L)
    
    st=zeros(length(anyon_basis(model)))
    st[1] = 1.0

    measure_outcome = bulk_evolution(model, st, measure_config)
    sample_measured_states, samples, sample_free_energy = measure_outcome.states, measure_outcome.samples, measure_outcome.free_energys
    EElis = [anyon_eelis(model, state_t)[5] for state_t in sample_measured_states]
    @test size(samples) == (2t, 5)
    @test EElis[1] ≈ 0.6895721925700435 atol = 1e-4
    @test EElis[end] > 0.7
    @test sample_free_energy[end] ≈  3.371812546192422 atol=1e-6
end

@testset "bulk_post_selection" begin
    L = 10
    τ = 0.1
    D = 15L
    model = AnyonModel(FibonacciAnyon(), L)
    measure_config = MeasureConfig(τ=τ, t₂=div(D,2), mode=:sample)
    st=zeros(length(anyon_basis(model)))
    st[1] = 1.0
    average_EElis=zeros(L-1)

    EE_tlis = zeros(D)
    samples = BitMatrix(zeros(Int8, D, div(L,2)))
    measure_outcome = bulk_evolution(model, st, measure_config, samples)
    sample_measured_states, samples, sample_free_energy = measure_outcome.states, measure_outcome.samples, measure_outcome.free_energys
    state_t = sample_measured_states[end]
    EE = anyon_eelis(model, state_t)[5]
    @test samples[end,:] == fill(0, div(L,2))
    @test EE ≈ 0.8098675501545762 atol = 1e-4
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

function fitCCEntEntScal(
    SvN_list::Vector{Float64};
    err::Vector{Float64}=0.0SvN_list,
    mincut::Int=1,
    pbc::Bool=false)

    # log of chord length / 6 for open boundary
    logChord(l, L) = @. log(sin(π * l /L))/6
    
    L = length(SvN_list) + 1

    # fit scaling
    lm(x,p) = @. p[1] * x + p[2]
    xdata = logChord([1:L-1;],L); #log.(sin.(π .* [1:L-1;] ./L))./6
    fit = curve_fit(lm, xdata[mincut:L-mincut], SvN_list[mincut:L-mincut], [0.5, 0.0])
    fitparam = fit.param
    cent = fitparam[1]
    cent_err = stderror(fit)[1]
    if pbc
        cent /= 2.0
        cent_err/= 2.0
    end
    println("cent ± cent_err is $(cent) ± $(cent_err)")

    return (cent, cent_err)
end

@testset "central_charge" begin
    N = 16
    τ = atanh(0.8) # critical point for IsingX
    model = AnyonModel(FibonacciAnyon(), N;)

    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0
    t = 6N
    samples = BitMatrix(zeros(Int8, 2t, div(N,2)))
    measure_config = MeasureConfig(τ=τ, t₂=t, mode=:sample)

    measure_outcome = bulk_evolution(model, st, measure_config, samples)
    generated_statelis, F = measure_outcome.states, measure_outcome.free_energys
    final_st = generated_statelis[end]
    EE = anyon_eelis(model, final_st)
    @test fitCCEntEntScal(EE, mincut=2, pbc=true)[1][1] ≈ 0.8 atol=1e-1

    samples = BitMatrix(ones(Int8, 15N, div(N,2)))
    measure_config = MeasureConfig(τ=τ, t₂=div(15N, 2), mode=:sample)
    ψ0, sites = initial_mps(N)
    measure_outcome_mps = bulk_evolution(model, sites, ψ0, measure_config, samples)
    generated_statelis_mps, F = measure_outcome_mps.states, measure_outcome_mps.free_energys
    final_mps = generated_statelis_mps[end]
    EE_mps = anyon_eelis(model, final_mps)
    c = fitCCEntEntScal(EE_mps, mincut=4, pbc =true)[1]
    @test c ≈ 0.7 atol=1e-1
end 