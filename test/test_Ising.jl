using FibonacciChain
using Test
using LinearAlgebra
using BitBasis
using Arpack
using Random
using LsqFit

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
    model = AnyonModel(IsingAnyon(), N, pbc=false)
    T = BitStr{N, Int}
    output1 = FibonacciChain.actingHam(model, T, bit"000", J=2.0, h=1.0) 
    states, weights = keys(output1), values(output1)
    @test [states...]== T.([bit"000", bit"100", bit"010", bit"001"])
    @test [weights...] ≈ [-4.0, -1.0, -1.0, -1.0]

    output2 = FibonacciChain.actingHam(model, T, bit"010") 
    states, weights = keys(output2), values(output2)
    @test [states...]== T.([bit"000", bit"110", bit"010", bit"011"])
    @test [weights...] ≈ [-1.0, -1.0, 2.0, -1.0]

    output3 = FibonacciChain.actingHam(model, T, bit"001") 
    states, weights = keys(output3), values(output3)
    @test [states...]== T.([bit"101", bit"000", bit"011", bit"001"])
    @test [weights...] ≈ [-1.0, -1.0, -1.0, 0.0]

    output4 = FibonacciChain.actingHam(model, T, bit"100") 
    states, weights = keys(output4), values(output4)
    @test [states...]== T.([bit"000", bit"100", bit"110", bit"101"])
    @test [weights...] ≈ [-1.0, 0.0, -1.0, -1.0]

    output = FibonacciChain.actingHam(model, T, bit"101")
    states, weights = keys(output), values(output)
    @test [states...]== T.([bit"101", bit"100", bit"111", bit"001"])
    @test [weights...] ≈ [2.0, -1.0, -1.0, -1.0]
end

@testset "actingHampbc" begin
    N = 3
    T = BitStr{N, Int}
    model = AnyonModel(IsingAnyon(), N, pbc=true)
    output1 = FibonacciChain.actingHam(model, T, bit"000", J=2.0, h=1.0) 
    states, weights = keys(output1), values(output1)
    @test [states...]== T.([bit"000",bit"100", bit"010", bit"001"])
    @test [weights...] ≈ [-6.0, -1.0, -1.0, -1.0]

    output2 = FibonacciChain.actingHam(model, T, bit"010") 
    states, weights = keys(output2), values(output2)
    @test [states...]== T.([bit"000", bit"110", bit"010", bit"011"])
    @test [weights...] ≈ [-1.0, -1.0, 1.0, -1.0]

    output3 = FibonacciChain.actingHam(model, T, bit"001") 
    states, weights = keys(output3), values(output3)
    @test [states...]== T.([bit"101", bit"000", bit"011", bit"001"])
    @test [weights...] ≈ [-1.0, -1.0, -1.0, 1.0]

    output4 = FibonacciChain.actingHam(model, T, bit"100") 
    states, weights = keys(output4), values(output4)
    @test [states...]== T.([bit"000",bit"100", bit"110", bit"101"])
    @test [weights...] ≈ [-1.0, 1.0, -1.0, -1.0]

    output = FibonacciChain.actingHam(AnyonModel(IsingAnyon(), 10, pbc=true), BitStr{10}, bit"1000010000")
    states, weights = keys(output), values(output)
    @test [states...] == BitStr{10}.([16, 560, 656, 532, 529, 784, 512, 536, 528, 530, 592])
    @test [weights...] ≈ vcat(-ones(8),[-2.0, -1.0, -1.0])
end

@testset "basis.jl" begin
    # Test the Ising basis creation
    fib_basis = anyon_basis(AnyonModel(IsingAnyon(), 5, pbc=false))
    @test length(fib_basis) == 32
    fib_basis = anyon_basis(AnyonModel(IsingAnyon(), 5, pbc=true))
    @test length(fib_basis) == 32
    # Test the Ising Hamiltonian
    fib_ham = anyon_ham(AnyonModel(IsingAnyon(), 5, pbc=true))
    @test size(fib_ham) == (32, 32)
    @test ishermitian(fib_ham)
    # Test the ground state energy
    @test eigvals(fib_ham)[1] ≈ -1/(2*sin(π/10))*4 atol=1e-10

    X= Float64[0 1; 1 0]
	Z= Float64[1 0; 0 -1]
	Id= Float64[1 0; 0 1]


    ⊗(A::AbstractArray, B::AbstractArray) = kron(A, B)
    H_temp = - (Z ⊗ Z ⊗ Id + Id ⊗ Z ⊗ Z  + X ⊗ Id ⊗ Id + Id ⊗ Id ⊗ X + Id ⊗ X ⊗ Id)

    @test anyon_ham(AnyonModel(IsingAnyon(), 3, pbc=false)) == H_temp

    H = anyon_ham(AnyonModel(IsingAnyon(), 3, pbc=true))
    @test H == H_temp - Z ⊗ Id ⊗ Z

    gs = eigvecs(fib_ham)[:,1]
    # Test the reduced density matrix function
    rdm = anyon_rdm(AnyonModel(IsingAnyon(), 5, pbc=true), collect(1:1), gs)
    @test size(rdm) == (2, 2)
    @test rdm ≈ [0.49999999999995276 0.3236067977499789; 0.3236067977499789 0.5000000000000473]
end

@testset "measure_basismap_IsingX" begin
    N = 3
    T = BitStr{N, Int}
    idx = 2
    pbc = false
    τ =0.0
    cstτ = cosh(τ/2) / √(2cosh(τ))
    sign = false
    basis0 = [T(0b000), T(0b001), T(0b010), T(0b100), T(0b101)]
    model = AnyonModel(IsingAnyon(), N, pbc=pbc, measure_operator=:X)

    output = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, 0.0)
    @test output[2] == (T(bit"001"), T(bit"011"), cstτ, 0.0)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, 0.0)
    @test output[4] == (T(bit"100"), T(bit"110"), cstτ, 0.0)
    @test output[5] == (T(bit"101"), T(bit"111"), cstτ, 0.0)

    sign = true
    output2 = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test output2 == output

    τ = 1.0
    sign = false
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, coef)
    @test output[2] == (T(bit"001"), T(bit"011"), cstτ, coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, coef)
    @test output[4] == (T(bit"100"), T(bit"110"), cstτ, coef)
    @test output[5] == (T(bit"101"), T(bit"111"), cstτ, coef)

    sign = true
    coef = -sinh(τ/2) / √(2cosh(τ))
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, coef)
    @test output[2] == (T(bit"001"), T(bit"011"), cstτ, coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, coef)
    @test output[4] == (T(bit"100"), T(bit"110"), cstτ, coef)
    @test output[5] == (T(bit"101"), T(bit"111"), cstτ, coef)

    τ = 1e3
    sign = false
    cstτ = 1/2
    coef = 1/2
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, coef)
    @test output[2] == (T(bit"001"), T(bit"011"), cstτ, coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, coef)
    @test output[4] == (T(bit"100"), T(bit"110"), cstτ, coef)
    @test output[5] == (T(bit"101"), T(bit"111"), cstτ, coef)


    idx = 3
    output = measure_basismap.(Ref(AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)), τ, basis0, idx, sign)
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
    sign = false
    basis0 = [T(0b000), T(0b001), T(0b010), T(0b100), T(0b101)]
    model = AnyonModel(IsingAnyon(), N, pbc=pbc, measure_operator=:ZZ)

    output = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"000"), cstτ, 0.0)
    @test output[2] == (T(bit"001"), T(bit"001"), cstτ, 0.0)
    @test output[3] == (T(bit"010"), T(bit"010"), cstτ, 0.0)
    @test output[4] == (T(bit"100"), T(bit"100"), cstτ, 0.0)
    @test output[5] == (T(bit"101"), T(bit"101"), cstτ, 0.0)

    sign = true
    output2 = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test output2 == output

    τ = 1.0
    sign = false
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"000"), cstτ+coef, 0.0)
    @test output[2] == (T(bit"001"), T(bit"001"), cstτ-coef, 0.0)
    @test output[3] == (T(bit"010"), T(bit"010"), cstτ-coef, 0.0)
    @test output[4] == (T(bit"100"), T(bit"100"), cstτ+coef, 0.0)
    @test output[5] == (T(bit"101"), T(bit"101"), cstτ-coef, 0.0)

    sign = true
    coef = -sinh(τ/2) / √(2cosh(τ))
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"000"), cstτ+coef, 0.0)
    @test output[2] == (T(bit"001"), T(bit"001"), cstτ-coef, 0.0)
    @test output[3] == (T(bit"010"), T(bit"010"), cstτ-coef, 0.0)
    @test output[4] == (T(bit"100"), T(bit"100"), cstτ+coef, 0.0)
    @test output[5] == (T(bit"101"), T(bit"101"), cstτ-coef, 0.0)

    idx=3
    τ = 1e3
    sign = false
    cstτ = 1/2
    coef = 1/2
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:ZZ)
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"000"), 1.0, 0.0)
    @test output[2] == (T(bit"001"), T(bit"001"), 0.0, 0.0)
    @test output[3] == (T(bit"010"), T(bit"010"), 1.0, 0.0)
    @test output[4] == (T(bit"100"), T(bit"100"), 0.0, 0.0)
    @test output[5] == (T(bit"101"), T(bit"101"), 1.0, 0.0)


    sign = true
    output = measure_basismap.(Ref(model), τ, basis0, idx, sign)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"000"), 0.0, 0.0)
    @test output[2] == (T(bit"001"), T(bit"001"), 1.0, 0.0)
    @test output[3] == (T(bit"010"), T(bit"010"), 0.0, 0.0)
    @test output[4] == (T(bit"100"), T(bit"100"), 1.0, 0.0)
    @test output[5] == (T(bit"101"), T(bit"101"), 0.0, 0.0)
end

@testset "measure_matrix" begin
    N = 3
    
    ⊗(a,b) = kron(a, b)
    τ = 1.0
    idx = 2
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))

    σx = [0.0 1.0; 1.0 0.0]
    σz = [1.0 0.0; 0.0 -1.0]
    # measuring X
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σx ⊗ I(2)

    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)
    Mpobc = FibonacciChain.measure_matrix(model, τ, idx, false)
    @test Mpobc == expected_matrix 

    coef = -sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σx ⊗ I(2)
    Mmobc = FibonacciChain.measure_matrix(model, τ, idx, true)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(8) 

    # measuring ZZ
    coef = sinh(τ/2) / √(2cosh(τ))
    model_ZZ = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:ZZ)
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σz ⊗ σz
    Mpobc = FibonacciChain.measure_matrix(model_ZZ, τ, idx, false)
    @test Mpobc == expected_matrix

    coef = -sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σz ⊗ σz
    Mmobc = FibonacciChain.measure_matrix(model_ZZ, τ, idx, true)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(8) 

    # Test with a different τ, idx
    idx = 3
    τ = 1000.0   
    cstτ = 0.5
    coef = 0.5
    expected_matrix = cstτ* I(8) + coef * σz ⊗ I(2) ⊗ σz
    Mpobc = FibonacciChain.measure_matrix(model_ZZ, τ, idx, false)
    @test Mpobc == expected_matrix

    coef = -0.5
    expected_matrix = cstτ* I(8) + coef * σz ⊗ I(2) ⊗ σz
    Mmobc = FibonacciChain.measure_matrix(model_ZZ, τ, idx, true)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(8) 
end

@testset "measuremap_IsingX" begin
    N = 3
    τ = 1.0
    idx = 2
    sign = false
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))
    model = AnyonModel(IsingAnyon(), N, measure_operator=:X)
    state = fill(1.0,2^N)
    output = measuremap(model, τ, state, idx, sign)
    @test output == (cstτ+coef) .* ones(2^N)
    
    sign = true
    coef = -sinh(τ/2) / √(2cosh(τ))
    output = measuremap(model, τ, state, idx, sign)
    @test output == (cstτ+coef) .* ones(2^N)
end

@testset "measuremap_IsingZZ" begin
    N = 3
    model = AnyonModel(IsingAnyon(), N, measure_operator=:ZZ)
    τ = 1.0
    idx = 2
    sign = false
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))

    state = fill(1.0,2^N)
    output = measuremap(model, τ, state, idx, sign)
    @test output == [cstτ+coef, cstτ-coef,cstτ-coef, cstτ+coef, cstτ+coef, cstτ-coef,cstτ-coef, cstτ+coef]
    
    sign = true
    coef = -sinh(τ/2) / √(2cosh(τ))
    output = measuremap(model, τ, state, idx, sign)
    @test output == [cstτ+coef, cstτ-coef,cstτ-coef, cstτ+coef, cstτ+coef, cstτ-coef,cstτ-coef, cstτ+coef]


    # Try with idx=3 pbc
    state = ones(2^N)
    output = measuremap(model, 1000.0, state, idx, sign)
    @test output == [0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0]
end


@testset "measurement_enumeration" begin
    N=6
    model = AnyonModel(IsingAnyon(), N, measure_operator=:X)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0 
    τ = 0.0
    measurement_sites = collect(2:2:N)

    final_states, trajectories, probabilities = measurement_enumeration(model, τ, st, measurement_sites)

    num_final_states = length(final_states)
    @test num_final_states == 2^length(measurement_sites)

    # all final states should be equally probable
    total_prob = sum(probabilities)
    @test isapprox(total_prob, 1.0, atol=1e-6)
    @test probabilities ≈ 1/8 .* ones(2^length(measurement_sites))
    @test sum(map(x->-x*log(x)/length(measurement_sites), probabilities)) ≈ log(2) # Shannon entropy non-measurement state

    model = AnyonModel(IsingAnyon(), N, measure_operator=:ZZ)
    final_states, trajectories, probabilities = measurement_enumeration(model, τ, st, measurement_sites)

    num_final_states = length(final_states)
    @test num_final_states == 2^length(measurement_sites)

    total_prob = sum(probabilities)
    @test isapprox(total_prob, 1.0, atol=1e-6)
    @test probabilities ≈ 1/8 .* ones(2^length(measurement_sites))
    @test sum(map(x->-x*log(x)/length(measurement_sites), probabilities)) ≈ log(2) # Shannon entropy non-measurement state
end

@testset "Boundary_Born" begin
    N=6
    model = AnyonModel(IsingAnyon(), N, measure_operator=:X)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0 
    τ = 3.802
    measure_config = MeasureConfig(τ=τ, t₂=1,rng = MersenneTwister(90), mode=:Born)
    measure_outcome = boundary_evolution(model, st, measure_config)
    sample_measured_state, sample, sample_free_energy = measure_outcome.state, measure_outcome.sample, measure_outcome.free_energy

    @test sample == BitVector([0,0,0,0,1,0])
    @test sample_free_energy ≈ 4.158883083359672 atol = 1e-6
end

@testset "Boundary_post_selection" begin
    N = 10
    τ = 1e3
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0

    measurement_sites = collect(1:N)
    config = MeasureConfig(τ=τ, t₂=1, mode=:sample)
    samples = BitVector(zeros(Int8, length(measurement_sites)))
    measure_outcome = boundary_evolution(model, st, config,samples)
    final_state_p, final_sequence_p, total_free_energy_p = measure_outcome.state, measure_outcome.sample, measure_outcome.free_energy

    # all final states should be equally probable, will give Nlog(2) free energy
    @test total_free_energy_p[1] /length(measurement_sites) ≈ log(2) atol = 1e-6

    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:ZZ)
    samples = BitVector(zeros(Int8, length(measurement_sites)))
    measure_outcome = boundary_evolution(model, final_state_p, config,samples, layer_idx=2)
    final_state_p, final_sequence_p, total_free_energy_p = measure_outcome.state, measure_outcome.sample, measure_outcome.free_energy

    @test total_free_energy_p[1] /length(measurement_sites) ≈ 0.6238324625039509 atol = 1e-6
    @test final_state_p[1] ≈ final_state_p[end] ≈ 1/√2
end

@testset "Bulk_Born" begin
    L = 6
    D = 2L
    τ = 1e3
    model = AnyonModel(IsingAnyon(), L, measure_operator=:X)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0

    config = MeasureConfig(τ=τ, t₂=D, rng = MersenneTwister(100), mode=:Born)
    measure_outcome = bulk_evolution(model, st, config)
    sample_measured_states, samples, sample_free_energy = measure_outcome.states, measure_outcome.samples, measure_outcome.free_energys
    EElis = [anyon_eelis(model, state_t)[div(L,2)] for state_t in sample_measured_states]
    @test size(samples) == (2D, L)
    # Each layer will erase previous info.
    @test EElis ≈ [log(2) for i in 1:D] atol = 1e-6
end

@testset "Bulk_post_selection" begin
    L = 6
    τ = 1000.0
    D = 10L
    model = AnyonModel(IsingAnyon(), L, pbc=true, measure_operator=:X)
    st=zeros(length(anyon_basis(model)))
    st[1] = 1.0
    average_EElis=zeros(L-1)

    EE_tlis = zeros(D)
    config = MeasureConfig(τ=τ, t₂=D, mode=:sample)
    samples= BitMatrix(zeros(Int, 2D, L))
    measure_outcome = bulk_evolution(model, st, config, samples)
    sample_measured_states, samples, sample_free_energy = measure_outcome.states, measure_outcome.samples, measure_outcome.free_energys
    state_t = sample_measured_states[end]
    EE = anyon_eelis(model, state_t)
    @test samples[end,:] == fill(0, L)
    @test EE ≈ log(2)*ones(L-1) atol = 1e-4
end

@testset "_apply_measurement_layer" begin
    N = 6
    τ = 1e3
    model = AnyonModel(IsingAnyon(), N, measure_operator=:X)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0

    config = MeasureConfig(τ=τ, t₂=N, rng = MersenneTwister(100), mode=:Born)
    measure_outcome = bulk_evolution(model, st, config)
    sample_measured_states, samples, sample_free_energy = measure_outcome.states, measure_outcome.samples, measure_outcome.free_energys

    measure_outcome_layer1 = FibonacciChain._apply_measurement_layer(model, τ, st, samples[1,:], layer_idx = 1)
    new_state, F1 = measure_outcome_layer1.state, measure_outcome_layer1.free_energy
    measure_outcome_layer2 = FibonacciChain._apply_measurement_layer(model, τ, new_state, samples[2,:], layer_idx = 2)
    new_state, F2 = measure_outcome_layer2.state, measure_outcome_layer2.free_energy
    @test new_state ≈ sample_measured_states[1]
    @test F1 ≈ sample_free_energy[1] atol=1e-6
    @test F2 ≈ sample_free_energy[2] atol=1e-6
end

@testset "boundary_evolution, bulk_evolution" begin
    N = 6
    τ = 1e3
    model = AnyonModel(IsingAnyon(), N, measure_operator=:X)
    st=zeros(length(anyon_basis(model)))
    st[1] = 1.0

    config = MeasureConfig(τ=τ, t₂=1, rng = MersenneTwister(100), mode=:Born)
    measure_outcome = boundary_evolution(model, st, config) # default layer_idx=1
    sample_measured_states, sample, sample_free_energy = measure_outcome.state, measure_outcome.sample, measure_outcome.free_energy
    config_post = MeasureConfig(τ=τ, t₂=1, mode=:sample)
    measure_outcome_post = boundary_evolution(model, st, config_post, sample)
    state, F = measure_outcome_post.state, measure_outcome_post.free_energy
    @test state ≈ sample_measured_states
    @test F ≈ sample_free_energy atol=1e-6

    config = MeasureConfig(τ=τ, t₂=N, rng = MersenneTwister(100), mode=:Born)
    measure_outcome = bulk_evolution(model, st, config)
    sample_measured_states, samples, sample_free_energy = measure_outcome.states, measure_outcome.samples, measure_outcome.free_energys
    config_generate = MeasureConfig(τ=τ, t₂=N, mode=:sample)
    measure_outcome_post = bulk_evolution(model, st, config_generate, samples)
    statelis, F = measure_outcome_post.states, measure_outcome_post.free_energys
    @test statelis ≈ sample_measured_states
    @test F ≈ sample_free_energy atol=1e-6  
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
    N = 12
    τ = log(1+√2) # critical point for IsingX
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)
    st=zeros(length(anyon_basis(model)))
    st[1] = 1.0

    t = 5N
    samples = BitMatrix(zeros(Int, 2t, N))
    config = MeasureConfig(τ=τ, t₂=t, mode=:sample)

    measure_outcome = bulk_evolution(model, st, config, samples)
    generated_statelis, F = measure_outcome.states, measure_outcome.free_energys
    final_st = generated_statelis[end]
    EE = anyon_eelis(model, final_st)
    @test fitCCEntEntScal(EE, mincut=2, pbc =true)[1][1] ≈ 0.5 atol=1e-2
    
    ψ0, sites = initial_mps(N)
    measure_outcome_mps = bulk_evolution(model, sites, ψ0, config, samples, cutoff=1e-10, maxdim=1000)
    generated_statelis_mps, F = measure_outcome_mps.states, measure_outcome_mps.free_energys
    final_mps = generated_statelis_mps[end]
    EE_mps = anyon_eelis(model, final_mps)
    @test EE_mps ≈ EE atol = 1e-5
    c = fitCCEntEntScal(EE_mps, mincut=2, pbc =true)[1]
    @test c ≈ 0.5 atol=1e-2
end 