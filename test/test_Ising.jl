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

    model = AnyonModel(IsingAnyon(), N, pbc=false, measure_operator=:X)
    Mpobc = FibonacciChain.measure_matrix(model, τ, idx, false)
    @test Mpobc == expected_matrix 

    coef = -sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σx ⊗ I(2)
    Mmobc = FibonacciChain.measure_matrix(model, τ, idx, true)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(8) 

    # measuring ZZ
    coef = sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σz ⊗ σz
    Mpobc = FibonacciChain.measure_matrix(model, τ, idx, false)
    @test Mpobc == expected_matrix

    coef = -sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ* I(8) + coef * I(2) ⊗ σz ⊗ σz
    Mmobc = FibonacciChain.measure_matrix(model, τ, idx, true)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 ≈ I(8) 
w
    # Test with a different τ, idx
    idx = 3
    τ = 1000.0   
    cstτ = 0.5
    coef = 0.5      
    expected_matrix = cstτ* I(8) + coef * σz ⊗ I(2) ⊗ σz 
    Mpobc = FibonacciChain.measure_matrix(model, τ, idx, false)
    @test Mpobc == expected_matrix 

    coef = -0.5
    expected_matrix = cstτ* I(8) + coef * σz ⊗ I(2) ⊗ σz
    Mmobc = FibonacciChain.measure_matrix(model, τ, idx, true)
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

@testset "Boundary_measure" begin
    N=6
    model = AnyonModel(IsingAnyon(), N, measure_operator=:X)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0 
    τ = 3.802
    measure_config = MeasureConfig(:random, nothing)
    sample_measured_states, samples, sample_free_energy = boundary_measure(N, τ, st, 1, 500, measure_config)

    num_final_states = length(sample_measured_states)
    @test num_final_states == 500

    @test samples[476,:] == fill(0, N)
    
end

@testset "Boundarypost_selection" begin
    N = 10
    τ = 1e3
    st = zeros(length(anyon_basis(N, anyon_type=:IsingX)))
    st[1] = 1.0 
    measurement_sites = collect(1:N)
    final_state_p, final_sequence_p, total_free_energy_p = boundary_post_selection(N, τ, st, 1, 0, anyon_type=:IsingX)

    # all final states should be equally probable, will give Nlog(2) free energy
    @test total_free_energy_p /length(measurement_sites) ≈ log(2) atol = 1e-6

    final_state_p, final_sequence_p, total_free_energy_p = boundary_post_selection(N, τ, final_state_p, 2, 0, anyon_type=:IsingZZ)

    @test total_free_energy_p /length(measurement_sites) ≈ 0.6238324625039509 atol = 1e-6
    @test final_state_p[1] ≈ final_state_p[end] ≈ 1/√2
end

@testset "Bulkmeasure" begin
    L = 6
    D = 2L
    τ = 1e3
    st = zeros(length(anyon_basis(L, anyon_type=:IsingX)))
    st[1] = 1.0

    sample_measured_states, samples, sample_free_energy = bulk_measure(L, 1000.0, st, D, MersenneTwister(100), anyon_type=:IsingX) 
    EElis = [anyon_eelis(L, state_t, anyon_type=:IsingX)[div(L,2)] for state_t in sample_measured_states]
    @test size(samples) == (2D, L)
    # Each layer will erase previous info.
    @test EElis ≈ [log(2) for i in 1:D] atol = 1e-6
end

@testset "Bulkpost_selection" begin
    L = 6
    τ = 1000.0
    D = 10L
    pbc = true
    st=zeros(length(anyon_basis(L, anyon_type=:IsingX)))
    st[1] = 1.0
    average_EElis=zeros(L-1)

    EE_tlis = zeros(D)
    sample_measured_states, samples, sample_free_energy = bulk_post_selection(L, τ, st, D, 0, pbc, anyon_type=:IsingX)
    state_t = sample_measured_states[end]
    EE = anyon_eelis(L, state_t, anyon_type=:IsingX)
    @test samples[end,:] == fill(0, L)
    @test EE ≈ log(2)*ones(L-1) atol = 1e-4
end

@testset "_apply_measurement_layer" begin
    N = 6
    τ = 1e3
    st = zeros(length(anyon_basis(N, anyon_type=:IsingX)))
    st[1] = 1.0

    sample_measured_states, samples, sample_free_energy = bulk_measure(N, τ, st, N, MersenneTwister(100), anyon_type=:IsingX)
    state_t = sample_measured_states[end]

    new_state, F1 = _apply_measurement_layer!(N, τ, st, samples[1,:], 1, anyon_type=:IsingX)
    new_state, F2 = _apply_measurement_layer!(N, τ, new_state, samples[2,:], 2, anyon_type=:IsingX)
    @test new_state ≈ sample_measured_states[1]
    @test F1 ≈ sample_free_energy[1] atol=1e-6
    @test F2 ≈ sample_free_energy[2] atol=1e-6
end

@testset "generate_state" begin
    N = 6
    τ = 1e3
    st=zeros(length(anyon_basis(N, anyon_type=:IsingX)))
    st[1] = 1.0

    sample_measured_states, samples, sample_free_energy = boundary_measure(N, τ, st, 1, 10, anyon_type=:IsingX)
    state, F = generate_state(τ, st, samples[1,:], anyon_type=:IsingX)
    @test state ≈ sample_measured_states[1]
    @test F[1] ≈ sample_free_energy[1] atol=1e-6

    sample_measured_states, samples, sample_free_energy = bulk_measure(N, τ, st, N, MersenneTwister(100), anyon_type=:IsingX)
    statelis, F = generate_state(τ, st, samples, anyon_type=:IsingX)
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
    st=zeros(length(anyon_basis(N, anyon_type=:IsingX)))
    st[1] = 1.0

    samples = zeros(Int, 10N, N)

    generated_statelis, F = generate_state(τ, st, samples, anyon_type=:IsingX)
    final_st = generated_statelis[end]
    EE = anyon_eelis(N, final_st, anyon_type=:IsingX)
    @test fitCCEntEntScal(EE, mincut=2, pbc =true)[1][1] ≈ 0.5 atol=1e-2
    
    ψ0, sites = initial_mps(N)
    generated_statelis_mps, F = generate_state_mps(τ, sites, ψ0, samples; pbc= true, anyon_type=:IsingX)
    final_mps = generated_statelis_mps[end]
    EE_mps = anyon_eelis(N, final_mps)
    @test EE_mps ≈ EE atol = 1e-6
    c = fitCCEntEntScal(EE_mps, mincut=2, pbc =true)[1]
    @test c ≈ 0.5 atol=1e-2
end 