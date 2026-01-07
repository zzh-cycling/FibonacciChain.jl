using FibonacciChain
using Test
using BitBasis
using LinearAlgebra 
using Random

@testset "antimap" begin
    ϕ = (1+√5)/2
    @test FibonacciChain.antimap(bit"000", 2) == (bit"000", bit"010", -ϕ^(-1), -ϕ^(-3/2))
    @test FibonacciChain.antimap(bit"010", 2) == (bit"010", bit"000",  -ϕ^(-2), -ϕ^(-3/2))
end

@testset "count_subBitStr" begin
    @test FibonacciChain.count_subBitStr(bit"00000") == 0
    @test FibonacciChain.count_subBitStr(bit"10000") == 0
    @test FibonacciChain.count_subBitStr(bit"10100") == 1
    @test FibonacciChain.count_subBitStr(bit"00100") == 0
    @test FibonacciChain.count_subBitStr(bit"10101") == 2
    @test FibonacciChain.count_subBitStr(bit"00101") == 1
    @test FibonacciChain.count_subBitStr(bit"010101") == 2 # Such config will be added additionally in PBC
end

@testset "actingHamobc" begin
    ϕ = (1+√5)/2
    model = AnyonModel(FibonacciAnyon(), 3, pbc=false)
    output1 = FibonacciChain.actingHam(model, bit"000") 
    states, weights = keys(output1), values(output1)
    @test [states...]== BitStr{3}.([bit"000", bit"010"])
    @test [weights...] ≈ [-ϕ^(-1), -ϕ^(-3/2)]
    output2 = FibonacciChain.actingHam(model, bit"010") 
    states, weights = keys(output2), values(output2)
    @test [states...]== BitStr{3}.([bit"000", bit"010"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output3 = FibonacciChain.actingHam(model, bit"001") 
    states, weights = keys(output3), values(output3)
    @test [states...]== BitStr{3}.([bit"001"])
    @test [weights...] ≈ [0.0]
    output4 = FibonacciChain.actingHam(model, bit"100") 
    states, weights = keys(output4), values(output4)
    @test [states...]== BitStr{3}.([bit"100"])
    @test [weights...] ≈ [0.0]
    output = FibonacciChain.actingHam(model, bit"101")
    states, weights = keys(output), values(output)
    @test [states...]== BitStr{3}.([bit"101"])
    @test [weights...] ≈ [-1.0]
end

@testset "Fsymmetry_coef" begin
    ϕ = (1+√5)/2
    N = 3
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    T = BitStr{N}
    output1 = FibonacciChain.Fsymmetry_coef(model, T(bit"000"), T(bit"010"))
    @test output1 ≈ -ϕ^(-3/2)
    output2 = FibonacciChain.Fsymmetry_coef(model, T(bit"010"), T(bit"000"))
    @test output2 ≈ -ϕ^(-3/2)
    output3 = FibonacciChain.Fsymmetry_coef(model, T(bit"000"), T(bit"000"))
    @test output3 ≈ -ϕ^(-3)
end

@testset "topological_symmetry_basismap" begin
    N = 4
    T = BitStr{N}
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    ϕ = (1+√5)/2
    @test FibonacciChain.topological_symmetry_basismap(model, T(bit"0000")) ≈ [ϕ^(-4), ϕ^(-5/2), ϕ^(-5/2), ϕ^(-5/2), ϕ^(-1), ϕ^(-5/2), ϕ^(-1)]
    @test FibonacciChain.topological_symmetry_basismap(model, T(bit"0100")) ≈ [ϕ^(-5/2), ϕ^(-1), -ϕ^(-2), ϕ^(-2), ϕ^(-1/2), -ϕ^(-2), ϕ^(-3/2)]
    @test FibonacciChain.topological_symmetry_basismap(model, T(bit"1010")) ≈ [ϕ^(-1), ϕ^(-3/2), ϕ^(-1/2), ϕ^(-3/2), ϕ^(-2), ϕ^(-1/2), 1]
end


@testset "topological_charge_operator" begin
    ϕ = (1+√5)/2
    N = 4
    T = BitStr{N}
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    Y =FibonacciChain.topological_charge_operator(model, T)
    Y = (Y + Y') / 2
    vals = eigvals(Y)
    @test vals ≈ [-1.0799610383969367, -0.23606797749978994, -0.23606797749978972, -0.23606797749978964, 0.4778136965285674, 1.9041523147215358, 3.079961038396939]
    @test vals[end]/vals[end-1] ≈ ϕ atol = 1e-3
end

@testset "actingHampbc" begin
    ϕ = (1+√5)/2
    model = AnyonModel(FibonacciAnyon(), 3, pbc=true)
    output1 = FibonacciChain.actingHam(model, bit"000") 
    states, weights = keys(output1), values(output1)
    @test [states...]== BitStr{3}.([bit"000",bit"100", bit"010", bit"001"])
    @test [weights...] ≈ [-3ϕ^(-1), -ϕ^(-3/2), -ϕ^(-3/2), -ϕ^(-3/2)]
    output2 = FibonacciChain.actingHam(model, bit"010") 
    states, weights = keys(output2), values(output2)
    @test [states...]== BitStr{3}.([bit"000", bit"010"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output3 = FibonacciChain.actingHam(model, bit"001") 
    states, weights = keys(output3), values(output3)
    @test [states...]== BitStr{3}.([bit"000", bit"001"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output4 = FibonacciChain.actingHam(model, bit"100") 
    states, weights = keys(output4), values(output4)
    @test [states...]== BitStr{3}.([bit"000",bit"100"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output = FibonacciChain.actingHam(AnyonModel(FibonacciAnyon(), 10), bit"1000010000")
    states, weights = keys(output), values(output)
    @test [states...] == BitStr{10}.([bit"1000010000", bit"0000010000",bit"1010010000", bit"1000010010", bit"1000010100", bit"1000000000", bit"1001010000"])
    @test [weights...] ≈ vcat([-(4ϕ^(-1)+2ϕ^(-2))],fill(-ϕ^(-3/2),6))
end

@testset "basis.jl" begin
    # Test the Fibonacci basis creation
    fib_basis = anyon_basis(AnyonModel(FibonacciAnyon(), 5, pbc=true))
    @test length(fib_basis) == 11
    fib_basis = anyon_basis(AnyonModel(FibonacciAnyon(), 5, pbc=false))
    @test length(fib_basis) == 13
    # Test the Fibonacci Hamiltonian
    fib_ham = anyon_ham(AnyonModel(FibonacciAnyon(), 5, pbc=true))
    @test size(fib_ham) == (11, 11)
    @test ishermitian(fib_ham)

    @test anyon_ham(AnyonModel(FibonacciAnyon(), 3, pbc=false)) ≈ [-0.6180339887498948 0.0 -0.48586827175664565 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; -0.48586827175664565 0.0 -0.3819660112501051 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0]
    @test anyon_ham(AnyonModel(FibonacciAnyon(), 3, pbc=true)) ≈ [-1.8541019662496843 -0.48586827175664565 -0.48586827175664565 -0.48586827175664565; -0.48586827175664565 -0.3819660112501051 0.0 0.0; -0.48586827175664565 0.0 -0.3819660112501051 0.0; -0.48586827175664565 0.0 0.0 -0.3819660112501051]
end

@testset "process_join" begin
    # [[000 ₍₂₎, 001 ₍₂₎, 010 ₍₂₎, 100 ₍₂₎, 101 ₍₂₎], [000 ₍₂₎, 001 ₍₂₎, 010 ₍₂₎, 100 ₍₂₎, 101 ₍₂₎]]
    lis1 = BitStr{2}[0, 1, 2]
    lis2 = BitStr{3}[0, 1, 2, 4, 5]
    res = FibonacciChain.process_join(lis1, lis2) 
    @test res == [join(l1, l2) for l1 in lis1 for l2 in lis2]

    # joint_basis
    res = FibonacciChain.joint_basis(FibonacciAnyon(), [2, 3])
    @test res == [join(l1, l2) for l1 in lis1 for l2 in lis2]

    lis1 = anyon_basis(AnyonModel(FibonacciAnyon(), 1, pbc=false))
    lis2 = anyon_basis(AnyonModel(FibonacciAnyon(), 2, pbc=false))
    res = FibonacciChain.joint_basis(FibonacciAnyon(), [1, 2])
    # Ensureing the order is 2*1, not 1*2
    # kron(st1*st2') is in order 1*2, while reshape(st1*st2',9) is in order 2*1
    @test res == vec([join(l1, l2) for l1 in lis1 for l2 in lis2])
    @test res != [join(l2, l1) for l1 in lis1, l2 in lis2]


    # Test disjoint system joint_basis
    res = FibonacciChain.joint_basis(FibonacciAnyon(), FibonacciAnyon(), [1], [2])
    @test res == [join(l1, l2) for l1 in anyon_basis(AnyonModel(FibonacciAnyon(), 1, pbc=false)) for l2 in anyon_basis(AnyonModel(FibonacciAnyon(), 2, pbc = false))]
    res = FibonacciChain.joint_basis(IsingAnyon(), IsingAnyon(), [1], [2])
    @test res == vec([join(l1, l2) for l1 in anyon_basis(AnyonModel(IsingAnyon(), 1, pbc=false)) for l2 in anyon_basis(AnyonModel(IsingAnyon(), 2, pbc = false))])
    res = FibonacciChain.joint_basis(IsingAnyon(), FibonacciAnyon(), [1], [2])
    @test res == vec([join(l1, l2) for l1 in anyon_basis(AnyonModel(IsingAnyon(), 1, pbc=false)) for l2 in anyon_basis(AnyonModel(FibonacciAnyon(), 2, pbc = false))])
    res = FibonacciChain.joint_basis(FibonacciAnyon(), IsingAnyon(), [1], [2])

    # move_subsystem
    res = FibonacciChain.move_subsystem(BitStr{5, Int}, bit"111", [1, 2, 5])
    @test res == BitStr{5}(bit"11001")

    # takeenviron
    bs, mask = BitStr{5}(0b11001), BitStr{5}(0b10001)
    env = FibonacciChain.takeenviron(bs, mask)
    sys = FibonacciChain.takesystem(bs, mask)
    @test env == BitStr{5}(0b01000)
    @test sys == BitStr{5}(0b10001)
end

@testset "connected components" begin
    v = [1, 2, 4, 5, 7]
    @test FibonacciChain.connected_components(v) == [[1, 2], [4, 5], [7]]
    @test FibonacciChain.connected_components([1,2,3,7,8,9]) == [[1, 2, 3], [7, 8, 9]]
end

@testset "anyon_rdm" begin
    N = 3
    model = AnyonModel(FibonacciAnyon(), N)
    st = ones(length(anyon_basis(model))); st /= norm(st)  # Normalize the state
    # The empty subsystem
    rdm = anyon_rdm(model, Int64[], st)
    @test rdm ≈ ones(Float64, 1,1)

    # The total system
    rdm = anyon_rdm(model, collect(1:N), st)
    @test rdm ≈ st*st'

    rdm = anyon_rdm(model, [1], st)
    @test rdm ≈ [0.75 0.25; 0.25 0.25]
end

@testset "anyon_rdm_Ising" begin
    N = 3
    st = zeros(2^N); st[1]=1; st[end]=1; st /= norm(st)  # Normalize the state
    model = AnyonModel(IsingAnyon(), N)
    rdm = anyon_rdm(model, Int[], st)
    @test rdm ≈ ones(Float64, 1,1)

    rdm = anyon_rdm(model, collect(1:N), st)
    @test rdm ≈ st*st'

    rdm = anyon_rdm(model, [1], st)
    @test rdm ≈ 0.5*I(2)
end

@testset "anyon_rdm_matrix" begin
    N = 4
    model = AnyonModel(FibonacciAnyon(), N)
    st = zeros(length(anyon_basis(model)));st[5]=1; st[end]=1; st /= norm(st)  # Normalize the state

    rdm = anyon_rdm(model, collect(1:2), st*st')
    @test rdm ≈ diagm([0.0, 0.5, 0.5])
    rdm = anyon_rdm(model, collect(1:N), st*st')
    @test rdm ≈ st*st'

    model2 = AnyonModel(IsingAnyon(), N)
    st = zeros(2^N); st[1]=1; st[end]=1; st /= norm(st)  # Normalize the state

    rdm = anyon_rdm(model2, [1,2], st*st') ≈ diagm([0.5, 0.0, 0.0, 0.5])

    st = zeros(2^2); st[1]=1; st[end]=1; st /= norm(st)  # Normalize the state
    rdm = anyon_rdm(AnyonModel(IsingAnyon(), 2), [1], st*st') ≈ 0.5 * I(2)
end

@testset "anyon_basis_K" begin
    N = 8
    model = AnyonModel(FibonacciAnyon(), N)
    k = 0
    fib_basis_k = anyon_basis(model, k)[1]
    @test length(fib_basis_k) == 8  # Check the length of the basis
    @test [i.buf for i in fib_basis_k] == [0, 1, 5, 9, 17, 21, 37, 85]
end

@testset "anyon_ham_K" begin
    N = 8
    model = AnyonModel(FibonacciAnyon(), N)
    k = 0
    fib_ham_k = anyon_ham(model, k)
    H = anyon_ham(model)
    @test size(fib_ham_k) == (8, 8)  # Check the size of the Hamiltonian matrix
    @test ishermitian(fib_ham_k)  # Check if the Hamiltonian is Hermitian
    @test eigvals(H)[1] ≈ eigvals(fib_ham_k)[1]  # Check if the ground state energy matches
end

@testset "mapst_sec2tot, anyon_rdm_sec" begin
    N = 8
    k = 0
    model = AnyonModel(FibonacciAnyon(), N)
    fib_ham_k = anyon_ham(model, k)
    H = anyon_ham(model)

    sec_gs = eigvecs(fib_ham_k)[:, 1]
    gs = eigvecs(H)[:, 1]
    mapped_st = FibonacciChain.mapst_sec2tot(model, sec_gs, k)
    @test isapprox(abs.(mapped_st), abs.(gs), atol=1e-9) # Check if the mapped state matches the ground state

    rdm_sec = anyon_rdm_sec(model, collect(1:div(N,2)), sec_gs, k)
    rdm = anyon_rdm(model, collect(1:div(N,2)), gs)
    @test rdm_sec ≈ rdm  # Check if the reduced density matrix matches
end

@testset "disjoint_rdm" begin
    N1 = 4
    N2 = 4
    model1 = AnyonModel(FibonacciAnyon(), N1)
    model2 = AnyonModel(FibonacciAnyon(), N2)
    state = zeros(length(anyon_basis(model1)) * length(anyon_basis(model2))); state[1] = 1; state[end] = 1
    state = state ./ norm(state)  # Normalize the state
    subsystemsA = [1, 2]
    subsystemsB = [1, 2]

    rdm_result = disjoint_rdm(model1, model2, subsystemsA, subsystemsB, state)
    @test size(rdm_result) == (9, 9)  # Check the size of the reduced density matrix
    @test rdm_result[1,1] ≈ rdm_result[end, end] ≈ 0.5
    
    # Test for one subsystem is empty
    rdm_result_empty1 = disjoint_rdm(model1, model2, Int64[], subsystemsB, state)
    @test diag(rdm_result_empty1) ≈ [0.5, 0.0, 0.5]

    rdm_result_empty2 = disjoint_rdm(model1, model2, subsystemsA, Int64[], state)
    @test diag(rdm_result_empty2) ≈ [0.5, 0.0, 0.5]
end

@testset "disjoint_rdm_Ising" begin
    N1 = 4
    N2 = 4
    model1 = AnyonModel(IsingAnyon(), N1)
    model2 = AnyonModel(IsingAnyon(), N2)
    state = zeros(length(anyon_basis(model1)) * length(anyon_basis(model2))); state[1] = 1; state[end] = 1
    state = state ./ norm(state)  # Normalize the state
    subsystemsA = [1, 2]
    subsystemsB = [1, 2]

    rdm_result = disjoint_rdm(model1, model2, subsystemsA, subsystemsB, state)
    @test size(rdm_result) == (16, 16)  # Check the size of the reduced density matrix
    @test rdm_result[1,1] ≈ rdm_result[end, end] ≈ 0.5

    rdm_result_empty1 = disjoint_rdm(model1, model2, Int64[], subsystemsB, state)
    @test all(diag(rdm_result_empty1) ≈ [0.5, 0.0, 0.0, 0.5])

    rdm_result_empty2 = disjoint_rdm(model1, model2, subsystemsA, Int64[], state)
    @test all(diag(rdm_result_empty2) ≈ [0.5, 0.0, 0.0, 0.5])

end

@testset "disjoint_rdm_mixed" begin
    N1 = 4
    N2 = 4
    model1 = AnyonModel(FibonacciAnyon(), N1)
    model2 = AnyonModel(IsingAnyon(), N2)
    state = zeros(length(anyon_basis(model1)) * length(anyon_basis(model2))); state[1] = 1; state[end] = 1
    state = state ./ norm(state)  # Normalize the state
    subsystemsA = [1, 2]
    subsystemsB = [1, 2]

    rdm_result = disjoint_rdm(model1, model2, subsystemsA, subsystemsB, state)
    @test size(rdm_result) == (12, 12)  # Check the size of the reduced density matrix
    @test rdm_result[1,1] ≈ rdm_result[end, end] ≈ 0.5
    
    # Test for one subsystem is empty
    rdm_result_empty1 = disjoint_rdm(model1, model2, Int64[], subsystemsB, state)
    @test diag(rdm_result_empty1) ≈ [0.5, 0.0, 0.0, 0.5]

    rdm_result_empty2 = disjoint_rdm(model1, model2, subsystemsA, Int64[], state)
    @test diag(rdm_result_empty2) ≈ [0.5, 0.0, 0.5]
end

@testset "disjoint_rdm and ref qubit" begin
    L=6
    model = AnyonModel(IsingAnyon(), L, pbc=true, measure_operator = :X)
    t₂ = 10
    sample = BitMatrix(ones(Int, 2t₂, L))
    τ = log(1 + √2)
    config = MeasurementConfig(τ = τ, t₂ = t₂, rng=MersenneTwister(1234), mode =:sample)
    
    initial_state = zeros(length(anyon_basis(model)))
    initial_state[1] = 1.0
    mo = bulk_evolution(model, initial_state, config,sample)
    statelis = mo.states
    
    ref_model = AnyonModel(IsingAnyon(), 2, pbc=false, measure_operator=:X)
    ref_congfig = MeasurementConfig(τ = τ, t₂ = t₂, rng=MersenneTwister(1234), mode =:sample)
    stlis, samples, sample_free_energy = reference_evolution(model, statelis, ref_config, sample)
    st=stlis[end]
    ρ1 = disjoint_rdm(ref_model, model, Int64[], collect(1:div(L,2)), st)
    ρ2 = disjoint_rdm(ref_model, model, Int64[], collect(1:L), st)
    ρ2r = anyon_rdm(model, collect(1:div(L,2)), ρ2)
    S1 = ee(ρ1)
    S2 = ee(ρ2r)
    @test S1 ≈ S2
end