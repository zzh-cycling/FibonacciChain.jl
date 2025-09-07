using FibonacciChain
using Test
using BitBasis
using LinearAlgebra 

@testset "antimap" begin
    ϕ = (1+√5)/2
    @test FibonacciChain.antimap(BitStr{3}, bit"000", 2) == (bit"000", bit"010", -ϕ^(-1), -ϕ^(-3/2))
    @test FibonacciChain.antimap(BitStr{3}, bit"010", 2) == (bit"010", bit"000",  -ϕ^(-2), -ϕ^(-3/2))
end

@testset "count_subBitStr" begin
    @test FibonacciChain.count_subBitStr(BitStr{5}, bit"00000") == 0
    @test FibonacciChain.count_subBitStr(BitStr{5}, bit"10000") == 0
    @test FibonacciChain.count_subBitStr(BitStr{5}, bit"10100") == 1
    @test FibonacciChain.count_subBitStr(BitStr{5}, bit"00100") == 0
    @test FibonacciChain.count_subBitStr(BitStr{5}, bit"10101") == 2
    @test FibonacciChain.count_subBitStr(BitStr{5}, bit"00101") == 1
    @test FibonacciChain.count_subBitStr(BitStr{6}, bit"010101") == 2 # Such config will be added additionally in PBC
end

@testset "actingHamobc" begin
    ϕ = (1+√5)/2
    output1 = FibonacciChain.actingHam(BitStr{3}, bit"000",false) 
    states, weights = keys(output1), values(output1)
    @test [states...]== BitStr{3}.([bit"000", bit"010"])
    @test [weights...] ≈ [-ϕ^(-1), -ϕ^(-3/2)]
    output2 = FibonacciChain.actingHam(BitStr{3}, bit"010",false) 
    states, weights = keys(output2), values(output2)
    @test [states...]== BitStr{3}.([bit"000", bit"010"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output3 = FibonacciChain.actingHam(BitStr{3}, bit"001",false) 
    states, weights = keys(output3), values(output3)
    @test [states...]== BitStr{3}.([bit"001"])
    @test [weights...] ≈ [0.0]
    output4 = FibonacciChain.actingHam(BitStr{3}, bit"100",false) 
    states, weights = keys(output4), values(output4)
    @test [states...]== BitStr{3}.([bit"100"])
    @test [weights...] ≈ [0.0]
    output = FibonacciChain.actingHam(BitStr{3}, bit"101",false)
    states, weights = keys(output), values(output)
    @test [states...]== BitStr{3}.([bit"101"])
    @test [weights...] ≈ [-1.0]
end

@testset "Fsymmetry_coef" begin
    ϕ = (1+√5)/2
    N=3
    T = BitStr{N}
    output1 = FibonacciChain.Fsymmetry_coef(T(bit"000"), T(bit"010"))
    @test output1 ≈ -ϕ^(-3/2)
    output2 = FibonacciChain.Fsymmetry_coef(T(bit"010"), T(bit"000"))
    @test output2 ≈ -ϕ^(-3/2)
    output3 = FibonacciChain.Fsymmetry_coef(T(bit"000"), T(bit"000"))
    @test output3 ≈ -ϕ^(-3)
end

@testset "topological_symmetry_basismap" begin
    N = 4
    T = BitStr{N}
    ϕ = (1+√5)/2
    @test FibonacciChain.topological_symmetry_basismap(T(bit"0000")) ≈ [ϕ^(-4), ϕ^(-5/2), ϕ^(-5/2), ϕ^(-5/2), ϕ^(-1), ϕ^(-5/2), ϕ^(-1)]
    @test FibonacciChain.topological_symmetry_basismap(T(bit"0100")) ≈ [ϕ^(-5/2), ϕ^(-1), -ϕ^(-2), ϕ^(-2), ϕ^(-1/2), -ϕ^(-2), ϕ^(-3/2)]
    @test FibonacciChain.topological_symmetry_basismap(T(bit"1010")) ≈ [ϕ^(-1), ϕ^(-3/2), ϕ^(-1/2), ϕ^(-3/2), ϕ^(-2), ϕ^(-1/2), 1]
end


@testset "topological_charge_operator" begin
    ϕ = (1+√5)/2
    N = 4
    T = BitStr{N}
    Y =FibonacciChain.topological_charge_operator(T)
    Y = (Y + Y') / 2
    vals = eigvals(Y)
    @test vals ≈ [-1.0799610383969367, -0.23606797749978994, -0.23606797749978972, -0.23606797749978964, 0.4778136965285674, 1.9041523147215358, 3.079961038396939]
    @test vals[end]/vals[end-1] ≈ ϕ atol = 1e-3
end

@testset "actingHampbc" begin
    ϕ = (1+√5)/2
    output1 = FibonacciChain.actingHam(BitStr{3}, bit"000") 
    states, weights = keys(output1), values(output1)
    @test [states...]== BitStr{3}.([bit"000",bit"100", bit"010", bit"001"])
    @test [weights...] ≈ [-3ϕ^(-1), -ϕ^(-3/2), -ϕ^(-3/2), -ϕ^(-3/2)]
    output2 = FibonacciChain.actingHam(BitStr{3}, bit"010") 
    states, weights = keys(output2), values(output2)
    @test [states...]== BitStr{3}.([bit"000", bit"010"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output3 = FibonacciChain.actingHam(BitStr{3}, bit"001") 
    states, weights = keys(output3), values(output3)
    @test [states...]== BitStr{3}.([bit"000", bit"001"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output4 = FibonacciChain.actingHam(BitStr{3}, bit"100") 
    states, weights = keys(output4), values(output4)
    @test [states...]== BitStr{3}.([bit"000",bit"100"])
    @test [weights...] ≈ [-ϕ^(-3/2), -ϕ^(-2)]
    output = FibonacciChain.actingHam(BitStr{10}, bit"1000010000")
    states, weights = keys(output), values(output)
    @test [states...] == BitStr{10}.([bit"1000010000", bit"0000010000",bit"1010010000", bit"1000010010", bit"1000010100", bit"1000000000", bit"1001010000"])
    @test [weights...] ≈ vcat([-(4ϕ^(-1)+2ϕ^(-2))],fill(-ϕ^(-3/2),6))
end

@testset "basis.jl" begin
    # Test the Fibonacci basis creation
    fib_basis = anyon_basis(5)
    @test length(fib_basis) == 11
    fib_basis = anyon_basis(5,false)
    @test length(fib_basis) == 13
    # Test the Fibonacci Hamiltonian
    fib_ham = anyon_ham(5)
    @test size(fib_ham) == (11, 11)
    @test ishermitian(fib_ham)

    @test anyon_ham(3,false) == [-0.6180339887498948 0.0 -0.48586827175664565 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; -0.48586827175664565 0.0 -0.3819660112501051 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0]
    @test anyon_ham(3) == [-1.8541019662496843 -0.48586827175664565 -0.48586827175664565 -0.48586827175664565; -0.48586827175664565 -0.3819660112501051 0.0 0.0; -0.48586827175664565 0.0 -0.3819660112501051 0.0; -0.48586827175664565 0.0 0.0 -0.3819660112501051]
end

@testset "process_join" begin
    # [[000 ₍₂₎, 001 ₍₂₎, 010 ₍₂₎, 100 ₍₂₎, 101 ₍₂₎], [000 ₍₂₎, 001 ₍₂₎, 010 ₍₂₎, 100 ₍₂₎, 101 ₍₂₎]]
    lis1 = BitStr{2}[0, 1, 2]
    lis2 = BitStr{3}[0, 1, 2, 4, 5]
    res = FibonacciChain.process_join(lis1, lis2) 
    @test res == [join(l1, l2) for l1 in lis1 for l2 in lis2]

    # joint_basis
    res = FibonacciChain.joint_basis([2, 3])
    @test res == [join(l1, l2) for l1 in lis1 for l2 in lis2]

    lis1 = anyon_basis(1,false)
    lis2 = anyon_basis(2,false)
    res = FibonacciChain.joint_basis([1, 2])  
    # Ensureing the order is 2*1, not 1*2
    # kron(st1*st2') is in order 1*2, while reshape(st1*st2',9) is in order 2*1
    @test res == vec([join(l1, l2) for l1 in lis1 for l2 in lis2])
    @test res != [join(l2, l1) for l1 in lis1, l2 in lis2]


    # Test disjoint system joint_basis
    res = FibonacciChain.joint_basis([1], [2], anyon_typeA=:Fibo, anyon_typeB=:Fibo)
    @test res == [join(l1, l2) for l1 in anyon_basis(1) for l2 in anyon_basis(2,false)]
    res = FibonacciChain.joint_basis([1], [2], anyon_typeA=:IsingX, anyon_typeB=:IsingX)
    @test res == vec([join(l1, l2) for l1 in anyon_basis(1, anyon_type=:IsingX) for l2 in anyon_basis(2, anyon_type=:IsingX, false)])
    res = FibonacciChain.joint_basis([1], [2], anyon_typeA=:IsingX, anyon_typeB=:Fibo)
    @test res == vec([join(l1, l2) for l1 in anyon_basis(1, anyon_type=:IsingX) for l2 in anyon_basis(2, false)])
    res = FibonacciChain.joint_basis([1], [2], anyon_typeA=:Fibo, anyon_typeB=:IsingX)

    # move_subsystem
    res = FibonacciChain.move_subsystem(BitStr{5, Int}, BitStr{3, Int}(bit"111"), [1, 2, 5])
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
    st = ones(length(anyon_basis(N))); st /= norm(st)  # Normalize the state
    # The empty subsystem
    rdm = anyon_rdm(N, Int64[], st)
    @test rdm ≈ ones(Float64, 1,1)

    # The total system
    rdm = anyon_rdm(N, collect(1:N), st)
    @test rdm ≈ st*st'

    rdm = anyon_rdm(N, [1], st)
    @test rdm ≈ [0.75 0.25; 0.25 0.25]
end

@testset "anyon_rdm_Ising" begin
    N = 3
    st = zeros(2^N); st[1]=1; st[end]=1; st /= norm(st)  # Normalize the state
    rdm = anyon_rdm(N, Int[], st)
    @test rdm ≈ ones(Float64, 1,1)

    rdm = anyon_rdm(N, collect(1:N), st, anyon_type=:IsingX)
    @test rdm ≈ st*st'
    
    rdm = anyon_rdm(N, [1], st, anyon_type=:IsingX) ≈ 0.5*I(2)
end

@testset "anyon_rdm_matrix" begin
    N = 4
    st = zeros(length(anyon_basis(N)));st[5]=1; st[end]=1; st /= norm(st)  # Normalize the state

    rdm = anyon_rdm(N, collect(1:2), st*st')
    @test rdm ≈ diagm([0.0, 0.5, 0.5])
    rdm = anyon_rdm(N, collect(1:N), st*st')
    @test rdm ≈ st*st'

    st = zeros(2^N); st[1]=1; st[end]=1; st /= norm(st)  # Normalize the state

    rdm = anyon_rdm(N, [1,2], st*st', anyon_type=:IsingX) ≈ diagm([0.5, 0.0, 0.0, 0.5])

    st = zeros(2^2); st[1]=1; st[end]=1; st /= norm(st)  # 
    rdm = anyon_rdm(2, [1], st*st', anyon_type=:IsingX) ≈ 0.5 * I(2)
end

@testset "anyon_basis_K" begin
    N = 8
    T = BitStr{N}
    k = 0
    fib_basis_k = anyon_basis(T, k)[1]
    @test length(fib_basis_k) == 8  # Check the length of the basis
    @test [i.buf for i in fib_basis_k] == [0, 1, 5, 9, 17, 21, 37, 85]
end

@testset "anyon_ham_K" begin
    N = 8
    T = BitStr{N}
    k = 0
    fib_ham_k = anyon_ham(T, k)
    H = anyon_ham(N)
    @test size(fib_ham_k) == (8, 8)  # Check the size of the Hamiltonian matrix
    @test ishermitian(fib_ham_k)  # Check if the Hamiltonian is Hermitian
    @test eigvals(H)[1] ≈ eigvals(fib_ham_k)[1]  # Check if the ground state energy matches
end

@testset "mapst_sec2tot, anyon_rdm_sec" begin
    N = 8
    k = 0
    T = BitStr{N}
    fib_ham_k = anyon_ham(T, k)
    H = anyon_ham(N)
    
    sec_gs = eigvecs(fib_ham_k)[:, 1]
    gs = eigvecs(H)[:, 1]
    mapped_st = FibonacciChain.mapst_sec2tot(T, sec_gs, k)
    @test mapped_st ≈ gs  # Check if the mapped state matches the ground state

    rdm_sec = anyon_rdm_sec(N, collect(1:div(N,2)), sec_gs, k)
    rdm = anyon_rdm(N, collect(1:div(N,2)), gs)
    @test rdm_sec ≈ rdm  # Check if the reduced density matrix matches
end

@testset "disjoint_rdm" begin
    N1 = 4
    N2 = 4
    T1 = BitStr{N1}
    T2 = BitStr{N2}
    state = zeros(length(anyon_basis(N1)) * length(anyon_basis(N2))); state[1] = 1; state[end] = 1
    state = state ./ norm(state)  # Normalize the state
    subsystemsA = [1, 2]
    subsystemsB = [1, 2]
    
    rdm_result = disjoint_rdm(T1, T2, subsystemsA, subsystemsB, state)
    @test size(rdm_result) == (9, 9)  # Check the size of the reduced density matrix
    @test rdm_result[1,1] ≈ rdm_result[end, end] ≈ 0.5
    
    # Test for one subsystem is empty
    rdm_result_empty1 = disjoint_rdm(T1, T2, Int64[], subsystemsB, state)
    @test diag(rdm_result_empty1) ≈ [0.5, 0.0, 0.5]

    rdm_result_empty2 = disjoint_rdm(T1, T2, subsystemsA, Int64[], state)
    @test diag(rdm_result_empty2) ≈ [0.5, 0.0, 0.5]
end

@testset "disjoint_rdm_Ising" begin
    N1 = 4
    N2 = 4
    T1 = BitStr{N1}
    T2 = BitStr{N2}
    state = zeros(length(anyon_basis(N1, anyon_type = :IsingX)) * length(anyon_basis(N2, anyon_type=:IsingX))); state[1] = 1; state[end] = 1
    state = state ./ norm(state)  # Normalize the state
    subsystemsA = [1, 2]
    subsystemsB = [1, 2]
    
    rdm_result = disjoint_rdm(N1, N2, subsystemsA, subsystemsB, state, anyon_typeA=:IsingX, anyon_typeB=:IsingX)
    @test size(rdm_result) == (16, 16)  # Check the size of the reduced density matrix
    @test rdm_result[1,1] ≈ rdm_result[end, end] ≈ 0.5

    rdm_result_empty1 = disjoint_rdm(N1, N2, Int64[], subsystemsB, state, anyon_typeA=:IsingX, anyon_typeB=:IsingX)
    @test all(diag(rdm_result_empty1) ≈ [0.5, 0.0, 0.0, 0.5]) 

    rdm_result_empty2 = disjoint_rdm(N1, N2, subsystemsA, Int64[], state, anyon_typeA=:IsingX, anyon_typeB=:IsingX)
    @test all(diag(rdm_result_empty2) ≈ [0.5, 0.0, 0.0, 0.5])

end

@testset "disjoint_rdm_mixed" begin
    N1 = 4
    N2 = 4
    T1 = BitStr{N1}
    T2 = BitStr{N2}
    state = zeros(length(anyon_basis(N1)) * length(anyon_basis(N2, anyon_type = :IsingX))); state[1] = 1; state[end] = 1
    state = state ./ norm(state)  # Normalize the state
    subsystemsA = [1, 2]
    subsystemsB = [1, 2]
    
    rdm_result = disjoint_rdm(T1, T2, subsystemsA, subsystemsB, state, anyon_typeB=:IsingX)
    @test size(rdm_result) == (12, 12)  # Check the size of the reduced density matrix
    @test rdm_result[1,1] ≈ rdm_result[end, end] ≈ 0.5
    
    # Test for one subsystem is empty
    rdm_result_empty1 = disjoint_rdm(T1, T2, Int64[], subsystemsB, state, anyon_typeB=:IsingX)
    @test diag(rdm_result_empty1) ≈ [0.5, 0.0, 0.0, 0.5]

    rdm_result_empty2 = disjoint_rdm(T1, T2, subsystemsA, Int64[], state, anyon_typeB=:IsingX)
    @test diag(rdm_result_empty2) ≈ [0.5, 0.0, 0.5]
end

@testset "disjoint_rdm and ref qubit" begin
    L=6
    pbc=true
    anyon_type = :IsingX
    sample = ones(Int, 20, L)
    τ = log(1 + √2)
    initial_state = zeros(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
    initial_state[1] = 1.0
    statelis = generate_state(τ, initial_state, sample, temp= true, anyon_type=anyon_type)


    st = reference_evolution(τ, statelis, sample, 1, 10, 16, anyon_type=:IsingX)
    ρ1 = disjoint_rdm(2, L, Int64[], collect(1:div(L,2)), st, anyon_typeA=:IsingX, anyon_typeB=:IsingX)
    ρ2 = disjoint_rdm(2, L, Int64[], collect(1:L), st, anyon_typeA=:IsingX, anyon_typeB=:IsingX)
    ρ2r = anyon_rdm(L, collect(1:div(L,2)), ρ2, anyon_type=:IsingX)
    S1 = ee(ρ1)
    S2 = ee(ρ2r)
    @test S1 ≈ S2
end