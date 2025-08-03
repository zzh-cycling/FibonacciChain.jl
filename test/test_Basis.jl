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
    @test isapprox(Y, Y')
    vals = eigvals(Y)
    @test vals == [-1.0799610383969367, -0.23606797749978994, -0.23606797749978972, -0.23606797749978964, 0.4778136965285674, 1.9041523147215358, 3.079961038396939]
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
    fib_basis = Fibonacci_basis(5)
    @test length(fib_basis) == 11
    fib_basis = Fibonacci_basis(5,false)
    @test length(fib_basis) == 13
    # Test the Fibonacci Hamiltonian
    fib_ham = Fibonacci_Ham(5)
    @test size(fib_ham) == (11, 11)
    @test ishermitian(fib_ham)

    @test Fibonacci_Ham(3,false) == [-0.6180339887498948 0.0 -0.48586827175664565 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; -0.48586827175664565 0.0 -0.3819660112501051 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0]
    @test Fibonacci_Ham(3) == [-1.8541019662496843 -0.48586827175664565 -0.48586827175664565 -0.48586827175664565; -0.48586827175664565 -0.3819660112501051 0.0 0.0; -0.48586827175664565 0.0 -0.3819660112501051 0.0; -0.48586827175664565 0.0 0.0 -0.3819660112501051]
    # Test the reduced density matrix function
    # rdm = FibonacciChain.rdm_Fibo(fib_basis, 2)
    # @test size(rdm) == (2, 2)

end

@testset "process_join" begin
    # [[000 ₍₂₎, 001 ₍₂₎, 010 ₍₂₎, 100 ₍₂₎, 101 ₍₂₎], [000 ₍₂₎, 001 ₍₂₎, 010 ₍₂₎, 100 ₍₂₎, 101 ₍₂₎]]
    lis1 = BitStr{2}[0, 1, 2]
    lis2 = BitStr{3}[0, 1, 2, 4, 5]
    res = FibonacciChain.process_join(lis1, lis2) 
    @test res == vec([join(l2, l1) for l1 in lis1, l2 in lis2])

    # joint_basis
    res = FibonacciChain.joint_basis([2, 3])
    @test res == vec([join(l2, l1) for l1 in lis1, l2 in lis2])
    
    lis1 = Fibonacci_basis(1,false)
    lis2 = Fibonacci_basis(2,false)
    res = FibonacciChain.joint_basis([1, 2])  
    @test res == vec([join(l2, l1) for l1 in lis1, l2 in lis2])
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

@testset "reference_rdm" begin
    N = 3
    st = ones(4)/2;
    site = 1
    add_site1 = FibonacciChain.add_reference_qubits!(N, st, site)
    
    rdm = reference_rdm(N, [1], add_site1)
    @test rdm == [0.75 0.0; 0.0 0.25]

    full_st = zeros(2^(N+1))
    inds = [1, 3, 5, 10]
    full_st[inds] .= 0.5
    rdm_Fibo(4, [1], full_st, measure_class=:IsingX) == [0.75 0.0; 0.0 0.25]

    add_st2 = FibonacciChain.add_reference_qubits!(N, add_site1, site)
    rdm2 = reference_rdm(N, [1], add_st2)
    rdm3 = reference_rdm(N, [2], add_st2)
    @test rdm == rdm2
    @test rdm2 == rdm3

    rdm4 = reference_rdm(N, [1, 2], add_st2)
    @test rdm4 == [0.75 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.25]
end

@testset "reference_rdm_Ising" begin
    N = 3
    st = ones(2^N); st /= norm(st)  # Normalize the state
    site = 1
    measure_class = :IsingX
    add_site1 = FibonacciChain.add_reference_qubits!(N, st, site, measure_class = measure_class)
    
    rdm = reference_rdm(N, [1], add_site1, measure_class = measure_class)
    @test rdm ≈ [0.5 0.0; 0.0 0.5]

    add_st2 = FibonacciChain.add_reference_qubits!(N, add_site1, site, measure_class = measure_class)
    rdm2 = reference_rdm(N, [1], add_st2, measure_class = measure_class)
    rdm3 = reference_rdm(N, [2], add_st2, measure_class = measure_class)
    @test rdm == rdm2
    @test rdm2 == rdm3

    rdm4 = reference_rdm(N, [1, 2], add_st2, measure_class = measure_class)
    @test rdm4 ≈ [0.4999999999999999 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.4999999999999999]
end

@testset "Fibonacci_basis_K" begin
    N = 8
    T = BitStr{N}
    k = 0
    fib_basis_k = Fibonacci_basis(T, k)[1]
    @test length(fib_basis_k) == 8  # Check the length of the basis
    @test [i.buf for i in fib_basis_k] == [0, 1, 5, 9, 17, 21, 37, 85]
end

@testset "Fibonacci_Ham_K" begin
    N = 8
    T = BitStr{N}
    k = 0
    fib_ham_k = Fibonacci_Ham(T, k)
    H = Fibonacci_Ham(N)
    @test size(fib_ham_k) == (8, 8)  # Check the size of the Hamiltonian matrix
    @test ishermitian(fib_ham_k)  # Check if the Hamiltonian is Hermitian
    @test eigvals(H)[1] ≈ eigvals(fib_ham_k)[1]  # Check if the ground state energy matches
end

@testset "mapst_sec2tot, rdm_Fibo_sec" begin
    N = 8
    k = 0
    T = BitStr{N}
    fib_ham_k = Fibonacci_Ham(T, k)
    H = Fibonacci_Ham(N)
    
    sec_gs = eigvecs(fib_ham_k)[:, 1]
    gs = eigvecs(H)[:, 1]
    mapped_st = FibonacciChain.mapst_sec2tot(T, sec_gs, k)
    @test mapped_st ≈ gs  # Check if the mapped state matches the ground state

    rdm_sec = rdm_Fibo_sec(N, collect(1:div(N,2)), sec_gs, k)
    rdm = rdm_Fibo(N, collect(1:div(N,2)), gs)
    @test rdm_sec ≈ rdm  # Check if the reduced density matrix matches
end

@testset "disjoint_rdm" begin
    N1 = 4
    N2 = 4
    T1 = BitStr{N1}
    T2 = BitStr{N2}
    state = zeros(length(Fibonacci_basis(N1)) * length(Fibonacci_basis(N2))); state[1] = 1; state[end] = 1
    state = state ./ norm(state)  # Normalize the state
    subsystemsA = [1, 2]
    subsystemsB = [1, 2]
    
    rdm_result = disjoint_rdm(T1, T2, subsystemsA, subsystemsB, state)
    @test size(rdm_result) == (9, 9)  # Check the size of the reduced density matrix
    @test rdm_result[1,1] ≈ rdm_result[end, end] ≈ 0.5
end