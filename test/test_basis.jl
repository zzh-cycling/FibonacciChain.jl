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

    # joint_pxp_basis
    res = FibonacciChain.joint_Fibo_basis([2, 3])
    @test res == vec([join(l2, l1) for l1 in lis1, l2 in lis2])

    # move_subsystem
    res = FibonacciChain.move_subsystem(BitStr{5, Int}, BitStr{3, Int}(0b101), [1, 2, 5])
    @test res == BitStr{5}(0b10001)

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

@testset "ladderrdm" begin
    # Test the ladder_rdm function
    N=2
    T = BitStr{N, Int}
    state = collect(1:3)
    vec = reshape(state*state', 9)
    rdm = ladderrdm(T, Int[1], vec)
    @test size(rdm) == (4, 4)
    @test rdm == [100 20 20 4; 20 40 4 8; 20 4 40 8; 4 8 8 16]

    # Test the ladder_rdm with a different basis
    vec/=norm(vec)
    rdm2 = FibonacciChain.ladderrdm(T, Int[1], vec)
    @test ishermitian(rdm2)
    @test tr(rdm2) ≈ 1.0 atol=1e-6
    @test sum(eigvals(rdm2)) ≈ 1.0 atol=1e-6
end

@testset "ladderChoi" begin
    N=3
    T = BitStr{N, Int}
    state = collect(1:4)
    state = reshape(state*state', 4^2)
    @test Float64.(ladderChoi(N, 0.0, state)) ≈ state/norm(state)
    
    ϕ = (1+√5)/2
    onechain_state = [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+3(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+3(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 4exp(-6im*π/5)]
    st = reshape(onechain_state*transpose(onechain_state), 16)
    @test ladderChoi(N, 1.0, state) ≈ st/norm(st)
 
    # At least for N = 4, the ladder Choi state is invariant under the ladder translation map twice.
    N=4
    energy, states = eigen(Fibonacci_Ham(N))
    antiGS= states[:, 1]
    len= length(antiGS)
    vecGS = reshape(antiGS*antiGS', len^2)
    plis = collect(0.0:0.1:1.0)
    for p in plis
        state = ladderChoi(N, p, vecGS)
        @show p
        @test isapprox(reduce((x, _) -> laddertranslationmap(N, x), 1:2; init=state),state, atol=1e-10)
        # @show isapprox(laddertranslationmap(N, laddertranslationmap(N, state)), state, atol=1e-10)
    end
    
    N=6
    energy, states = eigen(Fibonacci_Ham(N))
    antiGS= states[:, 1]
    len= length(antiGS)
    state = reshape(antiGS*antiGS', len^2)

    p=0.5
    for i in 2:2:N
        state=(1-p)*state+p*ladderbraidingmap(N, state, i)
        state/=norm(state) # normalize the state after each braiding
    end
end