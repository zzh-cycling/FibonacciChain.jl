using FibonacciChain
using Test
using LinearAlgebra 

@testset "Qmap" begin
    ϕ = (1+√5)/2
    @test FibonacciChain.Qmap(BitStr{3}, bit"000", 2) == (bit"000", bit"010", 1/ϕ, -1/ϕ^(3/2))
    @test FibonacciChain.Qmap(BitStr{3}, bit"010", 2) == (bit"010", bit"000", 1/ϕ^2, -1/ϕ^(3/2))
end

@testset "actingHam" begin
    ϕ = (1+√5)/2
    states, weights = FibonacciChain.actingHam(BitStr{3}, bit"000") 
    @test states== BitStr{3}.([bit"000", bit"010", bit"000", bit"100", bit"000", bit"001"])
    @test weights ≈ [1/ϕ, -1/ϕ^(3/2), 1/ϕ, -1/ϕ^(3/2), 1/ϕ, -1/ϕ^(3/2)]
    states, weights = FibonacciChain.actingHam(BitStr{3}, bit"010")
    @test states == BitStr{3}.([bit"010", bit"000"])
    @test weights ≈  [1/ϕ^2, -1/ϕ^(3/2)]
end

@testset "basis.jl" begin
    # Test the Fibonacci basis creation
    fib_basis = FibonacciChain.Fibonacci_basis(5)
    @test length(fib_basis) == 11
    fib_basis = FibonacciChain.Fibonacci_basis(5,false)
    @test length(fib_basis) == 13
    # Test the Fibonacci Hamiltonian
    fib_ham = FibonacciChain.Fibonacci_Ham(5)
    @test size(fib_ham) == (11, 11)
    @test ishermitian(fib_ham)

    # Test the reduced density matrix function
    # rdm = FibonacciChain.rdm_Fibo(fib_basis, 2)
    # @test size(rdm) == (2, 2)

    # Additional tests can be added here
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