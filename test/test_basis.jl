using FibonacciChain
using Test

@testset "basis.jl" begin
    # Test the Fibonacci basis creation
    fib_basis = FibonacciChain.Fibonacci_basis(5)
    @test length(fib_basis) == 5

    # Test the Fibonacci Hamiltonian
    fib_ham = FibonacciChain.Fibonacci_Ham(5)
    @test size(fib_ham) == (5, 5)

    # Test the reduced density matrix function
    rdm = FibonacciChain.rdm_Fibo(fib_basis, 2)
    @test size(rdm) == (2, 2)

    # Additional tests can be added here
end