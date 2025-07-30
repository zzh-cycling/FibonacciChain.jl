using FibonacciChain
using Test 
using SparseArrays

@testset "FibonacciChain Tests" begin
    N = 8
    H = Fibonacci_Ham(N)
    H_sparse = Fibonacci_Ham_sparse(N)

    @test Matrix(H_sparse) ≈ H
end

