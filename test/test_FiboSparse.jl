using FibonacciChain
using Test 
using SparseArrays
using Arpack

@testset "FibonacciChain Tests" begin
    N = 8
    H = Fibonacci_Ham(N)
    H_sparse = Fibonacci_Ham_sparse(N)

    @test Matrix(H_sparse) ≈ H
end

@testset "Fibonacci_Ham_sparse_K" begin
    N = 8
    H = Fibonacci_Ham(N)
    H_sparse = Fibonacci_Ham_sparse(N, 0)
    energy = eigvals(H)[1]
    energy_sparse, state_sparse = eigs(H_sparse, nev=1, which=:SR)
    @test energy ≈ energy_sparse[1]
end