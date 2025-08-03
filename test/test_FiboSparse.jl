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

@testset "iso_total2sec_sparse" begin
    N = 8
    k = 0
    iso_sparse = iso_total2sec_sparse(N, k)
    @test size(iso_sparse) == (Fibonacci_basis(N, k), Fibonacci_basis(N, 0))
    
    # Check if the isometry is correct
    basis = Fibonacci_basis(N, k)
    for state in basis
        @test norm(iso_sparse * state) ≈ 1.0
    end
end