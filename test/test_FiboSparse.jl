using FibonacciChain
using Test 
using SparseArrays
using Arpack

@testset "FibonacciChain Tests" begin
    N = 8
    H = anyon_ham(N)
    H_sparse = anyon_ham_sparse(N)

    @test Matrix(H_sparse) ≈ H
end

@testset "anyon_ham_sparse_K" begin
    N = 8
    H = anyon_ham(N)
    H_sparse = anyon_ham_sparse(N, 0)
    energy = eigvals(H)[1]
    energy_sparse, state_sparse = Arpack.eigs(H_sparse, nev=1, which=:SR)
    @test energy ≈ energy_sparse[1]
end

@testset "anyon_ham_sparse_J" begin
    N = 8
    H = anyon_ham_sparse(N, J=1.0, h=1.0, anyon_type=:IsingX)
    energy = Arpack.eigs(H, nev=1, which=:SR)[1][1]
    @test energy ≈ -1/(2*sin(π/2N))*4
end
