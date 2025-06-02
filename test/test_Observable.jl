using FibonacciChain
using Test
using BitBasis
using LinearAlgebra 

@testset "ee" begin
    N=6
    state=eigvecs(Fibonacci_Ham(N))[:,1]
    rdm=rdm_Fibo(N, collect(1:div(N,2)), state)
    @test size(rdm)==(5,5)
    @test isapprox(FibonacciChain.ee(rdm), 0.7619577865215983
    , atol=1e-5)
end

@testset "eelis" begin
    N=6
    state=eigvecs(Fibonacci_Ham(N))[:,1]
    EE_lis=eelis_Fibo_state(N,collect(1:N-1),state)
    @test length(EE_lis)==length(collect(1:N-1))
    @test all(EE_lis .> 0)
end

@testset "translation_matrix" begin
    N=8
    Mat=translation_matrix(N)
    @test size(Mat)==(47,47)
    @test isapprox(Mat*Mat',I,atol=1e-5) # Check if the matrix is unitary
    @test isapprox(Mat^N,I,atol=1e-5)
end

@testset "inversion_matrix" begin
    N=8
    Mat=inversion_matrix(N)
    @test size(Mat)==(47,47)
    @test isapprox(Mat*Mat',I,atol=1e-5) # Check if the matrix is unitary
    @test isapprox(Mat*Mat,I,atol=1e-5) 
end 

@testset "braiding_basismap" begin
    N=3
    T = BitStr{N, Int}

    state=T(bit"010")
    @test braiding_basismap(T, state, 1) == (T(bit"010"),exp(-6im*π/5))
    @test braiding_basismap(T, state, 2) == (T(bit"010"), T(bit"000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test braiding_basismap(T, state, 3) == (T(bit"010"), exp(-6im*π/5))

    state =T(bit"001")
    @test braiding_basismap(T, state, 1) == (T(bit"001"), exp(-6im*π/5))
    @test braiding_basismap(T, state, 2) == (T(bit"001"), exp(-6im*π/5))
    @test braiding_basismap(T, state, 3) == (T(bit"001"), T(bit"000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    state =T(bit"100")
    @test braiding_basismap(T, state, 1) == (T(bit"100"), T(bit"000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test braiding_basismap(T, state, 2) == (T(bit"100"), exp(-6im*π/5))
    @test braiding_basismap(T, state, 3) == (T(bit"100"), exp(-6im*π/5))

    state =T(bit"101") # Not in PBC basis
    @test braiding_basismap(T, state, 2) == (T(bit"101"), exp(-2im*π/5))
    @test braiding_basismap(T, state, 2, false) == (T(bit"101"), exp(-2im*π/5))

    state =T(bit"000")
    @test braiding_basismap(T, state, 1) == (T(bit"000"), T(bit"100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test braiding_basismap(T, state, 2) == (T(bit"000"), T(bit"010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test braiding_basismap(T, state, 3) == (T(bit"000"), T(bit"001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
end

@testset "braiding_matrix" begin
    N=3
    T = BitStr{N, Int}
    @test braiding_matrix(T, 2, false) ≈  ComplexF64[
    exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2) 0.0 (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2) 0.0 0.0;
    0.0 exp(-6im*π/5) 0.0 0.0 0.0; 
    (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2) 0.0 exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1) 0.0 0.0; 
    0.0 0.0 0.0 exp(-6im*π/5) 0.0; 0.0 0.0 0.0 0.0 exp(-2im*π/5)]
end

