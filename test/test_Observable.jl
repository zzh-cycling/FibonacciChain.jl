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

@testset "braidingsq_basismapN3" begin
    N=3
    T = BitStr{N, Int}

    state=T(bit"010")
    ϕ = (1+√5)/2
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"010"),exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"010"), T(bit"000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"010"), exp(-6im*π/5))

    state =T(bit"001")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"001"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"001"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"001"), T(bit"000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    state =T(bit"100")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"100"), T(bit"000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"100"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"100"), exp(-6im*π/5))

    state =T(bit"101") # Not in PBC basis
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"101"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2, false) == (T(bit"101"), exp(-2im*π/5))

    state =T(bit"000")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"000"), T(bit"100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"000"), T(bit"010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"000"), T(bit"001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
end

@testset "braidingsq_basismapN6" begin
    N=6
    T = BitStr{N, Int}
    ϕ = (1+√5)/2
    ## state 1 ###
    state =T(bit"000000")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"000000"), T(bit"100000"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"000000"), T(bit"010000"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"000000"), T(bit"001000"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"000000"), T(bit"000100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"000000"), T(bit"000010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"000000"), T(bit"000001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 2 ###
    state =T(bit"000001")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"000001"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"000001"), T(bit"010001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"000001"), T(bit"001001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"000001"),  T(bit"000101"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"000001"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"000001"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 3 ###
    state=T(bit"000010")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"000010"), T(bit"100010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"000010"), T(bit"010010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"000010"), T(bit"001010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"000010"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"000010"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"000010"), exp(-6im*π/5))

    ## state 4 ###
    state =T(bit"000100")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"000100"), T(bit"100100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"000100"), T(bit"010100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"000100"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"000100"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"000100"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"000100"), T(bit"000101"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 5 ###
    state =T(bit"000101")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"000101"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"000101"), T(bit"010101"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"000101"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"000101"), T(bit"000001"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"000101"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"000101"), T(bit"000100"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    
    ## state 6 ###
    state =T(bit"001000")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"001000"), T(bit"101000"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"001000"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"001000"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"001000"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"001000"), T(bit"001010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"001000"), T(bit"001001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 7 ###
    state =T(bit"001001")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"001001"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"001001"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"001001"), T(bit"000001"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"001001"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"001001"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"001001"), T(bit"001000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 8 ###
    state =T(bit"001010")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"001010"), T(bit"101010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"001010"),exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"001010"), T(bit"000010"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"001010"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"001010"), T(bit"001000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"001010"), exp(-6im*π/5))

    ## state 9 ###
    state =T(bit"010000")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"010000"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"010000"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"010000"),exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"010000"), T(bit"010100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"010000"), T(bit"010010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"010000"),(bit"010001"),  exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 10 ###
    state =T(bit"010001")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"010001"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"010001"), T(bit"000001"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"010001"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"010001"), T(bit"010101"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"010001"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"010001"), T(bit"010000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 11 ###
    state =T(bit"010010")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"010010"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"010010"), T(bit"000010"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"010010"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"010010"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"010010"), T(bit"010000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"010010"),exp(-6im*π/5))

    ## state 12 ###
    state =T(bit"010100")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"010100"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"010100"), T(bit"000100"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"010100"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"010100"), T(bit"010000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"010100"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"010100"), T(bit"010101"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    
    ## state 13 ###
    state =T(bit"010101")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"010101"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"010101"), T(bit"000101"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"010101"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"010101"), T(bit"010001"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"010101"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"010101"), T(bit"010100"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 14 ###
    state =T(bit"100000")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"100000"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"100000"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"100000"), T(bit"101000"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"100000"), T(bit"100100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"100000"), T(bit"100010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"100000"), exp(-6im*π/5))

    ## state 15 ###
    state =T(bit"100010")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"100010"), T(bit"000010"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"100010"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"100010"), T(bit"101010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"100010"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"100010"), T(bit"100000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"100010"), exp(-2im*π/5))

    ## state 16 ###
    state =T(bit"100100")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"100100"), T(bit"000100"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"100100"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"100100"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"100100"), T(bit"100000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"100100"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"100100"),exp(-6im*π/5))

    ## state 17 ###
    state =T(bit"101000")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"101000"), T(bit"001000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"101000"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"101000"), T(bit"100000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"101000"), exp(-6im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"101000"), T(bit"101010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"101000"),exp(-6im*π/5))

    ## state 18 ###
    state =T(bit"101010")
    @test FibonacciChain.braidingsq_basismap(T, state, 1) == (T(bit"101010"), T(bit"001010"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 2) == (T(bit"101010"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 3) == (T(bit"101010"), T(bit"100010"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 4) == (T(bit"101010"), exp(-2im*π/5))
    @test FibonacciChain.braidingsq_basismap(T, state, 5) == (T(bit"101010"), T(bit"101000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain.braidingsq_basismap(T, state, 6) == (T(bit"101010"), exp(-2im*π/5))

end

@testset "braidingsq_matrix" begin
    N=3
    T = BitStr{N, Int}
    ϕ = (1+√5)/2
    @test FibonacciChain.braidingsq_matrix(T, 2, false) ≈  ComplexF64[
    exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2) 0.0 (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2) 0.0 0.0;
    0.0 exp(-6im*π/5) 0.0 0.0 0.0; 
    (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2) 0.0 exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1) 0.0 0.0; 
    0.0 0.0 0.0 exp(-6im*π/5) 0.0; 0.0 0.0 0.0 0.0 exp(-2im*π/5)]

    B=[exp(4*π*im/5) 0 0 0 0;
    0 exp(-3*π*im/5) 0 0 0;
    0 0 exp(-3*π*im/5) 0 0;
    0 0 0 exp(-4*π*im/5)*ϕ^(-1) -exp(-2*π*im/5)*ϕ^(-1/2);
    0 0 0 -exp(-2*π*im/5)*ϕ^(-1/2) -ϕ^(-1)]

    U = [0 0 0 0 1;
        0 1 0 0 0;
        0 0 0 1 0;
        0 0 1 0 0;
        1 0 0 0 0]

    @test FibonacciChain.braidingsq_matrix(T, 2, false) ≈ U' * B^2 * U
    @test FibonacciChain.braidingsq_matrix(T, 2) ≈ (U' * B^2 * U)[1:4, 1:4]

    @test FibonacciChain.braidingsq_matrix(T, 1) ≈ ComplexF64[
    exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2) 0.0 0.0 (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2);
    0.0 exp(-6im*π/5) 0.0 0.0; 
    0.0 0.0  exp(-6im*π/5) 0.0; 
    (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2) 0.0 0.0 exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)]

    @test FibonacciChain.braidingsq_matrix(T, 3) ≈ ComplexF64[
    exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2) (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2) 0.0 0.0;
    (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2) exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1) 0.0 0.0;
    0.0 0.0 exp(-6im*π/5) 0.0; 
    0.0 0.0 0.0 exp(-6im*π/5)]

    ⊗(A,B) = kron(A,B)
    idx = [i.buf+1 for i in Fibonacci_basis(3,false)]
    Z=[1 0;0 -1]
    X=[0 1;1 0]
    P0=[1 0;0 0]
    P1=[0 0;0 1]
    σ2t = (1-2ϕ^(-1)) * Z +2ϕ^(-3/2) * X
    σt = ϕ^(-1)*Z +2ϕ^(-1/2) * X

    # F = P0 ⊗ I(2) ⊗ P1 + P1 ⊗ I(2) ⊗ P0 + P1 ⊗ X ⊗ P1 + P0 ⊗ σt ⊗ P0
    # Fc = U'*F[idx, idx]*U
    # B1 = exp(1im*π/10)*(P0 ⊗ exp(-7im*π/10 .* Z) ⊗ P1 + P1 ⊗ exp(-7im*π/10 .* Z) ⊗ P0 + P1 ⊗ exp(+7im*π/10 .* Z) ⊗ P1 + P0 ⊗ exp(-7im*π/10 .* σ2t) ⊗ P0)
    # Bc=B1[idx, idx]
    #  (U'*(exp(1im*π/10)*F'*(I(2)⊗exp(-7im*π/10 .* Z)⊗I(2))*F)[idx,idx])*U
    # exp(1im*π/10)*F'*(I(2)⊗exp(-7im*π/10 .* Z)⊗I(2))*F ≈ Bc
    # @test U'*Bc^2*U ≈ FibonacciChain.braidingsq_matrix(T, 2, false)
    # @test U'*Bc*U ≈ B
    # Bc^2
end

@testset "braidingsqmap" begin
    N=3
    state = collect(1:4)
    ϕ = (1+√5)/2
    @test braidingsqmap(N, state, 2) ≈ [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+3(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+3(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 4exp(-6im*π/5)]
    @test braidingsqmap(N, state, 1) ≈ [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+4(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), 3exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+4(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1))]
    @test braidingsqmap(N, state, 3) ≈ [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+2(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+2(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 3exp(-6im*π/5), 4exp(-6im*π/5)]
end


