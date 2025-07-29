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
    EE_lis=eelis_Fibo_state(N,state)
    @test length(EE_lis)==length(collect(1:N-1))
    @test all(EE_lis .> 0)
end

@testset "ee_Fiboladder_lis" begin
    N=3
    state=eigvecs(Fibonacci_Ham(N))[:,1]
    EE_lis=eelis_Fiboladder_state(N, kron(state, state))
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
    # Test braidingsq_matrix = the formulas of B^2 or not
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

    # Test that the braidingsq_matrix is equal to U' * B^2 * U
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
    σ2t = (1-2ϕ^(-1)) * Z - 2ϕ^(-3/2) * X
    σt = -ϕ^(-1)*Z +ϕ^(-1/2) * X

    pesudoF = [
        0 0 0 0 0
        0 1 0 0 0;
        0 0 1 0 0;
        0 0 0 ϕ^(-1) ϕ^(-1/2);
        0 0 0 ϕ^(-1/2) -ϕ^(-1)]
    F = [
        1 0 0 0 0
        0 1 0 0 0;
        0 0 1 0 0;
        0 0 0 ϕ^(-1) ϕ^(-1/2);
        0 0 0 ϕ^(-1/2) -ϕ^(-1)]
    # Zhu's formulas
    Fz = P0 ⊗ I(2) ⊗ P1 + P1 ⊗ I(2) ⊗ P0 + P1 ⊗ X ⊗ P1 + P0 ⊗ σt ⊗ P0
    Fc = U'*Fz[idx, idx]*U
    @test pesudoF ≈ Fc
    Bz = exp(1im*π/10)*(P0 ⊗ exp(-7im*π/10 .* Z) ⊗ P1 + P1 ⊗ exp(-7im*π/10 .* Z) ⊗ P0 + P1 ⊗ exp(+7im*π/10 .* Z) ⊗ P1 + P0 ⊗ exp(-7im*π/10 .* σ2t) ⊗ P0)
    Bc=U'*Bz[idx, idx]*U
    @test Bc ≈ B

end

@testset "braidingsqmap" begin
    N=3
    state = collect(1:4)
    ϕ = (1+√5)/2
    @test braidingsqmap(N, state, 2) ≈ [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+3(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+3(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 4exp(-6im*π/5)]
    @test braidingsqmap(N, state, 1) ≈ [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+4(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), 3exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+4(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1))]
    @test braidingsqmap(N, state, 3) ≈ [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+2(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+2(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 3exp(-6im*π/5), 4exp(-6im*π/5)]
end

@testset "build_extended_basis" begin
    N = 3
    k_total = 2
    T = BitStr{N+k_total, Int}
    basis = Fibonacci_basis(N)
    extended_basis = FibonacciChain.build_extended_basis(k_total,basis)

    @test length(extended_basis) == 4 * length(basis)
    @test extended_basis == T.([0, 1, 2, 4, 8, 9, 10, 12, 16, 17, 18, 20, 24, 25, 26, 28])

    k_total = 1
    T = BitStr{N+k_total, Int}
    extended_basis = FibonacciChain.build_extended_basis(k_total, basis)

    @test length(extended_basis) == 2 * length(basis)
    @test extended_basis == T.([0, 1, 2, 4, 8, 9, 10, 12])

    k_total = 0
    extended_basis = FibonacciChain.build_extended_basis(k_total, basis)
    @test length(extended_basis) == length(basis)
    @test extended_basis == basis
end

@testset "add_reference_qubits" begin
    N = 3
    st =ones(4)/2; 
    # Test add_reference_qubit
    add_st1 = FibonacciChain.add_reference_qubits!(N, st, 1)
    add_st2 = FibonacciChain.add_reference_qubits!(N, st, 2)
    add_st3 = FibonacciChain.add_reference_qubits!(N, st, 3)
    @test add_st1 == [0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5]
    @test add_st2 == [0.5, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0]
    @test add_st3 == [0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.0]

    add2_st1 = FibonacciChain.add_reference_qubits!(N, add_st1, 1)
    add2_st2 = FibonacciChain.add_reference_qubits!(N, add_st2, 2)
    add2_st3 = FibonacciChain.add_reference_qubits!(N, add_st3, 3)

    @test add2_st1 == [0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5]
    @test add2_st2 == [0.5, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0]
    @test add2_st3 == [0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0]

end

@testset "spatial_correlation" begin
    N=12
    mes = zeros(length(Fibonacci_basis(12, true, measure_class=:Fibo)))
    mes[233] = 1/√2 
    mes[end] = 1/√2 # set the last two qubits to be in the Bell state

    sclis = [spatial_correlation(N, mes, i, j, pbc) for i in 1:N-1, j in i+1:N]
    @test sclis ≈ log(2)*ones(11,11)
end

@testset "reference_measure_basismap" begin
    N = 3
    τ = 1.0
    sign = 0
    pbc = true
    k_old = 1
    T = BitStr{N, Int}
    basislis = Fibonacci_basis(N, pbc, measure_class=:Fibo)
    l = length(basislis)
    output11 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 1, sign, pbc, k_old=0)
    output12 = FibonacciChain.measure_basismap.(T, τ, basislis, 1, sign, pbc)
    output21 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 2, sign, pbc, k_old=0)
    output22 = FibonacciChain.measure_basismap.(T, τ, basislis, 2, sign, pbc)
    output31 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 3, sign, pbc, k_old=0)
    output32 = FibonacciChain.measure_basismap.(T, τ, basislis, 3, sign, pbc)
    @test all([all(output11[i] .≈ output12[i]) for i in 1:l])
    @test all([all(output21[i] .≈ output22[i]) for i in 1:l])
    @test all([all(output31[i] .≈ output32[i]) for i in 1:l])

    extended_basis = FibonacciChain.build_extended_basis(1, basislis)
    output13 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 1, sign, pbc, k_old=1)
    output23 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 2, sign, pbc, k_old=1)
    output33 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 3, sign, pbc, k_old=1)

    @test all([all([all(output13[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output23[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output33[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:1])

    extended_basis2 = FibonacciChain.build_extended_basis(2, basislis)
    output14 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 1, sign, pbc, k_old=2)
    output24 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 2, sign, pbc, k_old=2)
    output34 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 3, sign, pbc, k_old=2)

    @test all([all([all(output14[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output24[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output34[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:3])
end

@testset "reference_measuremap" begin
    N = 3
    τ = 1.0
    sign = 0
    pbc = true
    k_old = 1
    T = BitStr{N, Int}
    st = ones(4)/2;
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1)

    output13 = FibonacciChain.reference_measuremap(T, τ, add_st, 1, sign, pbc, k_old=1)
    output23 = FibonacciChain.reference_measuremap(T, τ, add_st, 2, sign, pbc, k_old=1)
    output33 = FibonacciChain.reference_measuremap(T, τ, add_st, 3, sign, pbc, k_old=1)
    @test output13 == [0.2859295753144778, 0.4692539498975694, 0.4692539498975694, -0.14412070965501542, -0.14412070965501542, 0.0, 0.0, 0.35595325543890144]
    @test output23 == [0.1418088656594624, 0.4692539498975694, 0.21183254578388602, 0.0, 0.0, 0.0, 0.0, 0.4692539498975694]
    @test output33 == [0.1418088656594624, 0.21183254578388602, 0.4692539498975694, 0.0, 0.0, 0.0, 0.0, 0.4692539498975694]
end

@testset "reference_apply_measurement_layer" begin
    N=8
    τ = 1000.0
    sign = 1
    pbc = true
    k_old = 1
    st = zeros(length(Fibonacci_basis(N, pbc))); st[1] = 1
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1)
    sample = zeros(Int, length(2:2:N))

    output1 = FibonacciChain.reference_apply_measurement_layer!(N, add_st, τ, sample, 1, pbc, k_old=1)
    output2 = FibonacciChain.apply_measurement_layer!(N, st, τ, sample, 1, pbc)
    @test output1[1:div(length(output1), 2)] == output2

    output3 = FibonacciChain.reference_apply_measurement_layer!(N, add_st, τ, sample, 2, pbc, k_old=1)
    output4 = FibonacciChain.apply_measurement_layer!(N, st, τ, sample, 2, pbc)
    @test output3[1:div(length(output3), 2)] == output4
    output5 = FibonacciChain.reference_apply_measurement_layer!(N, add_st, τ, sample, 3, pbc, k_old=1)
    output6 = FibonacciChain.apply_measurement_layer!(N, st, τ, sample, 3, pbc)
    @test output5[1:div(length(output5), 2)] == output6

    add_st2 = FibonacciChain.add_reference_qubits!(N, add_st, 1)
    output7 = FibonacciChain.reference_apply_measurement_layer!(N, add_st2, τ, sample, 4, pbc, k_old=2)
    output8 = FibonacciChain.apply_measurement_layer!(N, st, τ, sample, 4, pbc)
    @test output7[1:div(length(output7),4)] == output8
end

@testset "reference_generate_state" begin
    N = 8
    τ = 1000.0
    sign = 1
    pbc = true
    k_old = 1
    st = zeros(length(Fibonacci_basis(N, pbc))); st[1] = 1
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1)
    sample = zeros(Int, (3, length(2:2:N)))

    output1 = reference_generate_state(τ, add_st, sample, pbc, k_old=1)
    output2 = generate_state(τ, st, sample, pbc, temp= true)
    @test [i[1:47] for i in output1] == output2

    output3 = reference_generate_state(τ, add_st, sample, pbc, k_old=1)
    output4 = generate_state(τ, st, sample, pbc, temp= true)
    @test [i[1:47] for i in output3] == output4

    output5 = reference_generate_state(τ, add_st, sample, pbc, k_old=1)
    output6 = generate_state(τ, st, sample, pbc, temp= true)
    @test [i[1:47] for i in output5] == output6
end

@testset "temporal_correlation" begin
    N=12
    mes = zeros(length(Fibonacci_basis(12, true, measure_class=:Fibo)))
    mes[233] = 1/√2 
    mes[end] = 1/√2 # set the last two qubits to be in the Bell state
    
    sample = [0 1 1 1 1 1; 1 0 0 0 0 1; 1 1 1 1 1 1; 0 0 1 1 1 1; 1 0 1 1 1 1; 0 1 1 0 0 1; 1 1 1 1 1 1; 0 1 1 1 0 0; 0 1 1 1 1 0; 1 0 1 1 0 1]
    D= 10 
    # sample = zeros(Int, D, length(2:2:N))
    tc = temporal_correlation(τ, mes, sample, div(N,2), 5, 8)
    
    # tclis = [temporal_correlation(τ, mes, sample, div(N,2), i, j) for i in 1:D-1 for j in i+1:D]
    @test tc ≈ 0.03098628964295691
end