using FibonacciChain
using Test
using BitBasis
using LinearAlgebra 
using Random

@testset "ee" begin
    N=6
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    state=eigvecs(anyon_ham(model))[:,1]
    rdm=anyon_rdm(model, collect(1:div(N,2)), state)
    @test size(rdm)==(5,5)
    @test isapprox(FibonacciChain.ee(rdm), 0.7619577865215983, atol=1e-5)
end

@testset "eelis" begin
    N=6
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    state=eigvecs(anyon_ham(model))[:,1]
    EE_lis=anyon_eelis(model, state)
    @test length(EE_lis)==length(collect(1:N-1))
    @test all(EE_lis .> 0)
end

@testset "ee_Fibonacci_ladder_lis" begin
    N=3
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    state=eigvecs(anyon_ham(model))[:,1]
    EE_lis=anyonladder_eelis(model, kron(state, state))
    @test length(EE_lis)==length(collect(1:N-1))
    @test all(EE_lis .> 0)
end

@testset "translation_matrix" begin
    N=8
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    Mat=translation_matrix(model)
    @test size(Mat)==(47,47)
    @test isapprox(Mat*Mat',I,atol=1e-5) # Check if the matrix is unitary
    @test isapprox(Mat^N,I,atol=1e-5)
end

@testset "inversion_matrix" begin
    N=8
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    Mat=inversion_matrix(model)
    @test size(Mat)==(47,47)
    @test isapprox(Mat*Mat',I,atol=1e-5) # Check if the matrix is unitary
    @test isapprox(Mat*Mat,I,atol=1e-5) 
end 

@testset "braidingsq_basismapN3" begin
    N=3
    T = BitStr{N, Int}

    state=T(bit"010")
    ϕ = (1+√5)/2
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"010"), T(bit"010"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"010"), T(bit"000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"010"), T(bit"010"), exp(-6im*π/5), 0.0)

    state =T(bit"001")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"001"), T(bit"001"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"001"), T(bit"001"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"001"), T(bit"000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    state =T(bit"100")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"100"), T(bit"000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"100"), T(bit"100"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"100"), T(bit"100"), exp(-6im*π/5), 0.0)

    
    state =T(bit"000")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"000"), T(bit"100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"000"), T(bit"010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"000"), T(bit"001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    
    state =T(bit"101") # Not in PBC basis
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"101"), T(bit"101"), exp(-2im*π/5), 0.0)
    model_obc = AnyonModel(FibonacciAnyon(), N, pbc=false)
    @test FibonacciChain._braidingsq_apply(model_obc, state, 2) == (T(bit"101"), T(bit"101"), exp(-2im*π/5), 0.0)
end

@testset "braidingsq_basismapN6" begin
    N=6
    T = BitStr{N, Int}
    ϕ = (1+√5)/2
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    ## state 1 ###
    state =T(bit"000000")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"000000"), T(bit"100000"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"000000"), T(bit"010000"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"000000"), T(bit"001000"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"000000"), T(bit"000100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"000000"), T(bit"000010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"000000"), T(bit"000001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 2 ###
    state =T(bit"000001")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"000001"), T(bit"000001"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"000001"), T(bit"010001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"000001"), T(bit"001001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"000001"),  T(bit"000101"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"000001"), T(bit"000001"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"000001"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 3 ###
    state=T(bit"000010")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"000010"), T(bit"100010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"000010"), T(bit"010010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"000010"), T(bit"001010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"000010"), T(bit"000010"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"000010"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"000010"), T(bit"000010"), exp(-6im*π/5), 0.0)

    ## state 4 ###
    state =T(bit"000100")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"000100"), T(bit"100100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"000100"), T(bit"010100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"000100"), T(bit"000100"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"000100"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"000100"), T(bit"000100"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"000100"), T(bit"000101"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 5 ###
    state =T(bit"000101")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"000101"), T(bit"000101"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"000101"), T(bit"010101"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"000101"), T(bit"000101"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"000101"), T(bit"000001"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"000101"), T(bit"000101"), exp(-2im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"000101"), T(bit"000100"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 6 ###
    state =T(bit"001000")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"001000"), T(bit"101000"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"001000"), T(bit"001000"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"001000"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"001000"), T(bit"001000"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"001000"), T(bit"001010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"001000"), T(bit"001001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 7 ###
    state =T(bit"001001")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"001001"), T(bit"001001"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"001001"), T(bit"001001"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"001001"), T(bit"000001"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"001001"), T(bit"001001"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"001001"), T(bit"001001"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"001001"), T(bit"001000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 8 ###
    state =T(bit"001010")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"001010"), T(bit"101010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"001010"), T(bit"001010"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"001010"), T(bit"000010"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"001010"), T(bit"001010"), exp(-2im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"001010"), T(bit"001000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"001010"), T(bit"001010"), exp(-6im*π/5), 0.0)

    ## state 9 ###
    state =T(bit"010000")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"010000"), T(bit"010000"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"010000"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"010000"), T(bit"010000"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"010000"), T(bit"010100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"010000"), T(bit"010010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"010000"), T(bit"010001"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 10 ###
    state =T(bit"010001")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"010001"), T(bit"010001"), exp(-2im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"010001"), T(bit"000001"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"010001"), T(bit"010001"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"010001"), T(bit"010101"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"010001"), T(bit"010001"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"010001"), T(bit"010000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 11 ###
    state =T(bit"010010")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"010010"), T(bit"010010"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"010010"), T(bit"000010"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"010010"), T(bit"010010"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"010010"), T(bit"010010"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"010010"), T(bit"010000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"010010"), T(bit"010010"), exp(-6im*π/5), 0.0)

    ## state 12 ###
    state =T(bit"010100")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"010100"), T(bit"010100"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"010100"), T(bit"000100"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"010100"), T(bit"010100"), exp(-2im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"010100"), T(bit"010000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"010100"), T(bit"010100"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"010100"), T(bit"010101"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    
    ## state 13 ###
    state =T(bit"010101")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"010101"), T(bit"010101"), exp(-2im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"010101"), T(bit"000101"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"010101"), T(bit"010101"), exp(-2im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"010101"), T(bit"010001"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"010101"), T(bit"010101"), exp(-2im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"010101"), T(bit"010100"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))

    ## state 14 ###
    state =T(bit"100000")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"100000"), T(bit"000000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"100000"), T(bit"100000"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"100000"), T(bit"101000"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"100000"), T(bit"100100"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"100000"), T(bit"100010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"100000"), T(bit"100000"), exp(-6im*π/5), 0.0)

    ## state 15 ###
    state =T(bit"100010")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"100010"), T(bit"000010"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"100010"), T(bit"100010"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"100010"), T(bit"101010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"100010"), T(bit"100010"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"100010"), T(bit"100000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"100010"), T(bit"100010"), exp(-2im*π/5), 0.0)

    ## state 16 ###
    state =T(bit"100100")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"100100"), T(bit"000100"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"100100"), T(bit"100100"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"100100"), T(bit"100100"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"100100"), T(bit"100000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"100100"), T(bit"100100"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"100100"), T(bit"100100"), exp(-6im*π/5), 0.0)

    ## state 17 ###
    state =T(bit"101000")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"101000"), T(bit"001000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"101000"), T(bit"101000"), exp(-2im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"101000"), T(bit"100000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"101000"), T(bit"101000"), exp(-6im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"101000"), T(bit"101010"), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"101000"), T(bit"101000"), exp(-6im*π/5), 0.0)

    ## state 18 ###
    state =T(bit"101010")
    @test FibonacciChain._braidingsq_apply(model, state, 1) == (T(bit"101010"), T(bit"001010"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 2) == (T(bit"101010"), T(bit"101010"), exp(-2im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 3) == (T(bit"101010"), T(bit"100010"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 4) == (T(bit"101010"), T(bit"101010"), exp(-2im*π/5), 0.0)
    @test FibonacciChain._braidingsq_apply(model, state, 5) == (T(bit"101010"), T(bit"101000"), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2))
    @test FibonacciChain._braidingsq_apply(model, state, 6) == (T(bit"101010"), T(bit"101010"), exp(-2im*π/5), 0.0)

end

@testset "braidingsq_matrix" begin
    N=3
    model_obc = AnyonModel(FibonacciAnyon(), N, pbc=false)
    ϕ = (1+√5)/2
    # Test braidingsq_matrix = the formulas of B^2 or not
    @test FibonacciChain.braidingsq_matrix(model_obc, 2) ≈  ComplexF64[
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
    @test FibonacciChain.braidingsq_matrix(model_obc, 2) ≈ U' * B^2 * U

    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    @test FibonacciChain.braidingsq_matrix(model, 2) ≈ (U' * B^2 * U)[1:4, 1:4]

    @test FibonacciChain.braidingsq_matrix(model, 1) ≈ ComplexF64[
    exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2) 0.0 0.0 (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2);
    0.0 exp(-6im*π/5) 0.0 0.0; 
    0.0 0.0  exp(-6im*π/5) 0.0; 
    (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2) 0.0 0.0 exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)]

    @test FibonacciChain.braidingsq_matrix(model, 3) ≈ ComplexF64[
    exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2) (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2) 0.0 0.0;
    (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2) exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1) 0.0 0.0;
    0.0 0.0 exp(-6im*π/5) 0.0; 
    0.0 0.0 0.0 exp(-6im*π/5)]

    ⊗(A,B) = kron(A,B)
    idx = [i.buf+1 for i in anyon_basis(model_obc)]
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
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    @test braidingsqmap(model, state, 2) ≈ [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+3(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+3(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 4exp(-6im*π/5)]
    @test braidingsqmap(model, state, 1) ≈ [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+4(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), 3exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+4(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1))]
    @test braidingsqmap(model, state, 3) ≈ [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+2(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+2(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 3exp(-6im*π/5), 4exp(-6im*π/5)]
end

@testset "spatial_correlation" begin
    N=12
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    mes = zeros(length(anyon_basis(model)))
    mes[233] = 1/√2 
    mes[end] = 1/√2 # set the last two qubits to be in the Bell state

    sclis = [spatial_correlation(model, mes, i, j) for i in 1:N for j in 1:N if j!=i]
    @test sclis ≈ log(2)*ones(12*11)

    N=6
    model_ising = AnyonModel(IsingAnyon(), N, pbc=true)
    mes = zeros(length(anyon_basis(model_ising)))
    mes[22] = 1/√2 
    mes[43] = 1/√2
    sclis = [spatial_correlation(model_ising, mes, i, j) for i in 1:N for j in 1:N if j!=i]
    @test sclis ≈ log(2)*ones(6*5)
end

@testset "temporal_correlation" begin
    # N=12
    # τ = 1000.0
    # mes = zeros(length(anyon_basis(N, true, anyon_type=:Fibo)))
    # mes[233] = 1/√2 
    # mes[end] = 1/√2 # set the last two qubits to be in the Bell state
    
    # sample = [0 1 1 1 1 1; 1 0 0 0 0 1; 1 1 1 1 1 1; 0 0 1 1 1 1; 1 0 1 1 1 1; 0 1 1 0 0 1; 1 1 1 1 1 1; 0 1 1 1 0 0; 0 1 1 1 1 0; 1 0 1 1 0 1]
    # D= 10 
    # sample = zeros(Int, D, length(2:2:N))
    # tc = temporal_correlation(τ, mes, sample, div(N,2), 5, 8)
    
    # tclis = [temporal_correlation(τ, mes, sample, div(N,2), i, j) for i in 1:D-1 for j in i+1:D]
    # @test tc ≈ 0.03098628964295691

    N=4
    τ = 1000.0
    model = AnyonModel(IsingAnyon(), N, pbc=true)
    mes = zeros(length(anyon_basis(model)))
    mes[1]=1.0 # set the last two qubits to be in the Bell state
    
    D= 10 

    sample = BitMatrix(zeros(Int, D, N))
    config = MeasureConfig(τ=τ, t₂ =div(D,2), rng = MersenneTwister(1234), mode=:sample)
    measure_outcome = generate_state_by_measurement(model, mes, sample, config)
    statelis, free_energy = measure_outcome.state, measure_outcome.free_energy
    # Noting that the first state of statelis is not mes.
    spatial_corr_lis = [spatial_correlation(model, st, 1, div(N,2)) for st in statelis]
    @test spatial_corr_lis ≈ log(2)*ones(div(D,2))

    ref_config = MeasureConfig(rng = MersenneTwister(1234))
    final_stlis, samples, sample_free_energy = reference_evolution(model, statelis, sample, div(N,2), 2, 4, ref_config)
    final_st = final_stlis[end]
    tc = temporal_correlation(model, final_st)
    # tclis = [temporal_correlation(τ, mes, sample, div(N,2), i, j, anyon_type=:IsingX) for i in 1:D-1 for j in i+1:D]
    @test isapprox(tc, 0.0, atol=1e-6)
end

@testset "mutual_information" begin
    N=6
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    mes = zeros(length(anyon_basis(model)))
    mes[1]=1.0 # set the last two qubits to be in the Bell state
    
    mi = mutual_information(model, ([2], [5]), mes)
    @test mi ≈ 0.0

    model_ising = AnyonModel(IsingAnyon(), N, pbc=true)
    mes = zeros(length(anyon_basis(model_ising)))
    mes[1]=1.0 # set the last two qubits to be in the Bell state
    
    mi = mutual_information(model_ising, ([2], [5]), mes)
    @test mi ≈ 0.0
end

@testset "tri_mutual_information" begin
    N=6
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    mes = zeros(length(anyon_basis(model)))
    mes[1]=1.0 # set the last two qubits to be in the Bell state
    
    tmi = tri_mutual_information(model, ([1], [3], [5]), mes)
    @test tmi ≈ 0.0

    model_ising = AnyonModel(IsingAnyon(), N, pbc=true)
    mes = zeros(length(anyon_basis(model_ising)))
    mes[1]=1.0 # set the last two qubits to be in the Bell state
    
    tmi = tri_mutual_information(model_ising, ([1], [3], [5]), mes)
    @test tmi ≈ 0.0
end