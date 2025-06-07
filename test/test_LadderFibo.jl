using FibonacciChain
using Test
using LinearAlgebra 

@testset "ladderrdm" begin
    # Test the ladder_rdm function
    N=2
    T = BitStr{N, Int}
    state = collect(1:3)
    vec1 = kron(state, state)
    rdm = ladderrdm(T, Int[1], vec1)
    @test size(rdm) == (4, 4)
    @test rdm == [100 20 20 4; 20 40 4 8; 20 4 40 8; 4 8 8 16]

    # Test the ladder_rdm with a different basis
    vec1/=norm(vec1)
    rdm2 = FibonacciChain.ladderrdm(T, Int[1], vec1)
    @test ishermitian(rdm2)
    @test tr(rdm2) ≈ 1.0 atol=1e-6
    @test sum(eigvals(rdm2)) ≈ 1.0 atol=1e-6

    N = 4
    state = collect(1:49)
    rdm = ladderrdm(N, Int[1, 2], state)

    N = 6
    st1 = collect(1:18); st2 = collect(19:36)
    vec1 = kron(st1, st2)
    rdm = ladderrdm(6, Int[1,2,3], vec1)
    len = length(Fibonacci_basis(div(6, 2), false))^2
    @test size(rdm) == (len, len)
    @test rdm == [100 20 20 4; 20 40 4 8; 20 4 40 8; 4 8 8 16]

    st1 = zeros(18); st1[end] = 1.0 # |101010> state
    st2 = zeros(18); st2[13] = 1.0 # |010101> state
    splitlis = collect(1:3)
    rdm1 = rdm_Fibo(N, splitlis, st1) # |101><101|
    rdm2 = rdm_Fibo(N, splitlis, st2)
    # |010><010|
    rdm = ladderrdm(N, splitlis, kron(st1, st2)) # |101010><101010|
    @test kron(rdm1,rdm2) == rdm
     
    # rdm basis 15 is 010101, 23 is 101010
    st1 = zeros(18); st1[3] = 1.0; st1[5] = 1.0; st1/=norm(st1) # |000101>+|000010> state
    st2 = zeros(18); st2[end-1] = 1.0; st2[end] = 1.0; st2/=norm(st2)
    splitlis = collect(1:3) # |101000>+|101010> state
    rdm1 = rdm_Fibo(N, splitlis, st1)
    # |000><000|
    rdm2 = rdm_Fibo(N, splitlis, st2)
    # |101><101|
    rdm = ladderrdm(N, splitlis, kron(st1, st2))
    # |000101><000101|
    @test kron(rdm1,rdm2) == rdm
end 

@testset "ladderChoi" begin
    N=3
    T = BitStr{N, Int}
    state = collect(1:4)
    state = kron(state, state)
    @test Float64.(ladderChoi(N, 0.0, state)) ≈ state/norm(state)
    
    ϕ = (1+√5)/2
    onechain_state = [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+3(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+3(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 4exp(-6im*π/5)]
    st = kron(onechain_state, onechain_state)
    @test ladderChoi(N, 1.0, state) ≈ st/norm(st)
 
    # At least for N = 4, the ladder Choi state is invariant under the ladder translation map twice.
    N=4
    energy, states = eigen(Fibonacci_Ham(N))
    antiGS= states[:, 1]
    len= length(antiGS)
    vecGS = kron(antiGS, antiGS)
    plis = collect(0.0:0.1:1.0)
    for p in plis
        state = ladderChoi(N, p, vecGS)
        @test isapprox(reduce((x, _) -> laddertranslationmap(N, x), 1:2; init=state),state, atol=1e-10)
    end
    
    N=6
    energy, states = eigen(Fibonacci_Ham(N))
    antiGS= states[:, 1]
    len= length(antiGS)
    vecGS = kron(antiGS, antiGS)
    plis = collect(0.0:0.1:1.0)
    for p in plis
        state = ladderChoi(N, p, vecGS)
        @test isapprox(reduce((x, _) -> laddertranslationmap(N, x), 1:2; init=state),state, atol=1e-10)
    end
end


@testset "ladderbraidingmap" begin
    N=3
    state = collect(1:4)
    ϕ = (1+√5)/2
    onechain_state = [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+3(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+3(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 4exp(-6im*π/5)]
    @test ladderbraidingmap(N, kron(state, state), 2) ≈ kron(onechain_state, onechain_state)

    st1= collect(1:4);st2= collect(5:8)
    st1map = [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+3(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+3(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 4exp(-6im*π/5)]
    st2map = [5*(exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2))+7(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 6exp(-6im*π/5), 5(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+7(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 8exp(-6im*π/5)]
    @test braidingmap(N, st1, 2) ≈ st1map
    @test braidingmap(N, st2, 2) ≈  st2map
    @test ladderbraidingmap(N, kron(st1, st2), 2) ≈ kron(st1map, st2map)

    # translation invariance test
    N = 4
    state = fill(1,7)
    order = [1, 3, 4, 6, 7, 2, 5]
    onechain_state = braidingmap(N, braidingmap(N, state, 2),4)
    @test onechain_state[order][order] ≈ onechain_state
    
    twochain_state = ladderbraidingmap(N, ladderbraidingmap(N, kron(state, state), 2),4)
    @test laddertranslationmap(N, laddertranslationmap(N, twochain_state)) ≈ twochain_state

end

@testset "laddertranslationmap" begin
    N=3
    state = collect(1:4)
    order = [1, 3, 4, 2]
    ϕ = (1+√5)/2
    @test Int64.(laddertranslationmap(N, kron(state, state))) ≈ kron(state[order], state[order])

    st1= collect(1:4);st2= collect(5:8)
    @test Int64.(laddertranslationmap(N, kron(st1, st2))) ≈ kron(st1[order], st2[order])
end