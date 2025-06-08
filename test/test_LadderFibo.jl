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
    @test rdm == [25 15 15 9; 15 45 9 27; 15 9 45 27; 9 27 27 81]

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
    st2 = zeros(18); st2[end-1] = 1.0; st2[end] = 1.0; st2/=norm(st2) # |101000>+|101010> state
    splitlis = collect(1:3) 
    rdm1 = rdm_Fibo(N, splitlis, st1)
    # |000><000|
    rdm2 = rdm_Fibo(N, splitlis, st2)
    # |101><101|
    rdm = ladderrdm(N, splitlis, kron(st1, st2))
    # |000101><000101|
    @test kron(rdm1,rdm2) == rdm

    st1 = zeros(18); st1[3] = 1.0; st1[5] = 1.0; st1/=norm(st1) # |000101>+|000010> state
    st2 = zeros(18); st2[end-1] = 1.0; st2[end] = 1.0; st2/=norm(st2) # |101000>+|101010> state
    splitlis = collect(4:6) 
    rdm1 = rdm_Fibo(N, splitlis, st1)
    # 1/2(|101><101|+|010><010|)
    rdm2 = rdm_Fibo(N, splitlis, st2)
    # 1/2(|000><000|+|010><010|)
    rdm = ladderrdm(N, splitlis, kron(st1, st2))
    # 1/4（|101000><101000|+|101010><101010|+|010000><010000|+|010010><010010|） 
    # 11, 13, 21, 23 is 010000, 010010, 101000, 101010
    @test kron(rdm1,rdm2) == rdm

    # test whether the ladderChoi is product state of two mapped anti-GS states
    N=6
    energy, states = eigen(Fibonacci_Ham(N))
    antiGS= states[:, 1]
    len= length(antiGS)
    vecGS = kron(antiGS, antiGS)
    for i in 2:2:N
        antiGS = braidingmap(N, antiGS, i)
        antiGS/= norm(antiGS)
    end
    Choistate = ladderChoi(N, 1.0, vecGS)
    @test Choistate ≈ kron(antiGS, antiGS) 
    

    splitlis = collect(1:N-1)
    for m in eachindex(splitlis)
        subrho=ladderrdm(N, collect(1:splitlis[m]), Choistate)
        rdm_antiGS = rdm_Fibo(N, collect(1:splitlis[m]), antiGS)
        @test subrho ≈ kron(rdm_antiGS, rdm_antiGS)
    end
    
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