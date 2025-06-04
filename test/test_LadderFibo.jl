using FibonacciChain
using Test
using LinearAlgebra 

@testset "ladderrdm" begin
    # Test the ladder_rdm function
    N=2
    T = BitStr{N, Int}
    state = collect(1:3)
    vec = reshape(state*state', 9)
    rdm = ladderrdm(T, Int[1], vec)
    @test size(rdm) == (4, 4)
    @test rdm == [100 20 20 4; 20 40 4 8; 20 4 40 8; 4 8 8 16]

    # Test the ladder_rdm with a different basis
    vec/=norm(vec)
    rdm2 = FibonacciChain.ladderrdm(T, Int[1], vec)
    @test ishermitian(rdm2)
    @test tr(rdm2) ≈ 1.0 atol=1e-6
    @test sum(eigvals(rdm2)) ≈ 1.0 atol=1e-6
end

@testset "ladderChoi" begin
    N=3
    T = BitStr{N, Int}
    state = collect(1:4)
    state = reshape(state*state', 4^2)
    @test Float64.(ladderChoi(N, 0.0, state)) ≈ state/norm(state)
    
    ϕ = (1+√5)/2
    onechain_state = [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+3(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+3(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 4exp(-6im*π/5)]
    st = reshape(onechain_state*transpose(onechain_state), 16)
    @test ladderChoi(N, 1.0, state) ≈ st/norm(st)
 
    # At least for N = 4, the ladder Choi state is invariant under the ladder translation map twice.
    N=4
    energy, states = eigen(Fibonacci_Ham(N))
    antiGS= states[:, 1]
    len= length(antiGS)
    vecGS = reshape(antiGS*antiGS', len^2)
    plis = collect(0.0:0.1:1.0)
    for p in plis
        state = ladderChoi(N, p, vecGS)
        @show p
        @test isapprox(reduce((x, _) -> laddertranslationmap(N, x), 1:2; init=state),state, atol=1e-10)
        # @show isapprox(laddertranslationmap(N, laddertranslationmap(N, state)), state, atol=1e-10)
    end
    
    N=6
    energy, states = eigen(Fibonacci_Ham(N))
    antiGS= states[:, 1]
    len= length(antiGS)
    state = reshape(antiGS*antiGS', len^2)

    p=0.5
    for i in 2:2:N
        state=(1-p)*state+p*ladderbraidingmap(N, state, i)
        state/=norm(state) # normalize the state after each braiding
    end
end


@testset "ladderbraidingmap" begin
    N=3
    state = collect(1:4)
    ϕ = (1+√5)/2
    onechain_state = [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+3(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+3(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 4exp(-6im*π/5)]
    @test ladderbraidingmap(N, reshape(state*state',16), 2) ≈ reshape(onechain_state*transpose(onechain_state), 16)

    # translation invariance test
    N = 4
    state = fill(1,7)
    order = [1, 3, 4, 6, 7, 2, 5]
    onechain_state = braidingmap(N, braidingmap(N, state, 2),4)
    onechain_state[order][order] == onechain_state
    
    twochain_state = ladderbraidingmap(N, ladderbraidingmap(N, reshape(state*state',49), 2),4)
    @test laddertranslationmap(N, laddertranslationmap(N, twochain_state)) ≈ twochain_state
end

@testset "laddertranslationmap" begin
    N=3
    state = collect(1:4)
    order = [1, 3, 4, 2]
    ϕ = (1+√5)/2
    @test Int64.(laddertranslationmap(N, reshape(state*state',16))) ≈ reshape(state[order]*transpose(state[order]), 16)
end