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

    # st1 = collect(1:3); st2 = collect(4:6)
    # vec = reshape(st1*st2', 9)
    # rdm = ladderrdm(T, Int[1], vec)
    # @test size(rdm) == (4, 4)
    # @test rdm == [100 20 20 4; 20 40 4 8; 20 4 40 8; 4 8 8 16]
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
        @test isapprox(reduce((x, _) -> laddertranslationmap(N, x), 1:2; init=state),state, atol=1e-10)
    end
    
    N=6
    energy, states = eigen(Fibonacci_Ham(N))
    antiGS= states[:, 1]
    len= length(antiGS)
    vecGS = reshape(antiGS*antiGS', len^2)
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
    @test ladderbraidingmap(N, reshape(state*state',16), 2) ≈ reshape(onechain_state*transpose(onechain_state), 16)

    st1= collect(1:4);st2= collect(5:8)
    st1map = [exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2)+3(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 2exp(-6im*π/5), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+3(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 4exp(-6im*π/5)]
    st2map = [5*(exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2))+7(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2), 6exp(-6im*π/5), 5(exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)+7(exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1)), 8exp(-6im*π/5)]
    @test braidingmap(N, st1, 2) ≈ st1map
    @test braidingmap(N, st2, 2) ≈  st2map
    @test ladderbraidingmap(N, reshape(st1*st2',16), 2) ≈ reshape(st1map*transpose(st2map), 16)

    # translation invariance test
    N = 4
    state = fill(1,7)
    order = [1, 3, 4, 6, 7, 2, 5]
    onechain_state = braidingmap(N, braidingmap(N, state, 2),4)
    @test onechain_state[order][order] ≈ onechain_state
    
    twochain_state = ladderbraidingmap(N, ladderbraidingmap(N, reshape(state*state',49), 2),4)
    @test laddertranslationmap(N, laddertranslationmap(N, twochain_state)) ≈ twochain_state

end

@testset "laddertranslationmap" begin
    N=3
    state = collect(1:4)
    order = [1, 3, 4, 2]
    ϕ = (1+√5)/2
    @test Int64.(laddertranslationmap(N, reshape(state*state',16))) ≈ reshape(state[order]*transpose(state[order]), 16)

    st1= collect(1:4);st2= collect(5:8)
    @test Int64.(laddertranslationmap(N, reshape(st1*st2',16))) ≈ reshape(st1[order]*transpose(st2[order]), 16)
end