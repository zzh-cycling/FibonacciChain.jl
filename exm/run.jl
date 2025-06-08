using FibonacciChain
using LinearAlgebra
using BitBasis
using JLD

N=12
energy, states = eigen(Fibonacci_Ham(N))
antiGS= states[:, 1]
vecGS = kron(antiGS, antiGS)
splitlis=Vector(1:N-1)

###==varied p==###
probabilitylis=collect(0.0:0.05:1.0)
centlis=similar(probabilitylis)

for (idx, i) in enumerate(probabilitylis)
    @show i
    state = ladderChoi(N, i, vecGS)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        @show m
        subrho=ladderrdm(N, collect(1:splitlis[m]), state)
        @time EE_lis[m]=ee(subrho)
    end
    save("./exm/data/double_Fibo_ee_scaling_$(N)_prob_$(i).jld", "state", state, "EE_lis", EE_lis)
end