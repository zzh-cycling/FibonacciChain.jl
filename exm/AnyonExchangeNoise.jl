using FibonacciChain
using LinearAlgebra
using BitBasis
using JLD
using Plots
include("FitEntEntScal.jl")

# one pure Fibonacci chain #

N=18
energy, states = eigen(Fibonacci_Ham(N))
antiGS= states[:, 1]
vecGS = kron(antiGS, antiGS)
splitlis=Vector(1:N-1)

EE_lis=zeros(length(splitlis))
for m in eachindex(EE_lis)
    subrho=rdm_Fibo(N, collect(1:splitlis[m]), antiGS)
    @time EE_lis[m]=ee(subrho)
end



cent, fig = fitCCEntEntScal(EE_lis; mincut=2,pbc=true)
# savefig(fig, "./exm/fig/antiferro_Fibo_ee_scaling_$(N).pdf")
display(fig)

### two noisy Fibonacci chain with varied p

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

## one Noisy chain##

N=12
energy, states = eigen(Fibonacci_Ham(N))
antiGS= states[:, 1]
splitlis=Vector(1:N-1)

for i in 2:2:N
    antiGS= braidingsqmap(N, antiGS, i)
    antiGS/= norm(antiGS)
end

EE_lis=zeros(length(splitlis))
for m in eachindex(EE_lis)
    subrho=rdm_Fibo(N, collect(1:splitlis[m]), antiGS)
    EE_lis[m]=ee(subrho)
end
# save("./exm/data/single_Fibo_ee_scaling_$(N).jld", "state", antiGS, "EE_lis", EE_lis)
cent, fig = fitCCEntEntScal(EE_lis; mincut=2,pbc=true)
savefig(fig, "./exm/fig/single_Fibo_ee_scaling_$(N).pdf")
display(fig)

