using FibonacciChain
using LinearAlgebra
using BitBasis
using JLD
# using Plots
using Arpack
using SparseArrays
include("FitEntEntScal.jl")

# one pure Fibonacci chain #

N=18
energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
antiGS= states[:, 1]
EE_lis=eelis_Fibo_state(N, state)

cent, fig = fitCCEntEntScal(EE_lis; mincut=2,pbc=true)
# savefig(fig, "./exm/fig/antiferro_Fibo_ee_scaling_$(N).pdf")
display(fig)

## two noisy Fibonacci chain with varied p

N=16
energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
antiGS= states[:, 1]
vecGS = kron(antiGS, antiGS)

###==varied p==###
probabilitylis=collect(0.0:0.05:1.0)
centlis=similar(probabilitylis)

for (idx, i) in enumerate(probabilitylis)
    @show i
    state = ladderChoi(N, i, vecGS)
    
    @time EE_lis=eelis_Fiboladder_state(N, state)
    save("./exm/data/double_Fibo_ee_scaling_$(N)_prob_$(i).jld", "state", state, "EE_lis", EE_lis)
end

# one Noisy chain##

N=32
energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
antiGS= states[:, 1]
splitlis=Vector(1:N-1)

for i in 2:2:N
    antiGS= braidingsqmap(N, antiGS, i)
    antiGS/= norm(antiGS)
end

EE_lis=eelis_Fibo_state(N, antiGS)
# save("./exm/data/single_Fibo_ee_scaling_$(N).jld", "state", antiGS, "EE_lis", EE_lis)
cent, fig = fitCCEntEntScal(EE_lis; mincut=2,pbc=true)
save("./exm/data/single_Fibo_ee_scaling_$(N).jld", "EE_lis", EE_lis)
savefig(fig, "./exm/fig/single_Fibo_ee_scaling_$(N).pdf")
display(fig)

