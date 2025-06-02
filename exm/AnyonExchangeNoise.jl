using FibonacciChain
using LinearAlgebra
using BitBasis
include("FitEntEntScal.jl")

N=10
energy, states = eigen(Fibonacci_Ham(N))
antiGS= states[:, 1]
len= length(antiGS)
vecGS = reshape(antiGS*antiGS', len^2)
splitlis=Vector(1:N-1)

EE_lis=zeros(length(splitlis))
for m in eachindex(EE_lis)
    subrho=ladderrdm(N, collect(1:splitlis[m]), vecGS)
    EE_lis[m]=ee(subrho)
end

cent, fig = fitCCEntEntScal(EE_lis; mincut=2,pbc=true)
savefig(fig, "./exm/double_Fibo_ee_scaling_10.pdf")
display(fig)
