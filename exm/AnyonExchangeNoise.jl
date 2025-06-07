using FibonacciChain
using LinearAlgebra
using BitBasis
using JLD
include("FitEntEntScal.jl")

N=8
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
savefig(fig, "./exm/double_Fibo_ee_scaling_$(N).pdf")
display(fig)

###==varied p==###
probabilitylis=collect(0.0:0.05:1.0)
centlis=similar(probabilitylis)
for (idx, i) in enumerate(probabilitylis)
    state = ladderChoi(N, i, vecGS)
    @show laddertranslationmap(N, laddertranslationmap(N, state)) ≈ state
    EE_lis=zeros(length(splitlis))
    # save("./exm/data/double_Fibo_ee_scaling_10_prob_$(i).jld", "state", state, "EE_lis", EE_lis)
    for m in eachindex(EE_lis)
        subrho=ladderrdm(N, collect(1:splitlis[m]), state)
        EE_lis[m]=ee(subrho)
    end
    cent, fig = fitCCEntEntScal(EE_lis; mincut=2,pbc=true)
    centlis[idx]=cent
    # savefig(fig, "./exm/double_Fibo_ee_scaling_10_prob_$(i).pdf")
    display(fig)
end
using Plots
plot(probabilitylis, centlis, xlabel=L"p", ylabel=L"c_{cent}", label=false, marker=:circle)
savefig("./exm/fig/double_Fibo_10_pvscent.pdf")