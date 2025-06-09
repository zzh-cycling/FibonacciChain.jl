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

⊗(A,B) = kron(A,B)
ϕ = (1+√5)/2
Z=[1 0;0 -1]
X=[0 1;1 0]
P0=[1 0;0 0]
P1=[0 0;0 1]
σ2t = (1-2ϕ^(-1)) * Z +2ϕ^(-3/2) * X

R =(I(2)⊗(exp(1im*π/10)*exp(-7im*π/10 .* Z))⊗I(2))[idx, idx]
idx = [i.buf+1 for i in Fibonacci_basis(3,false)]
B = exp(1im*π/10)*(P0 ⊗ exp(-7im*π/10 .* Z) ⊗ P1 + P1 ⊗ exp(-7im*π/10 .* Z) ⊗ P0 + P1 ⊗ exp(+7im*π/10 .* Z) ⊗ P1 + P0 ⊗ exp(-7im*π/10 .* σ2t) ⊗ P0)

@show Bc=B[idx, idx]


Bc^2

## one chain##

N=10
energy, states = eigen(Fibonacci_Ham(N))
antiGS= states[:, 1]
splitlis=Vector(1:N-1)

for i in 2:2:N
    antiGS= braidingmap(N, antiGS, i)
    antiGS/= norm(antiGS)
end

EE_lis=zeros(length(splitlis))
for m in eachindex(EE_lis)
    subrho=rdm_Fibo(N, collect(1:splitlis[m]), antiGS)
    EE_lis[m]=ee(subrho)
end
cent, fig = fitCCEntEntScal(EE_lis; mincut=2,pbc=true)
display(fig)