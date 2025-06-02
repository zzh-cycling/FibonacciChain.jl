using FibonacciChain
using LinearAlgebra
using JLD
include("FitEntEntScal.jl")

function ee_Fibo_scaling_fig(N::Int64, state::Vector{ET},fit::String, mincut::Int64=1, pbc::Bool=true) where {ET}
    splitlis=Vector(1:N-1)
    EElis=eelis_Fibo_state(N, splitlis, state, pbc)

    if fit=="CC" 
        cent, fig=fitCCEntEntScal(EElis; mincut=mincut, pbc)
    end

    if fit=="Page"
        cent, fig=fitpage_curve(EElis; mincut=mincut)
    end

    if fit=="L+lnL"
        cent, fig=fitLpluslnL(EElis; mincut=mincut)
    end
    return cent, fig
end


N=20
energy, states = eigen(Fibonacci_Ham(20))
antiGS= states[:, 1]
save("./exm/Fibo_antiGS_20.jld", "antiGS", antiGS)
splitlis=Vector(1:N-1)
EElis=eelis_Fibo_state(N, splitlis, antiGS)
cent, fig = fitCCEntEntScal(EElis; mincut=4,pbc=true)
savefig(fig, "./exm/Fibo_ee_scaling_20.pdf")
display(fig)

energy, states = eigen(Fibonacci_ferroHam(20))
ferroGS= states[:, 1]
save("./exm/Fibo_ferroGS_20.jld", "ferroGS", ferroGS)
splitlis=Vector(1:N-1)
EElis=eelis_Fibo_state(N, splitlis, ferroGS)
cent, fig = fitCCEntEntScal(EElis; mincut=2,pbc=true)
savefig(fig, "./exm/ferroFibo_ee_scaling_20.pdf")
display(fig)
