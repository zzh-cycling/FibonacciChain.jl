using FibonacciChain
using LinearAlgebra
using JLD
using Arpack
# include("FitEntEntScal.jl")

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

N=34
energy, states =  eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
antiGS= states[:, 1]
EElis=eelis_Fibo_state(N, antiGS)
save("./exm/Fibo_antiGS_N$(N)_EElis.jld", "antiGS", antiGS, "EElis", EElis)
# cent, fig = fitCCEntEntScal(EElis; mincut=4,pbc=true)
# savefig(fig, "./exm/Fibo_ee_scaling_N$(N).pdf")
# display(fig)

# energy, states =  eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR) # Noting here sparse is for anti, not ferro
# ferroGS= states[:, 1]
# EElis=eelis_Fibo_state(N, ferroGS)
# save("./exm/Fibo_ferroGS_N$(N)_EElis.jld", "ferroGS", ferroGS, "EElis", EElis)
# cent, fig = fitCCEntEntScal(EElis; mincut=2,pbc=true)
# savefig(fig, "./exm/ferroFibo_ee_scaling_N$(N).pdf")
# display(fig)
