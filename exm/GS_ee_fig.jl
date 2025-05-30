using FibonacciChain
using Plots
using LaTeXStrings
using LsqFit
include("FitEntEntScal.jl")

function ee_Fibo_scaling_fig(N::Int64, state::Vector{ET},fit::String, pbc::Bool=true) where {ET}
    splitlis=Vector(1:N-1)
    EElis=eelis_Fibo_state(N, splitlis, state, pbc)

    if fit=="CC" 
        cent, fig=fitCCEntEntScal(EElis; mincut=1, pbc)
    end

    if fit=="Page"
        cent, fig=fitpage_curve(EElis; mincut=1)
    end

    if fit=="L+lnL"
        cent, fig=fitLpluslnL(EElis; mincut=1)
    end
    return cent, fig
end

energy, states = eigen(Fibonacci_Ham(18))
GS= states[:, end]
cent, fig = ee_Fibo_scaling_fig(18, GS, "CC")

display(fig)

energy, states = eigen(Fibonacci_ferroHam(16))
GS= states[:, end]
cent, fig = ee_Fibo_scaling_fig(16, GS, "CC")

display(fig)