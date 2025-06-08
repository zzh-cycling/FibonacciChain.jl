using JLD
include("FitEntEntScal.jl")
using Plots
using LaTeXStrings

probabilitylis=collect(0.0:0.05:1.0)
centlis=similar(probabilitylis)
centlisL=similar(probabilitylis)
centlisR=similar(probabilitylis)

for i in eachindex(probabilitylis)
    @show i
    state, EElis = load("./exm/data/double_Fibo_ee_scaling_10_prob_$(probabilitylis[i]).jld", "state", "EE_lis")
    cent, fig = fitCCEntEntScal(EElis; mincut=2, pbc=true)
    centlis[i] = cent

    centL, figL = fitpart(EElis; mincut=2, pbc=true, part=:L)
    centR, figR = fitpart(EElis; mincut=2, pbc=true, part=:R)
    centlisL[i] = centL
    centlisR[i] = centR

    # savefig(fig, "./exm/fig/double_Fibo_ee_scaling_10_prob_$(probabilitylis[i]).pdf")
    # savefig(figL, "./exm/fig/double_Fibo_ee_scaling_10_prob_$(probabilitylis[i])_L.pdf")
    # savefig(figR, "./exm/fig/double_Fibo_ee_scaling_10_prob_$(probabilitylis[i])_R.pdf")
end

fig = plot(probabilitylis, centlis, xlabel=L"p", ylabel=L"c_{cent}", label=L"c_{total}", marker=:circle)
plot!(fig, probabilitylis, centlisL, label=L"c_{L}", marker=:circle)
plot!(fig, probabilitylis, centlisR, label=L"c_{R}", marker=:circle)
savefig(fig, "./exm/fig/double_Fibo_10_pvscent.pdf")