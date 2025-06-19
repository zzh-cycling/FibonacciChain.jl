using FibonacciChain
using LinearAlgebra
using JLD
using Arpack
include("../FitEntEntScal.jl")

L_list = collect(10:2:20)
γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999])
τlis = atanh.(γlis)

ShannonEnt_Llis = Vector{Vector{Float64}}(undef, length(L_list))

for (id, N) in enumerate(L_list)
    @show N
    energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
    antiGS= states[:, 1]

    SEntlis = zeros(length(τlis))

    for (idx, τ) in enumerate(τlis)
        println("τ = $τ")
        measurement_sites = collect(2:2:N)
        final_states, trajectories, probabilities = measurement_enumeration(N, τ, antiGS, measurement_sites)

        # save("./exm/data/Born_eelis_N$(N)_τ$(τ).jld", "trajectories", trajectories, "probabilities", probabilities, "entropies", entropies)

        ShannonEnt = sum(x->-x*log(x), probabilities)
        SEntlis[idx] = ShannonEnt
    end

    ShannonEnt_Llis[id] = SEntlis
end

save("./exm/data/Born_enumedmeasure_ShannonEnt_N1020.jld", "ShannonEnt_Llis", ShannonEnt_Llis)