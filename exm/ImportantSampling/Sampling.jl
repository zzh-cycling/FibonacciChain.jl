using FibonacciChain
using LinearAlgebra
using JLD
include("../FitEntEntScal.jl")


function Born_eelis(N::Int64, τ::Float64, initial_state::Vector{ET}, measurement_sites::Vector{Int}, pbc::Bool=true) where {ET}

    final_states, trajectories, probabilities = measurement_enumeration(
        N, τ, initial_state, measurement_sites, pbc)
    
    num_final_states = length(final_states)
    println("Total number of final states: $(num_final_states)")
    println("Expected number: $(2^length(measurement_sites))")
    
    total_prob = sum(probabilities)
    println("Total probability: $(total_prob)")
    
    
    entropies = Vector{Vector{Float64}}(undef, num_final_states)
    for i in 1:num_final_states
        entropies[i] = eelis_Fibo_state(N, splitlis, final_states[i], pbc)
    end
    
    # avg_eelis = sum(probabilities .* entropies)
    # println("Born-averaged entanglement entropy: $(avg_entropy)")
    

    println("\nSome representative trajectories:")
    sorted_indices = sortperm(probabilities, rev=true)
    for i in 1:min(2, num_final_states)
        idx = sorted_indices[i]
        println("Trajectory $(trajectories[idx]): probability = $(probabilities[idx]/total_prob)")
        println("  Entanglement entropy = $(entropies[idx])")
    end
    
    return final_states, trajectories, probabilities, entropies
end

N=18
energy, states = eigen(Fibonacci_Ham(N))
antiGS= states[:, 1]
splitlis=Vector(1:N-1)

ee_lis=zeros(length(splitlis))
for m in eachindex(ee_lis)
    subrho=rdm_Fibo(N, collect(1:splitlis[m]), antiGS)
    @time ee_lis[m]=ee(subrho)
end
cent, fig = fitCCEntEntScal(ee_lis; mincut=2, pbc=true)
display(fig)

γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999])
τlis = atanh.(γlis)
centlis = zeros(length(τlis))


for (idx, τ) in enumerate(τlis)
    println("τ = $τ")
    measurement_sites = collect(2:2:N)
    @time final_states, trajectories, probabilities, entropies = Born_eelis(N, τ, antiGS, measurement_sites)
    
    save("./exm/data/Born_eelis_N$(N)_τ$(τ).jld", "final_states", final_states, 
         "trajectories", trajectories, "probabilities", probabilities, "entropies", entropies)

    EE_lis= sum(probabilities.*entropies)
    cent, fig = fitCCEntEntScal(EE_lis; mincut=2, pbc=true)
    centlis[idx] = cent
    display(fig)
    # savefig(fig, "./exm/fig/Born_eeliscc_N$(N)_τ$(τ).pdf")
end

save("./exm/data/Born_eeliscc_N$(N).jld", "centlis", centlis, "τlis", τlis)
fig=plot(γlis, centlis, xlabel=L"\gamma=\tanh(\tau)", ylabel=L"c_{cent}", marker=:o, legend=false)
savefig(fig, "./exm/fig/Born_eeliscc_N$(N).pdf")

centlis12=load("./exm/data/Born_eeliscc_N12.jld", "centlis")
centlis14=load("./exm/data/Born_eeliscc_N14.jld", "centlis")
centlis16=load("./exm/data/Born_eeliscc_N16.jld", "centlis")

fig = plot(γlis, centlis12, label=L"L=12", marker=:o, xlabel=L"\gamma=\tanh(\tau)", ylabel=L"c_{cent}",   legend_background_color=nothing,
        legend_foreground_color=nothing,)
plot!(fig, γlis, centlis14, label=L"L=14", marker=:x)
plot!(fig, γlis, centlis16, label=L"L=16", marker=:diamond)