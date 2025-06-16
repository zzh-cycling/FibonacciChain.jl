using FibonacciChain
using LinearAlgebra
using JLD
using Arpack
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
        entropies[i] = eelis_Fibo_state(N, final_states[i], pbc)
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

N=20
energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
antiGS= states[:, 1]

ee_lis=eelis_Fibo_state(N, antiGS)
cent, fig = fitCCEntEntScal(ee_lis; mincut=4, pbc=true)
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


L_list = [10, 12, 14, 16, 18, 20]
centlis10=load("./exm/data/Born_eeliscc_N10.jld", "centlis")
centlis12=load("./exm/data/Born_eeliscc_N12.jld", "centlis")
centlis14=load("./exm/data/Born_eeliscc_N14.jld", "centlis")
centlis16=load("./exm/data/Born_eeliscc_N16.jld", "centlis")
centlis18=load("./exm/data/Born_eeliscc_N18.jld", "centlis")
centlis20=load("./exm/data/Born_eeliscc_N20.jld", "centlis")
totalcentlis = hcat(centlis10, centlis12, centlis14, centlis16, centlis18, centlis20)

fig = plot(γlis, totalcentlis, marker=:o, xlabel=L"\gamma=\tanh(\tau)", ylabel=L"c_{cent}",   legend_background_color=nothing,
legend_foreground_color=nothing, 
label=L_list',c= cgrad(:blues, length(L_list)), marker_z=L_list',line_z=L_list',colorbar=false, lw=2)
annotate!(fig, [(0.115, 0.42, text(L"L=", 10, :black))])
savefig(fig, "./exm/fig/Born_eeliscc_scaling.pdf")


N=20
trajectorieslis=Vector{Vector{Vector{Symbol}}}(undef, length(τlis))
probabilitieslis=Vector{Vector{Float64}}(undef, length(τlis))
entropieslis=Vector{Vector{Vector{Float64}}}(undef, length(τlis))
for (idx,i) in enumerate(τlis)
    data = load("./exm/data/Born_eelis_N$(N)_τ$(i).jld")
    trajectories = data["trajectories"]
    probabilities = data["probabilities"]
    entropies = data["entropies"]
    
    trajectorieslis[idx] = trajectories
    probabilitieslis[idx] = probabilities
    entropieslis[idx] = entropies
end

save("./exm/data/Born_enumed_eelis_N$(N).jld", "trajectorieslis", trajectorieslis, 
     "probabilitieslis", probabilitieslis, "entropieslis", entropieslis)

L_cent_tau_lis = Vector{Vector{Float64}}(undef, length(L_list))
for (idx,i) in enumerate(L_list)
    data = load("./exm/data/Born_eeliscc_N$(i).jld", "centlis")
    L_cent_tau_lis[idx] = data
end

save("./exm/data/Born_enumedcc_tau_N1020.jld", "L_cent_tau_lis", L_cent_tau_lis, "τlis", τlis)