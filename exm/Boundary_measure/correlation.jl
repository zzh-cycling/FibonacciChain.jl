using FibonacciChain
using LinearAlgebra
using Statistics
using JLD

# Compute the boundary state quantum and classical correlations
function compute_correlations(L::Int; τ::Float64=0.0)
 
    data = load("exm/data/Boundary_measure/Born_Samples_N22/Born_Samples_N$(L)_τ$(τ).jld")
    sample_measured_states = data["sample_measured_states"]
    samples = data["samples"]
    
    model = AnyonModel(FibonacciAnyon(), L; pbc=true)
    n_samples = size(samples, 1)
    
 
    spins = 2 .* samples .- 1
    spin_avg = mean(spins, dims=1)[:]
    
 
    dlis = collect(1:div(L, 2))
    quantum_corr = zeros(length(dlis))
    quantum_corr_std = zeros(length(dlis))
    measure_corr = zeros(length(dlis))
    measure_corr_std = zeros(length(dlis))
    
    site1 = 1
    
    for (i, d) in enumerate(dlis)
        site2 = site1 + d
        
        # quantum correlation: I(R₁, R₂) =  S(R₁) + S(R₂) - S(R₁ ∪ R₂)
        q_corr_samples = [spatial_correlation(model, sample_measured_states[j], site1, site2) for j in 1:n_samples]
        quantum_corr[i] = mean(q_corr_samples)
        quantum_corr_std[i] = std(q_corr_samples) / sqrt(n_samples)
        
        # classical correlation: C(R₁, R₂) = ⟨σ_z(R₁) σ_z(R₂)⟩ - ⟨σ_z(R₁)⟩⟨σ_z(R₂)⟩
        m_corr_samples = spins[:, site1] .* spins[:, site2] .- spin_avg[site1] * spin_avg[site2]
        measure_corr[i] = mean(m_corr_samples)
        measure_corr_std[i] = std(m_corr_samples) / sqrt(n_samples)
    end
    
    return dlis, quantum_corr, quantum_corr_std, measure_corr, measure_corr_std
end


if abspath(PROGRAM_FILE) == @__FILE__
    L = 22
    dlis, q_corr, q_std, m_corr, m_std = compute_correlations(L)
    
    println("Distance | Quantum Corr ± std | Measure Corr ± std")
    println("-"^60)
    for (d, q, qs, m, ms) in zip(dlis, q_corr, q_std, m_corr, m_std)
        println("   $d     |   $(round(q, digits=4)) ± $(round(qs, digits=4))   |   $(round(m, digits=4)) ± $(round(ms, digits=4))")
    end
end

