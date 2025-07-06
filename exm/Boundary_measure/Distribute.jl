using JLD
using Plots
using Statistics
using LaTeXStrings
using Arpack
using FibonacciChain

function plot_inverse(τ::Float64)
    fig = plot(xlabel="1/ samples sizes", ylabel=L"\bar{S_c}", label=false, legend=:topright, legend_background_color=nothing, legend_foreground_color=nothing, )

    L_list = [22, 24, 26, 28, 30]
    c = cgrad(:blues, length(L_list), categorical=true)
    for (i, N) in enumerate(L_list)
        sample_measured_states, samples, sample_weights = load("./exm/data/Born_Samples_N$(N)/Born_Samples_N$(N)_τ$(τ).jld", "sample_measured_states", "samples", "sample_weights")
        FElis = -log.(sample_weights)
        sample_sizes = 1:length(FElis) 
        averages = [mean(FElis[1:i]) for i in sample_sizes]
        fig = plot!(fig, xscale=:log10, 1 ./sample_sizes, averages, label="$(N)", color = c[i])
    end
    savefig(fig, "./exm/fig/Sample_distribution/Inverse_average_Distribution_τ$(τ).pdf")

# savefig(fig, "./exm/fig/Fidelity_Entropy_Distribution_N22_τ3.8002.pdf")
end

function plot_compare(τ::Float64)
    fig = plot(xlabel="1/ samples sizes", ylabel=L"\bar{S_c}", label=false, legend=:topright, legend_background_color=nothing, legend_foreground_color=nothing, )

    L_list = collect(10:2:16)
    c = cgrad(:blues, length(L_list), categorical=true)
    
    γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999], 1.0)
    τlis = vcat(atanh.(vcat(collect(0.0:0.05:0.95), [0.99, 0.999])), 1e3)
    inds = findall(x -> x == τ, τlis)
    for (i, N) in enumerate(L_list)
        energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
        antiGS = states[:,1]
        measurement_sites = collect(2:2:N)
        @show (τ, i, N)
        sample_measured_states, samples, sample_weights = Sampling(N, τ, antiGS, measurement_sites, 10000)
        save("exm/data/Samples_Compared/Born_Samples_N$(N)_τ$(τ).jld", "sample_measured_states", sample_measured_states, "samples", samples, "sample_weights", sample_weights)
        # sample_measured_states,  samples, sample_weights = load("exm/data/Samples_Compared/Born_Samples_N$(N)_τ$(τ).jld", "sample_measured_states",  "samples","sample_weights")
        totalFE_enumed = load("exm/data/Boundary_measure/Born_enumed_FE_N1020.jld", "totalFE_enumed")
        FE = totalFE_enumed[i][inds]
        FElis = -log.(sample_weights)
        sample_sizes = 1:length(FElis) 
        averages = [mean(FElis[1:i]) for i in sample_sizes]
        fig = plot!(fig, xscale=:log10, 1 ./sample_sizes, averages, label="$(N)", color = c[i])
        plot!(fig, [1e-4, 1e-3], FE .*[1,1],  c=:Gray, label=false, linestyle=:dash, linewidth=2)
    end
    savefig(fig, "./exm/fig/Sample_distribution/Compare_Distribution_τ$(τ).pdf")

end

function plot_compare_distribution()
    γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999], 1.0)
    τlis = vcat(atanh.(vcat(collect(0.0:0.05:0.95), [0.99, 0.999])), 1e3)

    trajectorieslis, probabilitieslis = load("./exm/data/Boundary_measure/Born_enumed_eelis_N10.jld", "trajectorieslis", "probabilitieslis")
    for (idx, τ) in enumerate(τlis[end])
        sample_measured_states,  samples, sample_weights = load("exm/data/Samples_Compared/Born_Samples_N10_τ$(τ).jld", "sample_measured_states",  "samples","sample_weights")
        trajectories = trajectorieslis[idx]
        probabilities = probabilitieslis[idx]

        binary_digits_enum = map.(symbol -> symbol == :p ? 0 : 1, trajectories)
        decimal_value_enum = [sum(d * 2^(length(j) - i) for (i, d) in enumerate(j)) for j in binary_digits_enum]
        binary_digits_samples = map.(symbol -> symbol == :p ? 0 : 1, samples)
        decimal_value_samples = [sum(d * 2^(length(j) - i) for(i, d) in enumerate(j)) for j in binary_digits_samples]
        count_lis = [count(x->x==i, decimal_value_samples) for i in 0:31]
        @show decimal_value_enum, decimal_value_samples
        # fig = histogram(decimal_value_enum, alpha=0.3, bins=100, weights=probabilities, xlabel=L"s", ylabel=L"P", legend=false, title=latexstring("N=16, τ=$(τ)"))
        # histogram!(fig, decimal_value_samples, alpha=0.3, bins=100, normalize=:probability)
        # savefig(fig, "./exm/fig/Sample_distribution/N10_Probcompare_Distribution_τ$(τ).pdf")
    end
end

function plot_average(τ::Float64)
    fig = plot(xlabel="samples sizes", ylabel=L"\bar{S_c}", label=false, legend=:topright, legend_background_color=nothing, legend_foreground_color=nothing, )
    
    L_list = [22, 24, 26, 28, 30]
    c = cgrad(:blues, length(L_list), categorical=true)
    for (i, N) in enumerate(L_list)
        sample_measured_states, samples, sample_weights = load("./exm/data/Born_Samples_N$(N)/Born_Samples_N$(N)_τ$(τ).jld", "sample_measured_states", "samples", "sample_weights")
        FElis = -log.(sample_weights)
        sample_sizes = 1:length(FElis) 
        averages = [mean(FElis[1:i]) for i in sample_sizes]
        fig = plot!(fig, sample_sizes, averages, label="$(N)", color = c[i])
    end
    savefig(fig, "./exm/fig/average_ShannonEntropy_Distribution_τ$(τ).pdf")
end

function plot_sample(N::Int64, τ::Float64)
    # sample_measured_states, samples, sample_weights = load("./exm/data/Born_Samples_N$(N)/Born_Samples_N$(N)_τ$(τ).jld", "sample_measured_states", "samples", "sample_weights")
    energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
    antiGS = states[:,1]
    measurement_sites = collect(2:2:N)
    sample_measured_states, samples, sample_weights = Sampling(N, τ, antiGS, measurement_sites)
    num_samples=1000
    random_variable = -log.(sample_weights)
    Sc=mean(random_variable)
    stdvariance_FE = sqrt(var(random_variable))/sqrt(num_samples)

    binary_digits = map.(symbol -> symbol == :p ? 0 : 1, samples)
    decimal_value = [sum(d * 2^(length(j) - i) for (i, d) in enumerate(j)) for j in binary_digits]

    fig = histogram(decimal_value, bins=100, normalize=:probability, label="N=$(N), τ=$(τ)", xlabel=L"s", ylabel=L"P", legend=false, title=latexstring("N=$(N), \tau=$(τ)"))
    savefig(fig, "./exm/fig/Sample_Weights_Distribution_N$(N)_τ$(τ).pdf")
end

function plot_enum(N::Int64, τ::Float64)
    energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
    antiGS = states[:,1]
    measurement_sites = collect(2:2:N)
    final_states, trajectories, probabilities = measurement_enumeration(N, τ, antiGS, measurement_sites)
    Sc= sum(-log.(probabilities).*probabilities)
    binary_digits = map.(symbol -> symbol == :p ? 0 : 1, trajectories)
    decimal_value = [sum(d * 2^(length(j) - i) for (i, d) in enumerate(j)) for j in binary_digits]
    fig = histogram(decimal_value, bins=100, weights=probabilities, xlabel=L"s", ylabel=L"P", legend=false, title=latexstring("N=$(N), \tau=$(τ)"))
    savefig(fig, "./exm/fig/enumerate_distribution_N$(N)_τ$(τ).pdf")
end

# function save_boundary_postselection()
#     data = load("exm/data/post_selection/Boundary_postselection_probabilities.jld")
#     probability_m_tau_lis = data["probability_m_tau_lis"]
#     probability_p_tau_lis = data["probability_p_tau_lis"]
    
#     # Add N=30 data if it exists in the same format
#     N30_data = load("exm/data/post_selection/Boundary_postselection_N32.jld")
#     N30_m_probs = N30_data["total_weight_m_tau_lis"]
#     N30_p_probs = N30_data["total_weight_p_tau_lis"]
    
#     # Combine data
#     all_m_probs = vcat(probability_m_tau_lis, [N30_m_probs])
#     all_p_probs = vcat(probability_p_tau_lis, [N30_p_probs])
    
#     save("exm/data/post_selection/Boundary_postselection_probabilities.jld", "probability_m_tau_lis", all_m_probs, "probability_p_tau_lis", all_p_probs)
# end