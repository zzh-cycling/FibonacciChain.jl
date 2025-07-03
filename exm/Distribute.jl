using JLD
using Plots
using Statistics
using LaTeXStrings
using Arpack
using FibonacciChain

function plot()
    L_list = collect(22:2:28)
    c = cgrad(:blues, length(L_list), categorical=true)
    fig = plot(xlabel="1/ samples sizes", ylabel=L"\bar{S_c}", label=false, legend=:topright)

    for (i,N) in enumerate(L_list)
        sample_measured_states, samples, sample_weights = load("./exm/data/Born_Samples_N$(N)/Born_Samples_N$(N)_τ3.8002011672502.jld", "sample_measured_states", "samples", "sample_weights")
        FElis = -log.(sample_weights)
        sample_sizes = 1:length(FElis) 
        averages = [mean(FElis[1:i]) for i in sample_sizes]

        fig = scatter!(fig, 1 ./sample_sizes, averages, label="$(N)", color = c[i])
    end 
    savefig(fig, "./exm/fig/Inverse_Entropy_Distribution_N22_τ3.8002.pdf")

# savefig(fig, "./exm/fig/Fidelity_Entropy_Distribution_N22_τ3.8002.pdf")
end



function plot_sample(N::Int64, τ::Float64)
    # sample_measured_states, samples, sample_weights = load("./exm/data/Born_Samples_N$(N)/Born_Samples_N$(N)_τ$(τ).jld", "sample_measured_states", "samples", "sample_weights")
    energy, states = eigs(Fibonacci_Ham_sparse(N), nev=1, which=:SR)
    antiGS = states[:,1]
    measurement_sites = collect(2:2:N)
    sample_measured_states, samples, sample_weights = Sampling(N, τ, antiGS, measurement_sites)
    # num_samples=1000
    # Sc=sum(-log.(sample_weights))/num_samples
    # squared_diffs = (-log.(sample_weights) .- Sc).^2
    # variance_FE = sqrt(sum(squared_diffs) / num_samples)/sqrt(num_samples)
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