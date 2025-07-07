using FibonacciChain
using Plots
using LaTeXStrings
using JLD
using Statistics
# include("../FitEntEntScal.jl")

function samples_generate(L::Int64, τ::Float64, D::Int64=100L)
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    samples_num = 2000

    samples_lis = Vector{Matrix{Int}}(undef, samples_num)
    sample_weights_lis = Vector{Vector{Float64}}(undef, samples_num)
    # sample_measured_states_lis = Vector{Vector{Float64}}(undef, samples_num)
    for i in 1:samples_num
        @show i
        @time sample_measured_states, samples, sample_weights = Bulkmeasure(L, τ, st, D) 
        samples_lis[i] = samples
        sample_weights_lis[i] = sample_weights
        # sample_measured_states_lis[i] = sample_measured_states
    end
    
    save("exm/data/Bulk_measure/monitored_dynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld", "samples_lis", samples_lis, "sample_weights_lis", sample_weights_lis)
    # return sample_measured_states, samples, sample_weights
end

function monitored_dynamics(L::Int64, τ::Float64, D::Int64=500L)
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    average_EElis=zeros(L-1)
    
    samples_num = 2000

    all_EE_tlis= zeros(samples_num, D) 
    stderr_EElis = zeros(L-1)
    final_EElis = zeros(samples_num, L-1)

    all_FE_tlis = zeros(samples_num, D)
    final_FElis = zeros(samples_num)
    for i in 1:samples_num
        @show i
        sample_measured_states, samples, sample_weights = Bulkmeasure(L, τ, st, D) 
        all_EE_tlis[i, :] = [eelis_Fibo_state(L, j)[div(L,2)] for j in sample_measured_states]
        final_state = sample_measured_states[end]
        final_EElis[i, :] = eelis_Fibo_state(L, final_state)

        all_FE_tlis[i, :] = sample_weights
        final_FElis[i] = sample_weights[end]
    end

    average_EElis = mean(final_EElis, dims=1)[:]
    average_EE_tlis = mean(all_EE_tlis, dims=1)[:]
    stderr_EElis = (std(final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(all_EE_tlis, dims=1) ./ sqrt(samples_num))[:]
    
    average_FE = mean(final_FElis)
    stderr_FE = std(final_FElis) / sqrt(samples_num)
    average_FE_tlis = mean(all_FE_tlis, dims=1)[:]
    stderr_FE_tlis = (std(all_FE_tlis, dims=1) ./ sqrt(samples_num))[:]

    
    return average_EE_tlis, stderr_EE_tlis, average_EElis, stderr_EElis, 
           average_FE_tlis, stderr_FE_tlis, average_FE, stderr_FE
end

function sample_calculate(L::Int64, τ::Float64, D::Int64=100L)
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    average_EElis=zeros(L-1)
    
    samples_num = 2000

    all_EE_tlis= zeros(samples_num, D) 
    stderr_EElis = zeros(L-1)
    final_EElis = zeros(samples_num, L-1)


    samples_lis, sample_weights_lis = load("./exm/data/Bulk_measure/monitored_dynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld", "samples_lis", "sample_weights_lis")
    all_FE_tlis = hcat(sample_weights_lis...)'
    final_FElis = all_FE_tlis[:, end]
    for i in 1:samples_num
        @show i
        sample_measured_states = Generate_state(τ, st, samples_lis[i], true, true)
        all_EE_tlis[i, :] = [eelis_Fibo_state(L, j)[div(L,2)] for j in sample_measured_states]
        final_state = sample_measured_states[end]
        final_EElis[i, :] = eelis_Fibo_state(L, final_state)
    end

    average_EElis = mean(final_EElis, dims=1)[:]
    average_EE_tlis = mean(all_EE_tlis, dims=1)[:]
    stderr_EElis = (std(final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(all_EE_tlis, dims=1) ./ sqrt(samples_num))[:]

    average_FE = mean(final_FElis)
    stderr_FE = std(final_FElis) / sqrt(samples_num)
    average_FE_tlis = mean(all_FE_tlis, dims=1)[:]
    stderr_FE_tlis = (std(all_FE_tlis, dims=1) ./ sqrt(samples_num))[:]

    
    save("exm/data/Bulk_measure/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(D).jld", "average_EE_tlis", average_EE_tlis, "stderr_EE_tlis", stderr_EE_tlis, "average_EElis", average_EElis, "stderr_EElis", stderr_EElis, 
    "average_FE", average_FE, "stderr_FE", stderr_FE, 
    "average_FE_tlis", average_FE_tlis, "stderr_FE_tlis", stderr_FE_tlis)
end

# L = 10
# τ = 0.1
# D = 60L
# average_EElis, EE_tlis = post_selection(L, τ, D, :m)
# fig1=plot(range(0, D/2, D), EE_tlis, label=false, xlabel=L"t", ylabel=L"S_{vN}", title="Post-selection 0")
# # plot!([1, 30*10], 0.8933161189003952*[1,1],  c=:Gray, label=false, linestyle=:dash, linewidth=2)
# cent, fig = fitCCEntEntScal(average_EElis, mincut=2, pbc=true)
# display(fig)
# save("exm/data/monitored_dynamics.jld", "cent", cent, "fig", fig)
# save("exm/data/monitored_dynamics.jld", "average_EElis", average_EElis, "EE_tlis", EE_tlis)
# savefig(fig, "/Users/cycling/Documents/projects/NoisyFibonacciChain/figs/Bulk_measure/monitored_dynamics_N28_D15_eescaling.pdf")
# savefig(fig1, "/Users/cycling/Documents/projects/NoisyFibonacciChain/figs/Bulk_measure/monitored_dynamics_N28_D15_S.pdf")
## == monitored_dynamics == ##
L = 10
τ = log(1+√2)
D = 100L
inds = collect(1:10:D)
# inds = vcat(collect())
# inds = collect(1:6000)
all_EE_tlis, average_EE_tlis, stderr_EE_tlis, all_EElis, average_EElis, stderr_EElis = monitored_dynamics(L, τ, D)
# fig_monitored_N8=plot(range(0, D/2, D)[inds], average_EE_tlis[:][inds], yerror=stderr_EE_tlis[inds], label=false, xaxis=:log, xlabel=L"t", ylabel=L"S_{vN}", title=latexstring("N=$(L), τ=$(τ)"))
# fit_samples = scatter(range(0, D/2, D)[inds], all_EE_tlis'[inds,:], xlabel=L"t", ylabel=L"S_{vN}", legend = false, colors= cgrad(:blues, length(collect(1:100))), marker_z = collect(1:100)',  colorbar=false,  line_z = collect(1:100)', title=latexstring("N=$(L), τ=$(τ)"))
# cent, fig = fitCCEntEntScal(average_EElis, mincut=1, pbc=true)
# display(fig)

if length(ARGS) == 0
    println("No arguments provided.")
else
    for arg in ARGS
        println("Received argument: $arg")
        N=parse(Int64, arg)
        # samples_generate(N, τ)
        sample_calculate(N, τ)
    end
end