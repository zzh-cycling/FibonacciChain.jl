using FibonacciChain
using JLD
using Statistics
using BitBasis


τ = log(1+sqrt(2)) 

function compute_total(L::Int64, τ::Float64, index::Int64, D::Int64=35L)
    for index in index:index+99
        @time compute_ratio(L, τ, index, D)
    end
end

function compute_ratio(L::Int64, τ::Float64, index::Int64, D::Int64=35L)
    pbc = true
    measure_class = :Fibo
    sample = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$L/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "sample")

    initial_state = zeros(length(Fibonacci_basis(BitStr{L, Int}, pbc, measure_class=measure_class)))
    initial_state[1] = 1.0 # initial state is all zero state

    statelis = generate_state(τ, initial_state, sample, temp= true)

    final_st= statelis[end-20]
    spatial_corr = spatial_correlation(L, final_st, 1, div(L,2), pbc)
    
    timeslice1 = div(N, L)*8
    temporal_corr_lis = [temporal_correlation(τ, initial_state, sample, div(L,2), timeslice1, j) for j in timeslice1+5:D-20]

    save("exm/data/Bulk_measure/temporal_corr/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr", spatial_corr)
end

function corr_collect(L::Int64, τ::Float64, D::Int64=35L)
    samples_num = 10000
    temporal_corr_ensemble = Vector{Vector{Float64}}(undef, samples_num)
    spatial_corr_ensemble = Vector{Float64}(undef, samples_num)
     for i in 1:samples_num
        @show i
        temporal_corr_lis, spatial_corr = load("exm/data/Bulk_measure/temporal_corr/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld",  "temporal_corr_lis", "spatial_corr")
        temporal_corr_ensemble[i] = temporal_corr_lis
        spatial_corr_ensemble[i] = spatial_corr
    end

    average_temporal_corr = mean(temporal_corr_ensemble)
    average_spatial_corr = mean(spatial_corr_ensemble)
    temporal_corr_stderr = std(temporal_corr_ensemble) ./ sqrt(samples_num)
    spatial_corr_stderr = std(spatial_corr_ensemble) / sqrt(samples_num)

    save("exm/data/Bulk_measure/temporal_spatial_corr_L$(L)_τ$(τ)_D$(div(D,L)).jld", "average_temporal_corr", average_temporal_corr, 
    "temporal_corr_stderr", temporal_corr_stderr, 
    "average_spatial_corr", average_spatial_corr, 
    "spatial_corr_stderr", spatial_corr_stderr)
end
# 0.8s per run

# plot(collect(1:length(temporal_corr_lis))./L, temporal_corr_lis ./spatial_corr, xlabel=L"\Delta t /L", ylabel=L"g(0, \Delta t)/g_{space}", label=false)



if length(ARGS) == 0
    println("No arguments provided.")
else
    L=parse(Int64, ARGS[1])
    index=parse(Int64, ARGS[2])
    println("Received argument: $L, $index")
    # compute_ratio(L, τ, index)
    # compute_total(L, τ, index)
    corr_collect(L, τ)
end