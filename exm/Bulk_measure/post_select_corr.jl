using FibonacciChain
using JLD
using Statistics
using BitBasis

function compute_post_selection(L::Int64, τ::Float64, D::Int64=35L, start_point::Int64=24)
    pbc = true
    measure_class = :Fibo
    sample = ones(Int, D, length(2:2:L))

    initial_state = zeros(length(Fibonacci_basis(BitStr{L, Int}, pbc, measure_class=measure_class)))
    initial_state[1] = 1.0 # initial state is all zero state

    statelis = generate_state(τ, initial_state, sample, temp= true)

    final_st= statelis[end-3L]
    spatial_corr = spatial_correlation(L, final_st, 1, div(L,2), pbc)
    
    timeslice1 = L*start_point
    temporal_corr_lis = [temporal_correlation(τ, initial_state, sample, div(L,2), timeslice1, j) for j in timeslice1+L:D-3L]

    save("exm/data/Bulk_measure/temporal_corr/L$(L)/τ$(τ)/D$(div(D,L))_ps1.jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr", spatial_corr)
end

τ = log(1+sqrt(2)) 

if length(ARGS) == 0
    println("No arguments provided.")
else
    L=parse(Int64, ARGS[1])
    println("Received argument: $L")
    compute_post_selection(L, τ, D)
end

# temporal_corr_lis, spatial_corr = load("exm/data/Bulk_measure/temporal_corr/L$(L)/τ$(τ)/D$(div(D,L))_ps1.jld",  "temporal_corr_lis", "spatial_corr")

# plot(collect(1:length(temporal_corr_lis))./L, temporal_corr_lis ./spatial_corr, xlabel=L"\Delta t /L", ylabel=L"g(0, \Delta t)/g_{space}", label=false)

# inds = findall(x-> isapprox(x, 1.0, atol=0.1), temporal_corr_lis ./ spatial_corr)
# Δt = (collect(1:length(temporal_corr_lis))./L)[inds][end-1]
# α = log(1+√2)/π/Δt