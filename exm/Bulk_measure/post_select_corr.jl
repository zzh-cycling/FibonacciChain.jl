using FibonacciChain
using JLD
using Statistics
using BitBasis

function get_system_params(τ, L)
    if τ == log(1 + √2)
        D = 35L
        inds = collect(1:14:D)
        avg_range = 20L:D-5
    elseif τ == atanh(0.1)
        D = 2500L
        inds = collect(1:1000:D)
        avg_range = 1500L:D-5
    elseif τ == atanh(0.2)
        D = 500L
        inds = collect(1:100:D)
        avg_range = 250L:D-5
    elseif τ == atanh(0.3)
        D = 120L
        inds = collect(1:48:D)
        avg_range = 100L:D-5
    elseif τ == atanh(0.4)
        D = 100L
        inds = collect(1:40:D)
        avg_range = 80L:D-5
    elseif τ == atanh(0.5)
        D = 80L
        inds = collect(1:32:D)
        avg_range = 40L:D-5
    elseif τ == atanh(0.6)
        D = 45L
        inds = collect(1:20:D)
        avg_range = 30L:D-5
    elseif τ == atanh(0.8)
        D = 25L
        inds = collect(1:10:D)
        avg_range = 10L:D-5
    elseif τ == atanh(0.9) || τ == atanh(0.95)
        D = 8L
        inds = collect(1:4:D)
        avg_range = 4L:D-5
    elseif τ == atanh(0.999)
        D = 5L
        inds = collect(1:2:D)
        avg_range = 2L:D-5
    else
        D = 5L  # Default value for τ=1000.0
        inds = collect(1:2:D)
        avg_range = 2L:D-5
    end
    
    return D, inds, avg_range
end

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
    temporal_corr_lis = [temporal_correlation(τ, initial_state, sample, div(L,2), timeslice1, j) for j in timeslice1+5:timeslice1+2L]

    save("exm/data/Bulk_measure/temporal_corr/L$(L)/τ$(τ)/D$(div(D,L))_ps1.jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr", spatial_corr)
end

# function get_system_params_corr(τ)
#     if τ == log(1 + √2)
#         D = 35
#     elseif τ == atanh(0.1)
#         D = 2500
#     elseif τ == atanh(0.2)
#         D = 500
#     elseif τ == atanh(0.3)
#         D = 120
#     elseif τ == atanh(0.4)
#         D = 100
#     elseif τ == atanh(0.5)
#         D=80
#     elseif τ == atanh(0.6)
#         D=45
#     elseif τ == atanh(0.8)
#         D=25
#     elseif τ == atanh(0.9) || τ == atanh(0.95)
#         D=8
#     elseif τ == atanh(0.999)
#         D=5
#     else
#         D = 5  # Default value for τ=1000.0
#     end
#     return D
# end

# function plot_corr(L_list=collect(8:2:24))
#     c = cgrad(:blues, length(L_list), categorical=true)
    
#     fig = plot(
#         legend_background_color=nothing,
#         legend_foreground_color=nothing, 
#         xlabel=L"\Delta t /L",
#         ylabel=L"g(0, \Delta t)/g_{space}"
#     )
#     # annotate!(fig_monitored_N, [(335, 3.6, text(L"L=", 10, :black))])
#     for (idx, L) in enumerate(L_list)
#         D = get_system_params_corr(τ)
        
#         temporal_corr_lis, spatial_corr = load("exm/data/Bulk_measure/temporal_corr/L$(L)/τ$(τ)/D$(D)_ps1.jld",  "temporal_corr_lis", "spatial_corr")
        

#         plot!(fig ,collect(1:length(temporal_corr_lis))./L, temporal_corr_lis ./spatial_corr, label=latexstring("$(L)"), legendtitle=L"L", color=c[idx], linewidth=2)
#     end

#     return fig
# end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 


if length(ARGS) == 0
    println("No arguments provided.")
else
    L=parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    println("Received argument: $L, $inds")
    τ = τlis[inds]
    D, _, _ = get_system_params(τ, L)
    compute_post_selection(L, τ, D)
end



# τ = log(1 + √2)  
# fig = plot_corr(collect(8:4:12))
# inds = findall(x-> isapprox(x, 1.0, atol=0.1), temporal_corr_lis ./ spatial_corr)
# Δt = (collect(1:length(temporal_corr_lis))./L)[inds][end-1]
# α = log(1+√2)/π/Δt