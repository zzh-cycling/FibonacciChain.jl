using FibonacciChain
using LinearAlgebra
using JLD
using Statistics
using BitBasis
using LaTeXStrings
using Plots

function get_system_params(τ, L)
    table = Dict(
            atanh(0.1)  => (2500L, 1000, 1500L),
            atanh(0.2)  => (500L,  100, 250L),
            atanh(0.3)  => (120L,  48, 100L),
            atanh(0.4)  => (100L,  40, 80L),
            atanh(0.5)  => (80L,   32, 40L),
            atanh(0.6)  => (45L,   20, 30L),
            log(1 + √2) => (35L,   14, 20L),
            atanh(0.8)  => (25L,   10, 10L),
            atanh(0.9)  => (8L,    4, 4L),
            atanh(0.95) => (8L,    4, 4L),
            atanh(0.999)=> (5L,    2, 2L),
        )
    D, step, start = get(table, τ, (5L, 2, 2L))
    inds = collect(1:step:D)
    avg_range = start:D-5
    return D, inds, avg_range
end

function compute_post_selection_Ising(L::Int64, τ::Float64, D::Int64=20L, start_point::Int64=15, sign::Int64=0)
    pbc = true
    anyon_type = :IsingZ
    sample = (sign == 1) ? ones(Int, D, L) : zeros(Int, D, L)

    initial_state = zeros(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
    initial_state[1] = 1.0 # initial state is all zero state

    statelis = generate_state(τ, initial_state, sample, temp= true, anyon_type=anyon_type)
    
    timeslice1 = L*start_point
    final_st= statelis[L*start_point]
    spatial_corr = spatial_correlation(L, final_st, 1, div(L,2), pbc=pbc, anyon_type=anyon_type)

    temporal_corr_lis = [temporal_correlation(τ, initial_state, sample, div(L,2), timeslice1, j, anyon_type=:IsingX) for j in timeslice1+2:2:timeslice1+2L]

    save("exm/data/Bulk_measure/temporal_corr_Ising/L$(L)/τ$(τ)/D$(div(D,L))_ps$(sign).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr", spatial_corr)
end

function plot_ref_ee(eelis, gamma)
    # Plot the entanglement entropy evolution
    fig = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing, 
        xlabel=L"\Delta t /L",
        ylabel=L"S_{vN}",
        title=latexstring("γ= $(round(gamma, digits=3))"),
    )
    plot!(fig, collect(1:length(eelis))./L, eelis, color=:blues, linewidth=2, label=false)

    return fig
end

function spatial_temporal_corr_varying(L::Int64, τ::Float64, D::Int64=5L, block_size::Float64=0.3, sign::Int64=0, entangle_way::Symbol=:copy)
    # | ----> |____| ----> |
    # 0       D   D+δt   D+δt+t  
    # compute how the spatial and temporal correlation changes with t, the evolution time after add two ref qubits. block_size is the time interval δt between two ref qubits divided by L
    pbc = true
    anyon_type = :IsingX # Reset to |+> state, if IsingZ, reset to |0>; If copy entangle, irrelevant.
 
    # 1). First evolve to steady state with D time steps
    sample = (sign==0) ? zeros(Int, D, L) : ones(Int, D, L)
    initial_state = ones(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
    initial_state /= norm(initial_state) # initial state is all plus state
    δt = round(Int, block_size*L)
    δt = iseven(δt) ? δt : δt - 1
    statelis = generate_state(τ, initial_state, sample, temp= true, anyon_type=anyon_type)
    
    # tlis is the time list after adding two ref qubits.
    tlis = (anyon_type == :IsingX) ? collect(0:2:D) :  collect(1:2:D)
    # spatial_corr_lis = zeros(Float64, length(tlis))
    spatial_corr = spatial_correlation(L, statelis[end], 1, div(L,2),  pbc=pbc, anyon_type=anyon_type)
    temporal_corr_lis = zeros(Float64, length(tlis))
    # eelis = zeros(Float64, length(tlis))

    # 2). Then add two ref qubits at different time slices and evolve for δt time, and to final δt + t time.
    # The entangle_way is reset, need to do measurement, then form bell pair.
    if entangle_way == :reset
        if anyon_type == :IsingZ
            savesign = (sign == 0) ? 0 : 1
        elseif anyon_type == :IsingX
            savesign = (sign == 0) ? :p : :m
        end
    
        for (idx, t) in enumerate(tlis)
            t2 = t + δt    
            ref_sample = (sign == 0) ? zeros(Int, t+δt, L) : ones(Int, t+δt, L)
        
            ref2st = reference_evolution(τ, statelis, ref_sample, 1, t, t2, anyon_type=:IsingX, verbose=true)
            sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2st, anyon_type=:IsingX, traceref = false)
            eelis[idx] = ee(sysrdm)
            temporal_corr_lis[idx] = temporal_correlation(L, ref2st, anyon_type=:IsingX)
            # spatial_correlation.(L, final_st, 1, div(L,2),  pbc=pbc, anyon_type=anyon_type)
        end
        
        # save("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(div(D,L))_ps$(savesign)_$(div(δt,L)).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr_lis", spatial_corr_lis, "eelis", eelis)
        return temporal_corr_lis, spatial_corr_lis, eelis
    # entangle_way is copy, conditioned by the given site qubit.
    elseif entangle_way == :copy
        for (idx, t) in enumerate(tlis)
            ref_sample = (sign == 0) ? zeros(Int, t+δt + D, L) : ones(Int, t+δt + D, L)
            
            ref2st = reference_evolution(τ, statelis, ref_sample, 1, D, D+δt, anyon_type=:IsingX, verbose=true)
            sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2st, anyon_type=:IsingX, traceref = false)
            # eelis[idx] = ee(sysrdm)
            temporal_corr_lis[idx] = temporal_correlation(L, ref2st, anyon_type=:IsingX)
        end

        save("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(div(D,L))_$(div(δt,L)).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr", spatial_corr)
        # return temporal_corr_lis, spatial_corr, eelis
    else
        error("Unknown entanglement way")
    end
end

# temporal_corr_lis, spatial_corr_lis, eelis = load("./exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L10/τ0.8813735870195429/D10_ps0_8_seed100.jld", "temporal_corr_lis", "spatial_corr_lis", "eelis")

# fig = plot_ref_ee(eelis, 0.707)
# plot(temporal_corr_lis, ylim=(0, maximum(temporal_corr_lis)))
# plot!(spatial_corr_lis)

function get_system_params_corr(τ)
    D = Dict(
        atanh(0.1)  => 2500,
        atanh(0.2)  => 500,
        atanh(0.3)  => 120,
        atanh(0.4)  => 100,
        atanh(0.5)  => 80,
        atanh(0.6)  => 45,
        log(1 + √2) => 35,
        atanh(0.8)  => 25,
        atanh(0.9)  => 8,
        atanh(0.95) => 8,
        atanh(0.999)=> 5,
    )
    return get(D, τ, 5)   # 5 is the default value for τ=1000.0
end

function plot_corr(L_list=collect(8:2:24))
    c = cgrad(:blues, length(L_list), categorical=true)
    
    fig = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing, 
        xlabel=L"\Delta t /L",
        ylabel=L"g(0, \Delta t)/g_{space}",
        title=latexstring("γ= $(round(gamma, digits=3))"),
    )
    # annotate!(fig_monitored_N, [(335, 3.6, text(L"L=", 10, :black))])
    for (idx, L) in enumerate(L_list)
        D = get_system_params_corr(τ)
        
        temporal_corr_lis, spatial_corr = load("exm/data/Bulk_measure/temporal_corr_Ising/L$(L)/τ$(τ)/D$(D)_ps1.jld",  "temporal_corr_lis", "spatial_corr")
        

        # plot!(fig ,collect(1:length(temporal_corr_lis))./L, temporal_corr_lis ./spatial_corr, label=latexstring("$(L)"), legendtitle=L"L", color=c[idx], linewidth=2)
        plot!(fig ,collect(1:length(temporal_corr_lis))./L, temporal_corr_lis[1:end] ./spatial_corr, label=latexstring("s=1, L=$(L)"), color=c[idx], linewidth=2)

        temporal_corr_lis, spatial_corr = load("exm/data/Bulk_measure/temporal_corr_Ising/L$(L)/τ$(τ)/D$(D)_ps0.jld",  "temporal_corr_lis", "spatial_corr")
        plot!(fig ,collect(1:length(temporal_corr_lis))./L, temporal_corr_lis[1:end] ./spatial_corr, label=latexstring("s=0, L=$(L)"), color=c[idx], linewidth=2)
    end

    return fig
end

function plot_tc_tlis(L::Int64=10, D::Int64=10, τ::Float64=log(1+√2); anyon_type::Symbol=:IsingX)
    # Plot the temporal correlations vs t for different δt
    δtlis = collect(2:2:14)
    c = cgrad(:blues, length(δtlis)+1, categorical=true)

    if anyon_type == :IsingZ
        sign0 = 0
        sign1 = 1
    elseif anyon_type == :IsingX
        sign0 = :p
        sign1 = :m
    end

    fig = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing, 
        xlabel=L"t /L",
        ylabel=L"g(0, \Delta t)/g_{space}",
        title=latexstring("γ= $(round(tanh(τ), digits=3))"),
    )
    # annotate!(fig_monitored_N, [(335, 3.6, text(L"L=", 10, :black))])
  
    for (idx, δt) in enumerate(δtlis)

        temporal_corr_lis, spatial_corr = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_$(δt).jld",  "temporal_corr_lis", "spatial_corr")

        tlis = collect(1:length(temporal_corr_lis))./L
        plot!(fig, tlis, temporal_corr_lis ./ spatial_corr, label=latexstring("(δx,δt) = (0, $(δt/L)L)"), color=c[idx+1], linewidth=2)
    end

    return fig
end

function alpha_compute_corr(L, τ)
    D = get_system_params_corr(τ)

    temporal_corr_lis, spatial_corr = load("exm/data/Bulk_measure/temporal_corr_Ising/L$(L)/τ$(τ)/D$(D)_ps1.jld",  "temporal_corr_lis", "spatial_corr")

    inds = findall(x-> isapprox(x, 1.0, atol=0.1), temporal_corr_lis ./ spatial_corr)
    Δt = (collect(1:2:length(temporal_corr_lis))./L)[inds][end-1]
    α = 2*log(1+√2)/π/Δt
    return α
end

function plot_tc(L::Int, D::Int=10, τ::Float64=log(1+√2); anyon_type::Symbol=:IsingX)
    # Plot the temporal correlations vs δt

    if anyon_type == :IsingZ
        sign0 = 0
        sign1 = 1
    elseif anyon_type == :IsingX
        sign0 = :p
        sign1 = :m
    end

    fig = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing,
        xlabel=L"δt /L",
        ylabel=L"g(0, \Delta t)/g_{space}",
        title=latexstring("γ= $(round(tanh(τ), digits=3))"),
        ylim =(-0.2, 3.0),
    )
    
    δtlis = collect(2:2:14)

    tc0lis = Vector{Float64}(undef, length(δtlis))
    sc0 = 0
    for (idx, δt) in enumerate(δtlis)
        temporal_corr_lis, spatial_corr = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_$(δt).jld",  "temporal_corr_lis", "spatial_corr")
        tc0lis[idx] = temporal_corr_lis[end]
        sc0 = spatial_corr_lis[20]
    end

    tc1lis = Vector{Float64}(undef, length(δtlis))
    sc1 = 0
    for (idx, δt) in enumerate(δtlis)
        temporal_corr_lis, spatial_corr_lis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ0.8813735870195429/D10_ps$(sign1)_$(δt)_seed100.jld",  "temporal_corr_lis", "spatial_corr_lis")
        tc1lis[idx] = temporal_corr_lis[end]
        sc1 = spatial_corr_lis[20]
    end

    plot!(fig, δtlis, tc0lis./sc0, label=L"s=0", xticks=δtlis, color=:blues, linewidth=2, marker=:circle, markersize=4)
    plot!(fig, δtlis, tc1lis./sc1, label=L"s=1", xticks=δtlis, color=:reds, linewidth=2, marker=:circle, markersize=4)

    return fig, tc0lis./sc0, tc1lis./sc1
end


function alphalis_corr(γlis)
    αlis = [alpha_compute_corr(L, τ) for (L, τ) in zip(L_list[end], τlis)]
    return αlis
end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 
gamma=tanh(log(1 + √2))
# fig = plot_corr(collect(8:2:12))

fig_corr = plot_tc_tlis(8, anyon_type= :IsingZ)
# fig, t1lis, t2lis = plot_tc(10, :IsingX)

if length(ARGS) == 0
    println("No arguments provided.")
else
    L=parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    block_size = parse(Float64, ARGS[3])
    println("Received argument: $L, $inds")
    τ = τlis[inds]
    # D, _, _ = get_system_params(τ, L)
    # compute_post_selection(L, τ, D)
    spatial_temporal_corr_varying(L, τ, 10L, block_size, 100)
end

function trace_distance(ρ1, ρ2)
    diff = ρ1 - ρ2
    # 迹距离 = 1/2 * ||ρ1 - ρ2||₁
    return 0.5 * tr(sqrt(diff' * diff))
end

function fidelity(ρ1, ρ2)
    # 对于密度矩阵，保真度定义为 F(ρ1,ρ2) = tr(√(√ρ1 * ρ2 * √ρ1))²
    # 简化计算：
    sqrt_ρ1 = sqrt(ρ1)
    F = tr(sqrt(sqrt_ρ1 * ρ2 * sqrt_ρ1))^2
    return real(F)
end

# initial_state = zeros(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
# initial_state[1] = 1.0 # initial state is all zero state
# L=8; pbc=true; anyon_type=:IsingX; D=5L;
# τ =  log(1+√2)
# initial_state = ones(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
# initial_state /= norm(initial_state) # initial state is all plus state
# stlis1 = generate_state(τ, initial_state, ones(Int, D, L), temp= true, anyon_type=anyon_type)
# stlis0 = generate_state(τ, initial_state, zeros(Int, D, L), temp= true, anyon_type=anyon_type)
# fst0 = stlis0[end]
# fst1 = stlis1[end]
# spatial_corr0 = spatial_correlation(L, fst0, 1, div(L,2),  pbc=pbc, anyon_type=anyon_type)
# spatial_corr1 = spatial_correlation(L, fst1, 1, div(L,2),  pbc=pbc, anyon_type=anyon_type)

# ref_sample0 = zeros(Int, 20L, L)
# ref_sample1 = ones(Int, 20L, L)

# ref2st0 = reference_evolution(τ, stlis0, ref_sample0, 1, 4L, 10L, anyon_type=:IsingX, verbose=true)
# ref2st1 = reference_evolution(τ, stlis1, ref_sample1, 1, 4L, 10L, anyon_type=:IsingX, verbose=true)
# sys0 = reference_rdm(L, collect(1:L), ref2st0, anyon_type=:IsingX, traceref = false)
# sys1 = reference_rdm(L, collect(1:L), ref2st1, anyon_type=:IsingX, traceref = false)

# ρeelis0 = anyon_eelis(L, sys0, pbc, anyon_type=anyon_type)
# ρeelis1 = anyon_eelis(L, sys1, pbc, anyon_type=anyon_type)

# tc0= temporal_correlation(L, ref2st0, anyon_type=:IsingX)
# tc1= temporal_correlation(L, ref2st1, anyon_type=:IsingX)

# eelis0 = anyon_eelis(L, fst0, pbc, anyon_type=anyon_type)
# eelis1 = anyon_eelis(L, fst1, pbc, anyon_type=anyon_type)

# [ρeelis0 ρeelis1 eelis0 eelis1]

add1_st0 = add_reference_qubits!(L, fst0, site, pbc=pbc,
                                   anyon_type=anyon_type)
add1_st1 = add_reference_qubits!(L, fst1, site, pbc=pbc,
                                   anyon_type=anyon_type)
add2_st0 = add_reference_qubits!(L, add1_st0, site, pbc=pbc,
                                   anyon_type=anyon_type)
add2_st1 = add_reference_qubits!(L, add1_st1, site, pbc=pbc,
                                   anyon_type=anyon_type)
add3_st0 = add_reference_qubits!(L, add2_st0, site, pbc=pbc,
                                   anyon_type=anyon_type)
temporal_correlation(L, add3_st0, anyon_type=:IsingZ)
temporal_correlation(L, add2_st1, anyon_type=:IsingZ)

# ρ1 = reference_rdm(L, [2], add2_st1, pbc=pbc, anyon_type=anyon_type)

# sys0 = reference_rdm(L, collect(1:L), add2_st0, anyon_type=:IsingX, traceref = false)
# sys1 = reference_rdm(L, collect(1:L), add2_st1, anyon_type=:IsingX, traceref = false)

# ρeelis0 = anyon_eelis(L, sys0, pbc, anyon_type=anyon_type)
# ρeelis1 = anyon_eelis(L, sys1, pbc, anyon_type=anyon_type)