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

function compute_post_selection_Ising(L::Int64, τ::Float64, D::Int64=5L, δt::Int=2; sign::Int64=0, entangle_way::Symbol=:copy)
    pbc = true
    anyon_type = :IsingX
    sample = (sign == 1) ? ones(Int, D, L) : zeros(Int, D, L)

    initial_state = ones(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
    initial_state /= norm(initial_state) # initial state is all zero state

    statelis = generate_state(τ, initial_state, sample, temp= true, anyon_type=anyon_type)

    ref_sample = (sign == 0) ? zeros(Int, D+δt+D, L) : ones(Int, D+δt+D, L)
        
    if entangle_way == :copy
        ref3st = reference_evolution(τ, statelis, ref_sample, L÷2, D, D+δt, anyon_type=:IsingX, verbose=true)
        sysrdm = reference_rdm(L, collect(1:div(L,2)), ref3st, anyon_type=:IsingX, traceref = false)
        S = ee(sysrdm)
        spatio = (δt == 0) ? true : false
        spatial_corr, temporal_corr = ref_correlation(L, ref3st, anyon_type=:IsingX, spatio=spatio)
    end

    save("exm/data/Bulk_measure/temporal_corr_Ising/L$(L)/τ$(τ)/D$(div(D,L))_ps$(sign).jld", "temporal_corr", temporal_corr, "spatial_corr", spatial_corr, "S", S)
    return temporal_corr, spatial_corr, S
end

function spatial_temporal_corr_varyingt(L::Int64, τ::Float64, D::Int64=5L, δt::Int=2; sign::Int64=0, entangle_way::Symbol=:copy)
    # | ----> |____| ----> |
    # 0       D   D+δt   D+δt+t  
    # compute how the spatial and temporal correlation changes with t, the evolution time after add two ref qubits. block_size is the time interval δt between two ref qubits divided by L
    pbc = true
    anyon_type = :IsingX # Reset to |+> state, if IsingZ, reset to |0>; If copy entangle, irrelevant.
 
    # 1). First evolve to steady state with D time steps
    sample = (sign==0) ? zeros(Int, D, L) : ones(Int, D, L)
    initial_state = ones(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
    initial_state /= norm(initial_state) # initial state is all plus state
    
    δt = iseven(δt) ? δt : δt - 1
    statelis = generate_state(τ, initial_state, sample, temp= true, anyon_type=anyon_type)
    
    # tlis is the time list after adding two ref qubits.
    tlis = (anyon_type == :IsingX) ? collect(0:2:D) :  collect(1:2:D)
    spatial_corr_lis = zeros(Float64, length(tlis))
    temporal_corr_lis = zeros(Float64, length(tlis))
    eelis = zeros(Float64, length(tlis))

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
        
            ref3st = reference_evolution(τ, statelis, ref_sample, L÷2, D, D+δt, anyon_type=:IsingX, verbose=true)
            sysrdm = reference_rdm(L, collect(1:div(L,2)), ref3st, anyon_type=:IsingX, traceref = false)
            eelis[idx] = ee(sysrdm)
            spatio = (δt == 0) ? true : false
            spatial_corr, temporal_corr = ref_correlation(L, ref3st, anyon_type=:IsingX, spatio=spatio)
            temporal_corr_lis[idx] = temporal_corr
            spatial_corr_lis[idx] = spatial_corr
        end

        save("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(div(D,L))_ps$(sign)_$(δt).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr_lis", spatial_corr_lis, "eelis", eelis)
        return temporal_corr_lis, spatial_corr_lis, eelis
    else
        error("Unknown entanglement way")
    end
end


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

function plot_stc_tlis(L::Int64=10, D::Int64=10, τ::Float64=log(1+√2); anyon_type::Symbol=:IsingX, sign::Int=0)
    # Plot the spatio-temporal correlations vs t for different δt
    δtlis = collect(2:2:10)
    c = cgrad(:blues, length(δtlis)+1, categorical=true)
    
    fig = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing, 
        xlabel=L"t /L",
        ylabel=L"g(0, \Delta t), g_{space}",
        title=latexstring("γ= $(round(tanh(τ), digits=3))"),
    )
    # annotate!(fig_monitored_N, [(335, 3.6, text(L"L=", 10, :black))])
    temporal_corr_lis, spatial_corr_lis, eelis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_ps$(sign)_0.jld",  "temporal_corr_lis", "spatial_corr_lis", "eelis")

    tlis = collect(1:length(temporal_corr_lis))./L
    plot!(fig, tlis, spatial_corr_lis, label=latexstring("(δx,δt) = (L/2, 0)"), color=c[1], linestyle=:dash, linewidth=2)

    for (idx, δt) in enumerate(δtlis)

        temporal_corr_lis, spatial_corr_lis, eelis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_ps$(sign)_$(δt).jld",  "temporal_corr_lis", "spatial_corr_lis", "eelis")

        plot!(fig, tlis, temporal_corr_lis, label=latexstring("(δx,δt) = (0, $(δt/2L)L)"), color=c[idx+1], linewidth=2)
    end

    return fig
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

function plot_tc(L::Int, D::Int=10, τ::Float64=log(1+√2); anyon_type::Symbol=:IsingX)
    # Plot the temporal correlations vs δt
    
    fig = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing,
        xlabel=L"δt /L",
        ylabel=L"g(0, \Delta t)/g_{space}",
        title=latexstring("γ= $(round(tanh(τ), digits=3))"),
        # ylim =(-0.2, 3.0),
        )
        
        
    δtlis = collect(2:2:10)
    # tc0lis = Vector{Float64}(undef, length(δtlis))
    # spatial_corr_lis0 = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_ps0_0.jld", "spatial_corr_lis")
    # sc0 = spatial_corr_lis0[end]
    # for (idx, δt) in enumerate(δtlis)
    #     temporal_corr_lis, spatial_corr_lis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_ps0_$(δt).jld",  "temporal_corr_lis", "spatial_corr_lis")
    #     tc0lis[idx] = temporal_corr_lis[end]
    # end
    @show δtlis
    tc1lis = zeros(length(δtlis))
    spatial_corr_lis1 = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_ps1_0.jld", "spatial_corr_lis")
    sc1 = spatial_corr_lis1[end]
    for (idx, δt) in enumerate(δtlis)
        temporal_corr_lis, spatial_corr_lis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ0.8813735870195429/D10_ps1_$(δt).jld",  "temporal_corr_lis", "spatial_corr_lis")
        tc1lis[idx] = temporal_corr_lis[end]
    end

    # plot!(fig, δtlis./2L, tc0lis./sc0, label=L"s=0", xticks=δtlis./2L, color=:blues, linewidth=2, marker=:circle, markersize=4)
    plot!(fig, δtlis./2L, tc1lis./sc1, label=L"s=1", xticks=δtlis./2L, color=:reds, linewidth=2, marker=:circle, markersize=4)

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


function alphalis_corr(γlis)
    αlis = [alpha_compute_corr(L, τ) for (L, τ) in zip(L_list[end], τlis)]
    return αlis
end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 
gamma=tanh(log(1 + √2))


# fig_corr = plot_stc_tlis(10, anyon_type= :IsingZ, sign=0)
# fig= plot_tc(12, anyon_type= :IsingX)

if length(ARGS) == 0
    println("No arguments provided.")
else
    L=parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    δt = parse(Int, ARGS[3])
    println("Received argument: $L, $inds, $δt")
    τ = τlis[inds]
    # D, _, _ = get_system_params(τ, L)
    compute_post_selection_Ising(L, τ, 8L, δt, sign=0)
    # spatial_temporal_corr_varyingt(L, τ, 10L, δt, sign=1)
end

# L=8; pbc=true; anyon_type=:IsingX; D=5L;
# τ =  log(1+√2)
# initial_state = ones(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
# initial_state /= norm(initial_state) # initial state is all plus state
# stlis1 = generate_state(τ, initial_state, ones(Int, D, L), temp= true, anyon_type=anyon_type)
# stlis0 = generate_state(τ, initial_state, zeros(Int, D, L), temp= true, anyon_type=anyon_type)
# fst0 = stlis0[end]
# fst1 = stlis1[end]

# ref_sample0 = zeros(Int, D+50, L)
# ref_sample1 = ones(Int, D, L)
# using Profile
# @profile ref2st0 = reference_evolution(τ, stlis0, ref_sample0, L÷2, D, D+4, anyon_type=:IsingX, verbose=true)
# @code_warntype ref2st0 = reference_evolution(τ, stlis0, ref_sample0, L÷2, D, D+4, anyon_type=:IsingX, verbose=true)
# ref2st1 = reference_evolution(τ, stlis1, ref_sample1, L÷2, D, D, anyon_type=:IsingX, verbose=true)
# sys0 = reference_rdm(L, collect(1:L), ref2st0, anyon_type=:IsingX, traceref = false)
# sys1 = reference_rdm(L, collect(1:L), ref2st1, anyon_type=:IsingX, traceref = false)

# ρ1 = reference_rdm(L, [2], ref2st1, pbc=pbc, anyon_type=anyon_type)
# ρ2 = reference_rdm(L, [1], ref2st1, pbc=pbc, anyon_type=anyon_type)
# ρ12 = reference_rdm(L, [1, 2], ref2st1, pbc=pbc, anyon_type=anyon_type)
# ee(ρ1), ee(ρ2), ee(ρ12)
# I = ee(ρ1) + ee(ρ2) - ee(ρ12)

# ρeelis0 = anyon_eelis(L, sys0, pbc, anyon_type=anyon_type)
# ρeelis1 = anyon_eelis(L, sys1, pbc, anyon_type=anyon_type)

# tc0= temporal_correlation(L, ref2st0, anyon_type=:IsingX)
# tc1= temporal_correlation(L, ref2st1, anyon_type=:IsingX)

# eelis0 = anyon_eelis(L, fst0, pbc, anyon_type=anyon_type)
# eelis1 = anyon_eelis(L, fst1, pbc, anyon_type=anyon_type)

# [ρeelis0 ρeelis1 eelis0 eelis1]

# add1_st0 = add_reference_qubits!(L, fst0, L÷2, pbc=pbc,
#                                    anyon_type=anyon_type)
# add1_st1 = add_reference_qubits!(L, fst1, L÷2, pbc=pbc,
#                                    anyon_type=anyon_type)
# add2_st0 = add_reference_qubits!(L, add1_st0, L÷2, pbc=pbc,
#                                    anyon_type=anyon_type)
# add2_st1 = add_reference_qubits!(L, add1_st1, L÷2, pbc=pbc,
#                                    anyon_type=anyon_type)
# add3_st0 = add_reference_qubits!(L, add2_st0, L÷2, pbc=pbc,
#                                    anyon_type=anyon_type)
# temporal_correlation(L, add3_st0, anyon_type=:IsingZ)
# temporal_correlation(L, add2_st1, anyon_type=:IsingZ)

# ρ1 = reference_rdm(L, [2], add2_st1, pbc=pbc, anyon_type=anyon_type)

# sys0 = reference_rdm(L, collect(1:L), add2_st0, anyon_type=:IsingX, traceref = false)
# sys1 = reference_rdm(L, collect(1:L), add2_st1, anyon_type=:IsingX, traceref = false)

# ρeelis0 = anyon_eelis(L, sys0, pbc, anyon_type=anyon_type)
# ρeelis1 = anyon_eelis(L, sys1, pbc, anyon_type=anyon_type)

# state1 = add_reference_qubits!(L, fst0, 1, pbc=pbc,
#                                    anyon_type=anyon_type)
# state2 = add_reference_qubits!(L, state1, L÷2, pbc=pbc,
#                                    anyon_type=anyon_type)
# final_stlis2 = reference_generate_state(τ, state2, ref_sample0[D+1:end, :], pbc, anyon_type=anyon_type, temp=false)
# ρ = reference_rdm(L, collect(1:L), final_stlis2, anyon_type=:IsingX, traceref = false)
# temporal_correlation(L, final_stlis2, anyon_type=:IsingX)
# sc = spatial_correlation(L, ρ, 1, div(L,2),  pbc=pbc, anyon_type=anyon_type)

# ref_correlation(L, final_stlis2, anyon_type=:IsingX)

# state_addref3 = final_stlis2
# ρ1 = reference_rdm(N, [3], state_addref3, pbc=pbc, anyon_type=anyon_type)
# ρ2 = reference_rdm(N, [2], state_addref3, pbc=pbc, anyon_type=anyon_type) 
# ρ3 = reference_rdm(N, [1], state_addref3, pbc=pbc, anyon_type=anyon_type)
# ρ12 = reference_rdm(N, [2, 3], state_addref3, pbc=pbc, anyon_type=anyon_type)
# ρ23 = reference_rdm(N, [1, 2], state_addref3, pbc=pbc, anyon_type=anyon_type)
# ee(ρ1), ee(ρ2), ee(ρ3), ee(ρ12), ee(ρ23)