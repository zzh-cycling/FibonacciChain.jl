using FibonacciChain
using LinearAlgebra
using JLD
using Statistics
using BitBasis
using LaTeXStrings
using Plots
using Random

function get_system_params(τ, L)
    table = Dict(
            atanh(0.1)  => (1200L, 1000, 1500L),
            atanh(0.2)  => (250L,  100, 250L),
            atanh(0.3)  => (60L,  48, 100L),
            atanh(0.4)  => (50L,  40, 80L),
            atanh(0.5)  => (40L,   32, 40L),
            atanh(0.6)  => (22L,   20, 30L),
            log(1 + √2) => (16L,   14, 20L),
            atanh(0.8)  => (12L,   10, 10L),
            atanh(0.9)  => (4L,    4, 4L),
            atanh(0.95) => (4L,    4, 4L),
        )
    D, step, start = get(table, τ, (3L, 2, 2L))
    inds = collect(1:step:D)
    avg_range = start:D-5
    return D, inds, avg_range
end

# Here need to note that t is div(D,2), where D is the total evolution layer of get_system_params
function compute_post_selection_Ising(L::Int64, τ::Float64, t::Int64=5L, δt::Int=1; sign::Bool=false, entangle_way::Symbol=:copy)
    model = AnyonModel(IsingAnyon(), L; pbc=true, measure_operator=:X)
    sample = sign ? BitMatrix(ones(Bool, 2t, L)) : BitMatrix(zeros(Bool, 2t, L))

    initial_state = ones(length(anyon_basis(model)))
    initial_state /= norm(initial_state) # initial state is all zero state

    steady_config = MeasureConfig(τ=τ, mode=:sample, t₂=t, enable_τ_eff=false)
    mo = bulk_evolution(model, initial_state, steady_config, sample)
    statelis, Flis = mo.states, mo.free_energys

    ref_sample = sign ? BitMatrix(ones(Bool, 2*(t+δt+t), L)) : BitMatrix(zeros(Bool, 2*(t+δt+t), L))
       
    if entangle_way == :copy
        if δt == 0
            ref_config = MeasureConfig(τ=τ, mode=:sample, x₂ = L÷2+1, t₂=t, t₁ = t, verbose=true)
            ref2stlis, sample_layer, sample_free_energy = reference_evolution(model, statelis,ref_config, ref_sample) # to compute temporal correlation, add ref qubit at site L/2+1
            # ref2st, F= reference_apply_measurement_layer!(L, τ/2, ref2stlis[end], zeros(Int, L), D+1, pbc, anyon_type=:IsingX, extended_basis=ext_basis, k_old=2)
            ref2st = ref2stlis[end]
            spatial = true
            temporal = false
        else
            ref_config = MeasureConfig(τ=τ, mode=:sample, x₂ = L÷2+1, t₂=t + δt, t₁ = t, verbose=true, x₁ = L÷2+1)
            ref2stlis, sample_layer, sample_free_energy = reference_evolution(model, statelis,ref_config, ref_sample) # to compute temporal correlation, add ref qubit at site L/2+1
            # ref2st, F= reference_apply_measurement_layer!(L, τ/2, ref2stlis[end], zeros(Int, L), D+1, pbc, anyon_type=:IsingX, extended_basis=ext_basis, k_old=2)
            ref2st = ref2stlis[end]
            temporal = true
            spatial = false
        end
        spatial_corr, temporal_corr = ref_correlation(model, ref2st, spatial = spatial, temporal = temporal)
        sysrdm = reference_rdm(model, collect(1:div(L,2)), ref2st, traceref = false)
        S = ee(sysrdm)
    end

    save("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(div(2t,L))_ps$(sign)_$(δt).jld", "temporal_corr", temporal_corr, "spatial_corr", spatial_corr, "S", S)
    return temporal_corr, spatial_corr, S
end

function spatial_temporal_corr_varyingt(L::Int64, τ::Float64, t::Int64=5L, δt::Int=1; sign::Bool=false, entangle_way::Symbol=:copy)
    # | ----> |____| ----> |
    # 0       t   t+δt   t+δt+t  
    # compute how the spatial and temporal correlation changes with t, the evolution time after add two ref qubits. block_size is the time interval δt between two ref qubits divided by L
    model = AnyonModel(IsingAnyon(), L; pbc=true, measure_operator=:X) # Reset to |+> state, if IsingZ, reset to |0>; If copy entangle, irrelevant.
 
    # 1). First evolve to steady state with D time steps
    sample = sign ? BitMatrix(ones(Bool, 2t, L)) : BitMatrix(zeros(Bool, 2t, L))
    initial_state = ones(length(anyon_basis(model)))
    initial_state /= norm(initial_state) # initial state is all plus state
    
    config = MeasureConfig(τ=τ, mode=:sample, t₂=t, enable_τ_eff=false)
    mo = bulk_evolution(model, initial_state, config, sample)
    statelis, Flis = mo.states, mo.free_energys

    mp = model.measure_operator
    # tlis is the time list after adding two ref qubits.
    Dlis = (mp == :X) ? collect(0:2:2t) :  collect(1:2:2t)
    spatial_corr_lis = zeros(Float64, length(Dlis))
    temporal_corr_lis = zeros(Float64, length(Dlis))
    eelis = zeros(Float64, length(Dlis))

    # 2). Then add two ref qubits at different time slices and evolve for δt time, and to final δt + t time.
    # The entangle_way is reset, need to do measurement, then form bell pair.
    if entangle_way == :reset
        if mp == :Z
            savesign = sign ? 1 : 0
        elseif mp == :X
            savesign = sign ? :p : :m
        end
    
        for (idx, t) in enumerate(Dlis)
            t2 = t + δt    
            ref_sample = sign ? BitMatrix(ones(Int, t+δt, L)) : BitMatrix(zeros(Int, t+δt, L))

            ref_config = MeasureConfig(τ=τ, mode=:sample, x₂ = L÷2+1, t₂=t2, t₁ = t, verbose=true)
            ref_mo = reference_evolution(model, statelis,ref_config, ref_sample) # to compute temporal correlation, add ref qubit at site L/2+1
            ref2stlis, sample, sample_free_energy = ref_mo.states, ref_mo.samples, ref_mo.free_energys # to compute temporal correlation, add ref qubit at site 1
            sysrdm = reference_rdm(model, collect(1:div(L,2)), ref2stlis[end], traceref = false)
            eelis[idx] = ee(sysrdm)
            temporal_corr_lis[idx] = temporal_correlation(model, ref2stlis[end])
            # spatial_correlation.(L, final_st, 1, div(L,2),  pbc=pbc, anyon_type=anyon_type)

        end
        
        # save("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(div(D,L))_ps$(savesign)_$(div(δt,L)).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr_lis", spatial_corr_lis, "eelis", eelis)
        return temporal_corr_lis, spatial_corr_lis, eelis
    # entangle_way is copy, conditioned by the given site qubit.
    elseif entangle_way == :copy
        for (idx, t) in enumerate(Dlis)
            ref_sample = sign ? BitMatrix(ones(Int, t+δt + D, L)) : BitMatrix(zeros(Int, t+δt + D, L))

            if δt == 0
                ref_config = MeasureConfig(τ=τ, mode=:sample, x₂ = L÷2+1, t₂=t, t₁ = t, verbose=true)
                ref_mo = reference_evolution(model, statelis,ref_config, ref_sample) # to compute temporal correlation, add ref qubit at site L/2+1
                ref2stlis, sample_layer, sample_free_energy = ref_mo.states, ref_mo.samples, ref_mo.free_energys  # to compute temporal correlation, add ref qubit at site L/2+1
                # ref2st, F= reference_apply_measurement_layer!(L, τ/2, ref2stlis[end], zeros(Int, L), D+1, pbc, anyon_type=:IsingX, extended_basis=ext_basis, k_old=2)
                ref2st = ref2stlis[end]
                spatial = true
                temporal = false
            else
                ref_config = MeasureConfig(τ=τ, mode=:sample, x₂ = L÷2+1, t₂=t + δt, t₁ = t, verbose=true, x₁ = L÷2+1)
                ref_mo = reference_evolution(model, statelis,ref_config, ref_sample)
                ref2stlis, sample_layer, sample_free_energy = ref_mo.states, ref_mo.samples, ref_mo.free_energys  # to compute temporal correlation, add ref qubit at site L/2+1
                temporal = true
                spatial = false
            end
            sysrdm = reference_rdm(model, collect(1:div(L,2)), ref2stlis[end], traceref = false)
            eelis[idx] = ee(sysrdm)
            spatial_corr, temporal_corr = ref_correlation(model, ref2stlis[end], spatial=spatial, temporal=temporal)
            temporal_corr_lis[idx] = temporal_corr
            spatial_corr_lis[idx] = spatial_corr 
        end

        save("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(div(D,L))_ps$(sign)_$(δt).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr_lis", spatial_corr_lis, "eelis", eelis)
        return temporal_corr_lis, spatial_corr_lis, eelis
    end
end

function organize()
    Llis = collect(8:2:16)
    δtlis = collect(2:2:16)
    scLlis = zeros(Float64, length(Llis))
    tcLlis = zeros(Float64, length(Llis), length(δtlis))
    for (i, L) in enumerate(Llis)
        τ = log(1 + √2)
        D = 8
        spatial_corr = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_ps1_0.jld", "spatial_corr")
        scLlis[i] = spatial_corr
        for (j, δt) in enumerate(δtlis)
            temporal_corr = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_ps1_$(δt).jld",  "temporal_corr")
            tcLlis[i, j] = temporal_corr
        end
    end
    save("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/stc_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", scLlis, "tcLlis", tcLlis)
end

function get_system_params_corr(τ)
    D = Dict(
        atanh(0.1)  => 1200,
        atanh(0.2)  => 250,
        atanh(0.3)  => 60,
        atanh(0.4)  => 50,
        atanh(0.5)  => 40,
        atanh(0.6)  => 22,
        log(1 + √2) => 16,
        atanh(0.8)  => 12,
        atanh(0.9)  => 4,
        atanh(0.95) => 4,
    )
    return get(D, τ, 3)   # 5 is the default value for τ=1000.0
end

function plot_stc_tlis(L::Int64=10, D::Int64=10, τ::Float64=log(1+√2); sign::Int=0)
    # Plot the spatio-temporal correlations vs t for different δt
    δtlis = collect(1:6)
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

function plot_tc(L::Int, D::Int=10, τ::Float64=log(1+√2); sign::Int=1)
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
        
        
    δtlis = collect(1:6)
    @show δtlis

    tc1lis = zeros(length(δtlis))
    sc1 = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_ps$(sign)_0.jld", "spatial_corr")
    for (idx, δt) in enumerate(δtlis)
        temporal_corr= load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/L$(L)/τ$(τ)/D$(D)_ps$(sign)_$(δt).jld",  "temporal_corr")
        tc1lis[idx] = temporal_corr
    end

    plot!(fig, δtlis./L, tc1lis./sc1, label=L"s=1", xticks=δtlis./L, color=:reds, linewidth=2, marker=:circle, markersize=4)

    t_star = log(1 + √2)/ π
    plot!(fig, t_star*[1, 1],[minimum(tc1lis./sc1), maximum(tc1lis./sc1)], linestyle=:dash, color=:gray, label = false)
    plot!(fig, [0.05, 0.75], [1,1], linestyle=:dash, color=:gray, label = false) # horizontal line
    scatter!(fig, [t_star], [1], color=:black, marker=:star5, markersize=8, label=false) # 
    annotate!(fig, t_star+0.1, 1+0.15, text(L"(t^*=\frac{\log(1+\sqrt{2})}{\pi}, g/g=1)", 10, :black))

    return fig
end

function plot_stc_scaling(τ::Float64=log(1+√2))
    Llis = collect(8)
    δtlis = collect(1:8)

    scLlis, tcLlis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Ising/stc_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", "tcLlis")
    
    c = cgrad(:blues, length(Llis), categorical=true)
    fig = plot(
        legend_background_color=nothing,
        legend_foreground_color=nothing,
        xlabel=L"δt/L",
        ylabel=L"g(0, \delta t)/g(\delta x=L/2, 0)",
        title="spacetime self-dual",
        )

    for (i, L) in enumerate(Llis[1:end-1])    
        scatter!(fig, δtlis./(L), tcLlis[i, :]./scLlis[i], label=latexstring("L=$(L)"), color=c[i], marker=:circle, markersize=4)
    end

    plot!(fig, δtlis./(Llis[end]), tcLlis[end, :]./scLlis[end], label=latexstring("L=$(Llis[end])"), color=c[end], linewidth=2, marker=:circle, markersize=4)

    t_star = log(1 + √2)/ π
    plot!(fig, t_star*[1, 1],[minimum(tcLlis./scLlis), maximum(tcLlis./scLlis)], linestyle=:dash, color=:gray, label = false) # vertical line
    plot!(fig, [0.05, 0.75], [1,1], linestyle=:dash, color=:gray, label = false) # horizontal line
    scatter!(fig, [t_star], [1], color=:black, marker=:star5, markersize=8, label=false)
    annotate!(fig, t_star+0.15, 1+0.15, text(L"(t^*=\frac{\log(1+\sqrt{2})}{\pi}, g/g=1)", 8, :black))
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

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)


# fig_corr = plot_stc_tlis(10, anyon_type= :IsingZ, sign=0)
# fig= plot_tc(11, anyon_type= :IsingX)

if length(ARGS) == 0
    println("No arguments provided.")
else
    L=parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    δt = parse(Int, ARGS[3])
    println("Received argument: $L, $inds, $δt")
    τ = τlis[inds]
    # D, _, _ = get_system_params(τ, L)
    compute_post_selection_Ising(L, τ, 5L, δt, sign=1)
    # spatial_temporal_corr_varyingt(L, τ, 10L, δt, sign=1)
end

# L=8; pbc=true; anyon_type=:IsingX; D=5L;
# τ =  log(1+√2)
# initial_state = ones(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
# initial_state /= norm(initial_state) # initial state is all plus state
# stlis1, F = generate_state(τ, initial_state, ones(Int, D, L),  anyon_type=anyon_type)
# stlis0, F = generate_state(τ, initial_state, zeros(Int, D, L), anyon_type=anyon_type)
# fst0 = stlis0[end]
# fst1 = stlis1[end]

# ref_sample0 = zeros(Int, D+50, L)
# ref_sample1 = ones(Int, D, L)
# using Profile
# @profile ref2st0 = reference_evolution(L, τ, stlis0, ref_sample0, L÷2, D, D+4, anyon_type=:IsingX, verbose=true, rng= MersenneTwister(100))
# using BenchmarkTools
# @btime reference_evolution(L, τ, stlis0, ref_sample0, L÷2, D, D+4, anyon_type=:IsingX, verbose=false,rng= MersenneTwister(100))
# @code_warntype ref2st0 = reference_evolution(L, τ, stlis0, ref_sample0, L÷2, D, D+4, anyon_type=:IsingX, verbose=true)
# ref2st1 = reference_evolution(L, τ, stlis1, ref_sample1, L÷2, D, D, anyon_type=:IsingX, verbose=true)
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