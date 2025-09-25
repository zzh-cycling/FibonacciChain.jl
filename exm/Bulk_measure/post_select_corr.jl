using FibonacciChain
using LinearAlgebra
using JLD
using Statistics
using BitBasis
using LaTeXStrings
using Plots
using Random

function get_δtL(τ)
    table = Dict(
            atanh(0.1)  => (collect(25:35)),
            atanh(0.2)  => (collect(10:20)),
            atanh(0.3)  => (collect(5:15)),
            atanh(0.4)  => (collect(2:10)),)
    δtlis = get(table, τ, collect(1:8))
    return  δtlis
end

function organize(τ::Float64, sign::Int64=1)
    δtlis = get_δtL(τ)
    Llis = collect(8:4:20)
    scLlis = zeros(Float64, length(Llis))
    tcLlis = zeros(Float64, length(Llis), length(δtlis))
    for (i, L) in enumerate(Llis)
        D = get_system_params_corr(τ)
        spatial_corr = load("exm/data/Bulk_measure/spatial_temporal_corr/L$(L)/τ$(τ)/D$(D)_ps$(sign)_0.jld", "spatial_corr")
        scLlis[i] = spatial_corr
        for (j, δt) in enumerate(δtlis)
            temporal_corr = load("exm/data/Bulk_measure/spatial_temporal_corr/L$(L)/τ$(τ)/D$(D)_ps$(sign)_$(δt).jld",  "temporal_corr")
            tcLlis[i, j] = temporal_corr
        end
    end
    save("exm/data/Bulk_measure/spatial_temporal_corr/stc$(sign)_τ$(τ)_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", scLlis, "tcLlis", tcLlis)
end

function get_system_params(τ, L)
    table = Dict(
            atanh(0.1)  => (300L, 1000, 1500L),
            atanh(0.2)  => (60L,  100, 250L),
            atanh(0.3)  => (25L,  48, 100L),
            atanh(0.4)  => (20L,  40, 80L),
            atanh(0.5)  => (20L,   32, 40L),
            atanh(0.6)  => (15L,   20, 30L),
            log(1 + √2) => (10L,   14, 20L),
            atanh(0.8)  => (8L,   10, 10L),
            atanh(0.9)  => (5L,    4, 4L),
            atanh(0.95) => (5L,    4, 4L),
        )
    D, step, start = get(table, τ, (3L, 2, 2L))
    inds = collect(1:step:D)
    avg_range = start:D-5
    return D, inds, avg_range
end

function compute_post_selection(L::Int64, τ::Float64, D::Int64=5L, δt::Int64=2; sign::Int64=0, entangle_way::Symbol=:copy, rng = MersenneTwister(100))
    pbc = true
    sample = (sign == 1) ? ones(Int, 2*D, length(2:2:L)) : zeros(Int, 2*D, length(2:2:L))

    initial_state = zeros(length(anyon_basis(BitStr{L, Int}, pbc)))
    initial_state[1] = 1.0 # initial state is all zero state

    statelis, Flis = generate_state(τ, initial_state, sample, enable_τ_eff=false)
    ref_sample = (sign == 0) ? zeros(Int, 2*(D+δt+D), length(2:2:L)) : ones(Int, 2*(D+δt+D), length(2:2:L))

    if entangle_way == :copy
        if δt == 0
            ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D, verbose=true) # to compute temporal correlation, add ref qubit at site L/2+1
            spatial = true
            temporal = false
        else
            ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D+δt, x₁ = L÷2+1, verbose=true) # to compute temporal correlation, add ref qubit at site L/2+1
            temporal = true
            spatial = false
        end
        spatial_corr, temporal_corr = ref_correlation(L, ref2stlis[end], spatial = spatial, temporal = temporal)
        sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2stlis[end], traceref = false)
        S = ee(sysrdm)
    end

    save("exm/data/Bulk_measure/spatial_temporal_corr/L$(L)/τ$(τ)/D$(div(D,L))_ps$(sign)_$(δt).jld", "temporal_corr", temporal_corr, "spatial_corr", spatial_corr, "S", S)
    return temporal_corr, spatial_corr, S

end

function spatial_temporal_corr_varyingt(L::Int64, τ::Float64, D::Int64=5L, δt::Int=1; sign::Int64=0, entangle_way::Symbol=:copy)
    # | ----> |____| ----> |
    # 0       D   D+δt   D+δt+t  
    # compute how the spatial and temporal correlation changes with t, the evolution time after add two ref qubits. δt is the time interval between two ref qubits
    pbc = true
 
    # 1). First evolve to steady state with D time steps
    sample = (sign==0) ? zeros(Int, 2*D, length(2:2:L)) : ones(Int, 2*D, length(2:2:L))
    initial_state = zeros(length(anyon_basis(BitStr{L, Int}, pbc)))
    initial_state[1] = 1.0 # initial state is all zero state

    statelis, Flis = generate_state(τ, initial_state, sample, enable_τ_eff=false)

    # tlis is the time list after adding two ref qubits.
    tlis = collect(1:D)
    spatial_corr_lis = zeros(Float64, length(tlis))
    temporal_corr_lis = zeros(Float64, length(tlis))
    eelis = zeros(Float64, length(tlis))

    basis_F = anyon_basis(L, pbc)
    extended_basis = FibonacciChain.build_extended_basis(2, basis_F) 

    # 2). Then add two ref qubits at different time slices and evolve for δt time, and to final δt + t time.
    if entangle_way == :copy
        
        for (idx, t) in enumerate(tlis)
            @show "calculation time t:" t
            ref_sample = (sign == 0) ? zeros(Int, 2*(t + δt + D), length(2:2:L)) : ones(Int, 2*(t + δt + D), length(2:2:L))
            if δt == 0
                ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D) 
                spatial = true
                temporal = false
            else
                ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D+δt, x₁ = L÷2+1)
                temporal = true
                spatial = false
            end

            sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2stlis[end], traceref = false)
            eelis[idx] = ee(sysrdm)
            spatial_corr, temporal_corr = ref_correlation(L, ref2stlis[end], spatial=spatial, temporal=temporal)
            temporal_corr_lis[idx] = temporal_corr
            spatial_corr_lis[idx] = spatial_corr
        end

        save("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(div(D,L))_ps$(sign)_$(δt).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr_lis", spatial_corr_lis, "eelis", eelis)
        return temporal_corr_lis, spatial_corr_lis, eelis
    end
end

function get_system_params_corr(τ)
    D = Dict(
        atanh(0.1)  => 300,
        atanh(0.2)  => 60,
        atanh(0.3)  => 25,
        atanh(0.4)  => 20,
        atanh(0.5)  => 20,
        atanh(0.6)  => 15,
        log(1 + √2) => 10,
        atanh(0.8)  => 8,
        atanh(0.9)  => 5,
        atanh(0.95) => 5,
    )
    return get(D, τ, 3)   # 3 is the default value for τ=1000.0
end

function plot_stc_tlis(L::Int64=10, D::Int64=3, τ::Float64=log(1+√2); sign::Int=0)
    # Plot the spatio-temporal correlations vs t for different δt
    δtlis = get_δtL(τ)
    c = cgrad(:blues, length(δtlis)+1, categorical=true)
    
    fig = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing, 
        xlabel=L"t /L",
        ylabel=L"g(0, \Delta t), g_{space}",
        title=latexstring("γ= $(round(tanh(τ), digits=3))"),
    )
    
    temporal_corr_lis, spatial_corr_lis, eelis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(D)_ps$(sign)_0.jld",  "temporal_corr_lis", "spatial_corr_lis", "eelis")

    tlis = collect(1:length(temporal_corr_lis))./L
    plot!(fig, tlis, spatial_corr_lis, label=latexstring("(δx,δt) = (L/2, 0)"), color=c[1], linestyle=:dash, linewidth=2)
    
    for (idx, δt) in enumerate(δtlis)
        temporal_corr_lis, spatial_corr_lis, eelis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(D)_ps$(sign)_$(δt).jld",  "temporal_corr_lis", "spatial_corr_lis", "eelis")

        tlis = collect(1:length(temporal_corr_lis))./L
        plot!(fig, tlis, temporal_corr_lis, label=latexstring("(δx,δt) = (0, $(δt/L)L)"), color=c[idx+1], linewidth=2)
    end

    return fig
end
    
function plot_ref_ee(eelis, gamma, L)
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

function plot_tc(L::Int, D::Int=10, τ::Float64=log(1+√2))
    # Plot the temporal correlations vs δt
    δtlis = get_δtL(τ)

    fig = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing,
        xlabel=L"δt /L",
        ylabel=L"g(0, \Delta t)/g_{space}",
        title=latexstring("γ= $(round(tanh(τ), digits=3))"),
    )
    
    # tc0lis = Vector{Float64}(undef, length(δtlis))
    # spatial_corr_lis0 = load("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(D)_ps0_0.jld", "spatial_corr_lis")
    # sc0 = spatial_corr_lis0[end]
    # for (idx, δt) in enumerate(δtlis)
    #     temporal_corr_lis, spatial_corr_lis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(D)_ps0_$(δt).jld",  "temporal_corr_lis", "spatial_corr_lis")
    #     tc0lis[idx] = temporal_corr_lis[end]
    # end

    tc1lis = Vector{Float64}(undef, length(δtlis))
    spatial_corr_lis1 = load("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(D)_ps1_0.jld", "spatial_corr_lis")
    sc1 = spatial_corr_lis1[end]
    for (idx, δt) in enumerate(δtlis)
        temporal_corr_lis, spatial_corr_lis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(D)_ps1_$(δt).jld",  "temporal_corr_lis", "spatial_corr_lis")
        tc1lis[idx] = temporal_corr_lis[end]
    end

    # plot!(fig, δtlis./L, tc0lis./sc0, label=L"s=0", xticks=δtlis./L, color=:blues, linewidth=2, marker=:circle, markersize=4)
    plot!(fig, δtlis./L, tc1lis./sc1, label=L"s=1", xticks=δtlis./L, color=:reds, linewidth=2, marker=:circle, markersize=4)

    return fig, tc0lis./sc0, tc1lis./sc1
end

function plot_stc_scaling(τ::Float64=log(1+√2); sign::Int=1)
    Llis = collect(8:4:20)
    δtlis = get_δtL(τ)

    idx = findall(x->x==τ, τlis)
    tlis = tlis = [1.5622791153697357, 0.7780581224158777, 0.5083827557773009, 0.36709846127189905, 0.2760268365466788, 0.21191890015652085, 0.1568936137623848, 0.11830701021871359, 0.07607223477091314, 0.10374056220848625] # t obtained from α fitting
    scLlis, tcLlis = load("exm/data/Bulk_measure/spatial_temporal_corr/stc$(sign)_τ$(τ)_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", "tcLlis")

    c = cgrad(:blues, length(Llis), categorical=true)
    fig = plot(
        legend_background_color=nothing,
        legend_foreground_color=nothing,
        xlabel=L"δt/L",
        ylabel=L"g(0, \delta t)/g(\delta x=L/2, 0)",
        title=latexstring("γ= $(round(tanh(τ), digits=3)), s=$(sign)"),
        ylim = (0.97, 1.03),
        )

    for (i, L) in enumerate(Llis[1:end-1])    
        scatter!(fig, δtlis./(L), tcLlis[i, :]./scLlis[i], label=latexstring("L=$(L)"), color=c[i], marker=:circle, markersize=4)
    end

    plot!(fig, δtlis./(Llis[end]), tcLlis[end, :]./scLlis[end], label=latexstring("L=$(Llis[end])"), color=c[end], linewidth=2, marker=:circle, markersize=4)

    # for (i, L) in enumerate(Llis[1:end])    
    #     plot!(fig, δtlis./(2L), tcLlis[i, :]./scLlis[i], label=latexstring("L=$(L)"), color=c[i], marker=:circle, markersize=4, linewidth=2)
    # end
    scatter!(fig, tlis[idx],[1], label=false, color=:black, marker=:star5, markersize=6) # α fitting points
    plot!(fig, [0.05, 0.75], [1,1], linestyle=:dash, color=:gray, label = false) # horizontal line
    return fig, tcLlis./scLlis  
end

function alpha_compute_corr(τ; sign::Int=1)
    Llis = collect(8:4:20)
    δtlis= get_δtL(τ)

    scLlis, tcLlis = load("exm/data/Bulk_measure/spatial_temporal_corr/stc$(sign)_τ$(τ)_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", "tcLlis")
    
    ratio = tcLlis./scLlis
    inds = findall(x-> isapprox(x, 1.0, atol=0.1), ratio[end,:])
    
    linear_model(x,p) = p[1] * x .+ p[2]

    if !isempty(inds)
        fit = curve_fit(linear_model, δtlis[inds]./Llis[end], ratio[end, inds], [1.0, 1.0])
        a = fit.param[1]
        b = fit.param[2]
        t = (1-b)/a
        α = log(1+√2)/π/t
        return α,t
    else
        return NaN
    end
end

function alphalis_corr(τlis)
    αlis = [alpha_compute_corr(τ)[1] for (τ) in τlis]
    tlis = [alpha_compute_corr(τ)[2] for (τ) in τlis]
    @show tlis
    return αlis
end


# αlis = [0.1795773389079671, 0.36057708040946856, 0.551847840984768, 0.7642361812075126, 1.0163864125658897, 1.3238551444084465, 1.7881538925764218, 2.3713719554820867, 3.687941165583592]

# trueαlis = [0.17479998016347537, 0.35335053506014624, 0.5379878111388026, 0.7513193004900436, 0.9973854818280089, 1.2954457344112627, 1.7342095249297775, 2.324623740539361, 3.6859611811798443, 5.660898631724869, 70.46901458043763, -3.291284646075844e13]

# error = [2.7330430701558868, 2.0451491174598737, 2.5762720937165486, 1.7192265271295484, 1.9050739241818382, 2.1930220033566306, 3.11060266197238, 2.011001355938974, 0.05371690873624521] * 100%

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 


if length(ARGS) == 0
    println("No arguments provided.")
else
    L = parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    δt = parse(Int, ARGS[3])
    τ = τlis[inds]
    D, _, _ = get_system_params(τ, L)
    println("Computed spatial_temporal_corr_varyingt for L=$L, τ=$τ, D=$D, δt=$δt")
    spatial_temporal_corr_varyingt(L, τ, D, δt, sign=1)
    # compute_post_selection(L, τ, D, δt, sign=1)
end
