using FibonacciChain
using LinearAlgebra
using JLD
using Statistics
using BitBasis
using LaTeXStrings
using Plots

function get_δtL(τ)
    table = Dict(
            atanh(0.1)  => (collect(2:2:70), collect(8:4:20)),
            atanh(0.2)  => (collect(2:2:40), collect(8:4:20)),
            atanh(0.3)  => (collect(2:2:30), collect(8:4:20)),
            atanh(0.4)  => (collect(2:2:20), collect(8:4:20)),
            atanh(0.5)  => (collect(2:2:16), collect(8:4:20)),
            atanh(0.6)  => (collect(2:2:12), collect(8:4:20)),
        )
    δtlis, Llis = get(table, τ, (collect(2:2:10), collect(8:4:20)))
    return  δtlis, Llis
end

function organize(τ::Float64, sign::Int64=1)
    δtlis, Llis = get_δtL(τ)
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
            atanh(0.1)  => (600L, 1000, 1500L),
            atanh(0.2)  => (125L,  100, 250L),
            atanh(0.3)  => (50L,  48, 100L),
            atanh(0.4)  => (40L,  40, 80L),
            atanh(0.5)  => (40L,   32, 40L),
            atanh(0.6)  => (30L,   20, 30L),
            log(1 + √2) => (20L,   14, 20L),
            atanh(0.8)  => (15L,   10, 10L),
            atanh(0.9)  => (10L,    4, 4L),
            atanh(0.95) => (10L,    4, 4L),
            atanh(0.999)=> (5L,    2, 2L),
        )
    D, step, start = get(table, τ, (5L, 2, 2L))
    inds = collect(1:step:D)
    avg_range = start:D-5
    return D, inds, avg_range
end

function compute_post_selection(L::Int64, τ::Float64, D::Int64=10L, δt::Int64=2; sign::Int64=0, entangle_way::Symbol=:copy)
    pbc = true
    sample = (sign == 1) ? ones(Int, D, length(2:2:L)) : zeros(Int, D, length(2:2:L))

    initial_state = zeros(length(anyon_basis(BitStr{L, Int}, pbc)))
    initial_state[1] = 1.0 # initial state is all zero state

    statelis = generate_state(τ, initial_state, sample, temp= true)
    ref_sample = (sign == 0) ? zeros(Int, D+δt+D, length(2:2:L)) : ones(Int, D+δt+D, length(2:2:L))

    if entangle_way == :copy
        if δt == 0
            ref2st = reference_evolution(τ, statelis, ref_sample, L÷2+1, D, D, verbose=true) # to compute temporal correlation, add ref qubit at site L/2+1
            spatial = true
            temporal = false
        else
            ref2st = reference_evolution(τ, statelis, ref_sample, L÷2+1, D, D+δt, verbose=true, x₁ = L÷2+1) # to compute temporal correlation, add ref qubit at site L/2+1
            temporal = true
            spatial = false
        end
        spatial_corr, temporal_corr = ref_correlation(L, ref2st, spatial = spatial, temporal = temporal)
        sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2st, traceref = false)
        S = ee(sysrdm)
    end

    save("exm/data/Bulk_measure/spatial_temporal_corr/L$(L)/τ$(τ)/D$(div(D,L))_ps$(sign)_$(δt).jld", "temporal_corr", temporal_corr, "spatial_corr", spatial_corr, "S", S)
    return temporal_corr, spatial_corr, S

end

function spatial_temporal_corr_varyingt(L::Int64, τ::Float64, D::Int64=20L, δt::Int=2; sign::Int64=0, entangle_way::Symbol=:copy)
    # | ----> |____| ----> |
    # 0       D   D+δt   D+δt+t  
    # compute how the spatial and temporal correlation changes with t, the evolution time after add two ref qubits. δt is the time interval between two ref qubits
    pbc = true
 
    # 1). First evolve to steady state with D time steps
    sample = (sign==0) ? zeros(Int, D, length(2:2:L)) : ones(Int, D, length(2:2:L))
    initial_state = zeros(length(anyon_basis(BitStr{L, Int}, pbc)))
    initial_state[1] = 1.0 # initial state is all zero state
    
    δt = iseven(δt) ? δt : δt - 1
    statelis = generate_state(τ, initial_state, sample, temp= true)
    
    # tlis is the time list after adding two ref qubits.
    tlis = collect(0:2:D)
    spatial_corr_lis = zeros(Float64, length(tlis))
    temporal_corr_lis = zeros(Float64, length(tlis))
    eelis = zeros(Float64, length(tlis))

    # 2). Then add two ref qubits at different time slices and evolve for δt time, and to final δt + t time.
    if entangle_way == :copy
        for (idx, t) in enumerate(tlis)
            ref_sample = (sign == 0) ? zeros(Int, t+δt + D, length(2:2:L)) : ones(Int, t+δt + D, length(2:2:L))
        
            if δt == 0
                ref2st = reference_evolution(τ, statelis, ref_sample, L÷2+1, D, D, verbose=true) 
                spatial = true
                temporal = false
            else
                ref2st = reference_evolution(τ, statelis, ref_sample, L÷2+1, D, D+δt, verbose=true, x₁ = L÷2+1)
                temporal = true
                spatial = false
            end
            sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2st, traceref = false)
            eelis[idx] = ee(sysrdm)
            spatial_corr, temporal_corr = ref_correlation(L, ref2st, spatial=spatial, temporal=temporal)
            temporal_corr_lis[idx] = temporal_corr
            spatial_corr_lis[idx] = spatial_corr
        end

        save("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(div(D,L))_ps$(sign)_$(δt).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr_lis", spatial_corr_lis, "eelis", eelis)
        return temporal_corr_lis, spatial_corr_lis, eelis
    else
        error("Unknown entanglement way")
    end
end

function get_system_params_corr(τ)
    D = Dict(
        atanh(0.1)  => 600,
        atanh(0.2)  => 125,
        atanh(0.3)  => 50,
        atanh(0.4)  => 40,
        atanh(0.5)  => 40,
        atanh(0.6)  => 30,
        log(1 + √2) => 20,
        atanh(0.8)  => 15,
        atanh(0.9)  => 10,
        atanh(0.95) => 10,
        atanh(0.999)=> 5,
    )
    return get(D, τ, 5)   # 5 is the default value for τ=1000.0
end

function plot_stc_tlis(L::Int64=10, D::Int64=10, τ::Float64=log(1+√2); sign::Int=0)
    # Plot the spatio-temporal correlations vs t for different δt
    δtlis = collect(2:2:12)
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
    δtlis = collect(2:2:10)

    fig = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing,
        xlabel=L"δt /L",
        ylabel=L"g(0, \Delta t)/g_{space}",
        title=latexstring("γ= $(round(tanh(τ), digits=3))"),
    )
    
    tc0lis = Vector{Float64}(undef, length(δtlis))
    spatial_corr_lis0 = load("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(D)_ps0_0.jld", "spatial_corr_lis")
    sc0 = spatial_corr_lis0[end]
    for (idx, δt) in enumerate(δtlis)
        temporal_corr_lis, spatial_corr_lis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(D)_ps0_$(δt).jld",  "temporal_corr_lis", "spatial_corr_lis")
        tc0lis[idx] = temporal_corr_lis[end]
    end

    tc1lis = Vector{Float64}(undef, length(δtlis))
    spatial_corr_lis1 = load("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(D)_ps1_0.jld", "spatial_corr_lis")
    sc1 = spatial_corr_lis1[end]
    for (idx, δt) in enumerate(δtlis)
        temporal_corr_lis, spatial_corr_lis = load("exm/data/Bulk_measure/spatial_temporal_corr_varying/L$(L)/τ$(τ)/D$(D)_ps1_$(δt).jld",  "temporal_corr_lis", "spatial_corr_lis")
        tc1lis[idx] = temporal_corr_lis[end]
    end

    plot!(fig, δtlis./L, tc0lis./sc0, label=L"s=0", xticks=δtlis./L, color=:blues, linewidth=2, marker=:circle, markersize=4)
    plot!(fig, δtlis./L, tc1lis./sc1, label=L"s=1", xticks=δtlis./L, color=:reds, linewidth=2, marker=:circle, markersize=4)

    return fig, tc0lis./sc0, tc1lis./sc1
end

function plot_stc_scaling(τ::Float64=log(1+√2); sign::Int=1)
    Llis = collect(8:4:20)
    δtlis = get_δtL(τ)[1]

    scLlis, tcLlis = load("exm/data/Bulk_measure/spatial_temporal_corr/stc$(sign)_τ$(τ)_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", "tcLlis")

    c = cgrad(:blues, length(Llis), categorical=true)
    fig = plot(
        legend_background_color=nothing,
        legend_foreground_color=nothing,
        xlabel=L"δt/L",
        ylabel=L"g(0, \delta t)/g(\delta x=L/2, 0)",
        title=latexstring("γ= $(round(tanh(τ), digits=3)), s=$(sign)"),
        ylim=(0.95, 1.05),
        )

    for (i, L) in enumerate(Llis[1:end-1])    
        scatter!(fig, δtlis./(2L), tcLlis[i, :]./scLlis[i], label=latexstring("L=$(L)"), color=c[i], marker=:circle, markersize=4)
    end

    plot!(fig, δtlis./(2Llis[end]), tcLlis[end, :]./scLlis[end], label=latexstring("L=$(Llis[end])"), color=c[end], linewidth=2, marker=:circle, markersize=4)

    # for (i, L) in enumerate(Llis[1:end])    
    #     plot!(fig, δtlis./(2L), tcLlis[i, :]./scLlis[i], label=latexstring("L=$(L)"), color=c[i], marker=:circle, markersize=4, linewidth=2)
    # end

    # t_star = log(1 + √2)/ π
    # plot!(fig, t_star*[1, 1],[minimum(tcLlis./scLlis), maximum(tcLlis./scLlis)], linestyle=:dash, color=:gray, label = false) # vertical line
    plot!(fig, [0.05, 0.75], [1,1], linestyle=:dash, color=:gray, label = false) # horizontal line
    # scatter!(fig, [t_star], [1], color=:black, marker=:star5, markersize=8, label=false)
    # annotate!(fig, t_star+0.15, 1+0.15, text(L"(t^*=\frac{\log(1+\sqrt{2})}{\pi}, g/g=1)", 8, :black))
    return fig, tcLlis./scLlis  
end

function alpha_compute_corr(τ)
    δtlis, Llis = get_δtL(τ)

    scLlis, tcLlis = load("exm/data/Bulk_measure/spatial_temporal_corr/stc$(sign)_τ$(τ)_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", "tcLlis")
    
    ratio = tcLlis./scLlis
    inds = findall(x-> isapprox(x, 1.0, atol=0.1), ratio[end,:])
    
    linear_model(x,p) = p[1] * x .+ p[2]

    if !isempty(inds)
        fit = curve_fit(linear_model, δtlis[inds]./2/Llis[end], ratio[end, inds], [1.0, 1.0])
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

# 添加辅助函数

# αlis = [0.1727892970244856, 0.3400420163610316, 0.522784271836154, 0.6872600192202489, 0.9094895390374313, 1.1169076314891018, 1.410107539026646, 1.8083920281560268, 2.658491243484519, 3.6391966478869002]

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
    # spatial_temporal_corr_varyingt(L, τ, D, δt, sign=1)
    compute_post_selection(L, τ, D, δt, sign=1)
end
