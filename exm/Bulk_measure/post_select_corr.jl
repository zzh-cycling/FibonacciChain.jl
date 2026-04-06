using FibonacciChain
using LinearAlgebra
using JLD
using Statistics
using BitBasis
using LaTeXStrings
using Plots
using Random
using LsqFit

function get_δtL(τ, sign::Int64=1)
    if sign == 0
        table = Dict(
            atanh(0.1)  => (collect(10:30)),
            atanh(0.2)  => (collect(10:30)),
            atanh(0.3)  => (collect(10:30)),
            atanh(0.4)  => (collect(10:30)),)
        δtlis = get(table, τ, collect(1:10))
    else
        table = Dict(
                atanh(0.1)  => (collect(25:35)),
                atanh(0.2)  => (collect(10:20)),
                atanh(0.3)  => (collect(5:15)),
                atanh(0.4)  => (collect(2:10)),)
        δtlis = get(table, τ, collect(1:8))
    end
    return  δtlis
end

function organize(τ::Float64, sign::Int64=1)
    δtlis = get_δtL(τ, sign)
    Llis = (sign==1) ? collect(8:2:20) : collect(6:6:24)
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

function compute_post_selection(L::Int64, τ::Float64, t::Int64=5L, δt::Int64=2; sign::Bool=false, entangle_way::Symbol=:copy)
    model = AnyonModel(FibonacciAnyon(), L; pbc=true)
    sample = sign ? BitMatrix(ones(Int, 2t, length(2:2:L))) : BitMatrix(zeros(Int, 2t, length(2:2:L)))

    initial_state = zeros(length(anyon_basis(model)))
    initial_state[1] = 1.0 # initial state is all zero state

    pre_config = MeasureConfig(τ=τ, mode=:sample, t₂=div(t,2), enable_τ_eff=false)
    pre_mo = bulk_evolution(model, initial_state, pre_config, sample)
    Flis = pre_mo.free_energys
    ref_sample = sign ? BitMatrix(ones(Int, 2*(t + δt + t), length(2:2:L))) : BitMatrix(zeros(Int, 2*(t + δt + t), length(2:2:L)))

    if entangle_way == :copy
        if δt == 0
            ref_config = MeasureConfig(τ = τ, t₂ = t, t₁ = t, mode=:sample, x₂=L÷2+1, verbose = true)
            spatial = true
            temporal = false
        else
            ref_config = MeasureConfig(τ = τ, t₂ = t + δt, t₁ = t, mode=:sample, x₂=L÷2+1, x₁ = L÷2+1, verbose = true)
            temporal = true
            spatial = false
        end

        ref_mo = reference_evolution(model, pre_mo.state, ref_config, ref_sample)
        sample_layer, sample_free_energy = ref_mo.samples, ref_mo.free_energys # to compute temporal correlation, add ref qubit at site L/2+1

        spatial_corr, temporal_corr = ref_correlation(L, ref_mo.state, spatial = spatial, temporal = temporal)
        sysrdm = reference_rdm(L, collect(1:div(L,2)), ref_mo.state, traceref = false)
        S = ee(sysrdm)
    end

    save("exm/data/Bulk_measure/spatial_temporal_corr/L$(L)/τ$(τ)/D$(div(D,L))_ps$(sign)_$(δt).jld", "temporal_corr", temporal_corr, "spatial_corr", spatial_corr, "S", S)
    return temporal_corr, spatial_corr, S

end

function spatial_temporal_corr_varyingt(L::Int64, τ::Float64, t::Int64=5L, δt::Int=1; sign::Int64=0, entangle_way::Symbol=:copy)
    # | ----> |____| ----> |
    # 0       t   t+δt   t+δt+t  
    # compute how the spatial and temporal correlation changes with t, the evolution time after add two ref qubits. δt is the time interval between two ref qubits
    model = AnyonModel(FibonacciAnyon(), L; pbc=true)
 
    # 1). First evolve to steady state with D time steps
    sample = sign ? BitMatrix(ones(Int, 2t, length(2:2:L))) : BitMatrix(zeros(Int, 2t, length(2:2:L)))
    initial_state = zeros(length(anyon_basis(model)))
    initial_state[1] = 1.0 # initial state is all zero state

    pre_config = MeasureConfig(τ=τ, mode=:sample, t₂=div(t,2), enable_τ_eff=false)
    pre_mo = bulk_evolution(model, initial_state, pre_config, sample)
    Flis = pre_mo.free_energys

    # tlis is the time list after adding two ref qubits.
    tlis = collect(1:2t)
    spatial_corr_lis = zeros(Float64, length(tlis))
    temporal_corr_lis = zeros(Float64, length(tlis))
    eelis = zeros(Float64, length(tlis))
    
    # 2). Then add two ref qubits at different time slices and evolve for δt time, and to final δt + t time.
    if entangle_way == :copy
        
        for (idx, Δt) in enumerate(tlis)
            @show "calculation time t:" t
            ref_sample = sign ? BitMatrix(ones(Int, 2*(t+δt + Δt), L)) : BitMatrix(zeros(Int, 2*(Δt+δt + t), L))
            if δt == 0
                ref_config = MeasureConfig(τ=τ, mode=:sample, x₂ = L÷2+1, t₂=t, t₁ = t, verbose=true)
                spatial = true
                temporal = false
            else
                ref_config = MeasureConfig(τ=τ, mode=:sample, x₂ = L÷2+1, t₂=t + δt, t₁ = t, verbose=true, x₁ = L÷2+1)
                temporal = true
                spatial = false
            end
            
            ref_mo = reference_evolution(model, pre_mo.state,ref_config, ref_sample) # to compute temporal correlation, add ref qubit at site L/2+1
            sample_layer, sample_free_energy = ref_mo.samples, ref_mo.free_energys  # to compute temporal correlation, add ref qubit at site L/2+1

            sysrdm = reference_rdm(model, collect(1:div(L,2)), ref_mo.state, traceref = false)
            eelis[idx] = ee(sysrdm)
            spatial_corr, temporal_corr = ref_correlation(model, ref_mo.state, spatial=spatial, temporal=temporal)
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
    δtlis = get_δtL(τ, sign)
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
    δtlis = get_δtL(τ, sign)

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
    Llis = (sign == 1) ? collect(8:4:20) : collect(6:6:24)
    δtlis = get_δtL(τ, sign)

    idx = findall(x->x==τ, τlis)
    tlis = (sign==1) ? [1.555921499844019, 0.7656122752643163, 0.4989803597691102, 0.35908439421100435, 0.26702996240796666, 0.207503905626826, 0.15730972912766123, 0.1200023145127127, 0.07378287804757423, 0.04704355851034897, 0.0, 0.0] : [2.37017, 1.1654823245070127, 0.7631518032413759, 0.5555641984246366, 0.43075671034670565, 0.3429090132464811, 0.27246080615905355, 0.22173881849605775, 0.16639103643814843, 0.13907413514019382, 0.09403357048715472, 0.08707261034310568] # t obtained from α fitting
    scLlis, tcLlis = load("exm/data/Bulk_measure/spatial_temporal_corr/stc$(sign)_τ$(τ)_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", "tcLlis")

    c = cgrad(:blues, length(Llis), categorical=true)
    fig = plot(
        legend_background_color=nothing,
        legend_foreground_color=nothing,
        xlabel=L"δt/L",
        ylabel=L"g(0, \delta t)/g(\delta x=L/2, 0)",
        title=latexstring("γ= $(round(tanh(τ), digits=3)), s=$(sign)"),
        # ylim = (0.97, 1.03),
        )

    for (i, L) in enumerate(Llis)    
        scatter!(fig, δtlis./(L), tcLlis[i, :]./scLlis[i], label=latexstring("L=$(L)"), color=c[i], marker=:circle, markersize=4)
    end

    # index = 6
    # plot!(fig, δtlis./(Llis[index]), tcLlis[index, :]./scLlis[index], label=latexstring("L=$(Llis[index])"), color=c[index], linewidth=2, marker=:circle, markersize=4)

    # for (i, L) in enumerate(Llis[1:end])    
    #     plot!(fig, δtlis./(2L), tcLlis[i, :]./scLlis[i], label=latexstring("L=$(L)"), color=c[i], marker=:circle, markersize=4, linewidth=2)
    # end
    scatter!(fig, tlis[idx],[1], label=false, color=:black, marker=:star5, markersize=6) # α fitting points
    plot!(fig, [0.05, 0.75], [1,1], linestyle=:dash, color=:gray, label = false) # horizontal line
    return fig, tcLlis./scLlis  
end


function alpha_compute_corr(τ; sign::Int=1)
    Llis = (sign == 1) ? collect(8:4:20) : collect(6:6:24)
    δtlis= get_δtL(τ, sign)

    scLlis, tcLlis = load("exm/data/Bulk_measure/spatial_temporal_corr/stc$(sign)_τ$(τ)_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", "tcLlis")
    
    ratio = tcLlis./scLlis
    inds = findall(x-> isapprox(x, 1.0, atol=0.1), ratio[end,:])
    
    # linear_model(x,p) = p[1] * x .+ p[2]

    # if !isempty(inds)
    #     fit = curve_fit(linear_model, δtlis[inds]./Llis[end], ratio[end, inds], [1.0, 1.0])
    #     a = fit.param[1]
    #     b = fit.param[2]
    #     t = (1-b)/a
    #     α = log(1+√2)/π/t
    #     return α,t
    # else
    #     return NaN
    # end
    model(x, p) = @. p[1] .* exp.(-x ./p[2])  + p[3]
    fit = curve_fit(model, δtlis./Llis[end], ratio[end, :], [1.0, 1.0, 1.0])
    a, b, c = fit.param
    t = -b * log((1-c)/a)  # exponential fit

    α = log(1 + √2) / π / t

    return α, t
end

function alphalis_corr(τlis; sign::Int=1)
    αlis = [alpha_compute_corr(τ; sign=sign)[1] for (τ) in τlis]
    tlis = [alpha_compute_corr(τ; sign=sign)[2] for (τ) in τlis]
    @show tlis
    return αlis
end

# αlis = [0.1795773389079671, 0.36057708040946856, 0.551847840984768, 0.7642361812075126, 1.0163864125658897, 1.3238551444084465, 1.7881538925764218, 2.3713719554820867, 3.687941165583592]

# trueαlis1 = [0.18031110579660684, 0.3664386468630409, 0.5622464304996033, 0.7812924501662804, 1.0499308866604552, 1.352022388793634, 1.7834238716533293, 2.3978709594794476, 3.8023716828814287, 5.963620420165977, 0.0, 0.0]

# trueαlis0 = [0.118367005813756, 0.24071572796117682, 0.36762007896462434, 0.504982011017125, 0.6512955444008806, 0.8181468416752649, 1.02968911427875, 1.2652269371345004, 1.686087977905452, 2.0172688896233755, 2.9835081739017246, 3.2220226896161237]

# error = [2.7330430701558868, 2.0451491174598737, 2.5762720937165486, 1.7192265271295484, 1.9050739241818382, 2.1930220033566306, 3.11060266197238, 2.011001355938974, 0.05371690873624521] * 100%

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)


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
    compute_post_selection(L, τ, D, δt, sign=0)
end
