using FibonacciChain
using LinearAlgebra
using JLD
using Statistics
using BitBasis
using LaTeXStrings
using Plots

function organize(sign::Int64=1)
    Llis = collect(8:2:16)
    δtlis = collect(2:2:16)
    scLlis = zeros(Float64, length(Llis))
    tcLlis = zeros(Float64, length(Llis), length(δtlis))
    for (i, L) in enumerate(Llis)
        τ = log(1 + √2)
        D = 8
        spatial_corr = load("exm/data/Bulk_measure/spatial_temporal_corr/L$(L)/τ$(τ)/D$(D)_ps$(sign)_0.jld", "spatial_corr")
        scLlis[i] = spatial_corr
        for (j, δt) in enumerate(δtlis)
            temporal_corr = load("exm/data/Bulk_measure/spatial_temporal_corr/L$(L)/τ$(τ)/D$(D)_ps$(sign)_$(δt).jld",  "temporal_corr")
            tcLlis[i, j] = temporal_corr
        end
    end
    save("exm/data/Bulk_measure/spatial_temporal_corr/stc_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", scLlis, "tcLlis", tcLlis)
end

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

function compute_post_selection(L::Int64, τ::Float64, D::Int64=35L, start_point::Int64=24, sign::Int64=0)
    pbc = true
    sample = (sign == 1) ? ones(Int, D, length(2:2:L)) : zeros(Int, D, length(2:2:L))

    initial_state = zeros(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
    initial_state[1] = 1.0 # initial state is all zero state

    statelis = generate_state(τ, initial_state, sample, temp= true, anyon_type=anyon_type)
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
            spatial_corr, temporal_corr = ref_correlation(L, ref2st, spatio=spatial, temporal=temporal)
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

function plot_stc_tlis(L::Int64=10, D::Int64=10, τ::Float64=log(1+√2); anyon_type::Symbol=:Fibo, sign::Int=0)
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

function plot_tc(L::Int, D::Int=10, τ::Float64=log(1+√2); anyon_type::Symbol=:Fibo)
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

function alpha_compute_corr(L, τ)
    D = get_system_params_corr(τ)

    temporal_corr_lis, spatial_corr = load("exm/data/Bulk_measure/temporal_corr/L$(L)/τ$(τ)/D$(D)_ps1.jld",  "temporal_corr_lis", "spatial_corr")

    inds = findall(x-> isapprox(x, 1.0, atol=0.1), temporal_corr_lis ./ spatial_corr)
    if !isempty(inds)
        Δt = (collect(1:length(temporal_corr_lis))./L)[inds][end-1]
        α = log(1+√2)/π/Δt
        return α
    else
        return NaN
    end
end

function alphalis_corr(γlis, L_list)
    αlis = [alpha_compute_corr(L, τ) for (L, τ) in zip(L_list[end], τlis)]
    return αlis
end

# 添加辅助函数


γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 


if length(ARGS) == 0
    println("No arguments provided.")
else
    L = parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    println("Received argument: $L, $inds")
    τ = τlis[inds]
    D, _, _ = get_system_params(τ, L)
    
    if length(ARGS) >= 3
        δt = parse(Int, ARGS[3])
        spatial_temporal_corr_varying(L, τ, D, δt, sign=1)
    else
        compute_post_selection(L, τ, D)
    end
end

# 示例运行代码
for i in 8:2:12
    τ = τlis[findfirst(γlis .== 1.0)]
    D, _, _ = get_system_params(τ, i)
    compute_post_selection(i, τ, D, round(Int, 24/35*div(D,i)))
end

# 设置参数并生成图像
gamma = 1.0
τ = 1000.0
L_list = collect([8, 10, 12])
fig = plot_corr(L_list)

# savefig(fig, "exm/data/Bulk_measure/corr_plot_L$(L_list[1])$((L_list[end]))_τ$(τ).pdf")