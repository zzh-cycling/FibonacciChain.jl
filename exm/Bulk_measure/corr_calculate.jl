using FibonacciChain
using JLD
using Statistics
using Random
using LaTeXStrings
using Plots

function get_δtL(τ)
    table = Dict(
            atanh(0.1)  => (collect(25:35)),
            atanh(0.2)  => (collect(10:20)),
            atanh(0.3)  => (collect(5:15)),
            atanh(0.4)  => (collect(2:10)),)
    δtlis = get(table, τ, collect(1:8))
    return  δtlis
end

function organize( L::Int=8, τ::Float64= log(1 + √2))
    δtlis = collect(1:20)
    tcLlis = zeros(Float64, length(δtlis))
    tcstderrlis = zeros(Float64, length(δtlis))
    
    average_spatial_corr, spatial_corr_stderr = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt0_collect.jld", "average_spatial_corr", "spatial_corr_stderr")
    for (j, δt) in enumerate(δtlis)
        average_temporal_corr, temporal_corr_stderr = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)_collect.jld",  "average_temporal_corr", "temporal_corr_stderr")
        tcLlis[j] = average_temporal_corr
        tcstderrlis[j] = temporal_corr_stderr
    end
    
    save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/stc.jld", "average_spatial_corr", average_spatial_corr, "spatial_corr_stderr", spatial_corr_stderr, "δtlis", δtlis, "tcLlis", tcLlis, "tcstderrlis", tcstderrlis)
end

function get_system_params(τ, L)
    cfg = Dict(
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
    D, step, start = get(cfg, τ, (5L, 2, 2L))
    inds = collect(1:step:D)
    avg_range = start:D-5
    return D, inds, avg_range
end

function compute_total(L::Int64, τ::Float64, index::Int64, D::Int64=35L, δt::Int64=2)
    for index in index:index+99
        @show "Computing index=$index"
        @time compute_ratio(L, τ, index, D, δt)
    end
end

function compute_ratio(L::Int64, τ::Float64, index::Int64, D::Int64=16L, δt::Int64=2)
    pbc = true
    sample = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$L/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "sample")

    initial_state = zeros(length(anyon_basis(L, pbc)))
    initial_state[1] = 1.0 # initial state is all zero state

    rng = MersenneTwister(index)

    statelis, Flis = generate_state(τ, initial_state, sample, enable_τ_eff=false)
    D = div(D, 2)
    ref_sample = zeros(Int, 2*(D+δt+D), length(2:2:L))
    view(ref_sample, 1:2D, :) .= view(sample, :, :)

    if δt == 0
        ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D, verbose=false, rng = rng) # to compute temporal correlation, add ref qubit at site L/2+1
        spatial = true
        temporal = false
        view(sample_free_energy, 1:2D) .= view(Flis, :)
    else
        ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D+δt, x₁ = L÷2+1, verbose=false, rng = rng) # to compute temporal correlation, add ref qubit at site L/2+1
        temporal = true
        spatial = false
        view(sample_free_energy, 1:2D) .= view(Flis, :)
    end
    spatial_corr, temporal_corr = ref_correlation(L, ref2stlis[end], spatial = spatial, temporal = temporal)
    sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2stlis[end], traceref = false)
    S = ee(sysrdm)
    

    save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D,L))_Samples$(index).jld", "temporal_corr", temporal_corr, "spatial_corr", spatial_corr, "S", S, "sample_layer", sample_layer, "sample_free_energy", sample_free_energy)
    return temporal_corr, spatial_corr, S, sample, sample_free_energy
end


function corr_collect(L::Int64, τ::Float64, D::Int64=35L)
    samples_num = 2000
    δtlis = collect(0:20)  # Adjust this range based on the δt values you have used
    
    for δt in δtlis
        temporal_corr_ensemble = zeros(samples_num)
        spatial_corr_ensemble = zeros(samples_num)
        S_ensemble = zeros(samples_num)
        sample_free_energy_ensemble = zeros(samples_num, 2*(D+δt))

        for i in 1:samples_num
            @show i
            temporal_corr, spatial_corr, S, sample_free_energy = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D,2L))_Samples$(i).jld",  "temporal_corr", "spatial_corr", "S", "sample_free_energy")
            temporal_corr_ensemble[i] = temporal_corr
            spatial_corr_ensemble[i] = spatial_corr
            sample_free_energy_ensemble[i, :] = sample_free_energy
            S_ensemble[i] = S
        end
    
        average_temporal_corr = mean(temporal_corr_ensemble)
        average_spatial_corr = mean(spatial_corr_ensemble)
        temporal_corr_stderr = std(temporal_corr_ensemble) ./ sqrt(samples_num)
        spatial_corr_stderr = std(spatial_corr_ensemble) / sqrt(samples_num)
        average_EE = mean(S_ensemble)
        stderr_EE = std(S_ensemble) / sqrt(samples_num)
        average_free_energy_tlis = mean(sample_free_energy_ensemble, dims=1)[:]
        stderr_free_energy_tlis = (std(sample_free_energy_ensemble, dims=1) ./ sqrt(samples_num))[:]
    
        save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)_collect.jld", "average_temporal_corr", average_temporal_corr, 
        "temporal_corr_stderr", temporal_corr_stderr, 
        "average_spatial_corr", average_spatial_corr, 
        "spatial_corr_stderr", spatial_corr_stderr, 
        "average_EE", average_EE,
        "stderr_EE", stderr_EE,
        "average_free_energy_tlis", average_free_energy_tlis,
        "stderr_free_energy_tlis", stderr_free_energy_tlis)
    end
end

function alpha_compute_corr(τ)
    # Llis = collect(8:4:20)
    L = 8
    δtlis= collect(1:20)

    average_spatial_corr, spatial_corr_stderr, δtlis, tcLlis, tcstderrlis = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/stc.jld", "average_spatial_corr", "spatial_corr_stderr", "δtlis", "tcLlis", "tcstderrlis")

    ratio = tcLlis./average_spatial_corr
    inds = findall(x-> isapprox(x, 1.0, atol=0.1), ratio)
    
    linear_model(x,p) = p[1] * x .+ p[2]

    if !isempty(inds)
        fit = curve_fit(linear_model, δtlis[inds]./L, ratio[inds], [1.0, 1.0])
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


function plot_tc(inds::Int64)
    Llis = collect(8:2:10)
    τ = τlis[inds]

    plt = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing,
        xlabel=L"δt /L",
        ylabel=L"g(0, \Delta t)/g_{space}",
        title=latexstring("γ= $(round(tanh(τ), digits=3))"),
        ylim = (0.8, 1.2), 
    )
    
    c = cgrad(:reds, length(Llis), categorical=true)

    for (idx, L) in enumerate(Llis)
        average_spatial_corr, spatial_corr_stderr, δtlis, tcLlis, tcstderrlis = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/stc.jld", "average_spatial_corr", "spatial_corr_stderr", "δtlis", "tcLlis", "tcstderrlis")
    
        test_lis = [0.0, 1.371163364681693, 0.0, 0.0, 0.48786165394603076, 0.0, 0.0, 0.22, 0.0, 0.0, 0.0, 0.0]
    
    
        plot!(plt, δtlis ./ L, tcLlis./ average_spatial_corr, yerror = tcstderrlis ./ average_spatial_corr, label="L=$(L)", lw=2, marker=:o, ms=6, c = c[idx])
        scatter!(plt, [test_lis[inds]], [1.0], ms=8, label=false, mc=:red, m=:star5)
    end
    return plt
end

# αlis = [0.20460722142668827, 0.5750604170268026, 1.2752269371345004]

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)
seed_interval_lis = collect(1:100:2000)

if length(ARGS) == 0
    println("No arguments provided.")
else
    L = parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    δt = parse(Int, ARGS[3])
    seedinds = parse(Int, ARGS[4])
    τ = τlis[inds]
    D, _, _ = get_system_params(τ, L)
    seed = seed_interval_lis[seedinds]
    println("Computed spatial_temporal_corr_varyingt for L=$L, τ=$τ, D=$D, δt=$δt, seedlis=$(seed):$(seed+99)")
    # compute_ratio(L, τ, seed, D, δt)
    compute_total(L, τ, seed, D, δt)
    # corr_collect(L, τ)
end