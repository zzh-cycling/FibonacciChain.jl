using FibonacciChain
using JLD
using Statistics
using Random
using LaTeXStrings
using Plots
using LsqFit
using LinearAlgebra
using Measurements

function get_δtL_Born(τ, L)
    if L == 6
        table = Dict(
                atanh(0.1)  => vcat(collect(1:150), collect(152:2:620)),
                atanh(0.2)  => (collect(50:2:120)),
                atanh(0.3)  => (collect(1:40)),
                atanh(0.4)  => (collect(1:35)),
                atanh(0.5)  => (collect(1:25)),
                atanh(0.6)  => (collect(1:10)),)
        δtlis = get(table, τ, collect(1:10))
    elseif L == 8
        table = Dict(
                atanh(0.1)  => (collect(50:50:500)),
                atanh(0.2)  => (collect(65:2:145)),
                atanh(0.3)  => (collect(1:54)),
                atanh(0.4)  => (collect(1:45)),
                atanh(0.5)  => (collect(1:30)),
                atanh(0.6)  => (collect(1:16)),
                atanh(1/√2)  => (collect(1:12)),)
        δtlis = get(table, τ, collect(1:10))
    elseif L == 10
        table = Dict(
                atanh(0.1)  => (collect(100:50:500)),
                atanh(0.2)  => (collect(80:2:120)),
                atanh(0.3)  => (collect(30:60)),
                atanh(0.4)  => (collect(1:25)),
                atanh(0.5)  => (collect(1:16)),
                atanh(0.6)  => (collect(1:22)),
                atanh(1/√2)  => (collect(1:16)),)
        δtlis = get(table, τ, collect(1:10))
    elseif L == 12
        table = Dict(
                atanh(0.1)  => (collect(300:50:600)),
                atanh(0.2)  => (collect(100:2:126)),
                atanh(0.3)  => (collect(40:60)),
                atanh(0.4)  => (collect(15:30)),
                atanh(0.5)  => (collect(1:15)),
                atanh(0.6)  => (collect(1:25)),
                atanh(1/√2)  => (collect(1:18)),)
        δtlis = get(table, τ, collect(1:15))
    elseif L == 14
        table = Dict(
                atanh(0.1)  => sort(vcat([475,500, 525], collect(50:100:550))),
                atanh(0.2)  => (collect(120:2:140)),
                atanh(0.3)  => (collect(45:80)),
                atanh(0.4)  => (collect(25:39)),
                atanh(0.5)  => (collect(5:20)),
                atanh(0.6)  => (collect(1:28)),
                atanh(1/√2)  => (collect(1:8)),
                atanh(0.8)  => (collect(1:15)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 16
        table = Dict(
                atanh(0.1)  => vcat(collect(80:100:680), [740]),
                atanh(0.2)  => (collect(135:2:149)),
                atanh(0.3)  => (collect(55:65)),
                atanh(0.4)  => (collect(32:42)),
                atanh(0.5)  => (collect(4:22)),
                atanh(0.6)  => (collect(1:16)),
                atanh(1/√2)  => (collect(1:8)),
                atanh(0.8)  => (collect(1:16)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 18
        table = Dict(
                atanh(0.1)  => sort(vcat([630, 640, 650, 660, 680, 710 ,720], collect(600:50:800))),
                atanh(0.2)  => sort(vcat([130, 140, 158, 189, 190], collect(150:5:180))),
                atanh(0.3)  => sort(vcat(collect(71:74),collect(60:5:90))),
                atanh(0.4)  => sort(vcat([38,39,41,42], collect(25:5:50))),
                atanh(0.5)  => (collect(20:28)),
                atanh(0.6)  => (collect(5:15)),
                atanh(1/√2)  => (collect(1:10)),
                atanh(0.8)  => (collect(1:10)),
                atanh(0.9)  => (collect(1:10)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 20
        table = Dict(
                atanh(0.2)  => (collect(173:2:179)),
                atanh(0.3)  => (collect(72:4:80)),
                atanh(0.4)  => (collect(40:48)),
                atanh(0.5)  => (collect(26:34)),
                atanh(0.6)  => (collect(12:20)),
                atanh(1/√2)  => (collect(1:12)),)
        δtlis = get(table, τ, collect(1:10))
    else
        δtlis = collect(1:10)
    end
    return  δtlis
end

function organize(args::Tuple)
    L, τ = args
    δtlis = get_δtL_Born(τ, L)
    tcLlis = zeros(Float64, length(δtlis))
    tcstderrlis = zeros(Float64, length(δtlis))
    sample_numlis = zeros(Int, length(δtlis)+1)
    average_spatial_corr, spatial_corr_stderr, sample_num = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt0_collect.jld", "average_spatial_corr", "spatial_corr_stderr", "samples_num")
    sample_numlis[end] = sample_num
    for (j, δt) in enumerate(δtlis)
        average_temporal_corr, temporal_corr_stderr, sample_num = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)_collect.jld",  "average_temporal_corr", "temporal_corr_stderr", "samples_num")
        tcLlis[j] = average_temporal_corr
        tcstderrlis[j] = temporal_corr_stderr
        sample_numlis[j] = sample_num
    end

    save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/stc_L$(L)_τ$(τ).jld", "average_spatial_corr", average_spatial_corr, "spatial_corr_stderr", spatial_corr_stderr, "δtlis", δtlis, "tcLlis", tcLlis, "tcstderrlis", tcstderrlis, "sample_numlis", sample_numlis)
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
    D1 = (τ ∈ τlis[[3,4,5,6]]) ? D + 20L : D 
    ref_sample = zeros(Int, 2*(D+δt+D1), length(2:2:L))
    view(ref_sample, 1:2D, :) .= view(sample, :, :)

    if δt == 0
        ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D, verbose=false, rng = rng, mode=:Born) # to compute temporal correlation, add ref qubit at site L/2+1
        spatial = true
        temporal = false
        view(sample_free_energy, 1:2D) .= view(Flis, :)
        view(sample_layer, 1:2D, :) .= view(sample, :, :)
    else
        ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D+δt, x₁ = L÷2+1, verbose=false, rng = rng, mode=:Born) # to compute temporal correlation, add ref qubit at site L/2+1
        temporal = true
        spatial = false
        view(sample_free_energy, 1:2D) .= view(Flis, :)
        view(sample_layer, 1:2D, :) .= view(sample, :, :)
    end

    spatial_corr, temporal_corr = ref_correlation(L, ref2stlis[end], spatial = spatial, temporal = temporal)
    sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2stlis[end], traceref = false)
    S = ee(sysrdm)
    

    save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D1,L))_Samples$(index).jld", "temporal_corr", temporal_corr, "spatial_corr", spatial_corr, "S", S, "sample_layer", sample_layer, "sample_free_energy", sample_free_energy)
    return temporal_corr, spatial_corr, S, sample_layer, sample_free_energy
end

function spatial_temporal_corr_varyingt(L::Int64, τ::Float64, index::Int64, D::Int64=5L, δt::Int=1)
    # | ----> |____| ----> |
    # 0       D   D+δt   D+δt+t  
    # compute how the spatial and temporal correlation changes with t, the evolution time after add two ref qubits. δt is the time interval between two ref qubits
    pbc = true
 
    # 1). First evolve to steady state with D time steps
    sample = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$L/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "sample")
    initial_state = zeros(length(anyon_basis(L, pbc)))
    initial_state[1] = 1.0 # initial state is all zero state

    statelis, Flis = generate_state(τ, initial_state, sample, enable_τ_eff=false)
    D = div(D, 2)
    
    rng = MersenneTwister(index)
    # tlis is the time list after adding two ref qubits.
    D1 = (τ ∈ τlis[[3,4,5,6]]) ? D + 20L : D  # to adjust for longer evolution time for certain τ
    tlis = collect(1:D1)
    spatial_corr_lis = zeros(Float64, length(tlis))
    temporal_corr_lis = zeros(Float64, length(tlis))
    eelis = zeros(Float64, length(tlis))
    
    # 2). Then add two ref qubits at different time slices and evolve for δt time, and to final δt + t time.
    
    
    for (idx, t) in enumerate(tlis)
        @show "calculation time t:" t
        ref_sample = zeros(Int, 2*(D + δt + t), length(2:2:L)) 
        if δt == 0
            ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D, verbose=false, rng = rng, mode=:Born) # to compute temporal correlation, add ref qubit at site L/2+1
            spatial = true
            temporal = false
            view(sample_free_energy, 1:2D) .= view(Flis, :)
            view(sample_layer, 1:2D, :) .= view(sample, :, :)
        else
            ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D+δt, x₁ = L÷2+1, verbose=false, rng = rng, mode=:Born) # to compute temporal correlation, add ref qubit at site L/2+1
            temporal = true
            spatial = false
            view(sample_free_energy, 1:2D) .= view(Flis, :)
            view(sample_layer, 1:2D, :) .= view(sample, :, :)
        end
        sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2stlis[end], traceref = false)
        eelis[idx] = ee(sysrdm)
        spatial_corr, temporal_corr = ref_correlation(L, ref2stlis[end], spatial=spatial, temporal=temporal)
        temporal_corr_lis[idx] = temporal_corr
        spatial_corr_lis[idx] = spatial_corr
    end

    save("exm/data/Bulk_measure/spatial_temporal_corr_varying_Born/L$(L)/τ$(τ)/D$(div(D1,L))_$(δt).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr_lis", spatial_corr_lis, "eelis", eelis)
    return temporal_corr_lis, spatial_corr_lis, eelis
   
end

function corr_collect(arg::Tuple)
    L, τ, δt = arg
    D = get_system_params(τ, L)[1]
    D = div(D, 2) # true circuits depth
    D1 = D + get_correlation_dynamics_D(τ, L)
    samples_num = 1000
    println("Sample number: ", samples_num)
    
    success=0
        temporal_corr_ensemble = zeros(samples_num)
        spatial_corr_ensemble = zeros(samples_num)
        S_ensemble = zeros(samples_num)
        sample_free_energy_ensemble = zeros(samples_num, 2*(D1+δt+D))

        for i in 1:samples_num
            # @show i
            try
                temporal_corr, spatial_corr, S, sample_free_energy = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D1,L))_Samples$(i).jld",  "temporal_corr", "spatial_corr", "S", "sample_free_energy")
                temporal_corr_ensemble[i] = temporal_corr
                spatial_corr_ensemble[i] = spatial_corr
                sample_free_energy_ensemble[i, :] = sample_free_energy
                S_ensemble[i] = S
                success += 1
            catch e
                println("Error loading sample $(i) for L=$(L), τ=$(τ), δt=$(δt): ", e)
            end
        end
    
        average_temporal_corr = mean(temporal_corr_ensemble)
        average_spatial_corr = mean(spatial_corr_ensemble)
        temporal_corr_stderr = std(temporal_corr_ensemble) ./ sqrt(samples_num)
        spatial_corr_stderr = std(spatial_corr_ensemble) / sqrt(samples_num)
        average_EE = mean(S_ensemble)
        stderr_EE = std(S_ensemble) / sqrt(samples_num)
        average_free_energy_tlis = mean(sample_free_energy_ensemble, dims=1)[:]
        stderr_free_energy_tlis = (std(sample_free_energy_ensemble, dims=1) ./ sqrt(samples_num))[:]
    
        if success == samples_num
            save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)_collect.jld", "average_temporal_corr", average_temporal_corr, 
        "temporal_corr_stderr", temporal_corr_stderr, 
        "average_spatial_corr", average_spatial_corr, 
        "spatial_corr_stderr", spatial_corr_stderr, 
        "average_EE", average_EE,
        "stderr_EE", stderr_EE,
        "average_free_energy_tlis", average_free_energy_tlis,
        "stderr_free_energy_tlis", stderr_free_energy_tlis, 
        "samples_num", samples_num)
        println("Completed L=$(L), τ=$(τ), δt=$(δt)")
        else
            println("No successful samples loaded for L=$(L), τ=$(τ), δt=$(δt). Skipping save.")
        end
end

function alpha_with_error_wt(τ; L=16)
    # ---- read data ----
    file = "exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/stc_L$(L)_τ$(τ).jld"
    sc_μ = load(file, "average_spatial_corr")
    sc_σ = load(file, "spatial_corr_stderr")
    tc_μ = load(file, "tcLlis")
    tc_σ = load(file, "tcstderrlis")
    δtlis= load(file, "δtlis")

    ratio_μ = tc_μ ./ sc_μ
    ratio_σ = @. sqrt( (tc_σ/tc_μ)^2 + (sc_σ/sc_μ)^2 ) * ratio_μ   # error propagation formula
    wt = @. 1 / ratio_σ^2                                         # weights

    # ---- weighted fit y = a*x + b ----
    model(x, p) = @. p[1] .* exp.(-x ./p[2])  + p[3]
    fit = curve_fit(model, δtlis./L, ratio_μ, wt, [1.0, 1.0, 1.0])
    a, b, c = fit.param                                
    Σ = estimate_covar(fit)                            # 2×2 covariance matrix
    σa, σb, σc = sqrt.(diag(Σ))                            # a, b stderr
    a_m = a ± σa
    b_m = b ± σb
    c_m = c ± σc

    t = -b_m * log((1-c_m)/a_m)  # exponential fit
    # t = (1 - b_m) / a_m # first linear fit
    α = log(1 + √2) / π / t

    return α, t
end

function alphalis_corr(τlis)
    αlis = [alpha_with_error_wt(τ, L = 16)[1].val for (τ) in τlis]
    α_stderrlis = [alpha_with_error_wt(τ, L = 16)[1].err for (τ) in τlis]
    tlis = [alpha_with_error_wt(τ, L = 16)[2] for (τ) in τlis]
    @show tlis
    return αlis, α_stderrlis
end


function plot_tc(inds::Int64)
    # Llis = collect(6:2:14)
    Llis = collect(6:2:14)
    τ = τlis[inds]

    plt = plot(
        label=false,
        legend_background_color=nothing,
        legend_foreground_color=nothing,
        xlabel=L"δt /L",
        ylabel=L"g(0, \Delta t)/g_{space}",
        title=latexstring("γ= $(round(tanh(τ), digits=3))"),
        # ylim = (0.8, 1.2), 
    )
    
    c = cgrad(:blues, length(Llis), categorical=true)

    test_lis = [2.4, 1.1917091921566003, 0.783, 0.521, 0.45646685808273624, 0.3543349243759066, 0.282, 0.234, 0.196, 0.154, 0.12, 0.1]

    for (idx, L) in enumerate(Llis[1:end-1])
        average_spatial_corr, spatial_corr_stderr, δtlis, tcLlis, tcstderrlis = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/stc_L$(L)_τ$(τ).jld", "average_spatial_corr", "spatial_corr_stderr", "δtlis", "tcLlis", "tcstderrlis")


        scatter!(plt, δtlis ./ L, tcLlis./ average_spatial_corr, yerror = tcstderrlis ./ average_spatial_corr, label="L=$(L)", lw=2, marker=:circle, ms=6, c = c[idx])
    end
    scatter!(plt, [test_lis[inds]], [1.0], ms=8, label=false, mc=:red, m=:star5)

    average_spatial_corr, spatial_corr_stderr, δtlis, tcLlis, tcstderrlis = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(Llis[end])/stc_L$(Llis[end])_τ$(τ).jld", "average_spatial_corr", "spatial_corr_stderr", "δtlis", "tcLlis", "tcstderrlis")
    
    
    plot!(plt, δtlis ./ Llis[end], tcLlis./ average_spatial_corr, yerror = tcstderrlis ./ average_spatial_corr, label="L=$(Llis[end])", lw=2, marker=:circle, ms=6, c = c[end])
    
    model(x, p) = @. p[1] .* exp.(-x ./p[2])  + p[3]
    fit = curve_fit(model, δtlis./ Llis[end], tcLlis./ average_spatial_corr, [1.0, 1.0, 1.0])
    a,b,c = fit.param
    plot!(plt, δtlis ./ Llis[end], model(δtlis ./ Llis[end], fit.param), lw=2, ls=:dash, c=:black, label="$(round(a, digits=2))exp(-x/$(round(b, digits=2))) + $(round(c, digits=2))")
    return plt
end

# αlis = [0.11689580257066254, 0.23541811040484406, 0.3583013105614177, 0.5384835435116891, 0.6146118194603724, 0.7917648159117402, 0.9948578942184046, 1.1989313084170514, 1.431377174334643, 1.8217527673350005, 2.337916051413251, 2.8054992616959007]

# α_lis = [0.22996304034610712, 0.6025887408957163, 0.7737905089041682, 1.2185251168172877, 2.5759898698081787]
# α_stderrlis = [0.03939320246529536, 0.03727317853604598, 0.03256489438182987, 0.03328157129222018, 0.07449016381263754] # by c_{ent}/c_{Casimir}

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
    seed = parse(Int, ARGS[4])
    # seedinds = parse(Int, ARGS[4])
    τ = τlis[inds]
    D, _, _ = get_system_params(τ, L)
    # seed = seed_interval_lis[seedinds]
    println("Computed spatial_temporal_corr_varyingt for L=$L, τ=$τ, D=$D, δt=$δt, seedlis=$(seed):$(seed+99)")
    compute_ratio(L, τ, seed, D, δt)
    # compute_total(L, τ, seed, D, δt)
    # corr_collect(L, τ)
end