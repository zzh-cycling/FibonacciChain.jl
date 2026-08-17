using Distributed

const BULK_MEASURE_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BULK_MEASURE_CONFIG)
# using ClusterManagers

# const PROJECT_DIR = something(dirname(Base.active_project()), pwd())
# const NWORKERS = parse(Int, get(ENV, "SLURM_NTASKS", "512"))
# const CPUS_PER_TASK = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1"))
# addprocs(SlurmManager(NWORKERS), exeflags="--project=$(PROJECT_DIR) --threads=1")

@everywhere begin
    # using Pkg
    # Pkg.activate($PROJECT_DIR; io=devnull)
    using FibonacciChain
    using JLD
    using Statistics
    using Random
    using LaTeXStrings
    using Plots
    using LsqFit
    using LinearAlgebra
    using Measurements


function organize(args::Tuple)
    L, τ = args
    δtlis = get_δtL_Born(τ, L)
    tcLlis = zeros(Float64, length(δtlis))
    tcstderrlis = zeros(Float64, length(δtlis))
    sample_numlis = zeros(Int, length(δtlis)+1)
    average_spatial_corr, spatial_corr_stderr, sample_num = load(
        "exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt0_collect.jld",
        "average_spatial_corr",
        "spatial_corr_stderr",
        "samples_num",
    )
    sample_numlis[end] = sample_num
    for (j, δt) in enumerate(δtlis)
        average_temporal_corr, temporal_corr_stderr, sample_num = load(
            "exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)_collect.jld",
            "average_temporal_corr",
            "temporal_corr_stderr",
            "samples_num",
        )
        tcLlis[j] = average_temporal_corr
        tcstderrlis[j] = temporal_corr_stderr
        sample_numlis[j] = sample_num
    end

    save(
        "exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/stc_L$(L)_τ$(τ).jld",
        "average_spatial_corr",
        average_spatial_corr,
        "spatial_corr_stderr",
        spatial_corr_stderr,
        "δtlis",
        δtlis,
        "tcLlis",
        tcLlis,
        "tcstderrlis",
        tcstderrlis,
        "sample_numlis",
        sample_numlis,
    )
end

function get_correlation_dynamics_D(τ, L)
    cfg = Dict(atanh(0.4) => 20L, atanh(0.5) => 15L, atanh(0.6) => 10L)
    t = get(cfg, τ, 0)
    return t
end

function compute_ratio(L::Int64, τ_idx::Int64, index::Int64, δt::Int64 = 2)
    τ = τlis[τ_idx]
    D, _, _ = get_cfg_params_Born(τ_idx, L)
    model = fib_model(L)
    t = div(D, 2)
    sample = load(
        "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τ_idx)/t$(div(t,L))_samples$(index).jld",
        "sample",
    )
    # Here sample is of size D x (L/2), representing measurement outcomes at each layer for the monitored dynamics
    initial_state = zeros(length(anyon_basis(model)))
    initial_state[1] = 1.0 # initial state is all zero state

    t = div(D, 2)
    D1 = D + get_correlation_dynamics_D(τ, L)

    rng = MersenneTwister(index)
    config = MeasureConfig(τ = τ, mode = :sample, t₂ = t, enable_τ_eff = false)
    mo = bulk_evolution(model, initial_state, config, BitMatrix(sample))
    Flis = mo.free_energys

    t1 = t + get_correlation_dynamics_D(τ, L) # total evolution time after adding two ref qubits
    ref_sample = BitMatrix(zeros(Bool, 2*(t+δt+t1), length(2:2:L)))
    view(ref_sample, 1:D, :) .= Bool.(view(sample, :, :))

    if δt == 0
        ref_config =
            MeasureConfig(τ = τ, mode = :Born, rng = rng, x₂ = L÷2+1, t₂ = t, t₁ = t)
        spatial = true
        temporal = false
    else
        ref_config = MeasureConfig(
            τ = τ,
            mode = :Born,
            rng = rng,
            x₂ = L÷2+1,
            t₂ = t+δt,
            t₁ = t,
            x₁ = L÷2+1,
        )
        temporal = true
        spatial = false
    end

    ref_mo = reference_evolution(model, mo.state, ref_config, ref_sample)
    sample_layer, sample_free_energy = ref_mo.samples, ref_mo.free_energys # to compute temporal correlation, add ref qubit at site L/2+1
    view(sample_free_energy, 1:D) .= view(Flis, :)
    view(sample_layer, 1:D, :) .= view(sample, :, :)

    spatial_corr, temporal_corr =
        ref_correlation(model, ref_mo.state, spatial = spatial, temporal = temporal)
    sysrdm = reference_rdm(model, collect(1:div(L, 2)), ref_mo.state, traceref = false)
    S = ee(sysrdm)


    save(
        "exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D1,L))_Samples$(index).jld",
        "temporal_corr",
        temporal_corr,
        "spatial_corr",
        spatial_corr,
        "S",
        S,
        "sample_layer",
        sample_layer,
        "sample_free_energy",
        sample_free_energy,
    )
    return temporal_corr, spatial_corr, S, sample_layer, sample_free_energy
end

function spatial_temporal_corr_varyingt(
    L::Int64,
    τ::Float64,
    index::Int64,
    D::Int64 = 5L,
    δt::Int = 1,
)
    # | ----> |____| ----> |
    # 0       D   D+δt   D+δt+t  
    # compute how the spatial and temporal correlation changes with t, the evolution time after add two ref qubits. δt is the time interval between two ref qubits
    model = fib_model(L)

    # 1). First evolve to steady state with D time steps
    sample = load(
        "exm/data/Bulk_measure/Samples_monitored_dynamics/L$L/τ$(τ)/D$(div(D,L))_Samples$(index).jld",
        "sample",
    )
    initial_state = zeros(length(anyon_basis(model)))
    initial_state[1] = 1.0 # initial state is all zero state

    t = div(D, 2)
    D1 = D + get_correlation_dynamics_D(τ, L)
    pre_config = MeasureConfig(τ = τ, mode = :sample, t₂ = t, enable_τ_eff = false)
    pre_mo = bulk_evolution(model, initial_state, pre_config, BitMatrix(sample))
    Flis = pre_mo.free_energys

    rng = MersenneTwister(index)
    # tlis is the time list after adding two ref qubits.
    t1 = (τ ∈ τlis[[3, 4, 5, 6]]) ? t + 10L : t  # to adjust for longer evolution time for certain τ
    tlis = collect(1:2t1)
    spatial_corr_lis = zeros(Float64, length(tlis))
    temporal_corr_lis = zeros(Float64, length(tlis))
    eelis = zeros(Float64, length(tlis))

    # 2). Then add two ref qubits at different time slices and evolve for δt time, and to final δt + t time.


    for (idx, Δt) in enumerate(tlis)
        @show "calculation time t:" t
        ref_sample = BitMatrix(zeros(Int, 2*(t + δt + Δt), length(2:2:L)))
        if δt == 0
            ref_config =
                MeasureConfig(τ = τ, mode = :Born, rng = rng, x₂ = L÷2+1, t₂ = t, t₁ = t)
            spatial = true
            temporal = false
        else
            ref_config = MeasureConfig(
                τ = τ,
                mode = :Born,
                rng = rng,
                x₂ = L÷2+1,
                t₂ = t+δt,
                t₁ = t,
                x₁ = L÷2+1,
            )
            temporal = true
            spatial = false
        end

        ref_mo = reference_evolution(model, pre_mo.state, ref_config, ref_sample)
        sample_layer, sample_free_energy = ref_mo.samples, ref_mo.free_energys

        view(sample_free_energy, 1:D) .= view(Flis, :)
        view(sample_layer, 1:D, :) .= view(sample, :, :)
        sysrdm = reference_rdm(L, collect(1:div(L, 2)), ref_mo.state, traceref = false)
        eelis[idx] = ee(sysrdm)
        spatial_corr, temporal_corr =
            ref_correlation(L, ref_mo.state, spatial = spatial, temporal = temporal)
        temporal_corr_lis[idx] = temporal_corr
        spatial_corr_lis[idx] = spatial_corr
    end

    save(
        "exm/data/Bulk_measure/spatial_temporal_corr_varying_Born/L$(L)/τ$(τ)/D$(div(D1,L))_$(δt).jld",
        "temporal_corr_lis",
        temporal_corr_lis,
        "spatial_corr_lis",
        spatial_corr_lis,
        "eelis",
        eelis,
    )
    return temporal_corr_lis, spatial_corr_lis, eelis

end

function corr_collect(arg::Tuple)
    L, τ_idx, δt = arg
    τ = τlis[τ_idx]
    D = get_cfg_params_Born(τ_idx, L)[1]
    t = div(D, 2) # true circuits depth
    D1 = D + get_correlation_dynamics_D(τ, L)

    dir_path = "exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)"
    existing_files = filter(
            f -> startswith(f, "D$(div(D1,L))_Samples") && endswith(f, "jld"),
            readdir(dir_path),
        )
    samples_num = length(existing_files)
    println("Sample number: ", samples_num)

    success=0
    temporal_corr_ensemble = zeros(samples_num)
    spatial_corr_ensemble = zeros(samples_num)
    S_ensemble = zeros(samples_num)
    sample_free_energy_ensemble = zeros(samples_num, 2*(D1+δt))

    for (i, fname) in enumerate(existing_files)
        try
            temporal_corr, spatial_corr, S, sample_free_energy = load(joinpath(dir_path, fname),
                "temporal_corr",
                "spatial_corr",
                "S",
                "sample_free_energy",
            )
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
    average_free_energy_tlis = mean(sample_free_energy_ensemble, dims = 1)[:]
    stderr_free_energy_tlis =
        (std(sample_free_energy_ensemble, dims = 1) ./ sqrt(samples_num))[:]

    if success == samples_num
        save(
            "exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)_collect.jld",
            "average_temporal_corr",
            average_temporal_corr,
            "temporal_corr_stderr",
            temporal_corr_stderr,
            "average_spatial_corr",
            average_spatial_corr,
            "spatial_corr_stderr",
            spatial_corr_stderr,
            "average_EE",
            average_EE,
            "stderr_EE",
            stderr_EE,
            "average_free_energy_tlis",
            average_free_energy_tlis,
            "stderr_free_energy_tlis",
            stderr_free_energy_tlis,
            "samples_num",
            samples_num,
        )
        println("Completed L=$(L), τ=$(τ), δt=$(δt)")
    else
        println("No successful samples loaded for L=$(L), τ=$(τ), δt=$(δt). Skipping save.")
    end
end

function alpha_with_error_wt(τ; L = 16)
    # ---- read data ----
    file = "exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/stc_L$(L)_τ$(τ).jld"
    sc_μ = load(file, "average_spatial_corr")
    sc_σ = load(file, "spatial_corr_stderr")
    tc_μ = load(file, "tcLlis")
    tc_σ = load(file, "tcstderrlis")
    δtlis = load(file, "δtlis")

    ratio_μ = tc_μ ./ sc_μ
    ratio_σ = @. sqrt((tc_σ/tc_μ)^2 + (sc_σ/sc_μ)^2) * ratio_μ   # error propagation formula
    wt = @. 1 / ratio_σ^2                                         # weights

    # ---- weighted fit y = a*x + b ----
    model(x, p) = @. p[1] .* exp.(-x ./ p[2]) + p[3]
    fit = curve_fit(model, δtlis ./ L, ratio_μ, wt, [1.0, 1.0, 1.0])
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


seed_interval_lis = collect(1:100:2000)

function process_task(task)
    L, inds, index, δt = task
    try
        compute_ratio(L, inds, index, δt)
        return (L, inds, index, :success, nothing)
    catch e
        return (L, inds, index, :failed, e)
    end
end

end

if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia -p N corr_calculate.jl L τ_idx δt index_start index_end")
    println("Example: julia -p 16 corr_calculate.jl 10 7 2 1 1000")
else
    L = parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    δt = parse(Int, ARGS[3])
    index_start = parse(Int64, ARGS[4])
    index_end = parse(Int64, ARGS[5])

    println("=== Parallel Correlation Calculation ===")
    println("L = $L, τ_idx = $inds, δt = $δt")
    println("Sample index range: $index_start - $index_end")
    println("Total tasks: $(length(index_start:index_end))")
    println("Number of workers: $(nworkers())")

    taskslis = [(L, inds, i, δt) for i in index_start:index_end]

    println("\nStarting parallel processing...")
    results = pmap(process_task, taskslis; batch_size = 100)

    failed_tasks = [
        (L_res, inds, idx_res, error) for
        (L_res, inds, idx_res, status, error) in results if
        status != :success
    ]

    success_count = count(r -> r[4] == :success, results)
    failed_count = length(failed_tasks)

    println("\n=== Processing Complete ===")
    println("Total tasks: $(length(taskslis))")
    println("Successes: $success_count")
    println("Failures: $failed_count")

    if failed_count > 0
        println("\n=== Failed Task Details ===")
        for (i, (L_f, inds, idx_f, err)) in enumerate(failed_tasks)
            println("Failed $i: L=$L_f, τ_idx=$inds, index=$idx_f")
            println("  Error: $err")
        end

        failed_file = "failed_tasks_L$(L)_τidx$(inds)_δt$(δt)_batch$(index_start)_$(index_end).txt"
        open(failed_file, "w") do io
            println(io, "# Failed Task List")
            println(io, "# Format: L τ_idx δt sample_index")
            for (L_f, inds, idx_f, err) in failed_tasks
                println(io, "$L_f $inds $δt $idx_f  # Error: $err")
            end
        end
        println("\nFailed tasks saved to: $failed_file")
    end
end
