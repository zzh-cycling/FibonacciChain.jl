using Statistics
using Distributed

const BULK_MEASURE_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BULK_MEASURE_CONFIG)

println("requested workers: ", nworkers())
println("total procs:       ", nprocs())
@everywhere println("host: ", gethostname(), "  pid: ", getpid())

@everywhere begin
using FibonacciChain
using JLD
using Random
using LinearAlgebra

function get_correlation_dynamics_D(τ, L)
    cfg = Dict(atanh(0.3) => 25L, atanh(0.4) => 20L, atanh(0.5) => 15L, atanh(0.6) => 10L)
    t = get(cfg, τ, 0)
    return t
end

function spatial_temporal_corr_varyingt(args::Tuple)
    L, τind, index, δt = args
    try
        τ = τlis[τind]
        periods = get_cfg_params_Born(τind, L)[1]
        layers = 2 * periods
        # | ----> |____| ----> |
        # 0       layers   layers+δt   layers+δt+t
        # compute how the spatial and temporal correlation changes with t, the evolution time after add two ref qubits. δt is the time interval between two ref qubits

        model = fib_model(L)

        # 1). First evolve to steady state for the configured number of periods
        sample = load(
            "exm/data/Bulk_measure/monitored_dynamics/L$L/gammaind$(τind)/t$(div(periods,L))_samples$(index).jld",
            "sample",
        )
        println("Loaded sample for L=$(L), τ=$(τ), index=$(index)")
        initial_state = zeros(length(anyon_basis(model)))
        initial_state[1] = 1.0 # initial state is all zero state

        t = periods
        pre_config = MeasureConfig(τ = τ, mode = :sample, t₂ = t, enable_τ_eff = false)
        pre_mo = bulk_evolution(model, initial_state, pre_config, sample)
        Flis = pre_mo.free_energys

        rng = MersenneTwister(index)
        # tlis is the time list after adding two ref qubits.
        t1 = t + get_correlation_dynamics_D(τ, L) # total evolution time after adding two ref qubits
        tlis = collect(1:t1)
        spatial_corr_lis = zeros(Float64, length(tlis))
        temporal_corr_lis = zeros(Float64, length(tlis))
        eelis = zeros(Float64, length(tlis))

        # 2). Then add two ref qubits at different time slices and evolve for δt time, and to final δt + t time.


        for (idx, Δt) in enumerate(tlis)
            println("calculation time t: $t")
            ref_sample = BitMatrix(zeros(Int, 2*(t + δt + Δt), length(2:2:L)))
            if δt == 0
                ref_config = MeasureConfig(
                    τ = τ,
                    t₂ = t,
                    t₁ = t,
                    rng = rng,
                    mode = :Born,
                    x₂ = L÷2+1,
                )
                ref_mo = reference_evolution(model, pre_mo.state, ref_config, ref_sample)
                sample_layer, sample_free_energy = ref_mo.samples, ref_mo.free_energys  # to compute temporal correlation, add ref
                spatial = true
                temporal = false
                view(sample_free_energy, 1:layers) .= view(Flis, :)
                view(sample_layer, 1:layers, :) .= view(sample, :, :)
            else
                ref_config = MeasureConfig(
                    τ = τ,
                    t₂ = t + δt,
                    t₁ = t,
                    rng = rng,
                    mode = :Born,
                    x₂ = L÷2+1,
                    x₁ = L÷2+1,
                )
                ref_mo = reference_evolution(model, pre_mo.state, ref_config, ref_sample)
                sample_layer, sample_free_energy = ref_mo.samples, ref_mo.free_energys  # to compute temporal correlation, add ref qubit at site L/2+1
                temporal = true
                spatial = false
                view(sample_free_energy, 1:layers) .= view(Flis, :)
                view(sample_layer, 1:layers, :) .= view(sample, :, :)
            end
            sysrdm =
                reference_rdm(model, collect(1:div(L, 2)), ref_mo.state, traceref = false)
            eelis[idx] = ee(sysrdm)
            spatial_corr, temporal_corr =
                ref_correlation(model, ref_mo.state, spatial = spatial, temporal = temporal)
            temporal_corr_lis[idx] = temporal_corr
            spatial_corr_lis[idx] = spatial_corr
        end

        save(
            "exm/data/Bulk_measure/spatial_temporal_corr_varying_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(2t1,L))_Samples$(index).jld",
            "temporal_corr_lis",
            temporal_corr_lis,
            "spatial_corr_lis",
            spatial_corr_lis,
            "eelis",
            eelis,
        )
        return index, :success, nothing
    catch e
        return index, :error, e
    end
end

function collect_spatial_temporal_corr_varying_Born_data(
    L::Int64, τind::Int64, δt::Int64)
    
    τ = τlis[τind]
    periods = get_cfg_params_Born(τind, L)[1]
    src_dir = joinpath("exm/data/Bulk_measure/spatial_temporal_corr_varying_Born/L$(L)/τ$(τlis[τind])/dt$(δt)")
    dst_file = joinpath("exm/data/Bulk_measure/spatial_temporal_corr_varying_Born/L$(L)/τ$(τlis[τind])/dt$(δt)_collect.jld")
    mkpath(dirname(dst_file))
    files = sort(filter(
        path -> isfile(path) && endswith(path, ".jld") && !endswith(path, "_collect.jld"),
        readdir(src_dir; join = true),
    ))

    if isempty(files)
        error("No .jld files found in $(src_dir)")
    end

    spatial_corrs = Vector{Vector{Float64}}()
    temporal_corrs = Vector{Vector{Float64}}()
    ees = Vector{Vector{Float64}}()

    for src_file in files
        data = load(src_file)
        push!(spatial_corrs, Float64.(vec(data["spatial_corr_lis"])))
        push!(temporal_corrs, Float64.(vec(data["temporal_corr_lis"])))
        push!(ees, Float64.(vec(data["eelis"])))
    end

    spatial_matrix = hcat(spatial_corrs...)
    temporal_matrix = hcat(temporal_corrs...)
    ee_matrix = hcat(ees...)

    nfiles = length(files)
    scale = sqrt(nfiles)

    average_spatial_corr = vec(mean(spatial_matrix; dims = 2))
    average_temporal_corr = vec(mean(temporal_matrix; dims = 2))
    average_EE = vec(mean(ee_matrix; dims = 2))

    spatial_corr_stderr = vec(std(spatial_matrix; dims = 2, corrected = false)) ./ scale
    temporal_corr_stderr = vec(std(temporal_matrix; dims = 2, corrected = false)) ./ scale
    stderr_EE = vec(std(ee_matrix; dims = 2, corrected = false)) ./ scale

    mkpath(dirname(dst_file))
    save(
        dst_file,
        "average_spatial_corr", average_spatial_corr,
        "average_temporal_corr", average_temporal_corr,
        "average_EE", average_EE,
        "spatial_corr_stderr", spatial_corr_stderr,
        "temporal_corr_stderr", temporal_corr_stderr,
        "stderr_EE", stderr_EE,
    )

    println("Collected $(nfiles) files from $(src_dir) to $(dst_file)")
    return dst_file, nfiles
end


end

# main parallel function
function compute_parallel_batch(
    L::Int64,
    τ_idx::Int64,
    seed_range::Vector{Int} = collect(1:10000),
    )
    τ = τlis[τ_idx]
    println("Starting parallel computation with $(nprocs()) processes")
    println("Parameters: L=$L, τ=$τ")
    println("Computing seeds: $(seed_range[1]) to $(seed_range[end])")
    dtlis = get_δtL_Born(τ, L)
    # create parameter tuples for each task
    tasks = [(L, τ_idx, index, δt) for index in seed_range for δt in dtlis]

    # using pmap run
    println("Submitting $(length(tasks)) tasks to worker processes...")
    results = pmap(spatial_temporal_corr_varyingt, tasks, batch_size = 1)

    # outcome summary
    success_count = 0
    error_count = 0

    for (index, status, result) in results
        if status == :success
            success_count += 1
            println("✓ Task $index completed successfully")
        else
            error_count += 1
            println("✗ Task $index failed: $result")
        end
    end

    println("\nBatch computation completed:")
    println("  Successful: $success_count")
    println("  Failed: $error_count")
    println("  Total: $(success_count + error_count)")

    return results
end


if length(ARGS) == 0
    println("Usage examples:")
else
    L = parse(Int64, ARGS[1])
    τ_idx = parse(Int64, ARGS[2])
    index_start = parse(Int, ARGS[3])
    index_end = parse(Int, ARGS[4])
    compute_parallel_batch(L, τ_idx, collect(index_start:index_end))
end
