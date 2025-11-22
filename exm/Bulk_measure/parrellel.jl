using Distributed

@everywhere using FibonacciChain
@everywhere using JLD
@everywhere using Random
@everywhere using LinearAlgebra

@everywhere function get_system_params(τ, L)
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

@everywhere function compute_ratio(L::Int64, τ::Float64, index::Int64, D::Int64=16L, δt::Int64=2)
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
        ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D, verbose=false, rng = rng, mode=:Born)
        spatial = true
        temporal = false
        view(sample_free_energy, 1:2D) .= view(Flis, :)
    else
        ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D+δt, x₁ = L÷2+1, verbose=false, rng = rng, mode=:Born)
        temporal = true
        spatial = false
        view(sample_free_energy, 1:2D) .= view(Flis, :)
    end
    
    spatial_corr, temporal_corr = ref_correlation(L, ref2stlis[end], spatial = spatial, temporal = temporal)
    sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2stlis[end], traceref = false)
    S = ee(sysrdm)
    
    save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D,L))_Samples$(index).jld", 
         "temporal_corr", temporal_corr, "spatial_corr", spatial_corr, "S", S, 
         "sample_layer", sample_layer, "sample_free_energy", sample_free_energy)
    
    return temporal_corr, spatial_corr, S, sample_layer, sample_free_energy
end

# wrapper function
@everywhere function compute_single_task(task_params)
    L, τ, index, D, δt = task_params
    try
        result = compute_ratio(L, τ, index, D, δt)
        return (index, :success, result)
    catch e
        return (index, :error, string(e))
    end
end

# main parallel function
function compute_parallel_batch(L::Int64, τ::Float64, seed_range::UnitRange{Int}, D::Int64, δt::Int64)
    println("Starting parallel computation with $(nprocs()) processes")
    println("Parameters: L=$L, τ=$τ, D=$D, δt=$δt")
    println("Computing seeds: $(first(seed_range)) to $(last(seed_range))")
    
    # 创建任务参数列表
    tasks = [(L, τ, index, D, δt) for index in seed_range]
    
    # 使用pmap并行执行
    println("Submitting $(length(tasks)) tasks to worker processes...")
    results = pmap(compute_single_task, tasks, batch_size=1)
    
    # 处理结果
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

# 批量处理多个δt值
function compute_parallel_multiple_dt(L::Int64, τ::Float64, seed_range::UnitRange{Int}, D::Int64, δt_list::Vector{Int})
    total_tasks = length(seed_range) * length(δt_list)
    println("Starting computation for multiple δt values")
    println("Total tasks: $total_tasks")
    
    all_results = Dict()
    
    for δt in δt_list
        println("\n" * "="^50)
        println("Processing δt = $δt")
        println("="^50)
        
        results = compute_parallel_batch(L, τ, seed_range, D, δt)
        all_results[δt] = results
    end
    
    return all_results
end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0
τlis[findfirst(γlis .== 1/√2)] = log(1 + √2)

if length(ARGS) == 0
    println("Usage examples:")
    println("julia -p 4 parallel.jl compute L τ_idx δt start_seed end_seed")
    println("julia -p 4 parallel.jl collect L τ_idx")
    println("julia -p 4 parallel.jl batch L τ_idx start_seed end_seed")
    println("\nAvailable τ indices (1-$(length(τlis))):")
    for (i, (γ, τ)) in enumerate(zip(γlis, τlis))
        println("  $i: γ=$γ, τ=$τ")
    end
else
    mode = parse(Int64, ARGS[1])
    L = parse(Int64, ARGS[2])
    τ_idx = parse(Int64, ARGS[3])
    τ = τlis[τ_idx]
    D, _, _ = get_system_params(τ, L)
    
    if mode == 1
        # single δt value computation
        δt = parse(Int, ARGS[4])
        start_seed = parse(Int, ARGS[5])
        end_seed = parse(Int, ARGS[6])
        
        println("Single δt computation mode")
        compute_parallel_batch(L, τ, start_seed:end_seed, D, δt)
        
    elseif mode == 2
        # batch computation for multiple δt values
        start_seed = parse(Int, ARGS[4])
        end_seed = parse(Int, ARGS[5])
        δt_list = (L == 16) ? collect(0:10) : collect(0:20)
        
        println("Batch computation mode")
        compute_parallel_multiple_dt(L, τ, start_seed:end_seed, D, δt_list)

    end
end