using Distributed

@everywhere using FibonacciChain
@everywhere using JLD
@everywhere using Random
@everywhere using LinearAlgebra

println("requested workers: ", nworkers())
println("total procs:       ", nprocs())
@everywhere println("host: ", gethostname(), "  pid: ", getpid())

@everywhere function get_δtL(τ, L)
    if L == 6
        table = Dict(
                atanh(0.1)  => (collect(4:150)),
                atanh(0.2)  => (collect(0:20)),
                atanh(0.3)  => (collect(0:20)),
                atanh(0.4)  => (collect(0:20)),)
        δtlis = get(table, τ, collect(0:10))
    elseif L == 8
        if τ ∈ τlis[1:4]
            δtlis = collect(0:20)
        else
            δtlis = collect(0:10)
        end
    elseif L == 10
        table = Dict(
                atanh(0.1)  => vcat([0], collect(20:40)),
                atanh(0.2)  => vcat([0],collect(10:20)),
                atanh(0.3)  => (collect(0:20)),
                atanh(0.4)  => (collect(0:20)),)
        δtlis = get(table, τ, collect(0:10))
    elseif L == 12
        table = Dict(
                atanh(0.1)  => vcat([0], collect(80:100)),
                atanh(0.2)  => vcat([0], collect(30:40)),
                atanh(0.3)  => vcat([0], collect(15:25)),
                atanh(0.4)  => vcat([0], collect(15:25)),)
        δtlis = get(table, τ, collect(0:10))
    elseif L == 14
        table = Dict(
                atanh(0.1)  => vcat([0],collect(28:35)),
                atanh(0.2)  => vcat([0],collect(14:20)),
                atanh(0.3)  => vcat([0],collect(7:14)),
                atanh(0.4)  => vcat([0],collect(5:10)),
                atanh(0.5)  => vcat([0],collect(5:10)),)
        δtlis = get(table, τ, collect(0:8))
    elseif L == 16
        table = Dict(
                atanh(0.1)  => vcat([0],collect(32:42)),
                atanh(0.2)  => vcat([0],collect(16:24)),
                atanh(0.3)  => vcat([0],collect(8:16)),
                atanh(0.4)  => vcat([0],collect(4:12)),
                atanh(0.5)  => vcat([0],collect(4:12)),)
        δtlis = get(table, τ, collect(0:8))
    elseif L == 18
        table = Dict(
                atanh(0.1)  => vcat([0],collect(140:150)),
                atanh(0.2)  => vcat([0],collect(55:70)),
                atanh(0.3)  => vcat([0],collect(30:40)),
                atanh(0.4)  => vcat([0],collect(25:35)),
                atanh(0.5)  => vcat([0],collect(18:25)),
                atanh(0.6)  => vcat([0],collect(15:22)),
                atanh(1/√2)  => vcat([0],collect(10:18)),
                atanh(0.8)  => vcat([0],collect(10:18)),)
        δtlis = get(table, τ, collect(0:8))
    else
        δtlis = collect(1:10)
    end
    return  δtlis
end

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

function parse_filename(filename)
    # 解析 "missing_L12_inds1.txt" 格式
    # 使用正则表达式提取 L 和 inds
    m = match(r"missing_L(\d+)_inds(\d+)\.txt", filename)
    if m === nothing
        error("Cannot parse filename: $filename")
    end
    
    L = parse(Int, m.captures[1])
    inds = parse(Int, m.captures[2])
    return L, inds
end

# 新增：读取缺失任务文件
function read_missing_tasks(filename)
    if !isfile(filename)
        error("File not found: $filename")
    end
    
    tasks = Tuple{Int, Int}[]  # (δt, index)
    
    open(filename, "r") do file
        for line in eachline(file)
            line = strip(line)
            if isempty(line)
                continue
            end
            
            # 解析 "dt81 Sample1767" 格式
            parts = split(line)
            if length(parts) != 2
                @warn "Skipping invalid line: $line"
                continue
            end
            
            # 解析 dt 部分
            dt_str = parts[1]
            if !startswith(dt_str, "dt")
                @warn "Invalid dt format in line: $line"
                continue
            end
            δt = parse(Int, dt_str[3:end])  # 去掉 "dt" 前缀
            
            # 解析 sample 部分
            sample_str = parts[2]
            if !startswith(sample_str, "Sample")
                @warn "Invalid sample format in line: $line"
                continue
            end
            index = parse(Int, sample_str[7:end])  # 去掉 "Sample" 前缀
            
            push!(tasks, (δt, index))
        end
    end
    
    return tasks
end

# main parallel function
function compute_parallel_batch(L::Int64, τ::Float64, seed_range::UnitRange{Int}, D::Int64, δt::Int64)
    println("Starting parallel computation with $(nprocs()) processes")
    println("Parameters: L=$L, τ=$τ, D=$D, δt=$δt")
    println("Computing seeds: $(first(seed_range)) to $(last(seed_range))")
    
    # create parameter tuples for each task
    tasks = [(L, τ, index, D, δt) for index in seed_range]
    
    # using pmap run
    println("Submitting $(length(tasks)) tasks to worker processes...")
    results = pmap(compute_single_task, tasks, batch_size=1)
    
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

# batch mode for multiple δt values
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

function compute_missing_tasks_parallel(filename::String)
    println("Processing file: $filename")
    
    # 1. obtain L and inds from filename
    L, inds = parse_filename(basename(filename))
    println("Parsed parameters: L=$L, inds=$inds")
    
    # 2. obtain missing tasks
    missing_tasks = read_missing_tasks(filename)
    println("Found $(length(missing_tasks)) missing tasks")
    
    if isempty(missing_tasks)
        println("No missing tasks found!")
        return
    end
        
    τ = τlis[inds]
    D, _, _ = get_system_params(τ, L)
    
    println("Using τ=$τ, D=$D")
    
    # 4. 创建任务参数列表
    task_params = [(L, τ, index, D, δt) for (δt, index) in missing_tasks]
    
    # 5. 并行执行任务
    println("Starting parallel computation with $(nprocs()) processes...")
    println("Processing $(length(task_params)) tasks...")
    
    # 分批处理，避免内存问题
    batch_size = min(nworkers() * 4, 160)  # 每批最多160个任务
    total_success = 0
    total_failed = 0
    
    for batch_start in 1:batch_size:length(task_params)
        batch_end = min(batch_start + batch_size - 1, length(task_params))
        batch_tasks = task_params[batch_start:batch_end]
        
        println("Processing batch $(div(batch_start-1, batch_size) + 1): tasks $batch_start to $batch_end")
        
        # 使用 pmap 并行执行当前批次
        batch_results = pmap(compute_single_task, batch_tasks, batch_size=1)
        
        # 统计结果
        success_count = 0
        for (index, status, result) in batch_results
            if status == :success
                success_count += 1
                total_success += 1
            else
                total_failed += 1
                println("✗ Task failed - Index $index: $result")
            end
        end
        
        println("✓ Batch completed: $success_count/$(length(batch_tasks)) successful")
        
        # 强制垃圾回收
        @everywhere GC.gc()
    end
    
    println("\n" * "="^50)
    println("All batches completed!")
    println("Total successful: $total_success")
    println("Total failed: $total_failed")
    println("Success rate: $(round(total_success/(total_success + total_failed)*100, digits=1))%")
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
    
    if mode == 1
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])
        τ = τlis[τ_idx]
        D, _, _ = get_system_params(τ, L)
        # single δt value computation
        δt = parse(Int, ARGS[4])
        start_seed = parse(Int, ARGS[5])
        end_seed = parse(Int, ARGS[6])
        
        println("Single δt computation mode")
        compute_parallel_batch(L, τ, start_seed:end_seed, D, δt)
        
    elseif mode == 2
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])
        τ = τlis[τ_idx]
        D, _, _ = get_system_params(τ, L)
        # batch computation for multiple δt values
        start_seed = parse(Int, ARGS[4])
        end_seed = parse(Int, ARGS[5])
        δt_list = get_δt_list(L, τ)
        
        println("Batch computation mode")
        compute_parallel_multiple_dt(L, τ, start_seed:end_seed, D, δt_list)
    elseif mode == 3
        # compute missing tasks from file
        if length(ARGS) < 2
            error("Missing filename for mode 3")
        end
        
        filename = ARGS[2]
        compute_missing_tasks_parallel(filename)
    end
end