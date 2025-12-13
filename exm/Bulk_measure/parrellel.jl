using Distributed

@everywhere using FibonacciChain
@everywhere using JLD
@everywhere using Random
@everywhere using LinearAlgebra

println("requested workers: ", nworkers())
println("total procs:       ", nprocs())
@everywhere println("host: ", gethostname(), "  pid: ", getpid())

@everywhere function get_δtL_Born(τ, L)
    if L == 6
        table = Dict(
                atanh(0.1)  => vcat(collect(1:150), collect(152:2:620)),
                atanh(0.2)  => (collect(50:2:120)),
                atanh(0.3)  => (collect(1:45)),
                atanh(0.4)  => (collect(1:35)),
                atanh(0.5)  => (collect(1:25)),
                atanh(0.6)  => (collect(1:10)),)
        δtlis = get(table, τ, collect(1:10))
    elseif L == 8
        table = Dict(
                atanh(0.1)  => collect(50:25:500),
                atanh(0.2)  => (collect(65:2:145)),
                atanh(0.3)  => (collect(1:54)),
                atanh(0.4)  => (collect(1:45)),
                atanh(0.5)  => (collect(1:30)),
                atanh(0.6)  => (collect(1:16)),
                atanh(1/√2)  => (collect(1:12)),)
        δtlis = get(table, τ, collect(1:10))
    elseif L == 10
        table = Dict(
                atanh(0.1)  => collect(100:25:500),
                atanh(0.2)  => (collect(80:2:130)),
                atanh(0.3)  => (collect(30:60)),
                atanh(0.4)  => (collect(1:35)),
                atanh(0.5)  => (collect(1:24)),
                atanh(0.6)  => (collect(1:22)),
                atanh(1/√2)  => (collect(1:16)),)
        δtlis = get(table, τ, collect(1:10))
    elseif L == 12
        table = Dict(
                atanh(0.1)  => collect(300:25:600),
                atanh(0.2)  => (collect(100:2:136)),
                atanh(0.3)  => (collect(40:2:64)),
                atanh(0.4)  => (collect(15:35)),
                atanh(0.5)  => (collect(1:20)),
                atanh(0.6)  => (collect(1:25)),
                atanh(1/√2)  => (collect(1:18)),)
        δtlis = get(table, τ, collect(1:15))
    elseif L == 14
        table = Dict(
                atanh(0.1)  => (collect(50:25:550)),
                atanh(0.2)  => (collect(120:2:146)),
                atanh(0.3)  => (collect(45:80)),
                atanh(0.4)  => (collect(25:39)),
                atanh(0.5)  => (collect(5:24)),
                atanh(0.6)  => (collect(1:28)),
                atanh(1/√2)  => (collect(1:8)),
                atanh(0.8)  => (collect(1:15)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 16
        table = Dict(
                atanh(0.1)  => sort(vcat(collect(80:100:780), [640, 740])),
                atanh(0.2)  => (collect(135:2:160)),
                atanh(0.3)  => (vcat(collect(55:65), collect(67:2:75))),
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
                atanh(0.2)  => collect(173:2:193),
                atanh(0.3)  => (collect(68:4:88)),
                atanh(0.4)  => (collect(38:48)),
                atanh(0.5)  => (collect(24:34)),
                atanh(0.6)  => (collect(12:20)),
                atanh(1/√2)  => (collect(1:12)),)
        δtlis = get(table, τ, collect(1:10))
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

@everywhere function get_correlation_dynamics_D(τ, L)
    cfg = Dict(
        atanh(0.3)  => 50L,
        atanh(0.4)  => 40L,
        atanh(0.5)  => 30L,
        atanh(0.6)  => 20L
    )
    D = get(cfg, τ, 0)
    return D
end

@everywhere function compute_ratio(L::Int64, τ::Float64, index::Int64, D::Int64=16L, δt::Int64=2)
    pbc = true
    sample = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$L/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "sample")

    initial_state = zeros(length(anyon_basis(L, pbc)))
    initial_state[1] = 1.0 # initial state is all zero state

    rng = MersenneTwister(index)

    statelis, Flis = generate_state(τ, initial_state, sample, enable_τ_eff=false)
    D = div(D, 2) # true circuits depth
    D1 = D + get_correlation_dynamics_D(τ, L) # total evolution time after adding two ref qubits
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
    
    save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D1,L))_Samples$(index).jld", 
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
        @everywhere GC.gc()  # 强制垃圾回收
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
        δt_list = get_δtL_Born(τ, L)
        
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