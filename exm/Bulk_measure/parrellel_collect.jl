using Distributed

@everywhere using JLD
@everywhere using Statistics

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

@everywhere function corr_collect(arg::Tuple)
    L, τ, δt = arg
    D = get_system_params(τ, L)[1]
    D = div(D, 2) # true circuits depth
    D1 = D + get_correlation_dynamics_D(τ, L)
    samples_num = 10000
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
                atanh(0.2)  => collect(173:2:190),
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

@everywhere function dynamics_collect(arg::Tuple)
    L, τ, δt = arg
    D = get_system_params(τ, L)[1]
    D = div(D, 2) # true circuits depth
    D1 = D + get_correlation_dynamics_D(τ, L)
    samples_num = 1000
    println("Sample number: ", samples_num)
    success=0

   
    temporal_corr_ensemble = zeros(samples_num, div(D,2))
    spatial_corr_ensemble = zeros(samples_num, div(D,2))
    S_ensemble = zeros(samples_num, div(D,2))
        for i in 1:samples_num
            # @show i
            try
                temporal_corr_lis, spatial_corr_lis, eelis= load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D1,L))_Samples$(i).jld",  "temporal_corr_lis", "spatial_corr_lis", "eelis")
                temporal_corr_ensemble[i,:] += temporal_corr_lis
                spatial_corr_ensemble[i,:] += spatial_corr_lis
                S_ensemble[i,:] += eelis
                success += 1
            catch e
                println("Error loading sample $(i) for L=$(L), τ=$(τ), δt=$(δt): ", e)
            end
        end
    
        average_temporal_corr_tlis = mean(temporal_corr_ensemble, dims=1)[:]
        average_spatial_corr_tlis = mean(spatial_corr_ensemble, dims=1)[:]
        average_EE_tlis = mean(S_ensemble, dims=1)[:]
        temporal_corr_tlis_stderr = std(temporal_corr_ensemble, dims=1) ./ sqrt(samples_num)
        spatial_corr_tlis_stderr = std(spatial_corr_ensemble, dims=1) / sqrt(samples_num)
        stderr_EE_tlis = std(S_ensemble, dims=1) / sqrt(samples_num)
    
    if success == samples_num
        save("exm/data/Bulk_measure/spatial_temporal_corr_varying_Born/L$(L)/τ$(τ)/dt$(δt)_collect.jld", "average_temporal_corr", average_temporal_corr_tlis, 
        "temporal_corr_stderr", temporal_corr_tlis_stderr, 
        "average_spatial_corr", average_spatial_corr_tlis, 
        "spatial_corr_stderr", spatial_corr_tlis_stderr, 
        "average_EE", average_EE_tlis,
        "stderr_EE", stderr_EE_tlis,
        "samples_num", samples_num)
    println("Completed L=$(L), τ=$(τ), δt=$(δt)")
    else
        println("No successful samples loaded for L=$(L), τ=$(τ), δt=$(δt). Skipping save.")
    end
end

@everywhere function data_compress(arg::Tuple)
    L, τ, δt = arg
    D = get_system_params(τ, L)[1]
    D = div(D, 2) # true circuits depth
    D1 = D + get_correlation_dynamics_D(τ, L)
    # D1 = D
    samples_num = 10000
    println("Sample number: ", samples_num)
    
    success=0
        temporal_corr_ensemble = zeros(Float32, samples_num)
        spatial_corr_ensemble = zeros(Float32, samples_num)
        S_ensemble = zeros(Float32, samples_num)
        sample_free_energy_ensemble = zeros(Float32, samples_num, 2*(D1+δt+D))
        sample_ensemble = zeros(Bool, samples_num, 2*(D1+δt+D), div(L,2))

        for i in 1:samples_num
            try
                temporal_corr, spatial_corr, S, sample_free_energy, sample = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D1,L))_Samples$(i).jld",  "temporal_corr", "spatial_corr", "S", "sample_free_energy", "sample_layer")
                temporal_corr_ensemble[i] = Float32(temporal_corr)
                spatial_corr_ensemble[i] = Float32(spatial_corr)
                sample_free_energy_ensemble[i, :] = Float32.(sample_free_energy)
                sample_ensemble[i, :, :] = Bool.(sample)
                S_ensemble[i] = Float32(S)
                success += 1
            catch e
                println("Error loading sample $(i) for L=$(L), τ=$(τ), δt=$(δt): ", e)
            end
        end
    
        if success == samples_num
            save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/compressed_dt$(δt)_data.jld", 
        "temporal_corr_ensemble", temporal_corr_ensemble, 
        "spatial_corr_ensemble", spatial_corr_ensemble, 
        "S_ensemble", S_ensemble,
        "sample_free_energy_ensemble", sample_free_energy_ensemble,
        "sample_ensemble", sample_ensemble,
        "samples_num", samples_num)
        println("Completed L=$(L), τ=$(τ), δt=$(δt)")
        else
            println("No successful samples loaded for L=$(L), τ=$(τ), δt=$(δt). Skipping save.")
        end
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


@everywhere function compute_parallel_dt_GC(L::Int64, τ::Float64, δt_list::Vector{Int})
    total_tasks = length(δt_list)
    println("Starting computation for multiple δt values")
    println("Total tasks: $total_tasks")
    
    all_results = Dict()
    
    batch_size = 63  
    
    for δt_start in 1:batch_size:total_tasks
        jobs = []
        δt_end = min(δt_start + batch_size - 1, length(δt_list))
        println("\n" * "="^50)
        println("Processing δt = $δt_start")
        println("="^50)
        
        push!(jobs, [(L, τ, i) for i in δt_list[collect(δt_start:δt_end)]]...)
        @show jobs
        result = pmap(data_compress, jobs, batch_size=1)
        @everywhere GC.gc()  # 强制垃圾回收
    end
    
    return all_results
end

@everywhere function organize(args::Tuple)
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

@everywhere γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
@everywhere τlis = atanh.(γlis)
@everywhere τlis[end] = 1000.0 

Llis = collect(8:2:18)
jobs = []

for L in collect(18)
    for τ in τlis
        δt_list = get_δtL_Born(τ, L)
        for δt in δt_list
            push!(jobs, (L, τ, δt))
        end
    end
end


result = pmap(corr_collect, jobs, batch_size=1) 
# result = pmap(dynamics_collect, jobs, batch_size=1)  # Adjust batch_size as needed

println("All tasks completed.")