using Distributed

@everywhere using JLD
@everywhere using Statistics

@everywhere function corr_collect(arg::Tuple)
    L, τ, δt = arg
    D = get_system_params(τ, L)[1]
    samples_num = 1000
    println("Sample number: ", samples_num)
    
    success=0
        temporal_corr_ensemble = zeros(samples_num)
        spatial_corr_ensemble = zeros(samples_num)
        S_ensemble = zeros(samples_num)
        sample_free_energy_ensemble = zeros(samples_num, 2*(D+δt))

        for i in 1:samples_num
            # @show i
            try
                temporal_corr, spatial_corr, S, sample_free_energy = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D,2L))_Samples$(i).jld",  "temporal_corr", "spatial_corr", "S", "sample_free_energy")
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
        "stderr_free_energy_tlis", stderr_free_energy_tlis)
        println("Completed L=$(L), τ=$(τ), δt=$(δt)")
        else
            println("No successful samples loaded for L=$(L), τ=$(τ), δt=$(δt). Skipping save.")
        end
end



@everywhere function get_δtL(τ, L)
    if L == 6
        table = Dict(
                atanh(0.1)  => (collect(482:2:578)),
                atanh(0.2)  => (collect(100:2:120)),
                atanh(0.3)  => (collect(40:50)),
                atanh(0.4)  => (collect(25:35)),
                atanh(0.5)  => (collect(21:25)),
                atanh(0.6)  => (collect(5:15)),
                atanh(1/√2)  => (collect(5:10)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 8
        table = Dict(
                atanh(0.1)  => (collect(50:50:700)),
                atanh(0.2)  => (collect(135:2:145)),
                atanh(0.3)  => (collect(55:65)),
                atanh(0.4)  => (collect(35:45)),
                atanh(0.5)  => (collect(20:30)),
                atanh(0.6)  => (collect(8:16)),
                atanh(1/√2)  => (collect(5:12)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 10
        table = Dict(
                atanh(0.1)  => (collect(100:50:950)),
                atanh(0.2)  => (collect(170:2:190)),
                atanh(0.3)  => (collect(70:80)),
                atanh(0.4)  => (collect(45:55)),
                atanh(0.5)  => (collect(25:35)),
                atanh(0.6)  => (collect(15:22)),
                atanh(1/√2)  => (collect(11:16)),
                atanh(0.8)  => (collect(5:10)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 12
        table = Dict(
                atanh(0.1)  => (collect(300:50:1150)),
                atanh(0.2)  => (collect(228:2:230)),
                atanh(0.3)  => (collect(85:95)),
                atanh(0.4)  => (collect(55:65)),
                atanh(0.5)  => (collect(35:45)),
                atanh(0.6)  => (collect(18:25)),
                atanh(1/√2)  => (collect(10:18)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 14
        table = Dict(
                atanh(0.1)  => (collect(50:100:850)),
                atanh(0.2)  => (collect(253:2:265)),
                atanh(0.3)  => (collect(90:100)),
                atanh(0.4)  => (collect(65:75)),
                atanh(0.5)  => (collect(40:50)),
                atanh(0.6)  => (collect(20:28)),
                atanh(1/√2)  => (collect(15:23)),
                atanh(0.8)  => (collect(9:15)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 16
        table = Dict(
                atanh(0.1)  => (collect(80:100:580)),
                atanh(0.2)  => (collect(291:2:305)),
                atanh(0.3)  => (collect(115:125)),
                atanh(0.4)  => (collect(75:85)),
                atanh(0.5)  => (collect(45:55)),
                atanh(0.6)  => (collect(25:32)),
                atanh(1/√2)  => (collect(16:24)),
                atanh(0.8)  => (collect(9:16)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 18
        table = Dict(
                atanh(0.1)  => (collect(1780:2:1800)),
                atanh(0.2)  => (collect(320:2:340)),
                atanh(0.3)  => (collect(130:140)),
                atanh(0.4)  => (collect(86:95)),
                atanh(0.5)  => (collect(56:65)),
                atanh(0.6)  => (collect(29:35)),
                atanh(1/√2)  => (collect(20:25)),
                atanh(0.8)  => (collect(10:16)),
                atanh(0.9)  => (collect(10:10)),)
        δtlis = get(table, τ, collect(1:8))
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
    samples_num = 1000
    println("Sample number: ", samples_num)
    success=0

   
    temporal_corr_ensemble = zeros(samples_num, div(D,2))
    spatial_corr_ensemble = zeros(samples_num, div(D,2))
    S_ensemble = zeros(samples_num, div(D,2))

        for i in 1:samples_num
            # @show i
            try
                temporal_corr_lis, spatial_corr_lis, eelis= load("exm/data/Bulk_measure/spatial_temporal_corr_varying_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D,2L))_Samples$(i).jld",  "temporal_corr_lis", "spatial_corr_lis", "eelis")
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
        "stderr_EE", stderr_EE_tlis)
    println("Completed L=$(L), τ=$(τ), δt=$(δt)")
    else
        println("No successful samples loaded for L=$(L), τ=$(τ), δt=$(δt). Skipping save.")
    end
end

@everywhere function organize(args::Tuple)
    L, τ = args
    δtlis = get_δtL(τ, L)
    tcLlis = zeros(Float64, length(δtlis))
    tcstderrlis = zeros(Float64, length(δtlis))
    
    average_spatial_corr, spatial_corr_stderr = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt0_collect.jld", "average_spatial_corr", "spatial_corr_stderr")
    for (j, δt) in enumerate(δtlis)
        average_temporal_corr, temporal_corr_stderr = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)_collect.jld",  "average_temporal_corr", "temporal_corr_stderr")
        tcLlis[j] = average_temporal_corr
        tcstderrlis[j] = temporal_corr_stderr
    end

    save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/stc_L$(L)_τ$(τ).jld", "average_spatial_corr", average_spatial_corr, "spatial_corr_stderr", spatial_corr_stderr, "δtlis", δtlis, "tcLlis", tcLlis, "tcstderrlis", tcstderrlis)
end

@everywhere γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
@everywhere τlis = atanh.(γlis)
@everywhere τlis[end] = 1000.0 

Llis = collect(8:2:10)
jobs = []

for L in collect(18)
    for τ in τlis
        δt_list = get_δtL(τ, L)
        for δt in δt_list
            push!(jobs, (L, τ, δt))
        end
    end
end


result = pmap(corr_collect, jobs, batch_size=1) 
# result = pmap(dynamics_collect, jobs, batch_size=1)  # Adjust batch_size as needed

println("All tasks completed.")