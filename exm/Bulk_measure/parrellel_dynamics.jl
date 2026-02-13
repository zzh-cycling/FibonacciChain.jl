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
                atanh(0.1)  => (collect(580:2:620)),
                atanh(0.2)  => (collect(100:2:120)),
                atanh(0.3)  => (collect(40:50)),
                atanh(0.4)  => (collect(25:35)),
                atanh(0.5)  => (collect(15:25)),
                atanh(0.6)  => (collect(5:15)),
                atanh(1/√2)  => (collect(5:10)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 8
        table = Dict(
                atanh(0.1)  => (collect(785:2:805)),
                atanh(0.2)  => (collect(135:2:145)),
                atanh(0.3)  => (collect(55:65)),
                atanh(0.4)  => (collect(35:45)),
                atanh(0.5)  => (collect(20:30)),
                atanh(0.6)  => (collect(8:16)),
                atanh(1/√2)  => (collect(5:12)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 10
        table = Dict(
                atanh(0.1)  => (collect(985:2:1005)),
                atanh(0.2)  => (collect(170:2:190)),
                atanh(0.3)  => (collect(70:80)),
                atanh(0.4)  => (collect(45:55)),
                atanh(0.5)  => (collect(25:35)),
                atanh(0.6)  => (collect(15:22)),
                atanh(1/√2)  => (collect(8:16)),
                atanh(0.8)  => (collect(5:10)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 12
        table = Dict(
                atanh(0.1)  => (collect(1185:2:1205)),
                atanh(0.2)  => (collect(210:2:230)),
                atanh(0.3)  => (collect(85:95)),
                atanh(0.4)  => (collect(55:65)),
                atanh(0.5)  => (collect(35:45)),
                atanh(0.6)  => (collect(18:25)),
                atanh(1/√2)  => (collect(10:18)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 14
        table = Dict(
                atanh(0.1)  => (collect(1380:2:1400)),
                atanh(0.2)  => (collect(245:2:265)),
                atanh(0.3)  => (collect(90:100)),
                atanh(0.4)  => (collect(65:75)),
                atanh(0.5)  => (collect(40:50)),
                atanh(0.6)  => (collect(20:28)),
                atanh(1/√2)  => (collect(15:23)),
                atanh(0.8)  => (collect(5:15)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 16
        table = Dict(
                atanh(0.1)  => (collect(1580:2:1600)),
                atanh(0.2)  => (collect(285:2:305)),
                atanh(0.3)  => (collect(115:125)),
                atanh(0.4)  => (collect(75:85)),
                atanh(0.5)  => (collect(45:55)),
                atanh(0.6)  => (collect(25:32)),
                atanh(1/√2)  => (collect(16:24)),
                atanh(0.8)  => (collect(8:16)),)
        δtlis = get(table, τ, collect(1:8))
    elseif L == 18
        table = Dict(
                atanh(0.1)  => (collect(1780:2:1800)),
                atanh(0.2)  => (collect(320:2:340)),
                atanh(0.3)  => (collect(130:140)),
                atanh(0.4)  => (collect(85:95)),
                atanh(0.5)  => (collect(55:65)),
                atanh(0.6)  => (collect(28:35)),
                atanh(1/√2)  => (collect(18:25)),
                atanh(0.8)  => (collect(10:16)),
                atanh(0.9)  => (collect(5:10)),)
        δtlis = get(table, τ, collect(1:8))
    else
        δtlis = collect(1:10)
    end
    return  δtlis
end

@everywhere function get_correlation_dynamics_D(τ, L)
    cfg = Dict(
        atanh(0.3)  => 25L,
        atanh(0.4)  => 20L,
        atanh(0.5)  => 15L,
        atanh(0.6)  => 10L
    )
    t = get(cfg, τ, 0)
    return t
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

@everywhere function spatial_temporal_corr_varyingt(args::Tuple)
    try 
    L, τ, index, δt = args
    D = get_system_params(τ, L)[1]
    # | ----> |____| ----> |
    # 0       D   D+δt   D+δt+t  
    # compute how the spatial and temporal correlation changes with t, the evolution time after add two ref qubits. δt is the time interval between two ref qubits

    model = AnyonModel(FibonacciAnyon(), L; pbc=true)
 
    # 1). First evolve to steady state with D time steps
    sample = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$L/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "sample")
    println("Loaded sample for L=$(L), τ=$(τ), index=$(index)")
    initial_state = zeros(length(anyon_basis(model)))
    initial_state[1] = 1.0 # initial state is all zero state

    t = div(D, 2)
    pre_config = MeasureConfig(τ=τ, mode=:sample, t₂ = t, enable_τ_eff=false)
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
            ref_config = MeasureConfig(τ = τ, t₂ = t, t₁ = t, rng = rng, mode=:Born, x₂=L÷2+1)
            ref_mo = reference_evolution(model, pre_mo.state, ref_config, ref_sample)
            sample_layer, sample_free_energy = ref_mo.samples, ref_mo.free_energys  # to compute temporal correlation, add ref
            spatial = true
            temporal = false
            view(sample_free_energy, 1:D) .= view(Flis, :)
            view(sample_layer, 1:D, :) .= view(sample, :, :)
        else
            ref_config = MeasureConfig(τ = τ, t₂ = t + δt, t₁ = t, rng = rng, mode=:Born, x₂=L÷2+1, x₁ = L÷2+1)
            ref_mo = reference_evolution(model, pre_mo.state, ref_config, ref_sample)
            sample_layer, sample_free_energy = ref_mo.samples, ref_mo.free_energys  # to compute temporal correlation, add ref qubit at site L/2+1
            temporal = true
            spatial = false
            view(sample_free_energy, 1:D) .= view(Flis, :)
            view(sample_layer, 1:D, :) .= view(sample, :, :)
        end
        sysrdm = reference_rdm(model, collect(1:div(L,2)), ref_mo.state, traceref = false)
        eelis[idx] = ee(sysrdm)
        spatial_corr, temporal_corr = ref_correlation(model, ref_mo.state, spatial=spatial, temporal=temporal)
        temporal_corr_lis[idx] = temporal_corr
        spatial_corr_lis[idx] = spatial_corr
    end

    save("exm/data/Bulk_measure/spatial_temporal_corr_varying_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(2t1,L))_Samples$(index).jld", "temporal_corr_lis", temporal_corr_lis, "spatial_corr_lis", spatial_corr_lis, "eelis", eelis)
        return index, :success, nothing
    catch e
        return index, :error, e
    end
end

@everywhere γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
@everywhere τlis = atanh.(γlis)
@everywhere τlis[end] = 1000.0
@everywhere τlis[findfirst(γlis .== 1/√2)] = log(1 + √2)

# main parallel function
function compute_parallel_batch(L::Int64, τ::Float64, seed_range::Vector{Int}=collect(1:10000))
    println("Starting parallel computation with $(nprocs()) processes")
    println("Parameters: L=$L, τ=$τ")
    println("Computing seeds: $(seed_range[1]) to $(seed_range[end])")
    dtlis = get_δtL(τ, L)
    # create parameter tuples for each task
    tasks = [(L, τ, index, δt) for index in seed_range for δt in dtlis]
    
    # using pmap run
    println("Submitting $(length(tasks)) tasks to worker processes...")
    results = pmap(spatial_temporal_corr_varyingt, tasks, batch_size=1)
    
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
    τ = τlis[τ_idx]
    compute_parallel_batch(L, τ, collect(1:2))
end


