using Distributed
using FibonacciChain
using ITensorMPS, ITensors
using JLD
using Statistics
using Random

@everywhere begin
using FibonacciChain
using ITensorMPS, ITensors
using JLD
using Statistics
using Random

function samples_generate(L::Int64, τ::Float64, index::Int64, seed::Int64, D::Int64=5L)
    try
        rng = MersenneTwister(seed)
        
        model = AnyonModel(FibonacciAnyon(), L; pbc=true)
        ψ, sites = initial_mps(model)
        
        config = MeasureConfig(τ=τ, mode=:Born, t₂=div(D,2), rng=rng, cutoff=1e-12, maxdim=500)
        @time mps_mo = bulk_evolution(model, ψ, sites, config)
        sample_measured_states, sample, sample_free_energy = mps_mo.states, mps_mo.samples, mps_mo.free_energys
        
        halfchain_EE_tlis = [ee_mps(j, div(L,2)) for j in sample_measured_states]
        final_state = sample_measured_states[end]
        final_EElis = anyon_eelis(model, final_state)

        save("exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld", 
        "sample", sample, "sample_free_energy", sample_free_energy, "seed", seed, 
        "halfchain_EE_tlis", halfchain_EE_tlis, "final_EElis ", final_EElis)

        return (L, τ, index, seed, :success, nothing)
    catch e
        return (L, τ, index, seed, :failed, e)
    end
end


function samples_collect(L::Int64, τ::Float64, D::Int64=5L)
    samples_num = 5000
    ensemble = Vector{Matrix{Int}}(undef, samples_num)
    ensemble_free_energy = Vector{Vector{Float64}}(undef, samples_num)
    ensemble_seed = Vector{Int64}(undef, samples_num)
    ensemble_EE_dynamics= zeros(samples_num, D) 
    ensemble_final_EElis = zeros(samples_num, L-1)

     for i in 1:samples_num
        @show i
        @time sample, sample_free_energy, seed, halfchain_EE_tlis, final_EElis = load("exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(i).jld", "sample", "sample_free_energy", "seed", "halfchain_EE_tlis", "final_EElis")
        ensemble[i] = sample
        ensemble_free_energy[i] = sample_free_energy
        ensemble_seed[i] = seed
        ensemble_EE_dynamics[i, :] = halfchain_EE_tlis
        ensemble_final_EElis[i, :] = final_EElis
    end

    bulk_meanEElis = mean(ensemble_final_EElis, dims=1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims=1)[:]
    ensemble_stderr_EElis = (std(ensemble_final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims=1) ./ sqrt(samples_num))[:]

    save("exm/data/Bulk_measure/monitored_dynamics_mps/ensemble_L$(L)_τ$(τ)_D$(div(D,L)).jld", 
    "ensemble", ensemble, "ensemble_free_energy", ensemble_free_energy, "ensemble_seed", ensemble_seed,  
    "average_EE_tlis", average_EE_tlis, "stderr_EE_tlis", stderr_EE_tlis, 
    "bulk_meanEElis", bulk_meanEElis, "ensemble_stderr_EElis",ensemble_stderr_EElis)
end

function process_data(L::Int64, D::Int64=25L, τ::Float64=log(1+ √2))
    # timewindow = 8L:35L-10
    timewindow = 5L:D-5  # Adjusted time window for averaging
    load_data_path = "exm/data/Bulk_measure/monitored_dynamics_mps/ensemble_L$(L)_τ$(τ)_D$(div(D,L)).jld"
    data = load(load_data_path)
    
    average_EE_tlis, stderr_EE_tlis = data["average_EE_tlis"], data["stderr_EE_tlis"]
    bulk_meanEElis, ensemble_stderr_EElis = data["bulk_meanEElis"], data["ensemble_stderr_EElis"]
    ensemble_free_energy, ensemble_seed = data["ensemble_free_energy"], data["ensemble_seed"]
    
    function check_duplicates(seeds)
        if length(seeds) != length(unique(seeds))
            duplicates = findall(x -> count(==(x), seeds) > 1, unique(seeds))
            duplicate_values = unique(seeds)[duplicates]
            println("WARNING: Found duplicate seeds: $duplicate_values")
            return true
        else
            println("No duplicate seeds found in $(length(seeds)) seeds.")
            return false
        end
    end
    
    # Check if there are duplicates in the ensemble_seed
    has_duplicates = check_duplicates(ensemble_seed)
    
    temp = hcat(ensemble_free_energy...)
    time_average_free_energy = mean(temp[timewindow, :], dims=1) 
    bulk_FE = mean(time_average_free_energy)
    bulk_FE_stderr = std(time_average_free_energy) / sqrt(size(temp, 2))
    time_FEstderr = (std(temp, dims=2) ./ sqrt(size(temp, 2)))[:]
    time_FElis = mean(temp, dims=2)[:]
    
    save("exm/data/Bulk_measure/monitored_dynamics_mps/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld", 
        "average_EE_tlis", average_EE_tlis, 
        "stderr_EE_tlis", stderr_EE_tlis, 
        "bulk_meanEElis", bulk_meanEElis, 
        "ensemble_stderr_EElis", ensemble_stderr_EElis, 
        "time_average_free_energy", time_average_free_energy, 
        "bulk_FE", bulk_FE,
        "bulk_FE_stderr", bulk_FE_stderr, 
        "time_FEstderr", time_FEstderr, 
        "time_FElis", time_FElis, 
        "ensemble_seed", ensemble_seed)
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


γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)

# define a wrapper function for pmap
function process_task(task)
    L, τ, index, seed, D = task
    return samples_generate(L, τ, index, seed, D)
end
end


if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia -p N monitored_dynamics_mps.jl mode L τ_idx index_start index_end")
    println("Example: julia -p 16 monitored_dynamics_mps.jl 2 10 7 1 100")
else
    mode = parse(Int64, ARGS[1])
    if mode == 1
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])
        τ = τlis[τ_idx]
        D, _, _ = get_system_params(τ, L)
        samples_collect(L, τ, D)
        process_data(L, D, τ)
    elseif mode == 2
        L = parse(Int64, ARGS[2])
        inds = parse(Int64, ARGS[3])
        τ = τlis[inds]
        index_start = parse(Int64, ARGS[4])
        index_end = parse(Int64, ARGS[5])
        indexlis = collect(index_start:index_end)
        seedlis = -indexlis

        D, _, _ = get_system_params(τ, L)
        
        println("=== Parallel Sample Generation (MPS) ===")
        println("L = $L, τ_idx = $inds, τ = $τ, D = $D")
        println("Sample index range: $(indexlis[1]) - $(indexlis[end])")
        println("Total tasks: $(length(indexlis))")
        println("Number of workers: $(nworkers())")
        
        # create task list
        taskslis = [(L, τ, indexlis[i], seedlis[i], D) for i in eachindex(indexlis)]
        
        # use pmap for parallel processing
        println("\nStarting parallel processing...")
        results = pmap(process_task, taskslis; batch_size=50)
        
        # count successes and failures
        failed_tasks = [(L_res, τ_res, idx_res, seed_res, error) 
                        for (L_res, τ_res, idx_res, seed_res, status, error) in results 
                        if status != :success]
        
        success_count = count(r -> r[5] == :success, results)
        failed_count = length(failed_tasks)
        
        # summary report
        println("\n=== Processing Complete ===")
        println("Total tasks: $(length(taskslis))")
        println("Successes: $success_count")
        println("Failures: $failed_count")
        
        if failed_count > 0
            println("\n=== Failed Task Details ===")
            for (i, (L_f, τ_f, idx_f, seed_f, err)) in enumerate(failed_tasks)
                println("Failed $i: L=$L_f, τ=$τ_f, index=$idx_f, seed=$seed_f")
                println("  Error: $err")
            end
        end
    else
        error("Unknown mode: $mode")
    end
end