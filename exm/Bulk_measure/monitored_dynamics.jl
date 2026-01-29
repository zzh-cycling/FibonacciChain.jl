using Distributed
using FibonacciChain
using JLD
using Statistics
using Random

@everywhere begin
using FibonacciChain
using JLD
using Statistics
using Random

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0
τlis[findfirst(γlis .== 1/√2)] = log(1 + √2)

function get_system_params(τ, L)
    cfg = Dict(
        atanh(0.1)  => (2500L, 1000, 750L),
        atanh(0.2)  => (500L,  100, 120L),
        atanh(0.3)  => (120L,  48, 50L),
        atanh(0.4)  => (100L,  40, 40L),
        atanh(0.5)  => (80L,   32, 20L),
        atanh(0.6)  => (45L,   20, 15L),
        log(1 + √2) => (35L,   14, 10L),
        atanh(0.8)  => (25L,   10, 5L),
        atanh(0.9)  => (8L,    4, 2L),
        atanh(0.95) => (8L,    4, 2L),
        atanh(0.999)=> (5L,    2, 1L),
    )
    D, step, start = get(cfg, τ, (5L, 2, L))
    inds = collect(1:step:div(D,2))
    avg_range = start:div(D,2)-5
    return D, inds, avg_range
end

function samples_generate(L::Int64, τ::Float64, index::Int64, seed::Int64)
        try
            rng = MersenneTwister(seed)
            D, _, _ = get_system_params(τ, L)
            model = AnyonModel(FibonacciAnyon(), L; pbc=true)
            st = zeros(length(anyon_basis(model)))
            st[1] = 1.0
            
            config = MeasureConfig(τ=τ, mode=:Born, t₂=div(D,2), rng=rng)
            outcome = bulk_evolution(model, st, config)
            sample_measured_states = outcome.states
            sample = outcome.samples
            sample_free_energy = outcome.free_energys
            
            halfchain_EE_tlis = [ee(anyon_rdm(model, collect(1:div(L,2)), j)) for j in sample_measured_states]
            final_state = sample_measured_states[end]
            final_EElis = anyon_eelis(model, final_state)
            
            save("./exm/data/Bulk_measure/Observable_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "halfchain_EE_tlis", halfchain_EE_tlis, "final_EElis ", final_EElis, "seed", seed, "sample_free_energy", sample_free_energy)
            save("exm/data/Bulk_measure/Samples_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "sample", sample, "sample_free_energy", sample_free_energy, "seed", seed)
            
            return (L, τ, index, seed, :success, nothing)
        catch e
            return (L, τ, index, seed, :failed, e)
        end
end

# define a wrapper function for pmap
function process_task(task)
    L, τ, index, seed = task
    return samples_generate(L, τ, index, seed)
end

function samples_collect(L::Int64, τ::Float64)
    D = get_system_params(τ, L)[1]
    samples_num = 20000
    ensemble = Vector{BitMatrix}(undef, samples_num)
    ensemble_free_energy = Vector{Vector{Float32}}(undef, samples_num)
    ensemble_seed = zeros(Int64, samples_num)
     for i in 1:samples_num
        sample, sample_free_energy, seed = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(i).jld", "sample", "sample_free_energy", "seed")
        ensemble[i] = sample
        ensemble_free_energy[i] = sample_free_energy
        ensemble_seed[i] = seed
    end

    save("exm/data/Bulk_measure/Samples_monitored_dynamics/monitored_dynamics_ensemble_L$(L)_τ$(τ)_D$(div(D,L)).jld", "ensemble", ensemble, "ensemble_free_energy", ensemble_free_energy, "ensemble_seed", ensemble_seed)
end


function Observable_collect(L::Int64, τ::Float64)
    D = get_system_params(τ, L)[1]
    samples_num = 20000
    ensemble_free_energy = Vector{Vector{Float32}}(undef, samples_num)
    ensemble_seed = zeros(Int64, samples_num)
    ensemble_EE_dynamics= zeros(samples_num, div(D,2)) 
    ensemble_final_EElis = zeros(samples_num, L-1)

    for i in 1:samples_num
        halfchain_EE_tlis, final_EElis, seed, sample_free_energy = load("./exm/data/Bulk_measure/Observable_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(i).jld", "halfchain_EE_tlis", "final_EElis ", "seed",  "sample_free_energy")
        if length(halfchain_EE_tlis) == D
            halfchain_EE_tlis = halfchain_EE_tlis[2:2:end]
        end
        ensemble_EE_dynamics[i, :] = halfchain_EE_tlis
        ensemble_final_EElis[i, :] = final_EElis
        ensemble_seed[i] = seed
        ensemble_free_energy[i] = sample_free_energy
    end

    bulk_meanEElis = mean(ensemble_final_EElis, dims=1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims=1)[:]
    ensemble_stderr_EElis = (std(ensemble_final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims=1) ./ sqrt(samples_num))[:]

    
    save("exm/data/Bulk_measure/Observable_monitored_dynamics/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld", "average_EE_tlis", average_EE_tlis, "stderr_EE_tlis", stderr_EE_tlis, "bulk_meanEElis", bulk_meanEElis, "ensemble_stderr_EElis",ensemble_stderr_EElis, "ensemble_free_energy", ensemble_free_energy, "ensemble_seed", ensemble_seed)
end

function get_system_params(τ, L)
    cfg = Dict(
        atanh(0.1)  => (2500L, 1000, 750L),
        atanh(0.2)  => (500L,  100, 120L),
        atanh(0.3)  => (120L,  48, 50L),
        atanh(0.4)  => (100L,  40, 40L),
        atanh(0.5)  => (80L,   32, 20L),
        atanh(0.6)  => (45L,   20, 15L),
        log(1 + √2) => (35L,   14, 10L),
        atanh(0.8)  => (25L,   10, 5L),
        atanh(0.9)  => (8L,    4, 2L),
        atanh(0.95) => (8L,    4, 2L),
        atanh(0.999)=> (5L,    2, 1L),
    )
    D, step, start = get(cfg, τ, (5L, 2, L))
    inds = collect(1:step:div(D,2))
    avg_range = start:div(D,2)-5
    return D, inds, avg_range
end

function save_data_filename(L, τ, D)
    return "monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D, L)).jld"
end

function get_data_filename(L, τ, D)
    return "Born_Fibo_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D, L)).jld"
end

## == Process the data for entanglement entropy and free energy dynamics == ##
function process_data(L::Int64, τ::Float64=log(1+ √2))
    # timewindow = 8L:35L-10
    D, _, timewindow = get_system_params(τ, L)  # Adjusted time window for averaging
    DATA_DIR = "exm/data/Bulk_measure/Observable_monitored_dynamics/"
    load_data_path = joinpath(DATA_DIR, save_data_filename(L, τ, D))
    data = load(load_data_path)
    
    average_EE_tlis, stderr_EE_tlis = data["average_EE_tlis"], data["stderr_EE_tlis"] # S vs time
    bulk_meanEElis, ensemble_stderr_EElis = data["bulk_meanEElis"], data["ensemble_stderr_EElis"] # S(l) at final time slice
    ensemble_free_energy, ensemble_seed = data["ensemble_free_energy"], data["ensemble_seed"] # collection of free energy vs time for each sample and the corresponding seeds
    
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
    
    temp = hcat(ensemble_free_energy...) # fuse the free energy of each sample into a matrix 
    #  | -> sample
    #  | 
    #  ⬇️ time
    time_average_free_energy = mean(temp[timewindow, :], dims=1)  # each sample's time-averaged free energy
    bulk_FE = mean(time_average_free_energy) # average over samples
    bulk_FE_stderr = std(time_average_free_energy) / sqrt(size(temp, 2))
    time_FEstderr = (std(temp./2 , dims=2) ./ sqrt(size(temp, 2)))[:] 
    time_FElis = mean(temp, dims=2)[:] ./2 # average over samples vs time; S/2T

    save_data_path = joinpath(DATA_DIR, get_data_filename(L, τ, D))
    save(save_data_path, 
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
end

if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia -p N monitored_dynamics.jl L τ_idx index_start index_end")
    println("Example: julia -p 16 monitored_dynamics.jl 10 7 1 1000")
else
    mode = parse(Int64, ARGS[1])
    if mode == 1
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])
        τ = τlis[τ_idx]
        samples_collect(L, τ)
        Observable_collect(L, τ)
        # todoτlis = [τ]
        # taskslis = [(L, τ)]
        # results = pmap(samples_collect, [(L, τ)])
        # results = pmap(Observable_collect, [(L, τ)])
    elseif mode == 2
        L = parse(Int64, ARGS[2])
        inds = parse(Int64, ARGS[3])
        τ = τlis[inds]
        index_start = parse(Int64, ARGS[4])
        index_end = parse(Int64, ARGS[5])
        indexlis = collect(index_start:index_end)
        seedlis = -indexlis
        
        
        println("=== Parallel Sample Generation ===")
        println("L = $L, τ_idx = $inds, τ = $τ")
        println("Sample index range: $(indexlis[1]) - $(indexlis[end])")
        println("Total tasks: $(length(indexlis))")
        println("Number of workers: $(nworkers())")
        
        # create task list
        taskslis = [(L, τ, indexlis[i], seedlis[i]) for i in eachindex(indexlis)]
        
        # use pmap for parallel processing
        println("\nStarting parallel processing...")
        results = pmap(process_task, taskslis; batch_size=100)
        
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
            
            # save failed tasks to file
            failed_file = "failed_tasks_L$(L)_τidx$(inds)_batch$(index).txt"
            open(failed_file, "w") do io
                println(io, "# Failed Task List")
                println(io, "# Format: L τ_idx sample_index seed")
                for (L_f, τ_f, idx_f, seed_f, err) in failed_tasks
                    println(io, "$L_f $inds $idx_f $seed_f  # Error: $err")
                end
            end
            println("\nFailed tasks saved to: $failed_file")
        end
    end
end
