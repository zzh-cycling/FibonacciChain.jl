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

function get_system_params(τ_idx, L)
    cfg = Dict(
        1  => (1250L, 1000, 750L),
        2  => (250L,  100, 120L),
        3  => (100L,  40, 80L),
        4  => (50L,  40, 40L),
        5  => (40L,   32, 20L),
        6  => (22.5L,   20, 15L),
        7  => (17.5L,   14, 10L),
        8  => (12.5L,   10, 5L),
        9  => (4L,    4, 2L),
        10 => (4L,    4, 2L),
        11 => (2.5L,    2, 1L),
    )
    D, step, start = get(cfg, τ_idx, (2.5L, 2, L))
    inds = collect(1:step:div(D,2))
    avg_range = start:div(D,2)-5
    return D, inds, avg_range
end

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

function get_FE_avg_range(τind, L)
    # avoid Int(x) on non-integer Float64 (e.g. 1.2*48 = 57.6)
    toidx(x) = floor(Int, x * L)

    avg_table = Dict(
        2  => toidx(4):2:(40 * L - 10),
        3  => toidx(4):2:(40 * L - 10),
        4  => toidx(3):2:(28 * L - 10),
        5  => toidx(2):2:(20 * L - 10),
        6  => toidx(1.2):2:(12 * L - 10),
        7  => toidx(1.0):2:(10 * L - 10),
        8  => toidx(0.8):2:(7 * L - 10),
        9  => toidx(0.5):2:(3 * L - 10),
        10 => toidx(0.5):2:(2 * L - 10),
    )

    default_range = toidx(0.4):2:(2 * L - 10)
    return get(avg_table, τind, default_range)
end

function samples_generate(L::Int64, τ_idx::Int64, index::Int64)
        try
            τ = τlis[τ_idx]
            rng = MersenneTwister(index)
            t, _, _ = get_system_params(τ_idx, L)
            model = AnyonModel(FibonacciAnyon(), L; pbc=true)
            st = zeros(length(anyon_basis(model)))
            st[1] = 1.0
            
            config = MeasureConfig(τ=τ, mode=:Born, t₂=t, rng=rng)
            outcome = bulk_evolution(model, st, config)
            sample = outcome.samples
            sample_free_energy = outcome.free_energys
            
            halfchain_EE_tlis = outcome.entanglement_entropys
            final_state = outcome.state
            final_EElis = anyon_eelis(model, final_state)
            
            out_dir = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τ_idx)"
            save(joinpath(out_dir, "t$(div(t,L))_samples$(index).jld"), 
                "sample", sample, 
                "sample_free_energy", sample_free_energy, 
                "seed", index, 
                "halfchain_EE_tlis", halfchain_EE_tlis, 
                "final_EElis", final_EElis)
            return (L, τ_idx, index, :success, nothing)
        catch e
            return (L, τ_idx, index, :failed, e)
        end
end

# define a wrapper function for pmap
function process_task(task)
    L, τ_idx, index = task
    return samples_generate(L, τ_idx, index)
end

function samples_collect_process_data(L::Int64, τind::Int64)
           t, _, _= get_system_params(τind, L)
           dir_path = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τind)/"
            # dir_path = "exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/gammaind$(τind)"
           samples_num = length(filter(f -> startswith(f, "t$(t)_samples") && endswith(f, ".jld"), readdir(dir_path)))
           ensemble = Vector{BitMatrix}(undef, samples_num)
           ensemble_free_energy = Vector{Vector{Float32}}(undef, samples_num)
           ensemble_seed = zeros(samples_num)
           ensemble_EE_dynamics= zeros(samples_num, t * L) 
           ensemble_final_EElis = zeros(samples_num, L-1)

           existing_files = filter(f -> startswith(f, "t$(t)_samples") && endswith(f, ".jld2"), readdir(dir_path))
           for (i, fname) in enumerate(existing_files)
               sample, sample_free_energy, seed, halfchain_EE_tlis, final_EElis = load(joinpath(dir_path, fname), "sample", "sample_free_energy", "seed", "halfchain_EE_tlis", "final_EElis")
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

           check_duplicates(ensemble_seed)
           
           timewindow = get_FE_avg_range(τind, L)
           temp = hcat(ensemble_free_energy...)
           time_average_free_energy = mean(temp[2 .* timewindow, :], dims=1) 
           bulk_FE = mean(time_average_free_energy)
           bulk_FE_stderr = std(time_average_free_energy) / sqrt(size(temp, 2))
           time_FEstderr = (std(temp, dims=2) ./ sqrt(size(temp, 2)))[:]
           time_FElis = mean(temp, dims=2)[:]
           
           save("exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τind)/EE_FEdynamics_L$(L)_gamma$(τind)_t$(t).jld2", 
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
        process_data(L, τ)
        # todoτlis = [τ]
        # taskslis = [(L, τ)]
        # results = pmap(samples_collect, [(L, τ)])
        # results = pmap(Observable_collect, [(L, τ)])
    elseif mode == 2
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])
        index_start = parse(Int64, ARGS[4])
        index_end = parse(Int64, ARGS[5])
        indexlis = collect(index_start:index_end)
        seedlis = indexlis
        
        println("=== Parallel Sample Generation ===")
        println("L = $L, τ_idx = $τ_idx, τ = $τ")
        println("Sample index range: $(indexlis[1]) - $(indexlis[end])")
        println("Total tasks: $(length(indexlis))")
        println("Number of workers: $(nworkers())")
        
        # create task list
        taskslis = [(L, τ_idx, indexlis[i]) for i in eachindex(indexlis)]
        
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
