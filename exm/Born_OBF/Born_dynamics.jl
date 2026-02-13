using Distributed
using FibonacciChain
using LinearAlgebra
using JLD2
using Statistics
using Random

@everywhere begin 
using FibonacciChain
using LinearAlgebra
using JLD2
using Statistics
using Random
γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0
τlis[findfirst(γlis .== 1/√2)] = log(1 + √2)
λlis = vcat(collect(0.0:0.1:1.5), [11.0])

function get_dynamics_params(ind, λ)
    if ind == 1
            cfg = Dict(
                11.0 => (400, 14, 350),
            )
            t, step, start = get(cfg, λ, (200, 14, 180))
    elseif ind == 7
            cfg = Dict(
                12.0 => (10, 14, 10),
            )
            t, step, start = get(cfg, λ, (8, 14, 6))
    end
    inds = collect(1:step:t*step)
    avg_range = start:t
    return t, inds, avg_range
end

function born_dynamics_samples_generate(L::Int64, λ::Float64, ind::Int64, index::Int64)
        τ = τlis[ind]
        try
            rng = MersenneTwister(index)
            t, _, _ = get_dynamics_params(ind, λ)
            if λ >= 10.0
                model = AnyonModel(OBFAnyon(), L; λI=0.0, pbc=true)
            else
                model = AnyonModel(OBFAnyon(), L; λ=λ, pbc=true)
            end
            # NEED to set fermion parity sector here.
            st = ones(length(anyon_basis(model)))
            st ./= norm(st)
            
            config = MeasureConfig(τ=τ, mode=:Born, t₂=t*L, rng=rng)
            outcome = bulk_evolution(model, st, config)
            sample = outcome.samples
            sample_free_energy = outcome.free_energys
            
            halfchain_EE_tlis = outcome.entanglement_entropys
            final_state = outcome.state
            final_EElis = anyon_eelis(model, final_state)
            
            # Assume seed is the index
            save("./exm/data/OBF/Born_dynamics_records/L$(L)/τ$(τ)/λ$(λ)/t$(t)_samples$(index).jld2", "halfchain_EE_tlis", halfchain_EE_tlis, "final_EElis", final_EElis, "sample_free_energy", sample_free_energy, "sample", sample)
            
            return (L, λ, τ, index, :success, nothing)
        catch e
            return (L, λ, τ, index, :failed, e)
        end
end

# define a wrapper function for pmap
function process_task(task)
    L, λ, ind, index = task
    return born_dynamics_samples_generate(L, λ, ind, index)
end

function samples_collect(task)
    L, λ, ind = task

    t, _, timewindow = get_dynamics_params(ind, λ)
    τ = τlis[ind]
    samples_num = 30000
    measure_records_ensemble = Vector{BitMatrix}(undef, samples_num)
    ensemble_free_energy = Vector{Vector{Float32}}(undef, samples_num)
    ensemble_seed = zeros(samples_num)
    ensemble_EE_dynamics= zeros(samples_num, t*L) 
    ensemble_final_EElis = zeros(samples_num, L-1)
     for i in 1:samples_num
        sample, sample_free_energy, halfchain_EE_tlis, final_EElis = load("./exm/data/OBF/Born_dynamics_records/L$(L)/τ$(τ)/λ$(λ)/t$(t)_samples$(i).jld2", "sample", "sample_free_energy", "halfchain_EE_tlis", "final_EElis")
        measure_records_ensemble[i] = sample
        ensemble_free_energy[i] = sample_free_energy
        ensemble_EE_dynamics[i, :] = halfchain_EE_tlis
        ensemble_final_EElis[i, :] = final_EElis
        ensemble_seed[i] = i
    end

    bulk_meanEElis = mean(ensemble_final_EElis, dims=1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims=1)[:]
    ensemble_stderr_EElis = (std(ensemble_final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims=1) ./ sqrt(samples_num))[:]

    save("exm/data/OBF/Born_dynamics_records/L$(L)/τ$(τ)/ensemble_L$(L)_τ$(τ)_λ$(λ)_t$(t).jld2", "measure_records_ensemble", measure_records_ensemble, "ensemble_free_energy", ensemble_free_energy,     "ensemble_seed", ensemble_seed, "average_EE_tlis", average_EE_tlis, "stderr_EE_tlis", stderr_EE_tlis, "bulk_meanEElis", bulk_meanEElis, "ensemble_stderr_EElis",ensemble_stderr_EElis)
end

function data_process(L::Int, ind::Int64, λ::Float64)
    t, _, timewindow = get_dynamics_params(ind, λ)  # Adjusted time window for averaging
    τ = τlis[ind]
    data = load("exm/data/OBF/Born_dynamics_records/L$(L)/τ$(τ)/ensemble_L$(L)_τ$(τ)_λ$(λ)_t$(t).jld2")
    average_EE_tlis = data["average_EE_tlis"]
    stderr_EE_tlis = data["stderr_EE_tlis"]
    ensemble_free_energy = data["ensemble_free_energy"]
    bulk_meanEElis = data["bulk_meanEElis"]
    measure_records_ensemble = data["measure_records_ensemble"]
    ensemble_stderr_EElis = data["ensemble_stderr_EElis"]

    temp = hcat(ensemble_free_energy...) # fuse the free energy of each sample into a matrix 
    #  | -> sample
    #  | 
    #  ⬇️ time
    t1 = timewindow[1]
    t2 = timewindow[end]
    time_average_free_energy = mean(temp[collect(t1*L:t2*L-4), :], dims=1)  # each sample's time-averaged free energy
    bulk_FE = mean(time_average_free_energy) # average over samples
    bulk_FE_stderr = std(time_average_free_energy) / sqrt(size(temp, 2))
    time_FEstderr = (std(temp./2 , dims=2) ./ sqrt(size(temp, 2)))[:] 
    time_FElis = mean(temp, dims=2)[:] ./2 # average over samples vs time; S/2T
    
    save("exm/data/OBF/Born_dynamics_records/L$(L)/τ$(τ)/Observables_L$(L)_τ$(τ)_λ$(λ)_t$(t).jld2", 
        "average_EE_tlis", average_EE_tlis, 
        "stderr_EE_tlis", stderr_EE_tlis, 
        "bulk_meanEElis", bulk_meanEElis, 
        "ensemble_stderr_EElis", ensemble_stderr_EElis, 
        "bulk_FE", bulk_FE,
        "bulk_FE_stderr", bulk_FE_stderr, 
        "time_FEstderr", time_FEstderr, 
        "time_FElis", time_FElis)
end
end 


if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia -p N monitored_dynamics.jl L τ_idx λ_idx index_start index_end")
    println("Example: julia -p 16 monitored_dynamics.jl 12 7 11 1 1000")
else
    mode = parse(Int64, ARGS[1])
    if mode == 1
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])

        λlis = vcat(collect(0.0:0.1:1.5), [11.0])
        
        tasklis = [(L, λ, τ_idx) for λ in λlis]
        
        println("Total tasks: $(length(tasklis))")
        println("Number of workers: $(nworkers())")
        println("\nStarting parallel processing...")
        results = pmap(samples_collect, tasklis; batch_size=1)
        for (L, λ, ind) in tasklis
            data_process(L, ind, λ) 
        end

    elseif mode == 2
        L = parse(Int64, ARGS[2])
        τinds = parse(Int64, ARGS[3])
        λlis = vcat(collect(0.0:0.1:1.5), [11.0])
        index_start = parse(Int64, ARGS[4])
        index_end = parse(Int64, ARGS[5])
        indexlis = collect(index_start:index_end)
        
        # create task list
        taskslis = [(L, λ, τinds, indexlis[i]) for λ in λlis for i in eachindex(indexlis)]
        
        println("=== Parallel Sample Generation ===")
        println("L = $L, τ_idx = $τinds, λlis = $λlis")
        println("Sample index range: $(indexlis[1]) - $(indexlis[end])")
        println("Total tasks: $(length(taskslis))")
        println("Number of workers: $(nworkers())")
        
        
        # use pmap for parallel processing
        println("\nStarting parallel processing...")
        results = pmap(process_task, taskslis; batch_size=100)
        
        # count successes and failures
        failed_tasks = [(L_res, λ_res, τ_res, idx_res, error) 
                        for (L_res, λ_res, τ_res, idx_res, status, error) in results 
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
            for (i, (L_f, τ_f, λ_f, idx_f, err)) in enumerate(failed_tasks)
                println("Failed $i: L=$L_f, τ=$τ_f, λ=$λ_f, index=$idx_f")
                println("  Error: $err")
            end
            
            # save failed tasks to file
            failed_file = "failed_tasks_L$(L)_τidx$(τinds)_batch.txt"
            open(failed_file, "w") do io
                println(io, "# Failed Task List")
                println(io, "# Format: L τ_idx λ_idx sample_index  # Error Message")
                for (L_f, τ_f, λ_f, idx_f, err) in failed_tasks
                    println(io, "$L_f $τ_f $λ_f $idx_f  # Error: $err")
                end
            end
            println("\nFailed tasks saved to: $failed_file")
        end
    end
end
