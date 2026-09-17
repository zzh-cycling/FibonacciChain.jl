using Distributed

# Add workers if not already added
nprocs() == 1 && addprocs(Sys.CPU_THREADS - 1)

const BULK_MEASURE_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BULK_MEASURE_CONFIG)

@everywhere begin
    using JLD2
    using LinearAlgebra
    using Statistics
    using FibonacciChain

    # Broadcast parameters to all workers

    function data_path(inds, L, index)
        if L <= 30
            periods, _, _ = get_cfg_params_Born(inds, L)
            DATA_path = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(inds)/t$(div(periods,L))_samples$(index).jld"
        else
            t, _, _ = get_mps_params_Born(inds, L)
            DATA_path = "exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/gammaind$(inds)/t$(div(t,L))_samples$(index).jld2"
        end
        return DATA_path
    end

    function save_data_path(L, inds)
        τ = τlis[inds]
        if L <= 32
            periods, _, _ = get_cfg_params_Born(inds, L)
            t = div(periods, L)
        else
            t, _, _ = get_mps_params_Born(inds, L)
        end
        return "exm/data/Bulk_measure/defect_scaling/L$(L)/gammaind$(inds)"
    end

    function defect_FE(L::Int64, τ_idx::Int64, index::Int64)
        try
            τ = τlis[τ_idx]
            sample_path = data_path(τ_idx, L, index)
            sample = load(sample_path, "sample")
            sample[1:2:end, end] .= BitVector(1 .- sample[1:2:end, end]) # twist field, defect line (only odd layer, may correspond to superposition of many different primary field)
            parity = foldr(*, 2 .*sample[1, :] .-1)

            FElis=load(sample_path, "sample_free_energy")
            
            t, _, _ = get_cfg_params_Born(τ_idx, L)
            model = fib_model(L)
            st = zeros(length(anyon_basis(model)))
            st[1] = 1.0

            config = MeasureConfig(τ = τ, mode = :sample, t₂ = t)
            outcome = bulk_evolution(model, st, config, sample)
            sample_free_energy = outcome.free_energys

            save_path = save_data_path(L, τ_idx)
            save(joinpath(save_path, "defect_FE_t$(div(t,L))_samples$(index).jld2"), "defect_FElis", sample_free_energy, "FElis", FElis, "parity", parity)
            return (L, τ_idx, index, :success, nothing)
        catch e
            return (L, τ_idx, index, :failed, e)
        end
    end

     # define a wrapper function for pmap
    function process_task(task)
        L, τ_idx, index = task
        return defect_FE(L, τ_idx, index)
    end

    function defect_FE_collect(L::Int64, τind::Int64)
        try
            periods, _, _ = get_cfg_params_Born(τind, L)
            t = div(periods, L)
            dir_path = save_data_path(L, τind)
        
            samples_num = length(
                filter(
                    f -> startswith(f, "defect_FE_t$(t)") && endswith(f, ".jld2"),
                    readdir(dir_path),
                ),
            )
            println("collecting $(samples_num) sample files")
            layers = 2 * periods
            defect_FE_lis = zeros(samples_num, layers)
            FE_lis = zeros(samples_num, layers)
            parity_lis = zeros(samples_num)
            existing_files = filter(
                f -> startswith(f, "defect_FE_t$(t)") && endswith(f, ".jld2"),
                readdir(dir_path),
            )
            for (i, fname) in enumerate(existing_files)
                data = load(joinpath(dir_path, fname))
                defect_FElis = data["defect_FElis"]
                FElis = data["FElis"]
                parity = data["parity"]
                defect_FE_lis[i, :] = defect_FElis
                FE_lis[i, :] = FElis
                parity_lis[i] = parity
            end
            save("exm/data/Bulk_measure/defect_scaling/L$(L)/defect_FE_t$(t)_gamma$(τind)ensemble.jld2",
                "defect_FE_lis", defect_FE_lis,
                "FE_lis", FE_lis,
                "parity_lis", parity_lis,
            )
            avg_range = get_cfg_params_Born(τind, L)[3]
            defect_FElis_averaged = mean(defect_FE_lis, dims=2)
            FElis_averaged = mean(FE_lis, dims=2)
            defect_FE1 = mean(defect_FElis_averaged[avg_range][1:2:end])
            defect_FE2 = mean(defect_FElis_averaged[avg_range][2:2:end])
            delta_FEstd1 = std(defect_FElis_averaged[avg_range][1:2:end] .- FElis_averaged[avg_range][1:2:end]) / length(defect_FElis_averaged[avg_range][1:2:end])
            delta_FEstd2 = std(defect_FElis_averaged[avg_range][2:2:end] .- FElis_averaged[avg_range][2:2:end]) / length(defect_FElis_averaged[avg_range][2:2:end])
            FE = mean(FElis_averaged[avg_range])
            
            save("exm/data/Bulk_measure/defect_scaling/L$(L)/defect_FE_L$(L)_t$(t)_gamma$(τind)_collected.jld2",
                "defect_FE1", defect_FE1,
                "defect_FE2", defect_FE2,
                "delta_FEstd1", delta_FEstd1,
                "delta_FEstd2", delta_FEstd2,
                "FE", FE,
                "parity_lis", parity_lis,
            )
            return (L, τind, :success, nothing)
        catch e
            return (L, τind, :failed, e)
        end
    end
    
    function process_merge_task(task)
        L, τ_idx = task
        return defect_FE_collect(L, τ_idx)
    end
end


if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia -p N monitored_dynamics.jl L τ_idx index_start index_end")
    println("Example: julia -p 16 monitored_dynamics.jl 10 7 1 1000")
else
    mode = parse(Int64, ARGS[1])
    if mode == 1
        Llis = collect(12:2:14)
        τlis_idx = collect(10:10)
        merge_tasks = [(ll, tidx) for ll in Llis for tidx in τlis_idx]

        merge_results = pmap(process_merge_task, merge_tasks; batch_size = 1)

        failed_merges =
            [(Lr, τr, err) for (Lr, τr, status, err) in merge_results if status != :success]
        success_count = count(r -> r[3] == :success, merge_results)

        println("\n=== Merge Complete ===")
        println("Successes: $success_count")
        println("Failures: $(length(failed_merges))")

        if !isempty(failed_merges)
            println("\n=== Failed Merge Details ===")
            for (i, (Lf, τf, err)) in enumerate(failed_merges)
                println("Failed $i: L=$Lf, τ_idx=$τf")
                println("  Error: $err")
            end
        end
    elseif mode == 2
        τ_idx = parse(Int64, ARGS[2])
        index_start = parse(Int64, ARGS[3])
        index_end = parse(Int64, ARGS[4])
        indexlis = collect(index_start:index_end)
        seedlis = indexlis
        Llis = vcat(collect(12:2:14))

        println("=== Parallel Sample Generation ===")
        println("τ_idx = $τ_idx")
        println("Sample index range: $(indexlis[1]) - $(indexlis[end])")
        println("Total tasks: $(length(indexlis))")
        println("Number of workers: $(nworkers())")

        # create task list
        taskslis = [(L, τ_idx, indexlis[i]) for i in eachindex(indexlis) for L in Llis]
        
        # use pmap for parallel processing
        println("\nStarting parallel processing...")
        results = pmap(process_task, taskslis; batch_size = 100)

        # count successes and failures
        failed_tasks = [
            (L_res, τ_res, idx_res, status, error) for
            (L_res, τ_res, idx_res, status, error) in results if
            status != :success
        ]

        success_count = count(r -> r[4] == :success, results)
        failed_count = length(failed_tasks)

        # summary report
        println("\n=== Processing Complete ===")
        println("Total tasks: $(length(taskslis))")
        println("Successes: $success_count")
        println("Failures: $failed_count")

        if failed_count > 0
            println("\n=== Failed Task Details ===")
            for (i, (L_f, τ_f, idx_f, status_f, err)) in enumerate(failed_tasks)
                println("Failed $i: L=$L_f, τ=$τ_f, index=$idx_f, status=$status_f")
                println("  Error: $err")
            end

            # save failed tasks to file
            failed_file = "failed_tasks_τidx$(τ_idx)_batch_$(index_start)_$(index_end).txt"
            open(failed_file, "w") do io
                println(io, "# Failed Task List")
                println(io, "# Format: L τ_idx sample_index status error")
                for (L_f, τ_f, idx_f, status_f, err) in failed_tasks
                    println(io, "$L_f $τ_f $idx_f $status_f  # Error: $err")
                end
            end
            println("\nFailed tasks saved to: $failed_file")
        end
    end
end
