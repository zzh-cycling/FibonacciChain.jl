using Distributed
using FibonacciChain
using LinearAlgebra
using JLD2
using Statistics
using Random

const BORN_HEISENBERG_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BORN_HEISENBERG_CONFIG)

@everywhere begin
    using FibonacciChain
    using LinearAlgebra
    using JLD2
    using Statistics
    using Random

    function born_dynamics_samples_generate(L::Int64, Δ::Float64, ind::Int64, index::Int64)
        τ = τlis[ind]
        try
            rng = MersenneTwister(index)
            t, _, _ = get_born_dynamics_params(ind, L, Δ)
            model = heisenberg_model(L, Δ)
            st = ones(length(anyon_basis(model)))
            st ./= norm(st)

            config = MeasureConfig(τ = τ, mode = :Born, t₂ = t*L, rng = rng)
            outcome = bulk_evolution(model, st, config)
            sample = outcome.samples
            sample_free_energy = outcome.free_energys

            halfchain_EE_tlis = outcome.entanglement_entropys
            final_state = outcome.state
            final_EElis = anyon_eelis(model, final_state)

            # Assume seed is the index
            mkpath("./exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)/Δ$(Δ)")
            save(
                "./exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)/Δ$(Δ)/t$(t)_samples$(index).jld2",
                "halfchain_EE_tlis",
                halfchain_EE_tlis,
                "final_EElis",
                final_EElis,
                "sample_free_energy",
                sample_free_energy,
                "sample",
                sample,
            )

            return (L, Δ, τ, index, :success, nothing)
        catch e
            return (L, Δ, τ, index, :failed, e)
        end
    end

    # define a wrapper function for pmap
    function process_task(task)
        L, Δ, ind, index = task
        return born_dynamics_samples_generate(L, Δ, ind, index)
    end

    function samples_collect(task)
        L, Δ, ind = task

        t = get_born_dynamics_params(ind, L, Δ)[1]
        τ = τlis[ind]
        dir_path = "exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)/Δ$(Δ)"
        samples_num = length(
            filter(
                f -> startswith(f, "t$(t)_samples"),
                readdir(dir_path),
            ),
        )
        measure_records_ensemble = Vector{BitMatrix}(undef, samples_num)
        ensemble_free_energy = Vector{Vector{Float32}}(undef, samples_num)
        ensemble_seed = zeros(samples_num)
        ensemble_EE_dynamics = zeros(samples_num, t*L)
        ensemble_final_EElis = zeros(samples_num, L-1)

        existing_files = filter(
            f -> startswith(f, "t$(t)_samples"),
            readdir(dir_path),
        )
        for (i, fname) in enumerate(existing_files)
            sample, sample_free_energy, halfchain_EE_tlis, final_EElis = load(
                joinpath(dir_path, fname),
                "sample",
                "sample_free_energy",
                "halfchain_EE_tlis",
                "final_EElis",
            )
            measure_records_ensemble[i] = sample
            ensemble_free_energy[i] = sample_free_energy
            ensemble_EE_dynamics[i, :] = halfchain_EE_tlis
            ensemble_final_EElis[i, :] = final_EElis
            ensemble_seed[i] = i
        end

        bulk_meanEElis = mean(ensemble_final_EElis, dims = 1)[:]
        average_EE_tlis = mean(ensemble_EE_dynamics, dims = 1)[:]
        ensemble_stderr_EElis =
            (std(ensemble_final_EElis, dims = 1) ./ sqrt(samples_num))[:]
        stderr_EE_tlis = (std(ensemble_EE_dynamics, dims = 1) ./ sqrt(samples_num))[:]

        mkpath("exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)")
        save(
            "exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)/ensemble_L$(L)_τ$(τ)_Δ$(Δ)_t$(t).jld2",
            "measure_records_ensemble",
            measure_records_ensemble,
            "ensemble_free_energy",
            ensemble_free_energy,
            "ensemble_seed",
            ensemble_seed,
            "average_EE_tlis",
            average_EE_tlis,
            "stderr_EE_tlis",
            stderr_EE_tlis,
            "bulk_meanEElis",
            bulk_meanEElis,
            "ensemble_stderr_EElis",
            ensemble_stderr_EElis,
        )
    end

    function process_data(L::Int, ind::Int64, Δ::Float64)
        t, _, start = get_born_dynamics_params(ind, L, Δ)  # Adjusted time window for averaging
        τ = τlis[ind]
        data = load(
            "exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)/ensemble_L$(L)_τ$(τ)_Δ$(Δ)_t$(t).jld2",
        )
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

        timewindow = collect(start*L:L*t-2)

        time_FElis = mean(temp, dims = 2)[:]  # average over samples
        time_FEstderr = (std(temp, dims = 2) ./ sqrt(size(temp, 2)))[:]
        reshaped_time_FElis = sum(reshape(time_FElis, 2, L*t), dims=1)[:]
        time_FElis_window = reshaped_time_FElis[timewindow]
        bulk_FE = mean(time_FElis_window) ./2 # average over samples vs time; S/2T
        bulk_FE_stderr = std(time_FElis_window ./2) / sqrt(length(time_FElis_window))

        save(
            "exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)/Observables_L$(L)_τ$(τ)_Δ$(Δ)_t$(t).jld2",
            "average_EE_tlis",
            average_EE_tlis,
            "stderr_EE_tlis",
            stderr_EE_tlis,
            "bulk_meanEElis",
            bulk_meanEElis,
            "ensemble_stderr_EElis",
            ensemble_stderr_EElis,
            "bulk_FE",
            bulk_FE,
            "bulk_FE_stderr",
            bulk_FE_stderr,
            "time_FEstderr",
            time_FEstderr,
            "time_FElis",
            time_FElis,
        )
    end
end


if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia -p N Born_dynamics.jl mode L τ_idx index_start index_end")
    println("Example: julia -p 16 Born_dynamics.jl 2 12 7 1 1000")
else
    mode = parse(Int64, ARGS[1])
    if mode == 1
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])

        tasklis = [(L, Δ, τ_idx) for Δ in Δlis]

        println("Total tasks: $(length(tasklis))")
        println("Number of workers: $(nworkers())")
        println("\nStarting parallel processing...")
        results = pmap(samples_collect, tasklis; batch_size = 1)
        for (L, Δ, ind) in tasklis
            process_data(L, ind, Δ)
        end

    elseif mode == 2
        L = parse(Int64, ARGS[2])
        τinds = parse(Int64, ARGS[3])
        index_start = parse(Int64, ARGS[4])
        index_end = parse(Int64, ARGS[5])
        indexlis = collect(index_start:index_end)

        # create task list
        taskslis = [(L, Δ, τinds, indexlis[i]) for Δ in Δlis for i in eachindex(indexlis)]

        println("=== Parallel Sample Generation ===")
        println("L = $L, τ_idx = $τinds, Δlis = $Δlis")
        println("Sample index range: $(indexlis[1]) - $(indexlis[end])")
        println("Total tasks: $(length(taskslis))")
        println("Number of workers: $(nworkers())")


        # use pmap for parallel processing
        println("\nStarting parallel processing...")
        results = pmap(process_task, taskslis; batch_size = 100)

        # count successes and failures
        failed_tasks = [
            (L_res, Δ_res, τ_res, idx_res, error) for
            (L_res, Δ_res, τ_res, idx_res, status, error) in results if status != :success
        ]

        success_count = count(r -> r[5] == :success, results)
        failed_count = length(failed_tasks)

        # summary report
        println("\n=== Processing Complete ===")
        println("Total tasks: $(length(taskslis))")
        println("Successes: $success_count")
        println("Failures: $failed_count")

        if failed_count > 0
            println("\n=== Failed Task Details ===")
            for (i, (L_f, τ_f, Δ_f, idx_f, err)) in enumerate(failed_tasks)
                println("Failed $i: L=$L_f, τ=$τ_f, Δ=$Δ_f, index=$idx_f")
                println("  Error: $err")
            end

            # save failed tasks to file
            failed_file = "failed_tasks_L$(L)_τidx$(τinds)_batch.txt"
            open(failed_file, "w") do io
                println(io, "# Failed Task List")
                println(io, "# Format: L τ_idx Δ sample_index  # Error Message")
                for (L_f, τ_f, Δ_f, idx_f, err) in failed_tasks
                    println(io, "$L_f $τ_f $Δ_f $idx_f  # Error: $err")
                end
            end
            println("\nFailed tasks saved to: $failed_file")
        end
    end
end
