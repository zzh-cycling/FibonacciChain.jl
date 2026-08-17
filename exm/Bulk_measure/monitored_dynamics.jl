using Distributed
using FibonacciChain
using JLD
using JLD2
using Statistics
using Random
# using ClusterManagers

# const PROJECT_DIR = something(dirname(Base.active_project()), pwd())
# const NWORKERS = parse(Int, get(ENV, "SLURM_NTASKS", "512"))
# const CPUS_PER_TASK = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1"))
# addprocs(SlurmManager(NWORKERS), exeflags="--project=$(PROJECT_DIR) --threads=1")

const BULK_MEASURE_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BULK_MEASURE_CONFIG)

@everywhere begin
    # using Pkg
    # Pkg.activate($PROJECT_DIR; io=devnull)
    # const num_workers = nworkers()
    # @info("Number of workers: $num_workers")    
    using FibonacciChain
    using JLD
    using JLD2
    using Statistics
    using Random

    function process_merge_task(task)
        Lr, τr_idx = task
        try
            samples_collect_process_data(Lr, τr_idx)
            return (Lr, τr_idx, :success, nothing)
        catch e
            return (Lr, τr_idx, :failed, e)
        end
    end

    function get_FE_avg_range(τind, L)
        # avoid Int(x) on non-integer Float64 (e.g. 1.2*48 = 57.6)
        toidx(x) = floor(Int, x * L)

        avg_table = Dict(
            1 => toidx(100):2:(2500*L-10),
            2 => toidx(40):2:(500*L-10),
            3 => toidx(20):2:(200*L-10),
            4 => toidx(15):2:(100*L-10),
            5 => toidx(8):2:(80*L-10),
            6 => toidx(6):2:(45*L-10),
            7 => toidx(5):2:(35*L-10),
            8 => toidx(5):2:(25*L-10),
            9 => toidx(4):2:(8*L-10),
            10 => toidx(4):2:(8*L-10),
        )

        default_range = toidx(3):2:(5*L-10)
        return get(avg_table, τind, default_range)
    end

    function samples_generate(L::Int64, τ_idx::Int64, index::Int64)
        try
            τ = τlis_ext[τ_idx]
            rng = MersenneTwister(index)
            D, _, _ = get_cfg_params_Born(τ_idx, L)
            t = div(D, 2)
            model = fib_model(L)
            st = zeros(length(anyon_basis(model)))
            st[1] = 1.0

            config = MeasureConfig(τ = τ, mode = :Born, t₂ = t, rng = rng)
            outcome = bulk_evolution(model, st, config)
            sample = outcome.samples
            sample_free_energy = outcome.free_energys

            halfchain_EE_tlis = outcome.entanglement_entropys
            final_state = outcome.state
            final_EElis = anyon_eelis(model, final_state)

            out_dir = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τ_idx)"
            save(
                joinpath(out_dir, "t$(div(t,L))_samples$(index).jld"),
                "sample",
                sample,
                "sample_free_energy",
                sample_free_energy,
                "seed",
                index,
                "halfchain_EE_tlis",
                halfchain_EE_tlis,
                "final_EElis",
                final_EElis,
            )
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
        D, _, _ = get_cfg_params_Born(τind, L)
        t = div(D, 2L)
        dir_path = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τind)/"
        # dir_path = "exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/gammaind$(τind)"
        samples_num = length(
            filter(
                f -> startswith(f, "t$(t)_samples") && endswith(f, ".jld"),
                readdir(dir_path),
            ),
        )
        println("collecting $(samples_num) sample files")
        ensemble = Vector{BitMatrix}(undef, samples_num)
        ensemble_free_energy = Vector{Vector{Float32}}(undef, samples_num)
        ensemble_seed = zeros(samples_num)
        ensemble_EE_dynamics = zeros(samples_num, div(D, 2))
        ensemble_final_EElis = zeros(samples_num, L-1)

        existing_files = filter(
            f -> startswith(f, "t$(t)_samples") && endswith(f, ".jld"),
            readdir(dir_path),
        )
        for (i, fname) in enumerate(existing_files)
            sample, sample_free_energy, seed, halfchain_EE_tlis, final_EElis = load(
                joinpath(dir_path, fname),
                "sample",
                "sample_free_energy",
                "seed",
                "halfchain_EE_tlis",
                "final_EElis",
            )
            ensemble[i] = sample
            ensemble_free_energy[i] = sample_free_energy
            ensemble_seed[i] = seed
            if length(halfchain_EE_tlis) == div(D, 2)
                ensemble_EE_dynamics[i, :] = halfchain_EE_tlis
            else
                ensemble_EE_dynamics[i, :] = halfchain_EE_tlis[2:2:D]
            end
            #    ensemble_EE_dynamics[i, :] = halfchain_EE_tlis[2:2:D]
            ensemble_final_EElis[i, :] = final_EElis
        end

        bulk_meanEElis = mean(ensemble_final_EElis, dims = 1)[:]
        average_EE_tlis = mean(ensemble_EE_dynamics, dims = 1)[:]
        ensemble_stderr_EElis =
            (std(ensemble_final_EElis, dims = 1) ./ sqrt(samples_num))[:]
        stderr_EE_tlis = (std(ensemble_EE_dynamics, dims = 1) ./ sqrt(samples_num))[:]

        check_duplicates(ensemble_seed)

        timewindow = get_FE_avg_range(τind, L)
        temp = hcat(ensemble_free_energy...)
        time_average_free_energy = mean(temp[timewindow, :], dims = 1)
        bulk_FE = mean(time_average_free_energy)
        bulk_FE_stderr = std(time_average_free_energy) / sqrt(size(temp, 2))
        time_FEstderr = (std(temp, dims = 2) ./ sqrt(size(temp, 2)))[:]
        time_FElis = mean(temp, dims = 2)[:]

        save(
            "exm/data/Bulk_measure/monitored_dynamics/L$(L)/EE_FEdynamics_L$(L)_gamma$(τind)_t$(t).jld2",
            "average_EE_tlis",
            average_EE_tlis,
            "stderr_EE_tlis",
            stderr_EE_tlis,
            "bulk_meanEElis",
            bulk_meanEElis,
            "ensemble_stderr_EElis",
            ensemble_stderr_EElis,
            "time_average_free_energy",
            time_average_free_energy,
            "bulk_FE",
            bulk_FE,
            "bulk_FE_stderr",
            bulk_FE_stderr,
            "time_FEstderr",
            time_FEstderr,
            "time_FElis",
            time_FElis,
            "ensemble_seed",
            ensemble_seed,
        )
    end

    """
        collect_sample_y_expectations(L, τind)

    Replay every Born sample generated by `samples_generate(L, τind, index)` in
    fixed-sample (post-selection) mode and record the normalized topological
    charge expectation value `⟨Y⟩(t)` after every complete period.

    The original trajectory files are read-only. The returned per-trajectory
    dynamics, seeds, and source filenames are stored together in a separate
    JLD2 file so that every `⟨Y⟩(t)` remains traceable to its measurement
    record.
    """
    function collect_sample_y_expectations(L::Int64, τind::Int64)
        D, _, _ = get_cfg_params_Born(τind, L)
        periods = div(D, 2)
        time_in_L = div(periods, L)
        data_dir =
            "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τind)"
        isdir(data_dir) || error("Trajectory directory does not exist: $data_dir")

        sample_files = filter(
            file -> startswith(file, "t$(time_in_L)_samples") &&
                    endswith(file, ".jld"),
            readdir(data_dir),
        )
        isempty(sample_files) && error(
            "No t$(time_in_L) trajectory samples were found in $data_dir",
        )

        # Sort by the saved seed below; filename lexicographic order would place
        # samples100 before samples11.
        model = fib_model(L)
        initial_state = zeros(Float64, length(anyon_basis(model)))
        initial_state[1] = 1.0
        config = MeasureConfig(
            τ = τlis_ext[τind],
            mode = :sample,
            t₂ = periods,
            track_y_expectation = true,
        )

        samples_num = length(sample_files)
        seeds = Vector{Int}(undef, samples_num)
        y_expectation_dynamics = Matrix{Float32}(undef, samples_num, periods)

        for (i, file) in enumerate(sample_files)
            sample, seed = load(joinpath(data_dir, file), "sample", "seed")
            outcome = bulk_evolution(
                model,
                initial_state,
                config,
                BitMatrix(sample),
            )
            length(outcome.y_expectation_values) == periods || error(
                "Y-expectation dynamics length mismatch in $file: expected " *
                "$periods, got $(length(outcome.y_expectation_values))",
            )

            seeds[i] = Int(seed)
            y_expectation_dynamics[i, :] = outcome.y_expectation_values
        end

        length(unique(seeds)) == samples_num ||
            error("Duplicate trajectory seeds found in $data_dir")
        order = sortperm(seeds)
        seeds = seeds[order]
        y_expectation_dynamics = y_expectation_dynamics[order, :]
        sample_files = sample_files[order]

        average_y_expectation_values =
            vec(mean(y_expectation_dynamics; dims = 1))
        stderr_y_expectation_values =
            vec(std(y_expectation_dynamics; dims = 1, corrected = false)) ./
            sqrt(samples_num)
        final_y_expectation_values = y_expectation_dynamics[:, end]

        output_path = joinpath(
            data_dir,
            "Y_expectation_L$(L)_gamma$(τind)_t$(time_in_L).jld2",
        )
        JLD2.jldsave(
            output_path;
            L = L,
            τ_idx = τind,
            τ = τlis_ext[τind],
            gamma = tanh(τlis_ext[τind]),
            periods = periods,
            samples_num = samples_num,
            ensemble_seed = seeds,
            sample_files = sample_files,
            trajectory_y_expectation_values = y_expectation_dynamics,
            average_y_expectation_values = average_y_expectation_values,
            stderr_y_expectation_values = stderr_y_expectation_values,
            final_y_expectation_values = final_y_expectation_values,
        )
        return output_path
    end

end

if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage:")
    println("  mode 1: julia -p N monitored_dynamics.jl 1")
    println("  mode 2: julia -p N monitored_dynamics.jl 2 L τ_idx index_start index_end")
    println("  mode 3: julia monitored_dynamics.jl 3 L τ_idx")
else
    mode = parse(Int64, ARGS[1])
    if mode == 1
        Llis = collect(12:2:20)
        τlis_idx = collect(1:12)
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
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])
        index_start = parse(Int64, ARGS[4])
        index_end = parse(Int64, ARGS[5])
        indexlis = collect(index_start:index_end)
        seedlis = indexlis

        println("=== Parallel Sample Generation ===")
        println("L = $L, τ_idx = $τ_idx")
        println("Sample index range: $(indexlis[1]) - $(indexlis[end])")
        println("Total tasks: $(length(indexlis))")
        println("Number of workers: $(nworkers())")

        # create task list
        taskslis = [(L, τ_idx, indexlis[i]) for i in eachindex(indexlis)]

        # use pmap for parallel processing
        println("\nStarting parallel processing...")
        results = pmap(process_task, taskslis; batch_size = 100)

        # count successes and failures
        failed_tasks = [
            (L_res, τ_res, idx_res, error) for
            (L_res, τ_res, idx_res, status, error) in results if
            status != :success
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
            for (i, (L_f, τ_f, idx_f, err)) in enumerate(failed_tasks)
                println("Failed $i: L=$L_f, τ=$τ_f, index=$idx_f, seed=$seed_f")
                println("  Error: $err")
            end

            # save failed tasks to file
            failed_file = "failed_tasks_L$(L)_τidx$(inds)_batch$(index).txt"
            open(failed_file, "w") do io
                println(io, "# Failed Task List")
                println(io, "# Format: L τ_idx sample_index seed")
                for (L_f, τ_f, idx_f, err) in failed_tasks
                    println(io, "$L_f $inds $idx_f $seed_f  # Error: $err")
                end
            end
            println("\nFailed tasks saved to: $failed_file")
        end
    elseif mode == 3
        length(ARGS) == 3 || error("Mode 3 requires L and τ_idx")
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])
        output_path = collect_sample_y_expectations(L, τ_idx)
        println("Saved per-sample Y expectations to $output_path")
    end
end
