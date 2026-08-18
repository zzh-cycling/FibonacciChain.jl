using Distributed
using FibonacciChain
using JLD
using JLD2
using Statistics
using Random
using LinearAlgebra
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
    using LinearAlgebra

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

    function samples_collect_process_data(
        L::Int64,
        τind::Int64;
        sample_files::Union{Nothing,AbstractVector} = nothing,
        y_label::Union{Nothing,AbstractString} = nothing,
    )
        D, _, _ = get_cfg_params_Born(τind, L)
        t = div(D, 2L)
        dir_path = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τind)/"
        # dir_path = "exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/gammaind$(τind)"
        existing_files = isnothing(sample_files) ?
            filter(
                f -> startswith(f, "t$(t)_samples") && endswith(f, ".jld"),
                readdir(dir_path),
            ) : String.(sample_files)
        isempty(existing_files) && error("No sample files selected for L=$L, τ_idx=$τind")
        missing_files = filter(file -> !isfile(joinpath(dir_path, file)), existing_files)
        isempty(missing_files) || error("Missing sample files: $(join(missing_files, ", "))")

        samples_num = length(existing_files)
        println("collecting $(samples_num) sample files")
        ensemble = Vector{BitMatrix}(undef, samples_num)
        ensemble_free_energy = Vector{Vector{Float32}}(undef, samples_num)
        ensemble_seed = zeros(samples_num)
        ensemble_EE_dynamics = zeros(samples_num, div(D, 2))
        ensemble_final_EElis = zeros(samples_num, L-1)

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

        y_suffix = isnothing(y_label) ? "" : "_y$(y_label)"
        output_path = joinpath(
            "exm/data/Bulk_measure/monitored_dynamics/L$(L)",
            "EE_FEdynamics_L$(L)_gamma$(τind)$(y_suffix)_t$(t).jld2",
        )
        save(
            output_path,
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
        if !isnothing(y_label)
            JLD2.jldopen(output_path, "r+") do file
                file["y_label"] = String(y_label)
                file["source_sample_files"] = existing_files
            end
        end
        return output_path
    end

    """Collect EE/free-energy data separately in the two final-Y sectors."""
    function samples_collect_process_data_y_sectors(L::Int64, τind::Int64)
        D, _, _ = get_cfg_params_Born(τind, L)
        t = div(D, 2L)
        y_data_path = joinpath(
            "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τind)",
            "Y_expectation_L$(L)_gamma$(τind)_t$(t).jld2",
        )
        isfile(y_data_path) || error("Y-dynamics file does not exist: $y_data_path")
        y_data = JLD2.load(y_data_path)
        required = (
            "y_phi_sector_sample_files",
            "y_minus_invphi_sector_sample_files",
        )
        all(haskey(y_data, key) for key in required) ||
            error("Y-sector sample lists are missing; run mode 5 on $y_data_path")

        phi_path = samples_collect_process_data(
            L,
            τind;
            sample_files = y_data[required[1]],
            y_label = "phi",
        )
        minus_invphi_path = samples_collect_process_data(
            L,
            τind;
            sample_files = y_data[required[2]],
            y_label = "minus_invphi",
        )
        return (; phi = phi_path, minus_invphi = minus_invphi_path)
    end

    const Y_SECTOR_PHI = (1 + sqrt(5.0)) / 2
    const Y_SECTOR_MINUS_INVPHI = -inv(Y_SECTOR_PHI)
    const Y_REPLAY_CONTEXT_CACHE = Dict{NTuple{3,Int},Any}()
    const FINAL_Y_REPLAY_CONTEXT_CACHE = Dict{NTuple{3,Int},Any}()

    function _partition_y_sectors(final_y::AbstractVector)
        all(isfinite, final_y) || error("Y-sector classification requires finite values")
        positive = findall(>(0), final_y)
        negative = findall(<(0), final_y)
        zeros_num = count(iszero, final_y)
        zeros_num == 0 || error("Cannot assign $zeros_num trajectories with final Y = 0")
        isempty(positive) && error("No trajectories belong to the y = ϕ sector")
        isempty(negative) && error("No trajectories belong to the y = -1/ϕ sector")

        assignment = zeros(Int8, length(final_y))
        assignment[positive] .= 1
        assignment[negative] .= -1
        return positive, negative, assignment
    end

    function _mean_and_stderr(values::AbstractVector)
        isempty(values) && error("Cannot average an empty Y sector")
        return mean(values), std(values; corrected = false) / sqrt(length(values))
    end

    function _sector_metadata(final_y, seeds, sample_files)
        length(seeds) == length(final_y) || error("Seed count does not match Y values")
        length(sample_files) == length(final_y) ||
            error("Sample-file count does not match Y values")
        positive, negative, assignment = _partition_y_sectors(final_y)
        metadata = (
            y_sector_classification_rule = "sign(final_y_expectation_values)",
            y_sector_assignment = assignment,
            y_phi_sector_eigenvalue = Y_SECTOR_PHI,
            y_minus_invphi_sector_eigenvalue = Y_SECTOR_MINUS_INVPHI,
            y_phi_sector_samples_num = length(positive),
            y_minus_invphi_sector_samples_num = length(negative),
            y_phi_sector_ensemble_seed = seeds[positive],
            y_minus_invphi_sector_ensemble_seed = seeds[negative],
            y_phi_sector_sample_files = sample_files[positive],
            y_minus_invphi_sector_sample_files = sample_files[negative],
            y_phi_sector_final_y_expectation_values = final_y[positive],
            y_minus_invphi_sector_final_y_expectation_values = final_y[negative],
        )
        return metadata, positive, negative
    end

    function y_sector_dynamics_statistics(y_dynamics, seeds, sample_files)
        size(y_dynamics, 2) > 0 || error("Y dynamics must contain at least one period")
        metadata, positive, negative =
            _sector_metadata(y_dynamics[:, end], seeds, sample_files)
        positive_dynamics = view(y_dynamics, positive, :)
        negative_dynamics = view(y_dynamics, negative, :)
        return merge(metadata, (
            y_phi_sector_average_y_expectation_values =
                vec(mean(positive_dynamics; dims = 1)),
            y_minus_invphi_sector_average_y_expectation_values =
                vec(mean(negative_dynamics; dims = 1)),
            y_phi_sector_stderr_y_expectation_values =
                vec(std(positive_dynamics; dims = 1, corrected = false)) ./ sqrt(length(positive)),
            y_minus_invphi_sector_stderr_y_expectation_values =
                vec(std(negative_dynamics; dims = 1, corrected = false)) ./ sqrt(length(negative)),
        ))
    end

    function final_y_sector_statistics(final_y, seeds, sample_files)
        metadata, positive, negative = _sector_metadata(final_y, seeds, sample_files)
        positive_mean, positive_stderr = _mean_and_stderr(final_y[positive])
        negative_mean, negative_stderr = _mean_and_stderr(final_y[negative])
        return merge(metadata, (
            y_phi_sector_average_final_y_expectation = positive_mean,
            y_minus_invphi_sector_average_final_y_expectation = negative_mean,
            y_phi_sector_stderr_final_y_expectation = positive_stderr,
            y_minus_invphi_sector_stderr_final_y_expectation = negative_stderr,
        ))
    end

    function _replay_context(L, τind, periods, final_only)
        cache = final_only ? FINAL_Y_REPLAY_CONTEXT_CACHE : Y_REPLAY_CONTEXT_CACHE
        return get!(cache, (L, τind, periods)) do
            model = fib_model(L)
            initial_state = zeros(Float64, length(anyon_basis(model)))
            initial_state[1] = 1.0
            config = MeasureConfig(
                τ = τlis_ext[τind],
                mode = :sample,
                t₂ = periods,
                track_y_expectation = !final_only,
            )
            Y = final_only ? topological_charge_operator(model) : nothing
            (model, initial_state, config, Y)
        end
    end

    function _normalized_final_y_expectation(Y, state)
        norm_squared = real(dot(state, state))
        norm_squared > 0 || error("Cannot evaluate final Y for a zero-norm state")
        return Float32(real(dot(state, Y * state)) / norm_squared)
    end

    function _process_y_task(task, final_only)
        L, τind, periods, data_dir, file = task
        try
            model, initial_state, config, Y = _replay_context(L, τind, periods, final_only)
            sample, seed = load(joinpath(data_dir, file), "sample", "seed")
            outcome = bulk_evolution(model, initial_state, config, BitMatrix(sample))

            value = if final_only
                isempty(outcome.y_expectation_values) ||
                    error("Final-only replay unexpectedly recorded Y dynamics in $file")
                _normalized_final_y_expectation(Y, outcome.state)
            else
                length(outcome.y_expectation_values) == periods || error(
                    "Y-dynamics length mismatch in $file: expected $periods, " *
                    "got $(length(outcome.y_expectation_values))",
                )
                outcome.y_expectation_values
            end
            return (file, Int(seed), value, :success, nothing)
        catch e
            return (file, nothing, nothing, :failed, sprint(showerror, e))
        end
    end

    process_y_expectation_task(task) = _process_y_task(task, false)
    process_final_y_expectation_task(task) = _process_y_task(task, true)

    function _trajectory_inputs(L, τind)
        D, _, _ = get_cfg_params_Born(τind, L)
        periods = div(D, 2)
        time_in_L = div(periods, L)
        data_dir = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τind)"
        isdir(data_dir) || error("Trajectory directory does not exist: $data_dir")
        sample_files = filter(
            file -> startswith(file, "t$(time_in_L)_samples") && endswith(file, ".jld"),
            readdir(data_dir),
        )
        isempty(sample_files) && error("No t$(time_in_L) trajectory samples in $data_dir")
        return (; periods, time_in_L, data_dir, sample_files)
    end

    function _parallel_replay(L, τind, processor, heading, observable; batch_size)
        batch_size > 0 || error("batch_size must be positive")
        inputs = _trajectory_inputs(L, τind)
        tasks = [
            (L, τind, inputs.periods, inputs.data_dir, file)
            for file in inputs.sample_files
        ]
        samples_num = length(tasks)
        println("=== $heading ===")
        println("L = $L, τ_idx = $τind, trajectories = $samples_num")
        println("Number of workers: $(nworkers()), pmap batch size: $batch_size")

        results = pmap(processor, tasks; batch_size)
        failed = filter(result -> result[4] != :success, results)
        if !isempty(failed)
            details = join(["  $file: $err" for (file, _, _, _, err) in failed], "\n")
            error("$(length(failed)) of $samples_num $observable replays failed:\n$details")
        end

        seeds = Int[result[2] for result in results]
        length(unique(seeds)) == samples_num ||
            error("Duplicate trajectory seeds found in $(inputs.data_dir)")
        order = sortperm(seeds)
        return inputs, inputs.sample_files[order], seeds[order], results[order]
    end

    """Collect full per-period Y dynamics and both final-Y sectors."""
    function collect_sample_y_expectations(L::Int64, τind::Int64; batch_size::Int = 50)
        inputs, sample_files, seeds, results = _parallel_replay(
            L, τind, process_y_expectation_task,
            "Parallel Y-expectation replay", "Y-expectation"; batch_size,
        )
        y_dynamics = reduce(vcat, (permutedims(result[3]) for result in results))
        samples_num = length(seeds)
        average_y = vec(mean(y_dynamics; dims = 1))
        stderr_y = vec(std(y_dynamics; dims = 1, corrected = false)) ./ sqrt(samples_num)
        final_y = y_dynamics[:, end]
        sectors = y_sector_dynamics_statistics(y_dynamics, seeds, sample_files)
        output_path = joinpath(
            inputs.data_dir,
            "Y_expectation_L$(L)_gamma$(τind)_t$(inputs.time_in_L).jld2",
        )
        JLD2.jldsave(
            output_path;
            L,
            τ_idx = τind,
            τ = τlis_ext[τind],
            gamma = tanh(τlis_ext[τind]),
            periods = inputs.periods,
            samples_num,
            ensemble_seed = seeds,
            sample_files,
            trajectory_y_expectation_values = y_dynamics,
            average_y_expectation_values = average_y,
            stderr_y_expectation_values = stderr_y,
            final_y_expectation_values = final_y,
            sectors...,
        )
        return output_path
    end

    """Add or replace final-Y-sector statistics in an existing dynamics file."""
    function add_y_sector_averages!(output_path::AbstractString)
        isfile(output_path) || error("Y-expectation file does not exist: $output_path")
        data = JLD2.load(output_path)
        required = ("trajectory_y_expectation_values", "ensemble_seed", "sample_files")
        all(haskey(data, key) for key in required) || error("Missing datasets in $output_path")
        sectors = y_sector_dynamics_statistics(
            data[required[1]], data[required[2]], data[required[3]],
        )
        JLD2.jldopen(output_path, "r+") do file
            for (name, value) in pairs(sectors)
                key = String(name)
                haskey(file, key) && delete!(file, key)
                file[key] = value
            end
        end
        return output_path
    end

    """Collect only the final Y value, without recording per-period Y dynamics."""
    function collect_sample_final_y_expectations(
        L::Int64,
        τind::Int64;
        batch_size::Int = 50,
    )
        inputs, sample_files, seeds, results = _parallel_replay(
            L, τind, process_final_y_expectation_task,
            "Parallel final-Y replay", "final-Y"; batch_size,
        )
        final_y = Float32[result[3] for result in results]
        average_y, stderr_y = _mean_and_stderr(final_y)
        sectors = final_y_sector_statistics(final_y, seeds, sample_files)
        output_path = joinpath(
            inputs.data_dir,
            "Y_final_L$(L)_gamma$(τind)_t$(inputs.time_in_L).jld2",
        )
        JLD2.jldsave(
            output_path;
            L,
            τ_idx = τind,
            τ = τlis_ext[τind],
            gamma = tanh(τlis_ext[τind]),
            periods = inputs.periods,
            samples_num = length(seeds),
            ensemble_seed = seeds,
            sample_files,
            final_y_expectation_values = final_y,
            average_final_y_expectation = average_y,
            stderr_final_y_expectation = stderr_y,
            sectors...,
        )
        return output_path
    end

end

if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage:")
    println("  mode 1: julia -p N monitored_dynamics.jl 1")
    println("  mode 2: julia -p N monitored_dynamics.jl 2 L τ_idx index_start index_end")
    println("  mode 3: julia -p N monitored_dynamics.jl 3 L τ_idx")
    println("  mode 4: julia -p N monitored_dynamics.jl 4 L τ_idx")
    println("  mode 5: julia monitored_dynamics.jl 5 Y_dynamics_file")
    println("  mode 6: julia monitored_dynamics.jl 6 L τ_idx")
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
    elseif mode == 4
        length(ARGS) == 3 || error("Mode 4 requires L and τ_idx")
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])
        output_path = collect_sample_final_y_expectations(L, τ_idx)
        println("Saved per-sample final Y expectations to $output_path")
    elseif mode == 5
        length(ARGS) == 2 || error("Mode 5 requires a Y-dynamics JLD2 file")
        output_path = add_y_sector_averages!(ARGS[2])
        println("Added Y-sector averages to $output_path")
    elseif mode == 6
        length(ARGS) == 3 || error("Mode 6 requires L and τ_idx")
        L = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])
        output_paths = samples_collect_process_data_y_sectors(L, τ_idx)
        println("Saved y = ϕ sector data to $(output_paths.phi)")
        println("Saved y = -1/ϕ sector data to $(output_paths.minus_invphi)")
    else
        error("Unknown mode $mode; expected one of 1, 2, 3, 4, 5, or 6")
    end
end
