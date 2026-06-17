using Distributed
using FibonacciChain
using ITensorMPS, ITensors
using JLD2
using Statistics
using Random
using ClusterManagers

# const PROJECT_DIR = something(dirname(Base.active_project()), pwd())
# const NWORKERS = parse(Int, get(ENV, "SLURM_NTASKS", "512"))
# const CPUS_PER_TASK = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1"))
# addprocs(SlurmManager(NWORKERS), exeflags = "--project=$(PROJECT_DIR) --threads=1")

@everywhere begin
    # using Pkg
    # Pkg.activate($PROJECT_DIR; io = devnull)
    # const num_workers = nworkers()
    # @info("Number of workers: $num_workers")
    using FibonacciChain
    using ITensorMPS, ITensors
    using JLD2
    using Statistics
    using Random

    function samples_generate_Fibo(L::Int64, τind::Int64, index::Int64, χ::Int64 = 500)
        τ = τlis[τind]
        try
            t, _, _ = get_mps_params_Born(τind, L)
            rng = MersenneTwister(index)

            model = AnyonModel(FibonacciAnyon(), L; pbc = true)
            ψ, sites = initial_mps(L)
            config = MeasureConfig(
                τ = τ,
                mode = :Born,
                t₂ = t*L,
                rng = rng,
                cutoff = 1e-12,
                maxdim = χ,
            )

            @time mps_mo = bulk_evolution(model, sites, ψ, config)
            sample, sample_free_energy = mps_mo.samples, mps_mo.free_energys
            halfchain_EE_tlis = mps_mo.entanglement_entropys

            # Compute final state EE profile
            final_EElis = anyon_eelis(model, mps_mo.state)

            save(
                "exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/gammaind$(τind)/chi$(χ)/t$(t)_samples$(index)_chi$(χ).jld2",
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

            return (L, τ, index, :success, nothing)
        catch e
            return (L, τ, index, :failed, e)
        end
    end

    function samples_collect_process_data(L::Int64, τind::Int64, χ::Int64 = 500)
        t, _, _ = get_mps_params_Born(τind, L)
        dir_path = "exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/gammaind$(τind)/chi$(χ)"
        # dir_path = "exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/gammaind$(τind)"
        samples_num = length(
            filter(
                f -> startswith(f, "t$(t)_samples") && endswith(f, "_chi$(χ).jld2"),
                readdir(dir_path),
            ),
        )
        ensemble = Vector{BitMatrix}(undef, samples_num)
        ensemble_free_energy = Vector{Vector{Float32}}(undef, samples_num)
        ensemble_seed = zeros(samples_num)
        ensemble_EE_dynamics = zeros(samples_num, t * L)
        ensemble_final_EElis = zeros(samples_num, L-1)

        existing_files = filter(
            f -> startswith(f, "t$(t)_samples") && endswith(f, "_chi$(χ).jld2"),
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
            ensemble_EE_dynamics[i, :] = halfchain_EE_tlis
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
        time_average_free_energy = mean(temp[2 .* timewindow, :], dims = 1)
        bulk_FE = mean(time_average_free_energy)
        bulk_FE_stderr = std(time_average_free_energy) / sqrt(size(temp, 2))
        time_FEstderr = (std(temp, dims = 2) ./ sqrt(size(temp, 2)))[:]
        time_FElis = mean(temp, dims = 2)[:]

        save(
            "exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/gammaind$(τind)/EE_FEdynamics_L$(L)_gamma$(τind)_t$(t)_chi$(χ).jld2",
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


    function get_mps_params_Born(τind, L)
        cfg = if L <= 32
            Dict(
                1 => (1250, 1000, 600),
                2 => (250, 100, 150),
                3 => (40, 48, 30),
                4 => (28, 40, 30),
                5 => (40, 32, 24),
                6 => (22, 20, 15),
                7 => (10, 14, 7),
                8 => (12, 10, 8),
                9 => (3, 4, 3),
                10 => (4, 4, 2.5),
                11 => (3, 2, 2),
                12 => (2, 2, 1)
            )
        else
            Dict(
                1 => (700, 1000, 500),
                2 => (150, 100, 100),
                3 => (40, 48, 30),
                4 => (28, 40, 22),
                5 => (20, 32, 16),
                6 => (12, 20, 10),
                7 => (10, 14, 7),
                8 => (7, 10, 5.5),
                9 => (3, 4, 2.5),
                10 => (2, 4, 1.5),
                11 => (2, 2, 1.5),
                12 => (2, 2, 1)
            )
        end
        t, step, start = get(cfg, τind, (8, 2, 2))
        inds = collect(1:step:(t*L))
        avg_range = Int(start*L):2:(Int(t*L)-4)
        return t, inds, avg_range
    end


    function get_FE_avg_range(τind, L)
        # avoid Int(x) on non-integer Float64 (e.g. 1.2*48 = 57.6)
        toidx(x) = floor(Int, x * L)

        avg_table = if L<=32
            Dict(
                1 => toidx(100):(1250*L-10),
                2 => toidx(20):(250*L-10),
                3 => toidx(10):(40*L-10),
                4 => toidx(3):(28*L-10),
                5 => toidx(2):(40*L-10),
                6 => toidx(1.2):(22*L-10),
                7 => toidx(1.0):(10*L-10),
                8 => toidx(0.8):(12*L-10),
                9 => toidx(0.5):(3*L-10),
                10 => toidx(0.5):(4*L-10),
                11 => toidx(0.5):(3*L-10),
            )
        else
            Dict(
                1 => toidx(100):(1250*L-10),
                2 => toidx(20):(250*L-10),
                3 => toidx(10):(40*L-10),
                4 => toidx(3):(28*L-10),
                5 => toidx(2):(20*L-10),
                6 => toidx(1.2):(12*L-10),
                7 => toidx(1.0):(10*L-10),
                8 => toidx(0.8):(7*L-10),
                9 => toidx(0.5):(3*L-10),
                10 => toidx(0.5):(2*L-10),
                11 => toidx(0.5):(3*L-10),
            )
        end

        default_range = toidx(0.4):(2*L-10)
        return get(avg_table, τind, default_range)
    end


    γlis = vcat([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1], collect(0.77:0.01:0.79), collect(0.81:0.01:0.82), [0.825], collect(0.83:0.01:0.84))
    τlis = atanh.(γlis)
    τlis[7] = log(1 + √2)
    τlis[12] = 1000.0

    function get_default_chi_Born(ind, L)
        if L == 32
            chi64_table = Dict(3 => 150, 4 => 150, 7 => 150, 9 => 200)
            return get(chi64_table, ind, 80)
        elseif L == 48
            chi48_table = Dict(1 => 150)
            return get(chi48_table, ind, 200)
        elseif L == 128 && ind == 10
            return 300
        elseif L == 64
            chi64_table = Dict(
                3 => 250,
                4 => 250,
                5 => 300,
                6 => 175,
                7 => 250,
                8 => 300,
                9 => 200,
                10 => 250,
            )
            return get(chi64_table, ind, 110)
        end
    end

    function process_merge_task(task)
        Lr, τr_idx = task
        try
            chi = get_default_chi_Born(τr_idx, Lr)
            samples_collect_process_data(Lr, τr_idx, chi)
            return (Lr, τr_idx, :success, nothing)
        catch e
            return (Lr, τr_idx, :failed, e)
        end
    end

    # define a wrapper function for pmap
    function process_task(task)
        L, τ, index, χ = task
        return samples_generate_Fibo(L, τ, index, χ)
    end
end


if length(ARGS) == 0
    println("No arguments provided.")
    println(
        "Usage: julia -p N monitored_dynamics_mps.jl mode L τ_idx index_start index_end",
    )
    println("Example: julia -p 16 monitored_dynamics_mps.jl 2 10 7 1 100")
else
    mode = parse(Int64, ARGS[1])
    if mode == 1

        Llis = [32, 48, 64]
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
        χ = parse(Int64, ARGS[4])
        index_start = parse(Int64, ARGS[5])
        index_end = parse(Int64, ARGS[6])
        indexlis = collect(index_start:index_end)

        println("=== Parallel Sample Generation (MPS) ===")
        println("L = $L, τ_idx = $τ_idx, χ = $χ")
        # println("Total parallel workers: $(nworkers())")
        # println("CPUs per worker: $CPUS_PER_TASK")
        # println("Total system CPUs: $(nworkers() * CPUS_PER_TASK) cores")
        println("Sample index range: $(indexlis[1]) - $(indexlis[end])")
        println("Total tasks: $(length(indexlis))")
        println("Number of workers: $(nworkers())")
        println("Batch size: 1")
        println("Estimated batches: $(cld(length(indexlis), 1))")
        println()
        taskslis = [(L, τ_idx, indexlis[i], χ) for i in eachindex(indexlis)]

        # use pmap for parallel processing
        println("\nStarting parallel processing...")
        results = pmap(process_task, taskslis; batch_size = 1)

        # count successes and failures
        failed_tasks = [
            (L_res, τ_res, idx_res, error) for
            (L_res, τ_res, idx_res, status, error) in results if status != :success
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
            for (i, (L_f, τ_f, idx_f, err)) in enumerate(failed_tasks)
                println("Failed $i: L=$L_f, τ=$τ_f, index=$idx_f")
                println("Error: $err")
            end
        end
    else
        error("Unknown mode: $mode")
    end
end
