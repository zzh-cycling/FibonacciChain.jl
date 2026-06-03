using Distributed
using FibonacciChain
using LinearAlgebra
using JLD, JLD2
using Statistics
# using ClusterManagers

# const PROJECT_DIR = something(dirname(Base.active_project()), pwd())
# const NWORKERS = parse(Int, get(ENV, "SLURM_NTASKS", "512"))
# const CPUS_PER_TASK = parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1"))
# addprocs(SlurmManager(NWORKERS), exeflags="--project=$(PROJECT_DIR) --threads=1")

@everywhere begin
    # using Pkg
    # Pkg.activate($PROJECT_DIR; io=devnull)
    # const num_workers = nworkers()
    # @info("Number of workers: $num_workers")

    using FibonacciChain
    using LinearAlgebra
    using JLD, JLD2
    using Statistics

    γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
    τlis = atanh.(γlis)
    τlis[7] = log(1 + √2)  # atanh(1/√2) = log(1 + √2)
    τlis[end] = 1000.0     # Last value is for γ=1

    function get_cfg_params_Born(ind, L)
        cfg = Dict(
            1 => (2500L, 1000, 750L),
            2 => (500L, 100, 120L),
            3 => (200L, 40, 80L),
            4 => (100L, 40, 40L),
            5 => (80L, 32, 20L),
            6 => (45L, 20, 15L),
            7 => (35L, 14, 10L),
            8 => (25L, 10, 5L),
            9 => (8L, 4, 2L),
            10 => (8L, 4, 2L),
            11 => (5L, 2, 1L),
        )
        D, step, start = get(cfg, ind, (5L, 2, L))
        inds = collect(1:step:div(D, 2))
        avg_range = start:(div(D, 2)-5)
        return D, inds, avg_range
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
            )
        end
        t, step, start = get(cfg, τind, (2, 2, 1))
        inds = collect(1:step:(t*L))
        avg_range = Int(start*L):2:(Int(t*L)-4)
        return t, inds, avg_range
    end

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

    # ---------------------------------------------------------------------------
    # Build the transfer matrix for a fixed sample.
    # Column i = T(s) * e_i  (evolved basis vector, unnormalized).
    # ---------------------------------------------------------------------------
    function transfer_matrix(
        L::Int,
        τ_idx::Int,
        index::Int;
    )
        τ = τlis[τ_idx]

        D = get_cfg_params_Born(τ_idx, L)[1]
        t = div(D, 2L)
        path = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τ_idx)/t$(t)_samples$(index).jld"

        data = load(path)
        sample = data["sample"]

        model = AnyonModel(FibonacciAnyon(), L; pbc = true)
        model_basis = anyon_basis(model)
        l = length(model_basis)

        T = zeros(l, l)
        for i in 1:l
            st = zeros(l)
            st[i] = 1.0
            final_st = bulk_evolution_unnormalized(model, st, sample; τ = τ)
            T[:, i] = final_st
        end

        return T
    end
# Llis = collect(6:6:18)
# spectrum_lis = zeros(length(Llis), 10)
# τ = atanh(0.95) 
# for (i, L) in enumerate(Llis)
#     @show L
#     model = AnyonModel(FibonacciAnyon(), L)
#     sample = BitMatrix(zeros(Int8, 2, div(L, 2)))
#     T = transfer_matrix(model, τ, sample)
#     energy, states = eigen(T)
#     spectrum_lis[i, :] = -log.(real.(energy[end-9:end]))
# end
# save("transfer_matrix_spectrum_post_selection_fm.jld2", "Llis", Llis, "spectrum_lis", spectrum_lis)
    # ---------------------------------------------------------------------------
    # Parallel task wrappers
    # ---------------------------------------------------------------------------
    function compute_transfer_matrix_task(task)
        L, τ_idx, index = task
        try
            T = transfer_matrix(L, τ_idx, index)
            out_dir = "exm/data/Bulk_measure/transfer_matrix/L$(L)/gammaind$(τ_idx)"
            mkpath(out_dir)
            out_path = joinpath(out_dir, "transfer_matrix_samples$(index).jld")
            save(out_path, "T", T, "L", L, "τ_idx", τ_idx, "index", index)
            return (L, τ_idx, index, :success, nothing)
        catch e
            return (L, τ_idx, index, :failed, e)
        end
    end

    function compute_spectrum_task(task)
        L, τ_idx, index = task
        try
            T = transfer_matrix(L, τ_idx, index)
            energy = eigvals(T)
            λ = -log.(real(energy))
            return (L, τ_idx, index, :success, λ)
        catch e
            return (L, τ_idx, index, :failed, e)
        end
    end

    function compute_ps_spectrum_task(task)
        L, τ_idx, channel, λlength = task
        try
            T = transfer_matrix_post_selection(L, τ_idx, channel)
            D = get_cfg_params_Born(τ_idx, L)[1]
            energy = eigvals(T)
            # largest eigenvalues → smallest λ
            spectrum = -log.(real(energy[end-min(λlength, length(energy))+1:end])) / D * 2
            return (L, τ_idx, channel, :success, spectrum)
        catch e
            return (L, τ_idx, channel, :failed, e)
        end
    end

end

# =============================================================================
# Shell interface
# =============================================================================
if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia -p N transfer_matrix.jl mode [args...]")
    println("")
    println("Modes:")
    println("  1  L τ_idx index_start index_end")
    println("       Compute transfer matrices in parallel and save to file.")
    println("")
    println("  2  L τ_idx index_start index_end")
    println("       Compute eigenvalue spectra in parallel and save to file.")
    println("")
    println("  3  τ_idx channel L1 L2 ...")
    println("       Compute post-selection spectra for given L values.")
    println("       channel: true (AF) or false (FM).")
    println("")
    println("Examples:")
    println("  julia -p 16 transfer_matrix.jl 1 12 10 1 100")
    println("  julia -p 16 transfer_matrix.jl 2 12 10 1 100")
    println("  julia -p 4  transfer_matrix.jl 3 10 true 6 8 10 12")
else
    mode = parse(Int64, ARGS[1])

    if mode == 1
        # -------------------------------------------------------------------
        # Mode 1: parallel transfer matrix computation
        # -------------------------------------------------------------------
        L         = parse(Int64, ARGS[2])
        τ_idx     = parse(Int64, ARGS[3])
        idx_start = parse(Int64, ARGS[4])
        idx_end   = parse(Int64, ARGS[5])
        tasks = [(L, τ_idx, i) for i in idx_start:idx_end]

        println("=== Parallel Transfer Matrix Computation ===")
        println("L = $L, τ_idx = $τ_idx, τ = $(τlis[τ_idx])")
        println("Index range: $idx_start - $idx_end")
        println("Total tasks: $(length(tasks))")
        println("Number of workers: $(nworkers())")

        results = pmap(compute_transfer_matrix_task, tasks; batch_size=1)

        failed = [(Lr, τr, ir, err) for (Lr, τr, ir, status, err) in results if status != :success]
        success_count = count(r -> r[4] == :success, results)

        println("\n=== Complete ===")
        println("Successes: $success_count")
        println("Failures: $(length(failed))")

        if !isempty(failed)
            println("\n=== Failed Tasks ===")
            for (i, (Lf, τf, idxf, err)) in enumerate(failed)
                println("Failed $i: L=$Lf, τ_idx=$τf, index=$idxf")
                println("  Error: $err")
            end
            # save failed list
            failed_file = "failed_transfer_matrix_L$(L)_τ$(τ_idx)_$(idx_start)_$(idx_end).txt"
            open(failed_file, "w") do io
                println(io, "# Failed transfer matrix tasks")
                for (Lf, τf, idxf, err) in failed
                    println(io, "$Lf $τf $idxf  # $err")
                end
            end
            println("Failed list saved to: $failed_file")
        end

    elseif mode == 2
        # -------------------------------------------------------------------
        # Mode 2: parallel eigenvalue spectrum computation
        # -------------------------------------------------------------------
        L         = parse(Int64, ARGS[2])
        τ_idx     = parse(Int64, ARGS[3])
        idx_start = parse(Int64, ARGS[4])
        idx_end   = parse(Int64, ARGS[5])
        tasks = [(L, τ_idx, i) for i in idx_start:idx_end]

        println("=== Parallel Spectrum Computation ===")
        println("L = $L, τ_idx = $τ_idx, τ = $(τlis[τ_idx])")
        println("Index range: $idx_start - $idx_end")
        println("Total tasks: $(length(tasks))")
        println("Number of workers: $(nworkers())")

        results = pmap(compute_spectrum_task, tasks; batch_size=1)

        spectra  = Vector{Vector{Float64}}()
        indices  = Vector{Int}()
        for (Lr, τr, ir, status, result) in results
            if status == :success
                push!(spectra, result)
                push!(indices, ir)
            end
        end

        failed = [(Lr, τr, ir, err) for (Lr, τr, ir, status, err) in results if status != :success]
        success_count = length(spectra)

        println("\n=== Complete ===")
        println("Successes: $success_count")
        println("Failures: $(length(failed))")

        if success_count > 0
            out_dir = "exm/data/Bulk_measure/transfer_matrix/L$(L)/gammaind$(τ_idx)"
            mkpath(out_dir)
            out_path = joinpath(out_dir, "spectrum_$(idx_start)_$(idx_end).jld")
            save(out_path,
                 "spectra", spectra,
                 "indices", indices,
                 "L", L,
                 "τ_idx", τ_idx)
            println("Spectra saved to: $out_path")
        end

        if !isempty(failed)
            println("\n=== Failed Tasks ===")
            for (i, (Lf, τf, idxf, err)) in enumerate(failed)
                println("Failed $i: L=$Lf, τ_idx=$τf, index=$idxf")
                println("  Error: $err")
            end
        end

    elseif mode == 3
        # -------------------------------------------------------------------
        # Mode 3: post-selection spectra across L values
        # -------------------------------------------------------------------
        τ_idx   = parse(Int64, ARGS[2])
        channel = parse(Bool, ARGS[3])
        Llis    = [parse(Int64, arg) for arg in ARGS[4:end]]
        λlength = 10   # default number of eigenvalues to keep
        tasks   = [(L, τ_idx, channel, λlength) for L in Llis]

        println("=== Parallel Post-Selection Spectrum ===")
        println("τ_idx = $τ_idx, channel = $channel (", channel ? "AF" : "FM", ")")
        println("L values: $Llis")
        println("Total tasks: $(length(tasks))")
        println("Number of workers: $(nworkers())")

        results = pmap(compute_ps_spectrum_task, tasks; batch_size=1)

        L_success      = Vector{Int}()
        spectra_success = Vector{Vector{Float64}}()
        for (Lr, τr, ch, status, result) in results
            if status == :success
                push!(L_success, Lr)
                push!(spectra_success, result)
            end
        end

        failed = [(Lr, τr, ch, err) for (Lr, τr, ch, status, err) in results if status != :success]
        success_count = length(L_success)

        println("\n=== Complete ===")
        println("Successes: $success_count")
        println("Failures: $(length(failed))")

        if success_count > 0
            sign_str = channel ? "af" : "fm"
            out_name = "transfer_matrix_spectrum_post_selection_$(sign_str).jld"
            save(out_name,
                 "Llis", L_success,
                 "spectrum_lis", spectra_success,
                 "τ_idx", τ_idx)
            println("Saved to: $out_name")
        end

        if !isempty(failed)
            println("\n=== Failed Tasks ===")
            for (i, (Lf, τf, ch, err)) in enumerate(failed)
                println("Failed $i: L=$Lf, τ_idx=$τf, channel=$ch")
                println("  Error: $err")
            end
        end

    else
        println("Unknown mode: $mode")
        println("Valid modes are 1, 2, or 3.")
    end
end
