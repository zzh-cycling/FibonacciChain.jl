using Distributed
using FibonacciChain
using LinearAlgebra
using JLD2
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

    function scaling_dimension_Born(
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

        initial_stlis = Vector{Vector{Float64}}(undef, 10) # initialized 10 vectors
        for (idx, st) in enumerate(initial_stlis)
            st = zeros(Float64, l)
            st[idx] = 1.0
            initial_stlis[idx] = st
        end

        for step in 1:t
            sample_layer = sample[2step-1:2step, :]
            T = transfer_matrix(model, τ, sample_layer)
            for (idx, st) in enumerate(initial_stlis)
                initial_stlis[idx] = T * st
            end
        end

    
    end

    function transfer_matrix_Born(
        L::Int,
        τ_idx::Int,
        index::Int;
    )
        τ = τlis[τ_idx]

        D, inds, avg_range = get_cfg_params_Born(τ_idx, L)
        t = div(D, 2)
        # path = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τ_idx)/t$(t)_samples$(index).jld"
        # data = load(path)
        # sample = data["sample"]
        sample = BitMatrix(ones(Int8, D, div(L, 2)))
        
        model = AnyonModel(FibonacciAnyon(), L; pbc = true)
        basis = anyon_basis(model)
        l = length(basis)
        k = min(10, l)

        # Initialize k product states (basis vectors)
        states = zeros(Float64, l, k)
        for i in 1:k
            states[i, i] = 1.0
        end

        free_energies = zeros(k, t)

        for step in 1:t
            sample_layer = sample[2step-1:2step, :]

            # Apply transfer matrix to each state
            for i in 1:k
                states[:, i] = sample_evolution_unnormalized(
                    model, states[:, i], sample_layer; τ=τ, enable_τ_eff=false
                )
            end

            # Compute probabilities (squared norms) before normalization
            probs = vec(sum(abs2, states, dims=1))
            free_energies[:, step] = -log.(probs)

            # Normalize each state
            for i in 1:k
                norm_i = sqrt(probs[i])
                if norm_i > 0
                    states[:, i] ./= norm_i
                end
            end

            # QR orthogonalize
            Q = qr(states).Q
            states = Q[:, 1:k]
        end

        return free_energies
    end

    function transfer_matrix_parellel(
        model::AnyonModel{AT},
        τ::Float64,
        sign::Bool,
        index::Int
    ) where {AT}
        basis = anyon_basis(model)
        l = length(basis)
        sample = sign ? BitMatrix(ones(Int8, 2, div(L, 2))) : BitMatrix(zeros(Int8, 2, div(L, 2)))
        
        st = zeros(Float64, l)
        st[index] = 1.0
        # First layer
        col_indices1 = FibonacciChain._get_sample_column_indices(model, 1)
        layer1_sample = BitVector(sample[1, col_indices1])
        out1 = FibonacciChain._apply_measurement_layer(
            model, τ, st, layer1_sample;
            layer_idx = 1, normalized = false,
        )
        # Second layer
        col_indices2 = FibonacciChain._get_sample_column_indices(model, 2)
        layer2_sample = BitVector(sample[2, col_indices2])
        out2 = FibonacciChain._apply_measurement_layer(
            model, τ, out1.state, layer2_sample;
            layer_idx = 2, normalized = false,
        )
    
        label = sign ? "af" : "fm"
        save("exm/data/Bulk_measure/transfer_matrix_parellel_label_$(label)/L$(L)/gammaind$(τ_idx)/stlis_L$(L)_index$(index).jld2", "st", out2.state)
    end
    
    function collect_transfer_matrix(L, τ_idx, sign::Bool)
        model = AnyonModel(FibonacciAnyon(), L)
        basis = anyon_basis(model)
        l = length(basis)

        T = zeros(l, l)
        label = sign ? "af" : "fm"
        for index in 1:l
            st = load("exm/data/Bulk_measure/transfer_matrix_parellel_label_$(label)/L$(L)/gammaind$(τ_idx)/stlis_L$(L)_index$(index).jld2", "st")
            T[:, index] = st
        end
        save("exm/data/Bulk_measure/transfer_matrix_parellel_label_$(label)/transfer_matrix_L$(L)_gammaind$(τ_idx).jld2", "T", T)
    end

    function compute_ps_spectrum_task(sign::Bool=true)
        Llis = sign ? collect(8:2:20) : collect(6:6:18)
        spectrum_lis = zeros(length(Llis), 10)
        τ = atanh(0.95) 
        for (i, L) in enumerate(Llis)
            @show L
            model = AnyonModel(FibonacciAnyon(), L)
            sample = sign ? BitMatrix(ones(Int8, 2, div(L, 2))) : BitMatrix(zeros(Int8, 2, div(L, 2)))
            T = transfer_matrix(model, τ, sample)
            energy, states = eigen(T)
            spectrum_lis[i, :] = -log.(real.(energy[end-9:end]))
        end
        save_path = sign ? "transfer_matrix_spectrum_post_selection_af.jld2" : "transfer_matrix_spectrum_post_selection_fm.jld2"
        save(save_path, "Llis", Llis, "spectrum_lis", spectrum_lis)        
    end

    # ---------------------------------------------------------------------------
    # Parallel task wrappers
    # ---------------------------------------------------------------------------
    function compute_transfer_matrix_task(task)
        L, τ_idx, index = task
        model = AnyonModel(FibonacciAnyon(), L)
        try
            sign = true
            transfer_matrix_parellel(model, τlis[τ_idx], sign, index)
            return (L, τ_idx, index, :success, nothing)
        catch e
            return (L, τ_idx, index, :failed, e)
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

        results = pmap(compute_transfer_matrix_task, tasks; batch_size=10)

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

    elseif mode == 4
        L         = parse(Int64, ARGS[2])
        τ_idx     = parse(Int64, ARGS[3])
        collect_transfer_matrix(L, τ_idx, false)
    elseif mode == 5
        L         = parse(Int64, ARGS[2])
        τ_idx     = parse(Int64, ARGS[3])
        sign = false
        label = sign ? "af" : "fm"
        T = load("exm/data/Bulk_measure/transfer_matrix_parellel_label_$(label)/transfer_matrix_L$(L)_gammaind$(τ_idx).jld2", "T")
        energy, states = Arpack.eigs(T, nev=10, which=:LR)
        save("exm/data/Bulk_measure/transfer_matrix_parellel_label_$(label)/tf_spectrum_L$(L)_gammaind$(τ_idx).jld2", "energy", real.(energy))
    end
end
