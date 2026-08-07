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

const BULK_MEASURE_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BULK_MEASURE_CONFIG)

@everywhere begin
    # using Pkg
    # Pkg.activate($PROJECT_DIR; io=devnull)
    # const num_workers = nworkers()
    # @info("Number of workers: $num_workers")

    using FibonacciChain
    using LinearAlgebra
    using JLD, JLD2
    using Statistics
    using ITensors, ITensorMPS

    function get_default_chi_Born_tf(ind, L)
        if L == 32
            chi64_table = Dict(9 => 200)
            return get(chi64_table, ind, 150)
        elseif L == 48
            chi48_table = Dict(1 => 150)
            return get(chi48_table, ind, 200)
        elseif L == 30
            return 150
        elseif L == 34
            return 180
        elseif L == 36
            return 200
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
        else
            return 60
        end
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

    function compute_ps_spectrum_task(task)
        L, τ_idx, channel, λlength = task
        try
            τ = τlis[τ_idx]
            model = AnyonModel(FibonacciAnyon(), L)
            sample = channel ? BitMatrix(ones(Int8, 2, div(L, 2))) : BitMatrix(zeros(Int8, 2, div(L, 2)))
            T = transfer_matrix(model, τ, sample)
            energy = eigen(T).values
            # eigen() does NOT guarantee ordering for non-symmetric matrices. (But only due to the ABAB trotterization pattern, if T = AB could be written as T = A^(1/2) B A^(1/2), then eigen() would guarantee ordering.)
            # We must sort by absolute value to get the dominant eigenvalues.
            sorted_energy = sort(energy, by=abs, rev=true)
            spectrum = -log.(abs.(sorted_energy[1:λlength]))
            return (L, τ_idx, channel, :success, spectrum)
        catch e
            return (L, τ_idx, channel, :failed, e)
        end
    end

    
    function scaling_dimension_Born(
        L::Int,
        τ_idx::Int,
        index::Int;
        n_states::Int = 10,
    )   
    
        τ = τlis[τ_idx]
        D, inds, avg_range = get_cfg_params_Born(τ_idx, L)
        t = div(D, 2L)
        path = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τ_idx)/t$(t)_samples$(index).jld"
        data = load(path)
        sample = data["sample"]
        model = fib_model(L)
        FElis = transfer_matrix_subspace(model, τ, sample)
        # basis = anyon_basis(model)
        # l = length(basis)
        # function Tmap(v)        
        #     return sample_evolution_unnormalized(model, v, sample; τ=τ, enable_τ_eff=false)
        # end
        # st = zeros(l);st[1] = 1
        # tol_es = exp(-2*L^2)
        # krydim = 18L
        # iternum = 15000L
        # vals, vecs, info = eigsolve(Tmap, st, 10, :LM; tol=tol_es, krylovdim=krydim, ishermitian = false, maxiter=iternum);
        # FElis = -log.(vals)[1:10]/size(sample, 1)*2

        save("exm/data/Bulk_measure/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)/FElis_L$(L)_index$(index).jld2", "FElis", FElis)
        return FElis
    end

    function scaling_dimension_Born_mps(
        L::Int,
        τ_idx::Int,
        index::Int;
        n_states::Int = 10,
    )   
    
    
        χ = get_default_chi_Born_tf(τ_idx, L)
        τ = τlis[τ_idx]
        t, inds, avg_range = get_mps_params_Born(τ_idx, L)
        path = "exm/data/Bulk_measure/monitored_dynamics_mps/L$(L)/gammaind$(τ_idx)/chi$(χ)/t$(t)_samples$(index)_chi$(χ).jld2"
        data = load(path)
        sample = data["sample"]
        sites = siteinds("Qubit", L)
        model = fib_model(L)
        
        FElis = transfer_matrix_subspace_mps(
        model, sites, τ, sample; n_states = 10, cutoff = 1e-12, maxdim = χ)

        out_dir = "exm/data/Bulk_measure/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)/chi$(χ)"
        mkpath(out_dir)
        save(joinpath(out_dir, "FElis_L$(L)_index$(index).jld2"), "FElis", FElis)
        return FElis
    end

    function collect_transfer_matrix_Born(
        L::Int,
        τ_idx::Int;
        n_states::Int = 10,
        χ::Union{Int, Nothing} = nothing,
    )
        if isnothing(χ)
            # Exact diagonalization spectra
            D, _, _ = get_cfg_params_Born(τ_idx, L)
            t = div(D, 2)
            dir_path = "exm/data/Bulk_measure/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)"
            out_name = "ensemble_spectrum_L$(L)_gammaind$(τ_idx).jld2"
        else
            # MPS spectra
            t, _, _ = get_mps_params_Born(τ_idx, L)
            dir_path = "exm/data/Bulk_measure/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)/chi$(χ)"
            out_name = "ensemble_spectrum_L$(L)_gammaind$(τ_idx)_chi$(χ).jld2"
        end

        mkpath(dir_path)
        existing_files = filter(
            f -> startswith(f, "FElis_L$(L)_") && endswith(f, ".jld2"),
            readdir(dir_path),
        )
        samples_num = length(existing_files)
        println("collecting $(samples_num) sample files for L=$L, τ_idx=$τ_idx, χ=$χ")

        if samples_num == 0
            println("No sample files found in $(dir_path); skipping ensemble collection.")
            return nothing, nothing
        end

        spectra = zeros(n_states, t*L, samples_num)

        for (i, fname) in enumerate(existing_files)
            data = load(joinpath(dir_path, fname))
            FElis = data["FElis"]
            spectra[:, :, i] = FElis
        end

        # Ensemble statistics
        avg_spectrum = mean(spectra, dims = 3)[:, :, 1]
        stderr_spectrum = std(spectra, dims = 3)[:, :, 1] ./ sqrt(samples_num)

        # Save
        out_dir = "exm/data/Bulk_measure/tf_spectrum_Born/L$(L)"
        mkpath(out_dir)
        out_path = joinpath(out_dir, out_name)
        save(
            out_path,
            "avg_spectrum",
            avg_spectrum,
            "stderr_spectrum",
            stderr_spectrum,
            "n_samples",
            samples_num,
        )
        println("Ensemble spectrum saved to: $out_path")

        return avg_spectrum, stderr_spectrum
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
    
    function compute_tf_Born_task(task)
        L, τ_idx, index = task
        try
            result = scaling_dimension_Born(L, τ_idx, index)
            return (L, τ_idx, index, :success, result)
        catch e
            return (L, τ_idx, index, :failed, e)
        end
    end

    function compute_tf_Born_task_mps(task)
        L, τ_idx, index = task
        try
            result = scaling_dimension_Born_mps(L, τ_idx, index)
            return (L, τ_idx, index, :success, result)
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
    println("  6  L τ_idx [χ]")
    println("       Collect ensemble-averaged spectra from Born trajectories.")
    println("       If χ is omitted, collect exact spectra; otherwise collect MPS spectra.")
    println("")
    println("Examples:")
    println("  julia -p 16 transfer_matrix.jl 1 12 10 1 100")
    println("  julia -p 16 transfer_matrix.jl 2 12 10 1 100")
    println("  julia -p 4  transfer_matrix.jl 3 10 true 6 8 10 12")
    println("  julia         transfer_matrix.jl 6 8 10")
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

        results = pmap(compute_tf_Born_task, tasks; batch_size=50)

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

        results = pmap(compute_ps_spectrum_task, tasks; batch_size=50)

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
    elseif mode == 6
        # -------------------------------------------------------------------
        # Mode 6: collect ensemble-averaged Born spectra
        # -------------------------------------------------------------------
        L         = parse(Int64, ARGS[2])
        τ_idx     = parse(Int64, ARGS[3])
        χ         = length(ARGS) >= 4 ? parse(Int64, ARGS[4]) : nothing
        println("=== Collecting Ensemble-Averaged Born Spectrum ===")
        println("L = $L, τ_idx = $τ_idx, τ = $(τlis[τ_idx]), χ = $χ")
        avg, stderr = collect_transfer_matrix_Born(L, τ_idx; χ = χ)
        if isnothing(avg)
            println("\nNo sample files found; nothing to display.")
        else
            println("\n=== Average Spectrum (final time slice) ===")
            for i in 1:size(avg, 1)
                println("  λ_$i = $(avg[i, end]) ± $(stderr[i, end])")
            end
        end
     elseif mode == 7
        # -------------------------------------------------------------------
        # Mode 2: parallel eigenvalue spectrum computation
        # -------------------------------------------------------------------
        L         = parse(Int64, ARGS[2])
        τ_idx     = parse(Int64, ARGS[3])
        idx_start = parse(Int64, ARGS[4])
        idx_end   = parse(Int64, ARGS[5])
        tasks = [(L, τ_idx, i) for i in idx_start:idx_end]

        println("=== Parallel MPS Spectrum Computation ===")
        println("L = $L, τ_idx = $τ_idx, τ = $(τlis[τ_idx])")
        println("Index range: $idx_start - $idx_end")
        println("Total tasks: $(length(tasks))")
        println("Number of workers: $(nworkers())")

        results = pmap(compute_tf_Born_task_mps, tasks; batch_size=1)

        spectra  = Vector{Matrix{Float64}}()
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
            χ = get_default_chi_Born_tf(τ_idx, L)
            out_dir = "exm/data/Bulk_measure/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)/chi$(χ)"
            mkpath(out_dir)
            out_path = joinpath(out_dir, "spectrum_$(idx_start)_$(idx_end).jld2")
            save(out_path,
                 "spectra", spectra,
                 "indices", indices,
                 "L", L,
                 "τ_idx", τ_idx,
                 "χ", χ)
            println("Spectra saved to: $out_path")
        end

        if !isempty(failed)
            println("\n=== Failed Tasks ===")
            for (i, (Lf, τf, idxf, err)) in enumerate(failed)
                println("Failed $i: L=$Lf, τ_idx=$τf, index=$idxf")
                println("  Error: $err")
            end
        end
    end
end
