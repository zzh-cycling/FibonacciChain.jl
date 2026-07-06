using Distributed
using FibonacciChain
using LinearAlgebra
using ITensorMPS, ITensors
using JLD2
using Statistics
using Random

@everywhere begin
    using FibonacciChain
    using LinearAlgebra
    using ITensorMPS, ITensors
    using JLD2
    using Statistics
    using Random

    γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
    τlis = atanh.(γlis)
    τlis[end] = 1000.0
    τlis[findfirst(γlis .== 1/√2)] = log(1 + √2)

    function obf_model(L::Int, λ::Float64)
        if λ >= 10.0
            return AnyonModel(OBFAnyon(), L; λI = 0.0, pbc = true)
        else
            return AnyonModel(OBFAnyon(), L; λ = λ, pbc = true)
        end
    end

    function get_dynamics_params(ind, λ)
        if ind == 1
            cfg = Dict(11.0 => (400, 14, 350))
            t, step, start = get(cfg, λ, (200, 14, 180))
        elseif ind == 7
            cfg = Dict(11.0 => (18, 14, 4))
            t, step, start = get(cfg, λ, (18, 14, 4))
        end
        inds = collect(step:step:(t*step))
        avg_range = start*step:step:(t*step)-4
        return t, inds, avg_range
    end

    function get_dynamics_params_mps(ind, λ)
        if ind == 1
            cfg = Dict(11.0 => (400, 14, 350))
            t, step, start = get(cfg, λ, (200, 14, 180))
        elseif ind == 7
            cfg = Dict(11.0 => (10, 14, 10))
            t, step, start = get(cfg, λ, (8, 14, 6))
        else
            error("get_dynamics_params not configured for τ_idx=$ind and λ=$λ")
        end
        inds = collect(step:step:(t*step))
        avg_range = start:t
        return t, inds, avg_range
    end

    function scaling_dimension_Born(
        L::Int,
        τ_idx::Int,
        λ::Float64,
        index::Int;
        n_states::Int = 10,
    )
        τ = τlis[τ_idx]
        t, _, _ = get_dynamics_params(τ_idx, λ)
        model = obf_model(L, λ)

        path = "./exm/data/OBF/Born_dynamics_records/L$(L)/τ$(τ)/λ$(λ)/t$(t)_samples$(index).jld2"
        data = load(path)
        sample = data["sample"]

        FElis = transfer_matrix_subspace(model, τ, sample; n_states = n_states)

        out_dir = "exm/data/OBF/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)/λ$(λ)"
        mkpath(out_dir)
        save(joinpath(out_dir, "FElis_L$(L)_index$(index).jld2"), "FElis", FElis)
        return FElis
    end

    function scaling_dimension_Born_mps(
        L::Int,
        τ_idx::Int,
        λ::Float64,
        index::Int;
        n_states::Int = 10,
        χ::Int = 500,
        cutoff::Float64 = 1e-12,
    )
        τ = τlis[τ_idx]
        t, _, _ = get_dynamics_params_mps(τ_idx, λ)
        model = obf_model(L, λ)

        path = "exm/data/OBF/Born_dynamics_records_mps/L$(L)/gammaind$(τ_idx)/λ$(λ)/chi$(χ)/t$(t)_samples$(index)_chi$(χ).jld2"
        data = load(path)
        sample = data["sample"]

        sites = siteinds("Qubit", L)
        FElis = transfer_matrix_subspace_mps(
            model, sites, τ, sample;
            n_states = n_states,
            cutoff = cutoff,
            maxdim = χ,
        )

        out_dir = "exm/data/OBF/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)/λ$(λ)/chi$(χ)"
        mkpath(out_dir)
        save(joinpath(out_dir, "FElis_L$(L)_index$(index)_chi$(χ).jld2"), "FElis", FElis)
        return FElis
    end

    function collect_scaling_dimension_Born(
        L::Int,
        τ_idx::Int,
        λ::Float64;
        n_states::Int = 10,
        χ::Union{Int, Nothing} = nothing,
    )
        
        
        if isnothing(χ)
            dir_path = "exm/data/OBF/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)/λ$(λ)"
            prefix = "FElis_L$(L)_"
            out_name = "ensemble_spectrum_L$(L)_gammaind$(τ_idx)_λ$(λ).jld2"
            t, _, _ = get_dynamics_params(τ_idx, λ)
        else
            dir_path = "exm/data/OBF/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)/λ$(λ)/chi$(χ)"
            prefix = "FElis_L$(L)_"
            out_name = "ensemble_spectrum_L$(L)_gammaind$(τ_idx)_λ$(λ)_chi$(χ).jld2"
            t, _, _ = get_dynamics_params_mps(τ_idx, λ)
        end

        t_periods = t * L
        mkpath(dir_path)
        existing_files = filter(
            f -> startswith(f, prefix) && endswith(f, ".jld2"),
            readdir(dir_path),
        )
        samples_num = length(existing_files)
        println("collecting $(samples_num) sample files for L=$L, τ_idx=$τ_idx, λ=$λ")

        if samples_num == 0
            println("No sample files found in $(dir_path); skipping ensemble collection.")
            return
        end

        spectra = zeros(n_states, t_periods, samples_num)

        for (i, fname) in enumerate(existing_files)
            data = load(joinpath(dir_path, fname))
            FElis = data["FElis"]
            spectra[:, :, i] = FElis
        end

        avg_spectrum = mean(spectra, dims = 3)[:, :, 1]
        stderr_spectrum = std(spectra, dims = 3)[:, :, 1] ./ sqrt(samples_num)

        out_dir = "exm/data/OBF/tf_spectrum_Born/L$(L)"
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
    end

    function compute_tf_OBF_task(task)
        L, τ_idx, λ, index = task
        try
            result = scaling_dimension_Born(L, τ_idx, λ, index)
            return (L, τ_idx, λ, index, :success, result)
        catch e
            return (L, τ_idx, λ, index, :failed, e)
        end
    end

    function compute_tf_OBF_mps_task(task)
        L, τ_idx, λ, index, χ = task
        try
            result = scaling_dimension_Born_mps(L, τ_idx, λ, index; χ = χ)
            return (L, τ_idx, λ, index, :success, result)
        catch e
            return (L, τ_idx, λ, index, :failed, e)
        end
    end
end

# =============================================================================
# Shell interface
# =============================================================================
if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia -p N scaling_dimension.jl mode [args...]")
    println("")
    println("Modes:")
    println("  2  L τ_idx λ index_start index_end")
    println("       Compute exact OBF scaling dimensions in parallel.")
    println("")
    println("  7  L τ_idx λ χ index_start index_end")
    println("       Compute MPS OBF scaling dimensions in parallel.")
    println("")
    println("  6  L τ_idx λ [χ]")
    println("       Collect ensemble-averaged OBF spectra.")
    println("       If χ is omitted, collect exact spectra; otherwise collect MPS spectra.")
    println("")
    println("Examples:")
    println("  julia -p 16 scaling_dimension.jl 2 10 7 11.0 1 100")
    println("  julia -p 16 scaling_dimension.jl 7 12 1 0.856 500 1 100")
    println("  julia         scaling_dimension.jl 6 10 7 11.0")
    println("  julia         scaling_dimension.jl 6 12 1 0.856 500")
else
    mode = parse(Int64, ARGS[1])

    if mode == 1
        # -------------------------------------------------------------------
        # Mode 1: parallel exact OBF scaling dimension computation
        # -------------------------------------------------------------------
        L         = parse(Int64, ARGS[2])
        τ_idx     = parse(Int64, ARGS[3])
        λ         = parse(Float64, ARGS[4])
        idx_start = parse(Int64, ARGS[5])
        idx_end   = parse(Int64, ARGS[6])
        tasks = [(L, τ_idx, λ, i) for i in idx_start:idx_end]

        println("=== Parallel OBF Spectrum Computation (Exact) ===")
        println("L = $L, τ_idx = $τ_idx, τ = $(τlis[τ_idx]), λ = $λ")
        println("Index range: $idx_start - $idx_end")
        println("Total tasks: $(length(tasks))")
        println("Number of workers: $(nworkers())")

        results = pmap(compute_tf_OBF_task, tasks; batch_size = 50)

        spectra = Vector{Matrix{Float64}}()
        indices = Vector{Int}()
        for (Lr, τr, λr, ir, status, result) in results
            if status == :success
                push!(spectra, result)
                push!(indices, ir)
            end
        end

        failed = [
            (Lr, τr, λr, ir, err)
            for (Lr, τr, λr, ir, status, err) in results if status != :success
        ]
        success_count = length(spectra)

        println("\n=== Complete ===")
        println("Successes: $success_count")
        println("Failures: $(length(failed))")

        if success_count > 0
            out_dir = "exm/data/OBF/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)/λ$(λ)"
            mkpath(out_dir)
            out_path = joinpath(out_dir, "spectrum_$(idx_start)_$(idx_end).jld2")
            save(
                out_path,
                "spectra",
                spectra,
                "indices",
                indices,
                "L",
                L,
                "τ_idx",
                τ_idx,
                "λ",
                λ,
            )
            println("Spectra saved to: $out_path")
        end

        if !isempty(failed)
            println("\n=== Failed Tasks ===")
            for (i, (Lf, τf, λf, idxf, err)) in enumerate(failed)
                println("Failed $i: L=$Lf, τ_idx=$τf, λ=$λf, index=$idxf")
                println("  Error: $err")
            end
        end

    elseif mode == 2
        # -------------------------------------------------------------------
        # Mode 2: parallel MPS OBF scaling dimension computation
        # -------------------------------------------------------------------
        L         = parse(Int64, ARGS[2])
        τ_idx     = parse(Int64, ARGS[3])
        λ         = parse(Float64, ARGS[4])
        χ         = parse(Int64, ARGS[5])
        idx_start = parse(Int64, ARGS[6])
        idx_end   = parse(Int64, ARGS[7])
        tasks = [(L, τ_idx, λ, i, χ) for i in idx_start:idx_end]

        println("=== Parallel OBF Spectrum Computation (MPS) ===")
        println("L = $L, τ_idx = $τ_idx, τ = $(τlis[τ_idx]), λ = $λ, χ = $χ")
        println("Index range: $idx_start - $idx_end")
        println("Total tasks: $(length(tasks))")
        println("Number of workers: $(nworkers())")

        results = pmap(compute_tf_OBF_mps_task, tasks; batch_size = 1)

        spectra = Vector{Matrix{Float64}}()
        indices = Vector{Int}()
        for (Lr, τr, λr, ir, status, result) in results
            if status == :success
                push!(spectra, result)
                push!(indices, ir)
            end
        end

        failed = [
            (Lr, τr, λr, ir, err)
            for (Lr, τr, λr, ir, status, err) in results if status != :success
        ]
        success_count = length(spectra)

        println("\n=== Complete ===")
        println("Successes: $success_count")
        println("Failures: $(length(failed))")

        if success_count > 0
            out_dir = "exm/data/OBF/tf_spectrum_Born/L$(L)/gammaind$(τ_idx)/λ$(λ)/chi$(χ)"
            mkpath(out_dir)
            out_path = joinpath(out_dir, "spectrum_$(idx_start)_$(idx_end).jld2")
            save(
                out_path,
                "spectra",
                spectra,
                "indices",
                indices,
                "L",
                L,
                "τ_idx",
                τ_idx,
                "λ",
                λ,
                "χ",
                χ,
            )
            println("Spectra saved to: $out_path")
        end

        if !isempty(failed)
            println("\n=== Failed Tasks ===")
            for (i, (Lf, τf, λf, idxf, err)) in enumerate(failed)
                println("Failed $i: L=$Lf, τ_idx=$τf, λ=$λf, index=$idxf")
                println("  Error: $err")
            end
        end

    elseif mode == 3
        # -------------------------------------------------------------------
        # Mode 3: collect ensemble-averaged OBF spectra
        # -------------------------------------------------------------------
        L     = parse(Int64, ARGS[2])
        τ_idx = parse(Int64, ARGS[3])
        λ     = parse(Float64, ARGS[4])
        χ     = length(ARGS) >= 5 ? parse(Int64, ARGS[5]) : nothing

        println("=== Collecting Ensemble-Averaged OBF Spectrum ===")
        println("L = $L, τ_idx = $τ_idx, τ = $(τlis[τ_idx]), λ = $λ, χ = $χ")
        collect_scaling_dimension_Born(L, τ_idx, λ; χ = χ)
    end
end
