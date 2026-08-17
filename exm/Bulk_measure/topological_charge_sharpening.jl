using Distributed
using FibonacciChain
using JLD2
using LinearAlgebra
using Random
using Statistics

# Launch with `julia -p N` when parallel Born trajectories are needed.
const BULK_MEASURE_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BULK_MEASURE_CONFIG)

@everywhere begin
    using FibonacciChain
    using JLD2
    using LinearAlgebra
    using Random
    using Statistics

    function charge_sharpening_data_dir(
        mode::Symbol,
        L::Integer,
        τ_idx::Integer,
        t::Integer;
        sign::Bool = true,
    )
        mode ∈ (:Born, :sample) || error("mode must be one of :Born, :sample")
        branch = mode == :Born ? "Born" : sign ? "sample_AFM" : "sample_FM"
        return joinpath(
            "exm",
            "data",
            "Bulk_measure",
            "topological_charge_sharpening",
            branch,
            "L$(L)",
            "gammaind$(τ_idx)",
            "t$(t)",
        )
    end

    """
        simulate_topological_charge_sharpening(L, τ_idx, t; mode=:Born,
                                               seed=1, sign=true, samples=nothing)

    Simulate one topological-charge trajectory. In `:Born` mode, `seed` controls
    the measurement outcomes. In `:sample` mode, pass a fixed `samples` record or
    use `sign` to post-select every outcome to AFM (`true`) or FM (`false`).
    """
    function simulate_topological_charge_sharpening(
        L::Integer,
        τ_idx::Integer,
        t::Integer;
        mode::Symbol = :Born,
        seed::Integer = 1,
        sign::Bool = true,
        samples::Union{Nothing,BitMatrix} = nothing,
    )
        iseven(L) || error("L must be even, got $L")
        1 <= τ_idx <= length(τlis) || error("τ_idx must be in 1:$(length(τlis))")
        t >= 1 || error("t must be positive, got $t")
        mode ∈ (:Born, :sample) || error("mode must be one of :Born, :sample")

        τ = τlis[τ_idx]
        model = fib_model(L)
        initial_state = zeros(length(anyon_basis(model)))
        initial_state[1] = 1

        # `topological_charge_sharpening` records the reference entropy after
        # each full period. Compute the t = 0 value from the same Y-sector
        # decomposition before evolving the state.
        ϕ = (1 + √5) / 2
        yτ = -inv(ϕ)
        Ystate = topological_charge_operator(model) * initial_state
        state_y₁ = (Ystate .- yτ .* initial_state) ./ (ϕ - yτ)
        p_y₁ = real(dot(state_y₁, state_y₁))
        p_yτ = 1.0 - p_y₁
        initial_charge_entropy =
            -sum(p -> iszero(p) ? 0.0 : p * log(p), (p_y₁, p_yτ))

        if mode == :Born
            config = MeasureConfig(
                τ = τ,
                mode = :Born,
                t₂ = t,
                rng = MersenneTwister(seed),
                enable_τ_eff = false,
            )
            outcome = topological_charge_sharpening(model, initial_state, config)
        else
            if isnothing(samples)
                samples = BitMatrix(fill(sign, 2t, div(L, 2)))
            end
            config = MeasureConfig(
                τ = τ,
                mode = :sample,
                t₂ = t,
                enable_τ_eff = false,
            )
            outcome =
                topological_charge_sharpening(model, initial_state, config, samples)
        end

        reference_density_matrix = reference_rdm(model, [1], outcome.state)
        charge_entropy = Float64.(outcome.entanglement_entropys)

        return (
            L = Int(L),
            τ_idx = Int(τ_idx),
            τ = τ,
            t = Int(t),
            mode = mode,
            seed = Int(seed),
            sign = sign,
            initial_charge_entropy = initial_charge_entropy,
            charge_entropy = charge_entropy,
            reference_density_matrix = reference_density_matrix,
            samples = outcome.samples,
            sample_free_energy = Float64.(outcome.free_energys),
        )
    end

    function save_charge_sharpening_result(result)
        out_dir = charge_sharpening_data_dir(
            result.mode,
            result.L,
            result.τ_idx,
            result.t;
            sign = result.sign,
        )
        mkpath(out_dir)

        filename = if result.mode == :Born
            "trajectory_seed$(result.seed).jld2"
        else
            branch = result.sign ? "AFM" : "FM"
            "post_selection_$(branch).jld2"
        end
        output_path = joinpath(out_dir, filename)

        JLD2.jldsave(
            output_path;
            L = result.L,
            τ_idx = result.τ_idx,
            τ = result.τ,
            t = result.t,
            mode = String(result.mode),
            seed = result.seed,
            sign = result.sign,
            initial_charge_entropy = result.initial_charge_entropy,
            charge_entropy = result.charge_entropy,
            reference_density_matrix = result.reference_density_matrix,
            samples = result.samples,
            sample_free_energy = result.sample_free_energy,
        )
        return output_path
    end

    function run_charge_sharpening_task(task)
        L, τ_idx, t, seed = task
        try
            result = simulate_topological_charge_sharpening(
                L,
                τ_idx,
                t;
                mode = :Born,
                seed = seed,
            )
            output_path = save_charge_sharpening_result(result)
            return (L, τ_idx, t, seed, :success, output_path)
        catch e
            return (L, τ_idx, t, seed, :failed, e)
        end
    end

    function collect_charge_sharpening_born(
        L::Integer,
        τ_idx::Integer,
        t::Integer,
    )
        data_dir = charge_sharpening_data_dir(
            :Born,
            L,
            τ_idx,
            t,
        )
        files = sort(
            filter(
                path -> startswith(path, "trajectory_seed") && endswith(path, ".jld2"),
                readdir(data_dir),
            ),
        )
        isempty(files) && error("No Born trajectory files found in $data_dir")

        samples_num = length(files)
        ensemble_entropy = zeros(Float64, samples_num, t)
        seeds = zeros(Int, samples_num)
        initial_entropy = NaN

        for (i, filename) in enumerate(files)
            data = JLD2.load(joinpath(data_dir, filename))
            entropy = Float64.(data["charge_entropy"])
            length(entropy) == t || error(
                "charge-entropy length mismatch in $filename: expected $t, got $(length(entropy))",
            )
            ensemble_entropy[i, :] = entropy
            seeds[i] = data["seed"]
            if i == 1
                initial_entropy = data["initial_charge_entropy"]
            else
                isapprox(data["initial_charge_entropy"], initial_entropy; atol = 1e-12) ||
                    error("inconsistent initial charge entropy in $filename")
            end
        end
        length(unique(seeds)) == samples_num || error("duplicate trajectory seeds found")

        average_charge_entropy = vec(mean(ensemble_entropy; dims = 1))
        stderr_charge_entropy =
            vec(std(ensemble_entropy; dims = 1, corrected = false)) ./ √samples_num

        output_path = joinpath(data_dir, "charge_entropy_ensemble.jld2")
        JLD2.jldsave(
            output_path;
            L = Int(L),
            τ_idx = Int(τ_idx),
            τ = τlis[τ_idx],
            t = Int(t),
            initial_charge_entropy = initial_entropy,
            average_charge_entropy = average_charge_entropy,
            stderr_charge_entropy = stderr_charge_entropy,
            seeds = seeds,
            samples_num = samples_num,
        )
        return output_path, average_charge_entropy, stderr_charge_entropy
    end
end

function print_usage()
    println("Usage:")
    println(
        "  julia -p N exm/Bulk_measure/topological_charge_sharpening.jl born L τ_idx t seed_start seed_end",
    )
    println(
        "  julia exm/Bulk_measure/topological_charge_sharpening.jl sample L τ_idx t sign",
    )
    println(
        "  julia exm/Bulk_measure/topological_charge_sharpening.jl collect L τ_idx t",
    )
    println("Here sign=1 post-selects the AFM outcome and sign=0 the FM outcome.")
end

if isempty(ARGS)
    print_usage()
else
    action = lowercase(ARGS[1])

    if action == "born"
        length(ARGS) >= 6 || error("born mode requires L τ_idx t seed_start seed_end")
        L = parse(Int, ARGS[2])
        τ_idx = parse(Int, ARGS[3])
        t = parse(Int, ARGS[4])
        seed_start = parse(Int, ARGS[5])
        seed_end = parse(Int, ARGS[6])
        tasks = [(L, τ_idx, t, seed) for seed in seed_start:seed_end]

        println("=== Topological Charge Sharpening: Born Ensemble ===")
        println("L=$L, τ_idx=$τ_idx, t=$t, trajectories=$(length(tasks))")
        println("workers=$(nworkers())")
        results = pmap(run_charge_sharpening_task, tasks; batch_size = 1)
        failed = filter(result -> result[5] != :success, results)
        println("Successes: $(length(results) - length(failed))")
        println("Failures: $(length(failed))")
        for result in failed
            println("Failed seed $(result[4]): $(result[6])")
        end
    elseif action == "sample"
        length(ARGS) >= 5 || error("sample mode requires L τ_idx t sign")
        L = parse(Int, ARGS[2])
        τ_idx = parse(Int, ARGS[3])
        t = parse(Int, ARGS[4])
        sign = parse(Int, ARGS[5]) == 1

        result = simulate_topological_charge_sharpening(
            L,
            τ_idx,
            t;
            mode = :sample,
            sign = sign,
        )
        output_path = save_charge_sharpening_result(result)
        println("Saved fixed post-selection trajectory to $output_path")
    elseif action == "collect"
        length(ARGS) >= 4 || error("collect mode requires L τ_idx t")
        L = parse(Int, ARGS[2])
        τ_idx = parse(Int, ARGS[3])
        t = parse(Int, ARGS[4])
        output_path, _, _ = collect_charge_sharpening_born(
            L,
            τ_idx,
            t;
        )
        println("Saved Born ensemble statistics to $output_path")
    else
        print_usage()
        error("Unknown action: $(ARGS[1])")
    end
end
