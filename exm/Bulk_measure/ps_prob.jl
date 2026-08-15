using Distributed
using Statistics
using FibonacciChain
using JLD
using Random
using Printf

println("requested workers: ", nworkers())
println("total procs:       ", nprocs())
@everywhere println("host: ", gethostname(), "  pid: ", getpid())

const BULK_MEASURE_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BULK_MEASURE_CONFIG)

@everywhere begin
    using FibonacciChain
    using JLD
    using Random
    using Statistics
    include("../FitEntEntScal.jl")

    binary_distribution(p, rng) = rand(rng) < p ? 1 : 0

    function ps_prob_data_dir(L::Integer, τind::Integer, prob::Real)
        return joinpath(
            "exm/data/Bulk_measure/ps_prob_evolution",
            "prob$(prob)",
            "L$(L)",
            "gamma$(τind)",
        )
    end

    function ps_prob_sample_path(
        L::Integer,
        τind::Integer,
        D::Integer,
        prob::Real,
        seed::Integer,
    )
        return joinpath(
            ps_prob_data_dir(L, τind, prob),
            "L$(L)_t$(div(D,2L))_gamma$(τind)_prob$(prob)_sample$(seed).jld",
        )
    end

    function ps_prob_evolution(params)
        L, τind, prob, seed = params
        return ps_prob_evolution(L, τind, prob, seed)
    end

    function ps_prob_evolution(
        L::Integer,
        τind::Integer,
        prob::Real,
        seed::Integer,
    )
        D, _, _ = get_cfg_params_Born(τind, L)
        τ = τlis[τind]
        model = fib_model(L)
        initial_state = zeros(length(anyon_basis(model)))
        initial_state[1] = 1.0
        gate_num = div(D*L, 2)

        rng = MersenneTwister(seed)
        sample = BitMatrix(
            reshape([binary_distribution(prob, rng) for _ = 1:gate_num], D, div(L, 2)),
        )
        config = MeasureConfig(τ = τ, mode = :sample, t₂ = div(D, 2))
        mo = bulk_evolution(model, initial_state, config, sample)
        ee = Float64.(anyon_eelis(model, mo.state))
        ee_tlis = mo.entanglement_entropys
        sample_free_energy = Float32.(mo.free_energys)

        output_path = ps_prob_sample_path(L, τind, D, prob, seed)
        mkpath(dirname(output_path))
        save(
            output_path,
            "seed",
            seed,
            "ee",
            ee,
            "ee_tlis",
            ee_tlis,
            "sample_free_energy",
            sample_free_energy,
        )
    end
    # function ps_prob_evolution_Ising(params)
    # Try to reproduce the outcome of `Entanglement Transition in the Projective Transverse Field Ising Model`, but fail, due to it's not measurement at everywhere.
    function process_ps_prob_evolution(
        L::Integer,
        τind::Integer,
        prob::Real;
        samplelis = 1:10000,
    )
        D, _, avg_range = get_cfg_params_Born(τind, L)
        samplelis = collect(samplelis)
        isempty(samplelis) && throw(ArgumentError("samplelis must not be empty"))

        samples_num = length(samplelis)
        seedlis = zeros(Int64, samples_num)
        ensemble_ee = zeros(Float64, samples_num, L-1)
        ensemble_free_energy = Vector{Vector{Float32}}(undef, samples_num)
        ensemble_ee_tlis = zeros(Float64, samples_num, div(D, 2))

        for (i, sample) in enumerate(samplelis)
            sample_path = ps_prob_sample_path(L, τind, D, prob, sample)
            ee, ee_tlis, sample_free_energy, seed = load(
                sample_path,
                "ee",
                "ee_tlis",
                "sample_free_energy",
                "seed",
            )
            length(ee) == L-1 || error(
                "EE length mismatch in $sample_path: expected $(L-1), found $(length(ee))",
            )
            ensemble_ee[i, :] = ee
            if length(ee_tlis) == div(D, 2)
                ensemble_ee_tlis[i, :] = ee_tlis
            elseif length(ee_tlis) == D
                ensemble_ee_tlis[i, :] = ee_tlis[2:2:D]
            else
                error(
                    "EE dynamics length mismatch in $sample_path: expected $(div(D,2)) or $D, found $(length(ee_tlis))",
                )
            end
            ensemble_free_energy[i] = Float32.(sample_free_energy)
            seedlis[i] = seed
        end

        length(unique(seedlis)) == samples_num || error("duplicate seeds found in $seedlis")
        free_energy_lengths = unique(length.(ensemble_free_energy))
        length(free_energy_lengths) == 1 || error(
            "inconsistent free-energy lengths: $free_energy_lengths",
        )
        last(avg_range) <= only(free_energy_lengths) || error(
            "free-energy averaging range ends at $(last(avg_range)), but trajectories have length $(only(free_energy_lengths))",
        )

        average_ee = mean(ensemble_ee, dims = 1)[:]
        stderr_ee = (std(ensemble_ee, dims = 1) ./ sqrt(samples_num))[:]
        average_EE_tlis = mean(ensemble_ee_tlis, dims = 1)[:]
        stderr_EE_tlis = (std(ensemble_ee_tlis, dims = 1) ./ sqrt(samples_num))[:]

        temp = hcat(ensemble_free_energy...)
        time_average_free_energy = mean(temp[avg_range, :], dims = 1)
        bulk_FE = mean(time_average_free_energy)
        bulk_FE_stderr = std(time_average_free_energy) / sqrt(samples_num)
        time_FElis = mean(temp, dims = 2)[:]
        time_FEstderr = (std(temp, dims = 2) ./ sqrt(samples_num))[:]

        output_path = joinpath(
            "exm/data/Bulk_measure/ps_prob_evolution",
            "prob$(prob)",
            "L$(L)",
            "L$(L)_t$(div(D,2L))_gamma$(τind)_prob$(prob)_processed.jld",
        )
        save(
            output_path,
            "prob",
            prob,
            "average_ee",
            average_ee,
            "stderr_ee",
            stderr_ee,
            "average_EE_tlis",
            average_EE_tlis,
            "stderr_EE_tlis",
            stderr_EE_tlis,
            "time_average_free_energy",
            time_average_free_energy,
            "bulk_FE",
            bulk_FE,
            "bulk_FE_stderr",
            bulk_FE_stderr,
            "time_FElis",
            time_FElis,
            "time_FEstderr",
            time_FEstderr,
            "seedlis",
            seedlis,
            "avg_range",
            collect(avg_range),
        )

        return (
            seedlis = seedlis,
            bulk_FE = bulk_FE,
            bulk_FE_stderr = bulk_FE_stderr,
            average_EE_tlis = average_EE_tlis,
            stderr_EE_tlis = stderr_EE_tlis,
            time_FElis = time_FElis,
            time_FEstderr = time_FEstderr,
        )
    end

    function process_data(
        prob::Real;
        Llis = collect(8:2:20),
        τind_lis = [1, 3, 4, 7, 9, 10, 12],
    )
        summary_paths = String[]

        for τind in τind_lis
            τ = τlis[τind]
            cent_Llis = zeros(length(Llis))
            cent_stderrlis = zeros(length(Llis))
            bulk_FE_Llis = zeros(length(Llis))
            bulk_FE_stderrlis = zeros(length(Llis))

            for (id, L) in enumerate(Llis)
                D, _, _ = get_cfg_params_Born(τind, L)
                @show (L, τ, prob)
                cent, bulk_FE, bulk_FE_stderr = load(
                    joinpath(
                        "exm/data/Bulk_measure/ps_prob_evolution",
                        "prob$(prob)",
                        "L$(L)",
                        "L$(L)_t$(div(D,2L))_gamma$(τind)_prob$(prob)_processed.jld",
                    ),
                    "cent",
                    "bulk_FE",
                    "bulk_FE_stderr",
                )
                cent_Llis[id] = cent[1]
                cent_stderrlis[id] = cent[2]
                bulk_FE_Llis[id] = bulk_FE
                bulk_FE_stderrlis[id] = bulk_FE_stderr
            end

            summary_path = joinpath(
                "exm/data/Bulk_measure/ps_prob_evolution",
                "prob$(prob)",
                "cent_FE_L$(first(Llis))$(last(Llis))_gamma$(τind).jld",
            )
            save(
                summary_path,
                "prob",
                prob,
                "τind",
                τind,
                "Llis",
                collect(Llis),
                "cent_Llis",
                cent_Llis,
                "cent_stderrlis",
                cent_stderrlis,
                "bulk_FE_Llis",
                bulk_FE_Llis,
                "bulk_FE_stderrlis",
                bulk_FE_stderrlis,
            )
            push!(summary_paths, summary_path)
        end

        return summary_paths
    end

    function process_data(problis::AbstractVector; kwargs...)
        return [process_data(prob; kwargs...) for prob in problis]
    end
end


if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia -p N ps_prob.jl L_start L_end τ_idx prob seed_start seed_end")
else
    length(ARGS) == 6 || error("expected 6 arguments; received $(length(ARGS))")
    L_start = parse(Int64, ARGS[1])
    L_end = parse(Int64, ARGS[2])
    τind = parse(Int64, ARGS[3])
    prob = parse(Float64, ARGS[4])
    seed_start = parse(Int64, ARGS[5])
    seed_end = parse(Int64, ARGS[6])

    τ = τlis[τind]
    Llis = collect(L_start:2:L_end)
    seeds = collect(seed_start:seed_end)

    # Generate all (L, τind, prob, seed)
    tasks = [(L, τind, prob, seed) for L in Llis for seed in seeds]
    println(
        "Running $(length(tasks)) tasks: L=$Llis, τ=$τ, prob=$prob, seeds=$seed_start:$seed_end on $(nprocs()) workers",
    )

    results = @time pmap(tasks) do params
        try
            ps_prob_evolution(params)
            return (params, :success, nothing)
        catch e
            return (params, :failed, e)
        end
    end

    # statistics on results
    succeeded = filter(r -> r[2] == :success, results)
    failed = filter(r -> r[2] == :failed, results)
    println("\n=== Statistics ===")
    println("Success: $(length(succeeded)) / $(length(tasks))")
    println("Failed: $(length(failed)) / $(length(tasks))")

    if !isempty(failed)
        max_errors_to_print = min(length(failed), 20)
        println("\nFirst $max_errors_to_print failed tasks:")
        for (params, _, err) in failed[1:max_errors_to_print]
            print("  params=$params: ")
            showerror(stdout, err)
            println()
        end
        if length(failed) > max_errors_to_print
            println("  ... $(length(failed)-max_errors_to_print) additional failures omitted")
        end
    end
end
