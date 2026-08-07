using FibonacciChain
using JLD
using Statistics
using Random

include(joinpath(@__DIR__, "config.jl"))

function samples_generate(L::Int64, τ_idx::Int, index::Int64, D::Int64 = 8L)
    rng = MersenneTwister(index)

    model = ising_model(L)
    st = ones(length(anyon_basis(model)))
    st ./= norm(st)
    τ = τlis[τ_idx]
    t = div(D, 2)
    config = MeasureConfig(τ = τ, mode = :Born, t₂ = t, rng = rng)
    @time mo = bulk_evolution(model, st, config)
    sample, sample_free_energy = mo.samples, mo.free_energys

    halfchain_EE_tlis = mo.entanglement_entropys
    final_state = mo.state
    final_EElis = anyon_eelis(model, final_state)

    out_dir = "exm/data/Bulk_measure/Ising/L$(L)/gammaind$(τ_idx)"
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
    # return sample_measured_states, samples, sample_free_energy
end

function samples_collect_process_data(L::Int64, τind::Int64)
        D, _, _ = get_cfg_params_Born(τind, L)
        t = div(D, 2L)
        dir_path = "exm/data/Bulk_measure/Ising/L$(L)/gammaind$(τind)/"
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
            "exm/data/Bulk_measure/Ising/L$(L)/EE_FEdynamics_L$(L)_gamma$(τind)_t$(t).jld2",
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


if length(ARGS) == 0
    println("No arguments provided.")
else
    N=parse(Int64, ARGS[1])
    idx=parse(Int64, ARGS[2])
    index=parse(Int64, ARGS[3])
    println("Received argument: $N, $idx")
    samples_generate(N, idx, index)
    # samples_collect(N, τ)
end
