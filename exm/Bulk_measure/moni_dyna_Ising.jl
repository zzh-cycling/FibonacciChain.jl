using FibonacciChain
using JLD
using Statistics
using Random

function total_samples_generate(
    L::Int64,
    τ::Float64,
    index::Int64,
    seed::Int64,
    D::Int64 = 20L,
)
    for i = 0:99
        @show (index+i), (seed+i)
        samples_generate(L, τ, index+i, seed + i, D)
    end
end

function samples_generate(L::Int64, τ::Float64, index::Int64, seed::Int64, D::Int64 = 20L)
    rng = MersenneTwister(seed)

    model = AnyonModel(IsingAnyon(), L; pbc = true, anyon_type = :X)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0

    config = MeasureConfig(τ = τ, mode = :Born, t₂ = div(D, 2), rng = rng)
    mo = bulk_evolution(model, st, config)
    @time sample, sample_free_energy = mo.samples, mo.free_energys

    halfchain_EE_tlis = mo.entanglement_entropys
    final_state = mo.state
    final_EElis = anyon_eelis(model, final_state)


    save(
        "./exm/data/Bulk_measure/Ising/Observable_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld",
        "halfchain_EE_tlis",
        halfchain_EE_tlis,
        "final_EElis ",
        final_EElis,
        "seed",
        seed,
        "sample_free_energy",
        sample_free_energy,
    )

    save(
        "exm/data/Bulk_measure/Ising/Samples_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld",
        "sample",
        sample,
        "sample_free_energy",
        sample_free_energy,
        "seed",
        seed,
    )
    # return sample_measured_states, samples, sample_free_energy
end

function sample_continue_calculate(
    L::Int64,
    τ::Float64,
    index::Int64,
    seed::Int64,
    D::Int64 = 20L,
    additional_layers::Int64 = 15L,
)
    rng = MersenneTwister(seed)

    sample, sample_free_energy, seed = load(
        "exm/data/Bulk_measure/Ising/Samples_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld",
        "sample",
        "sample_free_energy",
        "seed",
    )
    model = AnyonModel(IsingAnyon(), L; pbc = true, anyon_type = :X)

    st = zeros(length(anyon_basis(model)));
    st[1] = 1.0
    pre_config = MeasureConfig(τ = τ, mode = :sample, t₂ = div(D, 2))
    pre_mo = bulk_evolution(model, st, pre_config, sample)
    final_st = pre_mo.state

    re_config =
        MeasureConfig(τ = τ, mode = :Born, t₂ = div(additional_layers, 2), rng = rng)
    after_mo = bulk_evolution(model, final_st, re_config)
    sample, sample_free_energy = after_mo.samples, after_mo.free_energys
    halfchain_EE_tlis = after_mo.entanglement_entropys
    final_state = after_mo.state
    final_EElis = anyon_eelis(model, final_state)


    save(
        "./exm/data/Bulk_measure/Ising/Observable_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D+additional_layers,L))_Samples$(index).jld",
        "halfchain_EE_tlis",
        halfchain_EE_tlis,
        "final_EElis ",
        final_EElis,
        "seed",
        seed,
        "sample_free_energy",
        sample_free_energy,
    )

    save(
        "exm/data/Bulk_measure/Ising/Samples_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D+additional_layers,L))_Samples$(index).jld",
        "sample",
        sample,
        "sample_free_energy",
        sample_free_energy,
        "seed",
        seed,
    )
end

function samples_collect(L::Int64, τ::Float64, D::Int64 = 20L)
    samples_num = 10000
    ensemble = Vector{Matrix{Int}}(undef, samples_num)
    ensemble_free_energy = Vector{Vector{Float64}}(undef, samples_num)
    ensemble_seed = Vector{Int64}(undef, samples_num)
    for i = 1:samples_num
        @show i
        @time sample, sample_free_energy, seed = load(
            "exm/data/Bulk_measure/Ising/Samples_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(i).jld",
            "sample",
            "sample_free_energy",
            "seed",
        )
        ensemble[i] = sample
        ensemble_free_energy[i] = sample_free_energy
        ensemble_seed[i] = seed
    end

    save(
        "exm/data/Bulk_measure/Ising/monitored_dynamics_ensemble_L$(L)_τ$(τ)_D$(div(D,L)).jld",
        "ensemble",
        ensemble,
        "ensemble_free_energy",
        ensemble_free_energy,
        "ensemble_seed",
        ensemble_seed,
    )
end


function Observable_collect(L::Int64, τ::Float64, D::Int64 = 20L)
    samples_num = 10000
    ensemble_free_energy = Vector{Vector{Float64}}(undef, samples_num)
    ensemble_seed = Vector{Int64}(undef, samples_num)
    ensemble_EE_dynamics = zeros(samples_num, D)
    ensemble_final_EElis = zeros(samples_num, L-1)

    for i = 1:samples_num
        @show i
        halfchain_EE_tlis, final_EElis, seed, sample_free_energy = load(
            "./exm/data/Bulk_measure/Ising/Observable_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(i).jld",
            "halfchain_EE_tlis",
            "final_EElis ",
            "seed",
            "sample_free_energy",
        )

        ensemble_EE_dynamics[i, :] = halfchain_EE_tlis
        ensemble_final_EElis[i, :] = final_EElis
        ensemble_seed[i] = seed
        ensemble_free_energy[i] = sample_free_energy
    end

    bulk_meanEElis = mean(ensemble_final_EElis, dims = 1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims = 1)[:]
    bulk_stderr_EElis = (std(ensemble_final_EElis, dims = 1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims = 1) ./ sqrt(samples_num))[:]


    save(
        "exm/data/Bulk_measure/Ising/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld",
        "average_EE_tlis",
        average_EE_tlis,
        "stderr_EE_tlis",
        stderr_EE_tlis,
        "bulk_meanEElis",
        bulk_meanEElis,
        "bulk_stderr_EElis",
        bulk_stderr_EElis,
        "ensemble_free_energy",
        ensemble_free_energy,
        "ensemble_seed",
        ensemble_seed,
    )
end

function monitored_dynamics(L::Int64, τ::Float64, D::Int64 = 20L, window = 5L:(D-5))
    model = AnyonModel(IsingAnyon(), L; pbc = true, anyon_type = :X)
    st=zeros(length(anyon_basis(model)))
    st[1] = 1.0
    bulk_meanEElis=zeros(L-1)

    samples_num = 10000

    ensemble_EE_dynamics = zeros(samples_num, D)
    bulk_stderr_EElis = zeros(L-1)
    final_EElis = zeros(samples_num, L-1)

    all_FE_tlis = zeros(samples_num, D)
    final_FElis = zeros(samples_num)
    seed_lis = zeros(Int64, samples_num)
    for i = 1:samples_num
        @show i
        config =
            MeasureConfig(τ = τ, mode = :Born, t₂ = div(D, 2), rng = MersenneTwister(i))
        @time mo = bulk_evolution(model, st, config)
        sample, sample_free_energy = mo.samples, mo.free_energys
        #
        ensemble_EE_dynamics[i, :] = mo.entanglement_entropys
        final_state = mo.state
        final_EElis[i, :] = anyon_eelis(model, final_state)

        all_FE_tlis[i, :] = sample_free_energy
        final_FElis[i] = sample_free_energy[end]
        seed_lis[i] = i
    end

    bulk_meanEElis = mean(final_EElis, dims = 1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims = 1)[:]
    bulk_stderr_EElis = (std(final_EElis, dims = 1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims = 1) ./ sqrt(samples_num))[:]

    time_average_free_energy = mean(all_FE_tlis[:, window], dims = 2)
    bulk_FE = mean(time_average_free_energy)
    bulk_FE_stderr = std(time_average_free_energy) / sqrt(size(all_FE_tlis, 1))
    time_FEstderr = (std(all_FE_tlis, dims = 1) ./ sqrt(size(all_FE_tlis, 1)))[:]
    time_FElis = mean(all_FE_tlis, dims = 1)[:]


    data_path = "exm/data/Bulk_measure/Ising/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld"

    save(
        data_path,
        "average_EE_tlis",
        average_EE_tlis,
        "stderr_EE_tlis",
        stderr_EE_tlis,
        "bulk_meanEElis",
        bulk_meanEElis,
        "bulk_stderr_EElis",
        bulk_stderr_EElis,
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
        seed_lis,
    )

    return average_EE_tlis,
    stderr_EE_tlis,
    bulk_meanEElis,
    bulk_stderr_EElis,
    time_FElis,
    time_FEstderr,
    bulk_FE,
    bulk_FE_stderr,
    seed_lis
end

τ = log(1 + sqrt(2)) # golden ratio
# N = 12
# average_EE_tlis, stderr_EE_tlis, bulk_meanEElis, bulk_stderr_EElis,
# time_FElis, time_FEstderr, bulk_FE, bulk_FE_stderr, seed_lis = monitored_dynamics(N, τ)
# fig = plot(collect(1:20*N), average_EE_tlis, yerror=stderr_EE_tlis, xlabel = L"t", ylabel=L"S_{vN}", label=false)

# plot(collect(1:20*N), time_FElis, yerror=time_FEstderr, xlabel = L"t", ylabel=L"S", label=false)

# cent, fig = fitCCEntEntScal(bulk_meanEElis, mincut=1, pbc=true)
# display(fig)
if length(ARGS) == 0
    println("No arguments provided.")
else
    N=parse(Int64, ARGS[1])
    index=parse(Int64, ARGS[2])
    seed=parse(Int64, ARGS[3])
    println("Received argument: $N, $index")
    total_samples_generate(N, τ, index, seed)
    # samples_generate(N, τ, index, seed)
    # Observable_collect(N, τ)
    # samples_collect(N, τ)
end
