using FibonacciChain
using JLD
using Statistics
using Random

function samples_generate(L::Int64, τ::Float64, index::Int64, seed::Int64, D::Int64=120L)
    rng = MersenneTwister(seed)
    
    st = zeros(length(anyon_basis(L)))
    st[1] = 1.0
    
    @time sample_measured_states, sample, sample_free_energy = bulk_measure(L, τ, st, div(D,2), rng, true) 
    
    halfchain_EE_tlis = [ee(anyon_rdm(L, collect(1:div(L,2)), j)) for j in sample_measured_states]
    final_state = sample_measured_states[end]
    final_EElis = anyon_eelis(L, final_state)

    
    save("./exm/data/Bulk_measure/Observable_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "halfchain_EE_tlis", halfchain_EE_tlis, "final_EElis ", final_EElis, "seed", seed, "sample_free_energy", sample_free_energy)

    save("exm/data/Bulk_measure/Samples_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "sample", sample, "sample_free_energy", sample_free_energy, "seed", seed)
    # return sample_measured_states, samples, sample_free_energy
end


function samples_collect(L::Int64, τ::Float64, D::Int64=120L)
    samples_num = 2000
    ensemble = Vector{Matrix{Int}}(undef, samples_num)
    ensemble_free_energy = Vector{Vector{Float64}}(undef, samples_num)
    ensemble_seed = Vector{Int64}(undef, samples_num)
     for i in 1:samples_num
        @show i
        @time sample, sample_free_energy, seed = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(i).jld", "sample", "sample_free_energy", "seed")
        ensemble[i] = sample
        ensemble_free_energy[i] = sample_free_energy
        ensemble_seed[i] = seed
    end

    save("exm/data/Bulk_measure/monitored_dynamics_ensemble_L$(L)_τ$(τ)_D$(div(D,L)).jld", "ensemble", ensemble, "ensemble_free_energy", ensemble_free_energy, "ensemble_seed", ensemble_seed)
end


function Observable_collect(L::Int64, τ::Float64, D::Int64=120L)
    samples_num = 2000
    ensemble_free_energy = Vector{Vector{Float64}}(undef, samples_num)
    ensemble_seed = Vector{Int64}(undef, samples_num)
    ensemble_EE_dynamics= zeros(samples_num, D) 
    ensemble_final_EElis = zeros(samples_num, L-1)

    for i in 1:samples_num
        @show i
        halfchain_EE_tlis, final_EElis, seed, sample_free_energy = load("./exm/data/Bulk_measure/Observable_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(i).jld", "halfchain_EE_tlis", "final_EElis ", "seed",  "sample_free_energy")

        ensemble_EE_dynamics[i, :] = halfchain_EE_tlis
        ensemble_final_EElis[i, :] = final_EElis
        ensemble_seed[i] = seed
        ensemble_free_energy[i] = sample_free_energy
    end

    bulk_meanEElis = mean(ensemble_final_EElis, dims=1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims=1)[:]
    ensemble_stderr_EElis = (std(ensemble_final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims=1) ./ sqrt(samples_num))[:]

    
    save("exm/data/Bulk_measure/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld", "average_EE_tlis", average_EE_tlis, "stderr_EE_tlis", stderr_EE_tlis, "bulk_meanEElis", bulk_meanEElis, "ensemble_stderr_EElis",ensemble_stderr_EElis, "ensemble_free_energy", ensemble_free_energy, "ensemble_seed", ensemble_seed)
end

function monitored_dynamics(L::Int64, τ::Float64, D::Int64=120L)
    st=zeros(length(anyon_basis(L)))
    st[1] = 1.0
    bulk_meanEElis=zeros(L-1)
    
    samples_num = 2000

    ensemble_EE_dynamics= zeros(samples_num, D) 
    ensemble_stderr_EElis = zeros(L-1)
    final_EElis = zeros(samples_num, L-1)

    all_FE_tlis = zeros(samples_num, D)
    final_FElis = zeros(samples_num)
    for i in 1:samples_num
        @show i
        sample_measured_states, sample, sample_free_energy = bulk_measure(L, τ, st, D) 
        ensemble_EE_dynamics[i, :] = [ee(anyon_rdm(L, collect(1:div(L,2)), j)) for j in sample_measured_states]
        final_state = sample_measured_states[end]
        final_EElis[i, :] = anyon_eelis(L, final_state)

        all_FE_tlis[i, :] = sample_free_energy
        final_FElis[i] = sample_free_energy[end]
    end

    bulk_meanEElis = mean(final_EElis, dims=1)[:]
    average_EE_tlis = mean(ensemble_EE_dynamics, dims=1)[:]
    ensemble_stderr_EElis = (std(final_EElis, dims=1) ./ sqrt(samples_num))[:]
    stderr_EE_tlis = (std(ensemble_EE_dynamics, dims=1) ./ sqrt(samples_num))[:]
    
    average_FE = mean(final_FElis)
    stderr_FE = std(final_FElis) / sqrt(samples_num)
    average_FE_tlis = mean(all_FE_tlis, dims=1)[:]
    stderr_FE_tlis = (std(all_FE_tlis, dims=1) ./ sqrt(samples_num))[:]

    
    return average_EE_tlis, stderr_EE_tlis, bulk_meanEElis, ensemble_stderr_EElis, 
           average_FE_tlis, stderr_FE_tlis, average_FE, stderr_FE
end

function process_data(L::Int64, D::Int64=25L, τ::Float64=log(1+ √2))
    # timewindow = 8L:35L-10
    timewindow = 5L:D-5  # Adjusted time window for averaging
    load_data_path = "exm/data/Bulk_measure/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld"
    data = load(load_data_path)
    
    average_EE_tlis, stderr_EE_tlis = data["average_EE_tlis"], data["stderr_EE_tlis"]
    bulk_meanEElis, ensemble_stderr_EElis = data["bulk_meanEElis"], data["ensemble_stderr_EElis"]
    ensemble_free_energy, ensemble_seed = data["ensemble_free_energy"], data["ensemble_seed"]
    
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
    
    # Check if there are duplicates in the ensemble_seed
    has_duplicates = check_duplicates(ensemble_seed)
    
    temp = hcat(ensemble_free_energy...)
    time_average_free_energy = mean(temp[timewindow, :], dims=1) 
    bulk_FE = mean(time_average_free_energy)
    bulk_FE_stderr = std(time_average_free_energy) / sqrt(size(temp, 2))
    time_FEstderr = (std(temp, dims=2) ./ sqrt(size(temp, 2)))[:]
    time_FElis = mean(temp, dims=2)[:]
    
    save("exm/data/Bulk_measure/simple/monitored_EE_FEdynamics_L$(L)_τ$(τ)_D$(div(D,L)).jld", 
        "average_EE_tlis", average_EE_tlis, 
        "stderr_EE_tlis", stderr_EE_tlis, 
        "bulk_meanEElis", bulk_meanEElis, 
        "ensemble_stderr_EElis", ensemble_stderr_EElis, 
        "time_average_free_energy", time_average_free_energy, 
        "bulk_FE", bulk_FE,
        "bulk_FE_stderr", bulk_FE_stderr, 
        "time_FEstderr", time_FEstderr, 
        "time_FElis", time_FElis, 
        "ensemble_seed", ensemble_seed)
end

function get_system_params(τ, L)
    cfg = Dict(
        atanh(0.1)  => (2500L, 1000, 1500L),
        atanh(0.2)  => (500L,  100, 250L),
        atanh(0.3)  => (120L,  48, 100L),
        atanh(0.4)  => (100L,  40, 80L),
        atanh(0.5)  => (80L,   32, 40L),
        atanh(0.6)  => (45L,   20, 30L),
        log(1 + √2) => (35L,   14, 20L),
        atanh(0.8)  => (25L,   10, 10L),
        atanh(0.9)  => (8L,    4, 4L),
        atanh(0.95) => (8L,    4, 4L),
        atanh(0.999)=> (5L,    2, 2L),
    )
    D, step, start = get(cfg, τ, (5L, 2, 2L))
    inds = collect(1:step:D)
    avg_range = start:D-5
    return D, inds, avg_range
end


γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)


if length(ARGS) == 0
    println("No arguments provided.")
else
    L = parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    τ = τlis[inds]
    index = parse(Int64, ARGS[3])
    seed = parse(Int64, ARGS[4])
    D, _, _ = get_system_params(τ, L)
    println("Computed Born Sample dynamics for L=$L, τ=$τ, D=$D, index=$index, seed=$seed")
    samples_generate(L, τ, index, seed, D)
    # Observable_collect(L, τ)
    # samples_collect(L, τ)
end
