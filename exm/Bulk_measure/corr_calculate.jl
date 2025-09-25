using FibonacciChain
using JLD
using Statistics
using Random

function get_δtL(τ)
    table = Dict(
            atanh(0.1)  => (collect(25:35)),
            atanh(0.2)  => (collect(10:20)),
            atanh(0.3)  => (collect(5:15)),
            atanh(0.4)  => (collect(2:10)),)
    δtlis = get(table, τ, collect(1:8))
    return  δtlis
end

function organize(τ::Float64, sign::Int64=1)
    δtlis = get_δtL(τ)
    Llis = collect(8:4:20)
    scLlis = zeros(Float64, length(Llis))
    tcLlis = zeros(Float64, length(Llis), length(δtlis))
    seedlis = collect(1:100)
    for (i, L) in enumerate(Llis)
        D, _, _ = get_system_params(τ, L)
        D = div(D, 2L)
        seed = seedlis[i]
        spatial_corr = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/D$(D)_Samples$(seed).jld", "spatial_corr")
        scLlis[i] = spatial_corr
        for (j, δt) in enumerate(δtlis)
            temporal_corr = load("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(D)_Samples$(seed).jld",  "temporal_corr")
            tcLlis[i, j] = temporal_corr
        end
    end
    save("exm/data/Bulk_measure/spatial_temporal_corr/stc$(sign)_τ$(τ)_L$(Llis[1])$(Llis[end])_t0$(δtlis[end]).jld", "scLlis", scLlis, "tcLlis", tcLlis)
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

function compute_total(L::Int64, τ::Float64, index::Int64, D::Int64=35L, δt::Int64=2)
    for index in index:index+99
        @show "Computing index=$index"
        @time compute_ratio(L, τ, index, D, δt)
    end
end

function compute_ratio(L::Int64, τ::Float64, index::Int64, D::Int64=16L, δt::Int64=2)
    pbc = true
    sample = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L$L/τ$(τ)/D$(div(D,L))_Samples$(index).jld", "sample")

    initial_state = zeros(length(anyon_basis(L, pbc)))
    initial_state[1] = 1.0 # initial state is all zero state

    rng = MersenneTwister(index)

    statelis, Flis = generate_state(τ, initial_state, sample, enable_τ_eff=false)
    D = div(D, 2)
    ref_sample = zeros(Int, 2*(D+δt+D), length(2:2:L))
    view(ref_sample, 1:2D, :) .= view(sample, :, :)

    if δt == 0
        ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D, verbose=false, rng = rng) # to compute temporal correlation, add ref qubit at site L/2+1
        spatial = true
        temporal = false
        view(sample_free_energy, 1:2D) .= view(Flis, :)
    else
        ref2stlis, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D+δt, x₁ = L÷2+1, verbose=false, rng = rng) # to compute temporal correlation, add ref qubit at site L/2+1
        temporal = true
        spatial = false
        view(sample_free_energy, 1:2D) .= view(Flis, :)
    end
    spatial_corr, temporal_corr = ref_correlation(L, ref2stlis[end], spatial = spatial, temporal = temporal)
    sysrdm = reference_rdm(L, collect(1:div(L,2)), ref2stlis[end], traceref = false)
    S = ee(sysrdm)
    

    save("exm/data/Bulk_measure/spatial_temporal_corr_Born/L$(L)/τ$(τ)/dt$(δt)/D$(div(D,L))_Samples$(index).jld", "temporal_corr", temporal_corr, "spatial_corr", spatial_corr, "S", S, "sample_layer", sample_layer, "sample_free_energy", sample_free_energy)
    return temporal_corr, spatial_corr, S, sample, sample_free_energy
end


function corr_collect(L::Int64, τ::Float64, D::Int64=35L)
    samples_num = 10000
    temporal_corr_ensemble = Vector{Vector{Float64}}(undef, samples_num)
    spatial_corr_ensemble = Vector{Float64}(undef, samples_num)
     for i in 1:samples_num
        @show i
        temporal_corr_lis, spatial_corr = load("exm/data/Bulk_measure/temporal_corr/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(i).jld",  "temporal_corr_lis", "spatial_corr")
        temporal_corr_ensemble[i] = temporal_corr_lis
        spatial_corr_ensemble[i] = spatial_corr
    end

    average_temporal_corr = mean(temporal_corr_ensemble)
    average_spatial_corr = mean(spatial_corr_ensemble)
    temporal_corr_stderr = std(temporal_corr_ensemble) ./ sqrt(samples_num)
    spatial_corr_stderr = std(spatial_corr_ensemble) / sqrt(samples_num)

    save("exm/data/Bulk_measure/temporal_spatial_corr_L$(L)_τ$(τ)_D$(div(D,L)).jld", "average_temporal_corr", average_temporal_corr, 
    "temporal_corr_stderr", temporal_corr_stderr, 
    "average_spatial_corr", average_spatial_corr, 
    "spatial_corr_stderr", spatial_corr_stderr)
end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 
seed_interval_lis = collect(1:100:2000)

if length(ARGS) == 0
    println("No arguments provided.")
else
    L = parse(Int64, ARGS[1])
    inds = parse(Int64, ARGS[2])
    δt = parse(Int, ARGS[3])
    seedinds = parse(Int, ARGS[4])
    τ = τlis[inds]
    D, _, _ = get_system_params(τ, L)
    seed = seed_interval_lis[seedinds]
    println("Computed spatial_temporal_corr_varyingt for L=$L, τ=$τ, D=$D, δt=$δt, seedlis=$(seed):$(seed+99)")
    # compute_ratio(L, τ, seed, D, δt)
    compute_total(L, τ, seed, D, δt)
    # corr_collect(L, τ)
end