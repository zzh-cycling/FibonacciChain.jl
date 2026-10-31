using FibonacciChain
using JLD
using JLD2
using LinearAlgebra
using Statistics

"""Equal-time C[x,y] = <O_x O_y†>, O_i = Πτ n_i Π1, with n_i the fusion-path bit."""
function projected_occupation_correlation(state, Y, occupations)
    norm_squared = sum(abs2, state)
    isfinite(norm_squared) && norm_squared > 0 || error("Invalid trajectory state")
    ϕ = (1 + sqrt(5)) / 2
    denominator = ϕ + inv(ϕ)
    tau_state = (ϕ .* state .- Y * state) ./ denominator
    # Columns are O_i†|ψ> = Π1 n_i Πτ|ψ>, so their Gram matrix has the requested order.
    frame = occupations .* tau_state
    frame = (Y * frame .+ inv(ϕ) .* frame) ./ denominator
    return (; correlation=frame' * frame / norm_squared,
        weight_tau=sum(abs2, tau_state) / norm_squared)
end

"""
Replay a saved monitored_dynamics.jl record once, measuring at complete periods.
Time zero is the original initial state. Only the original record's last layer
has half strength; intermediate time slices do not acquire a new endpoint layer.
"""
function replay_projected_occupation_correlation(model, sample, τ, times;
    sites=collect(1:model.N), Y=topological_charge_operator(model),
    initial_state=nothing,
)
    model isa AnyonModel{FibonacciAnyon} && model.pbc && iseven(model.N) ||
        throw(ArgumentError("Use an even periodic Fibonacci chain"))
    sample isa AbstractMatrix{Bool} || throw(ArgumentError("sample must be a Bool matrix"))
    size(sample, 2) == model.N ÷ 2 && iseven(size(sample, 1)) && size(sample, 1) > 0 ||
        throw(DimensionMismatch("Expected a (2 * periods) × (L ÷ 2) record"))
    periods = size(sample, 1) ÷ 2
    !isempty(times) && issorted(times) && allunique(times) &&
        all(t -> t isa Integer && 0 <= t <= periods, times) ||
        throw(ArgumentError("times must be sorted unique integer periods in 0:$periods"))
    !isempty(sites) && allunique(sites) && all(i -> i isa Integer && 1 <= i <= model.N, sites) ||
        throw(ArgumentError("sites must be distinct integers in 1:L"))
    τ >= 0 || throw(ArgumentError("τ must be nonnegative"))
    basis = anyon_basis(model)
    size(Y) == (length(basis), length(basis)) || throw(DimensionMismatch("Invalid Y matrix"))
    occupations = [Float64(path[model.N-site+1]) for path in basis, site in sites]
    state = initial_state === nothing ? [i == 1 ? 1.0 : 0.0 for i in eachindex(basis)] :
        copy(initial_state)
    length(state) == length(basis) && isfinite(norm(state)) && norm(state) > 0 ||
        throw(ArgumentError("Invalid initial state"))
    state ./= norm(state)
    correlations = zeros(ComplexF64, length(sites), length(sites), length(times))
    weight_tau = zeros(length(times))
    output_index = 1
    for period in 0:last(times)
        if period == times[output_index]
            result = projected_occupation_correlation(state, Y, occupations)
            correlations[:,:,output_index] = result.correlation
            weight_tau[output_index] = result.weight_tau
            output_index += 1
        end
        period == last(times) && break
        for layer in (2period+1):(2period+2)
            config = MeasureConfig(; τ=Float64(τ), t₂=1, mode=:sample,
                enable_τ_eff=layer == 2periods)
            state = boundary_evolution(model, state, config, BitVector(sample[layer,:]);
                layer_idx=layer).state
            all(isfinite, state) || error("Zero-probability or invalid record at layer $layer")
        end
    end
    return (; correlations, weight_tau, state)
end

"""Read original .jld samples, replay, and save per-seed C[x,y,time,seed] plus means/SEMs."""
function collect_projected_occupation_correlation(L, τ, times, sample_files, output)
    ispath(output) && error("Output already exists: $output")
    isempty(sample_files) && throw(ArgumentError("Provide at least one saved trajectory"))
    files = abspath.(sample_files)
    allunique(files) || error("Duplicate input files")
    model = AnyonModel(FibonacciAnyon(), L; pbc=true)
    # shortcut: the existing exact Y implementation is dense; large L needs an MPO implementation.
    Y = topological_charge_operator(model)
    sites = collect(1:L)
    seeds = Int[]
    source_periods = 0
    correlations = zeros(ComplexF64, L, L, length(times), length(files))
    weight_tau = zeros(length(times), length(files))
    for (index, path) in enumerate(files)
        data = JLD.load(path)
        sample = data["sample"]
        seed = data["seed"]
        seed isa Integer || error("Invalid seed in $path")
        seed in seeds && error("Duplicate seed $seed; do not combine overlapping ensembles")
        periods = size(sample, 1) ÷ 2
        if index == 1
            source_periods = periods
        else
            periods == source_periods || error("Source durations differ; endpoint conventions differ")
        end
        result = replay_projected_occupation_correlation(model, sample, τ, times; sites, Y)
        correlations[:,:,:,index] = result.correlations
        weight_tau[:,index] = result.weight_tau
        push!(seeds, seed)
        println("Replayed seed $seed ($index/$(length(files)))")
    end
    n = length(files)
    mean_correlation = dropdims(mean(correlations; dims=4); dims=4)
    stderr_real = n > 1 ? dropdims(std(real.(correlations); dims=4); dims=4) ./ sqrt(n) :
        fill(NaN, size(mean_correlation))
    stderr_imag = n > 1 ? dropdims(std(imag.(correlations); dims=4); dims=4) ./ sqrt(n) :
        fill(NaN, size(mean_correlation))
    mkpath(dirname(abspath(output)))
    JLD2.jldsave(output; L, τ, times=collect(times), sites, seeds, source_periods,
        source_sample_files=files, correlations, weight_tau, mean_correlation,
        stderr_real, stderr_imag, samples_num=n, time_unit="complete_two_layer_period",
        operator="O_i = Pi_tau n_i Pi_1; n_i = fusion-path bit 1 occupation",
        ordering="C[x,y,t,seed] = <psi(t)|O_x O_y^dagger|psi(t)>",
        initial_state="first_fusion_path_basis_vector", enable_τ_eff=true)
    return output
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 5 || error(
        "Usage: julia --project=. exm/Bulk_measure/projected_occupation_correlation.jl " *
        "L TAU_INDEX TIMES_CSV OUTPUT.jld2 SAMPLE1.jld [SAMPLE2.jld ...]",
    )
    include(joinpath(@__DIR__, "config.jl"))
    L = parse(Int, ARGS[1])
    τ_index = parse(Int, ARGS[2])
    times = parse.(Int, split(ARGS[3], ','))
    collect_projected_occupation_correlation(L, τlis_ext[τ_index], times, ARGS[5:end], ARGS[4])
end
