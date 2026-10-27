using FibonacciChain
using LinearAlgebra
using Random
using Statistics
using JLD2

"""
    reference_correlation_sample(model, initial_state; τ, seed, warmup, separation,
                                 final_periods, x1=1, x2=..., hybrid=nothing, record=nothing)

Exact two-reference copy protocol. Times count complete circuit periods. Insert
R1 after `warmup`, R2 after `warmup + separation`, then evolve `final_periods`.
For temporal correlations set x1=x2; spatial correlations use separation=0.
All reference branches share one outcome sampled from their joint Born weight.
`record` replays the returned mask/outcomes/angles without consuming the RNG.
No half-strength final layer is used. References are passive after insertion.
"""
function reference_correlation_sample(model, initial_state; τ::Real, seed::Int,
    warmup::Int, separation::Int, final_periods::Int, x1::Int = 1,
    x2::Int = separation == 0 ? model.N ÷ 2 + 1 : x1,
    hybrid::Union{Nothing,HybridConfig} = nothing, record = nothing,
)
    model.pbc && iseven(model.N) && model.N >= 4 ||
        throw(ArgumentError("Use an even periodic chain with L >= 4"))
    min(warmup, separation, final_periods, seed) >= 0 ||
        throw(ArgumentError("Times and seed must be nonnegative"))
    τ >= 0 || throw(ArgumentError("τ must be nonnegative"))
    1 <= x1 <= model.N && 1 <= x2 <= model.N || throw(ArgumentError("Invalid sites"))
    hybrid === nothing || model.basis isa FibonacciAnyon ||
        throw(ArgumentError("Hybrid unitaries require Fibonacci anyons"))
    dim = length(anyon_basis(model))
    length(initial_state) == dim || throw(DimensionMismatch("Expected a system-only state"))
    norm(initial_state) > 0 || throw(ArgumentError("Initial state must be nonzero"))
    current = reshape(normalize(ComplexF64.(initial_state)), dim, 1)
    buffer = similar(current)
    rng = MersenneTwister(seed)
    n_layers = FibonacciChain.layers_per_period(model)
    total = warmup + separation + final_periods
    dims = (n_layers * total, FibonacciChain._samples_per_layer(model))
    masks = falses(dims)
    active = falses(dims)
    outcomes = falses(dims)
    angles = fill(NaN, dims)
    free_energys = zeros(Float64, dims[1])
    if record !== nothing
        all(size(getproperty(record, key)) == dims for key in
            (:measurement_mask, :outcomes, :unitary_angles)) ||
            throw(DimensionMismatch("Record dimensions must be $dims"))
    end
    # Cache one period of local maps; OBF has 14 layers, Fibonacci has 2.
    layers = [FibonacciChain._obtain_measurement_config(model, layer, Float64(τ))
              for layer in 1:n_layers]
    for time in 0:total
        if time == warmup
            current = reshape(add_reference_qubits(model, vec(current), x1), dim, :)
            buffer = similar(current)
        end
        if time == warmup + separation
            current = reshape(add_reference_qubits(model, vec(current), x2), dim, :)
            buffer = similar(current)
        end
        time == total && break
        for layer in 1:n_layers
            row = time * n_layers + layer
            sites, measure_model, strength = layers[layer]
            basis = anyon_basis(measure_model)
            cols = FibonacciChain._get_sample_column_indices(model, layer)
            for (site, col) in zip(sites, cols)
                active[row, col] = true
                measured = record !== nothing ? record.measurement_mask[row, col] :
                    hybrid === nothing ? true : FibonacciChain._draw_hybrid_gate(rng, hybrid.p)
                masks[row, col] = measured
                if measured
                    FibonacciChain._measuremap_frame_impl!(
                        buffer, basis, measure_model, strength, current, site, false)
                    p0 = clamp(sum(abs2, buffer), 0.0, 1.0)
                    sign = record === nothing ? rand(rng) >= p0 : record.outcomes[row, col]
                    if sign
                        FibonacciChain._measuremap_frame_impl!(
                            buffer, basis, measure_model, strength, current, site, true)
                    end
                    probability = sum(abs2, buffer)
                    probability > 0 || error("Zero-probability outcome at layer $row, site $site")
                    outcomes[row, col] = sign
                    free_energys[row] -= log(probability)
                    buffer ./= sqrt(probability)
                else
                    hybrid === nothing && error("Measurement-only record contains an active unitary")
                    θ = record === nothing ?
                        (hybrid.random_angles ? 2π * rand(rng) : hybrid.θ) :
                        record.unitary_angles[row, col]
                    isfinite(θ) || error("Missing unitary angle at layer $row, site $site")
                    angles[row, col] = θ
                    FibonacciChain._measuremap_frame_impl!(
                        buffer, basis, measure_model, Inf, current, site, false)
                    buffer .= current .+ (cis(θ) - 1) .* buffer
                end
                current, buffer = buffer, current
            end
        end
    end
    state = vec(current)
    s1 = ee(reference_rdm(model, [1], state))
    s2 = ee(reference_rdm(model, [2], state))
    s12 = ee(reference_rdm(model, [1, 2], state))
    return (; L = model.N, τ = Float64(τ), seed, warmup, separation, final_periods,
        x1, x2, layers_per_period = n_layers, time_unit = "complete_period",
        insertion = "copy", enable_τ_eff = false, correlation = s1 + s2 - s12,
        reference_entropies = [s1, s2, s12], state,
        active_mask = active, measurement_mask = masks, outcomes,
        unitary_angles = angles, free_energys)
end

"""Bracketed linear interpolation of C_temporal(t*) = C_spatial. No extrapolation."""
function reference_alpha_crossing(L, separations, spatial, temporal)
    spatial > 0 && all(isfinite, temporal) || return (; alpha = NaN, crossing = NaN)
    roots = Float64[]
    for i in eachindex(separations)
        temporal[i] == spatial && push!(roots, Float64(separations[i]))
        if i < length(separations) && (temporal[i] - spatial) * (temporal[i+1] - spatial) < 0
            push!(roots, separations[i] + (separations[i+1] - separations[i]) *
                (spatial - temporal[i]) / (temporal[i+1] - temporal[i]))
        end
    end
    length(roots) == 1 && only(roots) > 0 || return (; alpha = NaN, crossing = NaN)
    crossing = only(roots)
    return (; alpha = asinh(1.0) * L / (π * crossing), crossing)
end

"""
    reference_correlation_summary(L, separations, spatial, temporal; ...)

Rows are paired seeds, columns are positive temporal separations. Average MI
before matching. Whole-row bootstrap retains shared spatial-denominator and
cross-separation covariance. Alpha assumes cylinder scaling and separation L/2:
alpha = asinh(1) L / (π t*). NaN means absent/ambiguous crossing, not alpha=0.
"""
function reference_correlation_summary(L, separations, spatial, temporal;
    bootstrap_samples::Int = 1000, bootstrap_seed::Int = 1234,
)
    n = length(spatial)
    n >= 2 || throw(ArgumentError("At least two trajectories are required for SEM"))
    size(temporal) == (n, length(separations)) || throw(DimensionMismatch("Invalid ensemble"))
    !isempty(separations) && all(>(0), separations) && issorted(separations) &&
        allunique(separations) || throw(ArgumentError("Use sorted unique positive separations"))
    bootstrap_samples >= 0 || throw(ArgumentError("bootstrap_samples must be nonnegative"))
    all(isfinite, spatial) && all(isfinite, temporal) || error("Nonfinite correlation samples")
    spatial_mean = mean(spatial)
    spatial_sem = std(spatial) / sqrt(n)
    temporal_mean = vec(mean(temporal; dims = 1))
    temporal_sem = vec(std(temporal; dims = 1)) / sqrt(n)
    estimate = reference_alpha_crossing(L, separations, spatial_mean, temporal_mean)
    rng = MersenneTwister(bootstrap_seed)
    bootstrap_alpha = fill(NaN, bootstrap_samples)
    for b in eachindex(bootstrap_alpha)
        rows = rand(rng, 1:n, n)
        bootstrap_alpha[b] = reference_alpha_crossing(L, separations,
            mean(spatial[rows]), vec(mean(temporal[rows, :]; dims = 1))).alpha
    end
    valid = filter(isfinite, bootstrap_alpha)
    # Do not quote a selected-subset interval if many resamples lack a crossing.
    reliable = isfinite(estimate.alpha) && length(valid) >= 2 &&
        length(valid) >= 0.95 * bootstrap_samples
    alpha_sem = reliable ? std(valid) : NaN
    alpha_ci = reliable ? quantile(valid, [0.025, 0.975]) : [NaN, NaN]
    return (; L, separations, spatial, temporal, spatial_mean, spatial_sem,
        temporal_mean, temporal_sem, estimate..., alpha_sem, alpha_ci,
        bootstrap_alpha, bootstrap_seed,
        bootstrap_valid = length(valid), bootstrap_samples,
        alpha_method = "bracketed_linear_crossing", time_unit = "complete_period")
end

"""Generate full-record trajectories and a paired-seed ensemble in a new directory."""
function run_reference_correlation_scan(generator, L, separations, seeds, output;
    bootstrap_samples::Int = 1000,
)
    seeds = collect(seeds)
    separations = collect(separations)
    length(seeds) >= 2 && allunique(seeds) || error("Use at least two distinct seeds")
    !isempty(separations) && all(>(0), separations) && issorted(separations) &&
        allunique(separations) || error("Use sorted unique positive temporal separations")
    ispath(output) && error("Output already exists: $output; use a new directory")
    mkpath(output)
    spatial = zeros(length(seeds))
    temporal = zeros(length(seeds), length(separations))
    for (i, seed) in enumerate(seeds)
        for (j, separation) in enumerate([0; separations])
            result = generator(separation, seed)
            result.L == L || error("Generator returned the wrong system size")
            result.x2 - result.x1 == (separation == 0 ? L ÷ 2 : 0) ||
                error("Alpha scan requires antipodal spatial and same-site temporal probes")
            j == 1 ? (spatial[i] = result.correlation) : (temporal[i, j-1] = result.correlation)
            # The state is reconstructible from metadata and record; avoid storing it.
            saved = (; (key => value for (key, value) in pairs(result) if key != :state)...)
            jldsave(joinpath(output, "dt$(separation)_seed$(seed).jld2"); saved...)
        end
        println("Completed seed $seed ($(i)/$(length(seeds)))")
    end
    summary = reference_correlation_summary(L, separations, spatial, temporal; bootstrap_samples)
    jldsave(joinpath(output, "averaged.jld2"); seeds, summary...)
    return summary
end
