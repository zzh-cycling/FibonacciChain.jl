raw"""
    HybridConfig(measurement::MeasureConfig; p=1.0, θ=π)
    HybridConfig(; p=1.0, θ=π, measure_config_keywords...)

Configuration for a Fibonacci monitored circuit in which every location of the
usual staggered spacetime lattice independently contains a measurement with
probability `p`, or the unitary

```math
U(z)=Π_1+e^{iθ}Π_\tau
```

with probability `1-p`. The embedded `MeasureConfig` controls the measurement
strength, time interval, RNG, and MPS truncation settings. In `:Born` mode the
gate locations and outcomes are sampled. In `:sample` mode a
`HybridGateSchedule` must be supplied for deterministic replay.
"""
struct HybridConfig
    measurement::MeasureConfig
    p::Float64
    θ::Float64

    function HybridConfig(measurement::MeasureConfig; p::Real = 1.0, θ::Real = π)
        0 <= p <= 1 || throw(ArgumentError("p must satisfy 0 ≤ p ≤ 1, got $p"))
        isfinite(θ) || throw(ArgumentError("θ must be finite, got $θ"))
        return new(measurement, Float64(p), Float64(θ))
    end
end

HybridConfig(; p::Real = 1.0, θ::Real = π, kwargs...) =
    HybridConfig(MeasureConfig(; kwargs...); p = p, θ = θ)

"""
    HybridGateSchedule(measurement_mask, outcomes)

A replayable spacetime realization. `measurement_mask[t,x] == true` denotes a
measurement and `false` denotes `U(θ)`. `outcomes` is meaningful only where the
mask is true. Both matrices use the same layer/column layout as
`Measurement_outcome_bulk.samples`.
"""
struct HybridGateSchedule
    measurement_mask::BitMatrix
    outcomes::BitMatrix

    function HybridGateSchedule(measurement_mask::AbstractMatrix{Bool}, outcomes::AbstractMatrix{Bool})
        size(measurement_mask) == size(outcomes) || throw(
            DimensionMismatch(
                "measurement_mask and outcomes must have the same size, got " *
                "$(size(measurement_mask)) and $(size(outcomes))",
            ),
        )
        return new(BitMatrix(measurement_mask), BitMatrix(outcomes))
    end
end

"""Result of exact-state or MPS hybrid monitored evolution."""
struct HybridMeasurementOutcome{ST}
    state::ST
    schedule::HybridGateSchedule
    free_energys::Vector{Float32}
    entanglement_entropys::Vector{Float32}
end

function _validate_hybrid_run(
    model::AnyonModel{FibonacciAnyon},
    config::HybridConfig,
    schedule::Union{Nothing,HybridGateSchedule},
)
    mc = config.measurement
    mc.mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")
    mc.mode == :Born && schedule !== nothing &&
        error("do not provide a HybridGateSchedule in mode=:Born")
    mc.mode == :sample && schedule === nothing &&
        error("mode=:sample requires a HybridGateSchedule")
    Δt = mc.t₂ - mc.t₁ + 1
    Δt >= 0 || error("t₂ must be >= t₁")
    dims = (Δt * layers_per_period(model), _samples_per_layer(model))
    schedule !== nothing && size(schedule.measurement_mask) != dims &&
        error("schedule size should be $dims, got $(size(schedule.measurement_mask))")
    return Δt, dims
end

# Apply Πτ without constructing a dense matrix. In the existing convention the
# projective `sign=false` Fibonacci Kraus operator is exactly Πτ.
function _apply_fibonacci_unitary_exact!(
    dest::Vector{T},
    basis,
    measure_model::AnyonModel{FibonacciAnyon},
    state::Vector{T},
    site::Int,
    θ::Float64,
) where {T<:Complex}
    θ == 0 && return copyto!(dest, state)
    fill!(dest, zero(T))
    _measuremap_impl!(dest, basis, measure_model, Inf, state, site, false)
    α = cis(θ) - 1
    @inbounds @simd for j in eachindex(dest, state)
        dest[j] = state[j] + α * dest[j]
    end
    return dest
end

function _hybrid_measure_event_exact!(
    dest::Vector{T},
    basis,
    measure_model,
    state::Vector{T},
    site::Int,
    τ::Float64,
    outcome::Bool;
    normalized::Bool,
) where {T}
    fill!(dest, zero(T))
    _measuremap_impl!(dest, basis, measure_model, τ, state, site, outcome)
    if normalized
        prob = sum(abs2, dest)
        prob > 0 || error("selected a zero-probability measurement outcome")
        dest .*= inv(sqrt(prob))
        return -log(prob)
    end
    return 0.0
end

function _draw_hybrid_gate(rng::MersenneTwister, p::Float64)
    p == 1 && return true
    p == 0 && return false
    return rand(rng) < p
end

"""
    bulk_evolution(model, state, config::HybridConfig, [schedule]; normalized=true)

Run the hybrid Fibonacci circuit using an exact state vector. The staggered
layer geometry is identical to the measurement-only evolution. Gate type is an
independent Bernoulli variable at every valid spacetime location.
"""
function bulk_evolution(
    model::AnyonModel{FibonacciAnyon},
    state::Vector,
    config::HybridConfig,
    schedule::Union{Nothing,HybridGateSchedule} = nothing;
    normalized::Bool = true,
)
    Δt, (D, ncols) = _validate_hybrid_run(model, config, schedule)
    mc = config.measurement
    born = mc.mode == :Born
    born && !normalized && error("mode=:Born requires normalized=true")
    masks = born ? falses(D, ncols) : copy(schedule.measurement_mask)
    outcomes = born ? falses(D, ncols) : copy(schedule.outcomes)
    free_energys = zeros(Float32, D)
    entropies = zeros(Float32, Δt)
    current = ComplexF64.(state)
    buffer = similar(current)
    n_layers = layers_per_period(model)

    for period in 1:Δt
        for layer in 1:n_layers
            glayer = (period - 1) * n_layers + layer
            τ_current =
                (period == Δt && layer == n_layers && mc.enable_τ_eff) ? mc.τ / 2 : mc.τ
            event_sites, measure_model, strength =
                _obtain_measurement_config(model, glayer, τ_current)
            basis = anyon_basis(measure_model)
            cols = _get_sample_column_indices(model, glayer)
            F = 0.0

            for (k, site) in enumerate(event_sites)
                col = cols[k]
                is_measurement = born ? _draw_hybrid_gate(mc.rng, config.p) : masks[glayer, col]
                masks[glayer, col] = is_measurement
                if is_measurement
                    outcome = if born
                        fill!(buffer, 0)
                        _measuremap_impl!(buffer, basis, measure_model, strength, current, site, false)
                        p0 = clamp(sum(abs2, buffer), 0.0, 1.0)
                        rand(mc.rng) >= p0
                    else
                        outcomes[glayer, col]
                    end
                    outcomes[glayer, col] = outcome
                    F += _hybrid_measure_event_exact!(
                        buffer,
                        basis,
                        measure_model,
                        current,
                        site,
                        strength,
                        outcome;
                        normalized = normalized,
                    )
                else
                    _apply_fibonacci_unitary_exact!(
                        buffer,
                        basis,
                        measure_model,
                        current,
                        site,
                        config.θ,
                    )
                end
                current, buffer = buffer, current
            end
            free_energys[glayer] = Float32(F)
        end
        entropies[period] = Float32(ee(anyon_rdm(model, collect(1:div(model.N, 2)), current)))
    end

    return HybridMeasurementOutcome(
        current,
        HybridGateSchedule(masks, outcomes),
        free_energys,
        entropies,
    )
end

function _fibonacci_unitary_operator_mps(measure_model, sites, site, θ)
    Pτ = _measurement_operator_mps_application(measure_model, sites, site, Inf, false)
    P1 = _measurement_operator_mps_application(measure_model, sites, site, Inf, true)
    return P1 + cis(θ) * Pτ
end

function _apply_hybrid_mps_operator(ψ::MPS, operator, config::MeasureConfig, truncate_now::Bool)
    if truncate_now
        ψ = apply(
            operator,
            ψ;
            cutoff = config.cutoff,
            mindim = config.mindim,
            maxdim = config.maxdim,
        )
    else
        ψ = apply(operator, ψ; cutoff = 0.0)
    end
    return ψ
end

"""
    bulk_evolution(model, sites, state::MPS, config::HybridConfig, [schedule]; normalized=true)

MPS implementation of the replayable hybrid Fibonacci circuit. `cutoff`,
`mindim`, `maxdim`, and `truncate_every_events` are inherited from the embedded
`MeasureConfig`.
"""
function bulk_evolution(
    model::AnyonModel{FibonacciAnyon},
    sites::Vector{<:Index},
    state::MPS,
    config::HybridConfig,
    schedule::Union{Nothing,HybridGateSchedule} = nothing;
    normalized::Bool = true,
)
    Δt, (D, ncols) = _validate_hybrid_run(model, config, schedule)
    mc = config.measurement
    mc.truncate_every_events >= 1 || error("truncate_every_events must be >= 1")
    born = mc.mode == :Born
    born && !normalized && error("mode=:Born requires normalized=true")
    masks = born ? falses(D, ncols) : copy(schedule.measurement_mask)
    outcomes = born ? falses(D, ncols) : copy(schedule.outcomes)
    free_energys = zeros(Float32, D)
    entropies = zeros(Float32, Δt)
    current = copy(state)
    n_layers = layers_per_period(model)
    constraint_projector = if mc.enforce_fibonacci_constraint
        model.pbc || error(
            "enforce_fibonacci_constraint currently requires periodic boundaries",
        )
        fibonacci_constraint_projector_mpo(sites)
    else
        nothing
    end

    for period in 1:Δt
        for layer in 1:n_layers
            glayer = (period - 1) * n_layers + layer
            τ_current =
                (period == Δt && layer == n_layers && mc.enable_τ_eff) ? mc.τ / 2 : mc.τ
            event_sites, measure_model, strength =
                _obtain_measurement_config(model, glayer, τ_current)
            cols = _get_sample_column_indices(model, glayer)
            F = 0.0
            nevents = length(event_sites)

            for (k, site) in enumerate(event_sites)
                col = cols[k]
                is_measurement = born ? _draw_hybrid_gate(mc.rng, config.p) : masks[glayer, col]
                masks[glayer, col] = is_measurement
                truncate_now = k % mc.truncate_every_events == 0 || k == nevents

                if is_measurement
                    outcome = outcomes[glayer, col]
                    if born
                        M0 = _measurement_operator_mps_application(measure_model, sites, site, strength, false)
                        ψ0 = _apply_hybrid_mps_operator(current, M0, mc, mc.truncate_every_events == 1)
                        p0 = clamp(real(inner(ψ0, ψ0)), 0.0, 1.0)
                        outcome = rand(mc.rng) >= p0
                        outcomes[glayer, col] = outcome
                    end
                    M = _measurement_operator_mps_application(measure_model, sites, site, strength, outcome)
                    next = _apply_hybrid_mps_operator(current, M, mc, truncate_now)
                    if normalized
                        prob = real(inner(next, next))
                        prob > 0 || error("selected a zero-probability measurement outcome")
                        F += -log(prob)
                        normalize!(next)
                    end
                    current = next
                else
                    if config.θ != 0
                        U = _fibonacci_unitary_operator_mps(measure_model, sites, site, config.θ)
                        current = _apply_hybrid_mps_operator(current, U, mc, truncate_now)
                        normalized && normalize!(current)
                    end
                end

                if mc.truncate_every_events > 1 && truncate_now
                    current = truncate(
                        current;
                        cutoff = mc.cutoff,
                        mindim = mc.mindim,
                        maxdim = mc.maxdim,
                    )
                    normalized && normalize!(current)
                end
            end
            if constraint_projector !== nothing
                current = _project_fibonacci_constraint(
                    constraint_projector,
                    current;
                    cutoff = mc.cutoff,
                    mindim = mc.mindim,
                    maxdim = mc.maxdim,
                )
            end
            free_energys[glayer] = Float32(F)
        end
        entropies[period] = Float32(ee_mps(current, div(length(sites), 2)))
    end

    return HybridMeasurementOutcome(
        current,
        HybridGateSchedule(masks, outcomes),
        free_energys,
        entropies,
    )
end
