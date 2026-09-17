raw"""
    HybridConfig(measurement::MeasureConfig; p=1.0, θ=π, random_angles=false)
    HybridConfig(; p=1.0, θ=π, random_angles=false, measure_config_keywords...)

Configuration for a Fibonacci monitored circuit in which every location of the
usual staggered spacetime lattice independently contains a measurement with
probability `p`, or the unitary

```math
U(θ)=Π_1+e^{iθ}Π_\tau
```

with probability `1-p`. The embedded `MeasureConfig` controls the measurement
strength, time interval, RNG, and MPS truncation settings. In `:Born` mode the
gate locations and outcomes are sampled. In `:sample` mode a
`HybridGateSchedule` must be supplied for deterministic replay. If
`random_angles=true`, every unitary brick independently draws
`θ[x,t] ~ Uniform(0, 2π)`; otherwise all unitary bricks use the fixed `θ`.
"""
struct HybridConfig
    measurement::MeasureConfig
    p::Float64
    θ::Float64
    random_angles::Bool

    function HybridConfig(
        measurement::MeasureConfig;
        p::Real = 1.0,
        θ::Real = π,
        random_angles::Bool = false,
    )
        0 <= p <= 1 || throw(ArgumentError("p must satisfy 0 ≤ p ≤ 1, got $p"))
        isfinite(θ) || throw(ArgumentError("θ must be finite, got $θ"))
        return new(measurement, Float64(p), Float64(θ), random_angles)
    end
end

HybridConfig(; p::Real = 1.0, θ::Real = π, random_angles::Bool = false, kwargs...) =
    HybridConfig(
        MeasureConfig(; kwargs...);
        p = p,
        θ = θ,
        random_angles = random_angles,
    )

"""
    HybridGateSchedule(measurement_mask, outcomes, unitary_angles)

A replayable spacetime realization. `measurement_mask[t,x] == true` denotes a
measurement and `false` denotes `U(θ)`. `outcomes` is meaningful only where the
mask is true. `unitary_angles[t,x]` stores the angle of a unitary brick and is
`NaN` at measurement locations. All arrays use the same layer/column layout as
`Measurement_outcome_bulk.samples`. The two-argument constructor remains
available for legacy fixed-angle schedules; its angle entries are `NaN` and
are resolved from `HybridConfig.θ` during replay.
"""
struct HybridGateSchedule
    measurement_mask::BitMatrix
    outcomes::BitMatrix
    unitary_angles::Matrix{Float64}

    function HybridGateSchedule(
        measurement_mask::AbstractMatrix{Bool},
        outcomes::AbstractMatrix{Bool},
        unitary_angles::AbstractMatrix{<:Real},
    )
        size(measurement_mask) == size(outcomes) || throw(
            DimensionMismatch(
                "measurement_mask and outcomes must have the same size, got " *
                "$(size(measurement_mask)) and $(size(outcomes))",
            ),
        )
        size(measurement_mask) == size(unitary_angles) || throw(
            DimensionMismatch(
                "measurement_mask and unitary_angles must have the same size, got " *
                "$(size(measurement_mask)) and $(size(unitary_angles))",
            ),
        )
        return new(
            BitMatrix(measurement_mask),
            BitMatrix(outcomes),
            Matrix{Float64}(unitary_angles),
        )
    end
end

HybridGateSchedule(measurement_mask::AbstractMatrix{Bool}, outcomes::AbstractMatrix{Bool}) =
    HybridGateSchedule(
        measurement_mask,
        outcomes,
        fill(NaN, size(measurement_mask)),
    )

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

_draw_hybrid_angle(rng::MersenneTwister, config::HybridConfig) =
    config.random_angles ? 2π * rand(rng) : config.θ

function _scheduled_hybrid_angle(
    schedule::HybridGateSchedule,
    layer::Int,
    column::Int,
    config::HybridConfig,
)
    θ = schedule.unitary_angles[layer, column]
    return isfinite(θ) ? θ : config.θ
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
    angles = born ? fill(NaN, D, ncols) : copy(schedule.unitary_angles)
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
                    θ_event = born ?
                              _draw_hybrid_angle(mc.rng, config) :
                              _scheduled_hybrid_angle(schedule, glayer, col, config)
                    angles[glayer, col] = θ_event
                    _apply_fibonacci_unitary_exact!(
                        buffer,
                        basis,
                        measure_model,
                        current,
                        site,
                        θ_event,
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
        HybridGateSchedule(masks, outcomes, angles),
        free_energys,
        entropies,
    )
end

function _hybrid_sector_projector(
    model::AnyonModel{FibonacciAnyon},
    sector::Union{Nothing,Symbol},
)
    l = length(anyon_basis(model))
    sector === nothing && return Matrix{Float64}(I, l, l), nothing, l, nothing
    model.pbc || error("A topological charge sector requires periodic boundaries")
    sector ∈ (:trivial, :tau) || throw(
        ArgumentError("sector must be nothing, :trivial, or :tau, got $sector"),
    )
    ϕ = (1 + √5) / 2
    y = sector == :trivial ? ϕ : -inv(ϕ)
    y_other = sector == :trivial ? -inv(ϕ) : ϕ
    Y = topological_charge_operator(model)
    P = Matrix(Symmetric((Y - y_other * I(l)) / (y - y_other)))
    return P, Y, round(Int, real(tr(P))), y
end

raw"""
    hybrid_lyapunov_spectrum(model, config, schedule;
                             n_states=10, sector=:trivial)

Compute the finite-time Lyapunov spectrum of an exact hybrid transfer-matrix
product, including every measurement projector and every stored random unitary
angle. The initial tangent frame is formed exactly as in `lyapunov_spectrum`:

```julia
states = zeros(l, k)
for i in 1:k
    states[i, i] = 1
end
```

and then projected into the requested topological sector before the initial QR.
Use `sector=nothing` for the full Hilbert space, `:trivial` for `y=1`, or `:tau`
for `y=τ`. For each period,

```math
T_tQ_{t-1}=Q_tR_t,
\qquad
\lambda_a(t)=\frac1t\sum_{s=1}^t\log|(R_s)_{aa}|.
```

The returned `free_energy_spectrum` is `-lyapunov_exponents`. No per-period
sorting is performed, since sorting would mix Oseledec directions.
"""
function hybrid_lyapunov_spectrum(
    model::AnyonModel{FibonacciAnyon},
    config::HybridConfig,
    schedule::HybridGateSchedule;
    n_states::Int = 10,
    sector::Union{Nothing,Symbol} = :trivial,
)
    n_states >= 1 || throw(ArgumentError("n_states must be positive"))
    n_layers = layers_per_period(model)
    D, ncols = size(schedule.measurement_mask)
    D % n_layers == 0 || error(
        "schedule row count $D must be divisible by $n_layers layers per period",
    )
    ncols == _samples_per_layer(model) || error(
        "schedule spatial dimension must be $(_samples_per_layer(model)), got $ncols",
    )
    periods = D ÷ n_layers
    P, Y, sector_dimension, y = _hybrid_sector_projector(model, sector)
    n_states <= sector_dimension || error(
        "requested n_states=$n_states, but the selected sector dimension is " *
        "$sector_dimension",
    )

    l = size(P, 1)
    states = zeros(ComplexF64, l, n_states)
    for i in 1:n_states
        states[i, i] = 1
    end
    factor = qr(P * states)
    minimum(abs.(diag(factor.R)[1:n_states])) > 1e-12 || error(
        "the first $n_states basis vectors are linearly dependent after sector projection; " *
        "reduce n_states",
    )
    states = Matrix(factor.Q)[:, 1:n_states]
    local_log_stretches = zeros(Float64, n_states, periods)
    sector_leakage = sector === nothing ? Float64[] : zeros(Float64, periods)
    mc = config.measurement
    step_config = HybridConfig(
        MeasureConfig(
            τ = mc.τ,
            t₂ = 1,
            mode = :sample,
            enable_τ_eff = false,
        );
        p = config.p,
        θ = config.θ,
        random_angles = config.random_angles,
    )

    for step in 1:periods
        rows = (step - 1) * n_layers + 1:step * n_layers
        step_schedule = HybridGateSchedule(
            schedule.measurement_mask[rows, :],
            schedule.outcomes[rows, :],
            schedule.unitary_angles[rows, :],
        )
        propagated = similar(states)
        for column in 1:n_states
            propagated[:, column] = bulk_evolution(
                model,
                states[:, column],
                step_config,
                step_schedule;
                normalized = false,
            ).state
        end

        factor = qr(propagated)
        stretches = abs.(diag(factor.R)[1:n_states])
        minimum(stretches) > 0 || error("Lyapunov frame collapsed at period $step")
        local_log_stretches[:, step] = log.(stretches)
        states = Matrix(factor.Q)[:, 1:n_states]
        if sector !== nothing
            states = Matrix(qr(P * states).Q)[:, 1:n_states]
            sector_leakage[step] = norm(Y * states - y * states) / norm(states)
        end
    end

    cumulative = cumsum(local_log_stretches; dims = 2)
    lyapunov_exponents = cumulative ./ reshape(collect(1:periods), 1, :)
    return (
        local_log_stretches = local_log_stretches,
        lyapunov_exponents = lyapunov_exponents,
        free_energy_spectrum = -lyapunov_exponents,
        sector = sector,
        sector_dimension = sector_dimension,
        sector_leakage = sector_leakage,
        final_frame = states,
    )
end

"""Exact common-record evolution diagnostics for topological-charge learnability."""
struct HybridLearnabilityOutcome{T}
    state_y1::Vector{T}
    state_ytau::Vector{T}
    schedule::HybridGateSchedule
    sampled_sector::Symbol
    log_record_probability_y1::Vector{Float64}
    log_record_probability_ytau::Vector{Float64}
    log_likelihood_ratio::Vector{Float64}
    posterior_y1::Vector{Float64}
    record_fidelity_estimator::Vector{Float64}
    conditional_entropy_estimator::Vector{Float64}
    mutual_information_estimator::Vector{Float64}
    bayes_error_estimator::Vector{Float64}
    y_expectation_y1::Vector{Float64}
    y_expectation_ytau::Vector{Float64}
end

function _normalized_hybrid_sector_state(
    model::AnyonModel{FibonacciAnyon},
    state::Vector,
    expected_y::Float64,
    label::String,
)
    ψ = ComplexF64.(state)
    norm2 = real(dot(ψ, ψ))
    norm2 > 0 || error("$label initial state has zero norm")
    ψ ./= sqrt(norm2)
    Y = topological_charge_operator(model)
    y = real(dot(ψ, Y * ψ))
    isapprox(y, expected_y; atol = 1e-8, rtol = 1e-8) || error(
        "$label initial state is not in the requested topological sector: ⟨Y⟩=$y",
    )
    return ψ, Y
end

function _hybrid_branch_probability!(
    dest::Vector{ComplexF64},
    basis,
    measure_model,
    state::Vector{ComplexF64},
    site::Int,
    strength::Float64,
    outcome::Bool,
)
    fill!(dest, 0)
    _measuremap_impl!(dest, basis, measure_model, strength, state, site, outcome)
    probability = sum(abs2, dest)
    if probability > 0
        dest .*= inv(sqrt(probability))
    end
    return probability
end

function _log_likelihood_ratio(logp_y1::Float64, logp_ytau::Float64)
    isfinite(logp_y1) && isfinite(logp_ytau) && return logp_y1 - logp_ytau
    logp_y1 == -Inf && isfinite(logp_ytau) && return -Inf
    isfinite(logp_y1) && logp_ytau == -Inf && return Inf
    error("the observed record has zero probability in both topological sectors")
end

function _posterior_from_log_likelihood(ell::Float64)
    ell == Inf && return 1.0
    ell == -Inf && return 0.0
    ell >= 0 && return inv(1 + exp(-ell))
    z = exp(ell)
    return z / (1 + z)
end

function _record_fidelity_estimator(ell::Float64)
    isfinite(ell) || return 0.0
    a = abs(ell)
    return 2exp(-a / 2) / (1 + exp(-a))
end

_binary_entropy(p::Float64) =
    (p == 0 || p == 1) ? 0.0 : -p * log(p) - (1 - p) * log1p(-p)

raw"""
    hybrid_bayesian_evolution(
        model, state_y1, state_ytau, config::HybridConfig,
        [schedule]; sampled_sector=:mixture, sector_rng=MersenneTwister(0)
    )

Evolve the two topological-sector hypotheses through the same exclusive hybrid
circuit and the same observed measurement record. `state_y1` must have
`Y` eigenvalue `ϕ`, and `state_ytau` eigenvalue `-1/ϕ`. In `:Born` mode,
`sampled_sector=:mixture` first chooses either hypothesis with equal prior and
samples every observed outcome from that hypothesis using the Born rule. A
fixed truth can be requested with `:trivial` or `:tau`. In `:sample` mode the
supplied `HybridGateSchedule` is replayed and `sampled_sector` is ignored.

At every measurement event, both normalized hypothesis states are updated and
their cumulative record likelihoods are accumulated. At the end of period
`t`, the log-likelihood ratio is

```math
\ell_t=\log\frac{P(\mathbf m_{\leq t}\mid y=1,\mathcal C)}
                 {P(\mathbf m_{\leq t}\mid y=\tau,\mathcal C)}.
```

For equal priors, the returned per-trajectory estimators obey

```math
F_t=\mathbb E_{P_{\rm mix}}[\operatorname{sech}(\ell_t/2)],
\qquad
H(Y\mid\mathbf M_t)=\mathbb E_{P_{\rm mix}}[h_2(\sigma(\ell_t))].
```

Consequently `mutual_information_estimator` is
`log(2) - conditional_entropy_estimator`, and `bayes_error_estimator` is
`min(posterior_y1, 1-posterior_y1)`. This routine is exact-state only.
"""
function hybrid_bayesian_evolution(
    model::AnyonModel{FibonacciAnyon},
    state_y1::Vector,
    state_ytau::Vector,
    config::HybridConfig,
    schedule::Union{Nothing,HybridGateSchedule} = nothing;
    sampled_sector::Symbol = :mixture,
    sector_rng::AbstractRNG = MersenneTwister(0),
)
    Δt, (D, ncols) = _validate_hybrid_run(model, config, schedule)
    model.pbc || error("Topological-sector learnability requires periodic boundaries")
    mc = config.measurement
    born = mc.mode == :Born
    born && sampled_sector ∉ (:mixture, :trivial, :tau) && throw(
        ArgumentError("sampled_sector must be :mixture, :trivial, or :tau"),
    )
    truth = if born && sampled_sector == :mixture
        rand(sector_rng, Bool) ? :trivial : :tau
    elseif born
        sampled_sector
    else
        :replay
    end

    ϕ = (1 + √5) / 2
    current_y1, Y = _normalized_hybrid_sector_state(model, state_y1, ϕ, "y=1")
    current_ytau, _ =
        _normalized_hybrid_sector_state(model, state_ytau, -inv(ϕ), "y=τ")
    buffer_y1 = similar(current_y1)
    buffer_ytau = similar(current_ytau)
    probe = similar(current_y1)
    alive_y1 = true
    alive_ytau = true

    masks = born ? falses(D, ncols) : copy(schedule.measurement_mask)
    outcomes = born ? falses(D, ncols) : copy(schedule.outcomes)
    angles = born ? fill(NaN, D, ncols) : copy(schedule.unitary_angles)
    logp_y1 = 0.0
    logp_ytau = 0.0
    log_record_y1 = zeros(Float64, Δt)
    log_record_ytau = zeros(Float64, Δt)
    log_likelihood = zeros(Float64, Δt)
    posterior_y1 = zeros(Float64, Δt)
    fidelity_estimator = zeros(Float64, Δt)
    entropy_estimator = zeros(Float64, Δt)
    mutual_information_estimator = zeros(Float64, Δt)
    bayes_error_estimator = zeros(Float64, Δt)
    y_expectation_y1 = fill(NaN, Δt)
    y_expectation_ytau = fill(NaN, Δt)
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

            for (k, site) in enumerate(event_sites)
                col = cols[k]
                is_measurement = born ? _draw_hybrid_gate(mc.rng, config.p) : masks[glayer, col]
                masks[glayer, col] = is_measurement

                if is_measurement
                    outcome = if born
                        truth_state = truth == :trivial ? current_y1 : current_ytau
                        fill!(probe, 0)
                        _measuremap_impl!(
                            probe,
                            basis,
                            measure_model,
                            strength,
                            truth_state,
                            site,
                            false,
                        )
                        p0 = clamp(sum(abs2, probe), 0.0, 1.0)
                        rand(mc.rng) >= p0
                    else
                        outcomes[glayer, col]
                    end
                    outcomes[glayer, col] = outcome

                    if alive_y1
                        q_y1 = _hybrid_branch_probability!(
                            buffer_y1,
                            basis,
                            measure_model,
                            current_y1,
                            site,
                            strength,
                            outcome,
                        )
                        if q_y1 > 0
                            logp_y1 += log(q_y1)
                            current_y1, buffer_y1 = buffer_y1, current_y1
                        else
                            alive_y1 = false
                            logp_y1 = -Inf
                        end
                    end
                    if alive_ytau
                        q_ytau = _hybrid_branch_probability!(
                            buffer_ytau,
                            basis,
                            measure_model,
                            current_ytau,
                            site,
                            strength,
                            outcome,
                        )
                        if q_ytau > 0
                            logp_ytau += log(q_ytau)
                            current_ytau, buffer_ytau = buffer_ytau, current_ytau
                        else
                            alive_ytau = false
                            logp_ytau = -Inf
                        end
                    end
                    truth == :trivial && !alive_y1 &&
                        error("sampled a zero-probability record in the y=1 truth sector")
                    truth == :tau && !alive_ytau &&
                        error("sampled a zero-probability record in the y=τ truth sector")
                else
                    θ_event = born ?
                              _draw_hybrid_angle(mc.rng, config) :
                              _scheduled_hybrid_angle(schedule, glayer, col, config)
                    angles[glayer, col] = θ_event
                    if alive_y1
                        _apply_fibonacci_unitary_exact!(
                            buffer_y1,
                            basis,
                            measure_model,
                            current_y1,
                            site,
                            θ_event,
                        )
                        current_y1, buffer_y1 = buffer_y1, current_y1
                    end
                    if alive_ytau
                        _apply_fibonacci_unitary_exact!(
                            buffer_ytau,
                            basis,
                            measure_model,
                            current_ytau,
                            site,
                            θ_event,
                        )
                        current_ytau, buffer_ytau = buffer_ytau, current_ytau
                    end
                end
            end
        end

        ell = _log_likelihood_ratio(logp_y1, logp_ytau)
        posterior = _posterior_from_log_likelihood(ell)
        entropy = _binary_entropy(posterior)
        log_record_y1[period] = logp_y1
        log_record_ytau[period] = logp_ytau
        log_likelihood[period] = ell
        posterior_y1[period] = posterior
        fidelity_estimator[period] = _record_fidelity_estimator(ell)
        entropy_estimator[period] = entropy
        mutual_information_estimator[period] = log(2) - entropy
        bayes_error_estimator[period] = min(posterior, 1 - posterior)
        alive_y1 &&
            (y_expectation_y1[period] = real(dot(current_y1, Y * current_y1)))
        alive_ytau &&
            (y_expectation_ytau[period] = real(dot(current_ytau, Y * current_ytau)))
    end

    return HybridLearnabilityOutcome(
        current_y1,
        current_ytau,
        HybridGateSchedule(masks, outcomes, angles),
        truth,
        log_record_y1,
        log_record_ytau,
        log_likelihood,
        posterior_y1,
        fidelity_estimator,
        entropy_estimator,
        mutual_information_estimator,
        bayes_error_estimator,
        y_expectation_y1,
        y_expectation_ytau,
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
    angles = born ? fill(NaN, D, ncols) : copy(schedule.unitary_angles)
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
                    θ_event = born ?
                              _draw_hybrid_angle(mc.rng, config) :
                              _scheduled_hybrid_angle(schedule, glayer, col, config)
                    angles[glayer, col] = θ_event
                    if θ_event != 0
                        U = _fibonacci_unitary_operator_mps(measure_model, sites, site, θ_event)
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
        HybridGateSchedule(masks, outcomes, angles),
        free_energys,
        entropies,
    )
end
