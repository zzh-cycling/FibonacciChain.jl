"""
    ising_reference_state(state; reference = [1, 1] / sqrt(2))

Attach one reference qubit to an Ising state: `kron(reference, state)`.
The first half of the resulting vector is the reference-0 (A) component;
the second half is the reference-1 (B) component. Normalize the joint state.
The reference labels the two staggered representations of the fusion chain.
"""
function ising_reference_state(state::AbstractVector;
    reference::AbstractVector = [1.0, 1.0] / sqrt(2))
    length(reference) == 2 || throw(DimensionMismatch("reference must have length 2"))
    ispow2(length(state)) || throw(ArgumentError("state length must be a power of two"))
    joint = float.(kron(reference, state))
    n = norm(joint)
    isfinite(n) && n > 0 || throw(ArgumentError("state must have finite, nonzero norm"))
    return joint / n
end

function _check_ising_reference_state(model, state)
    model.pbc || throw(ArgumentError("categorical KW requires periodic boundaries"))
    model.N >= 2 || throw(ArgumentError("at least two Ising sites are required"))
    length(state) == 2^(model.N + 1) ||
        throw(DimensionMismatch("expected a chain plus one reference qubit"))
    n2 = sum(abs2, state)
    isfinite(n2) && n2 > 0 || throw(ArgumentError("state must have finite, nonzero norm"))
    return n2
end

"""
    categorical_kw_expectation(model::AnyonModel{SpinHalf}, state)

Normalized categorical KW expectation on an Ising chain with one reference
qubit. With `state = [a; b]` and the reference ordered first,
`Y_sigma = [0 D'; D 0]` and the result is
`2real(dot(b, D*a)) / (dot(a,a) + dot(b,b))`.
This Hermitian operator has eigenvalues `sqrt(2), 0, -sqrt(2)`.
Unlike `kramers_wannier_expectation`, this requires `2^(L+1)` amplitudes.
"""
function categorical_kw_expectation(model::AnyonModel{SpinHalf}, state::AbstractVector)
    _check_ising_reference_state(model, state)
    return _categorical_kw_expectation(kramers_wannier_map(model), state)
end

function _categorical_kw_expectation(D::KramersWannierMap, state)
    dim = 2^D.N
    a = view(state, 1:dim)
    b = view(state, (dim + 1):(2dim))
    return 2real(dot(b, D * a)) / sum(abs2, state)
end

# Apply one physical fusion outcome to both reference branches. Existing
# Ising kernels fix both bit ordering and false -> positive outcome convention.
function _ising_reference_measure!(out, state, basis, model_x, model_zz,
    tau, layer, site, outcome)
    dim = length(basis)
    model_a, model_b = isodd(layer) ? (model_x, model_zz) : (model_zz, model_x)
    site_b = isodd(layer) ? site : mod1(site + 1, model_x.N)
    fill!(out, zero(eltype(out)))
    for (offset, model, idx) in ((0, model_a, site), (dim, model_b, site_b))
        @inbounds for i in eachindex(basis)
            result = _apply_result(model, tau, basis[i], idx, outcome)
            amplitude = state[offset + i]
            out[offset + i] += result.w1 * amplitude
            if result.w2 != 0
                j = searchsortedfirst(basis, result.s2)
                out[offset + j] += result.w2 * amplitude
            end
        end
    end
    return out
end

"""
    ising_reference_evolution(model, state, config, samples=nothing)

Exact local fusion measurements on a periodic Ising chain plus a reference
qubit, ordered as `state = [a; b]` (see `ising_reference_state`).
Each physical outcome is shared by the two reference branches:

* odd layer: X_j in A, Z_j Z_{j+1} in B;
* even layer: Z_j Z_{j+1} in A, X_{j+1} in B.

`:Born` samples using the full joint norm and normalizes the entire vector.
`:sample` replays a Boolean matrix of size `(2periods, L)` in the same
physical-bond ordering. `false` is the positive measurement outcome.
`periods = config.t₂ - config.t₁ + 1`; rows of `samples` cover this interval.
As in `bulk_evolution`, `enable_τ_eff` halves only the final second layer.
KW is always tracked, independently of `track_y_expectation`.

Returns `(state, samples, kw_expectation_values)`. The KW vector includes
the initial value at index 1, followed by one value after each period.
No posterior, entropy, free energy, or sector classification is computed.
"""
function ising_reference_evolution(model::AnyonModel{SpinHalf,:Ising},
    state::AbstractVector, config::MeasureConfig,
    samples::Union{Nothing,AbstractMatrix{Bool}} = nothing)
    n2 = _check_ising_reference_state(model, state)
    periods = config.t₂ - config.t₁ + 1
    periods >= 1 || throw(ArgumentError("at least one period is required"))
    isfinite(config.τ) && config.τ >= 0 ||
        throw(ArgumentError("measurement strength must be finite and nonnegative"))
    config.mode in (:Born, :sample) || throw(ArgumentError("mode must be :Born or :sample"))
    if config.mode == :sample
        samples !== nothing && size(samples) == (2periods, model.N) ||
            throw(DimensionMismatch("samples must have size ($(2periods), $(model.N))"))
    elseif samples !== nothing
        throw(ArgumentError("do not supply a record in :Born mode"))
    end
    record = config.mode == :sample ? BitMatrix(samples) : falses(2periods, model.N)
    current = Vector{promote_type(Float64, eltype(state))}(state) / sqrt(n2)
    buffer = similar(current)
    model_x = AnyonModel(SpinHalf(), model.N; model_type = :Ising, pbc = true, measure_operator = :X)
    model_zz = AnyonModel(SpinHalf(), model.N; model_type = :Ising, pbc = true, measure_operator = :ZZ)
    basis = anyon_basis(model_x)
    D = kramers_wannier_map(model)
    kw = Vector{Float64}(undef, periods + 1)
    kw[1] = _categorical_kw_expectation(D, current)
    for layer in 1:(2periods)
        tau = layer == 2periods && config.enable_τ_eff ? config.τ / 2 : config.τ
        for site in 1:model.N
            outcome = config.mode == :sample ? record[layer, site] : false
            _ising_reference_measure!(buffer, current, basis, model_x, model_zz,
                tau, layer, site, outcome)
            probability = sum(abs2, buffer)
            if config.mode == :Born && rand(config.rng) >= clamp(probability, 0.0, 1.0)
                outcome = true
                _ising_reference_measure!(buffer, current, basis, model_x, model_zz,
                    tau, layer, site, outcome)
                probability = sum(abs2, buffer)
            end
            isfinite(probability) && probability > 0 ||
                throw(ArgumentError("zero-probability outcome at layer $layer, site $site"))
            record[layer, site] = outcome
            buffer ./= sqrt(probability)
            current, buffer = buffer, current
        end
        if iseven(layer)
            kw[div(layer, 2) + 1] = _categorical_kw_expectation(D, current)
        end
    end
    return (state = current, samples = record, kw_expectation_values = kw)
end
