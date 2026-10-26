"""
    mpskit_sites(N)

Physical spaces for the MPSKit backend in the repository's fusion-label/qubit
encoding. These are ordinary ComplexSpaces, NOT physical Fibonacci objects.
This preserves the existing PBC Hilbert space and encoded entropy exactly.
"""
function mpskit_sites(N::Int)
    N >= 2 || throw(ArgumentError("MPSKit backend requires at least two sites"))
    return fill(TensorKit.ComplexSpace(2), N)
end

const _KitSites = Vector{TensorKit.ComplexSpace}

function _initial_mpskit(N::Int, plus::Bool)
    sites = mpskit_sites(N)
    unit = TensorKit.ComplexSpace(1)
    v = plus ? ComplexF64[1, 1] / sqrt(2) : ComplexF64[1, 0]
    tensors = [TensorKit.TensorMap(reshape(copy(v), 1, 2, 1),
        unit ⊗ s, unit) for s in sites]
    return MPSKit.FiniteMPS(tensors), sites
end

# Only import the tensor-product notation, avoiding collisions with this
# package's FibonacciAnyon, IsingAnyon, MPS, MPO and ITensor operations.
using TensorKit: ⊗, ←

"""
    mpskit_state(amplitudes, N)
    mps_statevector(state::MPSKit.FiniteMPS)

Small-system conversion using BitBasis ordering (site 1 is the most significant
bit). The vector has length 2^N, including forbidden fusion paths. These explicit
conversions are for validation and input/output, never used by the evolution.
"""
function mpskit_state(amplitudes::AbstractVector, N::Int)
    sites = mpskit_sites(N)
    length(amplitudes) == 2^N || throw(DimensionMismatch("expected 2^N amplitudes"))
    data = permutedims(reshape(ComplexF64.(amplitudes), ntuple(_ -> 2, N)), N:-1:1)
    return MPSKit.FiniteMPS(TensorKit.Tensor(data, prod(sites)))
end

"""
    mpskit_state(state::ITensorMPS.MPS)

Convert a qubit MPS tensor by tensor, without constructing its dense state.
The input is not mutated and its amplitude normalization is preserved.
"""
function mpskit_state(state::MPS)
    N = length(state)
    spaces = mpskit_sites(N)
    sites = siteinds(state)
    all(s -> dim(s) == 2, sites) || throw(ArgumentError("only qubit site indices are supported"))
    tensors = map(1:N) do i
        indices = i == 1 ? [sites[i], linkind(state,i)] : i == N ?
            [linkind(state,i-1), sites[i]] : [linkind(state,i-1), sites[i], linkind(state,i)]
        dl = i == 1 ? 1 : dim(linkind(state,i-1))
        dr = i == N ? 1 : dim(linkind(state,i))
        data = reshape(ComplexF64.(Array(state[i], indices...)), dl, 2, dr)
        TensorKit.TensorMap(data, TensorKit.ComplexSpace(dl) ⊗ spaces[i], TensorKit.ComplexSpace(dr))
    end
    return MPSKit.FiniteMPS(tensors)
end

"""Return encoded qubit amplitudes in BitBasis order; costs exponential memory."""
function mps_statevector(state::MPSKit.FiniteMPS)
    all(s -> s == TensorKit.ComplexSpace(2), MPSKit.physicalspace.(Ref(state), 1:length(state))) ||
        throw(ArgumentError("mps_statevector requires the encoded qubit backend"))
    data = convert(Array, convert(TensorKit.TensorMap, state))
    return vec(permutedims(reshape(data, ntuple(_ -> 2, length(state))), length(state):-1:1))
end

const _kit_I = ComplexF64[1 0; 0 1]
const _kit_X = ComplexF64[0 1; 1 0]
const _kit_Z = ComplexF64[1 0; 0 -1]
const _kit_P0 = ComplexF64[1 0; 0 0]
const _kit_P1 = ComplexF64[0 0; 0 1]

# A sum of product operators, with no full-Hilbert-space intermediate.
function _kit_product_mpo(sites::_KitSites, factors)
    unit = TensorKit.ComplexSpace(1)
    tensors = map(eachindex(sites)) do i
        a = get(factors, i, _kit_I)
        TensorKit.TensorMap(reshape(a, 1, 2, 2, 1), unit ⊗ sites[i], sites[i] ⊗ unit)
    end
    return MPSKit.FiniteMPO(tensors)
end

function _kit_sum_mpo(sites::_KitSites, terms)
    return reduce(+, (c * _kit_product_mpo(sites, Dict(factors)) for (c, factors) in terms))
end

function _kit_fibonacci_terms(model, i)
    N = model.N
    N >= 3 || throw(ArgumentError("Fibonacci local operators require N >= 3"))
    1 <= i <= N || throw(BoundsError(1:N, i))
    model.pbc || 2 <= i <= N-1 || throw(ArgumentError("OBC Fibonacci centers must be interior"))
    l, r = mod1(i - 1, N), mod1(i + 1, N)
    φ = (1 + sqrt(5.0)) / 2
    model.measure_operator in (:Antiferro, :Ferro) || throw(ArgumentError("invalid Fibonacci convention"))
    σ = model.measure_operator == :Antiferro ? 1.0 : -1.0
    return [
        (σ, [l => _kit_P0, i => _kit_Z, r => _kit_P1]),
        (σ, [l => _kit_P1, i => _kit_Z, r => _kit_P0]),
        (-σ, [l => _kit_P1, i => _kit_Z, r => _kit_P1]),
        (σ * (1 - 2 / φ), [l => _kit_P0, i => _kit_Z, r => _kit_P0]),
        (-2σ * φ^(-3/2), [l => _kit_P0, i => _kit_X, r => _kit_P0]),
    ]
end

function _kit_q_terms(model::AnyonModel{FibonacciAnyon}, i)
    return _kit_fibonacci_terms(model, i)
end

function _kit_q_terms(model::Union{AnyonModel{SpinHalf,:Ising},AnyonModel{SpinHalf,:OBF}}, i)
    op = model.measure_operator
    names = op == :X ? [_kit_X] : op == :ZZ ? [_kit_Z, _kit_Z] :
        op == :XZZ ? [_kit_X, _kit_Z, _kit_Z] :
        op == :ZZX ? [_kit_Z, _kit_Z, _kit_X] : throw(ArgumentError("unsupported measurement $op"))
    1 <= i <= model.N || throw(BoundsError(1:model.N, i))
    model.pbc || i + length(names) - 1 <= model.N || throw(ArgumentError("operator crosses OBC edge"))
    σ = op in (:XZZ, :ZZX) ? -1.0 : 1.0
    return [(σ, [mod1(i + k - 1, model.N) => a for (k, a) in enumerate(names)])]
end

function _kit_measurement_coefficients(τ::Real, outcome::Bool)
    τ >= 0 || throw(ArgumentError("measurement strength must be nonnegative"))
    # Stable also for τ=Inf. Match the legacy projective threshold.
    r = τ >= 100 ? 0.0 : exp(-τ)
    den = 2sqrt(1 + r^2)
    return (1 + r) / den, (outcome ? -1 : 1) * (1 - r) / den
end

function measurement_operator_mps(model::AnyonModel, sites::_KitSites, i::Int, τ::Real, outcome::Bool)
    length(sites) == model.N || throw(DimensionMismatch("model/site count mismatch"))
    a, b = _kit_measurement_coefficients(τ, outcome)
    return a * _kit_product_mpo(sites, Dict{Int,Matrix{ComplexF64}}()) +
        _kit_sum_mpo(sites, [(b * c, factors) for (c, factors) in _kit_q_terms(model, i)])
end

function _kit_local_term(c, factors)
    sorted = sort(factors; by = first)
    ops = [TensorKit.TensorMap(a, TensorKit.ComplexSpace(2), TensorKit.ComplexSpace(2)) for (_, a) in sorted]
    return Tuple(first.(sorted)) => c * reduce(⊗, ops)
end

function _kit_hamiltonian_terms(model::AnyonModel{FibonacciAnyon})
    terms = Pair[]
    for i in (model.pbc ? (1:model.N) : (2:model.N-1))
        append!(terms, [_kit_local_term(c/2, f) for (c, f) in _kit_fibonacci_terms(model, i)])
        shift = model.measure_operator == :Antiferro ? -0.5 : 0.5
        push!(terms, _kit_local_term(shift, [i => _kit_I]))
    end
    return terms
end

function _kit_hamiltonian_terms(model::Union{AnyonModel{SpinHalf,:Ising},AnyonModel{SpinHalf,:OBF}})
    obf = model isa AnyonModel{SpinHalf,:OBF}
    J = obf ? get_interaction_param(model, :λI, 1.0) : get_interaction_param(model, :J, 1.0)
    h = obf ? J : get_interaction_param(model, :h, 1.0)
    terms = Pair[_kit_local_term(-h, [i => _kit_X]) for i in 1:model.N]
    for i in 1:(model.pbc && model.N > 2 ? model.N : model.N-1)
        push!(terms, _kit_local_term(-J, [i => _kit_Z, mod1(i+1, model.N) => _kit_Z]))
    end
    if obf
        λ = get_interaction_param(model, :λ, 1.0)
        for i in 1:(model.pbc && model.N > 2 ? model.N : model.N-2)
            j, k = mod1(i+1, model.N), mod1(i+2, model.N)
            push!(terms, _kit_local_term(λ/2, [i => _kit_X, j => _kit_Z, k => _kit_Z]))
            push!(terms, _kit_local_term(λ/2, [i => _kit_Z, j => _kit_Z, k => _kit_X]))
        end
    end
    return terms
end

function anyon_ham(model::AnyonModel, sites::_KitSites)
    length(sites) == model.N || throw(DimensionMismatch("model/site count mismatch"))
    return MPSKit.FiniteMPOHamiltonian(sites, _kit_hamiltonian_terms(model))
end

function _kit_truncation(cutoff, mindim, maxdim)
    isfinite(cutoff) && cutoff >= 0 || throw(ArgumentError("cutoff must be finite and nonnegative"))
    1 <= mindim <= maxdim || throw(ArgumentError("require 1 <= mindim <= maxdim"))
    # ITensor cutoff is discarded squared weight; TensorKit truncerror uses
    # the 2-norm of discarded singular values. Algorithms need not choose the
    # same finite-χ approximation, but these tolerances have matching units.
    return (TensorKit.truncerror(; rtol = sqrt(cutoff)) |
        TensorKit.truncrank(mindim)) & TensorKit.truncrank(maxdim)
end

# Keep a single center when changing amplitudes. Scaling all cached bond
# centers can leave incompatible gauges when MPSKit lazily reconstructs them.
function _kit_scale(ψ::MPSKit.FiniteMPS, amplitude)
    state = copy(ψ)
    center = copy(state.AC[1])
    for i in 2:length(state)
        state.AR[i]
    end
    state.AC[1] = amplitude * center
    return state
end

# Direct-sum virtual bonds preserve relative phases independently of cached
# left/right gauges. Only local tensors are materialized, never a dense state.
function _kit_add(a::MPSKit.FiniteMPS, b::MPSKit.FiniteMPS)
    length(a) == length(b) || throw(DimensionMismatch("MPS lengths differ"))
    N = length(a)
    tensors = map(1:N) do i
        x = convert(Array, i == 1 ? a.AC[i] : a.AR[i])
        y = convert(Array, i == 1 ? b.AC[i] : b.AR[i])
        if N == 1
            data = x + y
        elseif i == 1
            data = cat(x, y; dims = 3)
        elseif i == N
            data = cat(x, y; dims = 1)
        else
            data = zeros(ComplexF64, size(x,1)+size(y,1), size(x,2), size(x,3)+size(y,3))
            data[1:size(x,1), :, 1:size(x,3)] .= x
            data[size(x,1)+1:end, :, size(x,3)+1:end] .= y
        end
        TensorKit.TensorMap(data, TensorKit.ComplexSpace(size(data,1)) ⊗
            TensorKit.ComplexSpace(size(data,2)), TensorKit.ComplexSpace(size(data,3)))
    end
    return MPSKit.FiniteMPS(tensors)
end

function _kit_compress(ψ; cutoff = 1e-10, mindim = 1, maxdim = 100)
    scheme = _kit_truncation(cutoff, mindim, maxdim)
    scale = norm(ψ)
    iszero(scale) && return copy(ψ)
    return _kit_scale(MPSKit.changebonds(_kit_scale(ψ, inv(scale)),
        MPSKit.SvdCut(; trscheme = scheme); normalize = false), scale)
end

# MPSKit's MPO multiplication removes singular values below machine epsilon
# using an absolute threshold. Scale the input first so unnormalized record
# amplitudes are not inadvertently erased just because their norm is tiny.
function _kit_apply(op, ψ::MPSKit.FiniteMPS)
    scale = norm(ψ)
    iszero(scale) && return copy(ψ)
    return _kit_scale(op * _kit_scale(ψ, inv(scale)), scale)
end

function measuremap(model::AnyonModel, ψ::MPSKit.FiniteMPS, sites::_KitSites, i::Int, τ::Real,
    outcome::Bool; cutoff::Float64 = 1e-10, mindim::Int = 1, maxdim::Int = 100,
    truncate_per_event::Bool = true)
    next = _kit_apply(measurement_operator_mps(model, sites, i, τ, outcome), ψ)
    truncate_per_event && (next = _kit_compress(next; cutoff, mindim, maxdim))
    p = norm(next)^2
    p > 0 || throw(ArgumentError("selected a zero-probability measurement outcome"))
    normalize!(next)
    return next, p
end

function ee_mps(ψ::MPSKit.FiniteMPS, b::Int)
    1 <= b < length(ψ) || throw(BoundsError(1:length(ψ)-1, b))
    # Preserve legacy semantics: this API is encoded entropy, never silently
    # switch to quantum-trace entropy for categorical tensors.
    all(i -> MPSKit.physicalspace(ψ, i) == TensorKit.ComplexSpace(2), 1:length(ψ)) ||
        throw(ArgumentError("ee_mps is encoded entropy; use MPSKit.entropy for anyonic tensors"))
    return MPSKit.entropy(ψ, b)
end

anyon_eelis(model::AnyonModel, ψ::MPSKit.FiniteMPS) = [ee_mps(ψ, b) for b in 1:model.N-1]

function fibonacci_constraint_violation(ψ::MPSKit.FiniteMPS; pbc::Bool = true)
    sites = mpskit_sites(length(ψ))
    terms = [_kit_local_term(1.0, [i => _kit_P1, mod1(i+1, length(ψ)) => _kit_P1])
        for i in 1:(pbc ? length(ψ) : length(ψ)-1)]
    # MPSKit expectation_value already divides by the state norm squared.
    return real(MPSKit.expectation_value(ψ, MPSKit.FiniteMPOHamiltonian(sites, terms)))
end

# Finite-state automaton implementing the exact no-adjacent-11 projector.
function fibonacci_constraint_projector_mpo(sites::_KitSites; pbc::Bool = true)
    N = length(sites)
    N >= 2 || throw(ArgumentError("at least two sites required"))
    channels = pbc ? [(a, b) for a in 0:1 for b in 0:1] : [(0, b) for b in 0:1]
    D = length(channels)
    tensors = map(1:N) do i
        dl, dr = i == 1 ? 1 : D, i == N ? 1 : D
        data = zeros(ComplexF64, dl, 2, 2, dr)
        for b in 0:1
            if i == 1
                q = findfirst(==((pbc ? b : 0, b)), channels)
                data[1, b+1, b+1, q] = 1
            else
                for (l, (firstbit, prev)) in enumerate(channels)
                    prev == b == 1 && continue
                    if i == N
                        pbc && firstbit == b == 1 && continue
                        data[l, b+1, b+1, 1] = 1
                    else
                        r = findfirst(==((firstbit, b)), channels)
                        data[l, b+1, b+1, r] = 1
                    end
                end
            end
        end
        TensorKit.TensorMap(data, TensorKit.ComplexSpace(dl) ⊗ sites[i], sites[i] ⊗ TensorKit.ComplexSpace(dr))
    end
    return MPSKit.FiniteMPO(tensors)
end

function _mpskit_groundstate(model; sweep_times = 20, maxdim = 100, cutoff = 1e-10,
    outputlevel = 1, seed = 1234, initial_linkdim = min(8, maxdim),
    constraint_penalty = 20.0, constraint_tolerance = 1e-6)
    sweep_times >= 1 || throw(ArgumentError("sweep_times must be positive"))
    1 <= initial_linkdim <= maxdim || throw(ArgumentError("invalid initial_linkdim"))
    sites = mpskit_sites(model.N)
    H = anyon_ham(model, sites)
    terms = _kit_hamiltonian_terms(model)
    fib = model isa AnyonModel{FibonacciAnyon}
    if fib
        constraint_penalty > 0 || throw(ArgumentError("constraint_penalty must be positive"))
        for i in 1:(model.pbc ? model.N : model.N-1)
            push!(terms, _kit_local_term(constraint_penalty,
                [i => _kit_P1, mod1(i+1, model.N) => _kit_P1]))
        end
    end
    penalized = MPSKit.FiniteMPOHamiltonian(sites, terms)
    rng = MersenneTwister(seed)
    exponent_limit = ceil(Int, log2(initial_linkdim))
    dims = [min(initial_linkdim, 2^min(i, model.N-i, exponent_limit)) for i in 0:model.N]
    tensors = [TensorKit.TensorMap(randn(rng, ComplexF64, dims[i], 2, dims[i+1]),
        TensorKit.ComplexSpace(dims[i]) ⊗ sites[i], TensorKit.ComplexSpace(dims[i+1]))
        for i in 1:model.N]
    ψ = MPSKit.FiniteMPS(tensors; normalize = true)
    if fib
        ψ = _kit_compress(fibonacci_constraint_projector_mpo(sites; pbc = model.pbc) * ψ; cutoff, maxdim)
        normalize!(ψ)
    end
    ψ, _, _ = MPSKit.find_groundstate(ψ, penalized, MPSKit.DMRG2(;
        maxiter = sweep_times, verbosity = outputlevel, tol = max(cutoff, 1e-12),
        trscheme = _kit_truncation(cutoff, 1, maxdim)))
    normalize!(ψ)
    if fib
        violation = fibonacci_constraint_violation(ψ; pbc = model.pbc)
        violation <= constraint_tolerance || error("DMRG constraint violation $violation exceeds tolerance")
    end
    return ψ, real(MPSKit.expectation_value(ψ, H))
end

"""
    fibonacci_category_projectors()

Native TensorKit morphisms on τ ⊗ τ, returned as `(vacuum, tau)`.
These are category objects, not the qubit-encoded operators used by the PBC
backend. Their quantum traces are 1 and φ; no dense qubit array is implied.
"""
function fibonacci_category_projectors()
    V = TensorKit.Vect[TensorKit.FibonacciAnyon](:τ => 1)
    P1 = zeros(ComplexF64, V ⊗ V ← V ⊗ V)
    Pτ = copy(P1)
    TensorKit.block(P1, TensorKit.FibonacciAnyon(:I)) .= 1
    TensorKit.block(Pτ, TensorKit.FibonacciAnyon(:τ)) .= 1
    return P1, Pτ
end

function topological_charge_mpo(sites::_KitSites; pbc::Bool = true)
    N = length(sites)
    N >= 2 || throw(ArgumentError("at least two sites required"))
    channels = pbc ? [(a, b, c, d) for a in 0:1 for b in 0:1 for c in 0:1 for d in 0:1] :
        [(0, 0, c, d) for c in 0:1 for d in 0:1]
    D = length(channels)
    weight(a, b, c, d) = _Fibonacci_Fsymbol(a, d, b, c)
    tensors = map(1:N) do i
        dl, dr = i == 1 ? 1 : D, i == N ? 1 : D
        data = zeros(ComplexF64, dl, 2, 2, dr)
        for out in 0:1, inp in 0:1
            if i == 1
                q = findfirst(==((pbc ? out : 0, pbc ? inp : 0, out, inp)), channels)
                data[1, out+1, inp+1, q] = 1
            else
                for (l, (fo, fi, po, pi)) in enumerate(channels)
                    w = weight(po, out, pi, inp)
                    if i == N
                        pbc && (w *= weight(out, fo, inp, fi))
                        data[l, out+1, inp+1, 1] = w
                    else
                        r = findfirst(==((fo, fi, out, inp)), channels)
                        data[l, out+1, inp+1, r] = w
                    end
                end
            end
        end
        TensorKit.TensorMap(data, TensorKit.ComplexSpace(dl) ⊗ sites[i], sites[i] ⊗ TensorKit.ComplexSpace(dr))
    end
    return MPSKit.FiniteMPO(tensors)
end

function _kit_validate(model, sites, state, mc; normalized = true)
    length(sites) == length(state) == model.N || throw(DimensionMismatch("model/site/state count mismatch"))
    mc.mode in (:Born, :sample) || throw(ArgumentError("mode must be :Born or :sample"))
    mc.mode == :Born && !normalized && throw(ArgumentError("Born mode requires normalized=true"))
    mc.truncate_every_events >= 1 || throw(ArgumentError("truncate_every_events must be positive"))
    mc.t₂ - mc.t₁ + 1 >= 0 || throw(ArgumentError("negative number of periods"))
    _kit_truncation(mc.cutoff, mc.mindim, mc.maxdim)
    mc.τ >= 0 || throw(ArgumentError("measurement strength must be nonnegative"))
    if mc.enforce_fibonacci_constraint || mc.track_y_expectation
        model isa AnyonModel{FibonacciAnyon} && model.pbc ||
            throw(ArgumentError("constraint projection / Y tracking require periodic Fibonacci models"))
    end
    if normalized
        isapprox(norm(state), 1; atol = 1e-10, rtol = 1e-10) ||
            throw(ArgumentError("normalized evolution requires a normalized initial state"))
    end
    return nothing
end

_kit_compress(ψ, mc::MeasureConfig) = _kit_compress(ψ; cutoff = mc.cutoff, mindim = mc.mindim, maxdim = mc.maxdim)

function _kit_normalize!(ψ)
    p = norm(ψ)^2
    p > 0 || throw(ArgumentError("selected a zero-probability measurement outcome"))
    normalize!(ψ)
    return p
end

# Preserve the legacy per-event Born convention: sample from the compressed
# outcome-0 norm and use 1-p0 for outcome 1. Do not mix a probability-policy
# change into a backend migration. Fixed-record replay records selected norms.
function _kit_layer(model, sites, ψ, mc, sample, layer, strength, cache;
    normalized = true, operator_transform = identity)
    positions, local_model, τ = _obtain_measurement_config(model, layer, strength)
    born = mc.mode == :Born
    !born && (sample === nothing || length(sample) != length(positions)) &&
        throw(DimensionMismatch("one outcome required per measurement event"))
    bits = born ? falses(length(positions)) : copy(sample)
    current = copy(ψ)
    F = 0.0
    for (k, i) in enumerate(positions)
        ops = get!(cache, (mod1(layer, layers_per_period(model)), i, strength)) do
            (operator_transform(measurement_operator_mps(local_model, sites, i, τ, false)),
                operator_transform(measurement_operator_mps(local_model, sites, i, τ, true)))
        end
        if born
            ψ0 = _kit_apply(ops[1], current)
            mc.truncate_every_events == 1 && (ψ0 = _kit_compress(ψ0, mc))
            p0 = clamp(norm(ψ0)^2, 0.0, 1.0)
            bits[k] = rand(mc.rng) >= p0
            if bits[k]
                current = _kit_apply(ops[2], current)
                mc.truncate_every_events == 1 && (current = _kit_compress(current, mc))
                F -= log1p(-p0)
            else
                current = ψ0
                F -= log(p0)
            end
            _kit_normalize!(current)
        else
            current = _kit_apply(ops[Int(bits[k])+1], current)
            mc.truncate_every_events == 1 && (current = _kit_compress(current, mc))
            normalized && (F -= log(_kit_normalize!(current)))
        end
        if mc.truncate_every_events > 1 && (k % mc.truncate_every_events == 0 || k == length(positions))
            current = _kit_compress(current, mc)
            normalized && _kit_normalize!(current)
        end
    end
    return Measurement_outcome_mps_boundary(current, bits, Float32(F))
end

function boundary_evolution(model::AnyonModel, sites::_KitSites, state::MPSKit.FiniteMPS,
    mc::MeasureConfig, sample::Union{Nothing,BitVector} = nothing; layer_idx::Int = 1)
    _kit_validate(model, sites, state, mc)
    return _kit_layer(model, sites, state, mc, sample, layer_idx, mc.τ, Dict())
end

function bulk_evolution(model::AnyonModel, sites::_KitSites, state::MPSKit.FiniteMPS,
    mc::MeasureConfig, samples::Union{Nothing,BitMatrix} = nothing; normalized::Bool = true,
    compute_entropy::Bool = true)
    _kit_validate(model, sites, state, mc; normalized)
    periods = mc.t₂ - mc.t₁ + 1
    layers = layers_per_period(model)
    dims = (periods * layers, _samples_per_layer(model))
    if mc.mode == :sample
        samples !== nothing && size(samples) == dims || throw(DimensionMismatch("expected samples of size $dims"))
    end
    bits = mc.mode == :Born ? falses(dims) : copy(samples)
    energies, entropies = zeros(Float32, dims[1]), zeros(Float32, periods)
    ys = mc.track_y_expectation ? zeros(Float32, periods) : Float32[]
    projector = mc.enforce_fibonacci_constraint ? fibonacci_constraint_projector_mpo(sites) : nothing
    Y = mc.track_y_expectation ? topological_charge_mpo(sites) : nothing
    current, cache = copy(state), Dict()
    for period in 1:periods
        for layer in 1:layers
            row = (period-1)*layers + layer
            strength = row == dims[1] && mc.enable_τ_eff ? mc.τ/2 : mc.τ
            cols = _get_sample_column_indices(model, row)
            out = _kit_layer(model, sites, current, mc,
                mc.mode == :Born ? nothing : BitVector(bits[row, cols]),
                row, strength, cache; normalized)
            current = out.state
            bits[row, cols] = out.sample
            energies[row] = out.free_energy
            if projector !== nothing
                current = _kit_compress(projector * current, mc)
                _kit_normalize!(current)
            end
        end
        compute_entropy && (entropies[period] = Float32(ee_mps(current, model.N÷2)))
        Y !== nothing && (ys[period] = Float32(real(MPSKit.expectation_value(current, Y))))
    end
    return Measurement_outcome_mps_bulk(current, bits, energies, entropies, ys)
end

function bulk_evolution(model::AnyonModel{FibonacciAnyon}, sites::_KitSites,
    state::MPSKit.FiniteMPS, config::HybridConfig,
    schedule::Union{Nothing,HybridGateSchedule} = nothing; normalized::Bool = true)
    periods, (D, ncols) = _validate_hybrid_run(model, config, schedule)
    mc = config.measurement
    _kit_validate(model, sites, state, mc; normalized)
    born = mc.mode == :Born
    masks = born ? falses(D, ncols) : copy(schedule.measurement_mask)
    bits = born ? falses(D, ncols) : copy(schedule.outcomes)
    angles = born ? fill(NaN, D, ncols) : copy(schedule.unitary_angles)
    energies, entropies = zeros(Float32, D), zeros(Float32, periods)
    ys = mc.track_y_expectation ? zeros(Float32, periods) : Float32[]
    Y = mc.track_y_expectation ? topological_charge_mpo(sites) : nothing
    projector = mc.enforce_fibonacci_constraint ? fibonacci_constraint_projector_mpo(sites) : nothing
    current, cache = copy(state), Dict()
    layers = layers_per_period(model)
    for period in 1:periods
        for layer in 1:layers
            row = (period-1)*layers + layer
            strength = row == D && mc.enable_τ_eff ? mc.τ/2 : mc.τ
            positions, local_model, τ = _obtain_measurement_config(model, row, strength)
            cols = _get_sample_column_indices(model, row)
            F = 0.0
            for (k, i) in enumerate(positions)
                col = cols[k]
                measurement = born ? _draw_hybrid_gate(mc.rng, config.p) : masks[row, col]
                masks[row, col] = measurement
                truncate_now = k % mc.truncate_every_events == 0 || k == length(positions)
                if measurement
                    ops = get!(cache, (layer, i, τ)) do
                        (measurement_operator_mps(local_model, sites, i, τ, false),
                            measurement_operator_mps(local_model, sites, i, τ, true))
                    end
                    if born
                        ψ0 = _kit_apply(ops[1], current)
                        mc.truncate_every_events == 1 && (ψ0 = _kit_compress(ψ0, mc))
                        bits[row, col] = rand(mc.rng) >= clamp(norm(ψ0)^2, 0.0, 1.0)
                    end
                    current = _kit_apply(ops[Int(bits[row, col])+1], current)
                    truncate_now && (current = _kit_compress(current, mc))
                    normalized && (F -= log(_kit_normalize!(current)))
                else
                    θ = born ? _draw_hybrid_angle(mc.rng, config) : _scheduled_hybrid_angle(schedule, row, col, config)
                    angles[row, col] = θ
                    if θ != 0
                        Pτ = measurement_operator_mps(local_model, sites, i, Inf, false)
                        P1 = measurement_operator_mps(local_model, sites, i, Inf, true)
                        current = _kit_apply(P1 + cis(θ)*Pτ, current)
                        truncate_now && (current = _kit_compress(current, mc))
                        normalized && _kit_normalize!(current)
                    end
                end
                if mc.truncate_every_events > 1 && truncate_now
                    current = _kit_compress(current, mc)
                    normalized && _kit_normalize!(current)
                end
            end
            if projector !== nothing
                current = _kit_compress(projector * current, mc)
                _kit_normalize!(current)
            end
            energies[row] = Float32(F)
        end
        entropies[period] = Float32(ee_mps(current, model.N÷2))
        Y !== nothing && (ys[period] = Float32(real(MPSKit.expectation_value(current, Y))))
    end
    return HybridMeasurementOutcome(current, HybridGateSchedule(masks, bits, angles), energies, entropies, ys)
end

function mps_measurement_enumeration(model::AnyonModel, ψ::MPSKit.FiniteMPS,
    sites::_KitSites, positions::Vector{Int}, τ::Float64; cutoff::Float64 = 1e-10, maxdim::Int = 100)
    states, records, probabilities = [copy(ψ)], [Bool[]], [1.0]
    for i in positions
        next_states = typeof(ψ)[]
        next_records, next_probabilities = Vector{Bool}[], Float64[]
        for (state, record, probability) in zip(states, records, probabilities), bit in (false, true)
            out = _kit_apply(measurement_operator_mps(model, sites, i, τ, bit), state)
            out = _kit_compress(out; cutoff, maxdim)
            p = norm(out)^2
            # Zero-weight branches retain a zero state, not NaN amplitudes.
            p > 0 && normalize!(out)
            push!(next_states, out)
            push!(next_records, [record; bit])
            push!(next_probabilities, probability*p)
        end
        states, records, probabilities = next_states, next_records, next_probabilities
    end
    return states, records, probabilities
end

function _kit_basis_state(bits)
    sites = mpskit_sites(length(bits))
    unit = TensorKit.ComplexSpace(1)
    tensors = map(eachindex(bits)) do i
        v = zeros(ComplexF64, 1, 2, 1)
        v[1, Int(bits[i])+1, 1] = 1
        TensorKit.TensorMap(v, unit ⊗ sites[i], unit)
    end
    return MPSKit.FiniteMPS(tensors)
end

function _kit_frame_qr(states; cutoff, maxdim)
    G = Hermitian([dot(a, b) for a in states, b in states])
    # Complex Gram matrices are essential for hybrid trajectories. Do not
    # discard imaginary overlaps as the historical real-only path did.
    factor = cholesky(G; check = false)
    if !issuccess(factor)
        throw(ArgumentError("MPS Lyapunov frame lost rank; reduce n_states or truncation. No arbitrary replacement direction is inserted."))
    end
    R = factor.U
    coefficients = inv(R)
    frame = map(eachindex(states)) do i
        terms = [coefficients[j, i] * states[j] for j in 1:i if !iszero(coefficients[j, i])]
        _kit_compress(reduce(_kit_add, terms); cutoff, maxdim)
    end
    return frame, abs.(diag(R))
end

function lyapunov_spectrum_mps(model::AnyonModel, sites::_KitSites, τ::Float64,
    sample::BitMatrix; n_states::Int = 10, initial_states = nothing,
    sector::Union{Nothing,Symbol} = nothing, cutoff::Float64 = 1e-10,
    maxdim::Int = 100, truncate_every_events::Int = 1)
    layers = layers_per_period(model)
    size(sample, 1) % layers == 0 || throw(DimensionMismatch("incomplete measurement period"))
    size(sample, 2) == _samples_per_layer(model) || throw(DimensionMismatch("sample columns"))
    n_states > 0 || throw(ArgumentError("n_states must be positive"))
    project = identity
    if sector !== nothing
        model isa AnyonModel{FibonacciAnyon} && model.pbc ||
            throw(ArgumentError("sector requires periodic Fibonacci chain"))
        sector in (:trivial, :tau) || throw(ArgumentError("sector must be :trivial or :tau"))
        φ = (1+sqrt(5.0))/2
        y, other = sector == :trivial ? (φ, -1/φ) : (-1/φ, φ)
        Y, P = topological_charge_mpo(sites), fibonacci_constraint_projector_mpo(sites)
        project = function (ψ)
            physical = _kit_compress(P * ψ; cutoff, maxdim)
            _kit_compress(inv(y-other) * (_kit_add(_kit_apply(Y, physical), _kit_scale(physical, -other))); cutoff, maxdim)
        end
    end
    if initial_states === nothing
        basis = anyon_basis(model)
        states = [_kit_basis_state([!iszero((b.buf >> (model.N-i)) & 1) for i in 1:model.N])
            for b in basis[1:min(n_states, length(basis))]]
    else
        isempty(initial_states) && throw(ArgumentError("empty initial frame"))
        states = copy.(initial_states)
        all(ψ -> ψ isa MPSKit.FiniteMPS && length(ψ) == model.N, states) ||
            throw(ArgumentError("initial frame must contain compatible MPSKit states"))
    end
    if sector !== nothing || initial_states !== nothing
        states = project.(states)
        states, _ = _kit_frame_qr(states; cutoff, maxdim)
    end
    periods = size(sample, 1) ÷ layers
    result = zeros(length(states), periods)
    mc = MeasureConfig(τ = τ, t₂ = 1, mode = :sample, enable_τ_eff = false,
        cutoff = cutoff, maxdim = maxdim, truncate_every_events = truncate_every_events)
    for t in 1:periods
        record = sample[(t-1)*layers+1:t*layers, :]
        states = [project(bulk_evolution(model, sites, ψ, mc, record;
            normalized = false, compute_entropy = false).state) for ψ in states]
        states, stretches = _kit_frame_qr(states; cutoff, maxdim)
        result[:, t] = -log.(stretches)
    end
    return result
end

function _kit_prepend_ancilla(ψ::MPSKit.FiniteMPS, bit::Bool)
    unit, site = TensorKit.ComplexSpace(1), TensorKit.ComplexSpace(2)
    data = zeros(ComplexF64, 1, 2, 1)
    data[1, Int(bit)+1, 1] = 1
    first = TensorKit.TensorMap(data, unit ⊗ site, unit)
    return MPSKit.FiniteMPS([first; ψ.AC[1]; ψ.AR[2:end]])
end

function add_reference_qubits(model::AnyonModel, ψ::MPSKit.FiniteMPS, sites::_KitSites,
    site_idx::Int = 1; k_new::Int = 1, verbose::Bool = false)
    length(sites) == length(ψ) || throw(DimensionMismatch("site/state count mismatch"))
    1 <= site_idx <= length(ψ) || throw(BoundsError(sites, site_idx))
    k_new in (0,1) || throw(ArgumentError("k_new must be 0 or 1"))
    k_new == 0 && return copy(ψ), copy(sites)
    a = _kit_product_mpo(sites, Dict(site_idx => _kit_P0)) * ψ
    b = _kit_product_mpo(sites, Dict(site_idx => _kit_P1)) * ψ
    # Direct controlled copy; no SWAP or exponential state conversion.
    joint = _kit_add(_kit_prepend_ancilla(a, false), _kit_prepend_ancilla(b, true))
    return joint, mpskit_sites(length(sites)+1)
end

function _kit_ancilla_identity(op::MPSKit.FiniteMPO)
    unit, site = TensorKit.ComplexSpace(1), TensorKit.ComplexSpace(2)
    first = TensorKit.TensorMap(reshape(_kit_I, 1, 2, 2, 1), unit ⊗ site, site ⊗ unit)
    return MPSKit.FiniteMPO([first; parent(op)])
end

function reference_rdm(model::AnyonModel{FibonacciAnyon}, ψ::MPSKit.FiniteMPS)
    length(ψ) == model.N+1 || throw(DimensionMismatch("expected ancilla plus N system sites"))
    data = convert(Array, ψ.AC[1])
    a = reshape(data, 2, :)
    return a * a'
end

function topological_charge_sharpening(model::AnyonModel{FibonacciAnyon}, sites::_KitSites,
    state::MPSKit.FiniteMPS, mc::MeasureConfig, samples::Union{Nothing,BitMatrix} = nothing)
    model.pbc || throw(ArgumentError("charge sharpening requires PBC"))
    ψ = copy(state)
    _kit_normalize!(ψ)
    _kit_validate(model, sites, ψ, mc)
    φ = (1+sqrt(5.0))/2
    Yψ = _kit_compress(topological_charge_mpo(sites) * ψ, mc)
    ψ1 = _kit_compress(inv(φ+1/φ) * (_kit_add(Yψ, _kit_scale(ψ, inv(φ)))), mc)
    ψτ = _kit_compress(inv(φ+1/φ) * (_kit_add(_kit_scale(ψ, φ), _kit_scale(Yψ, -1))), mc)
    current = _kit_compress(_kit_add(_kit_prepend_ancilla(ψτ, false), _kit_prepend_ancilla(ψ1, true)), mc)
    _kit_normalize!(current)
    periods, layers = mc.t₂-mc.t₁+1, layers_per_period(model)
    dims = (periods*layers, _samples_per_layer(model))
    born = mc.mode == :Born
    !born && (samples === nothing || size(samples) != dims) && throw(DimensionMismatch("sample shape must be $dims"))
    bits = born ? falses(dims) : copy(samples)
    energies, entropies = zeros(Float32, dims[1]), zeros(Float32, periods)
    P = mc.enforce_fibonacci_constraint ? _kit_ancilla_identity(fibonacci_constraint_projector_mpo(sites)) : nothing
    cache = Dict()
    for period in 1:periods
        for layer in 1:layers
            row = (period-1)*layers+layer
            strength = row == dims[1] && mc.enable_τ_eff ? mc.τ/2 : mc.τ
            cols = _get_sample_column_indices(model, row)
            out = _kit_layer(model, sites, current, mc, born ? nothing : BitVector(bits[row, cols]),
                row, strength, cache; operator_transform = _kit_ancilla_identity)
            current, bits[row, cols], energies[row] = out.state, out.sample, out.free_energy
            if P !== nothing
                current = _kit_compress(P * current, mc)
                _kit_normalize!(current)
            end
        end
        entropies[period] = Float32(ee_mps(current, 1))
    end
    return Measurement_outcome_mps_bulk(current, bits, energies, entropies)
end
