# Ising TY+ fusion category, arXiv:2308.00747v2, Sections 3.4.1–3.4.2.
# Tuple entries follow the paper's left-to-right link indexing.
function anyon_basis(model::AnyonModel{IsingAnyon}; symmetry_block = nothing)
    symmetry_block === nothing || throw(ArgumentError("Ising anyon symmetry blocks are not implemented"))
    L = model.N
    basis = NTuple{L,UInt8}[]
    isodd(L) && return basis
    n = L ÷ 2
    n < 8sizeof(Int)-2 || throw(ArgumentError("Chain too large to enumerate"))
    for sigma_parity in (0, 1), bits in 0:(2^n-1)
        push!(basis, ntuple(L) do i
            isodd(i) == (sigma_parity == 1) ? UInt8(2) :
                UInt8((bits >> (n-cld(i, 2))) & 1)
        end)
    end
    return sort!(basis)
end

function _ising_path_check(model::AnyonModel{IsingAnyon}, state::NTuple{L,UInt8}, i::Int) where {L}
    L == model.N || throw(DimensionMismatch("Fusion path length differs from model.N"))
    1 <= i <= L || throw(ArgumentError("Link index must lie in 1:model.N"))
    all(j -> state[j] <= 2 && xor(state[j] == 2, state[mod1(j+1,L)] == 2), 1:L) ||
        throw(ArgumentError("Invalid periodic Ising fusion path"))
    return nothing
end

# A_i = 2 P_i^I - I: X at an I/η link, ZZ across a σ link.
function _ising_involution(model::AnyonModel{IsingAnyon}, state::NTuple{L,UInt8}, i::Int) where {L}
    _ising_path_check(model, state, i)
    if state[i] == 2
        return state, state[mod1(i-1,L)] == state[mod1(i+1,L)] ? 1.0 : -1.0
    end
    return Base.setindex(state, UInt8(1)-state[i], i), 1.0
end

"""
    ising_fusion_projector(model::AnyonModel{IsingAnyon}, i; channel=:I)

Sparse projector onto the fusion channel `:I` or `:eta` of the two σ anyons
adjacent to link `i`. Uses the full, sorted `anyon_basis(model)`.
The TY+ F matrix on a σ–(I/η)–σ path is `[1 1; 1 -1]/√2`.
"""
function ising_fusion_projector(model::AnyonModel{IsingAnyon}, i::Int; channel::Symbol = :I)
    channel in (:I, :eta) || throw(ArgumentError("Fusion channel must be :I or :eta"))
    1 <= i <= model.N || throw(ArgumentError("Link index must lie in 1:model.N"))
    basis = anyon_basis(model)
    rows, cols, vals = Int[], Int[], Float64[]
    s = channel === :I ? 1.0 : -1.0
    for (col, state) in enumerate(basis)
        target, weight = _ising_involution(model, state, i)
        push!(rows, col); push!(cols, col); push!(vals, 0.5)
        push!(rows, searchsortedfirst(basis, target)); push!(cols, col); push!(vals, s*weight/2)
    end
    return sparse(rows, cols, vals, length(basis), length(basis))
end

# Exact projector normalization of Eq. (3.69); no dropped constant or factor.
function actingHam(model::AnyonModel{IsingAnyon}, state::NTuple{L,UInt8}) where {L}
    output = Dict{typeof(state),Float64}()
    for i in 1:model.N
        target, weight = _ising_involution(model, state, i)
        coupling = get_interaction_param(model, isodd(i) ? :J : :h, 1.0)
        output[state] = get(output, state, 0.0) - coupling/2
        output[target] = get(output, target, 0.0) - coupling*weight/2
    end
    return output
end

# Repo weak-measurement convention, extended to the fusion observable A_i.
# Eigen-amplitudes avoid overflow even for τ=Inf; true selects η as τ→Inf.
function _apply_result(model::AnyonModel{IsingAnyon}, τ::Float64,
    state::NTuple{L,UInt8}, i::Int, sign::Bool) where {L}
    τ >= 0 || throw(ArgumentError("Measurement strength must be nonnegative"))
    target, weight = _ising_involution(model, state, i)
    r = exp(-τ)
    a = inv(sqrt(1+r*r))
    b = r*a
    c, d = (a+b)/2, (sign ? -1 : 1)*(a-b)/2
    return target == state ? (s1=state, s2=state, w1=c+d*weight, w2=0.0) :
        (s1=state, s2=target, w1=c, w2=d*weight)
end

measure_basismap(model::AnyonModel{IsingAnyon}, τ::Float64,
    state::NTuple{L,UInt8}, i::Int, sign::Bool) where {L} =
    _apply_result(model, τ, state, i, sign)

function measure_matrix(model::AnyonModel{IsingAnyon}, τ::Float64, i::Int, sign::Bool)
    τ >= 0 || throw(ArgumentError("Measurement strength must be nonnegative"))
    P = ising_fusion_projector(model, i)
    r = exp(-τ)
    a = inv(sqrt(1+r*r))
    b = r*a
    return Matrix(sign ? b*P + a*(I-P) : a*P + b*(I-P))
end
