# =============================================================================
# Hamiltonian definitions
# (Basis.jl defines the Hilbert spaces/bases; local Hamiltonian maps, `actingHam`
# and `anyon_ham` live here together with the measurement operators.)
# =============================================================================
"""
    Fibomap(state::T, i::Int; ferro::Bool=false) where {N, T <: BitStr{N}}

Apply Fibonacci anyon projection term at site i (Temperley-Lieb generator).

# Fibonacci Hamiltonian structure:
The Hamiltonian H = ∑_i π_i acts on the fusion tree with local terms.
Each term π_i depends on the local fusion outcomes at sites i-1, i, i+1.

# Fusion rules for Fibonacci anyons:
- τ × τ = 1 + τ  (two τ's can fuse to vacuum 1 or τ)
- τ × 1 = τ     (τ with vacuum gives τ)
- 1 × 1 = 1     (two vacuums give vacuum)

# Returns
- `(state, X_state, diag_weight, off_diag_weight)`:
  - `state`: original state (diagonal contribution)
  - `X_state`: flipped state at site i (off-diagonal)
  - `diag_weight`, `off_diag_weight`: matrix element weights
"""
function Fibomap(state::T, i::Int; ferro::Bool = false) where {N,T<:BitStr{N}}
    ϕ = (1 + √5) / 2  # Golden ratio
    fl = bmask(T, N)
    X_state = flip(state, fl >> (i-1))

    # Read bit at site i: 0 = τ anyon
    # BitStr indexing: bit at site i is at position N - i + 1
    bit_i = readbit(state, N - i + 1)

    # Weights depend on local fusion outcome
    # Antiferro (ground state favors alternating): negative weights
    # Ferro (ground state favors aligned): positive weights
    sign = ferro ? 1 : -1

    if bit_i == 0  # site i is vacuum
        diag_weight = sign * ϕ^(-1)
    else           # site i is τ
        diag_weight = sign * ϕ^(-2)
    end
    off_diag_weight = sign * ϕ^(-3/2)

    return state, X_state, diag_weight, off_diag_weight
end

# Legacy aliases for backward compatibility
antimap(state::T, i::Int) where {N,T<:BitStr{N}} = Fibomap(state, i; ferro = false)
ferromap(state::T, i::Int) where {N,T<:BitStr{N}} = Fibomap(state, i; ferro = true)

function Isingmap(state::T, i::Int, pbc::Bool = true; kwargs...) where {N,T<:BitStr{N}}
    # H = - J ∑ Z_i Z_{i+1} - h ∑ X_i
    @assert 1 <= i <= N "i is expected to be in [1, $N], but got $i"

    fl = bmask(T, N)
    X(state, i) = flip(state, fl >> (i-1))
    J, h = get(kwargs, :J, 1.0), get(kwargs, :h, 1.0)

    # OBC special case: last site has no ZZ term
    if !pbc && i == N
        return X(state, i), -h
    end

    # Get site indices with PBC wrapping
    i1 = pbc ? mod1(i + 1, N) : i + 1

    bit_i = readbit(state, N - i + 1)
    bit_i1 = readbit(state, N - i1 + 1)

    zz_i1i2 = (bit_i == bit_i1) ? -J : J

    return state, X(state, i), zz_i1i2, -h
end

"""
    OBFmap(state::T, i::Int, pbc::Bool=true) where {N, T <: BitStr{N}}

Apply O'Brien-Fendley terms (X_i Z_{i} Z_{i+2} + Z_i Z_{i+1} X_{i+2}) at site i.

Returns (output_states, weights) for the XZZ + ZZX terms in the OBF Hamiltonian.
For OBC, valid range is 1 ≤ i ≤ N-2.
For PBC, valid range is 1 ≤ i ≤ N with periodic wrapping.
"""
function OBFmap(state::T, i::Int, pbc::Bool = true) where {N,T<:BitStr{N}}
    fl = bmask(T, N)
    X(state, j) = flip(state, fl >> (j-1))

    # Get site indices with PBC wrapping
    if pbc
        i1 = mod1(i + 1, N)
        i2 = mod1(i + 2, N)
    else
        @assert 1 <= i <= N - 2 "For OBC, i must be in [1, N-2], got $i"
        i1 = i + 1
        i2 = i + 2
    end

    # Read bits: note BitStr is 1-indexed from right, but we count from left
    # bit at site j is at position N - j + 1 in BitStr indexing
    bit_i = readbit(state, N - i + 1)
    bit_i1 = readbit(state, N - i1 + 1)
    bit_i2 = readbit(state, N - i2 + 1)

    # XZZ term: X_i Z_{i+1} Z_{i+2}
    # Z_{i+1} Z_{i+2} eigenvalue: +1 if same, -1 if different
    zz_i1i2 = (bit_i1 == bit_i2) ? 1 : -1
    xzz_state = X(state, i)
    xzz_weight = zz_i1i2  # λ for Hamiltonian (energy lowering)

    # ZZX term: Z_i Z_{i+1} X_{i+2}
    zz_ii1 = (bit_i == bit_i1) ? 1 : -1
    zzx_state = X(state, i2)
    zzx_weight = zz_ii1

    return xzz_state, zzx_state, xzz_weight, zzx_weight
end

function count_subBitStr(state::T) where {N,T<:BitStr{N}}
    n = length(state)
    n < 3 && return 0

    str100, str101, str001 = T(4), T(5), T(1) # 100, 101, 001
    num=0

    mask=bmask(T, 1, 2, 3)
    for i = 1:(n-2) # start from string right to left
        substr = state & (mask << (i-1))
        if substr == str101
            num += 1
        end
        str101 <<= 1
    end

    return num
end

"""
    actingHam(model::AnyonModel{FibonacciAnyon}, state::T) where {N, T <: BitStr{N}}

Act the Fibonacci anyon Hamiltonian on a given state.

# Physics
The Fibonacci anyon chain Hamiltonian consists of Temperley-Lieb generators e_i:
- H = -∑_i e_i (ferromagnetic) or H = +∑_i e_i (antiferromagnetic)
- Each e_i acts on the fusion space at site i with constraints from neighbors

# Fusion Rules
For Fibonacci anyons: τ × τ = 1 + τ
- Configuration 0x0 (neighbors are trivial): allows Fibomap operation
- Configuration 101, 100, 001: contributes diagonal energy from fusion constraints
- Configuration 111 (three consecutive τ): additional fusion contribution (PBC only)

# Returns
Dict{T, Float64}: mapping from output states to their coefficients
"""
function actingHam(model::AnyonModel{FibonacciAnyon}, state::T) where {N,T<:BitStr{N}}
    @assert num_digits(T) == N "The length of system is expected to be $N, but got $(num_digits(T))"

    pbc = model.pbc
    ferro = (model.measure_operator == :Ferro)
    sign = ferro ? 1 : -1

    # Helper to apply Fibomap and accumulate results
    function apply_fibomap!(output::Dict{T,Float64}, i::Int)
        state1, state2, weight1, weight2 = Fibomap(state, i; ferro = ferro)
        output[state1] = get(output, state1, 0.0) + weight1
        output[state2] = get(output, state2, 0.0) + weight2
    end

    output = Dict{T,Float64}()

    # Diagonal contribution from 101, 100, 001 patterns
    # Sign depends on ferro/antiferro: H = sign * ∑ (fusion constraints)
    output[state] = get(output, state, 0.0) + sign * count_subBitStr(state)

    # Bulk terms: sites 2 to N-1, apply Fibomap where 0x0 pattern exists
    mask = bmask(T, N, N-2)  # Mask to check neighbors
    for i = 2:(N-1)
        if state & (mask >> (i-2)) == 0  # Check 0x0 pattern
            apply_fibomap!(output, i)
        end
    end

    # Periodic boundary condition terms
    if pbc
        # Site 1: check if neighbors (site N-1 and site 2) form 0x0
        if state[1] == 0 && state[N-1] == 0
            apply_fibomap!(output, 1)
        end
        # Site N: check if neighbors (site 2 and site N) form 0x0
        if state[2] == 0 && state[N] == 0
            apply_fibomap!(output, N)
        end

        # 111 fusion contributions at boundaries
        mask1 = bmask(T, N, 2)      # Check pattern 1xxxx01
        mask2 = bmask(T, N-1, 1)    # Check pattern 01xxxx1
        if state & mask1 == mask1
            output[state] = get(output, state, 0.0) + sign
        end
        if state & mask2 == mask2
            output[state] = get(output, state, 0.0) + sign
        end
    end

    return output
end

function actingHam(model::AnyonModel{SpinHalf,:Ising}, state::T) where {N,T<:BitStr{N}}
    @assert num_digits(T) == N "The length of system is expected to be $N, but got $(num_digits(T))"

    pbc = model.pbc
    J = get_interaction_param(model, :J, 1.0)
    h = get_interaction_param(model, :h, 1.0)
    fl = bmask(T, N)
    X(st, j) = flip(st, fl >> (j-1))

    output = Dict{T,Float64}()

    # ZZ terms: -J ∑ Z_i Z_{i+1}
    n_bonds = pbc ? N : N - 1
    for i = 1:n_bonds
        i1 = mod1(i + 1, N)
        zz_val = (readbit(state, N-i+1) == readbit(state, N-i1+1)) ? -J : J
        output[state] = get(output, state, 0.0) + zz_val
    end

    # X terms: -h ∑ X_i
    for i = 1:N
        output[X(state, i)] = get(output, X(state, i), 0.0) - h
    end

    return output
end

function actingHam(model::AnyonModel{SpinHalf,:OBF}, state::T) where {N,T<:BitStr{N}}
    @assert num_digits(T) == N "The length of system is expected to be $N, but got $(num_digits(T))"

    pbc = model.pbc
    λ = get_interaction_param(model, :λ, 1.0)  # OBF coupling strength
    λI = get_interaction_param(model, :λI, 1.0)  # Ising coupling strength

    # Generate OBF model Hamiltonian: H = λ ∑ (X_i Z_{i+1} Z_{i+2} + Z_i Z_{i+1} X_{i+2}) - X - ZZ
    output = Dict{T,Float64}()
    for i = 1:N
        s1, s2, w1, w2 = Isingmap(state, i, pbc)
        output[s1] = get(output, s1, 0.0) + λI * w1
        output[s2] = get(output, s2, 0.0) + λI * w2
        state1, state2, weight1, weight2 = OBFmap(state, i, pbc)
        output[state1] = get(output, state1, 0.0) + λ/2 * weight1
        output[state2] = get(output, state2, 0.0) + λ/2 * weight2
    end

    return output
end

"""
    anyon_ham(model::AnyonModel)

Construct the Hamiltonian matrix for a 1D anyon chain.

# Arguments
- `model::AnyonModel`: An `AnyonModel` object specifying the anyon type, system size, 
  boundary conditions, interaction type, and model parameters (J, h, λ, etc.).

# Returns
- `Matrix{Float64}`: The Hamiltonian matrix constructed in the chosen basis.

# Supported Models
- **Fibonacci Anyons**: Supports `:Antiferro` and `:Ferro` interaction terms.
- **Ising Anyons**: Transverse field Ising model with parameters `J` and `h`.
- **OBF Anyons**: O'Brien-Fendley model with parameter `λ`.
- **Heisenberg chain**: Spin-1/2 XXZ chain `H = J ∑ (XX + YY + Δ ZZ)` (J = +1 AFM, J = -1 FM) with parameters `J` and `Δ`.

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; model_fibo = AnyonModel(FibonacciAnyon(), N; pbc=true, measure_operator=:Antiferro);

julia> H_fibo = anyon_ham(model_fibo); size(H_fibo)
(7, 7)

julia> model_ising = AnyonModel(SpinHalf(), N; model_type=:Ising, pbc=true, J=1.0, h=1.0); H_ising = anyon_ham(model_ising); size(H_ising)
(16, 16)
```
"""
function anyon_ham(model::AnyonModel{AT}) where {AT<:AbstractAnyonBasis}
    basis = anyon_basis(model)
    l = length(basis)
    H = zeros(Float64, (l, l))

    for i = 1:l
        output = actingHam(model, basis[i])
        for (m, weight) in output
            j = searchsortedfirst(basis, m)
            H[i, j] += weight
        end
    end

    return H
end
# Another method to write Fibonacci Hamiltonian is using the Measurement operator sum. For example, H = -∑ X_i, where X_i is the Temperley-Lieb generator acting on site i-1, i, and i+1. Pilis = [FibonacciChain.measure_matrix(BitStr{16, Int}, 1000.0, idx, 0) for idx in 1:N]. H = -sum(Pilis). This two Hamiltonian difference is not a constant, but like a arc in conformal energy spectrum below arc, but they have the same eigenstates.

"""
    anyon_ham(model::AnyonModel, ::Type{T}, k::Int; symmetry_block=nothing) where {N, T <: BitStr{N}}
    anyon_ham(model::AnyonModel, k::Int; symmetry_block=nothing)

Construct Hamiltonian matrix in specific momentum sector for 1D anyon chain.
    
# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (including J, h, λ, etc.)
- `T::Type`: BitStr type specifying chain length N (optional)
- `k::Int`: Momentum sector (0 ≤ k ≤ N-1)
- `symmetry_block`: Topological charge sector (optional)

# Returns
- `Matrix`: Hamiltonian matrix in chosen momentum sector (ComplexF64, or real if k=0 or k=N/2)

# Examples
```jldoctest
julia> using FibonacciChain

julia> model = AnyonModel(FibonacciAnyon(), 6; pbc=true);

julia> H_k0 = anyon_ham(model, 0);

julia> size(H_k0)[1] > 0
true
```
"""
function anyon_ham(
    model::AnyonModel{AT},
    ::Type{T},
    k::Int;
    symmetry_block = nothing,
) where {N,T<:BitStr{N},AT<:AbstractAnyonBasis}
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert symmetry_block === nothing || symmetry_block in [0, 1, :tau, :trivial] "symmetry_block is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $symmetry_block"

    basisK, basis_dic = anyon_basis(model, k, symmetry_block = symmetry_block)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    H = zeros(ComplexF64, (l, l))

    for i = 1:l
        n = basisK[i]
        output = actingHam(model, n)
        for (m, weight) in output
            mbar, d = get_representative(m)
            if mbar ∈ basisK
                j = searchsortedfirst(basisK, mbar)
                Yn = sqrt(length(basis_dic[n])) / N
                Ym = sqrt(length(basis_dic[mbar])) / N
                H[i, j] += Yn/Ym * omegak^d * weight
            end
        end
    end

    if k == 0 || k == div(N, 2)
        H = real(H)
    end
    H = (H + H') / 2
    return H
end

anyon_ham(
    model::AnyonModel{AT},
    k::Int;
    symmetry_block = nothing,
) where {AT<:AbstractAnyonBasis} =
    anyon_ham(model, BitStr{model.N,Int}, k; symmetry_block = symmetry_block)

# ============================================================================
# Spin-1/2 Heisenberg (XXZ) chain
# H = J ∑_i (X_i X_{i+1} + Y_i Y_{i+1} + Δ Z_i Z_{i+1}), J = +1 AFM (default), J = -1 FM
# It lives on the full 2^N spin-1/2 product basis (see `spin_half_basis` in Basis.jl),
# shared with the Ising and OBF chains.
# ============================================================================

"""
    Heisenbergmap(state::T, i::Int, pbc::Bool=true; J::Float64=1.0, Δ::Float64=1.0) where {N, T <: BitStr{N}}

Apply the Heisenberg (XXZ) bond term on sites `(i, i+1)`.

The Hamiltonian is H = J ∑_i (X_i X_{i+1} + Y_i Y_{i+1} + Δ Z_i Z_{i+1}),
with J = +1 antiferromagnetic (default) and J = -1 ferromagnetic. There is no field term.
For OBC the last site carries no bond, so zero weights are returned at `i == N`.

# Returns
- `(state, swapped_state, diag_weight, off_diag_weight)`:
  - `state`: original state (diagonal contribution `±J Δ`)
  - `swapped_state`: state with the anti-aligned pair `(i, i+1)` exchanged (from `J(XX + YY)`);
    equals `state` when the pair is aligned, since `XX + YY` annihilates `|00⟩` and `|11⟩`
  - `diag_weight`, `off_diag_weight`: matrix element weights (`2J` for the exchange)
"""
function Heisenbergmap(
    state::T,
    i::Int,
    pbc::Bool = true;
    J::Float64 = 1.0,
    Δ::Float64 = 1.0,
) where {N,T<:BitStr{N}}
    @assert 1 <= i <= N "i is expected to be in [1, $N], but got $i"

    # OBC special case: last site has no bond
    if !pbc && i == N
        return state, state, 0.0, 0.0
    end

    i1 = pbc ? mod1(i + 1, N) : i + 1
    bit_i = readbit(state, N - i + 1)
    bit_i1 = readbit(state, N - i1 + 1)

    if bit_i == bit_i1
        # Aligned pair: ZZ contributes +JΔ, XX + YY vanishes
        return state, state, J * Δ, 0.0
    else
        # Anti-aligned pair: ZZ contributes -JΔ, XX + YY exchanges the pair with amplitude 2J
        swapped = flip(state, bmask(T, N - i + 1, N - i1 + 1))
        return state, swapped, -J * Δ, 2J
    end
end

"""
    actingHam(model::AnyonModel{SpinHalf,:Heisenberg}, state::T) where {N, T <: BitStr{N}}

Act the spin-1/2 Heisenberg (XXZ) Hamiltonian on a given state:

    H = J ∑_i (X_i X_{i+1} + Y_i Y_{i+1} + Δ Z_i Z_{i+1})

with J = +1 antiferromagnetic (default) and J = -1 ferromagnetic. There is no field term.

# Parameters (passed as `kwargs` to `AnyonModel`)
- `J`: overall coupling strength (default 1.0)
- `Δ`: Ising (Z) anisotropy (default 1.0, i.e. the isotropic XXX chain)

# Returns
`Dict{T, Float64}`: mapping from output states to their coefficients
"""
function actingHam(model::AnyonModel{SpinHalf,:Heisenberg}, state::T) where {N,T<:BitStr{N}}
    @assert num_digits(T) == N "The length of system is expected to be $N, but got $(num_digits(T))"

    pbc = model.pbc
    J = get_interaction_param(model, :J, 1.0)
    Δ = get_interaction_param(model, :Δ, 1.0)

    output = Dict{T,Float64}()
    n_bonds = pbc ? N : N - 1
    for i = 1:n_bonds
        s1, s2, w1, w2 = Heisenbergmap(state, i, pbc; J = J, Δ = Δ)
        output[s1] = get(output, s1, 0.0) + w1
        w2 == 0 || (output[s2] = get(output, s2, 0.0) + w2)
    end

    return output
end

"""
    measure_basismap(model::AnyonModel, τ::Float64, state::T, i::Int, sign::Bool) where {T}

Map single basis state under measurement operation at site i.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, measure_operator)
- `τ::Float64`: Measurement strength parameter
- `state::T`: Input basis state
- `i::Int`: Measurement site index (1 ≤ i ≤ N)
- `sign::Bool`: Measurement outcome (false for τ, true for 1)

# Returns
- `Tuple`: Either `(basis1, basis2, coeff1, coeff2)` for superposition output or `(basis, coeff)` for single output

Maps individual basis states according to measurement protocols and fusion rules.

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; T = BitStr{N, Int};

julia> model = AnyonModel(FibonacciAnyon(), N; pbc=true, measure_operator=:Antiferro);

julia> state = T(0b0100);  # Single τ at site 2

julia> τ = 1.0;  # Measurement strength

julia> result = measure_basismap(model, τ, state, 2, false);

julia> length(result) ∈ [2, 4]  # Returns 2 or 4 elements depending on configuration
true
```
"""
function measure_basismap(
    model::AnyonModel{AT},
    τ::Float64,
    state::T,
    i::Int,
    sign::Bool,
) where {T,AT<:AbstractAnyonBasis}
    # default for PBC system, map basis (not state!!!), and index count from the left.
    @assert num_digits(T) == model.N "State length mismatch: expected $(model.N), got $(num_digits(T))"
    return _apply_result(model, τ, state, i, sign)
end

function _apply_result(
    model::AnyonModel{FibonacciAnyon},
    τ::Float64,
    state::T,
    i::Int,
    sign::Bool,
)::@NamedTuple{s1::T, s2::T, w1::Float64, w2::Float64} where {T}
    measure_operator = model.measure_operator

    N = model.N
    fl=bmask(T, N)
    X(state, i) = flip(state, fl >> (i-1))
    ϕ = (1+√5)/2

    if measure_operator == :reset && τ >= 1e2
        cstτ = 0.5
        coef = sign ? -0.5 : 0.5
        value = (state[N-i+1] == 0) ? cstτ + coef : cstτ - coef
        return (s1 = state, s2 = state, w1 = value, w2 = 0.0)
    else
        if τ >= 1e2
            # true is 1/vaccuum, false is 0/τ, Fibonacci anyon
            cstτ = 0.5
            coef = sign ? -0.5 : 0.5
        else
            cstτ = (exp(τ) + 1) / (2 * √(exp(2τ) + 1))
            coef =
                sign ? (1 - exp(τ)) / (2 * √(exp(2τ) + 1)) :
                (exp(τ) - 1) / (2 * √(exp(2τ) + 1))
        end

        if 2 <= i <= N-1
            mask=bmask(T, 1, 2, 3) << (N-i-1)
            str100, str101, str010, str001, str000 = T(4) << (N-i-1),
            T(5) << (N-i-1),
            T(2) << (N-i-1),
            T(1) << (N-i-1),
            T(0) << (N-i-1)
            if state & mask == str000
                return (
                    s1 = state,
                    s2 = X(state, i),
                    w1 = cstτ+coef*(1-2ϕ^(-1)),
                    w2 = -2*coef*ϕ^(-3/2),
                )
            elseif state & mask == str010
                return (
                    s1 = state,
                    s2 = X(state, i),
                    w1 = cstτ+coef*(2ϕ^(-1)-1),
                    w2 = -2*coef*ϕ^(-3/2),
                )
            elseif state & mask == str001
                return (s1 = state, s2 = state, w1 = cstτ+coef, w2 = 0.0)
            elseif state & mask == str100
                return (s1 = state, s2 = state, w1 = cstτ+coef, w2 = 0.0)
            elseif state & mask == str101
                return (s1 = state, s2 = state, w1 = cstτ-coef, w2 = 0.0)
            end
        end

        if model.pbc
            if i == 1 #count from the left
                mask=bmask(T, N, N-1, 1)
                str100, str101, str010, str001, str000 =
                    bmask(T, 1), bmask(T, N-1, 1), bmask(T, N), bmask(T, N-1), T(0)
                if state & mask == str000
                    return (
                        s1 = state,
                        s2 = X(state, i),
                        w1 = cstτ+coef*(1-2ϕ^(-1)),
                        w2 = -2*coef*ϕ^(-3/2),
                    )
                elseif state & mask == str010
                    return (
                        s1 = state,
                        s2 = X(state, i),
                        w1 = cstτ+coef*(2ϕ^(-1)-1),
                        w2 = -2*coef*ϕ^(-3/2),
                    )
                elseif state & mask == str001
                    return (s1 = state, s2 = state, w1 = cstτ+coef, w2 = 0.0)
                elseif state & mask == str100
                    return (s1 = state, s2 = state, w1 = cstτ+coef, w2 = 0.0)
                elseif state & mask == str101
                    return (s1 = state, s2 = state, w1 = cstτ-coef, w2 = 0.0)
                end
            elseif i == N #count from the left
                mask=bmask(T, N, 2, 1)
                str100, str101, str010, str001, str000 =
                    bmask(T, 2), bmask(T, N, 2), bmask(T, 1), bmask(T, N), T(0)
                if state & mask == str000
                    return (
                        s1 = state,
                        s2 = X(state, i),
                        w1 = cstτ+coef*(1-2ϕ^(-1)),
                        w2 = -2*coef*ϕ^(-3/2),
                    )
                elseif state & mask == str010
                    return (
                        s1 = state,
                        s2 = X(state, i),
                        w1 = cstτ+coef*(2ϕ^(-1)-1),
                        w2 = -2*coef*ϕ^(-3/2),
                    )
                elseif state & mask == str001
                    return (s1 = state, s2 = state, w1 = cstτ+coef, w2 = 0.0)
                elseif state & mask == str100
                    return (s1 = state, s2 = state, w1 = cstτ+coef, w2 = 0.0)
                elseif state & mask == str101
                    return (s1 = state, s2 = state, w1 = cstτ-coef, w2 = 0.0)
                end
            end
        end
    end
end

# Standalone Ising measurement logic (no model allocation needed)
function _apply_result_ising(
    measure_operator::Symbol,
    N::Int,
    pbc::Bool,
    τ::Float64,
    state::T,
    i::Int,
    sign::Bool,
)::@NamedTuple{s1::T, s2::T, w1::Float64, w2::Float64} where {T}
    fl = bmask(T, N)
    X_flip(st, j) = flip(st, fl >> (j-1))

    # Common coefficients for all operators except :reset/:Z
    if τ >= 1e2
        cstτ = 0.5
        coef = sign ? -0.5 : 0.5
    else
        cstτ = cosh(τ/2) / √(2cosh(τ))
        coef = sign ? -sinh(τ/2) / √(2cosh(τ)) : sinh(τ/2) / √(2cosh(τ))
    end

    if measure_operator == :X
        return (s1 = state, s2 = X_flip(state, i), w1 = cstτ, w2 = coef)

    elseif measure_operator == :ZZ
        if 1 <= i <= N-1
            eigenvalue = ((state >> (N - i)) & 1) == ((state >> (N - i - 1)) & 1) ? 1 : -1
        elseif pbc && i == N
            eigenvalue = (state & 1) == (state >> (N-1) & 1) ? 1 : -1
        end
        return (s1 = state, s2 = state, w1 = cstτ + coef * eigenvalue, w2 = 0.0)
    elseif measure_operator ∈ (:reset, :Z)
        # :reset/:Z has flipped sign convention for coef
        coef_z = sign ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
        if τ >= 1e2
            coef_z = sign ? -0.5 : 0.5
        end
        eigenvalue = (state[N-i+1] == 0) ? 1 : -1
        return (s1 = state, s2 = state, w1 = cstτ + coef_z * eigenvalue, w2 = 0.0)
    end
end

function _apply_result(
    model::AnyonModel{SpinHalf,:Ising},
    τ::Float64,
    state::T,
    i::Int,
    sign::Bool,
)::@NamedTuple{s1::T, s2::T, w1::Float64, w2::Float64} where {T}
    return _apply_result_ising(
        model.measure_operator,
        model.N,
        model.pbc,
        τ,
        state,
        i,
        sign,
    )
end

function _apply_result(
    model::AnyonModel{SpinHalf,:Heisenberg},
    τ::Float64,
    state::T,
    i::Int,
    sign::Bool,
)::@NamedTuple{s1::T, s2::T, w1::Float64, w2::Float64} where {T}
    # SWAP-bond measurement: the layer applies exp(s·τ·h_bond) on bond (i, i+1),
    # with s = +1 (sign = false) / -1 (sign = true) and h_bond = J (XX + YY + Δ ZZ).
    measure_operator = model.measure_operator
    @assert measure_operator == :SWAP "measure_operator must be :SWAP for Heisenberg chain, but got $measure_operator"

    N = model.N
    pbc = model.pbc
    J = get_interaction_param(model, :J, 1.0)
    Δ = get_interaction_param(model, :Δ, 1.0)

    # OBC: the last site has no bond (same OBC convention as _apply_result_ising)
    @assert pbc || 1 <= i <= N - 1 "Index i must be in [1, N-1] for open BC (SWAP)"

    i1 = pbc ? mod1(i + 1, N) : i + 1
    bit_i = readbit(state, N - i + 1)
    bit_i1 = readbit(state, N - i1 + 1)

    s = sign ? -1.0 : 1.0

    if τ >= 1e2
        # Projective limit (γ = 1): exp(s·τ·h_bond), up to its operator norm, becomes
        # the projector onto the extremal-eigenvalue sector of h_bond (lowest for
        # s = -1, highest for s = +1). Bond eigenvalues: aligned |00⟩/|11⟩ have
        # E_al = JΔ; the anti-aligned block has E_lo = -JΔ - 2J (singlet) and
        # E_hi = -JΔ + 2J (triplet0). At Δ = ∓1 the aligned sector is degenerate
        # with the winning anti-aligned eigenstate, so both are kept.
        E_al = J * Δ
        E_lo = -J * Δ - 2J
        E_hi = -J * Δ + 2J
        if bit_i == bit_i1
            keep = sign ? (E_al <= E_lo) : (E_al >= E_hi)
            return (s1 = state, s2 = state, w1 = keep ? 1.0 : 0.0, w2 = 0.0)
        elseif sign ? (E_lo <= E_al) : (E_hi >= E_al)
            # anti-aligned pair: project onto the singlet (s = -1) or triplet0 (s = +1)
            swapped = flip(state, bmask(T, N - i + 1, N - i1 + 1))
            return (s1 = state, s2 = swapped, w1 = 0.5, w2 = sign ? -0.5 : 0.5)
        else
            # whole anti-aligned block is not extremal (e.g. s = -1 with Δ < -1): killed
            return (s1 = state, s2 = state, w1 = 0.0, w2 = 0.0)
        end
    end

    if bit_i == bit_i1
        # Aligned pair: h_bond is diagonal with eigenvalue JΔ
        return (s1 = state, s2 = state, w1 = exp(s * τ * J * Δ), w2 = 0.0)
    else
        # Anti-aligned pair: on the {|01⟩, |10⟩} block h_bond = -JΔ·I + 2J·σ_x,
        # so exp(s·τ·h_bond) = exp(-sτJΔ) [cosh(2τJ)·I + s·sinh(2τJ)·σ_x]
        swapped = flip(state, bmask(T, N - i + 1, N - i1 + 1))
        return (
            s1 = state,
            s2 = swapped,
            w1 = exp(-s * τ * J * Δ) * cosh(2τ * J),
            w2 = s * exp(-s * τ * J * Δ) * sinh(2τ * J),
        )
    end
end

function _apply_result(
    model::AnyonModel{SpinHalf,:OBF},
    τ::Float64,
    state::T,
    i::Int,
    sign::Bool,
)::@NamedTuple{s1::T, s2::T, w1::Float64, w2::Float64} where {T}
    measure_operator = model.measure_operator

    N = model.N
    fl = bmask(T, N)
    X(state, i) = flip(state, fl >> (i-1))

    # Common coefficients for all operators, here in constrast to Ising case the sign convention is inverse.
    if τ >= 1e2
        cstτ = 0.5
        coef = sign ? 0.5 : -0.5
    else
        cstτ = cosh(τ/2) / √(2cosh(τ))
        coef = sign ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
    end

    # Helper: get ZZ eigenvalue for sites (j1, j2)
    zz_eigen(j1, j2) = ((state >> (N - j1)) & 1) == ((state >> (N - j2)) & 1) ? 1 : -1

    if model.pbc
        i1, i2 = mod1(i + 1, N), mod1(i + 2, N)
    else
        @assert 1 <= i <= N - 2 "Index i must be in [1, N-2] for OBC (OBF)"
        i1, i2 = i + 1, i + 2
    end

    if measure_operator ∈ (:ZZ, :X)
        return _apply_result_ising(measure_operator, N, model.pbc, τ, state, i, sign)
    elseif measure_operator == :XZZ
        # return XZZ components, and coefficients
        return (s1 = state, s2 = X(state, i), w1 = cstτ, w2 = coef * zz_eigen(i1, i2))
    elseif measure_operator == :ZZX
        # return ZZX components, and coefficients
        return (s1 = state, s2 = X(state, i2), w1 = cstτ, w2 = coef * zz_eigen(i, i1))
    end
end

function measure_matrix(
    model::AnyonModel{AT},
    τ::Float64,
    idx::Int,
    sign::Bool,
) where {AT<:AbstractAnyonBasis}

    if model.measure_operator ∈ [:Ferro, :Antiferro]
        @assert model.pbc || (2 <= idx <= model.N-1) "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif model.measure_operator == :ZZ
        @assert model.pbc || (1 <= idx <= model.N-1) "Index idx must be in [1, N-1] for open BC (IsingZZ)"
    elseif model.measure_operator ∈ (:X, :Z, :reset, :resetFibo)
        @assert model.pbc || (1 <= idx <= model.N) "Index idx must be in [1, N] for open BC (IsingX)"
    elseif model.measure_operator ∈ (:XZZ, :ZZX)
        @assert model.pbc || (1 <= idx <= model.N-2) "Index idx must be in [1, N-2] for open BC (OBF)"
    elseif model.measure_operator == :SWAP
        @assert model.pbc || (1 <= idx <= model.N-1) "Index idx must be in [1, N-1] for open BC (SWAP)"
    else
        error("Unknown measure class: $(model.basis)")
    end

    basis = anyon_basis(model)
    l = length(basis)
    Bmatrix = zeros(l, l)

    for i = 1:l
        s1, s2, w1, w2 = measure_basismap(model, τ, basis[i], idx, sign)

        if w2 == 0
            Bmatrix[i, i] += w1
        else
            j2 = searchsortedfirst(basis, s2)
            Bmatrix[i, i] += w1
            Bmatrix[i, j2] += w2
        end
    end

    return Bmatrix
end

"""
    measuremap(model::AnyonModel, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool)
    measuremap(model::AnyonModel, ψ::MPS, sites, i::Int, τ::Float64, sign::Bool; 
               cutoff::Float64=1e-10, maxdim::Int=100)

Apply measurement operator to quantum state and return post-measurement state.

# Methods

## Vector version (from `Measurement.jl`)
Apply measurement to a state vector.

### Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Input quantum state vector
- `idx::Int`: Measurement site index
- `sign::Bool`: Measurement outcome (false for τ, true for 1)

### Returns
- `Vector{ET}`: Post-measurement quantum state (unnormalized)

## MPS version (from `MPSMeasurement.jl`)
Apply measurement to an MPS state.

### Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `ψ::MPS`: Input quantum state
- `sites`: ITensor site indices
- `i::Int`: Measurement site
- `τ::Float64`: Measurement strength parameter
- `sign::Bool`: Measurement outcome (false for τ, true for 1)
- `cutoff::Float64=1e-10`: MPS truncation cutoff
- `maxdim::Int=100`: Maximum bond dimension

### Returns
- `MPS`: Post-measurement quantum state (normalized)
- `Float64`: Measurement probability
"""
function measuremap(
    model::AnyonModel{AT},
    τ::Float64,
    state::Vector{ET},
    idx::Int,
    sign::Bool,
) where {ET,AT<:AbstractAnyonBasis}
    basis = anyon_basis(model)
    mapped_state = zeros(ET, length(basis))
    return _measuremap_impl!(mapped_state, basis, model, τ, state, idx, sign)
end

"""
    measuremap!(mapped_state, model, τ, state, idx, sign, basis)

In-place version of `measuremap` that writes result into pre-allocated `mapped_state` buffer.
The `basis` argument should be obtained from `anyon_basis(model)` and cached for reuse.
"""
function measuremap!(
    mapped_state::Vector{ET},
    model::AnyonModel{AT},
    τ::Float64,
    state::Vector{ET},
    idx::Int,
    sign::Bool,
    basis,
) where {ET,AT<:AbstractAnyonBasis}
    fill!(mapped_state, zero(ET))
    return _measuremap_impl!(mapped_state, basis, model, τ, state, idx, sign)
end

# Type-stable inner implementation using function barrier pattern.
# The concrete type of `basis` is known here, making the loop type-stable.
function _measuremap_impl!(
    mapped_state::Vector{ET},
    basis::Vector{BT},
    model::AnyonModel{AT},
    τ::Float64,
    state::Vector{ET},
    idx::Int,
    sign::Bool,
) where {ET,BT,AT<:AbstractAnyonBasis}
    l = length(basis)
    @assert length(state) == l "state length is expected to be $l, but got $(length(state))"

    @inbounds for i = 1:l
        result = _apply_result(model, τ, basis[i], idx, sign)
        mapped_state[i] += result.w1 * state[i]
        if result.w2 != 0
            j2 = searchsortedfirst(basis, result.s2)
            mapped_state[j2] += result.w2 * state[i]
        end
    end

    return mapped_state
end

function laddermeasuremap(
    model::AnyonModel{AT},
    τ::Float64,
    state::Vector{ET},
    idx::Int,
    sign::Bool,
) where {ET,AT<:AbstractAnyonBasis}
    # input a superposition state, and output the braided state
    @assert model.pbc || (2 <= idx <= model.N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    basis=anyon_basis(model)
    l=length(basis)
    @assert l^2 == length(state) "state length is expected to be $(l^2), but got $(length(state))"
    mapped_state = zeros(ET, length(state))
    @inbounds for i = 1:l
        @inbounds for j = 1:l
            output1 = measure_basismap(model, τ, basis[i], idx, sign)
            output2 = measure_basismap(model, τ, basis[j], idx, sign)
            if length(output1) == 4 && length(output2) == 4
                basisi1, basisi2, coefi1, coefi2=output1
                basisj1, basisj2, coefj1, coefj2=output2
                i2=searchsortedfirst(basis, basisi2)
                j2=searchsortedfirst(basis, basisj2)
                # Here noting that the state is a vertorizing density matrix, so the index is i+(j-1)*len, not state[i], state[j]
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi1*coefj1
                mapped_state[(i-1)*l+j2]+=state[(i-1)*l+j]*coefi1*coefj2
                mapped_state[(i2-1)*l+j]+=state[(i-1)*l+j]*coefi2*coefj1
                mapped_state[(i2-1)*l+j2]+=state[(i-1)*l+j]*coefi2*coefj2
            elseif length(output1) == 4 && length(output2) == 2
                basisi1, basisi2, coefi1, coefi2=output1
                basisj, coefj=output2
                i2=searchsortedfirst(basis, basisi2)
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi1*coefj
                mapped_state[(i2-1)*l+j]+=state[(i-1)*l+j]*coefi2*coefj
            elseif length(output1) == 2 && length(output2) == 4
                basisi, coefi=output1
                basisj1, basisj2, coefj1, coefj2=output2
                j2=searchsortedfirst(basis, basisj2)
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi*coefj1
                mapped_state[(i-1)*l+j2]+=state[(i-1)*l+j]*coefi*coefj2
            else
                basisi, coefi=output1
                basisj, coefj=output2
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi*coefj
            end
        end
    end

    return mapped_state
end

"""
    measurement_enumeration(model::AnyonModel, τ::Float64, initial_state::Vector{ET}, measurement_sites::Vector{Int}) where {ET}

Enumerate all trajectories of measurements on a given initial state.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `initial_state::Vector{ET}`: Initial quantum state vector
- `measurement_sites::Vector{Int}`: Sites to measure

# Returns
- `Tuple{Vector{Vector{ET}}, Vector{Vector{Int64}}, Vector{Float64}}`:
  (final_states, trajectories, probabilities)

# Examples
```jldoctest
julia> using FibonacciChain, LinearAlgebra

julia> model = AnyonModel(FibonacciAnyon(), 4; pbc=true, measure_operator=:Antiferro);

julia> basis = anyon_basis(model);

julia> state = normalize(ones(length(basis)));

julia> states, trajs, probs = measurement_enumeration(model, 1.0, state, [2]);

julia> length(states) == 2  # Two possible outcomes
true
```
"""
function measurement_enumeration(
    model::AnyonModel{AT},
    τ::Float64,
    initial_state::Vector{ET},
    measurement_sites::Vector{Int},
) where {ET,AT<:AbstractAnyonBasis}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    # Initialize, only one initial state
    current_level_states = [copy(initial_state)]
    current_level_trajectories = [Bool[]]
    current_level_probabilities = [1.0]

    for (measurement_idx, site) in enumerate(measurement_sites)
        next_level_states = Vector{Vector{ET}}()
        next_level_trajectories = Vector{Vector{Int64}}()
        next_level_probabilities = Vector{Float64}()

        # Branching for each current state
        for (state_idx, state) in enumerate(current_level_states)
            current_trajectory = current_level_trajectories[state_idx]
            current_prob = current_level_probabilities[state_idx]

            state_after_p = measuremap(model, τ, state, site, false)
            prob_p = real(dot(state_after_p, state_after_p))

            normalized_state_p = state_after_p / sqrt(prob_p)
            new_trajectory_p = [current_trajectory; false]
            new_prob_p = current_prob * prob_p

            push!(next_level_states, normalized_state_p)
            push!(next_level_trajectories, new_trajectory_p)
            push!(next_level_probabilities, new_prob_p)

            state_after_m = measuremap(model, τ, state, site, true)
            prob_m = real(dot(state_after_m, state_after_m))

            normalized_state_m = state_after_m / sqrt(prob_m)
            new_trajectory_m = [current_trajectory; true]
            new_prob_m = current_prob * prob_m

            push!(next_level_states, normalized_state_m)
            push!(next_level_trajectories, new_trajectory_m)
            push!(next_level_probabilities, new_prob_m)

        end

        current_level_states = next_level_states
        current_level_trajectories = next_level_trajectories
        current_level_probabilities = next_level_probabilities

    end


    return current_level_states, current_level_trajectories, current_level_probabilities
end


"""
    measurement_tree_visualization(sample::BitMatrix)

Visualize a fixed measurement trajectory with staggered odd/even layers.

# Arguments
- `sample::BitMatrix`: Measurement outcome sequences (rows = layers, cols = sites).

# Behavior
- Prints each layer with odd/even rows staggered to show the alternating measurement pattern.
- Uses `●` for outcome 1 (true) and `○` for outcome 0 (false).
"""
function measurement_tree_visualization(sample::BitMatrix)
    println("Measurement Tree Visualization")
    println("=" ^ 40)

    n_layers, n_sites = size(sample)

    for layer in 1:n_layers
        # Convert bits to symbols
        symbols = [bit ? "●" : "○" for bit in sample[layer, :]]

        # Join with spacing
        line = join(symbols, "     ")

        # Stagger: even layers shifted right by half spacing
        indent = iseven(layer) ? "     " : "  "

        println("Layer $(lpad(layer, 2)): $(indent)$(line)")
    end
end


function _obtain_measurement_config(
    model::AnyonModel{FibonacciAnyon},
    layer_idx::Int,
    τ::Float64 = 1.0,
)
    measurement_sites = iseven(layer_idx) ? collect(1:2:model.N) : collect(2:2:model.N)
    measure_operator = :Antiferro
    measure_anyon_model = AnyonModel(
        FibonacciAnyon(),
        model.N;
        pbc = model.pbc,
        measure_operator = measure_operator,
    )
    measure_strength = τ
    return measurement_sites, measure_anyon_model, measure_strength
end

function _obtain_measurement_config(
    model::AnyonModel{SpinHalf,:Ising},
    layer_idx::Int,
    τ::Float64 = 1.0,
)
    measurement_sites = collect(1:model.N)
    measure_operator = iseven(layer_idx) ? :ZZ : :X
    measure_anyon_model = AnyonModel(
        model.basis,
        model.N;
        model_type = :Ising,
        pbc = model.pbc,
        measure_operator = measure_operator,
    )
    measure_strength = τ
    return measurement_sites, measure_anyon_model, measure_strength
end

function _obtain_measurement_config(
    model::AnyonModel{SpinHalf,:Heisenberg},
    layer_idx::Int,
    τ::Float64 = 1.0,
)
    # SWAP-bond measurement-only circuit: each layer applies exp(s·τ·h_bond) on a
    # staggered set of bonds, with h_bond = J (XX + YY + Δ ZZ).
    # Odd layers act on bonds (1,2), (3,4), ..., even layers on (2,3), (4,5), ...
    # (the last bond wraps around to site 1 for PBC).
    N = model.N
    J = get_interaction_param(model, :J, 1.0)
    Δ = get_interaction_param(model, :Δ, 1.0)

    measurement_sites = isodd(layer_idx) ? collect(1:2:N) : collect(2:2:N)
    if !model.pbc
        # OBC: bond (i, i+1) requires i <= N-1, so the dangling bond start i = N is excluded
        filter!(<(N), measurement_sites)
    end
    measure_anyon_model = AnyonModel(
        model.basis,
        N;
        model_type = :Heisenberg,
        pbc = model.pbc,
        measure_operator = :SWAP,
        J = J,
        Δ = Δ,
    )
    measure_strength = τ
    return measurement_sites, measure_anyon_model, measure_strength
end

function _obtain_measurement_config(
    model::AnyonModel{SpinHalf,:OBF},
    layer_idx::Int,
    τ::Float64 = 1.0,
)
    # OBF 8-layer period structure:
    # Layer 1, 13: √XZZ (sites 1,4,7...)
    # Layer 2, 12: √ZZX (sites 1,4,7...)
    # Layer 3, 11: √XZZ (sites 2,5,8...)
    # Layer 4, 10: √ZZX (sites 2,5,8...)
    # Layer 5, 9:  √XZZ (sites 3,6,9...)
    # Layer 6, 8:  √ZZX (sites 3,6,9...)
    # Layer 7:     X (all sites)
    # Final 14: ZZ
    # theoretically we use √ZZ, √XZZ₁, √ZZX₁, √XZZ₂, √ZZX₂, √XZZ₃, √ZZX₃, X, √ZZX₃, √XZZ₃, √ZZX₂, √XZZ₂, √ZZX₁, √XZZ₁, √ZZ, which is symmetricly trotterized and can be
    # √XZZ₁, √ZZX₁, √XZZ₂, √ZZX₂, √XZZ₃, √ZZX₃, X, √ZZX₃, √XZZ₃, √ZZX₂, √XZZ₂, √ZZX₁, √XZZ₁, ZZ
    # 1        2      3      4      5      6    7    8      9     10      11    12     13    14 

    phase = mod1(layer_idx, 14)
    λ = get_interaction_param(model, :λ, 1.0)  # OBF coupling strength
    λI = get_interaction_param(model, :λI, 1.0)  # Ising coupling strength

    N = model.N
    if phase == 1 || phase == 13
        # √XZZ measurement at (1, 4, 7)
        measurement_sites = collect(1:3:N)
        measure_operator = :XZZ
        measure_strength = λ*τ/4
    elseif phase == 2 || phase == 12
        # √ZZX measurement at (1, 4, 7)
        measurement_sites = collect(1:3:N)
        measure_operator = :ZZX
        measure_strength = λ*τ/4
    elseif phase == 3 || phase == 11
        # √XZZ measurement at (2, 5, 8)
        measurement_sites = collect(2:3:N)
        measure_operator = :XZZ
        measure_strength = λ*τ/4
    elseif phase == 4 || phase == 10
        # √ZZX measurement at (2, 5, 8)
        measurement_sites = collect(2:3:N)
        measure_operator = :ZZX
        measure_strength = λ*τ/4
    elseif phase == 5 || phase == 9
        # √XZZ measurement at (3, 6, 9)
        measurement_sites = collect(3:3:N)
        measure_operator = :XZZ
        measure_strength = λ*τ/4
    elseif phase == 6 || phase == 8
        # √ZZX measurement at (3, 6, 9)
        measurement_sites = collect(3:3:N)
        measure_operator = :ZZX
        measure_strength = λ*τ/4
    elseif phase == 7
        # X measurement at all sites
        measurement_sites = collect(1:N)
        measure_operator = :X
        measure_strength = λI * τ
    elseif phase == 14
        # ZZ measurement at all sites
        measurement_sites = collect(1:N)
        measure_operator = :ZZ
        measure_strength = λI * τ
    end

    measure_anyon_model = AnyonModel(
        model.basis,
        N;
        model_type = :OBF,
        pbc = model.pbc,
        measure_operator = measure_operator,
        λ = λ,
        λI = λI,
    )
    return measurement_sites, measure_anyon_model, measure_strength
end

"""
    _get_sample_column_indices(model::AnyonModel, layer_idx::Int) -> Vector{Int}

Get the column indices in the samples BitMatrix for a given layer.

For models where different layers have different numbers of measurement sites,
this maps the measurement sites to fixed column positions in the samples matrix.

# Returns
- Vector of column indices where this layer's samples should be stored/read
"""
function _get_sample_column_indices(model::AnyonModel{FibonacciAnyon}, layer_idx::Int)
    # Fibonacci: alternating even/odd sites, always N÷2 measurements
    # Columns 1:(N÷2) are used for all layers
    return collect(1:(model.N÷2))
end

function _get_sample_column_indices(model::AnyonModel{SpinHalf,:Ising}, layer_idx::Int)
    # Ising: all N sites measured each layer
    return collect(1:model.N)
end

function _get_sample_column_indices(model::AnyonModel{SpinHalf,:Heisenberg}, layer_idx::Int)
    # Heisenberg: SWAP-bond layers carry one sample per measured bond. Odd layers
    # measure bonds (1,2), (3,4), ..., even layers (2,3), (4,5), ...; for OBC the
    # even layer skips the dangling bond start i = N, so it holds one sample fewer.
    n_bonds = model.N ÷ 2
    (!model.pbc && iseven(layer_idx)) && (n_bonds -= 1)
    return collect(1:n_bonds)
end

function _get_sample_column_indices(model::AnyonModel{SpinHalf,:OBF}, layer_idx::Int)
    # OBF: different layers measure different sites, but all map to columns 1:N
    # The column index equals the site index being measured
    phase = mod1(layer_idx, 14)
    N = model.N

    if phase == 1 || phase == 13
        # XZZ: sites 1,4,7,... → columns 1,4,7,...
        return collect(1:3:N)
    elseif phase == 2 || phase == 12
        # ZZX: sites 1,4,7,... → columns 1,4,7,...
        return collect(1:3:N)
    elseif phase == 3 || phase == 11
        # XZZ: sites 2,5,8,... → columns 2,5,8,...
        return collect(2:3:N)
    elseif phase == 4 || phase == 10
        # ZZX: sites 2,5,8,... → columns 2,5,8,...
        return collect(2:3:N)
    elseif phase == 5 || phase == 9
        # XZZ: sites 3,6,9,... → columns 3,6,9,...
        return collect(3:3:N)
    elseif phase == 6 || phase == 8
        # ZZX: sites 3,6,9,... → columns 3,6,9,...
        return collect(3:3:N)
    elseif phase == 7
        # X: all sites 1:N → columns 1:N
        return collect(1:N)
    elseif phase == 14
        # ZZ: all sites 1:N → columns 1:N
        return collect(1:N)
    end
end

"""
    _samples_per_layer(model::AnyonModel) -> Int

Return the number of sample columns needed per layer (maximum across all layer types).
"""
_samples_per_layer(model::AnyonModel{FibonacciAnyon}) = model.N ÷ 2
_samples_per_layer(model::AnyonModel{SpinHalf,:Ising}) = model.N
_samples_per_layer(model::AnyonModel{SpinHalf,:OBF}) = model.N  # Max of all layer types
_samples_per_layer(model::AnyonModel{SpinHalf,:Heisenberg}) = model.N ÷ 2  # One sample per measured bond


struct Measurement_outcome_bulk{ET}
    state::ET
    samples::BitMatrix
    free_energys::Vector{Float32}
    entanglement_entropys::Vector{Float32}
    y_expectation_values::Vector{Float32}
end

# Keep the previous four-argument constructor available for callers (such as the
# reference-probe routines) that do not measure the topological Y charge.
Measurement_outcome_bulk(state, samples, free_energys, entanglement_entropys) =
    Measurement_outcome_bulk(
        state,
        samples,
        free_energys,
        entanglement_entropys,
        Float32[],
    )

struct Measurement_outcome_boundary{T}
    state::Vector{T}
    sample::BitVector
    free_energy::Float32
end

"""
Configuration struct for measurement evolution parameters.

# Fields
- `τ::Float64`: Measurement strength parameter
- `t₂::Int`: 2*Number of measurement layers (time steps)
- `rng::MersenneTwister`: Random number generator (default: `MersenneTwister()`)
- `mode::Symbol`: Sampling mode, one of `:sample`, `:Born` (default: `:sample`)
- `t₁::Int`: Starting layer index for evolution (default: 1)
- `verbose::Bool`: Verbosity flag for detailed output (default: false)
- `enable_τ_eff::Bool`: Whether to enable half-strength measurement for the last layer (default: true)
- `track_y_expectation::Bool`: Record the Fibonacci topological-symmetry expectation
  value after every complete period (default: false)
- `cutoff::Float64`: MPS truncation cutoff (default: `1e-12`)
- `maxdim::Int`: Maximum MPS bond dimension (default: 1000)
- `truncate_every_events::Int`: Number of measurement events between MPS truncations
  (default: 1)
"""
Base.@kwdef struct MeasureConfig
    τ::Float64
    t₂::Int
    rng::MersenneTwister = MersenneTwister()
    mode::Symbol = :sample
    t₁::Int = 1
    verbose::Bool = false
    enable_τ_eff::Bool = true
    track_y_expectation::Bool = false
    x₂::Int = 1
    x₁::Int = 1
    cutoff::Float64 = 1e-12
    maxdim::Int = 1000
    truncate_every_events::Int = 1
end

"""
    layers_per_period(model::AnyonModel) -> Int

Return the number of measurement layers per evolution period.
- Fibonacci: 2 layers
- Ising: 2 layers (X, ZZ)
- OBF: 14 layers (√XZZ₁, √ZZX₁, √XZZ₂, √ZZX₂, √XZZ₃, √ZZX₃, X, √ZZX₃, √XZZ₃, √ZZX₂, √XZZ₂, √ZZX₁, √XZZ₁, ZZ), here OBF represents XZZ + ZZX. At the end plus a final √ZZ layer.
- Heisenberg: 2 layers of SWAP-bond gates, staggered even/odd bonds
"""
layers_per_period(model::AnyonModel{FibonacciAnyon}) = 2
layers_per_period(model::AnyonModel{SpinHalf,:Ising}) = 2
layers_per_period(model::AnyonModel{SpinHalf,:OBF}) = 14
layers_per_period(model::AnyonModel{SpinHalf,:Heisenberg}) = 2

"""
    boundary_evolution(model::AnyonModel, state::Vector{T}, measure_config::MeasureConfig, 
                       sample::Union{Nothing, BitVector}=nothing; layer_idx::Int=1)
    boundary_evolution(model::AnyonModel, sites, state::MPS, sample::BitVector, 
                       measure_config::MeasureConfig; 
                       cutoff::Float64=1e-12, maxdim::Int=1000, layer_idx::Int=1)

Evolve a quantum state under a single layer of boundary measurements with a given trajectory.

# Methods

## Vector version (from `Measurement.jl`)
Evolve a state vector under boundary measurements.

### Arguments
- `model::AnyonModel`: The anyon model specifying system parameters.
- `state::Vector{T}`: The initial state vector.
- `measure_config::MeasureConfig`: Configuration containing measurement strength `τ` and mode.
- `sample::Union{Nothing, BitVector}=nothing`: The measurement outcomes for the layer.
- `layer_idx::Int=1`: The layer index for measurement configuration.

### Returns
- `Measurement_outcome_boundary`: A struct containing `state`, `sample`, and `free_energy`.

## MPS version (from `MPSMeasurement.jl`)
Evolve an MPS state under boundary measurements.

### Arguments
- `model::AnyonModel`: The anyon model specifying system parameters.
- `sites`: The ITensor site indices.
- `state::MPS`: The initial MPS state.
- `sample::BitVector`: The measurement outcomes for the layer.
- `measure_config::MeasureConfig`: Configuration containing measurement strength `τ` and mode (`:sample` or `:Born`).
- `cutoff::Float64=1e-12`: MPS truncation cutoff.
- `maxdim::Int=1000`: Maximum bond dimension for MPS operations.
- `layer_idx::Int=1`: The layer index for measurement configuration.

### Returns
- `Measurement_outcome_mps_boundary`: A struct containing:
  - `state::MPS`: The evolved MPS state.
  - `samples::BitVector`: The measurement outcomes for the layer.
  - `free_energy::Float64`: The free energy associated with the measurement layer.
"""
function boundary_evolution(
    anyon_model::AnyonModel{AT},
    state::Vector{T},
    measure_config::MeasureConfig,
    sample::Union{Nothing,BitVector} = nothing;
    layer_idx::Int = 1,
) where {T,AT<:AbstractAnyonBasis}

    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    τ_eff = measure_config.enable_τ_eff ? measure_config.τ / 2 : measure_config.τ
    if measure_config.mode == :sample
        N = anyon_model.N
        size(sample, 1) == _samples_per_layer(anyon_model) ||
            error("sample size mismatch with anyon_model $(N)")
        return _apply_measurement_layer(
            anyon_model,
            τ_eff,
            state,
            sample;
            layer_idx = layer_idx,
        )
    elseif measure_config.mode == :Born
        return _stochastic_measurement_layer(
            anyon_model,
            τ_eff,
            state;
            layer_idx = layer_idx,
            rng = measure_config.rng,
            verbose = measure_config.verbose,
        )
    end
end

"""
    _apply_measurement_layer(model::AnyonModel, τ::Float64, state::Vector{T}, 
                              layer_sample::BitVector; layer_idx::Int=1) where {T}

Apply deterministic measurements to a layer with given measurement outcomes.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `state::Vector{T}`: Quantum state vector
- `layer_sample::BitVector`: Measurement outcomes for the layer
- `layer_idx::Int=1`: Layer index (1-based) to determine measurement pattern

# Returns
- `Measurement_outcome_boundary`: A struct containing the post-measurement state, sample, and total free energy.
"""
function _apply_measurement_layer(
    anyon_model::AnyonModel{AT},
    τ::Float64,
    state::Vector{T},
    layer_sample::BitVector;
    layer_idx::Int64 = 1,
    normalized::Bool = true,
) where {T,AT<:AbstractAnyonBasis}
    # Helper function to apply deterministic measurements to a layer, connect measure on each site together.

    measurement_sites, measure_anyon_model, measurement_strength =
        _obtain_measurement_config(anyon_model, layer_idx, τ)

    mop = anyon_model.measure_operator
    N = anyon_model.N
    pbc = anyon_model.pbc
    if mop ∈ (:Ferro, :Antiferro)
        @assert pbc || 2 <= measurement_sites[1] || measurement_sites[end] <= N-1 "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif mop == :ZZ
        @assert pbc || 1 <= measurement_sites[1] || measurement_sites[end] <= N-1 "Index idx must be in [1, N-1] for open BC (ZZ)"
    elseif mop ∈ (:X, :Z, :reset, :resetFibo)
        @assert pbc || 1 <= measurement_sites[1] || measurement_sites[end] <= N "Index idx must be in [1, N] for open BC (X)"
    elseif mop ∈ (:XZZ, :ZZX)
        @assert pbc || 1 <= measurement_sites[1] || measurement_sites[end] <= N-2 "Index idx must be in [1, N-2] for open BC (OBF)"
    end

    # Cache basis and pre-allocate buffer to avoid per-site allocations
    basis = anyon_basis(measure_anyon_model)
    l = length(basis)
    buf = Vector{T}(undef, l)
    current_state = copy(state)

    if normalized
        total_free_energy = zero(real(T))
        for (idx, sign) in enumerate(layer_sample)
            # Apply measurement into pre-allocated buffer
            fill!(buf, zero(T))
            _measuremap_impl!(
                buf,
                basis,
                measure_anyon_model,
                measurement_strength,
                current_state,
                measurement_sites[idx],
                sign,
            )
            prob = sum(abs2, buf)
            total_free_energy += -log(prob)
            buf .*= inv(sqrt(prob))
            current_state, buf = buf, current_state  # swap buffers
        end

        return Measurement_outcome_boundary(
            current_state,
            layer_sample,
            Float32(total_free_energy),
        )
    else
        for (idx, sign) in enumerate(layer_sample)
            # Apply measurement into pre-allocated buffer (no normalization)
            fill!(buf, zero(T))
            _measuremap_impl!(
                buf,
                basis,
                measure_anyon_model,
                measurement_strength,
                current_state,
                measurement_sites[idx],
                sign,
            )
            current_state, buf = buf, current_state  # swap buffers
        end

        return Measurement_outcome_boundary(
            current_state,
            layer_sample,
            Float32(0.0),
        )
    end
end

"""
    _stochastic_measurement_layer(model::AnyonModel, τ::Float64, state::Vector{T};
                   layer_idx::Int=1, rng::MersenneTwister=MersenneTwister(), 
                   verbose::Bool=false) where {T}

Perform random measurement on a layer using Born rule sampling.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `state::Vector{T}`: Quantum state vector
- `layer_idx::Int=1`: Layer index (1-based) to determine measurement pattern
- `rng::MersenneTwister`: Random number generator
- `verbose::Bool=false`: Whether to print debug information

# Returns
- `Measurement_outcome_boundary`: A struct containing the post-measurement state, sample outcomes, and free energy.
"""
function _stochastic_measurement_layer(
    anyon_model::AnyonModel{AT},
    τ::Float64,
    state::Vector{T};
    layer_idx::Int64 = 1,
    rng::MersenneTwister = MersenneTwister(),
    verbose::Bool = false,
) where {T,AT<:AbstractAnyonBasis}

    measurement_sites, measure_anyon_model, measurement_strength =
        _obtain_measurement_config(anyon_model, layer_idx, τ)

    mop = anyon_model.measure_operator
    N = anyon_model.N
    pbc = anyon_model.pbc
    if mop ∈ (:Ferro, :Antiferro)
        @assert pbc || 2 <= measurement_sites[1] || measurement_sites[end] <= N-1 "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif mop == :ZZ
        @assert pbc || 1 <= measurement_sites[1] || measurement_sites[end] <= N-1 "Index idx must be in [1, N-1] for open BC (ZZ)"
    elseif mop ∈ (:X, :Z, :reset, :resetFibo)
        @assert pbc || 1 <= measurement_sites[1] || measurement_sites[end] <= N "Index idx must be in [1, N] for open BC (X)"
    elseif mop ∈ (:XZZ, :ZZX)
        @assert pbc || 1 <= measurement_sites[1] || measurement_sites[end] <= N-2 "Index idx must be in [1, N-2] for open BC (OBF)"
    end

    n = length(measurement_sites)
    sample = BitVector(undef, n)
    F_layer = 0.0

    # Cache basis and pre-allocate buffers
    basis = anyon_basis(measure_anyon_model)
    l = length(basis)
    buf0 = Vector{T}(undef, l)
    buf1 = Vector{T}(undef, l)
    current_state = copy(state)

    for (i, site) in enumerate(measurement_sites)
        # Compute 0-branch
        fill!(buf0, zero(T))
        _measuremap_impl!(
            buf0,
            basis,
            measure_anyon_model,
            measurement_strength,
            current_state,
            site,
            false,
        )
        p0 = sum(abs2, buf0)
        p1 = 1 - p0

        randomNumber = rand(rng)
        verbose && @show randomNumber
        if randomNumber < p0
            sample[i] = false
            buf0 .*= inv(sqrt(p0))
            current_state, buf0 = buf0, current_state
            F_layer += -log(p0)
            verbose && @show -log(p0)
        else
            # Compute 1-branch only when needed
            fill!(buf1, zero(T))
            _measuremap_impl!(
                buf1,
                basis,
                measure_anyon_model,
                measurement_strength,
                current_state,
                site,
                true,
            )
            sample[i] = true
            buf1 .*= inv(sqrt(p1))
            current_state, buf1 = buf1, current_state
            F_layer += -log(p1)
            verbose && @show -log(p1)
        end
    end
    return Measurement_outcome_boundary(current_state, sample, Float32(F_layer))
end

"""
    bulk_evolution(model::AnyonModel, state::Vector{ET}, measure_config::MeasureConfig,
                   samples::Union{Nothing,BitMatrix}=nothing)
    bulk_evolution(model::AnyonModel, sites, state::MPS, measure_config::MeasureConfig,
                   samples::Union{Nothing,BitMatrix}=nothing;
                   cutoff::Float64=1e-10, maxdim::Int=100)

Perform bulk measurement evolution from t₁ to t₂ on quantum state.

# Methods

## Vector version (from `Measurement.jl`)
Evolve a state vector under bulk measurements.

### Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `state::Vector{ET}`: Initial quantum state vector
- `measure_config::MeasureConfig`: Configuration struct containing τ, t₁, t₂, mode, rng, etc.
- `samples::Union{Nothing,BitMatrix}=nothing`: Predefined measurement sequences for `:sample` mode

### Returns
- `Measurement_outcome_bulk`: A struct containing:
  - `state::ET`: Final state after evolution
  - `samples::BitMatrix`: Measurement outcome sequences
  - `free_energys::Vector{Float32}`: Free energy for each layer
  - `entanglement_entropys::Vector{Float32}`: Half-chain EE at each period
  - `y_expectation_values::Vector{Float32}`: Normalized `Y` expectation after each
    period, or an empty vector when `track_y_expectation=false`

## MPS version (from `MPSMeasurement.jl`)
Evolve an MPS state under bulk measurements.

### Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `sites`: ITensor site indices
- `state::MPS`: Initial MPS quantum state
- `measure_config::MeasureConfig`: Configuration struct containing τ, t₁, t₂, mode, rng, etc.
- `samples::Union{Nothing,BitMatrix}=nothing`: Predefined measurement sequences for `:sample` mode
- `cutoff::Float64=1e-10`: MPS truncation cutoff
- `maxdim::Int=100`: Maximum bond dimension

### Returns
- `Measurement_outcome_mps_bulk`: A struct containing:
  - `state::MPS`: Final MPS state after evolution
  - `samples::BitMatrix`: Measurement outcome sequences
  - `free_energys::Vector{Float32}`: Free energy for each layer
  - `entanglement_entropys::Vector{Float32}`: Half-chain EE at each period

# Notes
- In `:Born` mode, samples are generated probabilistically via Born rule
- In `:sample` mode, `samples` must be provided as input
- `track_y_expectation=true` is currently supported by the exact-state Fibonacci
  evolution. The operator is constructed once per call and is not constructed when
  tracking is disabled.
- (2N+1) layers of measurements correspond to N time steps of evolution
"""
function bulk_evolution(
    anyon_model::AnyonModel{AT},   # DRY: don't repeat yourself.
    state::Vector{ET},
    measure_config::MeasureConfig,
    samples::Union{Nothing,BitMatrix} = nothing,
    normalized::Bool = true,
) where {ET,AT<:AbstractAnyonBasis}
    # ---------- Sample decided according to mode ----------
    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    current_state = copy(state)
    if mode == :Born
        return _born_measure(anyon_model, current_state, measure_config)
    else  # mode == :sample
        return _sample_measure(anyon_model, current_state, samples, measure_config, normalized)
    end
end

function _tracked_y_operator(model::AnyonModel, enabled::Bool)
    enabled || return nothing
    return _tracked_y_operator(model)
end

_tracked_y_operator(model::AnyonModel{FibonacciAnyon}) =
    topological_charge_operator(model)

function _tracked_y_operator(model::AnyonModel)
    error("Y expectation tracking is only supported for Fibonacci anyon models")
end

function _normalized_y_expectation(Y, state::AbstractVector)
    state_norm_squared = real(dot(state, state))
    state_norm_squared > 0 || error("cannot evaluate Y expectation for a zero-norm state")
    return Float32(real(dot(state, Y * state)) / state_norm_squared)
end

function _born_measure(
    model::AnyonModel{AT},
    current_state::Vector{ET},
    measure_config::MeasureConfig,
) where {AT,ET}

    n_cols = _samples_per_layer(model)  # Use max samples per layer
    Δt = measure_config.t₂ - measure_config.t₁ + 1
    τ = measure_config.τ
    enable_τ_eff = measure_config.enable_τ_eff
    rng = measure_config.rng
    verbose = measure_config.verbose
    Δt >= 0 || error("t₂ must be >= t₁")

    n_layers = layers_per_period(model)
    D = Δt * n_layers  # total number of layers

    # 1. Initialize sample matrix with max columns per layer
    samples = BitMatrix(zeros(Bool, D, n_cols))
    sample_free_energy = zeros(Float32, D)
    N = model.N
    entanglement_entropys = zeros(Float32, Δt)
    Y = _tracked_y_operator(model, measure_config.track_y_expectation)
    y_expectation_values =
        measure_config.track_y_expectation ? zeros(Float32, Δt) : Float32[]

    for period = 1:Δt
        # Apply all layers in this period
        for layer = 1:n_layers
            global_layer_idx = (period - 1) * n_layers + layer
            # Apply τ_eff only on the last layer of the last period
            τ_current = (period == Δt && layer == n_layers && enable_τ_eff) ? τ/2 : τ

            outcome = _stochastic_measurement_layer(
                model,
                τ_current,
                current_state;
                layer_idx = global_layer_idx,
                rng = rng,
                verbose = verbose,
            )
            current_state = outcome.state

            # Write samples to correct column indices for this layer
            col_indices = _get_sample_column_indices(model, global_layer_idx)
            samples[global_layer_idx, col_indices] = outcome.sample
            sample_free_energy[global_layer_idx] = outcome.free_energy
        end
        # Compute half-chain EE on-the-fly
        entanglement_entropys[period] =
            Float32(ee(anyon_rdm(model, collect(1:div(N, 2)), current_state)))
        if Y !== nothing
            y_expectation_values[period] = _normalized_y_expectation(Y, current_state)
        end
    end

    return Measurement_outcome_bulk(
        current_state,
        samples,
        sample_free_energy,
        entanglement_entropys,
        y_expectation_values,
    )
end

function _sample_measure(
    model::AnyonModel{AT},
    current_state::Vector{ET},
    samples::BitMatrix,
    measure_config::MeasureConfig,
    normalized::Bool = true,
) where {AT,ET}

    n_cols = _samples_per_layer(model)  # Use max samples per layer
    Δt = measure_config.t₂ - measure_config.t₁ + 1
    τ = measure_config.τ
    enable_τ_eff = measure_config.enable_τ_eff
    Δt >= 0 || error("t₂ must be >= t₁")

    n_layers = layers_per_period(model)
    D = Δt * n_layers  # total number of layers

    sample_free_energy = zeros(Float32, D)
    N = model.N
    entanglement_entropys = zeros(Float32, Δt)
    Y = _tracked_y_operator(model, measure_config.track_y_expectation)
    y_expectation_values =
        measure_config.track_y_expectation ? zeros(Float32, Δt) : Float32[]

    # 2. Validate sample matrix dimensions
    size(samples) == (D, n_cols) ||
        error("sample size should be ($D, $n_cols), got $(size(samples))")

    # 3. Deterministic trajectory for modes :sample
    #  Fibonacci: 2-layer period (even sites, odd sites)
    #  Ising (λ=0): 2-layer period (ZZ, X)
    #  OBF (λ≠0):   8-layer period (ZZ, X, OBF, X)

    #  If measure_operator is :Fibo, the measurement sites are half of N, circuits belike:
    #   1   1   1   1   1   1   1   1   1
    #     1   1   1   1   1   1   1   1    
    #   1   1   1   1   1   1   1   1   1
    #   -τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ- (head tail concatenation)

    # If measure_operator is :X or :ZZ, the measurement sites are N, circuits belike:
    #   Z₁Z₂ Z₂Z₃ Z₃Z₄ Z₄Z₅ Z₅Z₆ Z₆Z₇ Z₇Z₈ Z₈Z₁ (head tail concatenation)
    #  X    X    X    X    X    X    X    X
    #   Z₁Z₂ Z₂Z₃ Z₃Z₄ Z₄Z₅ Z₅Z₆ Z₆Z₇ Z₇Z₈ Z₈Z₁
    #  X    X    X    X    X    X    X    X
    #  ↑    ↑    ↑    ↑    ↑    ↑    ↑    ↑
    #  Or in majorana representation:
    # ---- --------  --------  --------  --------  --------  --------  --------  -----
    #  Z | |  ZZ  |  |  ZZ  |     ZZ  |  |  ZZ  |  |  ZZ  |  |  ZZ  |  |  ZZ  |  | Z
    # ---- --------  --------  --------  --------  --------  --------  --------  -----
    #  -------   -------   -------   -------   -------   -------   -------   ------
    #  |  X  |   |  X  |   |  X  |   |  X  |   |  X  |   |  X  |   |  X  |   |  X  |
    #  -------   -------   -------   -------   -------   -------   -------   ------
    #  γ₁   γ₂   γ₃   γ₄   γ₅   γ₆   γ₇   γ₈   γ₉  γ₁₀  γ₁₁  γ₁₂  γ₁₃  γ₁₄  γ₁₅  γ₁₆

    for period = 1:Δt
        # √M₁ᵒ M₁ᵉ √M₁ᵒ √M₁ᵒ M₁ᵉ √M₁ᵒ ⋯ √M₁ᵒ M₁ᵉ √M₁ᵒ→ M₁ᵉ M₁ᵒ M₁ᵉ M₁ᵒ ⋯ M₁ᵉ √M₁ᵒ. 
        # √X ZZ √X √X ZZ √X ⋯ √X ZZ √X→ √X ZZ X ZZ ⋯ X ZZ √X. To ensure each layer is hermitian, first layer doesn't matter.
        # Or √ZZ X √ZZ √ZZ X √ZZ ⋯ X √ZZ X √ZZ→ X ZZ X ZZ ⋯ X √ZZ, also works (we choose this one here).
        for layer = 1:n_layers
            global_layer_idx = (period - 1) * n_layers + layer
            # Apply τ_eff only on the last layer of the last period
            τ_current = (period == Δt && layer == n_layers && enable_τ_eff) ? τ/2 : τ

            # Read samples from correct column indices for this layer
            col_indices = _get_sample_column_indices(model, global_layer_idx)
            layer_sample = BitVector(samples[global_layer_idx, col_indices])

            outcome = _apply_measurement_layer(
                model,
                τ_current,
                current_state,
                layer_sample;
                layer_idx = global_layer_idx,
                # ---------------------------------------------------------------------------
                # Unnormalized sample evolution: apply a fixed sample without normalizing.
                # This gives a linear map T(s) acting on the state vector.
                # ---------------------------------------------------------------------------
                normalized,
            )
            current_state = outcome.state
            sample_free_energy[global_layer_idx] = outcome.free_energy
        end
        # Compute half-chain EE on-the-fly
        entanglement_entropys[period] =
            Float32(ee(anyon_rdm(model, collect(1:div(N, 2)), current_state)))
        if Y !== nothing
            y_expectation_values[period] = _normalized_y_expectation(Y, current_state)
        end
    end
    return Measurement_outcome_bulk(
        current_state,
        samples,
        sample_free_energy,
        entanglement_entropys,
        y_expectation_values,
    )
end


"""
    bayes_distort(γ::Float64, trajectories::Vector{Int64}, probabilities::Vector{Float64})

Distort measurement trajectories based on a Bayesian distortion factor γ.

This function implements the distortion process where each faithful sample s 
is converted to a distorted sample s̃ according to the conditional probability:
    P(s̃|s) = ∏ⱼ (1 + γ s̃ⱼ sⱼ)/2

Only works for single-layer measurements to generate new samples based on 
the projective limit measurement.

# Arguments
- `γ::Float64`: Distortion factor (readout fidelity parameter, 0 ≤ γ ≤ 1)
- `trajectories::Vector{Int64}`: Vector of measurement trajectories
- `probabilities::Vector{Float64}`: Corresponding probabilities for each trajectory

# Returns
- `Tuple{Vector{Vector{Int64}}, Vector{Float64}}`:
  (distorted_trajectories, distorted_probabilities) in corresponding order
"""
function bayes_distort(
    γ::Float64,
    trajectories::Vector{Int64},
    probabilities::Vector{Float64},
)

    # Dictionary to store the distorted trajectory probabilities
    distorted_prob_dict = Dict{Vector{Int64},Float64}()
    n_sites = length(trajectories)
    distorted_prob = Vector{Vector{Float64}}(undef, n_sites)
    transfer_matrix = [1 + γ 1 - γ; 1 - γ 1 + γ] / 2

    # For each original trajectory
    for (traj_idx, original_traj) in enumerate(trajectories)
        original_prob = probabilities[traj_idx]
        prob_distribution =
            (trajectories[traj_idx] == 1) ? [original_prob, 1 - original_prob] :
            [1 - original_prob, original_prob]
        distorted_prob[traj_idx] = transfer_matrix * prob_distribution
    end

    # Generate all possible distorted trajectories (2^n possibilities)
    for distorted_bits = 0:(2^n_sites-1)
        # Convert bit representation to ±1 trajectory
        prob = 1.0
        distorted_traj = Vector{Int64}(undef, n_sites)
        for j = 1:n_sites
            # Extract j-th bit and convert to ±1
            bit = (distorted_bits >> (j-1)) & 1
            distorted_traj[j] = bit
            prob *= distorted_prob[j][(bit==1) ? 1 : 2]  # bit + 1 because Julia is 1-indexed
        end



        if haskey(distorted_prob_dict, distorted_traj)
            distorted_prob_dict[distorted_traj] = prob
        else
            distorted_prob_dict[distorted_traj] = prob
        end
    end
    # Convert dictionary to vectors
    distorted_trajectories = collect(keys(distorted_prob_dict))
    distorted_probabilities = collect(values(distorted_prob_dict))

    return distorted_trajectories, distorted_probabilities
end

"""
    transfer_matrix(model::AnyonModel, τ::Float64, sample::BitMatrix)

Construct the transfer matrix for a fixed measurement trajectory.
Applies `_apply_measurement_layer` with `normalized=false` to each basis vector.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `sample::BitMatrix`: Measurement outcome sequences (rows = layers, cols = sites)

# Returns
- `Matrix{Float64}`: Transfer matrix where column i = T(sample) * e_i
"""
function transfer_matrix(
    model::AnyonModel{AT},
    τ::Float64,
    sample::BitMatrix,
) where {AT<:AbstractAnyonBasis}
    basis = anyon_basis(model)
    l = length(basis)
    TM = zeros(Float64, l, l)
    n_layers = layers_per_period(model)
    
    for i = 1:l
        st = zeros(Float64, l)
        st[i] = 1.0
        for layer_idx in 1:n_layers
            col_indices = _get_sample_column_indices(model, layer_idx)
            layer_sample = BitVector(sample[layer_idx, col_indices])
            out = _apply_measurement_layer(
                model, τ, st, layer_sample;
                layer_idx = layer_idx, normalized = false,
            )
            st = out.state
        end
        TM[:, i] = st
    end

    return TM
end

"""
    transfer_matrix_dynamics(model::AnyonModel, τ::Float64, sample::BitMatrix; n_spectrums::Int=10)

Compute the exact eigenvalue spectrum of the **cumulative** transfer matrix
T₁, T₁T₂, T₁T₂T₃, … at each time step via full exact diagonalization (ED).

At step `t` the cumulative transfer matrix is
    T_cum(t) = T_t ⋯ T₂ T₁
where each Tᵢ is the local transfer matrix for time slice `i`. The function
diagonalizes T_cum(t) exactly and records the dominant eigenvalues.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `sample::BitMatrix`: Measurement outcome sequences (rows = layers, cols = sites)
- `n_spectrums::Int=10`: Number of dominant eigenvalues to keep per step

# Returns
- `Matrix{ComplexF64}`: Array of shape `(n_spectrums, n_steps)` where column `t`
  contains the largest `n_spectrums` eigenvalues of the cumulative transfer
  matrix up to that step.
"""
function transfer_matrix_dynamics(
    model::AnyonModel{AT},
    τ::Float64,
    sample::BitMatrix;
    n_spectrums::Int = 10,
) where {AT<:AbstractAnyonBasis}

    n_layers = layers_per_period(model)
    D_layers, n_cols = size(sample)
    @assert D_layers % n_layers == 0 "Number of layers $D_layers must be divisible by $n_layers"
    t = D_layers ÷ n_layers
    n_cols == _samples_per_layer(model) ||
        error("sample size spatial dimension must be $(_samples_per_layer(model)), got $n_cols")

    basis = anyon_basis(model)
    l = length(basis)
    k = min(n_spectrums, l)

    spectrum_tlis = zeros(ComplexF64, k, t)
    TM_cum = Matrix{Float64}(I, l, l)

    for step in 1:t
        sample_layer = sample[(step - 1) * n_layers + 1 : step * n_layers, :]

        # Construct the local transfer matrix for this step
        TM_step = zeros(Float64, l, l)
        for i in 1:l
            st = zeros(Float64, l)
            st[i] = 1.0
            # Note that in the time end, we didn't use the effective τ, to ensure total transfer matrix is hermitian, as in Born case it born out to be non-hermitian.
            config = MeasureConfig(τ = τ, mode = :sample, t₂ = 1, enable_τ_eff = false)
            TM_step[:, i] = _sample_measure(
                model,
                st,
                sample_layer,
                config,
                false;
            ).state
        end

        # Cumulative product: T_cum = T_step * T_cum
        TM_cum = TM_step * TM_cum

        # Exact diagonalization of the cumulative transfer matrix
        energy = eigvals(TM_cum)
        sorted_energy = sort(energy, by = abs, rev = true)
        spectrum_tlis[:, step] = sorted_energy[1:k]
    end

    return spectrum_tlis
end


"""
    transfer_matrix_subspace(model::AnyonModel, τ::Float64, sample::BitMatrix; n_states::Int=10)

Compute the dominant spectrum of the transfer matrix via subspace iteration.

The algorithm initializes `n_states` product states (basis vectors), then
iteratively applies the transfer matrix for each time slice'''s measurement outcome.
At each step, states are normalized and QR-orthogonalized to keep the subspace
well-conditioned. Finally, the last transfer matrix is projected onto the
converged subspace (Rayleigh-Ritz) to obtain Ritz values, which converge to
the true eigenvalues.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `sample::BitMatrix`: Measurement outcome sequences (rows = layers, cols = sites).
  The number of rows must be divisible by `layers_per_period(model)`.
- `n_states::Int=10`: Number of initial basis vectors to propagate (and Ritz values to compute)

# Returns
- `Vector{Float64}`: Sorted `-log.(abs.(ritz_values))`.

# Examples
```jldoctest
julia> using FibonacciChain, LinearAlgebra

julia> L = 8; τ = atanh(0.95);

julia> model = AnyonModel(FibonacciAnyon(), L; pbc = true);

julia> sample = BitMatrix(ones(Int8, 2, div(L, 2)));

julia> spectrum = transfer_matrix_subspace(model, τ, sample; n_states = 5);

julia> length(spectrum) == 5
true
```
"""
function transfer_matrix_subspace(
    model::AnyonModel{AT},
    τ::Float64,
    sample::BitMatrix;
    n_states::Int = 10,
) where {AT<:AbstractAnyonBasis}
    # Here the transfer matrix is not hermitian, thus the Schur vector is not eigenvectors. We need to do Rayleigh-Ritz projection. But when the non-hermitian matrix is too ill-conditioned, this method fails.
    n_layers = layers_per_period(model)
    D_layers, n_cols = size(sample)
    @assert D_layers % n_layers == 0 "Number of layers $D_layers must be divisible by $n_layers"
    t = D_layers ÷ n_layers
    n_cols == _samples_per_layer(model) ||
        error("sample size spatial dimension must be $(_samples_per_layer(model)), got $n_cols")

    basis = anyon_basis(model)
    l = length(basis)
    k = min(n_states, l)

    # Initialize k product states (basis vectors)
    states = zeros(Float64, l, k)
    for i in 1:k
        states[i, i] = 1.0
    end
    spectrum_tlis = zeros(k, t)

    for step in 1:t
        sample_layer = sample[(step - 1) * n_layers + 1 : step * n_layers, :]
        for i in 1:k
            config = MeasureConfig(τ = τ, mode = :sample, t₂ = 1, enable_τ_eff = false)
            outcome = _sample_measure(
                model,
                states[:, i],
                sample_layer,
                config,
                false,
            )
            states[:, i] = outcome.state
        end
        
        Q, R = qr(states)
        states = Q[:, 1:k]
        # Note here do not sort, will distort the Lyapunov spectrum (singular eigenvalues, corresponds to the lnZ, not lnp)
        spectrum_tlis[:, step] = -log.(abs.(diag(R)))
    end

    return spectrum_tlis
end
