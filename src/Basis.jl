function Fibonacci_chain_OBC(::Type{T}) where {N, T <: BitStr{N}}
    # Generate Fibonacci chain for Fibonacci model with open boundary condition
    fib_chain=[[T(0), T(1)],[T(0), T(1), T(2)]]
    for i in 3:N
        push!(fib_chain,vcat([s << 1 for s in fib_chain[i-1]],[(s << 2 | T(1)) for s in fib_chain[i-2]]))
    end
    # each push we add a bit"0" or bit"01" to the end of the bit_string, and finally return the last element of the fib_chain
    return fib_chain[N]
end

function Fibonacci_chain_PBC(::Type{T}) where {N, T <: BitStr{N}}
    # Generate Fibonacci chain  for Fibonacci model with periodic boundary condition
    return filter(c -> iszero((c >> (N-1)) & (c & 1)), Fibonacci_chain_OBC(T))
end

abstract type AbstractAnyonType end
struct FibonacciAnyon <: AbstractAnyonType end
struct IsingAnyon <: AbstractAnyonType end
struct OBFAnyon <: AbstractAnyonType end

"""
    AnyonModel(anyon_type::AbstractAnyonType, N::Int; pbc::Bool=true, measure_operator::Symbol=:Antiferro, kwargs...)

Represents a 1D anyon chain model.

# Fields
- `anyon_type::AbstractAnyonType`: The type of anyon, e.g., `FibonacciAnyon()` or `IsingAnyon()`.
- `N::Int`: The number of sites in the chain.
- `pbc::Bool`: A boolean indicating whether periodic boundary conditions are applied (`true`) or not (`false`).
- `measure_operator::Symbol`: The operator type used for defining the Hamiltonian, e.g., `:Antiferro` or `:Ferro` for Fibonacci anyons.
- `params::Dict{Symbol, Any}`: Additional model parameters (e.g., `J`, `h` for Ising model).

# Examples
```julia
# Fibonacci anyon model
model_fibo = AnyonModel(FibonacciAnyon(), 6; pbc=true, measure_operator=:Antiferro)

# Ising model with custom couplings
model_ising = AnyonModel(IsingAnyon(), 10; pbc=true, measure_operator=:X, J=1.0, h=0.5)

# Access parameters
model_ising.params[:J]  # returns 1.0
```
"""
struct AnyonModel{AT<:AbstractAnyonType}
    anyon_type::AT
    N::Int
    pbc::Bool
    measure_operator::Symbol
    params::Dict{Symbol, Any}
    
    function AnyonModel(anyon_type::AT, N::Int; pbc::Bool=true, measure_operator::Symbol=:Antiferro, kwargs...) where {AT<:AbstractAnyonType}
        @assert N > 0 "N is expected to be greater than 0, but got $N"
        params = Dict{Symbol, Any}(kwargs...)
        return new{AT}(anyon_type, N, pbc, measure_operator, params)
    end
end

# Helper function to get parameter with default value
get_interaction_param(model::AnyonModel, key::Symbol, default) = get(model.params, key, default)

"""
    anyon_basis(model::AnyonModel; symmetry_block=nothing)
    anyon_basis(anyon_type::AbstractAnyonType, ::Type{T}; pbc::Bool=true, symmetry_block=nothing) where {N, T <: BitStr{N}}

Generate basis states for an anyon chain model.

This function provides two interfaces: a high-level one that takes an `AnyonModel` object, and a low-level one that takes the model parameters directly.

# Arguments
## High-level interface
- `model::AnyonModel`: An `AnyonModel` object containing the system parameters (anyon type, size, and boundary conditions).
- `symmetry_block`: Optional. The topological charge sector to filter the basis. Not yet implemented.

## Low-level interface
- `anyon_type::AbstractAnyonType`: The type of anyon, e.g., `FibonacciAnyon()` or `IsingAnyon()`.
- `T::Type`: The `BitStr` type specifying the chain length `N`.
- `pbc::Bool=true`: Specifies whether to use periodic (`true`) or open (`false`) boundary conditions.
- `symmetry_block`: Optional. The topological charge sector.

# Returns
- `Vector{T}`: A sorted vector of the basis states.

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> model = AnyonModel(FibonacciAnyon(), 4; pbc=true); basis = anyon_basis(model); length(basis)
7

julia> basis_fibo = anyon_basis(FibonacciAnyon(), BitStr{4, Int}; pbc=true); length(basis_fibo)
7

julia> model_ising = AnyonModel(IsingAnyon(), 4; pbc=true); basis_ising = anyon_basis(model_ising); length(basis_ising)
16
```
"""
function anyon_basis(::FibonacciAnyon, ::Type{T}; pbc::Bool=true, symmetry_block=nothing) where {N, T <: BitStr{N}}
    # Generate basis for Fibonacci model, return BitBasis form, which can be used as binary and decimal form. Here we both consider PBC and OBC
    @assert symmetry_block === nothing || symmetry_block in [0, 1, :tau, :trivial] "symmetry_block is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $symmetry_block"
    @assert T <: BitStr{N} "Type T must be a BitStr type"
    # If pbc is true, use Fibonacci_chain_PBC, otherwise use Fibonacci_chain_OBC
    if pbc
        if T == BitStr{1, Int}
            basis=Fibonacci_chain_OBC(T)
        else
            basis=Fibonacci_chain_PBC(T)
        end
    else
        basis=Fibonacci_chain_OBC(T)
    end
    sorted_basis=sort(basis)

    if symmetry_block !== nothing
        # Filter basis by topological charge
        # sorted_basis = filter(s -> topological_charge(s) == symmetry_block, sorted_basis)
        error("Filtering by topological charge is not implemented yet.")
    end

    return sorted_basis
end

function anyon_basis(::IsingAnyon, ::Type{T}; pbc::Bool=true, symmetry_block=nothing) where {N, T <: BitStr{N}}
    return [T(i) for i in 0:(2^N - 1)]
end

function anyon_basis(::OBFAnyon, ::Type{T}; pbc::Bool=true, symmetry_block=nothing) where {N, T <: BitStr{N}}
    return [T(i) for i in 0:(2^N - 1)]
end

anyon_basis(model::AnyonModel; symmetry_block=nothing) = anyon_basis(model.anyon_type, BitStr{model.N, Int}; pbc=model.pbc, symmetry_block=symmetry_block)

"""
    Fsymmetry_coef(::FibonacciAnyon, ::Type{T}, state::T, base::T; pbc::Bool=true) where {N, T <: BitStr{N}}
    Fsymmetry_coef(model::AnyonModel{FibonacciAnyon}, state::T, base::T)

Compute topological symmetry coefficient for state in given base configurations for Fibonacci anyon chain.

# Arguments
## Low-level interface
- `::FibonacciAnyon`: Fibonacci anyon type instance
- `::Type{T}`: BitStr type specifying chain length N
- `state::T`: Target state configuration
- `base::T`: Base state configuration  
- `pbc::Bool=true`: Periodic boundary conditions

## High-level interface
- `model::AnyonModel{FibonacciAnyon}`: Anyon model containing system parameters
- `state::T`: Target state configuration
- `base::T`: Base state configuration

# Returns
- `Float64`: Topological symmetry coefficient based on Fibonacci fusion rules

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; T = BitStr{N, Int};

julia> model = AnyonModel(FibonacciAnyon(), N; pbc=true);

julia> state = T(0b1010); ϕ = (1 + sqrt(5)) / 2;

julia> base = T(0b0101);

julia> coef = Fsymmetry_coef(model, state, base); coef ≈ 0.3819660112501051
true
```
"""
function Fsymmetry_coef(::FibonacciAnyon, state::T, base::T; pbc::Bool=true) where {N, T <: BitStr{N}}
    # Defined as, where idxin idxbond idxout ∈ state, idxbond' ∈ base, in Anyon basis, not in Fibonacci chain basis.
    #  %%%%%%%%%%%% τ, idxin, τ         idxbond
    #  %%
    #  %% 
    #  %%%%%%%%%%%%%
    #  %%
    #  %%
    #  %%%          idxout              idxbond'
    #  or effectively, the coefficient like: A_{x_1 x_2 x_2'}^{x_1'} 
    ϕ = (1+√5)/2
    prod=1

    @assert pbc || site != N "For OBC, site must not be $N, but got $site"
        for site in 1:N-1
            # Identify {x_2}, if x_2' is 1, return 1, otherwise, check x_3', if x_3' is 1, return 1.
            if state[N - site + 1] == 1
                prod *= 1
            else
                if base[N - site] == 1
                    prod *= 1
                else
                    # check x_3, if x_3 is 1
                    if state[N - site] == 0
                        # Determine coef according to x_2'
                        prod *= (base[N - site + 1] == 0) ? -ϕ^(-1) : ϕ^(-1/2)
                    else
                        prod *= (base[N - site + 1] == 0) ? ϕ^(-1/2) : ϕ^(-1)
                    end
                end
            end
        end

        if pbc
            # Check x_N', if x_N' is 1, return 1, otherwise, check x_1', if x_1' is 1, return 1.
            if state[1] == 1
                prod *= 1
            else
                if base[N] == 1
                    prod *= 1
                else
                    # check x_1, if x_1 is 1
                    if state[N] == 0
                        # Determine coef according to x_N'
                        prod *= (base[1] == 0) ? -ϕ^(-1) : ϕ^(-1/2)
                    else
                        prod *= (base[1] == 0) ? ϕ^(-1/2) : ϕ^(-1)
                    end
                end
            
                prod *= 1
            end
        end

        return prod
end
Fsymmetry_coef(model::AnyonModel{FibonacciAnyon}, state::T, base::T) where {N, T <: BitStr{N}} = 
    Fsymmetry_coef(model.anyon_type, state, base; pbc=model.pbc)

"""
    topological_symmetry_basismap(model::AnyonModel{FibonacciAnyon}, state::T) where {N, T <: BitStr{N}}

Compute topological symmetry coefficients for all basis states relative to given state.

# Arguments
- `model::AnyonModel{FibonacciAnyon}`: Fibonacci anyon model
- `state::T`: Reference state configuration

# Returns
- `Vector{Float64}`: Coefficients for each basis state

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; T = BitStr{N, Int};

julia> model = AnyonModel(FibonacciAnyon(), N; pbc=true);

julia> state = T(0b0000);  # Vacuum state

julia> coeffs = topological_symmetry_basismap(model, state);

julia> length(coeffs) == length(anyon_basis(model))
true
```
"""
function topological_symmetry_basismap(model::AnyonModel{AT}, state::T) where {N, T <: BitStr{N}, AT<:FibonacciAnyon}
    # Compute the topological symmetry map for a given state using the topological symmetry site map for all site
    # However, it seems not easy to separate the state into different topological symmetric sectors.
    basis = anyon_basis(model)
    coeflis = Vector{Float64}(undef, length(basis))
    
    # For each base in basis, check the state at each site
    for (idx, base) in enumerate(basis)
        coef = Fsymmetry_coef(model, state, base)
        coeflis[idx] = coef
    end
    return coeflis
end

"""
    topological_charge_operator(model::AnyonModel{FibonacciAnyon}, ::Type{T}) where {N, T <: BitStr{N}}

Compute the topological charge operator matrix for Fibonacci anyon model.

# Arguments
- `model::AnyonModel{FibonacciAnyon}`: Fibonacci anyon model
- `T::Type`: BitStr type specifying chain length N

# Returns
- `Matrix{Float64}`: Topological charge operator matrix (Y matrix)

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; T = BitStr{N, Int};

julia> model = AnyonModel(FibonacciAnyon(), N; pbc=true);

julia> Y = topological_charge_operator(model, T);

julia> size(Y) == (length(anyon_basis(model)), length(anyon_basis(model)))
true
```
"""
function topological_charge_operator(model::AnyonModel{AT}, ::Type{T}) where {N, T <: BitStr{N}, AT<:FibonacciAnyon}
    # compute the topological charge operator Yl in the Fibonacci model. default l=0, for tau. l=1, for vacuum.

    basis=anyon_basis(model)
    Ymatrix=hcat(topological_symmetry_basismap.(Ref(model), basis)...)

    return Ymatrix
end

"""
    anyon_basis(AT::AbstractAnyonType, ::Type{T}, k::Int; symmetry_block=nothing) where {N, T <: BitStr{N}}
    anyon_basis(model::AnyonModel, k::Int; symmetry_block=nothing)

Generate basis states in specific momentum sector `k` and topological sector `symmetry_block`.

# Arguments
## Low-level interface
- `AT::AbstractAnyonType`: Anyon type instance (e.g., `FibonacciAnyon()`, `IsingAnyon()`)
- `T::Type`: BitStr type specifying chain length N
- `k::Int`: Momentum sector (0 ≤ k ≤ N-1)
- `symmetry_block`: Topological charge sector (optional)

## High-level interface  
- `model::AnyonModel`: Anyon model containing system parameters
- `k::Int`: Momentum sector (0 ≤ k ≤ N-1)
- `symmetry_block`: Topological charge sector (optional)

# Returns
- `Vector{T}`: Basis states in momentum sector k
- `Dict{T, Vector{T}}`: Representative mapping for translation equivalence classes

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> model = AnyonModel(FibonacciAnyon(), 6; pbc=true);

julia> basisK, basis_dic = anyon_basis(model, 0);

julia> length(basisK) > 0
true
```
"""
function anyon_basis(AT::AbstractAnyonType, ::Type{T}, k::Int; symmetry_block=nothing) where {N, T <: BitStr{N}}
    #params: a int of lattice number, momentum of system, topological_charge symmetry_block, which default to be nothing
    #return: computational basis in given momentum kinetically constrained subspace with decimal int form in golden chain model
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert symmetry_block === nothing || symmetry_block in [0, 1, :tau, :trivial] "symmetry_block is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $symmetry_block"

    basisK = Vector{T}(undef, 0)
    basis = anyon_basis(AT, T; symmetry_block=symmetry_block)
    basis_dic = Dict{T, Vector{T}}()

    for i in basis
        category = get_representative(i)[1]
        if haskey(basis_dic, category)
            push!(basis_dic[category], i)
        else
            basis_dic[category] = [i]
        end
    end

    for j in eachindex(basis)
        n=basis[j]
        RS = get_representative(n)[1]
        if RS == n && (k * length(basis_dic[RS])) % N == 0
            push!(basisK, n)
        end
    end

    return basisK, basis_dic
end
anyon_basis(model::AnyonModel, k::Int; symmetry_block=nothing) = anyon_basis(model.anyon_type, BitStr{model.N, Int}, k; symmetry_block=symmetry_block)

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
function Fibomap(state::T, i::Int; ferro::Bool=false) where {N, T <: BitStr{N}}
    ϕ = (1 + √5) / 2  # Golden ratio
    fl = bmask(T, N)
    X_state = flip(state, fl >> (i-1))
    
    # Read bit at site i: 0 = vacuum (1), 1 = τ anyon
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
antimap(state::T, i::Int) where {N, T <: BitStr{N}} = Fibomap(state, i; ferro=false)
ferromap(state::T, i::Int) where {N, T <: BitStr{N}} = Fibomap(state, i; ferro=true)

function Isingmap(state::T, i::Int, pbc::Bool=true; kwargs...) where {N, T <: BitStr{N}}
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
    OBFmap(state::T, i::Int, pbc::Bool=true; λ=1.0) where {N, T <: BitStr{N}}

Apply O'Brien-Fendley terms (X_i Z_{i} Z_{i+2} + Z_i Z_{i+1} X_{i+2}) at site i.

Returns (output_states, weights) for the XZZ + ZZX terms in the OBF Hamiltonian.
For OBC, valid range is 1 ≤ i ≤ N-2.
For PBC, valid range is 1 ≤ i ≤ N with periodic wrapping.
"""
function OBFmap(state::T, i::Int, pbc::Bool=true; λ::Float64=1.0) where {N, T <: BitStr{N}}
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
    xzz_weight = λ * zz_i1i2  # λ for Hamiltonian (energy lowering)
    
    # ZZX term: Z_i Z_{i+1} X_{i+2}
    zz_ii1 = (bit_i == bit_i1) ? 1 : -1
    zzx_state = X(state, i2)
    zzx_weight = λ * zz_ii1
    
    return xzz_state, zzx_state, xzz_weight, zzx_weight
end

function count_subBitStr(state::T) where {N, T <: BitStr{N}}
    n = length(state)
    n < 3 && return 0 
    
    str100, str101, str001 = T(4), T(5), T(1) # 100, 101, 001
    num=0
    
    mask=bmask(T, 1, 2, 3)
    for i in 1:(n-2) # start from string right to left
        substr = state & (mask << (i-1))  
        if substr == str101
            num+= 1
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
function actingHam(model::AnyonModel{FibonacciAnyon}, state::T) where {N, T <: BitStr{N}}
    @assert num_digits(T) == N "The length of system is expected to be $N, but got $(num_digits(T))"
    
    pbc = model.pbc
    ferro = (model.measure_operator == :Ferro)
    sign = ferro ? 1 : -1
    
    # Helper to apply Fibomap and accumulate results
    function apply_fibomap!(output::Dict{T, Float64}, i::Int)
        state1, state2, weight1, weight2 = Fibomap(state, i; ferro=ferro)
        output[state1] = get(output, state1, 0.0) + weight1
        output[state2] = get(output, state2, 0.0) + weight2
    end
    
    output = Dict{T, Float64}()
    
    # Diagonal contribution from 101, 100, 001 patterns
    # Sign depends on ferro/antiferro: H = sign * ∑ (fusion constraints)
    output[state] = get(output, state, 0.0) + sign * count_subBitStr(state)
    
    # Bulk terms: sites 2 to N-1, apply Fibomap where 0x0 pattern exists
    mask = bmask(T, N, N-2)  # Mask to check neighbors
    for i in 2:N-1 
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

function actingHam(model::AnyonModel{IsingAnyon}, state::T) where {N, T <: BitStr{N}}
    @assert num_digits(T) == N "The length of system is expected to be $N, but got $(num_digits(T))"
    
    pbc = model.pbc
    J = get_interaction_param(model, :J, 1.0)
    h = get_interaction_param(model, :h, 1.0)
    fl = bmask(T, N)
    X(st, j) = flip(st, fl >> (j-1))
    
    output = Dict{T, Float64}()
    
    # ZZ terms: -J ∑ Z_i Z_{i+1}
    n_bonds = pbc ? N : N - 1
    for i in 1:n_bonds
        i1 = mod1(i + 1, N)
        zz_val = (readbit(state, N-i+1) == readbit(state, N-i1+1)) ? -J : J
        output[state] = get(output, state, 0.0) + zz_val
    end
    
    # X terms: -h ∑ X_i
    for i in 1:N
        output[X(state, i)] = get(output, X(state, i), 0.0) - h
    end

    return output
end

function actingHam(model::AnyonModel{OBFAnyon}, state::T) where {N, T <: BitStr{N}}
    @assert num_digits(T) == N "The length of system is expected to be $N, but got $(num_digits(T))"
    
    pbc = model.pbc
    λ = get_interaction_param(model, :λ, 1.0)  # OBF coupling strength

    # Generate OBF model Hamiltonian: H = λ ∑ (X_i Z_{i+1} Z_{i+2} + Z_i Z_{i+1} X_{i+2}) - X - ZZ
    output = Dict{T, Float64}()
    for i in 1:N
        s1, s2, w1, w2 = Isingmap(state, i, pbc)
        output[s1] = get(output, s1, 0.0) + w1
        output[s2] = get(output, s2, 0.0) + w2
        state1, state2, weight1, weight2 = OBFmap(state, i, pbc; λ=λ)
        output[state1] = get(output, state1, 0.0) + weight1
        output[state2] = get(output, state2, 0.0) + weight2
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

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; model_fibo = AnyonModel(FibonacciAnyon(), N; pbc=true, measure_operator=:Antiferro);

julia> H_fibo = anyon_ham(model_fibo); size(H_fibo)
(7, 7)

julia> model_ising = AnyonModel(IsingAnyon(), N; pbc=true, J=1.0, h=1.0); H_ising = anyon_ham(model_ising); size(H_ising)
(16, 16)
```
"""
function anyon_ham(model::AnyonModel{AT}) where {AT<:AbstractAnyonType}
    basis = anyon_basis(model)
    l = length(basis)
    H = zeros(Float64, (l, l))
    
    for i in 1:l
        output = actingHam(model, basis[i])
        for (m, weight) in output
            j = searchsortedfirst(basis, m)
            H[i, j] += weight
        end
    end

    return H
end
# Another method to write Fibonacci Hamiltonian is using the Measurement operator sum. For example, H = -∑ X_i, where X_i is the Temperley-Lieb generator acting on site i-1, i, and i+1. Pilis = [FibonacciChain.measure_matrix(BitStr{16, Int}, 1000.0, idx, 0) for idx in 1:N]. H = -sum(Pilis). This two Hamiltonian difference is not a constant, but like a arc in conformal energy spectrum below arc, but they have the same eigenstates.

function cyclebits(state::T) where {N, T <: BitStr{N}}
    #params: t is an integer, N is the length of the binary string
    #We also use this order: system size, state, circular shift bitstring 1 bit.
    # In case need to shift more than 1 bit, we can use a loop or recursion. or we leave a interface here  n_translations::Int
    mask = 1 << N - 1
    return ((state << 1) | (state >> (N - 1))) & mask
end

function get_representative(state::T) where {N, T <: BitStr{N}}
#Finds representative and representative translation for a state.
#State should be a decimal integer.

    representative = state
    translation = 0
    # cycle(bits) = (bits << 1) % (2^N - 1)  # Left shift the state by one position
    current = state
    for n_translation_sites in 1:N-1
        current = cyclebits(current)  # Cycle the bits
        if current < representative
            representative = current
            translation = n_translation_sites
        end
    end

    return representative, translation
end

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
function anyon_ham(model::AnyonModel{AT}, ::Type{T}, k::Int; symmetry_block=nothing) where {N, T <: BitStr{N}, AT <: AbstractAnyonType}
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert symmetry_block === nothing || symmetry_block in [0, 1, :tau, :trivial] "symmetry_block is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $symmetry_block"

    basisK, basis_dic = anyon_basis(model, k, symmetry_block=symmetry_block)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    H = zeros(ComplexF64, (l, l))

    for i in 1:l
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

anyon_ham(model::AnyonModel{AT}, k::Int; symmetry_block=nothing) where {AT <: AbstractAnyonType} = 
    anyon_ham(model, BitStr{model.N, Int}, k; symmetry_block=symmetry_block)

# join two lists of basis by make a product of two lists, b is placed after a (counts from left to right)
function process_join(a, b)
    # effectively, vec([join(a, b) for a in a for b in b]) # seems faster
    return [join(a, b) for a in a for b in b]
end

# create Fibonacci basis composed of multiple disjoint sub-chains
function joint_basis(AT::AbstractAnyonType, lengthlis::Vector{Int}, pbc::Bool=false;)
    # The element order in lengthlis doesn't matter, i.e., [2, 3] is only different with [3, 2] in basis order. Only your input state must consistent with the basis order.
    if isempty(lengthlis)
        return BitStr{0, Int}[]
    else
        return sort(mapreduce(len -> anyon_basis(AT, BitStr{len, Int}; pbc=pbc), process_join, lengthlis))
    end
end

function connected_components(v::Vector{Int})
    if isempty(v)
        return Int[]
    end

    sort!(v)

    result = []
    current_segment = [v[1]]

    for i in 2:length(v)
        if v[i] == v[i - 1] + 1
            push!(current_segment, v[i])
        else
            push!(result, current_segment)
            current_segment = [v[i]]
        end
    end

    push!(result, current_segment)

    return result
end

function move_subsystem(::Type{BitStr{M, INT}}, basis::BitStr{N, INT}, subsystems::Vector{Int}) where {M, N, INT}
    @assert length(subsystems) == N "subsystems length is expected to be $N, but got $(length(subsystems))"
    @assert M >= N "total length is expected to be greater than or equal to $N, but got $M"
    return sum(i -> BitStr{M}(readbit(basis.buf, i) << (M - subsystems[N-i+1])), 1:N)
end

# take environment part of a basis
takeenviron(x, mask::BitStr{l}) where {l} = x & (~mask)
# take system part of a basis
takesystem(x, mask::BitStr{l}) where {l} = (x & mask)

"""
    anyon_rdm(AT::AbstractAnyonType, ::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}; pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    anyon_rdm(AT::AbstractAnyonType, ::Type{T}, subsystems::Vector{Int64}, state::Matrix{ET}; pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    anyon_rdm(model::AnyonModel, subsystems::Vector{Int64}, state::Union{Vector{ET}, Matrix{ET}})

Compute the reduced density matrix (RDM) for a specified subsystem.

This function can operate on both state vectors and density matrices.

# Arguments
## High-level interface
- `model::AnyonModel`: An `AnyonModel` object containing the system parameters.
- `subsystems::Vector{Int64}`: A vector of site indices to keep in the reduced system. Indices are 1-based.
- `state::Union{Vector{ET}, Matrix{ET}}`: The quantum state, either as a state vector or a full density matrix.

## Low-level interface
- `AT::AbstractAnyonType`: The anyon type, e.g., `FibonacciAnyon()`.
- `T::Type`: The `BitStr` type specifying the total chain length `N`.
- `subsystems::Vector{Int64}`: The indices of the subsystem sites.
- `state::Union{Vector{ET}, Matrix{ET}}`: The quantum state vector or density matrix.
- `pbc::Bool=true`: Specifies whether the original system has periodic boundary conditions.

# Returns
- `Matrix{ET}`: The reduced density matrix for the specified subsystem.

# Details
The subsystem indices are 1-based and refer to the site positions in the chain. The function correctly handles both contiguous and non-contiguous subsystems.

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis, LinearAlgebra, Random

julia> Random.seed!(42); N = 4; model = AnyonModel(FibonacciAnyon(), N; pbc=true);

julia> basis = anyon_basis(model); state = randn(ComplexF64, length(basis)); state ./= norm(state);

julia> rdm = anyon_rdm(model, [1, 2], state); size(rdm)
(3, 3)

julia> ishermitian(rdm)
true

julia> isapprox(tr(rdm), 1.0; atol=1e-10)
true
```
"""
function anyon_rdm(AT::AbstractAnyonType, ::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true;) where {N,T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.

    unsorted_basis = anyon_basis(AT, T, pbc=pbc) # Input like: FibonacciAnyon(), IsingAnyon(), not the type, the instance. 
    subsystems=connected_components(subsystems)
    lengthlis=length.(subsystems)
    subsystems=vcat(subsystems...)
    # mask = bmask(T, subsystems...)
    mask = bmask(T, (N .-subsystems .+1)...)

    
    order = sortperm(unsorted_basis, by = x -> (takeenviron(x, mask), takesystem(x, mask))) #first sort by environment, then by system. The order of environment doesn't matter.

    if isempty(subsystems)
        return ones(ET, 1, 1) # Return empty matrix if no subsystems
    elseif subsystems == collect(1:N)
        # If subsystems is all sites, return the full density matrix
        return state * state' # Full density matrix
    end
    @assert length(unsorted_basis) == length(state) "state length is expected to be $(length(unsorted_basis)), but got $(length(state))"

    basis, state = unsorted_basis[order], state[order]

    
    reduced_basis = move_subsystem.(T, joint_basis(AT, lengthlis), Ref(subsystems))
    # The reduced_basis counting order doesn't matter, as long as it place subsystem part in basis correctly.
    # TIPS: mask must have common non-zero part with reduced_basis (thus for comparing)
    len = length(reduced_basis)
    # Initialize the reduced density matrix
    reduced_dm = zeros(ET, (len, len))

    # Keep track of indices where the key changes
    result_indices = Int[]
    current_key = -1
    for (idx, i) in enumerate(basis)
        key = takeenviron(i, mask)  # Get environment l bits
        if key != current_key
            @assert key > current_key "key is expected to be greater than $current_key, but got $key"
            push!(result_indices, idx)
            current_key = key
        end
    end

    # Add the final index to get complete ranges
    push!(result_indices, length(basis) + 1)

    for i in 1:length(result_indices)-1
        range = result_indices[i]:result_indices[i+1]-1         
        # Get indices in the reduced basis
        indices = searchsortedfirst.(Ref(reduced_basis), takesystem.(basis[range], mask))
        view(reduced_dm, indices, indices) .+= view(state, range) .* view(state, range)'
    end

    return reduced_dm
end
anyon_rdm(model::AnyonModel, subsystems::Vector{Int64}, state::Vector{ET}) where {ET} = anyon_rdm(model.anyon_type, BitStr{model.N, Int}, subsystems, state, model.pbc;)

function anyon_rdm(AT::AbstractAnyonType, ::Type{T}, subsystems::Vector{Int64}, state::Matrix{ET}, pbc::Bool=true;) where {N,T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.

    unsorted_basis = anyon_basis(AT, T, pbc=pbc)
    subsystems=connected_components(subsystems)
    lengthlis=length.(subsystems)
    subsystems=vcat(subsystems...)
    # mask = bmask(T, subsystems...)
    mask = bmask(T, (N .-subsystems .+1)...)

    
    order = sortperm(unsorted_basis, by = x -> (takeenviron(x, mask), takesystem(x, mask))) #first sort by environment, then by system. The order of environment doesn't matter.

    if isempty(subsystems)
        return ones(ET, 1, 1) # Return empty matrix if no subsystems
    elseif subsystems == collect(1:N)
        return state
    end
    @assert length(unsorted_basis) == size(state, 2) "state size is expected to be $(length(unsorted_basis)), but got $(size(state, 2))"

    basis, state = unsorted_basis[order], state[order, order]

    
    reduced_basis = move_subsystem.(T, joint_basis(AT, lengthlis), Ref(subsystems))
    # The reduced_basis counting order doesn't matter, as long as it place subsystem part in basis correctly.
    # TIPS: mask must have common non-zero part with reduced_basis (thus for comparing)
    len = length(reduced_basis)
    # Initialize the reduced density matrix
    reduced_dm = zeros(ET, (len, len))

    # Keep track of indices where the key changes
    result_indices = Int[]
    current_key = -1
    for (idx, i) in enumerate(basis)
        key = takeenviron(i, mask)  # Get environment l bits
        if key != current_key
            @assert key > current_key "key is expected to be greater than $current_key, but got $key"
            push!(result_indices, idx)
            current_key = key
        end
    end

    # Add the final index to get complete ranges
    push!(result_indices, length(basis) + 1)

    for i in 1:length(result_indices)-1
        range = result_indices[i]:result_indices[i+1]-1         
        # Get indices in the reduced basis
        indices = searchsortedfirst.(Ref(reduced_basis), takesystem.(basis[range], mask))
        view(reduced_dm, indices, indices) .+= view(state, range, range)
    end

    return reduced_dm
end
anyon_rdm(model::AnyonModel, subsystems::Vector{Int64}, state::Matrix{ET}) where {ET} = anyon_rdm(model.anyon_type, BitStr{model.N, Int}, subsystems, state, model.pbc;)

"""
    mapst_sec2tot(AT::AbstractAnyonType, ::Type{T}, state::Vector{ET}, k::Int; pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    mapst_sec2tot(model::AnyonModel, state::Vector{ET}, k::Int)

Map state in momentum sector to total Hilbert space.

# Arguments
## High-level interface
- `model::AnyonModel`: An `AnyonModel` object containing system parameters
- `state::Vector{ET}`: The state vector in the Hilbert space of the specified momentum sector
- `k::Int`: The momentum sector index (0 ≤ k ≤ N-1)

## Low-level interface
- `AT::AbstractAnyonType`: The anyon type.
- `T::Type`: The `BitStr` type specifying the chain length `N`.
- `state::Vector{ET}`: The state vector in the momentum sector.
- `k::Int`: The momentum sector index.
- `pbc::Bool=true`: Specifies periodic boundary conditions (required for momentum sectors).

# Returns
- `Vector{ET}`: The state vector represented in the total Hilbert space basis.

# Examples
```jldoctest
julia> using FibonacciChain, LinearAlgebra

julia> model = AnyonModel(FibonacciAnyon(), 6; pbc=true); H_k0 = anyon_ham(model, 0);

julia> eigvecs_k0 = eigvecs(H_k0); kstate = eigvecs_k0[:, 1];

julia> total_state = mapst_sec2tot(model, kstate, 0); length(total_state) == length(anyon_basis(model))
true
```
"""
function mapst_sec2tot(AT::AbstractAnyonType, ::Type{T}, state::Vector{ET}, k::Int; pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    # Map the symmetric sector Hilbert space state to total space state
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basis = anyon_basis(AT, T; pbc=pbc)
    k_dic = Dict{Int, Vector{Int64}}()
    basisK = Vector{T}(undef, 0)
    for i in eachindex(basis)
        base=basis[i]
        category = get_representative(base)[1]
        if haskey(k_dic, category)
            push!(k_dic[category], i)
        else
            k_dic[category] = [i]
        end
    end
    
    for j in eachindex(basis)
        n=basis[j]
        RS = get_representative(n)[1]
        if RS == n && (k * length(k_dic[RS])) % N == 0
            push!(basisK, n)
        end
    end
    
    total_state = zeros(ET, length(basis))
    for (i, basis) in enumerate(basisK)
        state_indices = k_dic[basis]  
        l = length(state_indices)
        total_state[state_indices] .+= 1/sqrt(l) * state[i]
    end

    return total_state
end
# High-level wrapper for mapst_sec2tot
mapst_sec2tot(model::AnyonModel, state::Vector{ET}, k::Int) where {ET} = mapst_sec2tot(model.anyon_type, BitStr{model.N, Int}, state, k; pbc=model.pbc)

"""
    anyon_rdm_sec(model::AnyonModel, subsystems::Vector{Int64}, kstate::Vector{ET}, k::Int) where {ET}

Compute reduced density matrix for specified subsystems from quantum state in specific momentum sector.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `subsystems::Vector{Int64}`: Indices of subsystem sites to keep
- `kstate::Vector{ET}`: State vector in momentum sector Hilbert space
- `k::Int`: Momentum sector (0 ≤ k ≤ N-1)

# Returns
- `Matrix{ET}`: Reduced density matrix in total Hilbert basis

# Examples
```jldoctest
julia> using FibonacciChain, LinearAlgebra

julia> model = AnyonModel(FibonacciAnyon(), 6; pbc=true);

julia> H_k0 = anyon_ham(model, 0);

julia> eigvecs_k0 = eigvecs(H_k0);

julia> kstate = eigvecs_k0[:, 1];

julia> rdm = anyon_rdm_sec(model, [1, 2, 3], kstate, 0);

julia> isapprox(tr(rdm), 1.0; atol=1e-10)
true
```
"""
function anyon_rdm_sec(model::AnyonModel{AT}, subsystems::Vector{Int64}, kstate::Vector{ET}, k::Int;) where {ET, AT <: AbstractAnyonType}
    l = length(anyon_basis(model, k)[1])
    @assert length(kstate) == l "state length is expected to be $(l), but got $(length(kstate))"
    state = mapst_sec2tot(model, kstate, k)
    reduced_dm = anyon_rdm(model, subsystems, state)
    return reduced_dm
end

# create Fibonacci basis composed of multiple disjoint chains with different basis type
function joint_basis(AT1::AbstractAnyonType, AT2::AbstractAnyonType, lengthlisA::Vector{Int}, lengthlisB::Vector{Int};subApbc::Bool=false, subBpbc::Bool=false)
    # subpbc is used to indicate whether the subsystem is periodic or not
    lisA = sort(mapreduce(len -> anyon_basis(AT1, BitStr{len, Int}, pbc= subApbc), process_join, lengthlisA))
    lisB = sort(mapreduce(len -> anyon_basis(AT2, BitStr{len, Int}, pbc= subBpbc), process_join, lengthlisB))
    return vec([join(i,j) for i in lisA for j in lisB])
end

"""
    disjoint_rdm(AT1::AbstractAnyonType, AT2::AbstractAnyonType, ::Type{T1}, ::Type{T2}, subsystemsA::Vector{Int64}, subsystemsB::Vector{Int64}, state::Vector{ET}; totalsubApbc::Bool=false, totalsubBpbc::Bool=false) where {N1, N2, T1 <: BitStr{N1}, T2 <: BitStr{N2}, ET}
    disjoint_rdm(model1::AnyonModel, model2::AnyonModel, subsystemsA::Vector{Int64}, subsystemsB::Vector{Int64}, state::Vector{ET})

Compute reduced density matrix for two joint different systems, or two parallel chains.

# Arguments
## Low-level interface
- `AT1::AbstractAnyonType`: Anyon type for subsystem A
- `AT2::AbstractAnyonType`: Anyon type for subsystem B
- `T1::Type`: BitStr type specifying chain length N1
- `T2::Type`: BitStr type specifying chain length N2
- `subsystemsA::Vector{Int64}`: Indices of subsystem A sites to keep
- `subsystemsB::Vector{Int64}`: Indices of subsystem B sites to keep
- `state::Vector{ET}`: Quantum state vector in total Hilbert space
- `totalsubApbc::Bool=false`: Whether subsystem A is periodic
- `totalsubBpbc::Bool=false`: Whether subsystem B is periodic

## High-level interface
- `model1::AnyonModel`: The model for the first subsystem (A).
- `model2::AnyonModel`: The model for the second subsystem (B).
- `subsystemsA::Vector{Int64}`: Indices of subsystem A sites to keep
- `subsystemsB::Vector{Int64}`: Indices of subsystem B sites to keep
- `state::Vector{ET}`: Quantum state vector in total Hilbert space

# Returns
- `Matrix{ET}`: Reduced density matrix

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis, LinearAlgebra

julia> N1, N2 = 4, 4;

julia> model1 = AnyonModel(FibonacciAnyon(), N1; pbc=true);

julia> model2 = AnyonModel(FibonacciAnyon(), N2; pbc=true);

julia> basisA = anyon_basis(model1);

julia> basisB = anyon_basis(model2);

julia> state = randn(ComplexF64, length(basisA) * length(basisB));

julia> state ./= norm(state);

julia> rdm = disjoint_rdm(model1, model2, [1,2], [1,2], state);

julia> size(rdm)          # reduced density matrix dimension
(9, 9)

julia> ishermitian(rdm)   # should be Hermitian
true

julia> isapprox(tr(rdm), 1.0; atol=1e-10)  # trace should be 1
true
```
"""
function disjoint_rdm(model1::AnyonModel, model2::AnyonModel, ::Type{newT}, subsystemsA::Vector{Int64}, subsystemsB::Vector{Int64}, state::Vector{ET}; totalsubApbc::Bool=false, totalsubBpbc::Bool=false) where {newT, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the two disjoint basis (two part in one chain), espeically, can be viewed as two parrelel chain. Given the system size, subsystem and basis type, to get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    # model.pbc is the basis boundary condition, while totalsubpbc is used to indicate whether the total subsystem is periodic or not. Roughly speaking, if partition is Int[] and collect(1:N2), thus the totalsubBpbc must be the same with model2.pbc.
    N1 = model1.N
    N2 = model2.N
    # newT should is BitStr{N1+N2, Int} 
    @assert all(subsystemsA .<= N1) "subsystemsA is expected to be in [1, $N1], but got $(subsystemsA)"
    @assert all(subsystemsB .<= N2) "subsystemsB is expected to be in [1, $N2], but got $(subsystemsB)"
    
    if isempty(subsystemsA) && isempty(subsystemsB)
        return ones(ET, 1, 1) # Return empty matrix if no subsystems
    elseif subsystemsA == collect(1:N1) && subsystemsB == collect(1:N2)
        return state * state'
    end

    unsorted_basisA = anyon_basis(model1)
    unsorted_basisB = anyon_basis(model2)
    lenubasisA = length(unsorted_basisA)
    lenubasisB = length(unsorted_basisB)
    doublebasis = vec([join(i,j) for i in unsorted_basisA for j in unsorted_basisB])
    # Align as [A, B] basis
    @assert lenubasisA*lenubasisB == length(state) "state length is expected to be $(lenubasisA*lenubasisB), but got $(length(state))"

    # Compute mask using original subsystems (before connected_components)
    # In the combined [A, B] basis, A occupies bits N1+N2 down to N2+1, B occupies bits N2 down to 1
    mask = bmask(newT, (N1 .- subsystemsA .+ 1) .+ N2..., (N2 .- subsystemsB .+ 1)...)
    
    # Process subsystemsA for reduced basis construction
    subsystemsA_cc = connected_components(subsystemsA)
    lengthlisA = length.(subsystemsA_cc)
    subsystemsA_flat = vcat(subsystemsA_cc...)
    
    # Process subsystemsB, shifted to combined basis indices
    subsystemsB_shifted = subsystemsB .+ N1
    subsystemsB_cc = connected_components(subsystemsB_shifted)
    lengthlisB = length.(subsystemsB_cc)
    subsystemsB_flat = vcat(subsystemsB_cc...)
    
    # Combined subsystems for move_subsystem
    combined_subsystems = vcat(subsystemsA_flat, subsystemsB_flat)

    # Construct reduced_basis
    if isempty(subsystemsA)
        # Only B subsystem, use lengthlis from B
        lengthlisB_orig = length.(connected_components(subsystemsB))
        reduced_basis = sort(move_subsystem.(newT, joint_basis(model2.anyon_type, lengthlisB_orig, totalsubBpbc), Ref(subsystemsB_flat)))
    elseif isempty(subsystemsB)
        # Only A subsystem
        reduced_basis = sort(move_subsystem.(newT, joint_basis(model1.anyon_type, lengthlisA, totalsubApbc), Ref(subsystemsA_flat)))
    else
        # Both subsystems present
        reduced_basis = sort(move_subsystem.(newT, joint_basis(model1.anyon_type, model2.anyon_type, lengthlisA, lengthlisB, subApbc=totalsubApbc, subBpbc=totalsubBpbc), Ref(combined_subsystems)))
    end

    order = sortperm(doublebasis, by = x -> (takeenviron(x, mask), takesystem(x, mask))) #first sort by environment, then by system. The order of environment doesn't matter.
    basis, state = doublebasis[order], state[order]
    
    len = length(reduced_basis)
    
    # Initialize the reduced density matrix
    reduced_dm = zeros(ET, (len, len))

    # Keep track of indices where the key changes
    result_indices = Int[]
    current_key = -1
    for (idx, i) in enumerate(basis)
        key = takeenviron(i, mask)  # Get environment l bits
        if key != current_key
            @assert key > current_key "key is expected to be greater than $current_key, but got $key"
            push!(result_indices, idx)
            current_key = key
        end
    end
    # Add the final index to get complete ranges
    push!(result_indices, length(basis) + 1)

    for i in 1:length(result_indices)-1
        range = result_indices[i]:result_indices[i+1]-1         
        # Get indices in the reduced basis
        indices = searchsortedfirst.(Ref(reduced_basis), takesystem.(basis[range], mask))
        view(reduced_dm, indices, indices) .+= view(state, range) .* view(state, range)'
    end

    return reduced_dm
end
disjoint_rdm(model1::AnyonModel, model2::AnyonModel, subsystemsA::Vector{Int64}, subsystemsB::Vector{Int64}, state::Vector{ET}; pbcA::Bool=false, pbcB::Bool=false) where {ET} =
disjoint_rdm(model1, model2, BitStr{model1.N + model2.N, Int}, subsystemsA, subsystemsB, state, totalsubApbc=pbcA, totalsubBpbc=pbcB
)