function Fibonacci_chain_OBC(::Type{T}) where {N,T<:BitStr{N}}
    # Generate Fibonacci chain for Fibonacci model with open boundary condition
    fib_chain=[[T(0), T(1)], [T(0), T(1), T(2)]]
    for i = 3:N
        push!(
            fib_chain,
            vcat(
                [s << 1 for s in fib_chain[i-1]],
                [(s << 2 | T(1)) for s in fib_chain[i-2]],
            ),
        )
    end
    # each push we add a bit"0" or bit"01" to the end of the bit_string, and finally return the last element of the fib_chain
    return fib_chain[N]
end

function Fibonacci_chain_PBC(::Type{T}) where {N,T<:BitStr{N}}
    # Generate Fibonacci chain  for Fibonacci model with periodic boundary condition
    return filter(c -> iszero((c >> (N-1)) & (c & 1)), Fibonacci_chain_OBC(T))
end

abstract type AbstractAnyonType end

# Each model type carries its own (possibly constrained) Hilbert space, defined
# once and for all by its `anyon_basis` method — this is the single extension
# point when adding a new chain:
struct FibonacciAnyon <: AbstractAnyonType end  # Fibonacci fusion space: bitstrings with no adjacent τ (1)s, dim ~ φ^N
struct IsingAnyon <: AbstractAnyonType end      # transverse-field Ising *spin* chain: full 2^N product basis
struct OBFAnyon <: AbstractAnyonType end        # O'Brien-Fendley spin-1/2 chain: full 2^N product basis
struct HeisenbergAnyon <: AbstractAnyonType end # spin-1/2 Heisenberg (XXZ) chain: full 2^N product basis

"""
    AnyonModel(anyon_type::AbstractAnyonType, N::Int; pbc::Bool=true, measure_operator::Symbol=:Antiferro, kwargs...)

Represents a 1D anyon/spin chain model.

# Hilbert space
Each anyon type defines its own (possibly constrained) Hilbert space via `anyon_basis`:
- `FibonacciAnyon()`: Fibonacci fusion space — bitstrings with no adjacent τ (`1`)s, dimension ~ φ^N.
- `IsingAnyon()`: transverse-field Ising *spin* chain — full 2^N product basis.
- `OBFAnyon()`: O'Brien-Fendley spin-1/2 chain — full 2^N product basis.
- `HeisenbergAnyon()`: spin-1/2 Heisenberg (XXZ) chain — full 2^N product basis.

# Fields
- `anyon_type::AbstractAnyonType`: The type of anyon, e.g., `FibonacciAnyon()` or `IsingAnyon()`.
- `N::Int`: The number of sites in the chain.
- `pbc::Bool`: A boolean indicating whether periodic boundary conditions are applied (`true`) or not (`false`).
- `measure_operator::Symbol`: The operator type used for defining the Hamiltonian, e.g., `:Antiferro` or `:Ferro` for Fibonacci anyons.
- `params::Dict{Symbol, Any}`: Additional model parameters (e.g., `J`, `h` for Ising model, `J`, `Jz`, `h` for Heisenberg chain).

# Examples
```julia
# Fibonacci anyon model
model_fibo = AnyonModel(FibonacciAnyon(), 6; pbc=true, measure_operator=:Antiferro)

# Ising model with custom couplings
model_ising = AnyonModel(IsingAnyon(), 10; pbc=true, measure_operator=:X, J=1.0, h=0.5)

# Heisenberg (XXZ) chain
model_heis = AnyonModel(HeisenbergAnyon(), 10; pbc=true, J=1.0, Jz=0.5, h=0.1)

# Access parameters
model_ising.params[:J]  # returns 1.0
```
"""
struct AnyonModel{AT<:AbstractAnyonType}
    anyon_type::AT
    N::Int
    pbc::Bool
    measure_operator::Symbol
    params::Dict{Symbol,Float64}
    function AnyonModel(
        anyon_type::FibonacciAnyon,
        N::Int;
        pbc::Bool = true,
        measure_operator::Symbol = :Antiferro,
        kwargs...,
    )
        @assert N > 0 "N is expected to be greater than 0, but got $N"
        @assert measure_operator ∈ [:Ferro, :Antiferro, :reset] "measure_operator must be :Ferro, :Antiferro, :reset for Fibonacci anyons"
        params = Dict{Symbol,Float64}(k => Float64(v) for (k, v) in kwargs)
        return new{FibonacciAnyon}(anyon_type, N, pbc, measure_operator, params)
    end
    function AnyonModel(
        anyon_type::IsingAnyon,
        N::Int;
        pbc::Bool = true,
        measure_operator::Symbol = :X,
        kwargs...,
    )
        @assert N > 0 "N is expected to be greater than 0, but got $N"
        @assert measure_operator in [:X, :ZZ, :Z, :reset] "measure_operator must be either :X, :ZZ, :Z, :reset for Ising anyons"
        params = Dict{Symbol,Float64}(k => Float64(v) for (k, v) in kwargs)
        return new{IsingAnyon}(anyon_type, N, pbc, measure_operator, params)
    end
    function AnyonModel(
        anyon_type::OBFAnyon,
        N::Int;
        pbc::Bool = true,
        measure_operator::Symbol = :X,
        kwargs...,
    )
        @assert N > 0 "N is expected to be greater than 0, but got $N"
        @assert measure_operator in [:XZZ, :ZZX, :ZZ, :X] "measure_operator must be :XZZ, :ZZX, :ZZ, :X for OBF anyons"
        params = Dict{Symbol,Float64}(k => Float64(v) for (k, v) in kwargs)
        return new{OBFAnyon}(anyon_type, N, pbc, measure_operator, params)
    end
    function AnyonModel(
        anyon_type::HeisenbergAnyon,
        N::Int;
        pbc::Bool = true,
        measure_operator::Symbol = :X,
        kwargs...,
    )
        @assert N > 0 "N is expected to be greater than 0, but got $N"
        @assert measure_operator in [:X, :ZZ, :Z, :reset] "measure_operator must be either :X, :ZZ, :Z, :reset for Heisenberg chain"
        params = Dict{Symbol,Float64}(k => Float64(v) for (k, v) in kwargs)
        return new{HeisenbergAnyon}(anyon_type, N, pbc, measure_operator, params)
    end
end

# Helper function to get parameter with default value
get_interaction_param(model::AnyonModel, key::Symbol, default) =
    get(model.params, key, default)

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
function anyon_basis(
    ::FibonacciAnyon,
    ::Type{T};
    pbc::Bool = true,
    symmetry_block = nothing,
) where {N,T<:BitStr{N}}
    # Generate basis for Fibonacci model, return BitBasis form, which can be used as binary and decimal form. Here we both consider PBC and OBC
    @assert symmetry_block === nothing || symmetry_block in [0, 1, :tau, :trivial] "symmetry_block is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $symmetry_block"
    @assert T <: BitStr{N} "Type T must be a BitStr type"
    # If pbc is true, use Fibonacci_chain_PBC, otherwise use Fibonacci_chain_OBC
    if pbc
        if T == BitStr{1,Int}
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

function anyon_basis(
    ::IsingAnyon,
    ::Type{T};
    pbc::Bool = true,
    symmetry_block = nothing,
) where {N,T<:BitStr{N}}
    return spin_half_basis(T)
end

function anyon_basis(
    ::OBFAnyon,
    ::Type{T};
    pbc::Bool = true,
    symmetry_block = nothing,
) where {N,T<:BitStr{N}}
    return spin_half_basis(T)
end

function anyon_basis(
    ::HeisenbergAnyon,
    ::Type{T};
    pbc::Bool = true,
    symmetry_block = nothing,
) where {N,T<:BitStr{N}}
    return spin_half_basis(T)
end

# Spin-1/2 chains (Ising, OBF, Heisenberg) share the full unconstrained product basis.
spin_half_basis(::Type{T}) where {N,T<:BitStr{N}} = T[T(i) for i = 0:(2^N-1)]

anyon_basis(model::AnyonModel; symmetry_block = nothing) = anyon_basis(
    model.anyon_type,
    BitStr{model.N,Int};
    pbc = model.pbc,
    symmetry_block = symmetry_block,
)


# Fibonacci F-symbol (F^{τστ}_δ)^μ_ν in the chain's bit convention:
# bit 1 = vacuum label, bit 0 = τ (fusion paths have no adjacent vacuums,
# i.e. no adjacent 1 bits). In the standard gauge the only nontrivial block is
# F^{τττ}_τ = [[ϕ⁻¹, ϕ^(-1/2)], [ϕ^(-1/2), -ϕ⁻¹]] (rows/cols ordered (vacuum, τ));
# all fusion-consistent symbols involving the vacuum equal 1, inconsistent ones vanish.
function _Fibonacci_Fsymbol(bσ::Int, bδ::Int, bμ::Int, bν::Int)
    ϕ = (1+√5)/2
    if bσ == 1 || bδ == 1
        # σ or δ is the vacuum: τ×vac = τ forces μ = ν = τ
        return (bμ == 0 && bν == 0) ? 1.0 : 0.0
    else
        # σ = δ = τ: F^{τττ}_τ block
        if bμ == 1 && bν == 1
            return ϕ^(-1)
        elseif bμ == bν
            return -ϕ^(-1)
        else
            return ϕ^(-1/2)
        end
    end
end

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
function Fsymmetry_coef(
    ::FibonacciAnyon,
    state::T,
    base::T;
    pbc::Bool = true,
) where {N,T<:BitStr{N}}
    # ⟨x₁x₂⋯x_N|Y|x'₁x'₂⋯x'_N⟩ = ∏ᵢ (F^{τ xᵢ τ}_{x'ᵢ₊₁})^{xᵢ₊₁}_{x'ᵢ},
    # where x (x') are the bits of `state` (`base`), read left to right,
    # i.e. xᵢ = state[N-i+1], and x_{N+1} ≡ x₁ for PBC.
    coef = 1.0
    nfactor = pbc ? N : N - 1
    for i = 1:nfactor
        j = mod1(i + 1, N)
        coef *= _Fibonacci_Fsymbol(
            state[N-i+1],
            base[N-j+1],
            state[N-j+1],
            base[N-i+1],
        )
        iszero(coef) && return 0.0
    end

    return coef
end
Fsymmetry_coef(
    model::AnyonModel{FibonacciAnyon},
    state::T,
    base::T,
) where {N,T<:BitStr{N}} = Fsymmetry_coef(model.anyon_type, state, base; pbc = model.pbc)

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
function topological_symmetry_basismap(
    model::AnyonModel{AT},
    state::T,
) where {N,T<:BitStr{N},AT<:FibonacciAnyon}
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
    topological_charge_operator(model::AnyonModel{AT},) where {AT<:FibonacciAnyon}

Compute the topological charge operator matrix for Fibonacci anyon model.

# Arguments
- `model::AnyonModel{FibonacciAnyon}`: Fibonacci anyon model

# Returns
- `Matrix{Float64}`: Topological charge operator matrix (Y matrix)

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; 

julia> model = AnyonModel(FibonacciAnyon(), N; pbc=true);

julia> Y = topological_charge_operator(model);

julia> size(Y) == (length(anyon_basis(model)), length(anyon_basis(model)))
true
```
"""
function topological_charge_operator(
    model::AnyonModel{AT},
) where {AT<:FibonacciAnyon}
    # compute the topological charge operator Yl in the Fibonacci model. default l=0, for tau. l=1, for vacuum.
    # Y[i, j] = ⟨basisᵢ|Y|basisⱼ⟩ = ∏ₖ (F^{τ xₖ τ}_{x'ₖ₊₁})^{xₖ₊₁}_{x'ₖ} with x = basisᵢ, x' = basisⱼ.
    basis=anyon_basis(model)
    Ymatrix=[Fsymmetry_coef(model, bra, ket) for bra in basis, ket in basis]

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
function anyon_basis(
    AT::AbstractAnyonType,
    ::Type{T},
    k::Int;
    symmetry_block = nothing,
) where {N,T<:BitStr{N}}
    #params: a int of lattice number, momentum of system, topological_charge symmetry_block, which default to be nothing
    #return: computational basis in given momentum kinetically constrained subspace with decimal int form in golden chain model
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert symmetry_block === nothing || symmetry_block in [0, 1, :tau, :trivial] "symmetry_block is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $symmetry_block"

    basisK = Vector{T}(undef, 0)
    basis = anyon_basis(AT, T; symmetry_block = symmetry_block)
    basis_dic = Dict{T,Vector{T}}()

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
anyon_basis(model::AnyonModel, k::Int; symmetry_block = nothing) =
    anyon_basis(model.anyon_type, BitStr{model.N,Int}, k; symmetry_block = symmetry_block)

function cyclebits(state::T) where {N,T<:BitStr{N}}
    #params: t is an integer, N is the length of the binary string
    #We also use this order: system size, state, circular shift bitstring 1 bit.
    # In case need to shift more than 1 bit, we can use a loop or recursion. or we leave a interface here  n_translations::Int
    mask = 1 << N - 1
    return ((state << 1) | (state >> (N - 1))) & mask
end

function get_representative(state::T) where {N,T<:BitStr{N}}
    #Finds representative and representative translation for a state.
    #State should be a decimal integer.

    representative = state
    translation = 0
    # cycle(bits) = (bits << 1) % (2^N - 1)  # Left shift the state by one position
    current = state
    for n_translation_sites = 1:(N-1)
        current = cyclebits(current)  # Cycle the bits
        if current < representative
            representative = current
            translation = n_translation_sites
        end
    end

    return representative, translation
end

# join two lists of basis by make a product of two lists, b is placed after a (counts from left to right)
function process_join(a, b)
    # effectively, vec([join(a, b) for a in a for b in b]) # seems faster
    return [join(a, b) for a in a for b in b]
end

# create Fibonacci basis composed of multiple disjoint sub-chains
function joint_basis(AT::AbstractAnyonType, lengthlis::Vector{Int}, pbc::Bool = false;)
    # The element order in lengthlis doesn't matter, i.e., [2, 3] is only different with [3, 2] in basis order. Only your input state must consistent with the basis order.
    if isempty(lengthlis)
        return BitStr{0,Int}[]
    else
        return sort(
            mapreduce(
                len -> anyon_basis(AT, BitStr{len,Int}; pbc = pbc),
                process_join,
                lengthlis,
            ),
        )
    end
end

function connected_components(v::Vector{Int})
    if isempty(v)
        return Int[]
    end

    sort!(v)

    result = []
    current_segment = [v[1]]

    for i = 2:length(v)
        if v[i] == v[i-1] + 1
            push!(current_segment, v[i])
        else
            push!(result, current_segment)
            current_segment = [v[i]]
        end
    end

    push!(result, current_segment)

    return result
end

function move_subsystem(
    ::Type{BitStr{M,INT}},
    basis::BitStr{N,INT},
    subsystems::Vector{Int},
) where {M,N,INT}
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
function anyon_rdm(
    AT::AbstractAnyonType,
    ::Type{T},
    subsystems::Vector{Int64},
    state::Vector{ET},
    pbc::Bool = true;
) where {N,T<:BitStr{N},ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.

    unsorted_basis = anyon_basis(AT, T, pbc = pbc) # Input like: FibonacciAnyon(), IsingAnyon(), not the type, the instance. 
    subsystems=connected_components(subsystems)
    lengthlis=length.(subsystems)
    subsystems=vcat(subsystems...)
    # mask = bmask(T, subsystems...)
    mask = bmask(T, (N .- subsystems .+ 1)...)


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

    for i = 1:(length(result_indices)-1)
        range = result_indices[i]:(result_indices[i+1]-1)
        # Get indices in the reduced basis
        indices = searchsortedfirst.(Ref(reduced_basis), takesystem.(basis[range], mask))
        view(reduced_dm, indices, indices) .+= view(state, range) .* view(state, range)'
    end

    return reduced_dm
end
anyon_rdm(model::AnyonModel, subsystems::Vector{Int64}, state::Vector{ET}) where {ET} =
    anyon_rdm(model.anyon_type, BitStr{model.N,Int}, subsystems, state, model.pbc;)

function anyon_rdm(
    AT::AbstractAnyonType,
    ::Type{T},
    subsystems::Vector{Int64},
    state::Matrix{ET},
    pbc::Bool = true;
) where {N,T<:BitStr{N},ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.

    unsorted_basis = anyon_basis(AT, T, pbc = pbc)
    subsystems=connected_components(subsystems)
    lengthlis=length.(subsystems)
    subsystems=vcat(subsystems...)
    # mask = bmask(T, subsystems...)
    mask = bmask(T, (N .- subsystems .+ 1)...)


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

    for i = 1:(length(result_indices)-1)
        range = result_indices[i]:(result_indices[i+1]-1)
        # Get indices in the reduced basis
        indices = searchsortedfirst.(Ref(reduced_basis), takesystem.(basis[range], mask))
        view(reduced_dm, indices, indices) .+= view(state, range, range)
    end

    return reduced_dm
end
anyon_rdm(model::AnyonModel, subsystems::Vector{Int64}, state::Matrix{ET}) where {ET} =
    anyon_rdm(model.anyon_type, BitStr{model.N,Int}, subsystems, state, model.pbc;)

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
function mapst_sec2tot(
    AT::AbstractAnyonType,
    ::Type{T},
    state::Vector{ET},
    k::Int;
    pbc::Bool = true,
) where {N,T<:BitStr{N},ET}
    # Map the symmetric sector Hilbert space state to total space state
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basis = anyon_basis(AT, T; pbc = pbc)
    k_dic = Dict{Int,Vector{Int64}}()
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
mapst_sec2tot(model::AnyonModel, state::Vector{ET}, k::Int) where {ET} =
    mapst_sec2tot(model.anyon_type, BitStr{model.N,Int}, state, k; pbc = model.pbc)

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
function anyon_rdm_sec(
    model::AnyonModel{AT},
    subsystems::Vector{Int64},
    kstate::Vector{ET},
    k::Int;
) where {ET,AT<:AbstractAnyonType}
    l = length(anyon_basis(model, k)[1])
    @assert length(kstate) == l "state length is expected to be $(l), but got $(length(kstate))"
    state = mapst_sec2tot(model, kstate, k)
    reduced_dm = anyon_rdm(model, subsystems, state)
    return reduced_dm
end

# create Fibonacci basis composed of multiple disjoint chains with different basis type
function joint_basis(
    AT1::AbstractAnyonType,
    AT2::AbstractAnyonType,
    lengthlisA::Vector{Int},
    lengthlisB::Vector{Int};
    subApbc::Bool = false,
    subBpbc::Bool = false,
)
    # subpbc is used to indicate whether the subsystem is periodic or not
    lisA = sort(
        mapreduce(
            len -> anyon_basis(AT1, BitStr{len,Int}, pbc = subApbc),
            process_join,
            lengthlisA,
        ),
    )
    lisB = sort(
        mapreduce(
            len -> anyon_basis(AT2, BitStr{len,Int}, pbc = subBpbc),
            process_join,
            lengthlisB,
        ),
    )
    return vec([join(i, j) for i in lisA for j in lisB])
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
function disjoint_rdm(
    model1::AnyonModel,
    model2::AnyonModel,
    ::Type{newT},
    subsystemsA::Vector{Int64},
    subsystemsB::Vector{Int64},
    state::Vector{ET};
    totalsubApbc::Bool = false,
    totalsubBpbc::Bool = false,
) where {newT,ET}
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
    doublebasis = vec([join(i, j) for i in unsorted_basisA for j in unsorted_basisB])
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
        reduced_basis = sort(
            move_subsystem.(
                newT,
                joint_basis(model2.anyon_type, lengthlisB_orig, totalsubBpbc),
                Ref(subsystemsB_flat),
            ),
        )
    elseif isempty(subsystemsB)
        # Only A subsystem
        reduced_basis = sort(
            move_subsystem.(
                newT,
                joint_basis(model1.anyon_type, lengthlisA, totalsubApbc),
                Ref(subsystemsA_flat),
            ),
        )
    else
        # Both subsystems present
        reduced_basis = sort(
            move_subsystem.(
                newT,
                joint_basis(
                    model1.anyon_type,
                    model2.anyon_type,
                    lengthlisA,
                    lengthlisB,
                    subApbc = totalsubApbc,
                    subBpbc = totalsubBpbc,
                ),
                Ref(combined_subsystems),
            ),
        )
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

    for i = 1:(length(result_indices)-1)
        range = result_indices[i]:(result_indices[i+1]-1)
        # Get indices in the reduced basis
        indices = searchsortedfirst.(Ref(reduced_basis), takesystem.(basis[range], mask))
        view(reduced_dm, indices, indices) .+= view(state, range) .* view(state, range)'
    end

    return reduced_dm
end
disjoint_rdm(
    model1::AnyonModel,
    model2::AnyonModel,
    subsystemsA::Vector{Int64},
    subsystemsB::Vector{Int64},
    state::Vector{ET};
    pbcA::Bool = false,
    pbcB::Bool = false,
) where {ET} = disjoint_rdm(
    model1,
    model2,
    BitStr{model1.N + model2.N,Int},
    subsystemsA,
    subsystemsB,
    state,
    totalsubApbc = pbcA,
    totalsubBpbc = pbcB,
)
