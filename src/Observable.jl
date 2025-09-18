"""
    ee(subrm::Matrix{ET}) where {ET}

Calculate entanglement entropy from reduced density matrix.

# Arguments
- `subrm::Matrix{ET}`: Hermitian reduced density matrix

# Returns
- `Float64`: von Neumann entanglement entropy

# Examples
```jldoctest
julia> using FibonacciChain, LinearAlgebra

julia> # Create a simple 2x2 density matrix
       ρ = [0.5 0.0; 0.0 0.5];  # Maximally mixed state

julia> entropy = ee(ρ);

julia> abs(entropy - log(2)) < 1e-10  # Should equal log(2) ≈ 0.693
true

julia> # Pure state has zero entropy
       ρ_pure = [1.0 0.0; 0.0 0.0];

julia> ee(ρ_pure) ≈ 0.0
true
```
"""
function ee(subrm::Matrix{ET}) where {ET}
    #  subrm=qi.ptrace(state*state',[2 for i in 1:N],[i for i in l+1:N])
    @assert ishermitian(subrm) "The reduced density matrix is not hermitian."
    spectrum=eigvals(subrm)
    EE=0
    for i in eachindex(spectrum)
        v=abs(spectrum[i])
            if v>1e-8
                EE+=-v*log(v)
            end
    end

    return EE
end

"""
    anyon_eelis(N::Int64, state::Union{Vector{ET}, Matrix{ET}}, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET}

Calculate entanglement entropy profile along the chain for given quantum state or density matrix.

# Arguments
- `N::Int64`: System size
- `state::Union{Vector{ET}, Matrix{ET}}`: Quantum state vector or density matrix in anyon basis
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Vector{Float64}`: Entanglement entropy at each bipartition from left to right

# Examples
```jldoctest
julia> using FibonacciChain, LinearAlgebra

julia> N = 6;

julia> # Create ground state of Fibonacci Hamiltonian
       H = anyon_ham(N, true);

julia> eigenvals, eigenvecs = eigen(H);

julia> ground_state = eigenvecs[:, 1];

julia> # Calculate entanglement entropy profile
       ee_profile = anyon_eelis(N, ground_state, true);

julia> length(ee_profile) == N - 1  # Profile has N-1 points
true

julia> all(x -> x ≥ 0, ee_profile)  # All entropies are non-negative
true
```
"""
function anyon_eelis(N::Int64,state::Union{Vector{ET}, Matrix{ET}},pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET}
    # Generate ee list for a given state from the left to the right
    splitlis=Vector(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        if m<= div(N,2)
            subrho=anyon_rdm(N, collect(1:m), state, pbc, anyon_type=anyon_type)
            EE_lis[m]=ee(subrho)
        else
            subrho=anyon_rdm(N, collect(m+1:N), state, pbc, anyon_type=anyon_type)
            EE_lis[m]=ee(subrho)
        end
    end
    return EE_lis
end

function anyonladder_eelis(N::Int64,state::Vector{ET},pbc::Bool=true) where {ET}
    # Generate ee list for a given state from the left to the right
    splitlis=Vector(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        if m<= div(N,2)
            subrho=ladderrdm(N, collect(1:m), state, pbc)
            EE_lis[m]=ee(subrho)
        else
            subrho=ladderrdm(N, collect(m+1:N), state, pbc)
            EE_lis[m]=ee(subrho)
        end
    end
    return EE_lis
end

"""
    translation_matrix(::Type{T}) where {N, T <: BitStr{N}}

Generate translation operator matrix for Fibonacci basis states.

# Arguments
- `T::Type`: BitStr type specifying chain length N

# Returns
- `Matrix{Float64}`: Translation matrix mapping each basis state to its translated version
"""
function translation_matrix(::Type{T}) where {N, T <: BitStr{N}}
    basis=anyon_basis(T) 
    l = length(basis) 
    Mat=zeros(Float64,(l,l))
    translated_basis = cyclebits.(basis) # Use broadcasting to apply cyclebits to each element in basis
    order = searchsortedfirst.(Ref(basis), translated_basis) # Find the indices of the translated basis in the original basis
    for i in 1:l
        Mat[i, order[i]] += 1.0
    end
    
    return Mat
end
translation_matrix(N::Int) = translation_matrix(BitStr{N, Int})

"""
    inversion_matrix(::Type{T}) where {N, T <: BitStr{N}}

Generate spatial inversion operator matrix for Fibonacci basis states.

# Arguments
- `T::Type`: BitStr type specifying chain length N

# Returns
- `Matrix{Float64}`: Inversion matrix mapping each basis state to its spatially reflected version
"""
function inversion_matrix(::Type{T}) where {N, T <: BitStr{N}}
    basis=anyon_basis(T)
    l=length(basis)
    Imatrix=zeros((l,l))
    # reversed_basis = map(breflect, basis) # The optimization try of using map function and broadcast
    reversed_basis=breflect.(basis)
    order = searchsortedfirst.(Ref(basis), reversed_basis) # Find the indices of the reversed basis in the original basis
    # Imatrix[CartesianIndex.(collect(1:length(basis)),searchsortedfirst.(Ref(basis), reversed_basis))].+=1.0
    for i in 1:l
        Imatrix[i,order[i]]+=1.0
    end
   
    return Imatrix
end
inversion_matrix(N::Int) = inversion_matrix(BitStr{N, Int})

"""
    braidingsqmap(::Type{T}, state::Vector{ET}, idx::Int, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}

Apply braiding squared operation to quantum state at specified site.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `state::Vector{ET}`: Quantum state vector in anyon basis
- `idx::Int`: Site index for braiding operation
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `Vector{ET}`: Transformed state after braiding operation

Braiding is a fundamental topological operation that exchanges adjacent anyons.
"""
function braidingsq_basismap(::Type{T}, state::T, i::Int, pbc::Bool=true) where {N, T <: BitStr{N}}
    # default for PBC system
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    ϕ = (1+√5)/2
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))
    
    if 2<= i <= N-1
        mask=bmask(T,1,2,3) << (N-i-1)
        str100, str101, str010, str001, str000 = T(4) << (N-i-1), T(5) << (N-i-1), T(2) << (N-i-1), T(1) << (N-i-1), T(0) << (N-i-1)
        if state & mask == str000
            return state, X(state,i), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
        elseif state & mask == str001
            return state, state, exp(-6im*π/5), 0
        elseif state & mask == str010
            return state, X(state,i), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
        elseif state & mask == str100
            return state, state, exp(-6im*π/5), 0
        elseif state & mask == str101
            return state, state, exp(-2im*π/5), 0
        end
    end
    if pbc
        if i == 1 #count from the left
        mask=bmask(T, N, N-1,1)
        str100, str101, str010, str001, str000 = bmask(T,1), bmask(T, N-1, 1), bmask(T, N), bmask(T, N-1), T(0)
            if state & mask == str000
                return state, X(state,i), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
            elseif state & mask == str001
                return state, state, exp(-6im*π/5), 0
            elseif state & mask == str010
                return state, X(state,i), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
            elseif state & mask == str100
                return state, state, exp(-6im*π/5), 0
            elseif state & mask == str101
                return state, state, exp(-2im*π/5), 0
            end
        elseif i == N #count from the left
        mask=bmask(T, N, 2, 1)
        str100, str101, str010, str001, str000 = bmask(T,2), bmask(T, N, 2), bmask(T, 1), bmask(T, N), T(0)
            if state & mask == str000
                return state, X(state,i), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
            elseif state & mask == str001
                return state, state, exp(-6im*π/5), 0
            elseif state & mask == str010
                return state, X(state,i), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
            elseif state & mask == str100
                return state, state, exp(-6im*π/5), 0
            elseif state & mask == str101
                return state, state, exp(-2im*π/5), 0
            end
        end
    end
end

function braidingsq_matrix(::Type{T}, idx::Int, pbc::Bool=true) where {N, T <: BitStr{N}}
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=anyon_basis(T, pbc)
    l=length(basis)
    Bmatrix=zeros(ComplexF64, (l,l))
    for i in 1:l
        outputstate1, outputstate2, output1, output2 = braidingsq_basismap(T, basis[i], idx, pbc)
        if output2 == 0
            Bmatrix[i,i]+=output1
        else
            j2=searchsortedfirst(basis, outputstate2)
            Bmatrix[i,i]+=output1
            Bmatrix[i,j2]+=output2
        end
    end
    
    return Bmatrix
end

"""
    braidingsqmap(::Type{T}, state::Vector{ET}, idx::Int, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET} 

Apply braiding squared operation to quantum state at specified site.
# Arguments
- `T::Type`: BitStr type specifying chain length N
- `state::Vector{ET}`: Quantum state vector in anyon basis
- `idx::Int`: Site index for braiding operation
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `Vector{ET}`: Transformed state after braiding operation
"""
function braidingsqmap(::Type{T}, state::Vector{ET}, idx::Int, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    # input a superposition state, and output the braided state
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=anyon_basis(T, pbc)
    l=length(basis)
    @assert l == length(state) "state length is expected to be $(l), but got $(length(state))"
    mapped_state = zeros(ComplexF64, length(state))
    for i in 1:l
        outputstate1, outputstate2, output1, output2 = braidingsq_basismap(T, basis[i], idx, pbc)
        if output2 == 0
            mapped_state[i]+=output1*state[i] # outputstate1 is the same as basis[i]
        else    
            j2=searchsortedfirst(basis, outputstate2)
            mapped_state[i]+=output1*state[i] # outputstate1 is the same as basis[i]
            mapped_state[j2]+=output2*state[i]
        end
    end
    
    return mapped_state
end
braidingsqmap(N::Int, state::Vector{ET}, idx::Int, pbc::Bool=true) where {ET} = braidingsqmap(BitStr{N, Int}, state, idx, pbc)

"""
    spatial_correlation(N::Int64, state::Union{Vector{ET}, Matrix{ET}}, site1::Int64, site2::Int64; pbc::Bool=true, anyon_type::Symbol=:Fibo) where {ET}

Calculate mutual information between two sites as spatial correlation measure.

# Arguments
- `N::Int64`: System size
- `state::Union{Vector{ET}, Matrix{ET}}`: Quantum state vector or density matrix
- `site1::Int64`: First site index (1 ≤ site1 ≤ N)
- `site2::Int64`: Second site index (1 ≤ site2 ≤ N, site2 ≠ site1)
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `Float64`: Mutual information I(1:2) = S(1) + S(2) - S(1,2)

Computes quantum mutual information as measure of spatial correlations.
"""
function spatial_correlation(N::Int64, state::Union{Vector{ET}, Matrix{ET}}, site1::Int64, site2::Int64; pbc::Bool=true, anyon_type::Symbol=:Fibo) where {ET}
    # Calculate the spatial correlation between two sites in a given state. For reference qubit added state, we need reference_rdm. For an initial state without reference qubit, we do not need anything.
    @assert 1 <= site1 <= N "Site1 index must be in the range [1, $(N)]"
    @assert 1 <= site2 <= N "Site2 index must be in the range [1, $(N)]"
    @assert site1 != site2 "Site1 and Site2 must be different"

    ρ1 = anyon_rdm(N, [site1], state, pbc, anyon_type=anyon_type)
    ρ2 = anyon_rdm(N, [site2], state, pbc, anyon_type=anyon_type)
    ρ12 = anyon_rdm(N, [site1, site2], state, pbc, anyon_type=anyon_type)
    
    correlation = ee(ρ1) + ee(ρ2) - ee(ρ12)

    return correlation
end

"""
    temporal_correlation(N::Int64, state_addref2::Vector{ET}; pbc::Bool=true, anyon_type::Symbol=:Fibo) where {ET}

Calculate temporal correlation using state with two reference qubits.

# Arguments
- `N::Int64`: System size
- `state_addref2::Vector{ET}`: Quantum state with two reference qubits added
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `Float64`: Temporal correlation measure between time slices

Uses reference qubit protocol to measure temporal correlations at single site.
"""
function temporal_correlation(N::Int64, state_addref2::Vector{ET}; pbc::Bool=true, anyon_type::Symbol=:Fibo) where {ET}
    # Calculate the temporal correlation between two time slices at one site in a given initial_state

    ρ1 = reference_rdm(N, [2], state_addref2, pbc=pbc, anyon_type=anyon_type)
    ρ2 = reference_rdm(N, [1], state_addref2, pbc=pbc, anyon_type=anyon_type) 
    ρ12 = reference_rdm(N, [1,2], state_addref2, pbc=pbc, anyon_type=anyon_type)
    correlation = ee(ρ1) + ee(ρ2) - ee(ρ12)

    return correlation
end

"""
    ref_correlation(N::Int64, state_addref3::Vector{ET}; pbc::Bool=true, anyon_type::Symbol=:Fibo) where {ET}

Calculate spatio-temporal correlation using state with three reference qubits.

# Arguments
- `N::Int64`: System size
- `state_addref3::Vector{ET}`: Quantum state with three reference qubits added
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: Model type
- `spatial::Bool=false`: If true, calculate only spatial correlation
- `temporal::Bool=false`: If true, calculate only temporal correlation

# Returns
- `Float64`: Spatio-temporal correlation measure between two any spacetime points.

Uses reference qubit protocol to measure spatio-temporal correlations at two any spacetime points.
"""
function ref_correlation(N::Int64, state_addref3::Vector{ET}; pbc::Bool=true, anyon_type::Symbol=:Fibo, spatial::Bool=false, temporal::Bool=false) where {ET}
    # Calculate the spatio-temporal correlation I(x₁, x₂, t₁, t₂) between two any spacetime points in a given initial_state
    # In basis, aligned as Ref3 Ref2 Ref1 |ψ_{1,2,...,N}>
    #                 Ref3  |   t₂
    #                  |
    #                  |
    #  Ref1 --------> Ref2      t₁
    #   x₁             x₂

    if spatial # pure spatial correlation, only 2 reference qubits
        @info "Only spatial correlation is calculated."
        spatial_corr = temporal_correlation(N, state_addref3, pbc=pbc, anyon_type=anyon_type)
        return spatial_corr, 0
    elseif temporal # pure temporal correlation, only 2 reference qubits
        @info "Only temporal correlation is calculated."
        temporal_corr = temporal_correlation(N, state_addref3, pbc=pbc, anyon_type=anyon_type)
        return 0, temporal_corr
    else # compute both spatial and temporal correlation, 3-point correlation
        ρ1 = reference_rdm(N, [3], state_addref3, pbc=pbc, anyon_type=anyon_type)
        ρ2 = reference_rdm(N, [2], state_addref3, pbc=pbc, anyon_type=anyon_type) 
        ρ3 = reference_rdm(N, [1], state_addref3, pbc=pbc, anyon_type=anyon_type)
        ρ12 = reference_rdm(N, [2, 3], state_addref3, pbc=pbc, anyon_type=anyon_type)
        ρ23 = reference_rdm(N, [1, 2], state_addref3, pbc=pbc, anyon_type=anyon_type)
    
        spatial_corr = ee(ρ1) + ee(ρ2) - ee(ρ12)
        temporal_corr = ee(ρ2) + ee(ρ3) - ee(ρ23)
    end

    return spatial_corr, temporal_corr
end

function trace_distance(ρ1::AbstractMatrix, ρ2::AbstractMatrix)
    diff = ρ1 - ρ2
    # trace distance = 1/2 * ||ρ1 - ρ2||₁
    return 0.5 * tr(sqrt(diff' * diff))
end

function fidelity(ρ1::AbstractMatrix, ρ2::AbstractMatrix)
    # For density matrix, fidelity defined as F(ρ1,ρ2) = tr(√(√ρ1 * ρ2 * √ρ1))²
    sqrt_ρ1 = sqrt(ρ1)
    F = tr(sqrt(sqrt_ρ1 * ρ2 * sqrt_ρ1))^2
    return real(F)
end

function fidelity(st1::AbstractVector, st2::AbstractVector)
    # For pure states, fidelity defined as F(ψ,φ) = |<ψ|φ>|²
    return abs(dot(st1, st2))^2
end