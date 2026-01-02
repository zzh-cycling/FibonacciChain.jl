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
    anyon_eelis(model::AnyonModel, state::Union{Vector{ET}, Matrix{ET}}) where {ET}

Calculate entanglement entropy profile along the chain for given quantum state or density matrix.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, anyon_type)
- `state::Union{Vector{ET}, Matrix{ET}}`: Quantum state vector or density matrix in anyon basis

# Returns
- `Vector{Float64}`: Entanglement entropy at each bipartition from left to right

# Examples
```jldoctest
julia> using FibonacciChain, LinearAlgebra

julia> N = 6;

julia> model = AnyonModel(FibonacciAnyon(), N, true);

julia> # Create ground state of Fibonacci Hamiltonian
       H = anyon_ham(model);

julia> eigenvals, eigenvecs = eigen(H);

julia> ground_state = eigenvecs[:, 1];

julia> # Calculate entanglement entropy profile
       ee_profile = anyon_eelis(model, ground_state);

julia> length(ee_profile) == N - 1  # Profile has N-1 points
true

julia> all(x -> x ≥ 0, ee_profile)  # All entropies are non-negative
true
```
"""
function anyon_eelis(model::AnyonModel, state::Union{Vector{ET}, Matrix{ET}}) where {ET}
    # Generate ee list for a given state from the left to the right
    N = model.N
    splitlis=Vector(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        if m<= div(N,2)
            subrho=anyon_rdm(model, collect(1:m), state)
            EE_lis[m]=ee(subrho)
        else
            subrho=anyon_rdm(model, collect(m+1:N), state)
            EE_lis[m]=ee(subrho)
        end
    end
    return EE_lis
end

function anyonladder_eelis(model::AnyonModel, state::Vector{ET}) where {ET}
    # Generate ee list for a given state from the left to the right
    N = model.N
    splitlis=Vector(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        if m<= div(N,2)
            subrho=ladderrdm(model, collect(1:m), state)
            EE_lis[m]=ee(subrho)
        else
            subrho=ladderrdm(model, collect(m+1:N), state)
            EE_lis[m]=ee(subrho)
        end
    end
    return EE_lis
end

"""
    translation_matrix(::Type{T}; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}

Generate translation operator matrix for Fibonacci basis states.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `Matrix{Float64}`: Translation matrix mapping each basis state to its translated version
"""
function translation_matrix(model::AnyonModel)
    basis=anyon_basis(model) 
    l = length(basis) 
    Mat=zeros(Float64,(l,l))
    translated_basis = cyclebits.(basis) # Use broadcasting to apply cyclebits to each element in basis
    order = searchsortedfirst.(Ref(basis), translated_basis) # Find the indices of the translated basis in the original basis
    for i in 1:l
        Mat[i, order[i]] += 1.0
    end
    
    return Mat
end

"""
    inversion_matrix(::Type{T}; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}

Generate spatial inversion operator matrix for Fibonacci basis states.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `Matrix{Float64}`: Inversion matrix mapping each basis state to its spatially reflected version
"""
function inversion_matrix(model::AnyonModel)
    basis=anyon_basis(model)
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
function _braidingsq_apply(model::AnyonModel{FibonacciAnyon}, state::T, i::Int) where {N, T <: BitStr{N}}
    # default for PBC system
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert num_digits(T) == N "State length mismatch: expected $(N), got $(num_digits(T))"
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
    if model.pbc
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

function braidingsq_matrix(model::AnyonModel{FibonacciAnyon}, idx::Int) 
    @assert model.pbc || (2 <= idx <= model.N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=anyon_basis(model)
    l=length(basis)
    Bmatrix=zeros(ComplexF64, (l,l))
    for i in 1:l
        outputstate1, outputstate2, output1, output2 = _braidingsq_apply(model, basis[i], idx)
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
function braidingsqmap(model::AnyonModel{FibonacciAnyon}, state::Vector{ET}, idx::Int) where {ET}
    # input a superposition state, and output the braided state
    @assert model.pbc || (2 <= idx <= model.N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=anyon_basis(model)
    l=length(basis)
    @assert l == length(state) "state length is expected to be $(l), but got $(length(state))"
    mapped_state = zeros(ComplexF64, length(state))
    for i in 1:l
        outputstate1, outputstate2, output1, output2 = _braidingsq_apply(model, basis[i], idx)
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

"""
    spatial_correlation(model::AnyonModel, state::Union{Vector{ET}, Matrix{ET}}, site1::Int64, site2::Int64) where {ET}

Calculate mutual information between two sites as spatial correlation measure.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, anyon_type)
- `state::Union{Vector{ET}, Matrix{ET}}`: Quantum state vector or density matrix
- `site1::Int64`: First site index (1 ≤ site1 ≤ N)
- `site2::Int64`: Second site index (1 ≤ site2 ≤ N, site2 ≠ site1)

# Returns
- `Float64`: Mutual information I(1:2) = S(1) + S(2) - S(1,2)

Computes quantum mutual information as measure of spatial correlations.
"""
function spatial_correlation(model::AnyonModel, state::Union{Vector{ET}, Matrix{ET}}, site1::Int64, site2::Int64) where {ET}
    # Calculate the spatial correlation between two sites in a given state. For reference qubit added state, we need reference_rdm. For an initial state without reference qubit, we do not need anything.
    N = model.N
    @assert 1 <= site1 <= N "Site1 index must be in the range [1, $(N)]"
    @assert 1 <= site2 <= N "Site2 index must be in the range [1, $(N)]"
    @assert site1 != site2 "Site1 and Site2 must be different"

    ρ1 = anyon_rdm(model, [site1], state)
    ρ2 = anyon_rdm(model, [site2], state)
    ρ12 = anyon_rdm(model, [site1, site2], state)
    
    correlation = ee(ρ1) + ee(ρ2) - ee(ρ12)

    return correlation
end

"""
    temporal_correlation(model::AnyonModel, state_addref2::Vector{ET}) where {ET}

Calculate temporal correlation using state with two reference qubits.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, anyon_type)
- `state_addref2::Vector{ET}`: Quantum state with two reference qubits added

# Returns
- `Float64`: Temporal correlation measure between time slices

Uses reference qubit protocol to measure temporal correlations at single site.
"""
function temporal_correlation(model::AnyonModel, state_addref2::Vector{ET}) where {ET}
    # Calculate the temporal correlation between two time slices at one site in a given initial_state
    # traceref default true, thus keep reference qubits
    ρ1 = reference_rdm(model, [2], state_addref2)
    ρ2 = reference_rdm(model, [1], state_addref2)
    ρ12 = reference_rdm(model, [1,2], state_addref2)
    correlation = ee(ρ1) + ee(ρ2) - ee(ρ12)

    return correlation
end

"""
    ref_correlation(model::AnyonModel, state_addref3::Vector{ET}; spatial::Bool=false, temporal::Bool=false) where {ET}

Calculate spatio-temporal correlation using state with three reference qubits.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, anyon_type)
- `state_addref3::Vector{ET}`: Quantum state with three reference qubits added
- `spatial::Bool=false`: If true, calculate only spatial correlation
- `temporal::Bool=false`: If true, calculate only temporal correlation

# Returns
- `Float64`: Spatio-temporal correlation measure between two any spacetime points.

Uses reference qubit protocol to measure spatio-temporal correlations at two any spacetime points.
"""
function ref_correlation(model::AnyonModel, state_addref3::Vector{ET}; spatial::Bool=false, temporal::Bool=false) where {ET}
    # Calculate the spatio-temporal correlation I(x₁, x₂, t₁, t₂) between two any spacetime points in a given initial_state
    # In basis, aligned as Ref3 Ref2 Ref1 |ψ_{1,2,...,N}>
    #                 Ref3  |   t₂
    #                  |
    #                  |
    #  Ref1 --------> Ref2      t₁
    #   x₁             x₂
    if spatial # pure spatial correlation, only 2 reference qubits, NEED to note that calculate spatial correlation also needs 2 reference qubits, so use `temporal_correlation` function, the only difference is where the reference qubits are added in the circuit.
        @info "Only spatial correlation is calculated."
        spatial_corr = temporal_correlation(model, state_addref3)
        return spatial_corr, 0
    elseif temporal # pure temporal correlation, only 2 reference qubits
        @info "Only temporal correlation is calculated."
        temporal_corr = temporal_correlation(model, state_addref3)
        return 0, temporal_corr
    else # compute both spatial and temporal correlation, 3-point correlation
        ρ1 = reference_rdm(model, [3], state_addref3)
        ρ2 = reference_rdm(model, [2], state_addref3)
        ρ3 = reference_rdm(model, [1], state_addref3)
        ρ12 = reference_rdm(model, [2, 3], state_addref3)
        ρ23 = reference_rdm(model, [1, 2], state_addref3)

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

function qfi(Ob::Vector{Float64}, state::Vector{T}) where T    
    # Calculate the quantum fisher information, espeically for diagonal operators.
    DeltaOb=state'*(Ob.^2 .*state)-(state'*(Ob.*state))^2
    # Calculate the Quantum Fisher Information
    # For spin 1/2, w/o 4
    F_Q = 4*DeltaOb

    return F_Q
end

function qfi(Ob::Matrix{Float64}, state::Vector{T}) where T
    rho=state*state'
    DeltaOb=tr(rho*Ob^2)-tr(rho*Ob)^2
    # Calculate the Quantum Fisher Information
    # For spin 1/2, w/o 4
    F_Q = 4*DeltaOb

    return F_Q
end

function anti_ferro_order(model::AnyonModel{AT}, ::Type{T}) where {N, T <: BitStr{N}, AT<:AbstractAnyonType}
#param N: Number of sites
#return:  antiferromagnetic order diagonal elements
#The eigenvectors of this operator are going from -N to N, increasing by 2, totally N+1 eigenvectors. Number of each eigenvalues is N choose k, 
#where k is the number of domain walls when we consider total Hilbert space. Defined as sum_i Z_i =1/2 (-1)^(i+1) * Z_i, we aim for spin systems.(S_Z= 1/2 Pauli Z)
    basis = anyon_basis(model)
    l=length(basis)
    anti_ferro = zeros(l)

    mask = bmask(T, collect(2:2:N)...)
    for (idx, str) in enumerate(basis)
        masked_str = flip(str, mask)
        Zi=sum([masked_str...].-1/2)
        anti_ferro[idx] = Zi
    end

    return anti_ferro
end
anti_ferro_order(model::AnyonModel{AT}) where {AT<:AbstractAnyonType} = anti_ferro_order(model, BitStr{model.N, Int})


"""
    mutual_information(model::AnyonModel, subsystems::Tuple{Vector{Int64}, Vector{Int64}}, state::Vector{ET}) where {ET}

Calculate the mutual information between two subsystems.

Computes I(A:B) = S_A + S_B - S_AB where S represents the Von Neumann entropy.
Mutual information quantifies the total correlation between subsystems A and B.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, anyon_type)
- `subsystems::Tuple{Vector{Int64}, Vector{Int64}}`: Tuple of (A_sites, B_sites)
- `state::Vector{ET}`: Quantum state vector

# Returns
- `Float64`: Mutual information I(A:B)

# Example
```julia
model = AnyonModel(FibonacciAnyon(), 10, true)
A_sites = [1, 2, 3]
B_sites = [7, 8, 9]
mi = mutual_information(model, (A_sites, B_sites), psi)
```
"""
function mutual_information(model::AnyonModel, subsystems::Tuple{Vector{Int64}, Vector{Int64}}, state::Vector{ET}) where {ET}
    A, B = subsystems
    # MI formula defined as: I(A:B) = S_A + S_B - S_AB
    # Calculate the reduced density matrices
    ρ_A = anyon_rdm(model, A, state)
    ρ_B = anyon_rdm(model, B, state)
    ρ_AB = anyon_rdm(model, vcat(A, B), state)
    # Calculate the Von Neumann entropies
    S_A = ee(ρ_A)
    S_B = ee(ρ_B)
    S_AB = ee(ρ_AB)
    # Calculate the mutual information
    I_AB = S_A + S_B - S_AB
    return I_AB
    
end



"""
    tri_mutual_information(model::AnyonModel, subsystems::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}}, state::Vector{ET}) where {ET}

Calculate the tripartite mutual information between three subsystems.

Computes I(A:B:C) = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC.
This measures genuine three-partie quantum correlations.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, anyon_type)
- `subsystems::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}}`: Tuple of (A_sites, B_sites, C_sites)
- `state::Vector{ET}`: Quantum state vector

# Returns
- `Float64`: Tripartite mutual information I(A:B:C)

# Example
```julia
model = AnyonModel(FibonacciAnyon(), 10, true)
A_sites = [1, 2]
B_sites = [5, 6]
C_sites = [9, 10]
tmi = tri_mutual_information(model, (A_sites, B_sites, C_sites), psi)
```
"""
function tri_mutual_information(model::AnyonModel, subsystems::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}}, state::Vector{ET}) where {ET}
    A, B, C = subsystems
    # TMI formula defined as: I(A:B:C) = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    
    ρ_A = anyon_rdm(model, A, state)
    ρ_B = anyon_rdm(model, B, state)
    ρ_C = anyon_rdm(model, C, state)

    ρ_AB = anyon_rdm(model, vcat(A,B), state)
    ρ_BC = anyon_rdm(model, vcat(B,C), state)
    ρ_AC = anyon_rdm(model, vcat(A,C), state)
    
    ρ_ABC = anyon_rdm(model, vcat(A,B,C), state)
    
    # Calculate the Von Neumann entropies
    
    S_A = ee(ρ_A)
    S_B = ee(ρ_B)
    S_C = ee(ρ_C)
    S_AB = ee(ρ_AB)
    S_BC = ee(ρ_BC)
    S_AC = ee(ρ_AC)
    S_ABC = ee(ρ_ABC)

    # Calculate the mutual information
    I_ABC = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    
    return I_ABC
end