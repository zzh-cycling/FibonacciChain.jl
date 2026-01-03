"""
    ladderbraidingsqmap(model::AnyonModel{FibonacciAnyon}, state::Vector{ET}, idx::Int) where {ET}

Apply braiding squared operation to density matrix state in vectorized form. Effectively, it describes the noise induced by inter-two layer chain braiding.

# Arguments
- `model::AnyonModel{FibonacciAnyon}`: Fibonacci anyon model containing system parameters (N, pbc)
- `state::Vector{ET}`: Density matrix state in vectorized form
- `idx::Int`: Site index for braiding operation

# Returns
- `Vector{ET}`: Transformed density matrix state after braiding

Operates on superposition states represented as vectorized density matrices.

# Examples
```jldoctest
julia> using FibonacciChain, LinearAlgebra, BitBasis

julia> N = 4;

julia> model = AnyonModel(FibonacciAnyon(), N, true);

julia> # Create PBC Fibonacci anyon basis 
       basis = anyon_basis(model);

julia> l = length(basis);

julia> ρ_vec = zeros(ComplexF64, l^2);

julia> # Initialize as identity matrix (vectorized)
       for i in 1:l
           ρ_vec[(i-1)*l + i] = 1.0/l;
       end

julia> # Apply braiding at site 2
       ρ_braided = ladderbraidingsqmap(model, ρ_vec, 2);

julia> norm(ρ_braided) > 0  # Should be non-zero
true
```
"""
function ladderbraidingsqmap(model::AnyonModel{FibonacciAnyon}, state::Vector{ET}, idx::Int) where {ET} 
    # input a superposition of basis, and output the braided state
    N = model.N
    @assert model.pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis = anyon_basis(model)
    l = length(basis)
    @assert l^2 == length(state) "state length is expected to be $(l^2), but got $(length(state))"
    
    mapped_state = zeros(ComplexF64, length(state))
    for i in 1:l
        # NOTING that in julia the matrix is column-major order, so we reshape a reduced density matrix to vector, its element will be like a0b0, a1b0, a2b0,,,
        for j in 1:l
            basisi1, basisi2, coefi1, coefi2 = _braidingsq_apply(model, basis[i], idx) 
            basisj1, basisj2, coefj1, coefj2 = _braidingsq_apply(model, basis[j], idx)

            if coefi2 == 0 # i has the same index as i1 i2
                j2=searchsortedfirst(basis, basisj2)
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi1*coefj1
                mapped_state[(i-1)*l+j2]+=state[(i-1)*l+j]*coefi1*coefj2
            elseif coefj2 == 0 # j has the same index as j1 j2
                i2=searchsortedfirst(basis, basisi2)  
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi1*coefj1
                mapped_state[(i2-1)*l+j]+=state[(i-1)*l+j]*coefi2*coefj1
            elseif coefi2 == 0 && coefj2 == 0
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi1*coefj1
            else
                i2=searchsortedfirst(basis, basisi2)
                j2=searchsortedfirst(basis, basisj2)
                # Here noting that the state is a vertorizing density matrix, so the index is i+(j-1)*len, not state[i], state[j]
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi1*coefj1
                mapped_state[(i-1)*l+j2]+=state[(i-1)*l+j]*coefi1*coefj2
                mapped_state[(i2-1)*l+j]+=state[(i-1)*l+j]*coefi2*coefj1
                mapped_state[(i2-1)*l+j2]+=state[(i-1)*l+j]*coefi2*coefj2
            end
        end
    end
    
    return mapped_state
end

"""
    ladderChoi(model::AnyonModel{FibonacciAnyon}, p::Float64, state::Vector{ET}) where {ET}

Apply probabilistic braiding noise channel to density matrix state in vec form.

# Arguments
- `model::AnyonModel{FibonacciAnyon}`: Fibonacci anyon model containing system parameters (N, pbc)
- `p::Float64`: Braiding probability (0 ≤ p ≤ 1)
- `state::Vector{ET}`: Density matrix state in vectorized form

# Returns
- `Vector{ET}`: State after applying Choi map with braiding noise

Implements noise channel: ρ → (1-p)ρ + p B(ρ) where B is braiding operation.
"""
function ladderChoi(model::AnyonModel{FibonacciAnyon}, p::Float64, state::Vector{ET}) where {ET}
    # The PBC anyon relation with basis like:
    #  _1 τ1 _2 τ2 _3 τ3 _4 τ4 _5(1), with _ representing the basis, if PBC, thus head tail _ are connected.
    @assert 0 <= p <= 1 "probability is expected to be in [0, 1], but got $p"
    N = model.N

    if model.pbc
        for i in 2:2:N
            state = (1-p)*state + p*ladderbraidingsqmap(model, state, i)
            state /= norm(state) # normalize the state after each braiding
        end
    else
        for i in 2:2:N-1
            state = (1-p)*state + p*ladderbraidingsqmap(model, state, i)
        end
    end

    return state
end

"""
    ladderrdm(model::AnyonModel, subsystems::Vector{Int64}, state::Vector{ET}) where {ET}

Compute reduced density matrix for specified subsystem from vectorized density matrix.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, anyon_type)
- `subsystems::Vector{Int64}`: Indices of subsystem sites to trace out
- `state::Vector{ET}`: Full density matrix state in vectorized form

# Returns
- `Matrix{Float64}`: Reduced density matrix of the subsystem

For ladder systems where both subsystems have equal dimensions.
"""
function ladderrdm(model::AnyonModel, subsystems::Vector{Int64}, state::Vector{ET}) where {ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    # The disjoin_rdm function need to be careful about the combing order of subsystems, as the order of subsystems in the disjoint basis matters. For example, if input state is 2*3 (counting from the left), the disjoint basis counts from the right, is 3*2. So must ensure the order of subsystems is consistent with the input state.

    # However, in this ladder_rdm function, the order of subsystems doesn't matter, because of two subsystems have the same length.
    return disjoint_rdm(model, model, subsystems, subsystems, state)
end

"""
    laddertranslationmap(model::AnyonModel, state::Vector{ET}) where {ET}

Apply translation operator to vectorized density matrix state.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, anyon_type)
- `state::Vector{ET}`: Density matrix state in vectorized form

# Returns
- `Vector{ET}`: Translated density matrix state

Translates both bra and ket parts of the density matrix consistently.
"""
function laddertranslationmap(model::AnyonModel, state::Vector{ET}) where {ET} 
    # input a superposition state, and output the translated state
    basis = anyon_basis(model)
    l = length(basis)
    @assert l^2 == length(state) "state length is expected to be $(l^2), but got $(length(state))"
    
    translated_basis = cyclebits.(basis) 
    order = searchsortedfirst.(Ref(basis), translated_basis) 
    
    mapped_state = zeros(ComplexF64, length(state))
    for i in 1:l
        for j in 1:l
           mapped_state[(i-1)*l+j] = state[(order[i]-1)*l+order[j]]
        end
    end
    
    return mapped_state
end