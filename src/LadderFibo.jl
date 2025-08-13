"""
    ladderbraidingsqmap(::Type{T}, state::Vector{ET}, idx::Int, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}, ET} -> Vector{ET}

# params: a type `T` of BitStr{N, Int}, `state` is the state vector in anyon basis, `idx` is the index of the site to be braided, `pbc` is a boolean value, `measure_class` is a symbol, which default to be :Fibo
# Return the braided state vector. In anyon basis, braiding is a fundamental operation that changes the state of the system by braiding the anyons at the specified index.
"""
function ladderbraidingsqmap(::Type{T}, state::Vector{ET}, idx::Int, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}, ET} 
    # input a superposition of basis, and output the braided state
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=Fibonacci_basis(T, pbc, measure_class=measure_class)
    l=length(basis)
    @assert l^2 == length(state) "state length is expected to be $(l^2), but got $(length(state))"
    
    mapped_state = zeros(ComplexF64, length(state))
    for i in 1:l
        # NOTING that in julia the matrix is column-major order, so we reshape a reduced density matrix to vector, its element will be like a0b0, a1b0, a2b0,,,
        for j in 1:l
            output1 = braidingsq_basismap(T, basis[i], idx, pbc)
            output2 = braidingsq_basismap(T, basis[j], idx, pbc)
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
ladderbraidingsqmap(N::Int, state::Vector{ET}, idx::Int, pbc::Bool=true; measure_class::Symbol=:Fibo) where {ET} = ladderbraidingsqmap(BitStr{N, Int}, state, idx, pbc, measure_class=measure_class)

"""
    ladderChoi(::Type{T}, p::Float64, state::Vector{ET}, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N,T <: BitStr{N}, ET} -> Vector{ET}    

# params: a type `T` of BitStr{N, Int}, `p` is the probability of braiding, `state` is the state vector in anyon basis, `pbc` is a boolean value, which default to be true, `measure_class` is a symbol, which default to be :Fibo
# The PBC anyon relation with basis like:
#  _1 τ1 _2 τ2 _3 τ3 _4 τ4 _5(1), with _ representing the basis, if PBC, thus head tail _ are connected.
# Return the state vector after braiding the anyons in the state with a given probability.
"""
function ladderChoi(::Type{T}, p::Float64, state::Vector{ET}, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N,T <: BitStr{N}, ET}
    # The PBC anyon relation with basis like:
    #  _1 τ1 _2 τ2 _3 τ3 _4 τ4 _5(1), with _ representing the basis, if PBC, thus head tail _ are connected.
    @assert 0 <= p <= 1 "probability is expected to be in [0, 1], but got $p"

    if pbc
        for i in 2:2:N
            state=(1-p)*state+p*ladderbraidingsqmap(T, state, i, pbc, measure_class=measure_class)
            state/=norm(state) # normalize the state after each braiding
        end
    else
        for i in 2:2:N-1
            state=(1-p)*state+p*ladderbraidingsqmap(T, state, i, pbc, measure_class=measure_class)
        end
    end

    return state
end
ladderChoi(N::Int, probability::Float64, state::Vector{ET}, pbc::Bool=true; measure_class::Symbol=:Fibo) where {ET} = ladderChoi(BitStr{N, Int}, probability, state, pbc, measure_class=measure_class)

"""
    ladderrdm(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N,T <: BitStr{N}, ET} -> Matrix{Float64}

# params: a type `T` of BitStr{N, Int}, `subsystems` is the indices of the subsystems to be reduced, `state` is the state vector in anyon basis, `pbc` is a boolean value, which default to be true, `measure_class` is a symbol, which default to be :Fibo
# Usually subsystem indices count from the right of binary string.
# The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.
# The disjoin_rdm function need to be careful about the combing order of subsystems, as the order of subsystems in the disjoint basis matters. For example, if input state is 2*3 (counting from the left), the disjoint basis counts from the right, is 3*2. So must ensure the order of subsystems is consistent with the input state.
# However, in this ladder_rdm function, the order of subsystems doesn't matter, because of two subsystems have the same length.
# Return the reduced density matrix of the subsystems in the state.
"""
function ladderrdm(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N,T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    # The disjoin_rdm function need to be careful about the combing order of subsystems, as the order of subsystems in the disjoint basis matters. For example, if input state is 2*3 (counting from the left), the disjoint basis counts from the right, is 3*2. So must ensure the order of subsystems is consistent with the input state.

    # However, in this ladder_rdm function, the order of subsystems doesn't matter, because of two subsystems have the same length.
    return disjoint_rdm(T, T, subsystems, subsystems, state, pbc)
end
ladderrdm(N::Int, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {ET} = ladderrdm(BitStr{N, Int}, subsystems, state, pbc)

"""
    laddertranslationmap(::Type{T}, state::Vector{ET}) where {N, T <: BitStr{N}, ET} -> Vector{ET}
# params: a type `T` of BitStr{N, Int}, `state` is the state vector in anyon basis
# input a superposition state, and output the translated state
"""
function laddertranslationmap(::Type{T}, state::Vector{ET}) where {N, T <: BitStr{N}, ET} 
    # input a superposition state, and output the translated state
    basis=Fibonacci_basis(T)
    l=length(basis)
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
laddertranslationmap(N::Int, state::Vector{ET}) where {ET} = laddertranslationmap(BitStr{N, Int}, state)