function measure_basismap(::Type{T}, τ:Float64, state::T, i::Int, pbc::Bool=true, sign::Symbol) where {N, T <: BitStr{N}}
    # default for PBC system
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (:p, :m) "sign must be either :p the plus or :m the minus"
    ϕ = (1+√5)/2
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))
    
    if 2<= i <= N-1
        mask=bmask(T,1,2,3) << (N-i-1)
        str100, str101, str010, str001, str000 = T(4) << (N-i-1), T(5) << (N-i-1), T(2) << (N-i-1), T(1) << (N-i-1), T(0) << (N-i-1)
        if state & mask == str000
            return state, X(state,i), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
        elseif state & mask == str001
            return state, exp(-6im*π/5)
        elseif state & mask == str010
            return state, X(state,i), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
        elseif state & mask == str100
            return state, exp(-6im*π/5)
        elseif state & mask == str101
            return state, exp(-2im*π/5)
        end
    end
    if pbc
        if i == 1 #count from the left
        mask=bmask(T, N, N-1,1)
        str100, str101, str010, str001, str000 = bmask(T,1), bmask(T, N-1, 1), bmask(T, N), bmask(T, N-1), T(0)
            if state & mask == str000
                return state, X(state,i), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
            elseif state & mask == str001
                return state, exp(-6im*π/5)
            elseif state & mask == str010
                return state, X(state,i), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
            elseif state & mask == str100
                return state, exp(-6im*π/5)
            elseif state & mask == str101
                return state, exp(-2im*π/5)
            end
        elseif i == N #count from the left
        mask=bmask(T, N, 2, 1)
        str100, str101, str010, str001, str000 = bmask(T,2), bmask(T, N, 2), bmask(T, 1), bmask(T, N), T(0)
            if state & mask == str000
                return state, X(state,i), exp(-2im*π/5)*ϕ^(-1)+exp(-6im*π/5)*ϕ^(-2), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
            elseif state & mask == str001
                return state, exp(-6im*π/5)
            elseif state & mask == str010
                return state, X(state,i), exp(-2im*π/5)*ϕ^(-2)+exp(-6im*π/5)*ϕ^(-1), (exp(-2im*π/5)-exp(-6im*π/5))*ϕ^(-3/2)
            elseif state & mask == str100
                return state, exp(-6im*π/5)
            elseif state & mask == str101
                return state, exp(-2im*π/5)
            end
        end
    end
end


function measure_matrix(::Type{T}, τ:Float64, idx::Int, pbc::Bool=true) where {N, T <: BitStr{N}}
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=Fibonacci_basis(T, pbc)
    l=length(basis)
    Bmatrix=zeros(ComplexF64, (l,l))
    for i in 1:l
        outcome = braidingsq_basismap(T, basis[i], idx, pbc)
        if length(outcome) == 4
            outputstate1, outputstate2, output1, output2=outcome
            j2=searchsortedfirst(basis, outputstate2)
            Bmatrix[i,i]+=output1
            Bmatrix[i,j2]+=output2
        else
            outputstate, output=outcome
            Bmatrix[i,i]+=output
        end
    end
    
    return Bmatrix
end

function measuremap(::Type{T}, τ:Float64, state::Vector{ET}, idx::Int, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    # input a superposition state, and output the braided state
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=Fibonacci_basis(T, pbc)
    l=length(basis)
    @assert l == length(state) "state length is expected to be $(l), but got $(length(state))"
    mapped_state = zeros(ComplexF64, length(state))
    for i in 1:l
        output = braidingsq_basismap(T, basis[i], idx, pbc)
        if length(output) == 4
            outputstate1, outputstate2, output1, output2=output
            j2=searchsortedfirst(basis, outputstate2)
            mapped_state[i]+=output1*state[i] # outputstate1 is the same as basis[i]
            mapped_state[j2]+=output2*state[i]
        else
            outputstate, output1=output # outputstate is the same as basis[i]
            mapped_state[i]+=output1*state[i]
        end
    end
    
    return mapped_state
end
