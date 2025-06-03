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

function eelis_Fibo_state(N::Int64,splitlis::Vector{Int64},state::Vector{ET},pbc::Bool=true) where {ET}
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        subrho=rdm_Fibo(N, collect(1:splitlis[m]), state, pbc)
        EE_lis[m]=ee(subrho)
    end
    return EE_lis
end

function translation_matrix(::Type{T}) where {N, T <: BitStr{N}}
    basis=Fibonacci_basis(T) 
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

function inversion_matrix(::Type{T}) where {N, T <: BitStr{N}}
    basis=Fibonacci_basis(T)
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

function braiding_basismap(::Type{T}, state::T, i::Int, pbc::Bool=true) where {N, T <: BitStr{N}}
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

function braiding_matrix(::Type{T}, idx::Int, pbc::Bool=true) where {N, T <: BitStr{N}}
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=Fibonacci_basis(T, pbc)
    l=length(basis)
    Bmatrix=zeros(ComplexF64, (l,l))
    for i in 1:l
        if length(braiding_basismap(T, basis[i], idx, pbc)) == 4
            outputstate1, outputstate2, output1, output2=braiding_basismap(T, basis[i], idx, pbc)
            j1=searchsortedfirst(basis, outputstate1)
            j2=searchsortedfirst(basis, outputstate2)
            Bmatrix[i,j2]+=output2
            Bmatrix[i,j1]+=output1
        else
            outputstate, output=braiding_basismap(T, basis[i], idx, pbc)
            j=searchsortedfirst(basis, outputstate)
            Bmatrix[i,j]+=output
        end
    end
    
    return Bmatrix
end

function braidingmap(::Type{T}, state::Vector{ET}, idx::Int, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    # input a superposition state, and output the braided state
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=Fibonacci_basis(T, pbc)
    l=length(basis)
    mapped_state = zeros(ComplexF64, length(state))
    for i in 1:l
        output = braiding_basismap(T, basis[i], idx, pbc)
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
braidingmap(N::Int, state::Vector{ET}, idx::Int, pbc::Bool=true) where {ET} = braidingmap(BitStr{N, Int}, state, idx, pbc)

function ladderbraidingmap(::Type{T}, state::Vector{ET}, idx::Int, pbc::Bool=true) where {N, T <: BitStr{N}, ET} 
    # input a superposition state, and output the braided state
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=Fibonacci_basis(T, pbc)
    len=length(basis)
    @assert len^2 == length(state) "state length is expected to be $(len^2), but got $(length(state))"
    
    mapped_state = zeros(ComplexF64, length(state))
    for i in 1:len
        for j in 1:len
            output1 = braiding_basismap(T, basis[i], idx, pbc)
            output2 = braiding_basismap(T, basis[j], idx, pbc)
            if length(output1) == 4 && length(output2) == 4
                basisi1, basisi2, coefi1, coefi2=output1
                basisj1, basisj2, coefj1, coefj2=output2
                i2=searchsortedfirst(basis, basisi2)
                j2=searchsortedfirst(basis, basisj2)
                mapped_state[(i-1)*len+j]+=state[i]*state[j]*coefi1*coefj1
                mapped_state[(i-1)*len+j2]+=state[i]*state[j]*coefi1*coefj2
                mapped_state[(i2-1)*len+j]+=state[i]*state[j]*coefi2*coefj1
                mapped_state[(i2-1)*len+j2]+=state[i]*state[j]*coefi2*coefj2
            elseif length(output1) == 4 && length(output2) == 2
                basisi1, basisi2, coefi1, coefi2=output1
                basisj, coefj=output2
                i2=searchsortedfirst(basis, basisi2)  
                mapped_state[(i-1)*len+j]+=state[i]*state[j]*coefi1*coefj
                mapped_state[(i2-1)*len+j]+=state[i]*state[j]*coefi2*coefj
            elseif length(output1) == 2 && length(output2) == 4
                basisi, coefi=output1
                basisj1, basisj2, coefj1, coefj2=output2
                j2=searchsortedfirst(basis, basisj2)
                mapped_state[(i-1)*len+j]+=state[i]*state[j]*coefi*coefj1
                mapped_state[(i-1)*len+j2]+=state[i]*state[j]*coefi*coefj2
            else
                basisi, coefi=output1
                basisj, coefj=output2
                mapped_state[(i-1)*len+j]+=state[i]*state[j]*coefi*coefj
            end
        end
    end
    
    return mapped_state
end
ladderbraidingmap(N::Int, state::Vector{ET}, idx::Int, pbc::Bool=true) where {ET} = ladderbraidingmap(BitStr{N, Int}, state, idx, pbc)

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