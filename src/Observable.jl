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
    Mat=zeros(Float64,(length(basis),length(basis)))
    for (i,n) in enumerate(basis)
        m=cyclebits(n)
        j=searchsortedfirst(basis, m)
        Mat[i,j]=1.0
    end
    
    return Mat
end
translation_matrix(N::Int) = translation_matrix(BitStr{N, Int})

function inversion_matrix(::Type{T}) where {N, T <: BitStr{N}}
    basis=Fibonacci_basis(T)
    l=length(basis)
    Imatrix=zeros((l,l))
    # reversed_basis = map(breflect, basis) # The optimization try of using map function and broadcast
    reversed_basis=similar(basis)
    for i in eachindex(basis)
        reversed_basis[i]=breflect(basis[i])
    end
    # Imatrix[CartesianIndex.(collect(1:length(basis)),searchsortedfirst.(Ref(basis), reversed_basis))].+=1.0
    for i in eachindex(basis)
        output=reversed_basis[i]
        j=searchsortedfirst(basis,output)
        Imatrix[i,j]+=1.0
    end
   
    return Imatrix
end
inversion_matrix(N::Int) = inversion_matrix(BitStr{N, Int})

function braiding(::Type{T}, idx::Int) where {N, T <: BitStr{N}}
    basis=Fibonacci_basis(T)
    l=length(basis)
    Bmatrix=zeros((l,l))
    for i in eachindex(basis)
        output=braidbits(basis[i], idx)
        j=searchsortedfirst(basis, output)
        Bmatrix[i,j]+=1.0
    end
    
    return Bmatrix
end

