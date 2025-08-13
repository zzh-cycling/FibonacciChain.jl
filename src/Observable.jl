```

    ee(::Matrix{ET}) where {ET} -> Float64

# params: a reduced density matrix `subrm` of type `ET`
# Return the entanglement entropy of the reduced density matrix.
```
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

```
    eelis_Fibo_state(N::Int64,state::Vector{ET},pbc::Bool=true; measure_class::Symbol=:Fibo) where {ET} -> Vector{Float64}  

# params: a int of lattice number, `state` is the state vector in anyon basis, `pbc` is a boolean value, `measure_class` is a symbol, which default to be :Fibo
# Return the entanglement entropy list for a given state from the left to the right.

```
function eelis_Fibo_state(N::Int64,state::Vector{ET},pbc::Bool=true; measure_class::Symbol=:Fibo) where {ET}
    # Generate ee list for a given state from the left to the right
    splitlis=Vector(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        if m<= div(N,2)
            subrho=rdm_Fibo(N, collect(1:m), state, pbc, measure_class=measure_class)
            EE_lis[m]=ee(subrho)
        else
            subrho=rdm_Fibo(N, collect(m+1:N), state, pbc, measure_class=measure_class)
            EE_lis[m]=ee(subrho)
        end
    end
    return EE_lis
end

function eelis_Fiboladder_state(N::Int64,state::Vector{ET},pbc::Bool=true) where {ET}
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

```
    translation_matrix(::Type{T}) where {N, T <: BitStr{N}} -> Matrix{Float64}

# params: a type `T` of BitStr{N, Int}
# Return the translation matrix for the Fibonacci basis, which is used to translate the state from one site to another.
```
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

```
    inversion_matrix(::Type{T}) where {N, T <: BitStr{N}} -> Matrix{Float64}
# params: a type `T` of BitStr{N, Int}
# Return the inversion matrix for the Fibonacci basis, which is used to invert the state from one site to another.
```
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

```
    braidingsqmap(::Type{T}, state::Vector{ET}, idx::Int, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}, ET} -> Vector{ET}
    
 # params: a type `T` of BitStr{N, Int}, `state` is the state vector in anyon basis, `idx` is the index of the site to be braided, `pbc` is a boolean value, `measure_class` is a symbol, which default to be :Fibo

# Return the braided state vector. In anyon basis, braiding is a fundamental operation that changes the state of the system by braiding the anyons at the specified index.
```
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

function braidingsq_matrix(::Type{T}, idx::Int, pbc::Bool=true) where {N, T <: BitStr{N}}
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

```
    braidingsqmap(::Type{T}, state::Vector{ET}, idx::Int, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}, ET} -> Vector{ET}

# params: a type `T` of BitStr{N, Int}, `state` is the state vector in anyon basis, `idx` is the index of the site to be braided, `pbc` is a boolean value, `measure_class` is a symbol, which default to be :Fibo
# Return the braided state vector. In anyon basis, braiding is a fundamental operation that changes the state of the system by braiding the anyons at the specified index.
```
function braidingsqmap(::Type{T}, state::Vector{ET}, idx::Int, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
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
braidingsqmap(N::Int, state::Vector{ET}, idx::Int, pbc::Bool=true) where {ET} = braidingsqmap(BitStr{N, Int}, state, idx, pbc)

```
    spatial_correlation(N::Int64, state::Union{Vector{ET}, Matrix{ET}}, site1::Int64, site2::Int64; pbc::Bool=true, measure_class::Symbol=:Fibo) where {ET} -> Float64

# Calculate the spatial correlation between two sites in a given state. For reference qubit added state, we need reference_rdm. For an initial state without reference qubit, we do not need anything.
# params: a int of lattice number, `state` is the state vector in anyon basis, `site1` and `site2` are the indices of the sites to be correlated, `pbc` is a boolean value, `measure_class` is a symbol, which default to be :Fibo
# Return the spatial correlation between the two sites.
```
function spatial_correlation(N::Int64, state::Union{Vector{ET}, Matrix{ET}}, site1::Int64, site2::Int64; pbc::Bool=true, measure_class::Symbol=:Fibo) where {ET}
    # Calculate the spatial correlation between two sites in a given state. For reference qubit added state, we need reference_rdm. For an initial state without reference qubit, we do not need anything.
    @assert 1 <= site1 <= N "Site1 index must be in the range [1, $(N)]"
    @assert 1 <= site2 <= N "Site2 index must be in the range [1, $(N)]"
    @assert site1 != site2 "Site1 and Site2 must be different"

    ρ1 = rdm_Fibo(N, [site1], state, pbc, measure_class=measure_class)
    ρ2 = rdm_Fibo(N, [site2], state, pbc, measure_class=measure_class)
    ρ12 = rdm_Fibo(N, [site1, site2], state, pbc, measure_class=measure_class)
    
    correlation = ee(ρ1) + ee(ρ2) - ee(ρ12)

    return correlation
end

```
    temporal_correlation(N::Int64, state_addref2::Vector{ET}; pbc::Bool=true, measure_class::Symbol=:Fibo) where {ET} -> Float64

# Calculate the temporal correlation between two time slices at one site in a given initial_state
# params: a int of lattice number, `state_addref2` is the state vector
# Return the temporal correlation between two time slices at one site in a given initial_state
```
function temporal_correlation(N::Int64, state_addref2::Vector{ET}; pbc::Bool=true, measure_class::Symbol=:Fibo) where {ET}
    # Calculate the temporal correlation between two time slices at one site in a given initial_state

    ρ1 = reference_rdm(N, [2], state_addref2, pbc=pbc, measure_class=measure_class)
    ρ2 = reference_rdm(N, [1], state_addref2, pbc=pbc, measure_class=measure_class) 
    ρ12 = reference_rdm(N, [1,2], state_addref2, pbc=pbc, measure_class=measure_class)
    correlation = ee(ρ1) + ee(ρ2) - ee(ρ12)

    return correlation
end