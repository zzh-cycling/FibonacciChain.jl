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

function build_extended_basis(k_total::Int, basis::Vector{ET}) where {ET}
    T = BitStr{k_total, Int}
    ref_strings = [T(i) for i in 0:(2^k_total-1)] 
    extended_basis = sort(
        mapreduce(suffix -> process_join(basis, [suffix]),
              vcat,
              ref_strings))

    return extended_basis
end
# process_join([suffix], basis) will give [0000, 0001, 0010, 0011], but process_join(basis, [suffix]) will give [0000, 0001, 0010, 0100]
function add_reference_qubits!(N::Int, state::Vector{ET}, site_idx::Int64; k_new::Int=1, pbc::Bool=true, measure_class::Symbol=:Fibo) where {ET}
    # Add k_new reference qubits to the state at the specified site_idx, and place them to the left part of basis (index N-site_idx+1)
    @assert 1 <= site_idx <= N "Site index must be in the range [1, N]"
    1>= k_new >= 0 || error("k_new must be in [0,1]")

    basis_F = Fibonacci_basis(N, pbc, measure_class=measure_class)
    len_F   = length(basis_F)
    l = length(state)
    

    # old reference qubit number, k_old
    k_old = round(Int, log2(length(state) ÷ len_F))
    length(state) == (2^k_old * len_F) ||
        error("state length is not compatible with (k_old, N), can not deduce k_old from state length")
    @assert k_old >= 0 "k_old must be non-negative, but got $(k_old)"
    
    # new reference prefix
    k_total = k_old + k_new
    new_dim = 2^k_total * len_F
    new_state = zeros(ET, new_dim)
    extended_basis = build_extended_basis(k_total, basis_F)

    inds = [i for i in 1:l*k_new if extended_basis[i][N- site_idx + 1]==0]

    co_inds = setdiff(1:l, inds)
    offset = length(state)   # new k_new qubit starting position

    new_state[inds] = state[inds]
    new_state[co_inds .+ offset] = state[co_inds]

    @show state, inds, co_inds, offset, extended_basis, N- site_idx + 1 + k_total, k_total
    return new_state
end

function spatial_correlation(N::Int64, state::Vector{ET}, site1::Int64, site2::Int64, pbc::Bool=true) where {ET}
    # Calculate the spatial correlation between two sites in a given state
    @assert 1 <= site1 <= length(state) "Site1 index must be in the range [1, length(state)]"
    @assert 1 <= site2 <= length(state) "Site2 index must be in the range [1, length(state)]"
    
    ρ1 = rdm_Fibo(N, [site1], state, pbc)
    ρ2 = rdm_Fibo(N, [site2], state, pbc)
    ρ12 = rdm_Fibo(N, [site1, site2], state, pbc)
    
    correlation = ee(ρ1) + ee(ρ2) - ee(ρ12)

    return correlation
end

function temporal_correlation(N::Int64, state::Vector{ET}, site::Int64, time_slice1::Int64, time_slice2::Int64, pbc::Bool=true) where {ET}
    # Calculate the temporal correlation between two sites in a given state
    @assert 1 <= site <= N "Site index must be in the range [1, N]"
    @assert 1 <= time_slice1 <= length(state) "Time slice 1 index must be in the range [1, length(state)]"
    @assert 1 <= time_slice2 <= length(state) "Time slice 2 index must be in the range [1, length(state)]"

    ρ1 = rdm_Fibo(N, [time_slice1], state, pbc)
    ρ2 = rdm_Fibo(N, [time_slice2], state, pbc) 
    ρ12 = rdm_Fibo(N, [time_slice1, time_slice2], state, pbc)
    correlation = ee(ρ1) + ee(ρ2) - ee(ρ12)

    return correlation
end