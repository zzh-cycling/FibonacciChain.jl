function Fibonacci_chain_OBC(::Type{T}) where {N, T <: BitStr{N}}
    # Generate Fibonacci chain for Fibonacci model with open boundary condition
    fib_chain=[[T(0), T(1)],[T(0), T(1), T(2)]]
    for i in 3:N
        push!(fib_chain,vcat([s << 1 for s in fib_chain[i-1]],[(s << 2 | T(1)) for s in fib_chain[i-2]]))
    end
    # each push we add a bit"0" or bit"01" to the end of the bit_string, and finally return the last element of the fib_chain
    return fib_chain[N]
end

function Fibonacci_chain_PBC(::Type{T}) where {N, T <: BitStr{N}}
    # Generate Fibonacci chain  for Fibonacci model with periodic boundary condition
    return filter(c -> iszero((c >> (N-1)) & (c & 1)), Fibonacci_chain_OBC(T))
end

function Fibonacci_basis(::Type{T},pbc::Bool=true) where {N, T <: BitStr{N}}
    # Generate basis for Fibonacci model, return both decimal and binary form, where we both consider PBC and OBC
    if pbc
        basis=Fibonacci_chain_PBC(T)
    else
        basis=Fibonacci_chain_OBC(T)
    end
    sorted_basis=sort(basis)
    return sorted_basis
end
Fibonacci_basis(N::Int, pbc::Bool=true) = Fibonacci_basis(BitStr{N, Int}, pbc)

function antimap(::Type{T}, state::T, i::Int) where {N, T <: BitStr{N}}
    # The type of n is DitStr{D, N, Int}, which is a binary string with length N in D-ary form.
    # Acting Hamiltonian on a given state in bitstr and return the output (states, weight) in bitstr
    # Here need to note that the order of the bitstr is from right to left, which is different from our counting order.
    ϕ = (1+√5)/2
    fl=bmask(T, N)

    X(state,i) = flip(state, fl >> (i-1))

    if (state & (1 << (N-i))) == 0
        return state, X(state,i), -ϕ^(-1), -ϕ^(-3/2)
    else
        return state, X(state,i), -ϕ^(-2), -ϕ^(-3/2)
    end
end

function ferromap(::Type{T}, state::T, i::Int) where {N, T <: BitStr{N}}
    ϕ = (1+√5)/2
    fl=bmask(T, N)

    X(state,i) = flip(state, fl >> (i-1))

    if (state & (1 << (N-i))) == 0
        return state, X(state,i), ϕ^(-1), ϕ^(-3/2)
    else
        return state, X(state,i), ϕ^(-2), ϕ^(-3/2)
    end
end

function count_subBitStr(::Type{T}, state::T) where {N, T <: BitStr{N}}
    n = length(state)
    n < 3 && return 0 
    
    str100, str101, str001 = T(4), T(5), T(1) # 100, 101, 001
    num=0
    
    mask=bmask(T, 1, 2, 3)
    for i in 1:(n-2) # start from string right to left
        substr = state & (mask << (i-1))  
        if substr == str101
            num+= 1
        end
        str101 <<= 1
    end
    
    return num
end

function actingHam(::Type{T}, state::T, pbc::Bool=true) where {N, T <: BitStr{N}}
    # The type of n is DitStr{D, N, Int}, which is a binary string with length N in D-ary form.
    # Acting Hamiltonian on a given state in bitstr and return the output states in bitstr
    # Here need to note that the order of the bitstr is from right to left, which is different from our counting order.
    mask=bmask(T, N, N-2)
    fl=bmask(T, N)
    ϕ = (1+√5)/2
    X(state,i) = flip(state, fl >> (i-1))

    output = Dict{T, Float64}()

    # count 101, 100, 001
    output[state] = get(output, state, 0.0) - count_subBitStr(T, state)

    # start from 2 site to N-1 site to count 0x0, because the first and last bits are not considered
    for i in 2:N-1 
        if state & (mask >> (i-2)) == 0
            state1, state2, weight1, weight2 = antimap(T, state, i)
            output[state1] = get(output, state1, 0.0) + weight1
            output[state2] = get(output, state2, 0.0) + weight2
        end
    end

    if pbc
        # 1 site antimap
        if state[1]==0 && state[N-1]==0
            state1, state2, weight1, weight2 = antimap(T, state, 1)
            output[state1] = get(output, state1, 0.0) + weight1
            output[state2] = get(output, state2, 0.0) + weight2
        end
        # N site antimap
        if state[2]==0 && state[N]==0
            state1, state2, weight1, weight2 = antimap(T, state, N)
            output[state1] = get(output, state1, 0.0) + weight1
            output[state2] = get(output, state2, 0.0) + weight2
        end
        mask1= bmask(T, N, 2)
        mask2= bmask(T, N-1, 1)
        # 1 site 111 fusion
        if state & mask1 == mask1
            output[state] = get(output, state, 0.0) - 1
        end
        # N site 111 fusion
        if state & mask2 == mask2
            output[state] = get(output, state, 0.0) - 1
        end
    end
    return output
end

function ferroactingHam(::Type{T}, state::T, pbc::Bool=true) where {N, T <: BitStr{N}}
    mask=bmask(T, N, N-2)
    fl=bmask(T, N)
    ϕ = (1+√5)/2
    X(state,i) = flip(state, fl >> (i-1))

    output = Dict{T, Float64}()

    # count 101, 100, 001
    output[state] = get(output, state, 0.0) + count_subBitStr(T, state)

    # start from 2 site to N-1 site to count 0x0, because the first and last bits are not considered
    for i in 2:N-1 
        if state & (mask >> (i-2)) == 0
            state1, state2, weight1, weight2 = ferromap(T, state, i)
            output[state1] = get(output, state1, 0.0) + weight1
            output[state2] = get(output, state2, 0.0) + weight2
        end
    end

    if pbc
        # 1 site ferromap
        if state[1]==0 && state[N-1]==0
            state1, state2, weight1, weight2 = ferromap(T, state, 1)
            output[state1] = get(output, state1, 0.0) + weight1
            output[state2] = get(output, state2, 0.0) + weight2
        end
        # N site ferromap
        if state[2]==0 && state[N]==0
            state1, state2, weight1, weight2 = ferromap(T, state, N)
            output[state1] = get(output, state1, 0.0) + weight1
            output[state2] = get(output, state2, 0.0) + weight2
        end
        mask1= bmask(T, N, 2)
        mask2= bmask(T, N-1, 1)
        # 1 site 111 fusion
        if state & mask1 == mask1
            output[state] = get(output, state, 0.0) + 1
        end
        # N site 111 fusion
        if state & mask2 == mask2
            output[state] = get(output, state, 0.0) + 1
        end
    end
    return output
end

function Fibonacci_Ham(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    # Generate Hamiltonian for Fibonacci model, automotically contain pbc or obc
    basis=Fibonacci_basis(T,pbc)

    l=length(basis)
    H=zeros(Float64,(l,l))
    for i in 1:l
        output=actingHam(T, basis[i], pbc) 
        states, weights = keys(output), values(output)
        for m in states
            j=searchsortedfirst(basis, m)
            H[i, j] += output[m]
        end
    end

    return H
end
Fibonacci_Ham(N::Int, pbc::Bool=true) = Fibonacci_Ham(BitStr{N, Int}, pbc)

function Fibonacci_ferroHam(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    # Generate Hamiltonian for Fibonacci model, automotically contain pbc or obc
    basis=Fibonacci_basis(T,pbc)

    l=length(basis)
    H=zeros(Float64,(l,l))
    for i in 1:l
        output=ferroactingHam(T, basis[i], pbc) 
        states, weights = keys(output), values(output)
        for m in states
            j=searchsortedfirst(basis, m)
            H[i, j] += output[m]
        end
    end

    return H
end
Fibonacci_ferroHam(N::Int, pbc::Bool=true) = Fibonacci_ferroHam(BitStr{N, Int}, pbc)

# join two lists of basis by make a product of two lists
function process_join(a, b)
    return vec([join(b, a) for a in a, b in b])
end

# create Fibonacci basis composed of multiple disjoint sub-chains
function joint_Fibo_basis(lengthlis::Vector{Int})
    return sort(mapreduce(len -> Fibonacci_basis(len, false), process_join, lengthlis))
end

function connected_components(v::Vector{Int})
    if isempty(v)
        return []
    end

    sort!(v)

    result = []
    current_segment = [v[1]]

    for i in 2:length(v)
        if v[i] == v[i - 1] + 1
            push!(current_segment, v[i])
        else
            push!(result, current_segment)
            current_segment = [v[i]]
        end
    end

    push!(result, current_segment)

    return result
end

function move_subsystem(::Type{BitStr{M, INT}}, basis::BitStr{N, INT}, subsystems::Vector{Int}) where {M, N, INT}
    @assert length(subsystems) == N "subsystems length is expected to be $N, but got $(length(subsystems))"
    @assert M >= N "total length is expected to be greater than or equal to $N, but got $M"
    return sum(i -> BitStr{M}(readbit(basis.buf, i) << (subsystems[i] - 1)), 1:N)
end

# take environment part of a basis
takeenviron(x, mask::BitStr{l}) where {l} = x & (~mask)
# take system part of a basis
takesystem(x, mask::BitStr{l}) where {l} = (x & mask)

function rdm_Fibo(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {N,T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    unsorted_basis = Fibonacci_basis(T, pbc)
    @assert length(unsorted_basis) == length(state) "state length is expected to be $(length(unsorted_basis)), but got $(length(state))"
    
    subsystems=connected_components(subsystems)
    lengthlis=length.(subsystems)
    subsystems=vcat(subsystems...)
    mask = bmask(T, subsystems...)

    
    order = sortperm(unsorted_basis, by = x -> (takeenviron(x, mask), takesystem(x, mask))) #first sort by environment, then by system. The order of environment doesn't matter.
    basis, state = unsorted_basis[order], state[order]
    
    reduced_basis = move_subsystem.(T, joint_Fibo_basis(lengthlis), Ref(subsystems))
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

    for i in 1:length(result_indices)-1
        range = result_indices[i]:result_indices[i+1]-1         
        # Get indices in the reduced basis
        indices = searchsortedfirst.(Ref(reduced_basis), takesystem.(basis[range], mask))
        view(reduced_dm, indices, indices) .+= view(state, range) .* view(state, range)'
    end

    return reduced_dm
end
rdm_Fibo(N::Int, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {ET} = rdm_Fibo(BitStr{N, Int}, subsystems, state, pbc)

function iso_total2K(::Type{T}, k::Int64) where {N, T <: BitStr{N}}
    #Function to map the total basis to the K space basis, actually is the isometry, defined as W'*W=I, W*W'=P, P^2=P
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basis = Fibonacci_basis(T)

    k_dic = Dict{Int, Vector{Int64}}()
    basisK = Vector{T}(undef, 0)
    for i in eachindex(basis)
        state=basis[i]
        category = get_representative(state)[1]
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

    iso = zeros((length(basis), length(keys(basisK))))
    
    for (i, state) in enumerate(basisK)
        state_indices = k_dic[state]  
        l = length(state_indices)
        iso[state_indices, i] .= 1/sqrt(l)
    end

    return iso
end
iso_total2K(N::Int, k::Int64) = iso_total2K(BitStr{N, Int}, k)

function mapstate_K2total(::Type{T}, state::Vector{ET}, k::Int64) where {N, T <: BitStr{N}, ET}
    # Map the K space state to total space state
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basis = Fibonacci_basis(T)
    k_dic = Dict{Int, Vector{Int64}}()
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
mapstate_K2total(N::Int, state::Vector{ET}, k::Int64) where {ET} = mapstate_K2total(BitStr{N, Int}, state, k)

function rdm_Fibo_K(::Type{T}, subsystems::Vector{Int64},kstate::Vector{ET}, k::Int64) where {N,T <: BitStr{N}, ET}
    @assert length(kstate) == length(Fibo_K_basis(T,k)[1]) "state length is expected to be $(length(Fibo_K_basis(T, k)[1])), but got $(length(kstate))"
    state = mapstate_K2total(T, kstate, k)
    reduced_dm = rdm_Fibo(T, subsystems, state)
    return reduced_dm
end
rdm_Fibo_K(N::Int, subsystems::Vector{Int64},state::Vector{ET}, k::Int64) where {ET} = rdm_Fibo_K(BitStr{N, Int}, subsystems, state, k)


function iso_K2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
#Function to map the MSS basis to the K space basis
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"
    basisK, k_dic = Fibo_K_basis(T, k)

    MSS_dic = Dict{Int, Vector{Int64}}()
    # MSS_dic is a dictionary, the key is the representative state of the inversion of n, and the value is the index of the state in the basisK. NOTE that MSS_dic is not sorted, so we need to sort it later.
    qlist = Vector{Int}(undef, 0)
    # Below procedure is to collapse the extra basis in K space that can be converted mutually to MSS space.
    if inv==1 && k==0 || inv==-1 && k==div(N,2)
        for i in eachindex(basisK)
            n = basisK[i]
            # here we calculate the representative state of the inversion of n
            nR = get_representative(breflect(n))[1]
            # For example, n = 41, nR=37, then we only need to keep n=37, and n=41 will be removed.
            if n <= min(nR, n)
                push!(qlist, length(Set([n, nR])))
            end
            n = min(nR, n)
            if haskey(MSS_dic, n)
                push!(MSS_dic[n], i)
            else
                MSS_dic[n] = [i]
            end
        end

    else
        for i in eachindex(basisK)
            n = basisK[i]
            nR = get_representative(breflect(n))[1]
            if n != nR
                n = min(nR, n)
                if haskey(MSS_dic, n)
                    push!(MSS_dic[n], i)
                else
                    MSS_dic[n] = [i]
                end
                push!(qlist, 2)
            end
            
        end    
        
    
    end

    iso = zeros((length(basisK), length(MSS_dic)))
    MSS_dic=sort(MSS_dic)
    for (i, state_index) in enumerate(values(MSS_dic))
        iso[state_index, i] .= 1/sqrt(qlist[i])
    end

    return iso
end
iso_K2MSS(N::Int, k::Int64, inv::Int64=1) = iso_K2MSS(BitStr{N, Int}, k, inv)

function mapstate_MSS2K(::Type{T}, state::Vector{ET}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}, ET}
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"

    basisK, k_dic = Fibo_K_basis(T, k)

    MSS_dic = Dict{Int, Vector{Int64}}()
    qlist = Vector{Int}(undef, 0)
   
    if inv==1
        for i in eachindex(basisK)
            n = basisK[i]
            nR = get_representative(breflect(n))[1]
            # For example, n = 41, nR=37, then we only need to keep n=37, and n=41 will be removed.
            if n <= min(nR, n)
                push!(qlist, length(Set([n, nR])))
            end
            n = min(nR, n)
            if haskey(MSS_dic, n)
                push!(MSS_dic[n], i)
            else
                MSS_dic[n] = [i]
            end
        end

    elseif inv==1 && k==div(N,2) || inv==-1 && k==0
        for i in eachindex(basisK)
            n = basisK[i]
            nR = get_representative(breflect(n))[1]
            if n != nR
                n = min(nR, n)
                if haskey(MSS_dic, n)
                    push!(MSS_dic[n], i)
                else
                    MSS_dic[n] = [i]
                end
                push!(qlist, 2)
            end
            
        end    
        
    
    end

    total_state = zeros(ET, length(basisK))
    MSS_dic=sort(MSS_dic)
    for (i, state_index) in enumerate(values(MSS_dic))
        total_state[state_index] .= 1/sqrt(qlist[i])*state[i]
    end

    return total_state
end
mapstate_MSS2K(N::Int, state::Vector{ET}, k::Int64, inv::Int64=1) where {ET} = mapstate_MSS2K(BitStr{N, Int}, state, k, inv)

mapstate_MSS2total(::Type{T}, state::Vector{ET}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}, ET} = mapstate_K2total(T, mapstate_MSS2K(T, state, k, inv), k)
mapstate_MSS2total(N::Int64, state::Vector{ET}, k::Int64, inv::Int64=1) where {ET} = mapstate_MSS2total(BitStr{N, Int}, state, k, inv)

function iso_total2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
    # Function to map the total basis to the MSS space basis, k can only equal to 0 or N/2(pi)
    iso = iso_total2K(T, k) * iso_K2MSS(T, k, inv)

    return iso
end
iso_total2MSS(N::Int, k::Int64, inv::Int64=1) = iso_total2MSS(BitStr{N, Int}, k, inv)

function rdm_Fibo_MSS(::Type{T}, subsystems::Vector{Int64}, mssstate::Vector{ET}, k::Int64, inv::Int64=1) where {N,T <: BitStr{N}, ET}
    @assert length(Fibo_MSS_basis(T, k, inv)[1]) == length(mssstate) "state length is expected to be $(length(Fibo_MSS_basis(T, k, inv)[1])), but got $(length(mssstate))"
    state=mapstate_MSS2total(T, mssstate, k, inv)
    reduced_dm = rdm_Fibo(T, subsystems, state)
    return reduced_dm
end
rdm_Fibo_MSS(N::Int64, subsystems::Vector{Int64},state::Vector{ET}, k::Int64, inv::Int64=1) where {ET} = rdm_Fibo_MSS(BitStr{N, Int}, subsystems, state, k, inv)


function cyclebits(state::T) where {N, T <: BitStr{N}}
    #params: t is an integer, N is the length of the binary string
    #We also use this order: system size, state, circular shift bitstring 1 bit.
    # In case need to shift more than 1 bit, we can use a loop or recursion. or we leave a interface here  n_translations::Int
    mask = 1 << N - 1
    return ((state << 1) | (state >> (N - 1))) & mask
end

function get_representative(state::T) where {N, T <: BitStr{N}}
#Finds representative and representative translation for a state.
#State should be a decimal integer.

    representative = state
    translation = 0
    # cycle(bits) = (bits << 1) % (2^N - 1)  # Left shift the state by one position
    current = state
    for n_translation_sites in 1:N-1
        current = cyclebits(current)  # Cycle the bits
        if current < representative
            representative = current
            translation = n_translation_sites
        end
    end

    return representative, translation
end

function Fibo_K_basis(::Type{T}, k::Int64) where {N, T <: BitStr{N}}
#params: a int of lattice number and momentum of system
#return: computational basis in given momentum kinetically constrained subspace with decimal int form in golden chain model
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basisK = Vector{T}(undef, 0)
    basis = Fibonacci_basis(T)


    basis_dic = Dict{T, Vector{T}}()
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
Fibo_K_basis(N::Int, k::Int64) = Fibo_K_basis(BitStr{N, Int}, k)


function Fib_MSS_basis(::Type{T}, k::Int64,inv::Int64=1) where {N, T <: BitStr{N}}
#params: a int of lattice number and momentum of system, we have considered the inversion symmetry
#return: computational basis in given momentum inversion symmetry subspace with decimal int form
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"
    # MSS is the list of states in the maximum symmetry sector
    MSS = Vector{T}(undef, 0)
    basisK, basis_dic = Fibo_K_basis(T, k)
    MSS_dic = Dict{T, Vector{T}}()

    # q is the number of states that are equivalent under inversion
    qlist = Vector{Int}(undef, 0)
    if inv==1 && k==0 || inv==-1 && k==div(N,2)
        for i in eachindex(basisK)
            n = basisK[i]
            # here we calculate the representative state of the inversion of n
            nR = get_representative(breflect(n))[1]
            if n <= min(nR, n)
                push!(MSS, n)
                MSS_dic[n] = basis_dic[n]
                push!(qlist, length(Set([n, nR])))
            end
        end

        return MSS, MSS_dic, qlist
        
    elseif inv==1 && k==div(N,2) || inv==-1 && k==0
        for i in eachindex(basisK)
                n = basisK[i]
                nR = get_representative(breflect(n))[1]
                if n <= min(nR, n)
                    push!(MSS, n)
                    MSS_dic[n] = basis_dic[n]
                    push!(qlist, length(Set([n, nR])))
                end
        end    
            index=findall(x -> x==2, qlist)
            new_MSS_dic = Dict(k => v for k in MSS[index] for v in [MSS_dic[k]])
            return MSS[index], new_MSS_dic, qlist
    end
          
end
Fibo_MSS_basis(N::Int, k::Int64, inv::Int64=1) = Fibo_MSS_basis(BitStr{N, Int}, k, inv)

function Fibo_K_Ham(::Type{T}, k::Int, Omega::Float64=1.0) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
#return: the Hamiltonian matrix in given K space

    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basisK, basis_dic = Fibo_K_basis(T, k)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    H = zeros(ComplexF64, (l, l))

    for i in 1:l
        n=basisK[i]
        output = actingHam(T, n, true)
        for m in output
            mbar, d = get_representative(m)
            if mbar ∈ basisK
                j=searchsortedfirst(basisK, mbar)
                Yn= sqrt(length(basis_dic[n])) / N
                Ym= sqrt(length(basis_dic[mbar])) / N
                H[i, j] += Yn/Ym * omegak^d
            end
        end
    end
    if k==0 || k==div(N,2)
        H=real(H)
    end
    H=(H+H')/2
    return H
end
Fibo_K_Ham(N::Int, k::Int, Omega::Float64=1.0) = Fibo_K_Ham(BitStr{N, Int}, k, Omega)

function Fibo_MSS_Ham(::Type{T}, k::Int, inv::Int64=1) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
#return: the Hamiltonian matrix in given maximum symmetry space 
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"

    omegak = exp(2im * π * k / N)
    
    if inv==1 && k==0 || inv==-1 && k==div(N,2)
        MSS, MSS_dic, qlist = Fibo_MSS_basis(T, k, inv)

        l = length(MSS)
        H = zeros(ComplexF64, (l, l))
        for i in 1:l
            n = MSS[i]
            Zn = sqrt(qlist[i]) / 4 * sqrt(length(MSS_dic[n])) / N
            output = actingHam(T, n, true)
            for m in output
                mbar, d = get_representative(m)
                inv_mbar = get_representative(breflect(mbar))[1]
                mtilde = min(mbar, inv_mbar)
                if mtilde ∈ MSS
                    j = searchsortedfirst(MSS, mtilde)
                    Zm = sqrt(qlist[j]) / 4 * sqrt(length(MSS_dic[mtilde])) / N
                    H[i, j] +=  Zn / Zm*omegak^d
                end
            end
        end
        H=real(H)
        H = (H + H') / 2

        return H
    elseif inv==1 && k==div(N,2) || inv==-1 && k==0
        MSS, MSS_dic, _ = Fibo_MSS_basis(T, k, inv)

        l = length(MSS)
        H = zeros(ComplexF64, (l, l))
        for i in 1:l
            n = MSS[i]
            Zn = 1 / 4 * sqrt(length(MSS_dic[n])) / N
            output = actingHam(T, n, true)
            for m in output
                mbar, d = get_representative(m)
                if mbar ∈ MSS
                    j=searchsortedfirst(MSS, mbar)
                    Zm = 1 / 4 * sqrt(length(MSS_dic[mbar])) / N
                    H[i, j] +=  Zn / Zm*omegak^d
                end
            end
        end
        H=real(H)
        H = (H + H') / 2  
        
        return H
    end
    
end
Fibo_MSS_Ham(N::Int, k::Int, inv::Int64=1) = Fibo_MSS_Ham(BitStr{N, Int}, k, inv) 

