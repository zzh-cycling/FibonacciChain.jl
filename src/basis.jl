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