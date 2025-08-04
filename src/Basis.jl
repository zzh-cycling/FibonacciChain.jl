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


function Fsymmetry_coef(state::T, base::T, pbc::Bool=true, measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}}
    # Defined as, where idxin idxbond idxout ∈ state, idxbond' ∈ base, in Anyon basis, not in Fibonacci chain basis.
    #  %%%%%%%%%%%% τ, idxin, τ         idxbond
    #  %%
    #  %% 
    #  %%%%%%%%%%%%%
    #  %%
    #  %%
    #  %%%          idxout              idxbond'
    #  or effectively, the coefficient like: A_{x_1 x_2 x_2'}^{x_1'} 
    ϕ = (1+√5)/2
    prod=1
    
    if measure_class == :Fibo
        @assert pbc || site != N "For OBC, site must not be $N, but got $site"
        for site in 1:N-1
            # Identify {x_2}, if x_2' is 1, return 1, otherwise, check x_3', if x_3' is 1, return 1.
            if state[N - site + 1] == 1
                prod *= 1
            else
                if base[N - site] == 1
                    prod *= 1
                else
                    # check x_3, if x_3 is 1
                    if state[N - site] == 0
                        # Determine coef according to x_2'
                        prod *= (base[N - site + 1] == 0) ? -ϕ^(-1) : ϕ^(-1/2)
                    else
                        prod *= (base[N - site + 1] == 0) ? ϕ^(-1/2) : ϕ^(-1)
                    end
                end
            end
        end

        if pbc
            # Check x_N', if x_N' is 1, return 1, otherwise, check x_1', if x_1' is 1, return 1.
            if state[1] == 1
                prod *= 1
            else
                if base[N] == 1
                    prod *= 1
                else
                    # check x_1, if x_1 is 1
                    if state[N] == 0
                        # Determine coef according to x_N'
                        prod *= (base[1] == 0) ? -ϕ^(-1) : ϕ^(-1/2)
                    else
                        prod *= (base[1] == 0) ? ϕ^(-1/2) : ϕ^(-1)
                    end
                end
            
                prod *= 1
            end
        end

        return prod

    else
        error("Unsupported measure_class: $measure_class")
    end
end

function topological_symmetry_basismap(state::T, pbc::Bool=true) where {N, T <: BitStr{N}}
    # Compute the topological symmetry map for a given state using the topological symmetry site map for all site

    basis = Fibonacci_basis(T, pbc, measure_class = :Fibo)
    coeflis = Vector{Float64}(undef, length(basis))
    
    # For each base in basis, check the state at each site
    for (idx, base) in enumerate(basis)
        coef = Fsymmetry_coef(state, base)
        coeflis[idx] = coef
    end
    return coeflis
end

function topological_charge_operator(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    # compute the topological charge operator Yl in the Fibonacci model. default l=0, for tau. l=1, for vacuum.
    
    basis=Fibonacci_basis(T, pbc, measure_class = :Fibo)
    Ymatrix=hcat(topological_symmetry_basismap.(basis)...)

    return Ymatrix
end

function Fibonacci_basis(::Type{T}, pbc::Bool=true; Y=nothing, measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}}
    # Generate basis for Fibonacci model, return BitBasis form, which can be used as binary and decimal form. Here we both consider PBC and OBC
    @assert N > 0 "N is expected to be greater than 0, but got $N"
    @assert Y === nothing || Y in [0, 1, :tau, :trivial] "Y is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $Y"
    @assert T <: BitStr{N} "Type T must be a BitStr type"
    if measure_class == :Fibo
        # Generate Fibonacci chain basis
        # If pbc is true, use Fibonacci_chain_PBC, otherwise use Fibonacci_chain_OBC
        if pbc
            if T == BitStr{1, Int}
                basis=Fibonacci_chain_OBC(T)
            else
                basis=Fibonacci_chain_PBC(T)
            end
        else
            basis=Fibonacci_chain_OBC(T)
        end
        sorted_basis=sort(basis)
    
        if Y !== nothing
            # Filter basis by topological charge
            # sorted_basis = filter(s -> topological_charge(s) == Y, sorted_basis)
            error("Filtering by topological charge is not implemented yet.")
        end
    
        return sorted_basis
    elseif measure_class == :IsingX || measure_class == :IsingZZ
        # Generate basis for Ising model
        return [T(i) for i in 0:(2^N - 1)]
    else
        error("Unsupported measure_class: $measure_class")
    end
end
Fibonacci_basis(N::Int, pbc::Bool=true; Y=nothing, measure_class::Symbol=:Fibo) = Fibonacci_basis(BitStr{N, Int}, pbc; Y=Y, measure_class=measure_class)

function Fibonacci_basis(::Type{T}, k::Int64; Y=nothing, measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system, topological_charge Y, which default to be nothing
#return: computational basis in given momentum kinetically constrained subspace with decimal int form in golden chain model
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert Y === nothing || Y in [0, 1, :tau, :trivial] "Y is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $Y"

    basisK = Vector{T}(undef, 0)
    basis = Fibonacci_basis(T, Y=Y, measure_class=measure_class)
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
Fibonacci_basis(N::Int64, k::Int64; Y=nothing, measure_class::Symbol=:Fibo) = Fibonacci_basis(BitStr{N, Int}, k, Y=Y, measure_class=measure_class)
    
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

function Isingmap(::Type{T}, state::T, i::Int, pbc::Bool=true) where {N, T <: BitStr{N}}
    @assert 1 <= i <= N "i is expected to be in [1, $N], but got $i"
    
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))

    if i == N
        if pbc
             if ((state >> (N - 1)) & 1) == (state & 1)
                return state, X(state,i), -1.0, -1.0 # If same, return -zz and -x
            else
                return state, X(state,i), 1.0, -1.0
            end
        else
            return X(state,i), -1.0
        end
    else
        if ((state >> (N - i - 1)) & 1) == ((state >> (N - i - 2)) & 1)
            return state, X(state,i), -1.0, -1.0 # If same, return -zz and -x
        else
            return state, X(state,i), 1.0, -1.0
        end
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

function actingHam(::Type{T}, state::T, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}}
    # The type of n is DitStr{D, N, Int}, which is a binary string with length N in D-ary form.
    # Acting Hamiltonian on a given state in bitstr and return the output states in bitstr
    # Here need to note that the order of the bitstr is from right to left, which is different from our counting order.
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))
    ϕ = (1+√5)/2

    if measure_class == :Fibo
        mask=bmask(T, N, N-2)
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
            # 1 site 111 fusion, check if is 1xxxx01
            if state & mask1 == mask1
                output[state] = get(output, state, 0.0) - 1
            end
            # N site 111 fusion, check if is 01xxxx1
            if state & mask2 == mask2
                output[state] = get(output, state, 0.0) - 1
            end
        end
        return output
    elseif measure_class == :IsingX || measure_class == :IsingZZ
        # Generate Ising model Hamiltonian
        output = Dict{T, Float64}()
        for i in 1:N-1
            state1, state2, weight1, weight2 = Isingmap(T, state, i, pbc)
            output[state1] = get(output, state1, 0.0) + weight1
            output[state2] = get(output, state2, 0.0) + weight2
        end

        if pbc 
            state1, state2, weight1, weight2 = Isingmap(T, state, N)
            output[state1] = get(output, state1, 0.0) + weight1
            output[state2] = get(output, state2, 0.0) + weight2
        end

        return output
    elseif measure_class == :Ferro
        mask=bmask(T, N, N-2)

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
    else
        error("Unsupported measure_class: $measure_class")
    end
end


function Fibonacci_Ham(::Type{T}, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}}
    # Generate Hamiltonian for Fibonacci model, automotically contain pbc or obc
    basis=Fibonacci_basis(T,pbc, measure_class=measure_class)

    l=length(basis)
    H=zeros(Float64,(l,l))
    for i in 1:l
        output=actingHam(T, basis[i], pbc; measure_class=measure_class) 
        states, weights = keys(output), values(output)
        for m in states
            j=searchsortedfirst(basis, m)
            H[i, j] += output[m]
        end
    end

    return H
end
Fibonacci_Ham(N::Int, pbc::Bool=true; measure_class::Symbol=:Fibo) = Fibonacci_Ham(BitStr{N, Int}, pbc; measure_class=measure_class)


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


function Fibonacci_Ham(::Type{T}, k::Int; Y=nothing, measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system and topological_charge of system
#return: the Hamiltonian matrix in given symmetric sector Hilbert space

    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert Y === nothing || Y in [0, 1, :tau, :trivial] "Y is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $Y"

    basisK, basis_dic =Fibonacci_basis(T, k, Y=Y, measure_class=measure_class)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    H = zeros(ComplexF64, (l, l))

    for i in 1:l
        n=basisK[i]
        output = actingHam(T, n, true)
        states, weights = keys(output), values(output)
        for m in states
            mbar, d = get_representative(m)
            if mbar ∈ basisK
                j=searchsortedfirst(basisK, mbar)
                Yn= sqrt(length(basis_dic[n])) / N
                Ym= sqrt(length(basis_dic[mbar])) / N
                H[i, j] += Yn/Ym * omegak^d*output[m]
            end
        end
    end
    if k==0 || k==div(N,2)
        H=real(H)
    end
    H=(H+H')/2
    return H
end
Fibonacci_Ham(N::Int, k::Int, Y=nothing) = Fibonacci_Ham(BitStr{N, Int}, k, Y)

# join two lists of basis by make a product of two lists, b is placed after a (counts from left to right)
function process_join(a, b)
    # effectively, vec([join(a, b) for a in a for b in b]) # seems faster
    return [join(a, b) for a in a for b in b]
end

# create Fibonacci basis composed of multiple disjoint sub-chains
function joint_basis(lengthlis::Vector{Int}, pbc::Bool=false;measure_class::Symbol=:Fibo)
    # The element order in lengthlis doesn't matter, i.e., [2, 3] is only different with [3, 2] in basis order. Only your input state must consistent with the basis order. 
    if isempty(lengthlis)
        return BitStr{0, Int}[]
    else
        return sort(mapreduce(len -> Fibonacci_basis(len, pbc, measure_class=measure_class), process_join, lengthlis))
    end
end

function connected_components(v::Vector{Int})
    if isempty(v)
        return Int[]
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
    return sum(i -> BitStr{M}(readbit(basis.buf, i) << (M - subsystems[N-i+1])), 1:N)
end

# take environment part of a basis
takeenviron(x, mask::BitStr{l}) where {l} = x & (~mask)
# take system part of a basis
takesystem(x, mask::BitStr{l}) where {l} = (x & mask)

function rdm_Fibo(::Type{T}, subsystems::Vector{Int64}, state::Union{Vector{ET}, Matrix{ET}}, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N,T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.

    unsorted_basis = Fibonacci_basis(T, pbc; measure_class=measure_class)
    subsystems=connected_components(subsystems)
    lengthlis=length.(subsystems)
    subsystems=vcat(subsystems...)
    # mask = bmask(T, subsystems...)
    mask = bmask(T, (N .-subsystems .+1)...)

    
    order = sortperm(unsorted_basis, by = x -> (takeenviron(x, mask), takesystem(x, mask))) #first sort by environment, then by system. The order of environment doesn't matter.

    if state isa Vector
        if isempty(subsystems)
            return ones(ET, 1, 1) # Return empty matrix if no subsystems
        elseif subsystems == collect(1:N)
            # If subsystems is all sites, return the full density matrix
            return state * state' # Full density matrix
        end
        @assert length(unsorted_basis) == length(state) "state length is expected to be $(length(unsorted_basis)), but got $(length(state))"

        basis, state = unsorted_basis[order], state[order]
    elseif state isa Matrix
        if isempty(subsystems)
            return ones(ET, 1, 1) # Return empty matrix if no subsystems
        elseif subsystems == collect(1:N)
            return state
        end
        @assert length(unsorted_basis) == size(state, 2) "state size is expected to be $(length(unsorted_basis)), but got $(size(state, 2))"

        basis, state = unsorted_basis[order], state[order, order]
    end

    
    reduced_basis = move_subsystem.(T, joint_basis(lengthlis, measure_class=measure_class), Ref(subsystems))
    # The reduced_basis counting order doesn't matter, as long as it place subsystem part in basis correctly.
    # TIPS: mask must have common non-zero part with reduced_basis (thus for comparing)
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
        view(reduced_dm, indices, indices) .+= (state isa Vector) ? view(state, range) .* view(state, range)' : view(state, range, range)
    end

    return reduced_dm
end
rdm_Fibo(N::Int, subsystems::Vector{Int64}, state::Union{Vector{ET}, Matrix{ET}}, pbc::Bool=true; measure_class::Symbol=:Fibo) where {ET} = rdm_Fibo(BitStr{N, Int}, subsystems, state, pbc; measure_class=measure_class)

function mapst_sec2tot(::Type{T}, state::Vector{ET}, k::Int64;measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}
    # Map the symmetric sector Hilbert space state to total space state
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basis = Fibonacci_basis(T, measure_class=measure_class)
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
mapst_sec2tot(N::Int, state::Vector{ET}, k::Int64; measure_class::Symbol=:Fibo) where {ET} = mapst_sec2tot(BitStr{N, Int}, state, k, measure_class=measure_class)

function rdm_Fibo_sec(::Type{T}, subsystems::Vector{Int64},kstate::Vector{ET}, k::Int64) where {N,T <: BitStr{N}, ET}
    @assert length(kstate) == length(Fibonacci_basis(T,k)[1]) "state length is expected to be $(length(Fibonacci_basis(T, k)[1])), but got $(length(kstate))"
    state = mapst_sec2tot(T, kstate, k)
    reduced_dm = rdm_Fibo(T, subsystems, state)
    return reduced_dm
end
rdm_Fibo_sec(N::Int, subsystems::Vector{Int64},state::Vector{ET}, k::Int64) where {ET} = rdm_Fibo_sec(BitStr{N, Int}, subsystems, state, k)

# create Fibonacci basis composed of multiple disjoint chains with different basis type
function joint_basis(lengthlisA::Vector{Int}, lengthlisB::Vector{Int};subApbc::Bool=false, subBpbc::Bool=false, measure_classA::Symbol=:Fibo, measure_classB::Symbol=:Fibo)
    # subpbc is used to indicate whether the subsystem is periodic or not
    lisA = sort(mapreduce(len -> Fibonacci_basis(len, subApbc, measure_class=measure_classA), process_join, lengthlisA))
    lisB = sort(mapreduce(len -> Fibonacci_basis(len, subBpbc, measure_class=measure_classB), process_join, lengthlisB))
    return vec([join(i, j) for i in lisA for j in lisB])
end

function disjoint_rdm(::Type{T1}, ::Type{T2}, subsystemsA::Vector{Int64}, subsystemsB::Vector{Int64}, state::Vector{ET}, pbc::Bool=true; totalsubApbc::Bool=false, totalsubBpbc::Bool=false, measure_classA::Symbol=:Fibo, measure_classB::Symbol=:Fibo) where {N1, N2,T1 <: BitStr{N1},T2 <: BitStr{N2}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the two disjoint basis (two part in one chain), espeically, can be viewed as two parrelel chain. Given the system size, subsystem and basis type, to get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    # pbc is the basis boundary condition, while totalsubpbc is used to indicate whether the total subsystem is periodic or not
    @assert all(subsystemsB .<= N2) "subsystemsB is expected to be in [1, $N2], but got $(subsystemsB)"
    if isempty(subsystemsA) && isempty(subsystemsB)
        return ones(ET, 1, 1) # Return empty matrix if no subsystems
    elseif subsystemsA == collect(1:N1) && subsystemsB == collect(1:N2)
        return state * state'
    end

    unsorted_basisA = Fibonacci_basis(T1, pbc, measure_class=measure_classA)
    unsorted_basisB = Fibonacci_basis(T2, pbc, measure_class=measure_classB)
    lenubasisA = length(unsorted_basisA)
    lenubasisB = length(unsorted_basisB)
    newT = BitStr{N1+N2, Int} # double the length of the basis
    doublebasis = reshape([join(i,j) for i in unsorted_basisA for j in unsorted_basisB], lenubasisA*lenubasisB)
    # Align as [A, B] basis
    @assert lenubasisA*lenubasisB == length(state) "state length is expected to be $(lenubasisA*lenubasisB), but got $(length(state))"

    mask = bmask(newT, (N1 .-subsystemsA .+1).+N2..., (N2 .-subsystemsB .+1)...)
    
    
    subsystems = vcat(subsystemsA, subsystemsB .+N1)
    subsystems = connected_components(subsystems)
    lengthlis = length.(subsystems)
    subsystems = vcat(subsystems...)
    
    subsystemsB = subsystemsB.+N1 # add the second half of the system to the subsystems
    subsystemsA = connected_components(subsystemsA)
    subsystemsB = connected_components(subsystemsB)
    lengthlisA = length.(subsystemsA)
    lengthlisB = length.(subsystemsB)
    subsystemsA = vcat(subsystemsA...)
    subsystemsB = vcat(subsystemsB...)

    if isempty(subsystemsA)
        # If subsystemsA is empty, we only consider subsystemsB, thus parameter is all about B, totalsubBpbc, measure_classB
        reduced_basis = move_subsystem.(newT, joint_basis(lengthlis, totalsubBpbc, measure_class = measure_classB), Ref(subsystems))
    elseif isempty(subsystemsB)
        reduced_basis = move_subsystem.(newT, joint_basis(lengthlis, totalsubApbc, measure_class = measure_classA), Ref(subsystems))
    else
        reduced_basis = move_subsystem.(newT, joint_basis(lengthlisA, lengthlisB, subApbc=totalsubApbc, subBpbc=totalsubBpbc, measure_classA = measure_classA, measure_classB = measure_classB), Ref(vcat(subsystemsA, subsystemsB)))
    end

    order = sortperm(doublebasis, by = x -> (takeenviron(x, mask), takesystem(x, mask))) #first sort by environment, then by system. The order of environment doesn't matter. Taking order starts from the left.
    basis, state = doublebasis[order], state[order]
    
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
disjoint_rdm(N1::Int64, N2::Int64, subsystemsA::Vector{Int64}, subsystemsB::Vector{Int64}, state::Vector{ET}, pbc::Bool=true, totalsubApbc::Bool=false, totalsubBpbc::Bool=false, measure_classA::Symbol=:Fibo, measure_classB::Symbol=:Fibo) where {ET} = 
disjoint_rdm(BitStr{N1, Int}, BitStr{N2, Int}, subsystemsA, subsystemsB, state, pbc, 
totalsubApbc=totalsubApbc, totalsubBpbc=totalsubBpbc,
measure_classA=measure_classA, measure_classB=measure_classB)

function reference_rdm(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}; pbc::Bool=true, measure_class::Symbol=:Fibo, traceref::Bool=true) where {N, T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    # N is the particle number of system, while k_old is the number of reference qubit, which is deduced from the state length.
    unsorted_basis = Fibonacci_basis(T, pbc; measure_class=measure_class)
    len_F   = length(unsorted_basis)
    k_old = round(Int, log2(length(state) ÷ len_F))
    
    length(state) == (2^k_old * len_F) || error("state length is not compatible with (k_old, N), can not deduce k_old from state length")
    @assert 2^k_old*length(unsorted_basis) == length(state) "state length is expected to be $(2^k_old*length(unsorted_basis)), but got $(length(state))"
    if traceref
        # If traceref is true, we need to trace out the reference qubit. otherwise, we trace out system.
        return disjoint_rdm(BitStr{k_old, Int}, T, subsystems, Int[], state, pbc; measure_classA=:IsingX, measure_classB=measure_class)
    else
        totalsubBpbc = (length(subsystems) == N) ? true : false
        return disjoint_rdm(BitStr{k_old, Int}, T, Int[], subsystems, state, pbc; measure_classA=:IsingX, totalsubBpbc=totalsubBpbc, measure_classB=measure_class)
    end
   
  
end
reference_rdm(N::Int, subsystems::Vector{Int}, state::Vector{ET}; pbc::Bool=true, measure_class::Symbol=:Fibo, traceref::Bool=false) where {ET} = reference_rdm(BitStr{N, Int}, subsystems, state, pbc=pbc, measure_class=measure_class, traceref=traceref)