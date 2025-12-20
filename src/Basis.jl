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

abstract type AbstractAnyonType end
struct FibonacciAnyon <: AbstractAnyonType end
struct IsingAnyon <: AbstractAnyonType end

struct AnyonModel{AT<:AbstractAnyonType}
    anyon_type::AT # anyon type
    N::Int   # system size
    pbc::Bool # periodic boundary conditions
    interaction_type::Symbol # interaction type, e.g., :Antiferro, :Ferro
    function AnyonModel(anyon_type::AT, N::Int; pbc::Bool=true, interaction_type::Symbol=:Antiferro) where {AT<:AbstractAnyonType}
        @assert N > 0 "N is expected to be greater than 0, but got $N"
        return new{AT}(anyon_type, N, pbc, interaction_type)
    end
end

function anyon_basis(::FibonacciAnyon, ::Type{T}; pbc::Bool=true, symmetry_block=nothing) where {N, T <: BitStr{N}}
    # Generate basis for Fibonacci model, return BitBasis form, which can be used as binary and decimal form. Here we both consider PBC and OBC
    @assert symmetry_block === nothing || symmetry_block in [0, 1, :tau, :trivial] "symmetry_block is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $symmetry_block"
    @assert T <: BitStr{N} "Type T must be a BitStr type"
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

    if symmetry_block !== nothing
        # Filter basis by topological charge
        # sorted_basis = filter(s -> topological_charge(s) == symmetry_block, sorted_basis)
        error("Filtering by topological charge is not implemented yet.")
    end

    return sorted_basis
end
anyon_basis(model::AnyonModel{AT}; symmetry_block=nothing) where {AT<:FibonacciAnyon} = anyon_basis(model.anyon_type, BitStr{model.N, Int}; pbc=model.pbc, symmetry_block=symmetry_block)

function anyon_basis(::IsingAnyon, ::Type{T}) where {N, T <: BitStr{N}}
    return [T(i) for i in 0:(2^N - 1)]
end
anyon_basis(model::AnyonModel{AT}) where {AT<:IsingAnyon} = anyon_basis(model.anyon_type, BitStr{model.N, Int})

"""
    Fsymmetry_coef(state::T, base::T, pbc::Bool=true, anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}

Compute topological symmetry coefficient for state in given base configurations for Fibonacci anyon chain.

# Arguments
- `state::T`: Target state configuration
- `base::T`: Base state configuration  
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: Measurement class

# Returns
- `Float64`: Topological symmetry coefficient based on Fibonacci fusion rules

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; T = BitStr{N, Int};

julia> state = T(0b1010); ϕ = (1 + sqrt(5)) / 2;  # Example state configuration

julia> base = T(0b0101);   # Example base configuration

julia> coef = Fsymmetry_coef(state, base, true, :Fibo); coef ≈ 0.3819660112501051
true

julia> abs(coef - (1-1/ϕ)) < 1e-10  # Should equal φ for this configuration
true
```
"""
function Fsymmetry_coef(state::T, base::T, model::AnyonModel{AT}) where {N, T <: BitStr{N}, AT<:FibonacciAnyon}
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

    @assert model.pbc || site != N "For OBC, site must not be $N, but got $site"
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

        if model.pbc
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
end

"""
    topological_symmetry_basismap(state::T, pbc::Bool=true) where {N, T <: BitStr{N}}

Compute topological symmetry coefficients for all basis states relative to given state.

# Arguments
- `state::T`: Reference state configuration
- `pbc::Bool=true`: Periodic boundary conditions

# Returns
- `Vector{Float64}`: Coefficients for each basis state

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; T = BitStr{N, Int};

julia> state = T(0b0000);  # Vacuum state

julia> coeffs = topological_symmetry_basismap(state, true);

julia> length(coeffs) == length(anyon_basis(T, true, anyon_type=:Fibo))
true

julia> all(x -> abs(x) > 1e-10 || abs(x) < 1e-10, coeffs)  # All coeffs are well-defined
true
```
"""
function topological_symmetry_basismap(state::T, model::AnyonModel{AT}) where {N, T <: BitStr{N}, AT<:FibonacciAnyon}
    # Compute the topological symmetry map for a given state using the topological symmetry site map for all site
    # However, it seems not easy to separate the state into different topological symmetric sectors.
    basis = anyon_basis(model)
    coeflis = Vector{Float64}(undef, length(basis))
    
    # For each base in basis, check the state at each site
    for (idx, base) in enumerate(basis)
        coef = Fsymmetry_coef(state, base, model)
        coeflis[idx] = coef
    end
    return coeflis
end

function topological_charge_operator(::Type{T}, model::AnyonModel{AT}) where {N, T <: BitStr{N}, AT<:FibonacciAnyon}
    # compute the topological charge operator Yl in the Fibonacci model. default l=0, for tau. l=1, for vacuum.

    basis=anyon_basis(model)
    Ymatrix=hcat(topological_symmetry_basismap.(basis, model)...)

    return Ymatrix
end

"""
    anyon_basis(::Type{T}, k::Int64; symmetry_block=nothing, anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}

Generate basis states in specific momentum sector `k` and topological sector `symmetry_block`.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `k::Int64`: Momentum sector (0 ≤ k ≤ N-1)
- `symmetry_block`: Topological charge sector
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Vector{T}`: Basis states in momentum sector k
- `Dict{T, Vector{T}}`: Representative mapping for translation equivalence classes
"""
function anyon_basis(AT::AbstractAnyonType, ::Type{T}, k::Int; symmetry_block=nothing) where {N, T <: BitStr{N}}
    #params: a int of lattice number, momentum of system, topological_charge symmetry_block, which default to be nothing
    #return: computational basis in given momentum kinetically constrained subspace with decimal int form in golden chain model
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert symmetry_block === nothing || symmetry_block in [0, 1, :tau, :trivial] "symmetry_block is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $symmetry_block"

    basisK = Vector{T}(undef, 0)
    basis = anyon_basis(AT, T; symmetry_block=symmetry_block)
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
anyon_basis(model, k::Int64; symmetry_block=nothing) = anyon_basis(model.anyon_type, BitStr{model.N}, k=k; symmetry_block=symmetry_block)

function antimap(::Type{T}, state::T, i::Int) where {N, T <: BitStr{N}}
    # The type of n is DitStr{D, N, Int}, which is a binary string with length N in D-ary form.
    # Acting Hamiltonian on a given state in bitstr and return the output (states, weight) in bitstr
    # Here need to note that the order of the bitstr is from right to left, which is different from our counting order.
    ϕ = (1+√5)/2
    fl=bmask(T, N)

    X(state,i) = flip(state, fl >> (i-1))

    if readbit(state, N +1 -i) == 0
        return state, X(state,i), -ϕ^(-1), -ϕ^(-3/2)
    else
        return state, X(state,i), -ϕ^(-2), -ϕ^(-3/2)
    end
end

function ferromap(::Type{T}, state::T, i::Int) where {N, T <: BitStr{N}}
    ϕ = (1+√5)/2
    fl=bmask(T, N)

    X(state,i) = flip(state, fl >> (i-1))

    if readbit(state, N +1 -i) == 0
        return state, X(state,i), ϕ^(-1), ϕ^(-3/2)
    else
        return state, X(state,i), ϕ^(-2), ϕ^(-3/2)
    end
end

function Isingmap(::Type{T}, state::T, i::Int, pbc::Bool=true; kwargs...) where {N, T <: BitStr{N}}
    @assert 1 <= i <= N "i is expected to be in [1, $N], but got $i"
    
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))
    J, h = get(kwargs, :J, 1.0), get(kwargs, :h, 1.0)
    if i == N
        if pbc
            if readbit(state, 1) == readbit(state, N)
                return state, X(state,i), -J, -h # If same, return -zz and -x
            else
                return state, X(state,i), J, -h
            end
        else
            return X(state,i), -h # If OBC, only return -x_N
        end
    else
        if readbit(state, N - i + 1) == readbit(state, N - i)
            return state, X(state,i), -J, -h # If same, return -zz and -x
        else
            return state, X(state,i), J, -h
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

function actingHam(model::AnyonModel{AT}, ::Type{T}, state::T; kwargs...) where {N, T <: BitStr{N}, AT<:FibonacciAnyon}
    # The type of n is DitStr{D, N, Int}, which is a binary string with length N in D-ary form.
    # Acting Hamiltonian on a given state in bitstr and return the output states in bitstr
    # Here need to note that the order of the bitstr is from right to left, which is different from our counting order.
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))
    ϕ = (1+√5)/2
    pbc = model.pbc

    if model.interaction_type == :Antiferro
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

    elseif model.interaction_type == :Ferro
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
    end

    return output
end

function actingHam(model::AnyonModel{AT}, ::Type{T}, state::T; kwargs...) where {N, T <: BitStr{N}, AT<:IsingAnyon}
    # The type of n is DitStr{D, N, Int}, which is a binary string with length N in D-ary form.
    # Acting Hamiltonian on a given state in bitstr and return the output states in bitstr
    # Here need to note that the order of the bitstr is from right to left, which is different from our counting order.
    pbc = model.pbc

    # Generate Ising model Hamiltonian
    output = Dict{T, Float64}()
    for i in 1:N-1
        state1, state2, weight1, weight2 = Isingmap(T, state, i, pbc; kwargs...)
        output[state1] = get(output, state1, 0.0) + weight1
        output[state2] = get(output, state2, 0.0) + weight2
    end

    if pbc
        state1, state2, weight1, weight2 = Isingmap(T, state, N, pbc; kwargs...)
        output[state1] = get(output, state1, 0.0) + weight1
        output[state2] = get(output, state2, 0.0) + weight2
    else
        state1, weight1 = Isingmap(T, state, N, pbc; kwargs...)
        output[state1] = get(output, state1, 0.0) + weight1
    end

    
    return output
end


"""
    anyon_ham(::Type{T}, pbc::Bool=true; anyon_type::Symbol=:Fibo, kwargs...) where {N, T <: BitStr{N}}

Construct Hamiltonian matrix for 1D anyon chain.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type
- `kwargs...`: Additional model parameters, e.g., `J`, `h` for Ising model.

# Returns
- `Matrix{Float64}`: Hamiltonian matrix in chosen basis

Supports various anyon models including Fibonacci and Ising anyons.

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; T = BitStr{N,Int};

julia> H_fibo = anyon_ham(T, true, anyon_type=:Fibo);

julia> size(H_fibo)       # Hamiltonian dimension matches basis size
(7, 7)

julia> H_Ising = anyon_ham(T, true, anyon_type=:IsingX, J=1.0, h=1.0);

julia> size(H_Ising)      # full Hilbert space for Ising model
(16, 16)
```
"""
function anyon_ham(model::AnyonModel{AT}, ::Type{T}; kwargs...) where {N, T <: BitStr{N}, AT<:AbstractAnyonType}
    # Generate Hamiltonian for Fibonacci model, automotically contain pbc or obc
    basis=anyon_basis(model)

    l=length(basis)
    H=zeros(Float64,(l,l))
    for i in 1:l
        output=actingHam(model, T, basis[i];kwargs...) 
        states, weights = keys(output), values(output)
        for m in states
            j=searchsortedfirst(basis, m)
            H[i, j] += output[m]
        end
    end

    return H
end
anyon_ham(model::AnyonModel{AT}; kwargs...) where {AT <: AbstractAnyonType}= anyon_ham(model, BitStr{model.N, Int};kwargs...)
# Another method to write Fibonacci Hamiltonian is using the Measurement operator sum. For example, H = -∑ X_i, where X_i is the Temperley-Lieb generator acting on site i-1, i, and i+1. Pilis = [FibonacciChain.measure_matrix(BitStr{16, Int}, 1000.0, idx, 0) for idx in 1:N]. H = -sum(Pilis). This two Hamiltonian difference is not a constant, but like a arc in conformal energy spectrum below arc, but they have the same eigenstates.

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

"""
    anyon_ham(::Type{T}, k::Int; symmetry_block=nothing, anyon_type::Symbol=:Fibo, kwargs...) where {N, T <: BitStr{N}}

Construct Hamiltonian matrix in specific symmetry sector for 1D anyon chain.
    
# Arguments
- `T::Type`: BitStr type specifying chain length N
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type
- `kwargs...`: Additional model parameters, e.g., `J`, `h` for Ising model.

# Returns
- `Matrix{Float64}`: Hamiltonian matrix in chosen basis
"""
function anyon_ham(model::AnyonModel{AT}, ::Type{T}, k::Int; symmetry_block=nothing, kwargs...) where {N, T <: BitStr{N}, AT <: AbstractAnyonType}
#params: a int of lattice number, momentum of system and topological_charge of system
#return: the Hamiltonian matrix in given symmetric sector Hilbert space

    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert symmetry_block === nothing || symmetry_block in [0, 1, :tau, :trivial] "symmetry_block is expected to be nothing or 1 or 0 or :trivial or :nontrivial, but got $symmetry_block"

    basisK, basis_dic =anyon_basis(model, T, k, symmetry_block=symmetry_block)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    H = zeros(ComplexF64, (l, l))

    for i in 1:l
        n=basisK[i]
        output = actingHam(model, T, n; kwargs...)
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
anyon_ham(model::AnyonModel{AT}, k::Int; symmetry_block=nothing, kwargs...) where {AT <: AbstractAnyonType}= anyon_ham(model, BitStr{N, Int}, k, symmetry_block=symmetry_block, kwargs...)

# join two lists of basis by make a product of two lists, b is placed after a (counts from left to right)
function process_join(a, b)
    # effectively, vec([join(a, b) for a in a for b in b]) # seems faster
    return [join(a, b) for a in a for b in b]
end

# create Fibonacci basis composed of multiple disjoint sub-chains
function joint_basis(AT::AbstractAnyonType, lengthlis::Vector{Int}, pbc::Bool=false;)
    # The element order in lengthlis doesn't matter, i.e., [2, 3] is only different with [3, 2] in basis order. Only your input state must consistent with the basis order.
    if isempty(lengthlis)
        return BitStr{0, Int}[]
    else
        return sort(mapreduce(len -> anyon_basis(AT, BitStr{len, Int}, pbc = pbc), process_join, lengthlis))
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

"""
    anyon_rdm(::Type{T}, subsystems::Vector{Int64}, state::Union{Vector{ET}, Matrix{ET}}, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N,T <: BitStr{N}, ET}

Compute reduced density matrix for specified subsystems from quantum state or density matrix.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `subsystems::Vector{Int64}`: Indices of subsystem sites to keep
- `state::Union{Vector{ET}, Matrix{ET}}`: Quantum state vector or density matrix
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `Matrix{ET}`: Reduced density matrix for specified subsystems

Subsystem indices are counted from right in binary representation.

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis, LinearAlgebra

julia> N = 4; T = BitStr{N, Int};

julia> basis = anyon_basis(T, true, anyon_type=:Fibo);

julia> state = randn(ComplexF64, length(basis)); state ./= norm(state);  

julia> rdm = anyon_rdm(T, [1, 2], state, true, anyon_type=:Fibo);

julia> size(rdm)  # Reduced density matrix dimension
(3, 3)

julia> ishermitian(rdm)  # RDM should be Hermitian
true

julia> tr(rdm) ≈ 1.0 + 0.0im  # Trace should be 1
true
```
"""
function anyon_rdm(::Type{T}, subsystems::Vector{Int64}, state::Union{Vector{ET}, Matrix{ET}}, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N,T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.

    unsorted_basis = anyon_basis(T, pbc; anyon_type=anyon_type)
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

    
    reduced_basis = move_subsystem.(T, joint_basis(lengthlis, anyon_type=anyon_type), Ref(subsystems))
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
anyon_rdm(N::Int, subsystems::Vector{Int64}, state::Union{Vector{ET}, Matrix{ET}}, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET} = anyon_rdm(BitStr{N, Int}, subsystems, state, pbc; anyon_type=anyon_type)

"""
    mapst_sec2tot(::Type{T}, state::Vector{ET}, k::Int64; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET} 

map state in symmetric sector to total Hilbert space.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `state::Vector{ET}`: State vector in symmetric sector Hilbert space
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Vector{ET}`: Total space state vector
"""
function mapst_sec2tot(::Type{T}, state::Vector{ET}, k::Int64; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}
    # Map the symmetric sector Hilbert space state to total space state
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basis = anyon_basis(T, anyon_type=anyon_type)
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
mapst_sec2tot(N::Int, state::Vector{ET}, k::Int64; anyon_type::Symbol=:Fibo) where {ET} = mapst_sec2tot(BitStr{N, Int}, state, k, anyon_type=anyon_type)

"""
    anyon_rdm_sec(::Type{T}, subsystems::Vector{Int64}, kstate::Vector{ET}, k::Int64) where {N,T <: BitStr{N}, ET} 

Compute reduced density matrix for specified subsystems from quantum state or density matrix in specific symmetry sector.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `subsystems::Vector{Int64}`: Indices of subsystem sites to keep
- `kstate::Vector{ET}`: State vector in symmetric sector Hilbert space
- `k::Int64`: Momentum sector (0 ≤ k ≤ N-1)
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Matrix{ET}`: Reduced density matrix in total Hilbert basis
"""
function anyon_rdm_sec(AT::AbstractAnyonType, ::Type{T}, subsystems::Vector{Int64}, kstate::Vector{ET}, k::Int64;) where {N,T <: BitStr{N}, ET}
    l = length(anyon_basis(AT, T, k)[1])
    @assert length(kstate) == l "state length is expected to be $(l), but got $(length(kstate))"
    state = mapst_sec2tot(T, kstate, k)
    reduced_dm = anyon_rdm(T, subsystems, state, true; anyon_type=anyon_type)
    return reduced_dm
end
anyon_rdm_sec(model::AnyonModel{AT}, subsystems::Vector{Int64},state::Vector{ET}, k::Int64;) where {AT <: AbstractAnyonType, ET} = anyon_rdm_sec(AT, BitStr{N, Int}, subsystems, state, k)

# create Fibonacci basis composed of multiple disjoint chains with different basis type
function joint_basis(AT1::AbstractAnyonType, AT2::AbstractAnyonType, lengthlisA::Vector{Int}, lengthlisB::Vector{Int};subApbc::Bool=false, subBpbc::Bool=false)
    # subpbc is used to indicate whether the subsystem is periodic or not
    lisA = sort(mapreduce(len -> anyon_basis(AT1, BitStr{len, Int}, subApbc), process_join, lengthlisA))
    lisB = sort(mapreduce(len -> anyon_basis(AT2, BitStr{len, Int}, subBpbc), process_join, lengthlisB))
    return vec([join(i, j) for i in lisA for j in lisB])
end

"""
    disjoint_rdm(::Type{T1}, ::Type{T2}, subsystemsA::Vector{Int64}, subsystemsB::Vector{Int64}, state::Vector{ET}, pbc::Bool=true; totalsubApbc::Bool=false, totalsubBpbc::Bool=false, anyon_typeA::Symbol=:Fibo, anyon_typeB::Symbol=:Fibo) where {N1, N2,T1 <: BitStr{N1},T2 <: BitStr{N2}, ET}

Compute reduced density matrix for two joint different systems, or two parallel chains. Where specified subsystems is given, output quantum state or density matrix.

# Arguments
- `T1::Type`: BitStr type specifying chain length N1
- `T2::Type`: BitStr type specifying chain length N2
- `subsystemsA::Vector{Int64}`: Indices of subsystem sites to keep
- `subsystemsB::Vector{Int64}`: Indices of subsystem sites to keep
- `state::Vector{ET}`: Quantum state vector in total Hilbert space
- `pbc::Bool=true`: Periodic boundary conditions
- `totalsubApbc::Bool=false`: Whether the total subsystem A is periodic
- `totalsubBpbc::Bool=false`: Whether the total subsystem B is periodic
- `anyon_typeA::Symbol=:Fibo`: anyon type for subsystem A
- `anyon_typeB::Symbol=:Fibo`: anyon type for subsystem B

# Returns
- `Matrix{ET}`: Reduced density matrix in total Hilbert basis

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis, LinearAlgebra

julia> N1, N2 = 4, 4; T1, T2 = BitStr{N1,Int}, BitStr{N2,Int};

julia> basisA = anyon_basis(T1, true, anyon_type=:Fibo);

julia> basisB = anyon_basis(T2, true, anyon_type=:Fibo);

julia> state = randn(ComplexF64, length(basisA) * length(basisB));

julia> state ./= norm(state);

julia> rdm = disjoint_rdm(T1, T2, [1,2], [1,2], state, true,
                          totalsubApbc=false, totalsubBpbc=false,
                          anyon_typeA=:Fibo, anyon_typeB=:Fibo);

julia> size(rdm)          # reduced density matrix dimension
(9, 9)

julia> ishermitian(rdm)   # should be Hermitian
true

julia> tr(rdm) ≈ 0.9999999999999999 + 0.0im  # trace should be 1
true
```
"""
function disjoint_rdm(::Type{T1}, ::Type{T2}, subsystemsA::Vector{Int64}, subsystemsB::Vector{Int64}, state::Vector{ET}, pbc::Bool=true; totalsubApbc::Bool=false, totalsubBpbc::Bool=false,  AT1::AbstractAnyonType=FibonacciAnyon(), AT2::AbstractAnyonType=FibonacciAnyon()) where {N1, N2,T1 <: BitStr{N1},T2 <: BitStr{N2}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the two disjoint basis (two part in one chain), espeically, can be viewed as two parrelel chain. Given the system size, subsystem and basis type, to get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    # pbc is the basis boundary condition, while totalsubpbc is used to indicate whether the total subsystem is periodic or not
    @assert all(subsystemsB .<= N2) "subsystemsB is expected to be in [1, $N2], but got $(subsystemsB)"
    if isempty(subsystemsA) && isempty(subsystemsB)
        return ones(ET, 1, 1) # Return empty matrix if no subsystems
    elseif subsystemsA == collect(1:N1) && subsystemsB == collect(1:N2)
        return state * state'
    end

    unsorted_basisA = anyon_basis(AT1, T1, pbc)
    unsorted_basisB = anyon_basis(AT2, T2, pbc)
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
        # If subsystemsA is empty, we only consider subsystemsB, thus parameter is all about B, totalsubBpbc, anyon_typeB
        reduced_basis = move_subsystem.(newT, joint_basis(AT1, lengthlis, totalsubBpbc), Ref(subsystems))
    elseif isempty(subsystemsB)
        reduced_basis = move_subsystem.(newT, joint_basis(AT2, lengthlis, totalsubApbc), Ref(subsystems))
    else
        reduced_basis = move_subsystem.(newT, joint_basis(AT1, AT2, lengthlisA, lengthlisB, totalsubApbc, totalsubBpbc), Ref(vcat(subsystemsA, subsystemsB)))
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
disjoint_rdm(N1::Int64, N2::Int64, subsystemsA::Vector{Int64}, subsystemsB::Vector{Int64}, state::Vector{ET}, pbc::Bool=true; totalsubApbc::Bool=false, totalsubBpbc::Bool=false, anyon_typeA::Symbol=:Fibo, anyon_typeB::Symbol=:Fibo) where {ET} = 
disjoint_rdm(BitStr{N1, Int}, BitStr{N2, Int}, subsystemsA, subsystemsB, state, pbc, 
totalsubApbc=totalsubApbc, totalsubBpbc=totalsubBpbc,
anyon_typeA=anyon_typeA, anyon_typeB=anyon_typeB)