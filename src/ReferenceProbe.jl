function build_extended_basis(k_total::Int, basis::Vector{ET}) where {ET}
    T = BitStr{k_total, Int}
    ref_strings = [T(i) for i in 0:(2^k_total-1)] 
    extended_basis = sort(
        mapreduce(suffix -> process_join([suffix], basis),
              vcat,
              ref_strings))

    return extended_basis
end
# process_join([suffix], basis) will give [0000, 0001, 0010, 0011], but process_join(basis, [suffix]) will give [0000, 0001, 0010, 0100]

function reference_measure_basismap(::Type{T}, τ::Float64, state::ET, i::Int, sign::Int64, pbc::Bool=true; k_old::Int64=1, measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}
    # default for PBC system, map basis
    @assert k_old >= 0 "k_old must be at least 0, but got $(k_old)"
 
    mask = bmask(BitStr{N+k_old, Int}, 1:N...)
    action_state = T(takesystem(state, mask))
    return measure_basismap(T, τ, action_state, i, sign, pbc, measure_class=measure_class)
   
end

function reference_measuremap(::Type{T}, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true;k_old::Int64=1,measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}
    # input a superposition state with reference qubit, and output the measured state. k_old is the number of reference qubits in the state.
    if measure_class == :Fibo
        @assert pbc || (2 <= idx <= N-1) "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif measure_class == :IsingZZ
        @assert pbc || (1 <= idx <= N-1) "Index idx must be in [1, N-1] for open BC (IsingZZ)"
    elseif measure_class ∈ (:IsingX, :IsingZ, :reset, :resetFibo)
        @assert pbc || (1 <= idx <= N) "Index idx must be in [1, N] for open BC (IsingX)"
    else
        error("Unknown measure class: $measure_class")
    end
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"
    # @assert k_old >= 1 "k_old must be at least 1, but got $(k_old)" # because join(bit"1", outputstate2), in contrast to reference_measure_basismap

    basis=Fibonacci_basis(T, pbc, measure_class=measure_class)
    
    # Noting that basis is not consisten with state, but extended_basis is.
    extended_basis = build_extended_basis(k_old, basis) 
    @assert length(basis)*2^k_old == length(state) == length(extended_basis) "state length is expected to be $(length(basis)*2^k_old), but got $(length(state))"
    l=length(extended_basis)

    mapped_state = zeros(ET, length(state))

    pretype = BitStr{k_old, Int}
    mask = bmask(BitStr{N+k_old, Int}, 1:N...)
    
    for (i, ext_basis_i) in enumerate(extended_basis)
        output = reference_measure_basismap(T, τ, ext_basis_i, idx, sign, pbc, k_old=k_old, measure_class=measure_class)
        
        prefix_i = pretype(takeenviron(ext_basis_i, mask) >> N)
        
        if length(output) == 4
            outputstate1, outputstate2, output1, output2=output
            j2=searchsortedfirst(extended_basis, join(prefix_i, outputstate2))
            mapped_state[i]+=output1*state[i] # outputstate1 is the same as basis[i]
            mapped_state[j2]+=output2*state[i]
        else
            outputstate, output1=output # outputstate is the same as basis[i]
            mapped_state[i]+=output1*state[i]
        end
    end

    return mapped_state
end
reference_measuremap(N::Int, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true; k_old::Int64=1, measure_class::Symbol=:Fibo) where {ET} = reference_measuremap(BitStr{N, Int}, τ, state, idx, sign, pbc, k_old=k_old, measure_class=measure_class)

function reference_apply_measurement_layer!(N::Int64, state::Vector{T}, τ::Float64, layer_sample::Vector{Int64}, layer_idx::Int64, pbc::Bool=true; k_old::Int64=1, measure_class::Symbol=:Fibo) where {T}
    if measure_class == :Fibo
        if layer_idx % 2 == 1
            measurement_sites = collect(2:2:N)  # odd sites anyons, even sites qubits
        else
            measurement_sites = collect(1:2:N)  # even sites anyons, odd sites qubits
        end
        for (idx, measurement_type) in enumerate(layer_sample)
            state = reference_measuremap(N, τ, state, measurement_sites[idx], measurement_type, pbc, k_old=k_old, measure_class = measure_class)
            normalize!(state)
        end
        return state
    
    elseif measure_class == :IsingX || measure_class == :IsingZZ
        measurement_sites = collect(1:N)
        if layer_idx % 2 == 1
            # odd layers: measure X
            for (idx, measurement_type) in enumerate(layer_sample)
                state = reference_measuremap(N, τ, state, measurement_sites[idx], measurement_type, pbc, k_old=k_old, measure_class = :IsingX)
                normalize!(state)
            end
            return state
        else
            # even layers: measure ZZ
            for (idx, measurement_type) in enumerate(layer_sample)
                state = reference_measuremap(N, τ, state, measurement_sites[idx], measurement_type, pbc, k_old= k_old, measure_class = :IsingZZ)
                normalize!(state)
            end
            return state
        end
    else
        error("Unknown measure class: $measure_class")
    end
end

function reference_generate_state(τ::Float64, state::Vector{T}, sample::ET, pbc::Bool=true; temp::Bool=true, k_old::Int64=1, measure_class::Symbol=:Fibo) where{T, ET}
    @assert ET == Matrix{Int} "ET must be Matrix{Int} for reference_generate_state"

    D, N = size(sample, 1), 2 * size(sample, 2) 
    statelis = temp ? Vector{Vector{T}}(undef, D) : nothing
 
    for layer in 1:D
        τ_eff = (layer == D) ? τ/2 : τ
        state = reference_apply_measurement_layer!(N, state, τ_eff, sample[layer, :], layer, pbc, k_old=k_old, measure_class = measure_class)

        if temp
            statelis[layer] = copy(state)
        end
    end
    
    return temp ? statelis : state
end

function add_reference_qubits!(N::Int, state::Vector{ET}, site_idx::Int64, rng::MersenneTwister=MersenneTwister(); k_new::Int=1, pbc::Bool=true, measure_class::Symbol=:Fibo) where {ET}
    # Add k_new reference qubits to the state at the specified site_idx, and place them to the left part of basis (index N-site_idx+1) to form a maximally entangled state.
    @assert 1 <= site_idx <= N "Site index must be in the range [1, N]"
    1 >= k_new >= 0 || error("k_new must be in [0,1]")
    # Because each qubit can only concat with one reference qubit, so k_new can only be 0 or 1. If need to add more reference qubits, use add_reference_qubits! multiple times at different site.
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
    extended_basis_old = build_extended_basis(k_old, basis_F)
    T = BitStr{N + k_old, Int}
    fl=bmask(T, N)
    mask = bmask(BitStr{N+k_total, Int}, 1:N...)
    X(state,i) = flip(state, fl >> (i-1))
    
    offset = length(state)   # new k_new qubit starting position
    # NEED to do RESET qubit to 0!!! then concat with the new reference qubit. 
    resettype = (measure_class ∈ (:Fibo,)) ? :resetFibo : :reset
    state_after_0 = reference_measuremap(N, 1000.0, state, site_idx, 0, pbc,k_old=k_old, measure_class = resettype)
    prob_sqrt0 = state_after_0' * state_after_0
    prob_sqrt1 = 1 - prob_sqrt0
    random_number = rand(rng)  
    @show prob_sqrt0, random_number
    # random_number = 0.1
    
    if random_number < prob_sqrt0
        current_state = state_after_0 ./ sqrt(prob_sqrt0)
        
        # inds is the indices of the basis that has 0 at site_idx. flipped_inds is the indices of the Bell pair corresponding basis in the extended_basis. existed_inds is the indices of the flipped basis that exists in the constraint Hilbert space.
        inds = [i for i in 1:l*k_new if extended_basis_old[i][N- site_idx + 1]==0]

        @show inds
        flipped_basis = X.(extended_basis_old[inds], site_idx)
        flipped_inds = indexin(flipped_basis, T.(takesystem.(extended_basis, mask)))
        existed_inds = findall(!isnothing, flipped_inds)
        
        flipped_inds = flipped_inds[existed_inds] .+ offset
        new_state[inds] = current_state[inds]
        new_state[flipped_inds] = current_state[existed_inds]
        new_state ./= norm(new_state)
        return new_state
    else
        state_after_1 = reference_measuremap(N, 1000.0, state, site_idx, 1, pbc, k_old=k_old, measure_class = resettype)
        current_state = state_after_1 ./ sqrt(prob_sqrt1)

        inds = [i for i in 1:l*k_new if extended_basis_old[i][N- site_idx + 1]==1]
    
        flipped_basis = X.(extended_basis_old[inds], site_idx)
        flipped_inds = indexin(flipped_basis, T.(takesystem.(extended_basis, mask)))
        existed_inds = findall(!isnothing, flipped_inds)

        flipped_inds = flipped_inds[existed_inds]
        new_state[inds .+ offset] = current_state[inds]
        new_state[flipped_inds] = current_state[inds]
        new_state ./= norm(new_state)

        return new_state
    end

end

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