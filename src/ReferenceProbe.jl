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

function reference_measure_basismap(::Type{T}, τ::Float64, state::ET, i::Int, sign::Int64, pbc::Bool=true; k_old::Int64=1, anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}
    # default for PBC system, map basis
    @assert k_old >= 0 "k_old must be at least 0, but got $(k_old)"
 
    mask = bmask(BitStr{N+k_old, Int}, 1:N...)
    action_state = T(takesystem(state, mask))
    return measure_basismap(T, τ, action_state, i, sign, pbc, anyon_type=anyon_type)
   
end

function reference_measuremap(::Type{T}, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true;k_old::Int64=1,anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}
    # input a superposition state with reference qubit, and output the measured state. k_old is the number of reference qubits in the state.
    if anyon_type == :Fibo
        @assert pbc || (2 <= idx <= N-1) "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif anyon_type == :IsingZZ
        @assert pbc || (1 <= idx <= N-1) "Index idx must be in [1, N-1] for open BC (IsingZZ)"
    elseif anyon_type ∈ (:IsingX, :IsingZ, :reset, :resetFibo)
        @assert pbc || (1 <= idx <= N) "Index idx must be in [1, N] for open BC (IsingX)"
    else
        error("Unknown measure class: $anyon_type")
    end
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    basis=anyon_basis(T, pbc, anyon_type=anyon_type)
    
    # Noting that basis is not consisten with state, but extended_basis is.
    extended_basis = build_extended_basis(k_old, basis) 
    @assert length(basis)*2^k_old == length(state) == length(extended_basis) "state length is expected to be $(length(basis)*2^k_old), but got $(length(state))"
    l=length(extended_basis)

    mapped_state = zeros(ET, length(state))

    pretype = BitStr{k_old, Int}
    mask = bmask(BitStr{N+k_old, Int}, 1:N...)
    
    for (i, ext_basis_i) in enumerate(extended_basis)
        output = reference_measure_basismap(T, τ, ext_basis_i, idx, sign, pbc, k_old=k_old, anyon_type=anyon_type)
        
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
reference_measuremap(N::Int, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true; k_old::Int64=1, anyon_type::Symbol=:Fibo) where {ET} = reference_measuremap(BitStr{N, Int}, τ, state, idx, sign, pbc, k_old=k_old, anyon_type=anyon_type)

function reference_apply_measurement_layer!(N::Int64, state::Vector{T}, τ::Float64, layer_sample::Vector{Int64}, layer_idx::Int64, pbc::Bool=true; k_old::Int64=1, anyon_type::Symbol=:Fibo) where {T}
    if anyon_type == :Fibo
        if layer_idx % 2 == 1
            measurement_sites = collect(2:2:N)  # odd sites anyons, even sites qubits
        else
            measurement_sites = collect(1:2:N)  # even sites anyons, odd sites qubits
        end
        for (idx, measurement_type) in enumerate(layer_sample)
            state = reference_measuremap(N, τ, state, measurement_sites[idx], measurement_type, pbc, k_old=k_old, anyon_type = anyon_type)
            normalize!(state)
        end
        return state
    
    elseif anyon_type == :IsingX || anyon_type == :IsingZZ
        measurement_sites = collect(1:N)
        if layer_idx % 2 == 1
            # odd layers: measure X
            for (idx, measurement_type) in enumerate(layer_sample)
                state = reference_measuremap(N, τ, state, measurement_sites[idx], measurement_type, pbc, k_old=k_old, anyon_type = :IsingX)
                normalize!(state)
            end
            return state
        else
            # even layers: measure ZZ
            for (idx, measurement_type) in enumerate(layer_sample)
                state = reference_measuremap(N, τ, state, measurement_sites[idx], measurement_type, pbc, k_old= k_old, anyon_type = :IsingZZ)
                normalize!(state)
            end
            return state
        end
    else
        error("Unknown measure class: $anyon_type")
    end
end


"""
    reference_generate_state(N::Int64, τ::Float64, state::Vector{ET}, sample; pbc::Bool=true, temp::Bool=true, k_old::Int64=1, anyon_type::Symbol=:Fibo) where {ET}

Generate quantum state evolution under measurement protocol with reference qubits.

# Arguments
- `N::Int64`: System size
- `τ::Float64`: Evolution time parameter
- `state::Vector{ET}`: Initial quantum state with reference qubits
- `sample`: Measurement sample configuration
- `pbc::Bool=true`: Periodic boundary conditions
- `temp::Bool=true`: Return trajectory if true, final state if false
- `k_old::Int64=1`: Number of existing reference qubits
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `Vector{ET}` or trajectory: Evolved state or time evolution trajectory

# Examples
```jldoctest
julia> using FibonacciChain, Random, LinearAlgebra

julia> N = 4;

julia> initial_state = normalize!(ones(Float64, length(anyon_basis(N))));

julia> add_ref = add_reference_qubits!(N, initial_state, 1, verbose=false);

julia> Random.seed!(42);

julia> sample = ones(Int, 2, 2);

julia> τ = 0.5;

julia> trajectory = reference_generate_state(τ, add_ref, sample, temp=true);

julia> length(trajectory) == size(sample, 1)
true
```
"""
function reference_generate_state(τ::Float64, state::Vector{T}, sample::ET, pbc::Bool=true; temp::Bool=true, k_old::Int64=1, anyon_type::Symbol=:Fibo, verbose=false) where{T, ET}
    @assert ET == Matrix{Int} "ET must be Matrix{Int} for reference_generate_state"

    D = size(sample, 1)
    N = (anyon_type == :Fibo) ? 2*round(Int, size(sample, 2)) : size(sample, 2)
    
    statelis = temp ? Vector{Vector{T}}(undef, D) : nothing
 
    for layer in 1:D
        τ_eff = (layer == D) ? τ/2 : τ
        state = reference_apply_measurement_layer!(N, state, τ_eff, sample[layer, :], layer, pbc, k_old=k_old, anyon_type = anyon_type)

        if temp
            statelis[layer] = copy(state)
        end
    end
    
    return temp ? statelis : state
end

"""
    add_reference_qubits!(N::Int, state::Vector{ET}, site_idx::Int64, rng::MersenneTwister=MersenneTwister(); k_new::Int=1, pbc::Bool=true, anyon_type::Symbol=:Fibo) where {ET}

Add reference qubits to quantum state at specified site for correlation measurements.

# Arguments
- `N::Int`: System size
- `state::Vector{ET}`: Quantum state vector
- `site_idx::Int64`: Site index for reference qubit insertion (1 ≤ site_idx ≤ N)
- `k_new::Int=1`: Number of new reference qubits to add (0 or 1)
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type
- `verbose::Bool=false`: Enable verbose output

# Returns
- `Vector{ET}`: New state with reference qubits in maximally entangled configuration

For multiple reference qubits, call this function multiple times at different sites.
"""
function add_reference_qubits!(N::Int, state::Vector{ET}, site_idx::Int64; k_new::Int=1, pbc::Bool=true, anyon_type::Symbol=:Fibo, verbose = false) where {ET}
    # Add k_new reference qubits to the state at the specified site_idx, and place them to the left part of basis (index N-site_idx+1) to form a maximally entangled state.
    @assert 1 <= site_idx <= N "Site index must be in the range [1, N]"
    1 >= k_new >= 0 || error("k_new must be in [0,1]")
    # Because each qubit can only concat with one reference qubit, so k_new can only be 0 or 1. If need to add more reference qubits, use add_reference_qubits! multiple times at different site.
    basis_F = anyon_basis(N, pbc, anyon_type=anyon_type)
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
    
    # NEED to do RESET qubit to 0 or +!!! then concat with the new reference qubit. 
    resettable = Dict(
            :Fibo  => :resetFibo,
            :IsingX => :IsingX,
            :IsingZ => :reset,
            :IsingZZ => :reset,
        )
    anyon_type ∈ keys(resettable) || error("Unknown anyon type: $anyon_type")
    resettype = resettable[anyon_type]
    
    state_after_0 = reference_measuremap(N, 1000.0, state, site_idx, 0, pbc,k_old=k_old, anyon_type = resettype)
    prob_sqrt0 = state_after_0' * state_after_0

    verbose && @show random_number, prob_sqrt0

    current_state = state_after_0 ./ sqrt(prob_sqrt0)
    
    # inds is the indices of the basis that has 0 at site_idx. flipped_inds is the indices of the Bell pair corresponding basis in the extended_basis. existed is the indices of the flipped basis that exists in the constraint Hilbert space.
    inds = findall(x -> x!=0.0, current_state)

    flipped_basis = X.(extended_basis_old[inds], site_idx)
    added_basis = join.(BitStr{1, Int}.(readbit.(extended_basis_old[inds], N-site_idx+1)), extended_basis_old[inds])
    added_flipped_basis = join.(BitStr{1, Int}.(readbit.(flipped_basis, N-site_idx+1)), flipped_basis)
    
    flipped_inds = indexin(added_flipped_basis, extended_basis)
    added_inds = indexin(added_basis, extended_basis)
    existed = findall(!isnothing, flipped_inds)
    existed_added = findall(!isnothing, added_inds)
    flipped_inds = flipped_inds[existed]
    added_inds = added_inds[existed_added]
    
    new_state[added_inds] += current_state[inds[existed_added]]
    new_state[flipped_inds] += current_state[inds[existed]]
    new_state ./= norm(new_state)
    return new_state
    
end


# subsystems is the system to keep, not throw!!!
function reference_rdm(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}; pbc::Bool=true, anyon_type::Symbol=:Fibo, traceref::Bool=true) where {N, T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    # N is the particle number of system, while k_old is the number of reference qubit, which is deduced from the state length.
    unsorted_basis = anyon_basis(T, pbc; anyon_type=anyon_type)
    len_F   = length(unsorted_basis)
    k_old = round(Int, log2(length(state) ÷ len_F))
    
    length(state) == (2^k_old * len_F) || error("state length is not compatible with (k_old, N), can not deduce k_old from state length")
    @assert 2^k_old*length(unsorted_basis) == length(state) "state length is expected to be $(2^k_old*length(unsorted_basis)), but got $(length(state))"

    if traceref
        # If traceref is true, we need to trace out the reference qubit. otherwise, we trace out system.
        return disjoint_rdm(BitStr{k_old, Int}, T, subsystems, Int[], state, pbc; anyon_typeA=:IsingX, anyon_typeB=anyon_type)
    else
        totalsubBpbc = (length(subsystems) == N) ? true : false
        return disjoint_rdm(BitStr{k_old, Int}, T, Int[], subsystems, state, pbc; anyon_typeA=:IsingX, totalsubBpbc=totalsubBpbc, anyon_typeB=anyon_type)
    end
   
  
end
reference_rdm(N::Int, subsystems::Vector{Int}, state::Vector{ET}; pbc::Bool=true, anyon_type::Symbol=:Fibo, traceref::Bool=true) where {ET} = reference_rdm(BitStr{N, Int}, subsystems, state, pbc=pbc, anyon_type=anyon_type, traceref=traceref)

"""
    reference_evolution(τ, forward, sample, site, time_slice1, time_slice2; pbc=true, seed::Int64=100, anyon_type::Symbol=:Fibo, temp::Bool=false)

Compute temporal correlation between two time slices using cached forward evolution.

# Arguments
- `τ`: Evolution time parameter
- `forward`: Cached forward state evolution trajectory
- `sample`: Measurement sample configuration
- `site`: Site index for reference qubit insertion
- `time_slice1`: First time slice index
- `time_slice2`: Second time slice index (must be > time_slice1)
- `pbc=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type
- `temp::Bool=false`: Return trajectory if true
- `verbose::Bool=false`: Enable verbose output

# Returns
- Reference qubit added state between specified time slices

Avoids redundant computation by reusing forward evolution results.
"""
function reference_evolution(τ, forward, sample, site, time_slice1, time_slice2;
pbc=true, anyon_type::Symbol=:Fibo, temp::Bool=false, verbose=false)
    # This function is used to compute the temporal_correlation at different time slices cache, avoiding the repeated calculation of the state evolution. INPUT the forward state evolution.
    # time_slice1 and time_slice2 are the indices of the time slices in the sample.
    N = (anyon_type == :Fibo) ? 2*round(Int, size(sample, 2)) : size(sample, 2)
    D = size(sample, 1)
    @assert 1 <= site <= N "Site index must be in the range [1, N]"    
    @assert 1 <= time_slice1 <= D "Time slice 1 index must be in the range [1, $(D)]"
    @assert 1 <= time_slice2 <= D "Time slice 2 index must be in the range [1, $(D)]"
    @assert time_slice1 < time_slice2 "Time slice 1 must before time slice 2"

    # 1) 0 → t₁
    state = forward[time_slice1]

    # 2) add reference qubit 1 at site
    state1 = add_reference_qubits!(N, state, site, pbc=pbc,
                                   anyon_type=anyon_type, verbose=verbose)
    
    # 3) t₁ → t₂ evolution
    final_stlis1 = reference_generate_state(τ, state1, sample[time_slice1+1:time_slice2, :], pbc, anyon_type=anyon_type, temp=temp)
    
    if temp
        statelis = Vector{eltype(forward)}(undef, D)
        # 4) add reference qubit 2 at site
        state2 = add_reference_qubits!(N, final_stlis1[end], site, pbc=pbc, anyon_type=anyon_type, verbose=verbose) 

        # 5) t₂ → D evolution
        final_stlis2 = reference_generate_state(τ, state2, sample[time_slice2+1:end, :], pbc, k_old=2, anyon_type=anyon_type, temp=temp)
        statelis[1:time_slice1] = forward[1:time_slice1]
        statelis[time_slice1+1:time_slice2] = final_stlis1
        statelis[time_slice2+1:end] = final_stlis2
        return statelis
    else
        state2 = add_reference_qubits!(N, final_stlis1, site, pbc=pbc, anyon_type=anyon_type)
        final_stlis2 = reference_generate_state(τ, state2, sample[time_slice2+1:end, :], pbc, k_old=2, anyon_type=anyon_type, temp=temp)
        return final_stlis2
    end
end
