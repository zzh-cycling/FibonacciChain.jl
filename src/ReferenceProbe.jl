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

function reference_measuremap(::Type{T}, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true; extended_basis::Vector{newT}, k_old::Int64=1, anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET, newT}
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
    @assert newT == BitStr{N + k_old, Int} "The extended_basis should be BitStr{$(N + k_old), Int}, but got $newT"

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
reference_measuremap(N::Int, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true; extended_basis::Vector{newT}, k_old::Int64=1, anyon_type::Symbol=:Fibo) where {ET, newT} = reference_measuremap(BitStr{N, Int}, τ, state, idx, sign, pbc, k_old=k_old, anyon_type=anyon_type, extended_basis=extended_basis)

function concat_bell_pair!(N::Int, current_state::Vector{ET}, extended_basis_old, extended_basis; site_idx::Int64, k_old::Int64=1, new_dim::Int64=2^(N+k_old), verbose =false) where {ET}
    # inds is the indices of the basis that has 0 at site_idx. flipped_inds is the indices of the Bell pair corresponding basis in the extended_basis. existed is the indices of the flipped basis that exists in the constraint Hilbert space.
    inds = findall(x -> x!=0.0, current_state)

    new_state = zeros(ET, new_dim)
    T = BitStr{N + k_old, Int}
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))

    basis0 = extended_basis_old[inds] .& bmask(BitStr{N + k_old, Int}, setdiff(1:N, N-site_idx+1)...)
    basis1 = X.(basis0, site_idx)
    basis00 = join.(Ref(BitStr{1, Int}(bit"0")), basis0)
    basis11 = join.(Ref(BitStr{1, Int}(bit"1")), basis1)
    inds00 = indexin(basis00, extended_basis)
    inds11 = indexin(basis11, extended_basis)
    existed00 = findall(!isnothing, inds00)
    existed11 = findall(!isnothing, inds11)
    inds00 = inds00[existed00]
    inds11 = inds11[existed11]
    
    verbose && @show current_state, basis0, basis1, basis00, basis11, inds00, inds11, new_state

    new_state[inds00] += current_state[inds[existed00]]
    new_state[inds11] += current_state[inds[existed11]]
    new_state ./= norm(new_state)

    verbose && @show new_state

    return new_state
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
- `entangle_way::Symbol=:copy`: Method to entangle reference qubit, either `:copy` or `:reset`

# Returns
- `Vector{ET}`: New state with reference qubits in maximally entangled configuration

For multiple reference qubits, call this function multiple times at different sites.
"""
function add_reference_qubits!(N::Int, state::Vector{ET}, site_idx::Int64; k_new::Int=1, pbc::Bool=true, anyon_type::Symbol=:Fibo, verbose = false, entangle_way::Symbol = :copy) where {ET}
    # Add k_new reference qubits to the state at the specified site_idx, and place them to the left part of basis (index N-site_idx+1) to form a maximally entangled state.
    @assert 1 <= site_idx <= N "Site index must be in the range [1, N]"
    1 >= k_new >= 0 || error("k_new must be in [0,1]")
    # Because each qubit can only concat with one reference qubit, so k_new can only be 0 or 1. If need to add more reference qubits, use add_reference_qubits! multiple times at different site.
    basis_F = anyon_basis(N, pbc, anyon_type=anyon_type)
    len_F   = length(basis_F)
        

    # old reference qubit number, k_old
    k_old = round(Int, log2(length(state) ÷ len_F))
    length(state) == (2^k_old * len_F) ||
        error("state length is not compatible with (k_old, N), can not deduce k_old from state length")
    @assert k_old >= 0 "k_old must be non-negative, but got $(k_old)"
    
    # new reference prefix
    k_total = k_old + k_new
    new_dim = 2^k_total * len_F
    extended_basis = build_extended_basis(k_total, basis_F)
    extended_basis_old = build_extended_basis(k_old, basis_F)
    T = BitStr{N + k_old, Int}
    fl=bmask(T, N)
    mask = bmask(BitStr{N+k_total, Int}, 1:N...)
    X(state,i) = flip(state, fl >> (i-1))
    
    # entangle the reference qubit with the system qubit at site_idx
    if entangle_way == :reset
        verbose && @info "Using reset way to add reference qubit"
        # NEED to do RESET qubit to 0 or +!!! then concat with the new reference qubit. 
        resettable = Dict(
                :Fibo  => :resetFibo,
                :IsingX => :IsingX,
                :IsingZ => :reset,
                :IsingZZ => :reset,
            )
        anyon_type ∈ keys(resettable) || error("Unknown anyon type: $anyon_type")
        resettype = resettable[anyon_type]

        state_after_0 = reference_measuremap(N, 1000.0, state, site_idx, 0, pbc, extended_basis=extended_basis_old, k_old=k_old, anyon_type=resettype)
        state_after_1 = reference_measuremap(N, 1000.0, state, site_idx, 1, pbc, extended_basis=extended_basis_old, k_old=k_old, anyon_type=resettype)

        prob_sqrt0 = state_after_0' * state_after_0
        prob_sqrt1 = 1 - prob_sqrt0

        verbose && @show prob_sqrt0, prob_sqrt1, state_after_1

        current_statep = state_after_0 ./ sqrt(prob_sqrt0)
        current_statem = state_after_1 ./ sqrt(prob_sqrt1)

        verbose && @show current_statem

        new_statep = concat_bell_pair!(N, current_statep, extended_basis_old, extended_basis, site_idx = site_idx, k_old = k_old, new_dim = new_dim)

        new_statem = concat_bell_pair!(N, current_statem, extended_basis_old, extended_basis, site_idx = site_idx, k_old = k_old, new_dim = new_dim)

        return prob_sqrt0, prob_sqrt1, new_statep, new_statem

    elseif entangle_way == :copy
        verbose && @info "Using copy way to add reference qubit"
      
        inds = findall(x -> x!=0.0, state)
        new_state = zeros(ET, new_dim)

        basis0_inds = findall(x -> (x .& bmask(BitStr{N + k_old, Int}, N-site_idx+1)) == 0, extended_basis_old[inds])
        basis1_inds = findall(x -> (x .& bmask(BitStr{N + k_old, Int}, N-site_idx+1)) == 1 << (N-site_idx), extended_basis_old[inds])
        basis00 = join.(Ref(BitStr{1, Int}(bit"0")), extended_basis_old[inds[basis0_inds]])
        basis11 = join.(Ref(BitStr{1, Int}(bit"1")), extended_basis_old[inds[basis1_inds]])

        inds00 = indexin(basis00, extended_basis)
        inds11 = indexin(basis11, extended_basis)
        existed00 = findall(!isnothing, inds00)
        existed11 = findall(!isnothing, inds11)
        inds00 = inds00[existed00]
        inds11 = inds11[existed11]

        new_state[inds00] = state[inds[basis0_inds[existed00]]]
        new_state[inds11] = state[inds[basis1_inds[existed11]]]
        new_state ./= norm(new_state)

        return new_state
    else
        error("Unknown entangle way: $entangle_way")
    end
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
    reference_apply_measurement_layer!(N::Int64, τ::Float64, state::Vector{T}, 
    layer_sample::Vector{Int64}, 
    layer_idx::Int64, 
    pbc::Bool=true; 
    anyon_type::Symbol=:Fibo, 
    return_free_energy::Bool=false) where {T}

Generate measurement samples at 1 layer with post_selection outcomes, i.e., with time step 1.

# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `layer_sample::Vector{Int64}`: Measurement outcomes for the layer
- `layer_idx::Int64`: Layer index (1-based) to determine measurement pattern
- `pbc::Bool=true`: Periodic boundary conditions
- `extended_basis::Vector{newT}`: Extended basis with reference qubits
- `k_old::Int64=1`: Number of reference qubits in the state
- `anyon_type::Symbol=:Fibo`: anyon type
- `return_free_energy::Bool=false`: Whether to return total free energy

# Returns
`Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
- `post-measurement states`
- `measurement outcomes` 
- `free energies`

Samples measurement outcomes with given measurement outcomes.
"""
function reference_apply_measurement_layer!(N::Int64, τ::Float64, state::Vector{ET},
    layer_sample::Vector{Int64}, layer_idx::Int64,
    pbc::Bool=true; 
    extended_basis::Vector{newT}, k_old::Int64=1, 
    anyon_type::Symbol=:Fibo,
    return_free_energy::Bool=false) where {ET, newT}

    total_free_energy = zero(real(ET))
    measurement_sites, measure_type = _obtain_measurement_config(N, layer_idx, anyon_type)  

    for (idx, sign) in enumerate(layer_sample)
        state = reference_measuremap(N, τ, state, measurement_sites[idx], sign, pbc, 
        k_old=k_old, anyon_type = measure_type, extended_basis=extended_basis)

        prob = real(dot(state, state))
        total_free_energy += -log(prob)
        state ./= sqrt(prob)
    end

    return return_free_energy ? (state, total_free_energy) :  state

end

"""
    reference_sample_layer!(N::Int64, τ_eff::Float64, state::Vector{T}, 
    rng::MersenneTwister = MersenneTwister(), 
    layer_idx::Int64=1, 
    pbc::Bool=true; 
    anyon_type::Symbol=:Fibo, 
    extended_basis::Vector{newT}, k_old::Int64=1, 
    return_free_energy::Bool=true) where {T, newT}

Do the random measurement on a layer, in contrast to _apply_measurement_layer! No sample input, but output the sample.

# Arguments
- `N::Int`: Chain length N
- `τ_eff::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `rng::MersenneTwister = MersenneTwister()`: Random number generator
- `layer_idx::Int64`: Layer index (1-based) to determine measurement pattern
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type
- `return_free_energy::Bool=false`: Whether to return total free energy
- `extended_basis::Vector{newT}`: Extended basis with reference qubits
- `k_old::Int64=1`: Number of reference qubits in the state

# Returns
`Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
- `state_Born_measured`, 
- `samples`, 
- `free_energy` of one layer, if `return_free_energy=true`, otherwise only return `(state_Born_measured, samples)`.

Samples measurement outcomes with Born rule driven trajectories.
"""
function reference_sample_layer!(N::Int64, τ_eff::Float64, state::Vector{T}, 
    rng::MersenneTwister = MersenneTwister(), 
    layer_idx::Int64=1, 
    pbc::Bool=true; 
    anyon_type::Symbol=:Fibo, 
    extended_basis::Vector{newT}, k_old::Int64=1, 
    return_free_energy::Bool=false) where {T, newT}

    measurement_sites, measure_type = _obtain_measurement_config(N, layer_idx, anyon_type)  
    n = length(measurement_sites)
    sample = zeros(Int, n)
    F_layer = 0.0


    for (i, site) in enumerate(measurement_sites)
        # first 0 branch
        ψ0 = reference_measuremap(N, τ_eff, state, site, 0, pbc, 
        k_old=k_old, anyon_type = measure_type, extended_basis=extended_basis)
        p0 = real(dot(ψ0, ψ0))
        p1 = 1 - p0

        if rand(rng) < p0
            sample[i] = 0
            state = ψ0 ./ sqrt(p0)
            F_layer += -log(p0)
        else
            # else 1 branch
            ψ1 = reference_measuremap(N, τ_eff, state, site, 1, pbc, 
            k_old=k_old, anyon_type = measure_type, extended_basis=extended_basis)
            sample[i] = 1
            state = ψ1 ./ sqrt(p1)
            F_layer += -log(p1)
        end
    end
    return return_free_energy ? (state, sample, F_layer) : (state, sample)
end


"""
    reference_generate_state(τ::Float64, state::Vector{T}, sample::ET, pbc::Bool=true;
    temp::Bool=true, anyon_type::Symbol=:Fibo, verbose=false, 
    rng::MersenneTwister=MersenneTwister(),
    mode::Symbol=:sample, return_free_energy::Bool=false) where{T, ET}

Generate quantum state evolution under measurement protocol with reference qubits.

# Arguments
- `N::Int64`: System size
- `τ::Float64`: Evolution time parameter
- `state::Vector{ET}`: Initial quantum state with reference qubits
- `sample`: Measurement sample configuration
- `pbc::Bool=true`: Periodic boundary conditions
- `temp::Bool=true`: Return trajectory if true, final state if false
- `anyon_type::Symbol=:Fibo`: Model type
- `verbose::Bool=false`: Enable verbose output
- `rng::MersenneTwister=MersenneTwister()`: Random number generator
- `mode::Symbol = :sample`: Evolution mode, either `:sample` for given trajectory or `:Born` for Born rule driven
- `return_free_energy::Bool=false`: Whether to return total free energy

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
function reference_generate_state(τ::Float64, state::Vector{T}, sample::ET, pbc::Bool=true;
    temp::Bool=true, anyon_type::Symbol=:Fibo, verbose=false, 
    rng::MersenneTwister=MersenneTwister(),
    mode::Symbol=:sample, return_free_energy::Bool=false) where{T, ET}

    @assert ET == Matrix{Int} "ET must be Matrix{Int} for reference_generate_state"

    D = size(sample, 1)
    N = (anyon_type == :Fibo) ? 2*round(Int, size(sample, 2)) : size(sample, 2)
    
    statelis = temp ? Vector{Vector{T}}(undef, D) : nothing
    sample_free_energy = zeros(Float64, D) 

    basis_F = anyon_basis(N, pbc, anyon_type=anyon_type)

    len_F   = length(basis_F)
    # old reference qubit number, k_old
    k_old = round(Int, log2(length(state) ÷ len_F))
    
    # Noting that basis is not consisten with state, but extended_basis is.
    extended_basis = build_extended_basis(k_old, basis_F) 

    @assert length(basis_F)*2^k_old == length(state) == length(extended_basis) "state length is expected to be $(length(basis_F)*2^k_old), but got $(length(state))"

    if mode == :Born
        verbose && @info "Using Born rule driven evolution mode"
        for layer in 1:D
            verbose && @info "Evolving layer $layer / $D"
            τ_eff = (layer == D) ? τ/2 : τ

            state, sample[layer, :], sample_free_energy[layer] = reference_sample_layer!(N, τ_eff, state, rng, layer, pbc, k_old=k_old, anyon_type = anyon_type, extended_basis=extended_basis, return_free_energy=true)

            temp && (statelis[layer] = copy(state))
        end

    elseif mode == :sample
        verbose && @info "Using given sample evolution mode"

        for layer in 1:D
            verbose && @info "Evolving layer $layer / $D"
            τ_eff = (layer == D) ? τ/2 : τ

            state, sample_free_energy[layer] = reference_apply_measurement_layer!(N, τ_eff, state, sample[layer, :], layer, pbc, k_old=k_old, anyon_type = anyon_type, extended_basis=extended_basis, return_free_energy=true)
     
            temp && (statelis[layer] = copy(state))
        end
    end  
    
    final_st = temp ? statelis : state

    return return_free_energy ? (final_st, sample, sample_free_energy) : (final_st, sample)
end


"""
    reference_evolution(N::Int, τ::Float64, forward::Vector{ET}, sample::Matrix{Int}, 
    x₂::Int, t₁, t₂; x₁::Int=1, 
    rng=Random.default_rng(), pbc=true, 
    anyon_type::Symbol=:Fibo, temp::Bool=false, verbose=false, 
mode::Symbol = :sample, return_free_energy::Bool=false) where {ET}

Compute temporal correlation between two time slices using cached forward evolution and given trajectory (doing the post selection).

# Arguments
- `N::Int`: System size
- `τ`: Evolution time parameter
- `forward::Vector{ET}`: Cached forward state evolution trajectory
- `sample::Matrix{Int}`: Measurement sample configuration
- `x₂`: Site index for reference qubit insertion
- `t₁`: First time slice index
- `t₂`: Second time slice index (must be >= t₁)
- `x₁::Int=1`: Spatial site index for first reference qubit
- `rng=Random.default_rng()`: Random number generator
- `pbc=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type
- `temp::Bool=false`: Return trajectory if true
- `verbose::Bool=false`: Enable verbose output
- `mode::Symbol = :sample`: Evolution mode, either `:sample` for given trajectory or `:Born` for Born rule driven
- `return_free_energy::Bool=false`: Whether to return total free energy
    
# Returns
- Reference qubit added state between specified spacetime slices.

Avoids redundant computation by reusing forward evolution results.
"""
function reference_evolution(N::Int, τ::Float64, forward::Vector{ET}, sample::Matrix{Int}, 
    x₂::Int, t₁, t₂; x₁::Int=1, 
    rng=Random.default_rng(), pbc=true, 
    anyon_type::Symbol=:Fibo, temp::Bool=false, verbose=false, 
mode::Symbol = :sample, return_free_energy::Bool=false) where {ET}
    # This function is used to evolve state to compute the spatial-temporal correlation at different space and time slices cache, avoiding the repeated calculation of the state evolution. INPUT the forward state evolution, given the trajectory.
    # t₁ and t₂ are the indices of the time slices in the sample. Spacetime is:
    # | ---->|
    # forward
    #       Ref1                    x₁ = 1            
    #        |
    #        |
    #       Ref2 --> Ref3               x₂ = N/2
    # | ----> |______| -----> |
    # 1      t₁   t₂=t₁+δt   D

    n_measure = (anyon_type == :Fibo) ? N÷2 : N
    D = size(sample, 1)  

    @assert size(sample, 2) == n_measure "Sample size spatial dimension must be $n_measure, but got $(size(sample, 2))"
    @assert 1 <= t₁ <= t₂ <= D "Time slice t₁ must before time slice t₂, both must be in the range [1, $(D)]"
    @assert 1 <= x₁ <= x₂ <= N "Site index x₁ must be smaller than site index x₂, both must be in the range [1, $(N)]"
    @assert mode ∈ [:sample, :Born] "mode must be either :sample or :Born, but got $mode"

    δt = t₂ - t₁ # if δt = 0, pure spatial correlation
    δx = abs(x₂ - x₁) # if δx = 0, pure temporal correlation
    
    # 1) 0 → t₁, the steady state, at the t₁ of forward evolution
    state = forward[t₁]
    temp && (statelis = Vector{ET}(undef, D) )
    temp && (view(statelis, 1:t₁) .= view(forward, 1:t₁))
    sample_layer = zeros(Int, size(sample, 1), n_measure)
    view(sample_layer, 1:t₁, :) .= view(sample, 1:t₁, :)
    sample_free_energy= zeros(Float64, D)
    view(sample_free_energy, 1:t₁) .= 0.0

    if δt > 0 && δx > 0 # 3 ref qubits, both spatial and temporal correlation, actually 3-point correlation.
        verbose && @info "t₁ = $(t₁), t₂ = $(t₂), x₁ = $(x₁), x₂ = $(x₂), 3 refs"

        # 2) add reference qubit 1 at x₁, and reference qubit 2 at x₂
        state1 = add_reference_qubits!(N, state, x₁, pbc=pbc,
                                       anyon_type=anyon_type, verbose=verbose)
        state2 = add_reference_qubits!(N, state1, x₂, pbc=pbc,
                                       anyon_type=anyon_type, verbose=verbose)
        
        # 3) t₁ → t₂ evolution, or δt
        
        final_stlis1, samples1, free_energy1 = reference_generate_state(τ, state2, sample[t₁+1:t₂, :], pbc, rng=rng, anyon_type=anyon_type, temp=true, verbose=verbose, mode=mode, return_free_energy=true)
        
       
        # 4) add reference qubit 3 at x₂
        state3 = add_reference_qubits!(N, final_stlis1[end], x₂, pbc=pbc, anyon_type=anyon_type, verbose=verbose) 

        # 5) t₂ → D evolution
        final_stlis2, samples2, free_energy2 = reference_generate_state(τ, state3, sample[t₂+1:end, :], pbc, rng=rng, anyon_type=anyon_type, temp=true, verbose=verbose, mode=mode, return_free_energy=true)

        temp && (view(statelis, t₁+1:t₂) .= view(final_stlis1, :))
        temp && (view(statelis, t₂+1:D, :) .= view(final_stlis2, :))

        sample_layer[t₁+1:t₂, :] .= samples1
        sample_layer[t₂+1:end, :] .= samples2
        sample_free_energy[t₁+1:t₂] .= free_energy1
        sample_free_energy[t₂+1:end] .= free_energy2
        final_stlis3 = temp ? statelis : final_stlis2[end]

    elseif δt == 0 # 2 ref qubits, pure 2-point spatial correlation
        verbose && @info "x₁ = $(x₁), x₂ = $(x₂), δx = $(δx), at time slice t₁ = t₂ = $(t₁), 2 refs"
    
        # 2) add reference qubit 1 at x₁, and reference qubit 2 at x₂
        state1 = add_reference_qubits!(N, state, x₁, pbc=pbc,
                                       anyon_type=anyon_type, verbose=verbose)
        state2 = add_reference_qubits!(N, state1, x₂, pbc=pbc,
                                       anyon_type=anyon_type, verbose=verbose)
        
        # 3) t₁ → D evolution
        final_stlis2, samples2, free_energy2 = reference_generate_state(τ, state2, sample[t₁+1:end, :], pbc, rng=rng, anyon_type=anyon_type, temp=true, verbose=verbose, mode=mode, return_free_energy=true)

        temp && (view(statelis, t₁+1:D) .= view(final_stlis2, :))
        sample_layer[t₁+1:end, :] .= samples2
        sample_free_energy[t₁+1:end] .= free_energy2
        final_stlis3 = temp ? statelis : final_stlis2[end]

    elseif δx == 0 # 2 ref qubits, pure 2-point temporal correlation
        verbose && @info "t₁ = $(t₁), t₂ = $(t₂), δt = $(δt), at site x₁ = x₂ = $(x₂), 2 refs"

        # 2) add reference qubit 1 at x₁
        state1 = add_reference_qubits!(N, state, x₂, pbc=pbc, anyon_type=anyon_type, verbose=verbose)

        # 3) t₁ → t₂ evolution, or δt
        final_stlis1, samples1, free_energy1 = reference_generate_state(τ, state1, sample[t₁+1:t₂, :], pbc, rng=rng, anyon_type=anyon_type, temp=true, verbose=verbose, mode=mode, return_free_energy=true)

        # 4) add reference qubit 2 at x₂
        state2 = add_reference_qubits!(N, final_stlis1[end], x₂, pbc=pbc, anyon_type=anyon_type, verbose=verbose)
        
        # 5) t₂ → D evolution
        final_stlis2, samples2, free_energy2 = reference_generate_state(τ, state2, sample[t₂+1:end, :], pbc, rng=rng, anyon_type=anyon_type, temp=true, verbose=verbose, mode=mode, return_free_energy=true)

        temp && (view(statelis, t₁+1:t₂) .= view(final_stlis2, :))
        temp && (view(statelis, t₂+1:D) .= view(final_stlis2, :))
        sample_layer[t₁+1:t₂, :] .= samples1
        sample_layer[t₂+1:end, :] .= samples2
        sample_free_energy[t₁+1:t₂] .= free_energy1
        sample_free_energy[t₂+1:end] .= free_energy2

        final_stlis3 = temp ? statelis : final_stlis2[end]
        
    end
    
    return return_free_energy ? (final_stlis3, sample_layer, sample_free_energy) : (final_stlis3, sample_layer)
end