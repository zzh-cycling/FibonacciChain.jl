function build_extended_basis(basis::Vector{ET}, ::Type{T}) where {k_total, ET, T <: BitStr{k_total, Int}}
    ref_strings = [T(i) for i in 0:(2^k_total-1)] 
    extended_basis = sort(
        mapreduce(suffix -> process_join([suffix], basis),
              vcat,
              ref_strings))

    return extended_basis
end
build_extended_basis(k_total::Int, basis::Vector{ET}) where {ET} = build_extended_basis(basis, BitStr{k_total, Int})
# process_join([suffix], basis) will give [0000, 0001, 0010, 0011], but process_join(basis, [suffix]) will give [0000, 0001, 0010, 0100]

function reference_measure_basismap(model::AnyonModel, ::Type{T}, ::Type{newT}, τ::Float64, state::ET, i::Int, sign::Bool; k_old::Int64=1) where {N, M, T <: BitStr{N}, ET, newT <: BitStr{M}}
    # default for PBC system, map basis
    @assert k_old >= 0 "k_old must be at least 0, but got $(k_old)"
    @assert M == N + k_old "The output basis should be with length $(N + k_old), but got $M"
    
    mask = bmask(newT, 1:N...)
    action_state = T(takesystem(state, mask))
    return measure_basismap(model, τ, action_state, i, sign)
    
end

num_digits(::Type{<:BitStr{N}}) where N = N
function reference_measuremap(model::AnyonModel, ::Type{T}, ::Type{pretype}, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool; extended_basis::Vector{newT}) where {N, k_old, T <: BitStr{N}, ET, newT <: BitStr, pretype <: BitStr{k_old, Int}}
    # input a superposition state with reference qubit, and output the measured state. k_old is the number of reference qubits in the state.
    if model.measure_operator ∈ [:Ferro, :Antiferro]
        @assert model.pbc || (2 <= idx <= model.N-1) "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif model.measure_operator == :ZZ
        @assert model.pbc || (1 <= idx <= model.N-1) "Index idx must be in [1, N-1] for open BC (ZZ)"
    elseif model.measure_operator ∈ (:X, :Z, :reset, :resetFibo)
        @assert model.pbc || (1 <= idx <= model.N) "Index idx must be in [1, N] for open BC (X)"
    else
        error("Unknown measure class: $(model.anyon_type)")
    end
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"
    @assert num_digits(newT) == N + k_old "The output basis should be with length $(N + k_old), but got $(num_digits(newT))"

    mapped_state = zeros(ET, length(state))

    mask = bmask(newT, 1:N...)
    
    for (i, ext_basis_i) in enumerate(extended_basis)
        outputstate1, outputstate2, output1, output2 = reference_measure_basismap(model, T, newT, τ, ext_basis_i, idx, sign, k_old=k_old)

        prefix_i = pretype(takeenviron(ext_basis_i, mask) >> N)

        if output2 ==0
            mapped_state[i]+=output1*state[i]
        else
            j2=searchsortedfirst(extended_basis, join(prefix_i, outputstate2))
            mapped_state[i]+=output1*state[i] # outputstate1 is the same as basis[i]
            mapped_state[j2]+=output2*state[i]
        end
    end

    return mapped_state
end

"""
    reference_measuremap(model::AnyonModel, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool;
                         extended_basis::Vector{newT}, k_old::Int64=1)
    reference_measuremap(model::AnyonModel, ::Type{T}, ::Type{pretype}, τ::Float64, state::Vector{ET},
                         idx::Int, sign::Bool; extended_basis::Vector{newT})

Apply a measurement operator to a quantum state with reference qubits.

This function extends `measuremap` to handle states that include reference qubits, 
preserving the reference qubit subspace while applying measurements to the system.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Quantum state vector with reference qubits
- `idx::Int`: Site index for measurement
- `sign::Bool`: Measurement outcome (false=0, true=1)
- `extended_basis::Vector{newT}`: Extended basis including reference qubits
- `k_old::Int64=1`: Number of reference qubits in the state

# Returns
- `Vector{ET}`: The post-measurement state (unnormalized)

# Notes
The extended basis must be constructed via `build_extended_basis` before calling this function.
"""
reference_measuremap(model::AnyonModel, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool; extended_basis::Vector{newT}, k_old::Int64=1) where {ET, newT} = reference_measuremap(model, BitStr{model.N, Int}, BitStr{k_old, Int}, τ, state, idx, sign, extended_basis=extended_basis)

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
function add_reference_qubits!(model::AnyonModel, state::Vector{ET}, site_idx::Int64; k_new::Int=1, verbose = false, entangle_way::Symbol = :copy) where {ET}
    # Add k_new reference qubits to the state at the specified site_idx, and place them to the left part of basis (index N-site_idx+1) to form a maximally entangled state.
    N = model.N
    @assert 1 <= site_idx <= N "Site index must be in the range [1, $(N)]"
    1 >= k_new >= 0 || error("k_new must be in [0,1]")
    # Because each qubit can only concat with one reference qubit, so k_new can only be 0 or 1. If need to add more reference qubits, use add_reference_qubits! multiple times at different site.
    basis_F = anyon_basis(model)
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

        state_after_0 = reference_measuremap(model, 1000.0, state, site_idx, false, extended_basis=extended_basis_old, k_old=k_old)
        state_after_1 = reference_measuremap(model, 1000.0, state, site_idx, true,extended_basis=extended_basis_old, k_old=k_old)

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


"""
    reference_rdm(model::AnyonModel, subsystems::Vector{Int64}, state::Vector{ET}; traceref::Bool=true)

Compute the reduced density matrix for a state with reference qubits.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `subsystems::Vector{Int64}`: Indices of system sites to keep (not trace out)
- `state::Vector{ET}`: Quantum state vector with reference qubits
- `traceref::Bool=true`: If `true`, trace out reference qubits and keep system subsystems;
                          if `false`, trace out system and keep reference qubits

# Returns
- `Matrix{ET}`: The reduced density matrix

# Notes
- Subsystem indices count from the right of the binary string representation
- The number of reference qubits `k_old` is automatically deduced from the state length
"""
function reference_rdm(model::AnyonModel, subsystems::Vector{Int64}, state::Vector{ET}; traceref::Bool=true) where {ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    # N is the particle number of system, while k_old is the number of reference qubit, which is deduced from the state length.
    unsorted_basis = anyon_basis(model)
    len_F   = length(unsorted_basis)
    k_old = round(Int, log2(length(state) ÷ len_F))
    
    length(state) == (2^k_old * len_F) || error("state length is not compatible with (k_old, N), can not deduce k_old from state length")
    @assert 2^k_old*length(unsorted_basis) == length(state) "state length is expected to be $(2^k_old*length(unsorted_basis)), but got $(length(state))"

    if traceref
        # If traceref is true, we need to trace out the reference qubit. otherwise, we trace out system.
        ref_model = AnyonModel(IsingAnyon(), k_old, pbc=false)
        return disjoint_rdm(ref_model, model, subsystems, Int[], state;)
    else
        ref_model = AnyonModel(IsingAnyon(), k_old, pbc=false)
        totalsubBpbc = (length(subsystems) == N) ? true : false
        model = AnyonModel(model.anyon_type, model.N, pbc=totalsubBpbc)
        return disjoint_rdm(ref_model, model, Int[], subsystems, state;)
    end
   
  
end


"""
    reference_boundary_evolution(model::AnyonModel, state::Vector{T}, measure_config::MeasureConfig,
                                 sample::Union{Nothing, BitVector}=nothing; layer_idx::Int64=1)

Perform a single measurement layer evolution on a state with reference qubits.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `state::Vector{T}`: Quantum state vector with reference qubits
- `measure_config::MeasureConfig`: Configuration containing τ, mode, rng, etc.
- `sample::Union{Nothing, BitVector}=nothing`: Predefined measurement outcomes for `:sample` mode
- `layer_idx::Int64=1`: Layer index (1-based) to determine measurement pattern

# Returns
- `Measurement_outcome_boundary`: A struct containing:
  - `state::Vector{T}`: Post-measurement state
  - `sample::BitVector`: Measurement outcomes for this layer
  - `free_energy::Float64`: Free energy contribution from this layer

# Notes
- In `:Born` mode, samples are generated probabilistically via Born rule
- In `:sample` mode, `sample` must be provided as input
"""
function reference_boundary_evolution(model::AnyonModel, state::Vector{T}, measure_config::MeasureConfig, sample::Union{Nothing, BitVector} = nothing; layer_idx::Int64=1) where{T} 
    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")


    basis_F = anyon_basis(model)

    len_F   = length(basis_F)
    # old reference qubit number, k_old
    k_old = round(Int, log2(length(state) ÷ len_F))
    
    # Noting that basis is not consisten with state, but extended_basis is.
    extended_basis = build_extended_basis(k_old, basis_F) 

    τ_eff = measure_config.enable_τ_eff ? τ/2 : τ
    if mode == :sample
        N = anyon_model.N
        size(sample, 1) == measurement_num(anyon_model.anyon_type)*(N ÷ 2) || error("sample size mismatch with anyon_model $(N)")
        verbose && @info "Using given sample evolution mode"
        _reference_apply_measurement_layer(model, τ_eff, state, sample, layer_idx; extended_basis=extended_basis, k_old=k_old)
    elseif mode == :Born
        verbose && @info "Using Born rule driven evolution mode"
        _reference_sample_layer(model, τ_eff, state, rng, layer_idx; extended_basis=extended_basis, k_old=k_old)
    end
end


"""
    _reference_apply_measurement_layer!(N::Int64, τ::Float64, state::Vector{ET},
    layer_sample::Vector{Int64}, layer_idx::Int64,
    pbc::Bool=true; 
    extended_basis::Vector{newT}, k_old::Int64=1, 
    anyon_type::Symbol=:Fibo) where {ET, newT}

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

# Returns
`Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}`: 
- `post-measurement states`
- `free energies`

Samples measurement outcomes with given measurement outcomes.
"""
function _reference_apply_measurement_layer(model::AnyonModel, τ::Float64, state::Vector{ET},
    layer_sample::Vector{Bool}, layer_idx::Int64;
    extended_basis::Vector{newT}, k_old::Int64=1) where {ET, newT}

    total_free_energy = zero(real(ET))
    measurement_sites, measure_anyon_model = _obtain_measurement_config(model, layer_idx)  

    for (idx, sign) in enumerate(layer_sample)
        state = reference_measuremap(measure_anyon_model, τ, state, measurement_sites[idx], sign, extended_basis=extended_basis, k_old=k_old)

        prob = real(dot(state, state))
        total_free_energy += -log(prob)
        state ./= sqrt(prob)
    end

    return Measurement_outcome_boundary(state, layer_sample, total_free_energy)
end


"""
    _reference_sample_layer!(N::Int64, τ_eff::Float64, state::Vector{T}, 
    rng::MersenneTwister = MersenneTwister(), 
    layer_idx::Int64=1, 
    pbc::Bool=true; 
    anyon_type::Symbol=:Fibo, 
    extended_basis::Vector{newT}, k_old::Int64=1) where {T, newT}

Do the random measurement on a layer, in contrast to _apply_measurement_layer! No sample input, but output the sample.

# Arguments
- `N::Int`: Chain length N
- `τ_eff::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `rng::MersenneTwister = MersenneTwister()`: Random number generator
- `layer_idx::Int64`: Layer index (1-based) to determine measurement pattern
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type
- `extended_basis::Vector{newT}`: Extended basis with reference qubits
- `k_old::Int64=1`: Number of reference qubits in the state

# Returns
`Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}`: 
- `state_Born_measured`, 
- `samples`, 
- `free_energy` of one layer.

Samples measurement outcomes with Born rule driven trajectories.
"""
function _reference_sample_layer(model::AnyonModel, τ_eff::Float64, state::Vector{T}, 
    rng::MersenneTwister = MersenneTwister(), 
    layer_idx::Int64=1; 
    extended_basis::Vector{newT}, k_old::Int64=1, verbose::Bool=false) where {T, newT}

    measurement_sites, measure_anyon_model = _obtain_measurement_config(model, layer_idx)
    n = length(measurement_sites)
    sample = BitVector(zeros(Bool, n))
    F_layer = 0.0


    for (i, site) in enumerate(measurement_sites)
        # first 0 branch
        ψ0 = reference_measuremap(measure_anyon_model, τ_eff, state, site, false, extended_basis=extended_basis, k_old=k_old)
        p0 = real(dot(ψ0, ψ0))
        p1 = 1 - p0

        random_number = rand(rng)
        verbose && @show random_number
        if random_number < p0
            sample[i] = 0
            state = ψ0 ./ sqrt(p0)
            F_layer += -log(p0)
            verbose && @show -log(p0)
        else
            # else 1 branch
            ψ1 = reference_measuremap(measure_anyon_model, τ_eff, state, site, true, extended_basis=extended_basis, k_old=k_old)
            sample[i] = 1
            state = ψ1 ./ sqrt(p1)
            F_layer += -log(p1)
            verbose && @show -log(p1)
        end
    end
    return Measurement_outcome_boundary(state, sample, F_layer)
end


"""
    reference_bulk_evolution(model::AnyonModel, state::Vector{T}, measure_config::MeasureConfig,
                             samples::Union{Nothing,BitMatrix}=nothing) where{T}

Perform bulk measurement evolution from t₁ to t₂ on quantum state with reference qubits.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `state::Vector{T}`: Initial quantum state vector with reference qubits
- `measure_config::MeasureConfig`: Configuration struct containing τ, t₁, t₂, mode, rng, etc.
- `samples::Union{Nothing,BitMatrix}=nothing`: Predefined measurement sequences for `:sample` mode

# Returns
- `Measurement_outcome_bulk`: A struct containing:
  - `states::Vector{Vector{T}}`: Intermediate states at each time step
  - `samples::BitMatrix`: Measurement outcome sequences
  - `free_energys::Vector{Float64}`: Free energy for each layer

# Notes
- In `:Born` mode, samples are generated probabilistically via Born rule
- In `:sample` mode, `samples` must be provided as input
- The state should already contain reference qubits added via `add_reference_qubits!`
"""
function reference_bulk_evolution(model::AnyonModel, state::Vector{T}, measure_config::MeasureConfig, samples::Union{Nothing,BitMatrix}=nothing) where{T}
    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    basis_F = anyon_basis(model)
    len_F = length(basis_F)
    
    # Deduce k_old from state length
    k_old = round(Int, log2(length(state) ÷ len_F))
    length(state) == (2^k_old * len_F) || error("state length is not compatible with (k_old, N)")
    
    # Build extended basis
    extended_basis = build_extended_basis(k_old, basis_F)
    
    current_state = copy(state)
    if mode == :Born
        return _reference_born_measure(model, current_state, measure_config; extended_basis=extended_basis, k_old=k_old)
    elseif mode == :sample
        return _reference_sample_measure(model, current_state, samples, measure_config; extended_basis=extended_basis, k_old=k_old)
    end
end

"""
    _reference_born_measure(model::AnyonModel, state::Vector{ET}, measure_config::MeasureConfig;
                            extended_basis::Vector{newT}, k_old::Int64=1) where {ET, newT}

Evolve a reference qubit state using probabilistic Born-rule sampling.

This internal helper function is called by `reference_bulk_evolution` when `mode` is `:Born`.

# Arguments
- `model::AnyonModel`: The anyon model.
- `state::Vector{ET}`: The initial state with reference qubits.
- `measure_config::MeasureConfig`: Configuration containing `τ`, `t₁`, `t₂`, `rng`, etc.
- `extended_basis::Vector{newT}`: Extended basis with reference qubits.
- `k_old::Int64=1`: Number of reference qubits in the state.

# Returns
- `Measurement_outcome_bulk`: A struct containing:
  - `states::Vector{Vector{ET}}`: Intermediate states at each full time step.
  - `samples::BitMatrix`: The generated measurement outcome sequences.
  - `free_energys::Vector{Float64}`: The free energy for each measurement layer.
"""
function _reference_born_measure(model::AnyonModel{AT}, current_state::Vector{ET}, measure_config::MeasureConfig; extended_basis::Vector{newT}, k_old::Int64=1) where {AT, ET, newT}
    n_measure = measurement_num(model.anyon_type)*(model.N÷2)
    τ = measure_config.τ
    t₁ = measure_config.t₁
    t₂ = measure_config.t₂
    rng = measure_config.rng
    enable_τ_eff = measure_config.enable_τ_eff
    verbose = measure_config.verbose

    Δt = t₂ - t₁ + 1
    Δt >= 0 || error("t₂ must be >= t₁")
    D = Δt * 2 # number of layers to evolve

    # Initialize sample matrix
    samples = BitMatrix(undef, (D, n_measure))
    sample_free_energy = zeros(D)
    states = Vector{Vector{ET}}(undef, Δt)

    for period in 1:Δt
        τ_eff = (period == Δt && enable_τ_eff) ? τ/2 : τ
        
        outcome1 = _reference_sample_layer(model, τ, current_state, rng, 2*period-1; 
                                           extended_basis=extended_basis, k_old=k_old, verbose=verbose)
        current_state = outcome1.state
        samples[2*period-1, :] = outcome1.sample
        sample_free_energy[2*period-1] = outcome1.free_energy

        outcome2 = _reference_sample_layer(model, τ_eff, current_state, rng, 2*period; 
                                           extended_basis=extended_basis, k_old=k_old, verbose=verbose)
        current_state = outcome2.state
        samples[2*period, :] = outcome2.sample
        sample_free_energy[2*period] = outcome2.free_energy

        states[period] = current_state
    end

    return Measurement_outcome_bulk(states, samples, sample_free_energy)
end

"""
    _reference_sample_measure(model::AnyonModel, state::Vector{ET}, samples::BitMatrix,
                              measure_config::MeasureConfig;
                              extended_basis::Vector{newT}, k_old::Int64=1) where {ET, newT}

Evolve a reference qubit state using a predefined measurement trajectory.

This internal helper function is called by `reference_bulk_evolution` when `mode` is `:sample`.

# Arguments
- `model::AnyonModel`: The anyon model.
- `state::Vector{ET}`: The initial state with reference qubits.
- `samples::BitMatrix`: The predefined measurement outcomes.
- `measure_config::MeasureConfig`: Configuration containing `τ`, `t₁`, `t₂`, etc.
- `extended_basis::Vector{newT}`: Extended basis with reference qubits.
- `k_old::Int64=1`: Number of reference qubits in the state.

# Returns
- `Measurement_outcome_bulk`: A struct containing:
  - `states::Vector{Vector{ET}}`: Intermediate states at each full time step.
  - `samples::BitMatrix`: The input measurement outcome sequences.
  - `free_energys::Vector{Float64}`: The free energy for each measurement layer.
"""
function _reference_sample_measure(model::AnyonModel{AT}, current_state::Vector{ET}, samples::BitMatrix, measure_config::MeasureConfig; extended_basis::Vector{newT}, k_old::Int64=1) where {AT, ET, newT}
    n_measure = measurement_num(model.anyon_type)*(model.N÷2)
    τ = measure_config.τ
    t₁ = measure_config.t₁
    t₂ = measure_config.t₂
    enable_τ_eff = measure_config.enable_τ_eff

    Δt = t₂ - t₁ + 1
    Δt >= 0 || error("t₂ must be >= t₁")
    D = Δt * 2 # number of layers to evolve

    # Validate sample matrix
    isnothing(samples) && error("When mode=:sample, samples must be provided as BitMatrix")
    size(samples) == (D, n_measure) || error("sample size should be ($D, $n_measure)")

    sample_free_energy = zeros(D)
    states = Vector{Vector{ET}}(undef, Δt)

    for period in 1:Δt
        τ_eff = (period == Δt && enable_τ_eff) ? τ/2 : τ

        outcome1 = _reference_apply_measurement_layer(model, τ, current_state, 
                                                      samples[2*period-1, :], 2*period-1;
                                                      extended_basis=extended_basis, k_old=k_old)
        current_state = outcome1.state
        sample_free_energy[2*period-1] = outcome1.free_energy

        outcome2 = _reference_apply_measurement_layer(model, τ_eff, current_state,
                                                      samples[2*period, :], 2*period;
                                                      extended_basis=extended_basis, k_old=k_old)
        current_state = outcome2.state
        sample_free_energy[2*period] = outcome2.free_energy

        states[period] = current_state
    end

    return Measurement_outcome_bulk(states, samples, sample_free_energy)
end


"""
    reference_evolution(N::Int, τ::Float64, forward::Vector{ET}, sample::Matrix{Bool}, 
                        x₂::Int, t₁, t₂; x₁::Int=1, rng=MersenneTwister(), pbc=true, 
                        anyon_type::Symbol=:Fibo, verbose=false, mode::Symbol=:sample)
    reference_evolution(model::AnyonModel, sites, forward::Vector{MPS}, sample::Matrix{Bool},
                        x₂::Int, t₁, t₂, measure_config::MeasureConfig; x₁::Int=1, verbose=false)

Compute temporal/spatial correlation between two spacetime points using cached forward evolution.

This function avoids recomputing the forward evolution up to time `t₁` by using the cached `forward` states.

# Methods

## Vector version (from `ReferenceProbe.jl`)
Compute correlation using state vectors.

### Arguments
- `N::Int`: System size
- `τ::Float64`: Measurement strength parameter
- `forward::Vector{ET}`: Cached forward state evolution trajectory
- `sample::Matrix{Bool}`: Measurement sample configuration
- `x₂::Int`: Site index for second reference qubit insertion
- `t₁::Int`: First time slice index
- `t₂::Int`: Second time slice index (must be >= t₁)
- `x₁::Int=1`: Spatial site index for first reference qubit
- `rng::MersenneTwister=MersenneTwister()`: Random number generator
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: Anyon type (`:Fibo` or `:Ising`)
- `verbose::Bool=false`: Enable verbose output
- `mode::Symbol=:sample`: Evolution mode (`:sample` or `:Born`)

### Returns
- Reference qubit added states between specified spacetime slices.

## MPS version (from `MPSMeasurement.jl`)
Compute correlation using MPS states.

### Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `sites`: ITensor site indices
- `forward::Vector{MPS}`: Cached forward state evolution trajectory from a previous `bulk_evolution` run.
- `sample::Matrix{Bool}`: The full measurement sample configuration for the entire evolution.
- `x₂::Int`: The site index for the second reference qubit insertion.
- `t₁::Int`: The first time slice index for correlation.
- `t₂::Int`: Second time slice index (must be >= t₁)
- `measure_config::MeasureConfig`: Configuration containing τ, rng, mode, etc.
- `x₁::Int=1`: Site index for first reference qubit
- `verbose::Bool=false`: Enable verbose output for debugging.

### Returns
- `Vector{MPS}`: A vector of MPS states for the specified spacetime slices.
- `BitMatrix`: The measurement samples used for the evolution.
- `Vector{Float64}`: The free energy calculated for each layer.
"""
function reference_evolution(model::AnyonModel, forward::Vector{ET}, measure_config::MeasureConfig, sample::BitMatrix) where {ET}
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

    N = model.N
    τ = measure_config.τ
    t₁ = measure_config.t₁
    t₂ = measure_config.t₂
    x₁ = measure_config.x₁
    x₂ = measure_config.x₂
    rng = measure_config.rng
    verbose = measure_config.verbose
    mode = measure_config.mode
    n_measure = measurement_num(model.anyon_type)*(N÷2)
    Δt = size(sample, 1) ÷ 2
    D = size(sample, 1)   # D is the number of layers, while Δt is the true time(# period)

    @assert size(sample, 2) == n_measure "Sample size spatial dimension must be $n_measure, but got $(size(sample, 2))"
    @assert 1 <= t₁ <= t₂ <= D÷2 "Time slice t₁ must before time slice t₂, both must be in the range [1, $(D÷2)]"
    @assert 1 <= x₁ <= x₂ <= N "Site index x₁ must be smaller than site index x₂, both must be in the range [1, $(N)]"
    @assert mode ∈ [:sample, :Born] "mode must be either :sample or :Born, but got $mode"

    # Build model based on anyon_type
    δt = t₂ - t₁ # if δt = 0, pure spatial correlation
    δx = abs(x₂ - x₁) # if δx = 0, pure temporal correlation
    
    # 1) 0 → t₁, the steady state, at the t₁ of forward evolution
    state = forward[t₁]
    statelis = Vector{ET}(undef, Δt) 
    view(statelis, 1:t₁) .= view(forward, 1:t₁)
    sample_layer = BitMatrix(undef, (size(sample, 1), n_measure))
    view(sample_layer, 1:t₁, :) .= Bool.(view(sample, 1:t₁, :))
    sample_free_energy = zeros(Float64, D)

    if δt > 0 && δx > 0 # 3 ref qubits, both spatial and temporal correlation, actually 3-point correlation.
        verbose && @info "t₁ = $(t₁), t₂ = $(t₂), x₁ = $(x₁), x₂ = $(x₂), 3 refs"

        # 2) add reference qubit 1 at x₁, and reference qubit 2 at x₂
        state1 = add_reference_qubits!(model, state, x₁; verbose=verbose)
        state2 = add_reference_qubits!(model, state1, x₂; verbose=verbose)
        
        # 3) t₁ → t₂ evolution, or δt
        config1 = MeasureConfig(τ=τ, t₂=(t₂-t₁), rng=rng, mode=mode, t₁=1, verbose=verbose, enable_τ_eff=false)
        outcome1 = reference_bulk_evolution(model, state2, config1, BitMatrix(sample[2*t₁+1:2*t₂, :]))

        # 4) add reference qubit 3 at x₂
        state3 = add_reference_qubits!(model, outcome1.states[end], x₂; verbose=verbose)

        # 5) t₂ → D evolution
        config2 = MeasureConfig(τ=τ, t₂=(Δt-t₂), rng=rng, mode=mode, t₁=1, verbose=verbose, enable_τ_eff=true)
        outcome2 = reference_bulk_evolution(model, state3, config2, BitMatrix(sample[2*t₂+1:end, :]))

        view(statelis, t₁+1:t₂) .= view(outcome1.states, :)
        view(statelis, t₂+1:Δt) .= view(outcome2.states, :)

        sample_layer[2*t₁+1:2*t₂, :] .= outcome1.samples
        sample_layer[2*t₂+1:end, :] .= outcome2.samples
        sample_free_energy[2*t₁+1:2*t₂] .= outcome1.free_energys
        sample_free_energy[2*t₂+1:end] .= outcome2.free_energys

    elseif δt == 0 # 2 ref qubits, pure 2-point spatial correlation
        verbose && @info "x₁ = $(x₁), x₂ = $(x₂), δx = $(δx), at time slice t₁ = t₂ = $(t₁), 2 refs"
    
        # 2) add reference qubit 1 at x₁, and reference qubit 2 at x₂
        state1 = add_reference_qubits!(model, state, x₁; verbose=verbose)
        state2 = add_reference_qubits!(model, state1, x₂; verbose=verbose)
        
        # 3) t₁ → D evolution
        config2 = MeasureConfig(τ=τ, t₂=(Δt-t₁), rng=rng, mode=mode, t₁=1, verbose=verbose, enable_τ_eff=true)
        outcome2 = reference_bulk_evolution(model, state2, config2, BitMatrix(sample[2*t₁+1:end, :]))

        view(statelis, t₁+1:Δt) .= view(outcome2.states, :)
        sample_layer[2*t₁+1:end, :] .= outcome2.samples
        sample_free_energy[2*t₁+1:end] .= outcome2.free_energys

    elseif δx == 0 # 2 ref qubits, pure 2-point temporal correlation
        verbose && @info "t₁ = $(t₁), t₂ = $(t₂), δt = $(δt), at site x₁ = x₂ = $(x₂), 2 refs"

        # 2) add reference qubit 1 at x₁
        state1 = add_reference_qubits!(model, state, x₂; verbose=verbose)

        # 3) t₁ → t₂ evolution, or δt
        config1 = MeasureConfig(τ=τ, t₂=(t₂-t₁), rng=rng, mode=mode, t₁=1, verbose=verbose, enable_τ_eff=false)
        outcome1 = reference_bulk_evolution(model, state1, config1, BitMatrix(sample[2*t₁+1:2*t₂, :]))

        # 4) add reference qubit 2 at x₂
        state2 = add_reference_qubits!(model, outcome1.states[end], x₂; verbose=verbose)
        
        # 5) t₂ → D evolution
        config2 = MeasureConfig(τ=τ, t₂=(Δt-t₂), rng=rng, mode=mode, t₁=1, verbose=verbose, enable_τ_eff=true)
        outcome2 = reference_bulk_evolution(model, state2, config2, BitMatrix(sample[2*t₂+1:end, :]))

        view(statelis, t₁+1:t₂) .= view(outcome1.states, :)
        view(statelis, t₂+1:Δt) .= view(outcome2.states, :)
        sample_layer[2*t₁+1:2*t₂, :] .= outcome1.samples
        sample_layer[2*t₂+1:end, :] .= outcome2.samples
        sample_free_energy[2*t₁+1:2*t₂] .= outcome1.free_energys
        sample_free_energy[2*t₂+1:end] .= outcome2.free_energys
    end

    return Measurement_outcome_bulk(statelis, sample_layer, sample_free_energy)
end