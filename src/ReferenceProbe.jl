"""
        build_extended_basis(k_total::Int, basis::Vector{T}) where {T<:BitStr}
        build_extended_basis(basis::Vector{ET}, ::Type{T}) where {k_total, ET, T<:BitStr{k_total, Int}}

Construct the extended basis that concatenates `k_total` reference qubits (as a prefix)
with the anyon constraint system basis.

# Arguments
- `k_total::Int`: Total number of reference qubits to include as prefix.
- `basis::Vector{T}`: Sorted anyon basis of the system (suffix).
- `::Type{T}`: BitStr type for the extended basis when using the low-level interface.

# Returns
- `Vector{T}`: Sorted extended basis with reference prefix concatenated to system basis.

# Notes
- Reference bits occupy the highest-order positions (leftmost) in the binary representation.
- Use high-level `build_extended_basis(k_total, basis)` for convenience; low-level variant
    accepts the extended `BitStr` type.
"""
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

function concat_bell_pair(N::Int, current_state::Vector{ET}, extended_basis_old, extended_basis; site_idx::Int64, k_old::Int64=1, new_dim::Int64=2^(N+k_old), verbose =false) where {ET}
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
    add_reference_qubits(model::AnyonModel, state::Vector{ET}, site_idx::Int64; 
                         k_new::Int=1, verbose::Bool=false) -> Vector{ET}

Add reference qubits to a quantum state at specified site using copy method for correlation measurements.

# Principle

Reference qubits are ancillary qubits used to probe spatial and temporal correlations in 
monitored quantum systems. The key idea is to create a maximally entangled state between 
the reference qubit and a system qubit at a specific site, enabling measurement of 
correlations through the reference qubit's reduced density matrix.

## Copy Method

Creates a controlled-copy (CNOT-like) entanglement:
- If the system qubit at `site_idx` is in state |0⟩, the reference qubit becomes |0⟩
- If the system qubit at `site_idx` is in state |1⟩, the reference qubit becomes |1⟩

Mathematically, for a state |ψ⟩ = α|0⟩ + β|1⟩ at `site_idx`:
```
    |ψ⟩ → α|00⟩ + β|11⟩  (system ⊗ reference)
    |ψ⟩ → α|A0⟩ + β|B1⟩  (system ⊗ reference), where A, B are orthogonal states
```

This preserves the original state's structure while creating classical correlation.

For reset method (measure system qubit first, then create Bell pair), use 
[`add_reference_qubits_reset`](@ref) instead.

## Basis Structure

The extended basis combines reference qubits (prefix) with the Fibonacci/Ising 
constraint basis (suffix):
```
    |extended_basis⟩ = |ref₁ ref₂ ... refₖ⟩ ⊗ |system basis⟩
```

Reference qubits are stored as the leftmost (highest) bits in the binary representation.

# Arguments
- `model::AnyonModel`: Anyon model containing system size N, boundary conditions, etc.
- `state::Vector{ET}`: Quantum state vector (may already contain k_old reference qubits)
- `site_idx::Int64`: Site index for reference qubit insertion (1 ≤ site_idx ≤ N)
- `k_new::Int=1`: Number of new reference qubits to add (must be 0 or 1)
- `verbose::Bool=false`: Enable verbose output for debugging

# Returns
- `Vector{ET}`: New normalized state with reference qubit added

# Notes
- Each system qubit can only be entangled with one reference qubit
- To add multiple reference qubits at different sites, call this function multiple times
- The number of existing reference qubits (k_old) is automatically deduced from state length
- State dimension grows as: new_dim = 2^(k_old + k_new) × len(Fibonacci_basis)

# Example
```julia
# Add first reference qubit at site 1
state1 = add_reference_qubits(model, initial_state, 1)

# Add second reference qubit at site N÷2
state2 = add_reference_qubits(model, state1, N÷2)

# Now state2 has 2 reference qubits for 2-point correlation measurement
```

See also: [`add_reference_qubits_reset`](@ref), [`reference_rdm`](@ref), [`reference_evolution`](@ref), [`build_extended_basis`](@ref)
"""
function add_reference_qubits(model::AnyonModel, state::Vector{ET}, site_idx::Int64; k_new::Int=1, verbose::Bool=false)::Vector{ET} where {ET}
    # Add k_new reference qubits to the state at the specified site_idx using :copy mode.
    # For :reset mode, use add_reference_qubits_reset instead.
    N = model.N
    @assert 1 <= site_idx <= N "Site index must be in the range [1, $(N)]"
    1 >= k_new >= 0 || error("k_new must be in [0,1]")
    # Because each qubit can only concat with one reference qubit, so k_new can only be 0 or 1. If need to add more reference qubits, use add_reference_qubits multiple times at different site.
    basis_F = anyon_basis(model)
    len_F   = length(basis_F)
        

    # old reference qubit number, k_old
    k_old = round(Int, log2(length(state) ÷ len_F))
    length(state) == (2^k_old * len_F) ||
        error("state length is not compatible with ($(k_old), $(N)), can not deduce k_old from state length")
    @assert k_old >= 0 "k_old must be non-negative, but got $(k_old)"
    
    # new reference prefix
    k_total = k_old + k_new
    new_dim = 2^k_total * len_F
    extended_basis = build_extended_basis(k_total, basis_F)
    extended_basis_old = build_extended_basis(k_old, basis_F)

    # :copy mode - Create classical correlation via copy operation
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
end

"""
    add_reference_qubits_reset(model::AnyonModel, state::Vector{ET}, site_idx::Int64; 
                                k_new::Int=1, verbose::Bool=false) -> Tuple{Float64, Float64, Vector{ET}, Vector{ET}}

Add reference qubits using reset method: first measures/resets the system qubit, then creates a Bell pair.

# Principle
1. Project the system qubit onto |0⟩ or |1⟩ basis (strong measurement)
2. Create a Bell pair between the reset qubit and reference qubit
3. Returns both branches with their probabilities

This method is useful when you want to prepare a known initial state before 
adding the reference qubit.

# Arguments
- `model::AnyonModel`: Anyon model containing system size N, boundary conditions, etc.
- `state::Vector{ET}`: Quantum state vector (may already contain k_old reference qubits)
- `site_idx::Int64`: Site index for reference qubit insertion (1 ≤ site_idx ≤ N)
- `k_new::Int=1`: Number of new reference qubits to add (must be 0 or 1)
- `verbose::Bool=false`: Enable verbose output for debugging

# Returns
- `Tuple{Float64, Float64, Vector{ET}, Vector{ET}}`: 
  (prob₀, prob₁, state_after_0, state_after_1) for both measurement branches

See also: [`add_reference_qubits`](@ref), [`reference_rdm`](@ref), [`reference_evolution`](@ref)
"""
function add_reference_qubits_reset(model::AnyonModel, state::Vector{ET}, site_idx::Int64; k_new::Int=1, verbose::Bool=false)::Tuple{Float64, Float64, Vector{ET}, Vector{ET}} where {ET}
    # Reset mode: measure/reset qubit, then concat with new reference qubit.
    N = model.N
    @assert 1 <= site_idx <= N "Site index must be in the range [1, $(N)]"
    1 >= k_new >= 0 || error("k_new must be in [0,1]")
    
    basis_F = anyon_basis(model)
    len_F   = length(basis_F)
    
    # old reference qubit number, k_old
    k_old = round(Int, log2(length(state) ÷ len_F))
    length(state) == (2^k_old * len_F) ||
        error("state length is not compatible with ($(k_old), $(N)), can not deduce k_old from state length")
    @assert k_old >= 0 "k_old must be non-negative, but got $(k_old)"
    
    # new reference prefix
    k_total = k_old + k_new
    new_dim = 2^k_total * len_F
    extended_basis = build_extended_basis(k_total, basis_F)
    extended_basis_old = build_extended_basis(k_old, basis_F)

    verbose && @info "Using reset way to add reference qubit"
    # NEED to do RESET qubit to 0 or +!!! then concat with the new reference qubit. 
  
    reset_model = AnyonModel(model.anyon_type, model.N, measure_operator=reset_type(model), pbc=model.pbc)
    state_after_0 = reference_measuremap(reset_model, 1000.0, state, site_idx, false, extended_basis=extended_basis_old, k_old=k_old)
    state_after_1 = reference_measuremap(reset_model, 1000.0, state, site_idx, true, extended_basis=extended_basis_old, k_old=k_old)

    prob_sqrt0 = state_after_0' * state_after_0
    prob_sqrt1 = 1 - prob_sqrt0

    verbose && @info "Reset probabilities: prob_0 = $(prob_sqrt0), prob_1 = $(prob_sqrt1)"

    current_statep = state_after_0 ./ sqrt(prob_sqrt0)
    current_statem = state_after_1 ./ sqrt(prob_sqrt1)

    new_statep = concat_bell_pair(N, current_statep, extended_basis_old, extended_basis, site_idx = site_idx, k_old = k_old, new_dim = new_dim)
    new_statem = concat_bell_pair(N, current_statem, extended_basis_old, extended_basis, site_idx = site_idx, k_old = k_old, new_dim = new_dim)

    return (prob_sqrt0, prob_sqrt1, new_statep, new_statem)
end

reset_type(model::AnyonModel{FibonacciAnyon}) = :reset
reset_type(model::AnyonModel{IsingAnyon}) = (model.measure_operator == :X) ? :X : :reset

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
    
    N = model.N
    length(state) == (2^k_old * len_F) || error("state length is not compatible with ($(k_old), $(N)), can not deduce k_old from state length")
    @assert 2^k_old*length(unsorted_basis) == length(state) "state length is expected to be $(2^k_old*length(unsorted_basis)), but got $(length(state))"

    if traceref
        # If traceref is true, we need to trace out the reference qubit. otherwise, we trace out system.
        ref_model = AnyonModel(IsingAnyon(), k_old, pbc=false)
        return disjoint_rdm(ref_model, model, subsystems, Int[], state;)
    else
        ref_model = AnyonModel(IsingAnyon(), k_old, pbc=false)
        totalsubBpbc = (length(subsystems) == model.N) ? true : false
        return disjoint_rdm(ref_model, model, Int[], subsystems, state; pbcB = totalsubBpbc)
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

    τ = measure_config.τ
    τ_eff = measure_config.enable_τ_eff ? τ/2 : τ
    verbose = measure_config.verbose
    if mode == :sample
        N = model.N
        size(sample, 1) == _samples_per_layer(model) || error("sample size mismatch with anyon_model $(N)")
        verbose && @info "Using given sample evolution mode"
        _reference_apply_measurement_layer(model, τ_eff, state, sample, layer_idx; extended_basis=extended_basis, k_old=k_old)
    elseif mode == :Born
        verbose && @info "Using Born rule driven evolution mode"
        _reference_sample_layer(model, τ_eff, state, measure_config.rng, layer_idx; extended_basis=extended_basis, k_old=k_old)
    end
end


"""
    _reference_apply_measurement_layer(model::AnyonModel, τ::Float64, state::Vector{ET},
                                       layer_sample::BitVector, layer_idx::Int64;
                                       extended_basis::Vector{newT}, k_old::Int64=1) where {ET, newT}

Apply deterministic measurements to a layer with given outcomes for states with reference qubits.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Quantum state vector with reference qubits
- `layer_sample::BitVector`: Measurement outcomes for the layer
- `layer_idx::Int64`: Layer index (1-based) to determine measurement pattern
- `extended_basis::Vector{newT}`: Extended basis with reference qubits
- `k_old::Int64=1`: Number of reference qubits in the state

# Returns
- `Measurement_outcome_boundary`: A struct containing the post-measurement state, sample, and total free energy.
"""
function _reference_apply_measurement_layer(model::AnyonModel, τ::Float64, state::Vector{ET},
    layer_sample::BitVector, layer_idx::Int64;
    extended_basis::Vector{newT}, k_old::Int64=1) where {ET, newT}

    total_free_energy = zero(real(ET))
    measurement_sites, measure_anyon_model = _obtain_measurement_config(model, layer_idx)  

    for (idx, sign) in enumerate(layer_sample)
        state = reference_measuremap(measure_anyon_model, τ, state, measurement_sites[idx], sign, extended_basis=extended_basis, k_old=k_old)

        prob = real(dot(state, state))
        total_free_energy += -log(prob)
        state ./= sqrt(prob)
    end

    return Measurement_outcome_boundary(state, layer_sample, Float32(total_free_energy))
end


"""
    _reference_sample_layer(model::AnyonModel, τ_eff::Float64, state::Vector{T}, 
                            rng::MersenneTwister=MersenneTwister(), layer_idx::Int64=1;
                            extended_basis::Vector{newT}, k_old::Int64=1, 
                            verbose::Bool=false) where {T, newT}

Perform random measurement on a layer using Born rule sampling for states with reference qubits.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ_eff::Float64`: Effective measurement strength parameter
- `state::Vector{T}`: Quantum state vector with reference qubits
- `rng::MersenneTwister=MersenneTwister()`: Random number generator
- `layer_idx::Int64=1`: Layer index (1-based) to determine measurement pattern
- `extended_basis::Vector{newT}`: Extended basis with reference qubits
- `k_old::Int64=1`: Number of reference qubits in the state
- `verbose::Bool=false`: Whether to print debug information

# Returns
- `Measurement_outcome_boundary`: A struct containing the post-measurement state, sample outcomes, and free energy.
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
    return Measurement_outcome_boundary(state, sample, Float32(F_layer))
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
- The state should already contain reference qubits added via `add_reference_qubits`
"""
function reference_bulk_evolution(model::AnyonModel, state::Vector{T}, measure_config::MeasureConfig, samples::Union{Nothing,BitMatrix}=nothing) where{T}
    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    basis_F = anyon_basis(model)
    len_F = length(basis_F)
    
    # Deduce k_old from state length
    k_old = round(Int, log2(length(state) ÷ len_F))
    length(state) == (2^k_old * len_F) || error("state length is not compatible with ($(k_old), $(N))")
    
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
    n_measure = _samples_per_layer(model)
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
    sample_free_energy = zeros(Float32, D)
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
    n_measure = _samples_per_layer(model)
    τ = measure_config.τ
    t₁ = measure_config.t₁
    t₂ = measure_config.t₂
    enable_τ_eff = measure_config.enable_τ_eff

    Δt = t₂ - t₁ + 1
    # Debug: verbose && @show size(samples), Δt, t₁, t₂ ,Δt * 2
    Δt >= 0 || error("t₂ must be >= t₁")
    D = Δt * 2 # number of layers to evolve

    # Validate sample matrix
    isnothing(samples) && error("When mode=:sample, samples must be provided as BitMatrix")
    size(samples) == (D, n_measure) || error("sample size should be ($D, $n_measure)")

    sample_free_energy = zeros(Float32, D)
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
    reference_evolution(model::AnyonModel, forward::Vector{ET}, measure_config::MeasureConfig, 
                        sample::Union{BitMatrix, Nothing}=nothing) -> Measurement_outcome_bulk

Compute temporal/spatial correlation between two spacetime points using cached forward evolution.

This function avoids recomputing the forward evolution up to time `t₁` by using the cached 
`forward` states, then adds reference qubits and continues evolution to measure correlations.

# Spacetime Diagram
```
    | ---->|
    forward
          Ref1                    x₁ = 1            
           |
           |
          Ref2 --> Ref3               x₂ = N/2
    | ----> |______| -----> |
    1      t₁   t₂=t₁+δt   Δt
    | --------------------->|
    sample
```

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, anyon_type, etc.)
- `forward::Vector{ET}`: Cached forward state evolution trajectory from a previous `bulk_evolution` run
- `measure_config::MeasureConfig`: Configuration struct containing:
  - `τ::Float64`: Measurement strength parameter
  - `t₁::Int`: First time slice index for correlation
  - `t₂::Int`: Second time slice index (must be >= t₁)
  - `x₁::Int`: Site index for first reference qubit
  - `x₂::Int`: Site index for second reference qubit (must be >= x₁)
  - `mode::Symbol`: Evolution mode (`:sample` or `:Born`)
  - `rng::MersenneTwister`: Random number generator (for `:Born` mode)
  - `verbose::Bool`: Enable verbose output
- `sample::Union{BitMatrix, Nothing}=nothing`: Full measurement sample configuration for the entire evolution

# Returns
- `Measurement_outcome_bulk`: A struct containing:
  - `states::Vector{ET}`: State evolution trajectory with reference qubits
  - `samples::BitMatrix`: Measurement samples used for the evolution
  - `free_energys::Vector{Float64}`: Free energy calculated for each layer

# Correlation Types
Based on `δt = t₂ - t₁` and `δx = |x₂ - x₁|`:
- `δt > 0 && δx > 0`: 3-point correlation (both spatial and temporal), uses 3 reference qubits
- `δt == 0`: Pure 2-point spatial correlation at fixed time, uses 2 reference qubits
- `δx == 0`: Pure 2-point temporal correlation at fixed site, uses 2 reference qubits

# Notes
- The `forward` states should come from a previous `bulk_evolution` call
- Reference qubits are added via [`add_reference_qubits`](@ref) at the specified spacetime points
- In `:sample` mode, `sample` must be provided; in `:Born` mode, samples are generated via Born rule

See also: [`add_reference_qubits`](@ref), [`reference_bulk_evolution`](@ref), [`ref_correlation`](@ref)
"""
function reference_evolution(model::AnyonModel, forward::Vector{ET}, measure_config::MeasureConfig, sample::BitMatrix) where {ET}
    # This function is used to evolve state to compute the spatial-temporal correlation at different space and time slices cache, avoiding the repeated calculation of the state evolution. INPUT the forward state evolution, given the trajectory.
    # The length of forward is to t₁, but the time length of sample is 0-> D, need to note that the input sample will be updated in :Born mode, and the smaple front part must be consistent with forward.
    # t₁ and t₂ are the indices of the time slices in the sample. Spacetime is:
    # | ---->|
    # forward
    #       Ref1                    x₁ = 1            
    #        |
    #        |
    #       Ref2 --> Ref3               x₂ = N/2
    # | ----> |______| -----> |
    # 1      t₁   t₂=t₁+δt   Δt
    # | --------------------->|
    # sample

    N = model.N
    τ = measure_config.τ
    t₁ = measure_config.t₁
    t₂ = measure_config.t₂
    x₁ = measure_config.x₁
    x₂ = measure_config.x₂
    rng = measure_config.rng
    verbose = measure_config.verbose
    mode = measure_config.mode
    n_measure = _samples_per_layer(model)
    Δt = size(sample, 1) ÷ 2
    D = size(sample, 1)   # D is the number of layers, while Δt is the true time(# period), each time slice have two layers.

    @assert size(sample, 2) == n_measure "Sample size spatial dimension must be $n_measure, but got $(size(sample, 2))"
    @assert 1 <= t₁ <= t₂ <= Δt "Time slice t₁ must before time slice t₂, both must be in the range [1, $(D÷2)]"
    @assert 1 <= x₁ <= x₂ <= N "Site index x₁ must be smaller than site index x₂, both must be in the range [1, $(N)]"
    @assert mode ∈ [:sample, :Born] "mode must be one of :sample, :Born" # It will only evolve according to the mode, even if sample is given.

    # Build model based on anyon_type
    δt = t₂ - t₁ # if δt = 0, pure spatial correlation
    δx = abs(x₂ - x₁) # if δx = 0, pure temporal correlation
    
    # 1) 1 → t₁, the steady state, at the t₁ of forward evolution
    state = forward[t₁]
    statelis = Vector{ET}(undef, Δt) 
    view(statelis, 1:t₁) .= view(forward, 1:t₁)
    sample_layer = BitMatrix(undef, (D, n_measure))
    view(sample_layer, 1:2t₁, :) .= view(sample, 1:2t₁, :)
    sample_free_energy = zeros(Float32, D)

    if δt > 0 && δx > 0 # 3 ref qubits, both spatial and temporal correlation, actually 3-point correlation.
        verbose && @info "t₁ = $(t₁), t₂ = $(t₂), x₁ = $(x₁), x₂ = $(x₂), 3 refs"

        # 2) add reference qubit 1 at x₁, and reference qubit 2 at x₂
        state1 = add_reference_qubits(model, state, x₁; verbose=verbose)
        state2 = add_reference_qubits(model, state1, x₂; verbose=verbose)
        
        # 3) t₁ → t₂ evolution, or δt
        config1 = MeasureConfig(τ=τ, t₂= t₁+δt, rng=rng, mode=mode, t₁=t₁+1, verbose=verbose, enable_τ_eff=false)
        outcome1 = reference_bulk_evolution(model, state2, config1, sample[2*t₁+1:2*t₂, :])

        # 4) add reference qubit 3 at x₂
        state3 = add_reference_qubits(model, outcome1.states[end], x₂; verbose=verbose)

        # 5) t₂ → D evolution
        config2 = MeasureConfig(τ=τ, t₂=Δt, rng=rng, mode=mode, t₁=t₂+1, verbose=verbose, enable_τ_eff=true)
        outcome2 = reference_bulk_evolution(model, state3, config2, sample[2*t₂+1:end, :])

        view(statelis, t₁+1:t₂) .= view(outcome1.states, :)
        view(statelis, t₂+1:Δt) .= view(outcome2.states, :)

        sample_layer[2*t₁+1:2*t₂, :] .= outcome1.samples
        sample_layer[2*t₂+1:end, :] .= outcome2.samples
        sample_free_energy[2*t₁+1:2*t₂] .= outcome1.free_energys
        sample_free_energy[2*t₂+1:end] .= outcome2.free_energys

    elseif δt == 0 # 2 ref qubits, pure 2-point spatial correlation
        verbose && @info "x₁ = $(x₁), x₂ = $(x₂), δx = $(δx), at time slice t₁ = t₂ = $(t₁), 2 refs"
    
        # 2) add reference qubit 1 at x₁, and reference qubit 2 at x₂
        state1 = add_reference_qubits(model, state, x₁; verbose=verbose)
        state2 = add_reference_qubits(model, state1, x₂; verbose=verbose)
        
        # 3) t₁ → D evolution
        config2 = MeasureConfig(τ=τ, t₂= Δt, rng=rng, mode=mode, t₁=t₁+1, verbose=verbose, enable_τ_eff=true)
        # Debug: verbose && @show size(sample),  Δt, t₁, t₂, size(sample[2*t₁+1:end, :])
        outcome2 = reference_bulk_evolution(model, state2, config2, sample[2*t₁+1:end, :])

        view(statelis, t₁+1:Δt) .= view(outcome2.states, :)
        sample_layer[2*t₁+1:end, :] .= outcome2.samples
        sample_free_energy[2*t₁+1:end] .= outcome2.free_energys

    elseif δx == 0 # 2 ref qubits, pure 2-point temporal correlation
        verbose && @info "t₁ = $(t₁), t₂ = $(t₂), δt = $(δt), at site x₁ = x₂ = $(x₂), 2 refs"

        # 2) add reference qubit 1 at x₂ at time slice t₁
        state1 = add_reference_qubits(model, state, x₂; verbose=verbose)

        # 3) t₁ → t₂ evolution, or δt
        config1 = MeasureConfig(τ=τ, t₂= t₁ +δt, rng=rng, mode=mode, t₁=t₁+1, verbose=verbose, enable_τ_eff=false)
        outcome1 = reference_bulk_evolution(model, state1, config1, sample[2*t₁+1:2*t₂, :])

        # 4) add reference qubit 2 at x₂ at time slice t₂
        state2 = add_reference_qubits(model, outcome1.states[end], x₂; verbose=verbose)
        
        # 5) t₂ → D evolution
        config2 = MeasureConfig(τ=τ, t₂=Δt, rng=rng, mode=mode, t₁=t₂+1, verbose=verbose, enable_τ_eff=true)
        outcome2 = reference_bulk_evolution(model, state2, config2, sample[2*t₂+1:end, :])

        view(statelis, t₁+1:t₂) .= view(outcome1.states, :)
        view(statelis, t₂+1:Δt) .= view(outcome2.states, :)
        sample_layer[2*t₁+1:2*t₂, :] .= outcome1.samples
        sample_layer[2*t₂+1:end, :] .= outcome2.samples
        sample_free_energy[2*t₁+1:2*t₂] .= outcome1.free_energys
        sample_free_energy[2*t₂+1:end] .= outcome2.free_energys
    end

    return Measurement_outcome_bulk(statelis, sample_layer, sample_free_energy)
end