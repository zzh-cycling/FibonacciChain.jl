"""
    measure_basismap(model::AnyonModel, τ::Float64, state::T, i::Int, sign::Bool) where {T}

Map single basis state under measurement operation at site i.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, measure_operator)
- `τ::Float64`: Measurement strength parameter
- `state::T`: Input basis state
- `i::Int`: Measurement site index (1 ≤ i ≤ N)
- `sign::Bool`: Measurement outcome (false for +, true for -)

# Returns
- `Tuple`: Either `(basis1, basis2, coeff1, coeff2)` for superposition output or `(basis, coeff)` for single output

Maps individual basis states according to measurement protocols and fusion rules.

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; T = BitStr{N, Int};

julia> model = AnyonModel(FibonacciAnyon(), N; pbc=true, measure_operator=:Antiferro);

julia> state = T(0b0100);  # Single τ at site 2

julia> τ = 1.0;  # Measurement strength

julia> result = measure_basismap(model, τ, state, 2, false);

julia> length(result) ∈ [2, 4]  # Returns 2 or 4 elements depending on configuration
true
```
"""
function measure_basismap(model::AnyonModel{AT}, τ::Float64, state::T, i::Int, sign::Bool) where {T, AT<:AbstractAnyonType}
    # default for PBC system, map basis (not state!!!), and index count from the left.
    @assert 1 <= i <= model.N "Index i must be in the range [1, N]"
    @assert sign in (0, 1) "sign must be either 0 the plus, 1 the minus"

    return _apply_result(model, τ, state, i, sign)
end

function _apply_result(model::AnyonModel{FibonacciAnyon}, τ::Float64, state::T, i::Int, sign::Bool) where {T}
    @assert model.measure_operator ∈ [:Ferro, :Antiferro] "measure_operator must be :Ferro or :Antiferro"
    N = model.N
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))
    ϕ = (1+√5)/2

    if τ >= 1e2
        # true is 1, false is 0
        cstτ = 0.5
        coef = sign ? -0.5 : 0.5 
    else
        cstτ = (exp(τ) + 1) / (2 * √(exp(2τ) + 1))
        coef = sign ? (1 - exp(τ)) / (2 * √(exp(2τ) + 1)) : (exp(τ) - 1) / (2 * √(exp(2τ) + 1)) 
    end
    
    if 2<= i <= N-1
        mask=bmask(T,1,2,3) << (N-i-1)
        str100, str101, str010, str001, str000 = T(4) << (N-i-1), T(5) << (N-i-1), T(2) << (N-i-1), T(1) << (N-i-1), T(0) << (N-i-1)
        if state & mask == str000
            return state, X(state,i), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2)
        elseif state & mask == str010
            return state, X(state,i), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2)
        elseif state & mask == str001
            return state, state, cstτ+coef, 0
        elseif state & mask == str100
            return state, state, cstτ+coef, 0
        elseif state & mask == str101
            return state, state, cstτ-coef, 0
        end
    end

    if model.pbc
        if i == 1 #count from the left
        mask=bmask(T, N, N-1,1)
        str100, str101, str010, str001, str000 = bmask(T,1), bmask(T, N-1, 1), bmask(T, N), bmask(T, N-1), T(0)
            if state & mask == str000
                return state, X(state,i), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2)
            elseif state & mask == str010
                return state, X(state,i), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2)
            elseif state & mask == str001
                return state, state, cstτ+coef, 0
            elseif state & mask == str100
                return state, state, cstτ+coef, 0
            elseif state & mask == str101
                return state, state, cstτ-coef, 0
            end
        elseif i == N #count from the left
        mask=bmask(T, N, 2, 1)
        str100, str101, str010, str001, str000 = bmask(T,2), bmask(T, N, 2), bmask(T, 1), bmask(T, N), T(0)
            if state & mask == str000
                return state, X(state,i), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2)
            elseif state & mask == str010
                return state, X(state,i), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2)
            elseif state & mask == str001
                return state, state, cstτ+coef, 0
            elseif state & mask == str100
                return state, state, cstτ+coef, 0
            elseif state & mask == str101
                return state, state, cstτ-coef, 0
            end
        end
    end
end

function _apply_result(model::AnyonModel{IsingAnyon}, τ::Float64, state::T, i::Int, sign::Bool) where {T}
    @assert model.measure_operator in [:X, :ZZ] "measure_operator must be either :X or :ZZ"

    N = model.N
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))

    if model.measure_operator == :X
        if τ >= 1e2
            cstτ = 0.5
            coef = sign ? -0.5 : 0.5 
        else
            cstτ = cosh(τ/2) / √(2cosh(τ))
            coef = sign ? -sinh(τ/2) / √(2cosh(τ)) : sinh(τ/2) / √(2cosh(τ)) 
        end

        return state, X(state,i), cstτ, coef

    elseif model.measure_operator == :ZZ
        if τ >= 1e2
            cstτ = 0.5
            coef = sign ? -0.5 : 0.5
        else
            cstτ = cosh(τ/2) / √(2cosh(τ))
            coef = sign ? -sinh(τ/2) / √(2cosh(τ)) : sinh(τ/2) / √(2cosh(τ))
        end

        if 1<= i <= N-1
            if ((state >> (N - i)) & 1) == ((state >> (N - i -1)) & 1)
                return state, state, cstτ+coef, 0
            else
                return state, state, cstτ-coef, 0
            end
        end

        if pbc && i == N
            if (state & 1) == (state >> (N-1) & 1)
                return state, state, cstτ+coef, 0
            else
                return state, state, cstτ-coef, 0
            end
        end
    end
end

# function _apply_result(anyontype::IsingAnyon, τ::Float64, state::T, i::Int, sign::Int64) where {T}
#     For :reset, :resetFibo, :IsingZ, not sure if needed yet.
#     if τ >= 1e2
#         cstτ = 0.5
#         coef = sign ? 0.5 : -0.5
#     else
#         cstτ = cosh(τ/2) / √(2cosh(τ))
#         coef = sign ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
#     end

#     return state, state, (state[N - i + 1] == 0) ? cstτ + coef : cstτ - coef, 0
# end


function measure_matrix(model::AnyonModel{AT}, τ::Float64, idx::Int, sign::Bool) where {AT<:AbstractAnyonType}

    if model.measure_operator ∈ [:Ferro, :Antiferro]
        @assert model.pbc || (2 <= idx <= model.N-1) "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif model.measure_operator == :IsingZZ
        @assert model.pbc || (1 <= idx <= model.N-1) "Index idx must be in [1, N-1] for open BC (IsingZZ)"
    elseif model.measure_operator ∈ (:IsingX, :IsingZ, :reset, :resetFibo)
        @assert model.pbc || (1 <= idx <= model.N) "Index idx must be in [1, N] for open BC (IsingX)"
    else
        error("Unknown measure class: $(model.anyon_type)")
    end

    basis = anyon_basis(model)
    l = length(basis)
    Bmatrix = zeros(l, l)

    for i in 1:l
        s1, s2, w1, w2 = measure_basismap(model, τ, basis[i], idx, sign)

        if w2 == 0
            Bmatrix[i, i] += w1
        else
            j2 = searchsortedfirst(basis, s2)
            Bmatrix[i, i] += w1
            Bmatrix[i, j2] += w2
        end
    end

    return Bmatrix
end

function measuremap(model::AnyonModel{AT}, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool) where {ET, AT<:AbstractAnyonType}
    # input a superposition state, and output the measured state (tedancy fusion to 0 or 1 in Fibonacci measure class, or X ZZ in Ising measure class)
    if model.measure_operator ∈ [:Ferro, :Antiferro]
        @assert model.pbc || (2 <= idx <= model.N-1) "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif model.measure_operator == :IsingZZ
        @assert model.pbc || (1 <= idx <= model.N-1) "Index idx must be in [1, N-1] for open BC (IsingZZ)"
    elseif model.measure_operator ∈ (:IsingX, :IsingZ, :reset, :resetFibo)
        @assert model.pbc || (1 <= idx <= model.N) "Index idx must be in [1, N] for open BC (IsingX)"
    else
        error("Unknown measure class: $(model.anyon_type)")
    end
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    basis=anyon_basis(model)
    l=length(basis)
    @assert l == length(state) "state length is expected to be $(l), but got $(length(state))"
    mapped_state = zeros(ET, length(state))
    for i in 1:l
        outputstate1, outputstate2, output1, output2 = measure_basismap(model, τ, basis[i], idx, sign)
        
        if output2 == 0
            mapped_state[i]+=output1*state[i] # outputstate is the same as basis[i]
        else
            j2=searchsortedfirst(basis, outputstate2)
            mapped_state[i]+=output1*state[i] # outputstate1 is the same as basis[i]
            mapped_state[j2]+=output2*state[i]
        end
    end
    
    return mapped_state
end

function laddermeasuremap(model::AnyonModel{AT}, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool) where {ET, AT<:AbstractAnyonType}
    # input a superposition state, and output the braided state
    @assert model.pbc || (2 <= idx <= model.N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    basis=anyon_basis(model)
    l=length(basis)
    @assert l^2 == length(state) "state length is expected to be $(l^2), but got $(length(state))"
    mapped_state = zeros(ET, length(state))
    for i in 1:l
        for j in 1:l
            output1 = measure_basismap(model, τ, basis[i], idx, sign)
            output2 = measure_basismap(model, τ, basis[j], idx, sign)
            if length(output1) == 4 && length(output2) == 4
                basisi1, basisi2, coefi1, coefi2=output1
                basisj1, basisj2, coefj1, coefj2=output2
                i2=searchsortedfirst(basis, basisi2)
                j2=searchsortedfirst(basis, basisj2)
                # Here noting that the state is a vertorizing density matrix, so the index is i+(j-1)*len, not state[i], state[j]
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi1*coefj1
                mapped_state[(i-1)*l+j2]+=state[(i-1)*l+j]*coefi1*coefj2
                mapped_state[(i2-1)*l+j]+=state[(i-1)*l+j]*coefi2*coefj1
                mapped_state[(i2-1)*l+j2]+=state[(i-1)*l+j]*coefi2*coefj2
            elseif length(output1) == 4 && length(output2) == 2
                basisi1, basisi2, coefi1, coefi2=output1
                basisj, coefj=output2
                i2=searchsortedfirst(basis, basisi2)  
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi1*coefj
                mapped_state[(i2-1)*l+j]+=state[(i-1)*l+j]*coefi2*coefj
            elseif length(output1) == 2 && length(output2) == 4
                basisi, coefi=output1
                basisj1, basisj2, coefj1, coefj2=output2
                j2=searchsortedfirst(basis, basisj2)
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi*coefj1
                mapped_state[(i-1)*l+j2]+=state[(i-1)*l+j]*coefi*coefj2
            else
                basisi, coefi=output1
                basisj, coefj=output2
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi*coefj
            end
        end
    end
    
    return mapped_state
end

"""
    measurement_enumeration(model::AnyonModel, τ::Float64, initial_state::Vector{ET}, measurement_sites::Vector{Int}) where {ET}

Enumerate all trajectories of measurements on a given initial state.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `initial_state::Vector{ET}`: Initial quantum state vector
- `measurement_sites::Vector{Int}`: Sites to measure

# Returns
- `Tuple{Vector{Vector{ET}}, Vector{Vector{Int64}}, Vector{Float64}}`:
  (final_states, trajectories, probabilities)

# Examples
```jldoctest
julia> using FibonacciChain, LinearAlgebra

julia> model = AnyonModel(FibonacciAnyon(), 4; pbc=true, measure_operator=:Antiferro);

julia> basis = anyon_basis(model);

julia> state = normalize(ones(length(basis)));

julia> states, trajs, probs = measurement_enumeration(model, 1.0, state, [2]);

julia> length(states) == 2  # Two possible outcomes
true
```
"""
function measurement_enumeration(model::AnyonModel{AT}, τ::Float64, initial_state::Vector{ET}, measurement_sites::Vector{Int}) where {ET, AT<:AbstractAnyonType}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"
    
    # Initialize, only one initial state
    current_level_states = [copy(initial_state)]
    current_level_trajectories = [Bool[]]
    current_level_probabilities = [1.0]
    
    for (measurement_idx, site) in enumerate(measurement_sites)
        next_level_states = Vector{Vector{ET}}()
        next_level_trajectories = Vector{Vector{Int64}}()
        next_level_probabilities = Vector{Float64}()
        
        # Branching for each current state
        for (state_idx, state) in enumerate(current_level_states)
            current_trajectory = current_level_trajectories[state_idx]
            current_prob = current_level_probabilities[state_idx]

            state_after_p = measuremap(model, τ, state, site, false)
            prob_p = real(dot(state_after_p, state_after_p))

            normalized_state_p = state_after_p / sqrt(prob_p)
            new_trajectory_p = [current_trajectory; 0]
            new_prob_p = current_prob * prob_p
            
            push!(next_level_states, normalized_state_p)
            push!(next_level_trajectories, new_trajectory_p)
            push!(next_level_probabilities, new_prob_p)

            state_after_m = measuremap(model, τ, state, site, true)
            prob_m = real(dot(state_after_m, state_after_m))

            normalized_state_m = state_after_m / sqrt(prob_m)
            new_trajectory_m = [current_trajectory; 1]
            new_prob_m = current_prob * prob_m
            
            push!(next_level_states, normalized_state_m)
            push!(next_level_trajectories, new_trajectory_m)
            push!(next_level_probabilities, new_prob_m)

        end
        
        current_level_states = next_level_states
        current_level_trajectories = next_level_trajectories
        current_level_probabilities = next_level_probabilities

    end
    

    return current_level_states, current_level_trajectories, current_level_probabilities
end


function measurement_tree_visualization(trajectories::Vector{Vector{Int64}}, probabilities::Vector{Float64})
    """
    visualize measurement tree
    """
    total_prob = sum(probabilities)
    normalized_probs = probabilities / total_prob
    
    println("Measurement Tree Visualization:")
    println("==============================")
    
    max_length = maximum(length.(trajectories))
    
    for traj_length in 0:max_length
        level_indices = findall(t -> length(t) == traj_length, trajectories)
        if !isempty(level_indices)
            println("Level $(traj_length):")
            for idx in level_indices
                traj = trajectories[idx]
                prob = normalized_probs[idx]
                indent = "  " ^ traj_length
                if traj_length == 0
                    println("$(indent)Initial state (prob: 1.0)")
                else
                    println("$(indent)$(traj) (prob: $(round(prob, digits=6)))")
                end
            end
            println()
        end
    end
end

"""
    transfer_matrix(model::AnyonModel, τ::Float64; sign::Bool=true)

Construct the transfer matrix for two measurement layers.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `sign::Bool=true`: Measurement outcome sign

# Returns
- `Matrix{Float64}`: Transfer matrix between two measurement layers
"""
function transfer_matrix(model::AnyonModel{AT}, τ::Float64; sign::Bool=true) where {AT<:AbstractAnyonType}
    measurement_sites1, measure_type = _obtain_measurement_config(model.N, 1)  
    measurement_sites2, measure_type = _obtain_measurement_config(model.N, 2)  

    basis = anyon_basis(model)
    l = length(basis)
    TM = zeros(Float64, l, l)   

    for i in 1:l
        # First layer
        for (idx, site) in enumerate(measurement_sites1)
            outputstate1, outputstate2, output1, output2 = measure_basismap(model, τ, basis[i], site, sign)

            if output2 == 0
                TM[i,i] += output1
            else
                j2 = searchsortedfirst(basis, outputstate2)
                TM[i,i] += output1
                TM[i,j2] += output2
            end
        end
    end

    for i in 1:l
        # Second layer
        for (idx, site) in enumerate(measurement_sites2)
            outputstate1, outputstate2, output1, output2 = measure_basismap(model, τ, basis[i], site, sign)

            if output2 == 0
                TM[i,i] += output1
            else
                j2 = searchsortedfirst(basis, outputstate2)
                TM[i,i] += output1
                TM[i,j2] += output2
            end
        end
    end

    return TM
end

function _obtain_measurement_config(model::AnyonModel{FibonacciAnyon}, layer_idx::Int)
        measurement_sites = iseven(layer_idx) ? collect(1:2:model.N) : collect(2:2:model.N) # measurement sites for Fibonacci different layer, even layers measure at odd sites, odd layers measure at even sites,
        measure_operator = :Antiferro
        measure_anyon_model = AnyonModel(FibonacciAnyon(), model.N, pbc = model.pbc, measure_operator = measure_operator)
    return measurement_sites, measure_anyon_model
end

function _obtain_measurement_config(model::AnyonModel{IsingAnyon}, layer_idx::Int)
    # measure at all sites!!!
    measurement_sites = collect(1:model.N)
    measure_operator = iseven(layer_idx) ? :IsingZZ : :IsingX # measurement sites for Ising each layer, odd layers: measure X, even layers: measure ZZ
    measure_anyon_model = AnyonModel(IsingAnyon(), model.N, pbc = model.pbc, measure_operator = measure_operator)
    return measurement_sites, measure_anyon_model
end

"""
    _apply_measurement_layer!(model::AnyonModel, τ::Float64, state::Vector{T}, 
                              layer_sample::BitVector; layer_idx::Int=1) where {T}

Apply deterministic measurements to a layer with given measurement outcomes.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `state::Vector{T}`: Quantum state vector (modified in place)
- `layer_sample::BitVector`: Measurement outcomes for the layer
- `layer_idx::Int=1`: Layer index (1-based) to determine measurement pattern

# Returns
- `Tuple{Vector{T}, Float64}`: (post-measurement state, total free energy)
"""
function _apply_measurement_layer!(anyon_model::AnyonModel{AT}, τ::Float64, state::Vector{T}, 
    layer_sample::BitVector; layer_idx::Int64=1) where {T, AT<:AbstractAnyonType}
    # Helper function to apply deterministic measurements to a layer, connect measure on each site together.

    total_free_energy = zero(real(T))

    measurement_sites, measure_anyon_model = _obtain_measurement_config(anyon_model, layer_idx)  

    for (idx, sign) in enumerate(layer_sample)
        # Apply measurement at site measurement_sites[idx] with outcome sign
        state = measuremap(measure_anyon_model, τ, state, measurement_sites[idx], sign)
        prob = real(dot(state, state))
        total_free_energy += -log(prob)
        state ./= sqrt(prob)
    end

    return state, total_free_energy
end

"""
    _sample_layer!(model::AnyonModel, τ_eff::Float64, state::Vector{T};
                   layer_idx::Int=1, rng::MersenneTwister=MersenneTwister(), 
                   verbose::Bool=false) where {T}

Perform random measurement on a layer using Born rule sampling.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ_eff::Float64`: Effective measurement strength parameter
- `state::Vector{T}`: Quantum state vector
- `layer_idx::Int=1`: Layer index (1-based) to determine measurement pattern
- `rng::MersenneTwister`: Random number generator
- `verbose::Bool=false`: Whether to print debug information

# Returns
- `Tuple{Vector{T}, BitVector, Float64}`: (post-measurement state, sample outcomes, free energy)
"""
function _sample_layer!(anyon_model::AnyonModel{AT}, τ_eff::Float64, state::Vector{T};
    layer_idx::Int64=1,
    rng::MersenneTwister = MersenneTwister(),
    verbose::Bool=false) where {T, AT<:AbstractAnyonType}

    measurement_sites, measure_anyon_model = _obtain_measurement_config(anyon_model, layer_idx)  
    n = length(measurement_sites)
    sample = BitVector(zeros(Bool, n))
    F_layer = 0.0


    for (i, site) in enumerate(measurement_sites)
        # first 0 branch
        ψ0 = measuremap(measure_anyon_model, τ_eff, state, site, false)
        p0 = real(dot(ψ0, ψ0))
        p1 = 1 - p0

        randomNumber = rand(rng)
        verbose && @show randomNumber
        if randomNumber < p0
            sample[i] = 0
            state = ψ0 ./ sqrt(p0)
            F_layer += -log(p0)
            verbose && @show -log(p0)
        else
            # else 1 branch
            ψ1 = measuremap(measure_anyon_model, τ_eff, state, site, true)
            sample[i] = 1
            state = ψ1 ./ sqrt(p1)
            F_layer += -log(p1)
            verbose && @show -log(p1)
        end
    end
    return state, sample, F_layer
end

"""
Configuration struct for measurement evolution parameters.

# Fields
- `τ::Float64`: Measurement strength parameter
- `t₂::Int`: Number of measurement layers (time steps)
- `rng::MersenneTwister`: Random number generator (default: `MersenneTwister()`)
- `mode::Symbol`: Sampling mode, one of `:sample`, `:Born` (default: `:sample`)
- `t₁::Int`: Starting layer index for evolution (default: 1)
- `verbose::Bool`: Verbosity flag for detailed output (default: false)
- `enable_τ_eff::Bool`: Whether to enable half-strength measurement for the last layer (default: true)
"""
Base.@kwdef struct MeasureConfig
    τ::Float64
    t₂::Int
    rng::MersenneTwister  = MersenneTwister()
    mode::Symbol = :sample
    t₁::Int = 1
    verbose::Bool = false
    enable_τ_eff::Bool = true
end

"""
    measure_evolution!(model::AnyonModel, state::Vector{ET}, 
                       samples::Union{Nothing,BitMatrix}, 
                       measure_config::MeasureConfig) where {ET}
                  
Perform measurement evolution from t₁ to t₂ on the initial state.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `state::Vector{ET}`: Initial quantum state vector
- `samples::Union{Nothing,BitMatrix}`: Predefined measurement sequences for `:sample` mode
- `measure_config::MeasureConfig`: Configuration struct containing τ, t₁, t₂, mode, etc.

# Returns
- `Tuple{Vector{Vector{ET}}, BitMatrix, Vector{Float64}}`:
  - `states`: Intermediate states at each time step
  - `samples`: Measurement outcome sequences
  - `sample_free_energy`: Free energy for each layer

# Notes
- In `:Born` mode, samples are generated probabilistically via Born rule
- In `:sample` mode, samples must be provided as input
- (2N+1) layers of measurements correspond to N time steps of evolution
"""
function measure_evolution!(anyon_model::AnyonModel{AT},   # DRY: don't repeat yourself.
                  state::Vector{ET},
                  samples::Union{Nothing,BitMatrix},
                  measure_config::MeasureConfig) where {ET, AT<:AbstractAnyonType}
    # ---------- Sample decided according to mode ----------
    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    current_state = copy(state)
    if measure_config.mode == :Born
        _born_measure(anyon_model, current_state, measure_config)
    elseif measure_config.mode == :sample
        _sample_measure(anyon_model, current_state, samples, measure_config)
    end
end

function _born_measure(model::AnyonModel{AT}, current_state::Vector{ET}, measure_config::MeasureConfig) where {AT, ET}

    n_measure = measurement_num(model.anyon_type)*(model.N÷2)
    Δt = measure_config.t₂ - measure_config.t₁ + 1
    Δt >= 0 || error("t₂ must be >= t₁")
    D = Δt * 2 # number of layers to evolve

    # 1. Initialize sample matrix
    samples = BitMatrix(undef, (D, n_measure))   # to be filled during sampling
    sample_free_energy = zeros(D)
    states = Vector{Vector{ET}}(undef, Δt)

    for period in 1:Δt
    
        # Random sampling for this period
        τ = measure_config.τ
        enable_τ_eff = measure_config.enable_τ_eff
        rng = measure_config.rng
        verbose = measure_config.verbose
        τ_eff = (period == Δt && enable_τ_eff) ? τ/2 : τ
        current_state, samples[2*period-1, :], sample_free_energy[2*period-1] = _sample_layer!(model, τ, current_state, layer_idx=2*period-1, rng=rng, verbose=verbose)
        current_state, samples[2*period, :], sample_free_energy[2*period] = _sample_layer!(model, τ_eff, current_state, layer_idx=2*period, rng=rng, verbose=verbose)

        states[period] = current_state
    end

    return states, samples, sample_free_energy
end

function _sample_measure(model::AnyonModel{AT}, current_state::Vector{ET}, samples::BitMatrix, measure_config::MeasureConfig) where {AT, ET}

        n_measure = measurement_num(model.anyon_type)*model.N÷2
        Δt = measure_config.t₂ - measure_config.t₁ + 1
        Δt >= 0 || error("t₂ must be >= t₁")
        D = Δt * 2 # number of layers to evolve

        samples = BitMatrix(undef, (D, n_measure))   # to be filled during sampling
        sample_free_energy = zeros(D)
        states = Vector{Vector{ET}}(undef, Δt)
        # 2. Validate sample matrix
        isnothing(samples) && error("When mode=:sample samples must be ::BitMatrix")
        size(samples) == (D, n_measure) ||
            error("sample size should be ($D, $n_measure)")

        # 3. Deterministic trajectory for modes :sample, directly deterministic trajectory
    
        #  If anyon_type is :Fibo, the measurement sites are half of N, circuits belike:
        #   1   1   1   1   1   1   1   1   1
        #     1   1   1   1   1   1   1   1    
        #   1   1   1   1   1   1   1   1   1
        #   -τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ- (head tail concatenation)

        # If anyon_type is :IsingX or :IsingZZ, the measurement sites are N, circuits belike:
        #   Z₁Z₂ Z₂Z₃ Z₃Z₄ Z₄Z₅ Z₅Z₆ Z₆Z₇ Z₇Z₈ Z₈Z₁ (head tail concatenation)
        #  X    X    X    X    X    X    X    X
        #   Z₁Z₂ Z₂Z₃ Z₃Z₄ Z₄Z₅ Z₅Z₆ Z₆Z₇ Z₇Z₈ Z₈Z₁
        #  X    X    X    X    X    X    X    X
        #  ↑    ↑    ↑    ↑    ↑    ↑    ↑    ↑
        #  Or in majorana representation:
        # ---- --------  --------  --------  --------  --------  --------  --------  -----
        #  Z | |  ZZ  |  |  ZZ  |     ZZ  |  |  ZZ  |  |  ZZ  |  |  ZZ  |  |  ZZ  |  | Z
        # ---- --------  --------  --------  --------  --------  --------  --------  -----
        #  -------   -------   -------   -------   -------   -------   -------   ------
        #  |  X  |   |  X  |   |  X  |   |  X  |   |  X  |   |  X  |   |  X  |   |  X  |
        #  -------   -------   -------   -------   -------   -------   -------   ------
        #  γ₁   γ₂   γ₃   γ₄   γ₅   γ₆   γ₇   γ₈   γ₉  γ₁₀  γ₁₁  γ₁₂  γ₁₃  γ₁₄  γ₁₅  γ₁₆
    
        for period in 1:Δt
            # √M₁ᵉ √M₁ᵒ √M₁ᵉ √M₁ᵉ √M₁ᵒ √M₁ᵉ ⋯ √M₁ᵉ √M₁ᵒ √M₁ᵉ→ √M₁ᵉ M₁ᵒ M₁ᵉ M₁ᵒ ⋯ M₁ᵉ M₁ᵒ √M₁ᵉ. 
            # √X √ZZ √X √X √ZZ √X ⋯ √X √ZZ √X→ √X ZZ X ZZ ⋯ X ZZ √X. To ensure each layer is hermitian, first layer doesn't matter.
            # Or √ZZ X √ZZ √ZZ X √ZZ ⋯ X √ZZ X √ZZ→ X ZZ X ZZ ⋯ X √ZZ, also works
            τ = measure_config.τ
            enable_τ_eff = measure_config.enable_τ_eff
            τ_eff = (period == Δt && enable_τ_eff) ? τ/2 : τ
            current_state, sample_free_energy[2*period-1] = _apply_measurement_layer!(
                            model, τ, current_state,
                            samples[2*period-1, :], layer_idx=2*period-1)
            current_state, sample_free_energy[2*period] = _apply_measurement_layer!(
                            model, τ_eff, current_state,
                            samples[2*period, :], layer_idx=2*period)

            states[period] = current_state
        end
    return states, samples, sample_free_energy
end

struct Measurement_outcome{T}
    state::Vector{T}
    samples::Matrix{Int}
    free_energy::Float64
end

function generate_state_by_measurement(anyon_model::AnyonModel{AT}, state::Vector{T}, samples, measure_config::MeasureConfig) where{T, AT<:AbstractAnyonType}
    N = measurement_num(anyon_model.anyon_type) * size(samples, 2)
    D = size(samples, 1) # number of layers
    t₂ = D ÷ 2 # number of time steps/ periods

    final_state, sample, free_energy = measure_evolution!(anyon_model, state, samples,  measure_config)
    return Measurement_outcome(final_state, sample, free_energy)
end
measurement_num(::FibonacciAnyon) = 1
measurement_num(::IsingAnyon) = 2

function generate_state_by_measurement(anyon_type::AbstractAnyonType, τ::Float64, state::Vector{T}, sample::BitVector; pbc::Bool=true, layer_idx::Int=1, enable_τ_eff::Bool=true) where{T} 
    N = (anyon_type == :Fibo) ? length(sample) * 2 : length(sample)
    D = size(sample, 1) # number of layers
    t₂ = D ÷ 2 # number of time steps/ periods

    final_state, sample, free_energy = measure_evolution!(N, τ, state, t₂; 
    pbc=pbc, anyon_type=anyon_type, mode=:sample, 
    sample=sample, enable_τ_eff=enable_τ_eff)
    return Measurement_outcome(final_state, sample, free_energy)
end

"""
    boundary_measure(model::AnyonModel, state::Vector{ET}, layer_idx::Int=1, 
                     num_samples::Int=1000, measure_config::MeasureConfig) where {ET}

Generate measurement samples at boundary sites with probabilistic outcomes (Born rule).

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `state::Vector{ET}`: Initial quantum state vector
- `layer_idx::Int=1`: Measurement layer index
- `num_samples::Int=1000`: Number of measurement samples to generate
- `measure_config::MeasureConfig`: Configuration containing rng and other parameters

# Returns
- `Tuple{Vector{Vector{Float64}}, BitMatrix, Vector{Float64}}`:
  (post-measurement states, measurement sequences, free energies)
"""
function boundary_measure(model::AnyonModel{AT}, state::Vector{ET}, layer_idx::Int=1, num_samples::Int=1000, measure_config::MeasureConfig) where {ET, AT<:AbstractAnyonType}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    n_measure = measurement_num(model.anyon_type)*(model.N ÷ 2)
    sample_measured_states = Vector{Vector{Float64}}(undef, num_samples)
    samples = BitMatrix(undef, num_samples, n_measure)
    sample_free_energy = Vector{Float64}(undef, num_samples)

    for sample_idx in 1:num_samples
        final_state, sample, free_energy = _sample_layer!(model, state, rng, layer_idx = layer_idx)

        sample_measured_states[sample_idx] = final_state
        samples[sample_idx, :] = sample
        sample_free_energy[sample_idx] = free_energy
    end
    
    return sample_measured_states, samples, sample_free_energy
end


"""
    boundary_post_selection(model::AnyonModel, state::Vector{ET}, 
                            layer_idx::Int=1, sign::Bool=true) where {ET}

Apply deterministic measurements at boundary sites with specified outcomes.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `state::Vector{ET}`: Initial quantum state vector
- `layer_idx::Int=1`: Measurement layer index
- `sign::Bool=true`: Measurement outcome (false for all +, true for all -)

# Returns
- `Tuple{Vector{ET}, BitVector, Float64}`:
  (post-measurement state, measurement outcomes, total free energy)
"""
function boundary_post_selection(model::AnyonModel{AT}, state::Vector{ET}, layer_idx::Int=1, sign::Bool=true) where {ET, AT<:AbstractAnyonType}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    n_measure = measurement_num(model.anyon_type)*(model.N ÷ 2)
    sample = (sign == 0) ? BitVector(zeros(Bool, n_measure)) : BitVector(ones(Bool, n_measure))

    state_measured, total_free_energy =  _apply_measurement_layer!(
        model, state, sample, layer_idx
    )  # return the final state after applying all measurements in the layer

    return state_measured, sample, total_free_energy
end


"""
    bulk_measure(model::AnyonModel, state::Vector{ET}, t::Int64, 
                 measure_config::MeasureConfig) where {ET}

Generate measurement samples at bulk sites with probabilistic outcomes.

Performs both time and spatial axis evolution for bulk measurements.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `state::Vector{ET}`: Initial quantum state vector
- `t::Int64`: Number of measurement layers (depth)
- `measure_config::MeasureConfig`: Measurement configuration parameters

# Returns
- `Tuple{Vector, BitMatrix, Float64}`:
  (final state, measurement samples, sample free energy)

Samples measurement outcomes probabilistically based on Born rule.
"""
function bulk_measure(model::AnyonModel{AT}, state::Vector{ET}, t::Int64, measure_config::MeasureConfig) where {ET, AT<:AbstractAnyonType}

    final_state, samples, sample_free_energy = measure_evolution!(
        model, state, t, measure_config;
    )

    return final_state, samples, sample_free_energy
end

"""
    bulk_post_selection(model::AnyonModel, state::Vector{ET}, sign::Bool,
                        measure_config::MeasureConfig) where {ET}

Generate measurement samples at bulk sites with deterministic (post-selected) outcomes.

Performs both time and spatial axis evolution with specified measurement outcomes.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `state::Vector{ET}`: Initial quantum state vector
- `sign::Bool`: Measurement outcome (false for all +, true for all -)
- `measure_config::MeasureConfig`: Measurement configuration parameters

# Returns
- `Tuple{Vector, BitMatrix, Float64}`:
  (post-measurement states, measurement samples, sample free energy)
"""
function bulk_post_selection(
        model::AnyonModel{AT},
        state::Vector{ET},
        sign::Bool,
        measure_config::MeasureConfig
    ) where {ET}

    # 1. Build the sample, all 0 or all 1, each layer N/2 or N measurements, depending on anyon_type, total D layers.
    n_measure = measurement_num(model.anyon_type)*(model.N ÷ 2)

    sample = (sign == 1) ? BitMatrix(undef, (2D, n_measure)) : BitMatrix(undef, (2D, n_measure))

    # 2. generate_state to run all the layers
    sample_measured_states, samples, sample_free_energy = measure_evolution!(
        model, state, sample, measure_config
    )        # need all the intermediate states, measure_evolution! return each layer free energy

    return sample_measured_states, sample, sample_free_energy
end



"""
    bayes_distort(γ::Float64, trajectories::Vector{Int64}, probabilities::Vector{Float64})

Distort measurement trajectories based on a Bayesian distortion factor γ.

This function implements the distortion process where each faithful sample s 
is converted to a distorted sample s̃ according to the conditional probability:
    P(s̃|s) = ∏ⱼ (1 + γ s̃ⱼ sⱼ)/2

Only works for single-layer measurements to generate new samples based on 
the projective limit measurement.

# Arguments
- `γ::Float64`: Distortion factor (readout fidelity parameter, 0 ≤ γ ≤ 1)
- `trajectories::Vector{Int64}`: Vector of measurement trajectories
- `probabilities::Vector{Float64}`: Corresponding probabilities for each trajectory

# Returns
- `Tuple{Vector{Vector{Int64}}, Vector{Float64}}`:
  (distorted_trajectories, distorted_probabilities) in corresponding order
"""
function bayes_distort(γ::Float64, trajectories::Vector{Int64}, probabilities::Vector{Float64})
    
    # Dictionary to store the distorted trajectory probabilities
    distorted_prob_dict = Dict{Vector{Int64}, Float64}()
    n_sites = length(trajectories)
    distorted_prob = Vector{Vector{Float64}}(undef, n_sites)
    transfer_matrix = [1 + γ 1 - γ; 1 - γ 1 + γ] / 2
    
    # For each original trajectory
    for (traj_idx, original_traj) in enumerate(trajectories)
        original_prob = probabilities[traj_idx]
        prob_distribution = (trajectories[traj_idx] == 1) ? [original_prob, 1 - original_prob] : [1 - original_prob, original_prob]
        distorted_prob[traj_idx] = transfer_matrix *prob_distribution
    end
    
    # Generate all possible distorted trajectories (2^n possibilities)
    for distorted_bits in 0:(2^n_sites - 1)
        # Convert bit representation to ±1 trajectory
        prob = 1.0
        distorted_traj = Vector{Int64}(undef, n_sites)
        for j in 1:n_sites
            # Extract j-th bit and convert to ±1
            bit = (distorted_bits >> (j-1)) & 1
            distorted_traj[j] = bit
            prob*= distorted_prob[j][(bit==1) ? 1 : 2]  # bit + 1 because Julia is 1-indexed
        end
        
        
        
        if haskey(distorted_prob_dict, distorted_traj)
            distorted_prob_dict[distorted_traj] = prob
        else
            distorted_prob_dict[distorted_traj] = prob
        end
    end
    # Convert dictionary to vectors
    distorted_trajectories = collect(keys(distorted_prob_dict))
    distorted_probabilities = collect(values(distorted_prob_dict))
    
    return distorted_trajectories, distorted_probabilities
end