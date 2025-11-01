struct AnyonModel{AT<:AbstractAnyonType}
    anyon_type::AT
    N::Int   # system size
    pbc::Bool
end

"""
    measure_basismap(anyon_model::AnyonModel, τ::Float64, state::T, i::Int, sign::Int64) where {T}

Map single basis state under measurement operation at site i.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `τ::Float64`: Measurement strength parameter
- `state::T`: Input basis state
- `i::Int`: Measurement site index (1 ≤ i ≤ N)
- `sign::Int64`: Measurement outcome (0 for +, 1 for -)
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- Basis-dependent output: Either `(basis, coefficient)` or `(basis1, basis2, coeff1, coeff2)`

Maps individual basis states according to measurement protocols and fusion rules. Here we choose Heisenberg-like preferring way, selecting specific fusion outcome.

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis

julia> N = 4; T = BitStr{N, Int};

julia> state = T(0b0100);  # Single τ at site 2

julia> τ = 1.0;  # Measurement strength

julia> # Measure at site 2 with outcome + (sign=0)
       result = measure_basismap(T, τ, state, 2, 0, true);

julia> length(result) ∈ [2, 4]  # Returns 2 or 4 elements depending on configuration
true

julia> # The first element is always a basis state
       typeof(result[1]) == T
true
```
"""
function measure_basismap(::Type{T}, τ::Float64, state::T, i::Int, sign::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}
    # default for PBC system, map basis (not state!!!), and index count from the left.
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (0, 1) "sign must be either 0 the plus, 1 the minus"
    
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))
    ϕ = (1+√5)/2

    return _apply_result(anyon_type, τ, state, i, sign, measure_operator)
end

function _apply_result(anyontype::FibonacciAnyon, τ::Float64, state::T, i::Int, sign::Int64, measure_operator::Symbol) where {T}
    @assert measure_operator ∈ [:reset, :Fibo] "measure_operator must be :Fibo or :reset"
    if τ >= 1e2
        cstτ = 0.5
        coef = sign == 0 ? 0.5 : -0.5
    else
        cstτ = (exp(τ) + 1) / (2 * √(exp(2τ) + 1))
        coef = sign == 0 ? (exp(τ) - 1) / (2 * √(exp(2τ) + 1)) : (1 - exp(τ)) / (2 * √(exp(2τ) + 1))
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

    if pbc
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

function _apply_result(anyontype::IsingType, τ::Float64, state::T, i::Int, sign::Int64, measure_operator::Symbol) where {T}
    @assert measure_operator in [:X, :ZZ] "measure_operator must be either :X or :ZZ"
    if measure_operator == :X
        if τ >= 1e2
            cstτ = 0.5
            coef = sign == 0 ? 0.5 : -0.5
        else
            cstτ = cosh(τ/2) / √(2cosh(τ))
            coef = sign == 0 ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
        end

        return state, X(state,i), cstτ, coef

    elseif measure_operator == :ZZ
        if τ >= 1e2
            cstτ = 0.5
            coef = sign == 0 ? 0.5 : -0.5
        else
            cstτ = cosh(τ/2) / √(2cosh(τ))
            coef = sign == 0 ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
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

function _apply_result(anyontype::IsingType, τ::Float64, state::T, i::Int, sign::Int64) where {T}
        if τ >= 1e2
            cstτ = 0.5
            coef = sign == 0 ? 0.5 : -0.5
        else
            cstτ = cosh(τ/2) / √(2cosh(τ))
            coef = sign == 0 ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
        end

        return state, state, (state[N - i + 1] == 0) ? cstτ + coef : cstτ - coef, 0
end


function measure_matrix(::Type{T}, τ::Float64, idx::Int, sign::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}

    if anyon_type == :Fibo
        @assert pbc || (2 <= idx <= N-1) "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif anyon_type == :IsingZZ
        @assert pbc || (1 <= idx <= N-1) "Index idx must be in [1, N-1] for open BC (IsingZZ)"
    elseif anyon_type ∈ (:IsingX, :IsingZ, :reset, :resetFibo)
        @assert pbc || (1 <= idx <= N) "Index idx must be in [1, N] for open BC (IsingX)"
    else
        error("Unknown measure class: $anyon_type")
    end

    basis = anyon_basis(T, pbc; anyon_type = anyon_type)
    l = length(basis)
    Bmatrix = zeros(l, l)

    for i in 1:l
        s1, s2, w1, w2 = measure_basismap(T, τ, basis[i], idx, sign, pbc; anyon_type = anyon_type)

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

function measuremap(::Type{T}, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}
    # input a superposition state, and output the measured state (tedancy fusion to 0 or 1 in Fibonacci measure class, or X ZZ in Ising measure class)
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
    l=length(basis)
    @assert l == length(state) "state length is expected to be $(l), but got $(length(state))"
    mapped_state = zeros(ET, length(state))
    for i in 1:l
        outputstate1, outputstate2, output1, output2 = measure_basismap(T, τ, basis[i], idx, sign, pbc, anyon_type=anyon_type)
        
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
measuremap(N::Int, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET} = measuremap(BitStr{N, Int}, τ, state, idx, sign, pbc, anyon_type=anyon_type)

function laddermeasuremap(::Type{T}, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true, anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}
    # input a superposition state, and output the braided state
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    basis=anyon_basis(T, pbc, anyon_type=anyon_type)
    l=length(basis)
    @assert l^2 == length(state) "state length is expected to be $(l^2), but got $(length(state))"
    mapped_state = zeros(ET, length(state))
    for i in 1:l
        for j in 1:l
            output1 = measure_basismap(T, τ, basis[i], idx, sign, pbc, anyon_type=anyon_type)
            output2 = measure_basismap(T, τ, basis[j], idx, sign, pbc, anyon_type=anyon_type)
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
laddermeasuremap(N::Int, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET} = laddermeasuremap(BitStr{N, Int}, τ, state, idx, sign, pbc, anyon_type=anyon_type)

"""
    measurement_enumeration(::Type{T}, τ::Float64, initial_state::Vector{ET}, measurement_sites::Vector{Int}, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}

Enumerating all trajectories of measurements on a given initial state.

Args:T, τ, initial_state, measurement_sites, pbc

Returns:
    final_states, trajectories, probabilities
"""
function measurement_enumeration(anyon_model::AnyonModel, τ::Float64, initial_state::Vector{ET}, measurement_sites::Vector{Int}) where {ET}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"
    
    # Initialize, only one initial state
    current_level_states = [copy(initial_state)]
    current_level_trajectories = [Int64[]]
    current_level_probabilities = [1.0]
    
    for (measurement_idx, site) in enumerate(measurement_sites)
        next_level_states = Vector{Vector{ET}}()
        next_level_trajectories = Vector{Vector{Int64}}()
        next_level_probabilities = Vector{Float64}()
        
        # Branching for each current state
        for (state_idx, state) in enumerate(current_level_states)
            current_trajectory = current_level_trajectories[state_idx]
            current_prob = current_level_probabilities[state_idx]
            
            state_after_p = measuremap(T, τ, state, site, 0, pbc, anyon_type=anyon_type)
            prob_p = real(dot(state_after_p, state_after_p))

            normalized_state_p = state_after_p / sqrt(prob_p)
            new_trajectory_p = [current_trajectory; 0]
            new_prob_p = current_prob * prob_p
            
            push!(next_level_states, normalized_state_p)
            push!(next_level_trajectories, new_trajectory_p)
            push!(next_level_probabilities, new_prob_p)
        
            
            state_after_m = measuremap(T, τ, state, site, 1, pbc, anyon_type=anyon_type)
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

measurement_enumeration(N::Int, τ::Float64, initial_state::Vector{ET}, measurement_sites::Vector{Int}, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET} = 
measurement_enumeration(BitStr{N, Int}, τ, initial_state, measurement_sites, pbc, anyon_type=anyon_type)


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

function _obtain_measurement_config(N::Int, layer_idx::Int, anyon_type::Symbol=:Fibo)
    if anyon_type == :Fibo
        measurement_sites = iseven(layer_idx) ? collect(1:2:N) : collect(2:2:N) # measurement sites for Fibonacci different layer, even layers measure at odd sites, odd layers measure at even sites, 
        measure_type = :Fibo
    elseif anyon_type ∈ (:IsingX, :IsingZ, :IsingZZ, :reset, :resetFibo)
        # measure at all sites!!!
        measurement_sites = collect(1:N)
        measure_type = iseven(layer_idx) ? :IsingZZ : :IsingX # measurement sites for Ising each layer, odd layers: measure X, even layers: measure ZZ
    else
        error("Unknown measure class: $anyon_type")
    end

    return measurement_sites, measure_type
end

"""
    _apply_measurement_layer!(N::Int64, τ::Float64, state::Vector{T}, 
    layer_sample::Vector{Int64}, 
    layer_idx::Int64, 
    pbc::Bool=true; 
    anyon_type::Symbol=:Fibo) where {T}

Generate measurement samples at 1 layer with post_selection outcomes, i.e., with time step 1.

# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `layer_sample::Vector{Int64}`: Measurement outcomes for the layer
- `layer_idx::Int64`: Layer index (1-based) to determine measurement pattern
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
`Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
- `post-measurement states`
- `free energies`

Samples measurement outcomes with given measurement outcomes.
"""
function _apply_measurement_layer!(anyon_model::AnyonModel, τ::Float64, state::Vector{T}, 
    layer_sample::Vector{Int64}, layer_idx::Int64) where {T}
    # Helper function to apply deterministic measurements to a layer, connect measure on each site together.

    total_free_energy = zero(real(T))

    measurement_sites, measure_type = _obtain_measurement_config(N, layer_idx, anyon_type)  

    for (idx, sign) in enumerate(layer_sample)
        # Apply measurement at site measurement_sites[idx] with outcome sign
        state = measuremap(N, τ, state, measurement_sites[idx], sign, pbc, anyon_type = measure_type)
        prob = real(dot(state, state))
        total_free_energy += -log(prob)

        state ./= sqrt(prob)
    end

    return state, total_free_energy
end

"""
    _sample_layer!(N::Int64, τ_eff::Float64, state::Vector{T}, 
    rng::MersenneTwister = MersenneTwister(), 
    layer_idx::Int64, 
    pbc::Bool=true; 
    anyon_type::Symbol=:Fibo) where {T}

Do the random measurement on a layer, in contrast to _apply_measurement_layer! No sample input, but output the sample.

# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `rng::MersenneTwister = MersenneTwister()`: Random number generator
- `layer_idx::Int64`: Layer index (1-based) to determine measurement pattern
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
`Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
- `state_Born_measured`, 
- `samples`, 
- `free_energy` of one layer.

Samples measurement outcomes with Born rule driven trajectories.
"""
function _sample_layer!(N::Int64, τ_eff::Float64, state::Vector{T}, 
    rng::MersenneTwister = MersenneTwister(), 
    layer_idx::Int64=1, 
    pbc::Bool=true; 
    anyon_type::Symbol=:Fibo, 
    verbose::Bool=false) where {T}

    measurement_sites, measure_type = _obtain_measurement_config(N, layer_idx, anyon_type)  
    n = length(measurement_sites)
    sample = zeros(Int, n)
    F_layer = 0.0


    for (i, site) in enumerate(measurement_sites)
        # first 0 branch
        ψ0 = measuremap(N, τ_eff, state, site, 0, pbc, anyon_type = measure_type)
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
            ψ1 = measuremap(N, τ_eff, state, site, 1, pbc, anyon_type = measure_type)
            sample[i] = 1
            state = ψ1 ./ sqrt(p1)
            F_layer += -log(p1)
            verbose && @show -log(p1)
        end
    end
    return state, sample, F_layer
end

"""
- `τ::Float64`: Measurement strength parameter
- `t₂::Int`: Number of measurement layers (time steps)
- `rng::MersenneTwister = MersenneTwister()`: Random number generator
- `mode::Symbol = :sample`: Sampling mode, one of `:sample`, `:Born`
- `t₁::Int = 1`: Starting layer index for evolution (default is 1)
- `verbose::Bool = false`: Verbosity flag for detailed output
- `enable_τ_eff::Bool = true`: Whether to enable half-strength measurement for the last layer
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
    measure_evolution!(N::Int,
                  τ::Float64,
                  state::Vector{ET},
                  t₂::Int = 1;
                  rng::MersenneTwister = MersenneTwister(),
                  pbc::Bool = true,
                  anyon_type::Symbol = :Fibo,
                  mode::Symbol = :sample,
                  t₁::Int = 1,
                  sample::Union{Nothing,Matrix{Int}}=nothing,
                  verbose::Bool = false, enable_τ_eff::Bool = true) where {ET}
                  
Perform measurement evolution from t₁ (default to be 1) to t₂ on the initial state with specified parameters, returning final states, measurement sequences, and free energies, with options for different sampling modes and temporal settings. KEEP IN MIND: In such definition, (2N+1) layers of measurements correspond to N time steps of evolution, with the last layer being half-strength. The first and end layer should not be √ZZ or √M₁ᵒ

# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `t₂::Int=1`: Number of measurement layers (time steps)
- `rng`: Random number generator (default: `MersenneTwister()`)
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type
- `mode::Symbol=:prob`: evolution mode, one of `:prob`, `:sample`, `:Born`, Born is random sampling, driven by Born rule, the other is deterministic post-selection evolution, with given measurement outcomes.
- `sample::Union{Nothing,Matrix{Int}}=nothing`: Predefined measurement sequences for `:sample` mode
- `t₁::Int=1`: Starting layer index for evolution (default is 1)
- `verbose::Bool=false`: Verbosity flag for detailed output
- `enable_τ_eff::Bool=true`: Whether to enable half-strength measurement for the last layer

# Returns
- `Tuple{Union{Vector{ET}, Vector{Vector{ET}}}, Matrix{Int}, Vector{Float64}}`: 
  (final state or all intermediate states, measurement sequences, free energies)
- states        : final state
- samples       : measurement sequences (Matrix{Int})
- sample_free_energy : cumulative free energy for each layer (or each sample)
"""
function measure_evolution!(anyon_model::AnyonModel,   # DRY: don't repeat yourself.
                  state::Vector{ET},
                  sample::Union{Nothing,Matrix{Int}},
                  measure_config::MeasureConfig) where {ET}

    n_measure = anyon_type == :Fibo ? N÷2 : N

    # ---------- Sample decided according to mode ----------
    Δt = t₂ - t₁ + 1 
    Δt >= 0 || error("t₂ must be >= t₁")
    D = Δt * 2 # number of layers to evolve
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    sample_free_energy = zeros(D) # free energy of each layer
    states = Vector{Vector{ET}}(undef, Δt)  # states of each layer
    current_state = copy(state)
    if measure_config.mode == :Born
        _born_measure(model, current_state, sample, measure_config)
    elseif measure_config.mode == :sample
        _sample_measure(model, current_state, sample, measure_config)
    end
end

function _born_measure(model, current_state, sample, measure_config)
         # 1. Initialize sample matrix
        sample = zeros(Int, D, n_measure)   # to be filled during sampling

        for period in 1:Δt
        
            # Random sampling for this period
            τ_eff = (period == Δt && enable_τ_eff) ? τ/2 : τ
            current_state, sample[2*period-1, :], sample_free_energy[2*period-1] = _sample_layer!(N, τ, current_state, rng, 2*period-1, pbc, anyon_type = anyon_type, verbose=verbose)
            current_state, sample[2*period, :], sample_free_energy[2*period] = _sample_layer!(N, τ_eff, current_state, rng, 2*period, pbc, anyon_type = anyon_type, verbose=verbose)

            states[period] = current_state
        end
        
    elseif mode == :sample
        isnothing(sample) && error("When mode=:sample sample must be ::Matrix{Int}")
        size(sample) == (D, n_measure) ||
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
            τ_eff = (period == Δt && enable_τ_eff) ? τ/2 : τ
            current_state, sample_free_energy[2*period-1] = _apply_measurement_layer!(
                            N, τ, current_state,
                            sample[2*period-1, :], 2*period-1, pbc;
                            anyon_type=anyon_type)
            current_state, sample_free_energy[2*period] = _apply_measurement_layer!(
                            N, τ_eff, current_state,
                            sample[2*period, :], 2*period, pbc;
                            anyon_type=anyon_type)

            states[period] = current_state
        end
    end
    
    return states, sample, sample_free_energy
end

struct StateByMeasurement{T}
    state::Vector{T}
    samples::Matrix{Int}
    free_energy::Float64
end
function generate_state_by_measurement(anyon_model::AnyonModel, state::Vector{T}, sample::Matrix{Int}; τ::Float64, enable_τ_eff::Bool=true) where{T}
    N = measurement_time(anyon_model.anyon_type) * size(sample, 2)
    D = size(sample, 1) # number of layers
    t₂ = D ÷ 2 # number of time steps/ periods

    final_state, sample, free_energy = measure_evolution!(anyon_model, τ, state, t₂; 
        mode=:sample, 
        sample=sample, enable_τ_eff=enable_τ_eff)
    return StateByMeasurement(final_state, sample, free_energy)
end
measurement_time(::FibonacciAnyon) = 2
measurement_time(::IsingAnyon) = 1

function generate_state_by_measurement(anyon_type::AbstractAnyonType, τ::Float64, state::Vector{T}, sample::Vector{Int}; pbc::Bool=true, layer_idx::Int=1, enable_τ_eff::Bool=true) where{T} 
    N = (anyon_type == :Fibo) ? length(sample) * 2 : length(sample)
    D = size(sample, 1) # number of layers
    t₂ = D ÷ 2 # number of time steps/ periods

    final_state, sample, free_energy = measure_evolution!(N, τ, state, t₂; 
    pbc=pbc, anyon_type=anyon_type, mode=:sample, 
    sample=sample, enable_τ_eff=enable_τ_eff)
    return StateByMeasurement(final_state, sample, free_energy)
end

"""
    boundary_measure(::Type{T}, τ::Float64, state::Vector{ET}, layer_idx::Int, num_samples::Int=1000, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}

Generate measurement samples at boundary sites with probabilistic outcomes, i.e., without time axis evolution. DO num_samples times sampling.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `layer_idx::Int`: measurement layer index, (layer_idx is even, measurement sites is collect(1:2:N), layer_idx is odd, measurement sites is collect(2:2:N) for Fibonacci anyon type, layer_idx is always 1 for Ising anyon type)
- `num_samples::Int=1000`: Number of measurement samples to generate
- `rng::MersenneTwister`: Random number generator
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Tuple{Vector{Vector{Float64}}, Matrix{Int64}, Vector{Float64}}`: 
  (post-measurement states, measurement sequences, free energies)

Samples measurement outcomes probabilistically based on Born rule.
"""
function boundary_measure(N::Int, τ::Float64, state::Vector{ET}, layer_idx::Int=1, num_samples::Int=1000, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"
    
    n_measure = anyon_type == :Fibo ? N ÷ 2 : N
    sample_measured_states = Vector{Vector{Float64}}(undef, num_samples)
    samples = zeros(Int, num_samples, n_measure)
    sample_free_energy = Vector{Float64}(undef, num_samples)

    for sample_idx in 1:num_samples
        final_state, sample, free_energy = _sample_layer!(N, τ, state,  rng, layer_idx, pbc; anyon_type=anyon_type)

        sample_measured_states[sample_idx] = final_state
        samples[sample_idx, :] = sample
        sample_free_energy[sample_idx] = free_energy
    end
    
    return sample_measured_states, samples, sample_free_energy
end


"""
    boundary_post_selection(N::Int64, τ::Float64, state::Vector{ET}, layer_idx::Int=1, sign::Int64=1, pbc::Bool=true; anyon_type::Symbol=:Fibo)

Generate measurement samples at boundary sites with post_selection outcomes, i.e., given spatial evolution without time axis evolution.

# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `layer_idx::Int=1`: measurement layer index, (layer_idx is even, measurement sites is collect(1:2:N), layer_idx is odd, measurement sites is collect(2:2:N) for Fibonacci anyon type, layer_idx is always 1 for Ising anyon type)
- `sign::Int64=1`: Measurement outcome sign (0 or 1)
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
  (post-measurement states, measurement sequences, free energies)

Samples measurement outcomes with given measurement outcomes.
"""
function boundary_post_selection(N::Int64, τ::Float64, state::Vector{ET}, layer_idx::Int=1, sign::Int64=1, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    n_measure = (anyon_type == :Fibo) ? N ÷ 2 : N
    sample = (sign == 0) ? zeros(Int, n_measure) : ones(Int, n_measure)

    state_measured, total_free_energy =  _apply_measurement_layer!(
        N, τ, state, sample, layer_idx, pbc;
        anyon_type=anyon_type
    )  # return the final state after applying all measurements in the layer

    return state_measured, sample, total_free_energy
end


"""
    bulk_measure(N::Int64, τ::Float64, state::Vector{ET}, t₂::Int64, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true; anyon_type::Symbol=:Fibo)

Generate measurement samples at bulk sites with probabilistic outcomes, i.e., with time axis evolution, together with spatial evolution axis as bulk.

# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `t₂::Int64`: Number of measurement layers (depth), or time step (over 2)
- `rng::MersenneTwister`: Random number generator
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
  (post-measurement states, measurement sequences, free energies)

Samples measurement outcomes probabilistically based on Born rule.
"""
function bulk_measure(N::Int64, τ::Float64, state::Vector{ET}, D::Int64, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true; anyon_type::Symbol=:Fibo, verbose::Bool=false) where {ET}

    final_state, samples, sample_free_energy = measure_evolution!(
        N, τ, state, D;
        rng=rng,
        pbc=pbc,
        anyon_type=anyon_type,
        mode=:Born,
        verbose=verbose
    )  # need all the intermediate states, measure_evolution! return each layer free energy
    return final_state, samples, sample_free_energy
end

"""
    bulk_post_selection(N::Int64, τ::Float64, state::Vector{ET}, t₂::Int64, sign::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo)

Generate measurement samples at bulk sites with post_selection outcomes, i.e., with time axis evolution, together with spatial evolution axis as bulk.

# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `t₂::Int64`: Number of measurement layers (depth), or time step (over 2)
- `sign::Int64`: Measurement sign, 0 for positive, 1 for negative
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
  (post-measurement states, measurement sequences, free energies)

Samples measurement outcomes with given measurement outcomes.
"""
function bulk_post_selection(
        N::Int64,
        τ::Float64,
        state::Vector{ET},
        D::Int64,
        sign::Int64,
        pbc::Bool=true;
        anyon_type::Symbol=:Fibo
    ) where {ET}

    # 1. Build the sample, all 0 or all 1, each layer N/2 or N measurements, depending on anyon_type, total D layers.
    n_measure = (anyon_type == :Fibo) ? N ÷ 2 : N

    sample = (sign == 1) ? ones(Int, 2D, n_measure) : zeros(Int, 2D, n_measure)

    # 2. generate_state to run all the layers
    sample_measured_states, samples, sample_free_energy = measure_evolution!(
        N, τ, state, D;
        pbc=pbc,
        anyon_type=anyon_type,
        mode=:sample,
        sample=sample
    )        # need all the intermediate states, measure_evolution! return each layer free energy

    return sample_measured_states, sample, sample_free_energy
end



"""
Distort the measurement trajectories based on a Bayesian distortion factor γ. Noting that it only works for one layer measurement to generate new sample for the other factor γ based on Projective limit measurement.

This function implements the distortion process where each faithful sample s is converted to a distorted sample s̃ according to the conditional probability:
P(s̃|s) = ∏ⱼ (1 + γ s̃ⱼ sⱼ)/2

Args:
    γ: Distortion factor (readout fidelity parameter, 0 ≤ γ ≤ 1).
    trajectories: Vector of measurement trajectories.
    probabilities: Corresponding probabilities for each trajectory.
    
Returns:
    Tuple of (distorted_trajectories, distorted_probabilities) where:
    - distorted_trajectories: All possible distorted trajectories
    - distorted_probabilities: Their corresponding probabilities after distortion
    in corresponding order.
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