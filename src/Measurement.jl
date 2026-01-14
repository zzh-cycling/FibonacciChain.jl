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
    measure_operator = model.measure_operator
    @assert measure_operator ∈ [:Ferro, :Antiferro, :reset] "measure_operator must be :Ferro, :Antiferro, :reset"
    N = model.N
    @assert num_digits(T) == N "State length mismatch: expected $(N), got $(num_digits(T))"

    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))
    ϕ = (1+√5)/2

    if measure_operator == :reset && τ >= 1e2
        cstτ = 0.5
        coef = sign ? -0.5 : 0.5
        
        return state, state, (state[N - i + 1] == 0) ? cstτ + coef : cstτ - coef, 0.0
    else
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
                return state, state, cstτ+coef, 0.0
            elseif state & mask == str100
                return state, state, cstτ+coef, 0.0
            elseif state & mask == str101
                return state, state, cstτ-coef, 0.0
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
                    return state, state, cstτ+coef, 0.0
                elseif state & mask == str100
                    return state, state, cstτ+coef, 0.0
                elseif state & mask == str101
                    return state, state, cstτ-coef, 0.0
                end
            elseif i == N #count from the left
            mask=bmask(T, N, 2, 1)
            str100, str101, str010, str001, str000 = bmask(T,2), bmask(T, N, 2), bmask(T, 1), bmask(T, N), T(0)
                if state & mask == str000
                    return state, X(state,i), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2)
                elseif state & mask == str010
                    return state, X(state,i), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2)
                elseif state & mask == str001
                    return state, state, cstτ+coef, 0.0
                elseif state & mask == str100
                    return state, state, cstτ+coef, 0.0
                elseif state & mask == str101
                    return state, state, cstτ-coef, 0.0
                end
            end
        end
    end
end

function _apply_result(model::AnyonModel{IsingAnyon}, τ::Float64, state::T, i::Int, sign::Bool) where {T}
    measure_operator = model.measure_operator
    @assert measure_operator in [:X, :ZZ, :Z, :reset] "measure_operator must be either :X, :ZZ, :Z, :reset"

    N = model.N
    fl = bmask(T, N)
    X(state, i) = flip(state, fl >> (i-1))

    # Common coefficients for all operators except :reset/:Z
    if τ >= 1e2
        cstτ = 0.5
        coef = sign ? -0.5 : 0.5
    else
        cstτ = cosh(τ/2) / √(2cosh(τ))
        coef = sign ? -sinh(τ/2) / √(2cosh(τ)) : sinh(τ/2) / √(2cosh(τ))
    end

    # Helper: get ZZ eigenvalue for sites (j1, j2)
    zz_eigen(j1, j2) = ((state >> (N - j1)) & 1) == ((state >> (N - j2)) & 1) ? 1 : -1

    if measure_operator == :X
        return state, X(state, i), cstτ, coef

    elseif measure_operator == :ZZ
        if 1 <= i <= N-1
            eigenvalue = zz_eigen(i, i+1)
        elseif model.pbc && i == N
            eigenvalue = (state & 1) == (state >> (N-1) & 1) ? 1 : -1
        end
        return state, state, cstτ + coef * eigenvalue, 0.0

    elseif measure_operator ∈ [:reset, :Z]
        # :reset/:Z has flipped sign convention for coef
        coef_z = sign ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
        if τ >= 1e2
            coef_z = sign ? -0.5 : 0.5
        end
        eigenvalue = (state[N - i + 1] == 0) ? 1 : -1
        return state, state, cstτ + coef_z * eigenvalue, 0.0
    end
end

function _apply_result(model::AnyonModel{OBFAnyon}, τ::Float64, state::T, i::Int, sign::Bool) where {T}
    measure_operator = model.measure_operator
    @assert measure_operator in [:XZZ, :ZZX, :ZZ, :X] "measure_operator must be :XZZ, :ZZX, :ZZ, :X"

    N = model.N
    fl = bmask(T, N)
    X(state, i) = flip(state, fl >> (i-1))
    
    # Common coefficients for all operators, here in constrast to Ising case the sign convention is inverse.
    if τ >= 1e2
        cstτ = 0.5
        coef = sign ? 0.5 : -0.5
    else
        cstτ = cosh(τ/2) / √(2cosh(τ))
        coef = sign ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
    end

    # Helper: get ZZ eigenvalue for sites (j1, j2)
    zz_eigen(j1, j2) = ((state >> (N - j1)) & 1) == ((state >> (N - j2)) & 1) ? 1 : -1
    
    if model.pbc
        i1, i2 = mod1(i + 1, N), mod1(i + 2, N)
    else
        @assert 1 <= i <= N - 2 "Index i must be in [1, N-2] for OBC (OBF)"
        i1, i2 = i + 1, i + 2
    end

    if measure_operator ∈ [:ZZ, :X]
        new_model = AnyonModel(IsingAnyon(), N; pbc=model.pbc, measure_operator=measure_operator)
        return _apply_result(new_model, τ, state, i, sign)
    elseif measure_operator == :XZZ
        # return XZZ components, and coefficients
        return state, X(state, i), cstτ, coef * zz_eigen(i1, i2)
    elseif measure_operator == :ZZX
        # return ZZX components, and coefficients
        return state, X(state, i2), cstτ, coef * zz_eigen(i, i1)
    end
end

function measure_matrix(model::AnyonModel{AT}, τ::Float64, idx::Int, sign::Bool) where {AT<:AbstractAnyonType}

    if model.measure_operator ∈ [:Ferro, :Antiferro]
        @assert model.pbc || (2 <= idx <= model.N-1) "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif model.measure_operator == :ZZ
        @assert model.pbc || (1 <= idx <= model.N-1) "Index idx must be in [1, N-1] for open BC (IsingZZ)"
    elseif model.measure_operator ∈ (:X, :Z, :reset, :resetFibo)
        @assert model.pbc || (1 <= idx <= model.N) "Index idx must be in [1, N] for open BC (IsingX)"
    elseif model.measure_operator ∈ (:XZZ, :ZZX)
        @assert model.pbc || (1 <= idx <= model.N-2) "Index idx must be in [1, N-2] for open BC (OBF)"
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

"""
    measuremap(model::AnyonModel, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool)
    measuremap(model::AnyonModel, ψ::MPS, sites, i::Int, τ::Float64, sign::Bool; 
               cutoff::Float64=1e-10, maxdim::Int=100)

Apply measurement operator to quantum state and return post-measurement state.

# Methods

## Vector version (from `Measurement.jl`)
Apply measurement to a state vector.

### Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Input quantum state vector
- `idx::Int`: Measurement site index
- `sign::Bool`: Measurement outcome (false for +, true for -)

### Returns
- `Vector{ET}`: Post-measurement quantum state (unnormalized)

## MPS version (from `MPSMeasurement.jl`)
Apply measurement to an MPS state.

### Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `ψ::MPS`: Input quantum state
- `sites`: ITensor site indices
- `i::Int`: Measurement site
- `τ::Float64`: Measurement strength parameter
- `sign::Bool`: Measurement outcome (false for +, true for -)
- `cutoff::Float64=1e-10`: MPS truncation cutoff
- `maxdim::Int=100`: Maximum bond dimension

### Returns
- `MPS`: Post-measurement quantum state (normalized)
- `Float64`: Measurement probability
"""
function measuremap(model::AnyonModel{AT}, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool) where {ET, AT<:AbstractAnyonType}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"
    @assert length(state) == length(anyon_basis(model)) "state length is expected to be $(length(anyon_basis(model))), but got $(length(state))"
    return _measuremap_generic(model, τ, state, idx, sign)
end

"""
    measuremap(model::AnyonModel{OBFAnyon}, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool)

Specialized measuremap for OBFAnyon that handles `:OBF` operator (XZZ + ZZX combined).

When `measure_operator == :OBF`, this applies both XZZ and ZZX measurements sequentially
at site `idx`, which is more efficient than two separate layer applications.
"""
function measuremap(model::AnyonModel{OBFAnyon}, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool) where {ET}
    if model.measure_operator == :OBF
        # OBF = XZZ followed by ZZX at the same site
        # This combines two measurements into one layer for efficiency, because e^(-i τ OBF) = e^(-i τ XZZ) e^(-i τ ZZX) when [XZZ, ZZX] = 0
        model_xzz = AnyonModel(OBFAnyon(), model.N; pbc=model.pbc, measure_operator=:XZZ, model.params...)
        model_zzx = AnyonModel(OBFAnyon(), model.N; pbc=model.pbc, measure_operator=:ZZX, model.params...)
        
        # Apply XZZ first, then ZZX
        state = measuremap(model_xzz, τ, state, idx, sign)
        state = measuremap(model_zzx, τ, state, idx, sign)
        return state
    else
        # For :XZZ, :ZZX, :X, :ZZ, use the generic implementation
        return _measuremap_generic(model, τ, state, idx, sign)
    end
end

# Generic implementation extracted for reuse
function _measuremap_generic(model::AnyonModel{AT}, τ::Float64, state::Vector{ET}, idx::Int, sign::Bool) where {ET, AT<:AbstractAnyonType}
    if model.measure_operator ∈ [:Ferro, :Antiferro]
        @assert model.pbc || (2 <= idx <= model.N-1) "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif model.measure_operator == :ZZ
        @assert model.pbc || (1 <= idx <= model.N-1) "Index idx must be in [1, N-1] for open BC (ZZ)"
    elseif model.measure_operator ∈ (:X, :Z, :reset, :resetFibo)
        @assert model.pbc || (1 <= idx <= model.N) "Index idx must be in [1, N] for open BC (X)"
    elseif model.measure_operator ∈ (:XZZ, :ZZX)
        @assert model.pbc || (1 <= idx <= model.N-2) "Index idx must be in [1, N-2] for open BC (OBF)"
    else
        error("Unknown measure operator: $(model.measure_operator)")
    end
    
    basis = anyon_basis(model)
    l = length(basis)
    mapped_state = zeros(ET, l)
    
    for i in 1:l
        outputstate1, outputstate2, output1, output2 = measure_basismap(model, τ, basis[i], idx, sign)
        
        if output2 == 0
            mapped_state[i] += output1 * state[i]
            # outputstate1 is the same as basis[i]
        else
            j2 = searchsortedfirst(basis, outputstate2)
            mapped_state[i] += output1 * state[i]
            # outputstate1 is the same as basis[i]
            mapped_state[j2] += output2 * state[i]
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
            new_trajectory_p = [current_trajectory; false]
            new_prob_p = current_prob * prob_p
            
            push!(next_level_states, normalized_state_p)
            push!(next_level_trajectories, new_trajectory_p)
            push!(next_level_probabilities, new_prob_p)

            state_after_m = measuremap(model, τ, state, site, true)
            prob_m = real(dot(state_after_m, state_after_m))

            normalized_state_m = state_after_m / sqrt(prob_m)
            new_trajectory_m = [current_trajectory; true]
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


"""
    measurement_tree_visualization(trajectories::Vector{Vector{Int64}}, probabilities::Vector{Float64})

Visualize a measurement tree given trajectories and their probabilities.

# Arguments
- `trajectories::Vector{Vector{Int64}}`: Measurement outcome sequences at each branch.
- `probabilities::Vector{Float64}`: Corresponding probabilities for each trajectory.

# Behavior
- Normalizes probabilities to sum to 1.
- Prints levels from root to leaves with indentation representing depth.
"""
function measurement_tree_visualization(trajectories::Vector{Vector{Int64}}, probabilities::Vector{Float64})
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

function _obtain_measurement_config(model::AnyonModel{FibonacciAnyon}, layer_idx::Int, τ::Float64=1.0)
        measurement_sites = iseven(layer_idx) ? collect(1:2:model.N) : collect(2:2:model.N) # measurement sites for Fibonacci different layer, even layers measure at odd sites, odd layers measure at even sites,
        measure_operator = :Antiferro
        measure_anyon_model = AnyonModel(FibonacciAnyon(), model.N; pbc = model.pbc, measure_operator = measure_operator, model.params...)
        measure_strength = τ
    return measurement_sites, measure_anyon_model, measure_strength
end

function _obtain_measurement_config(model::AnyonModel{IsingAnyon}, layer_idx::Int, τ::Float64=1.0)
    # measure at all sites!!!
    measurement_sites = collect(1:model.N)
    measure_operator = iseven(layer_idx) ? :ZZ : :X # measurement sites for Ising each layer, odd layers: measure X, even layers: measure ZZ, but actual Trotterization pattern is √ZZ X √ZZ X.
    measure_anyon_model = AnyonModel(IsingAnyon(), model.N; pbc = model.pbc, measure_operator = measure_operator, model.params...)
    measure_strength = τ
    return measurement_sites, measure_anyon_model, measure_strength
end

function _obtain_measurement_config(model::AnyonModel{OBFAnyon}, layer_idx::Int, τ::Float64=1.0)
    # OBF 8-layer period structure:
    # Layer 1: X (all sites)
    # Layer 2: √ZZ (all sites)  
    # Layer 3,4,5: √OBF on 3 groups (sites 1,4,7...), (sites 2,5,8...), (sites 3,6,9...)
    # Layer 6,7: √OBF reverse order
    # Layer 8: √ZZ (all sites)
    # Final: √X
    # √X, √OBF₁, √OBF₂, OBF₃, √OBF₂, √OBF₁, √X, ZZ

    phase = mod1(layer_idx, 8)
    λ = get_interaction_param(model, :λ, 1.0)  # OBF coupling strength
    if phase == 1 || phase == 7
        # X measurement at all bonds
        measurement_sites = collect(1:model.N)
        measure_operator = :X
        measure_strength = τ/2
    elseif phase == 2 || phase == 6
        # √OBF₁: group index = (1, 4, 7)
        measurement_sites = collect(1:3:model.N)
        measure_operator = :OBF
        measure_strength = λ*τ/2
    elseif phase == 3 || phase == 5
        # √OBF₂: group index = (2, 5, 8)
        measurement_sites = collect(2:3:model.N)
        measure_operator = :OBF
        measure_strength = λ*τ/2
    elseif phase == 4
        # OBF₃: group index = (3, 6, 9)
        measurement_sites = collect(3:3:model.N)
        measure_operator = :OBF
        measure_strength = λ*τ
    elseif phase == 8
        # ZZ measurement at all sites
        measurement_sites = collect(1:model.N)
        measure_operator = :ZZ
        measure_strength = τ
    end
    
    measure_anyon_model = AnyonModel(OBFAnyon(), model.N; pbc = model.pbc, measure_operator = measure_operator, model.params...)
    return measurement_sites, measure_anyon_model, measure_strength
end

"""
    _get_sample_column_indices(model::AnyonModel, layer_idx::Int) -> Vector{Int}

Get the column indices in the samples BitMatrix for a given layer.

For models where different layers have different numbers of measurement sites,
this maps the measurement sites to fixed column positions in the samples matrix.

# Returns
- Vector of column indices where this layer's samples should be stored/read
"""
function _get_sample_column_indices(model::AnyonModel{FibonacciAnyon}, layer_idx::Int)
    # Fibonacci: alternating even/odd sites, always N÷2 measurements
    # Columns 1:(N÷2) are used for all layers
    return collect(1:(model.N ÷ 2))
end

function _get_sample_column_indices(model::AnyonModel{IsingAnyon}, layer_idx::Int)
    # Ising: all N sites measured each layer
    return collect(1:model.N)
end

function _get_sample_column_indices(model::AnyonModel{OBFAnyon}, layer_idx::Int)
    # OBF: different layers measure different sites, but all map to columns 1:N
    # The column index equals the site index being measured
    phase = mod1(layer_idx, 8)
    N = model.N
    
    if phase == 1 || phase == 7
        # X: all sites 1:N → columns 1:N
        return collect(1:N)
    elseif phase == 2 || phase == 6
        # OBF group 1: sites 1,4,7,... → columns 1,4,7,...
        return collect(1:3:N)
    elseif phase == 3 || phase == 5
        # OBF group 2: sites 2,5,8,... → columns 2,5,8,...
        return collect(2:3:N)
    elseif phase == 4
        # OBF group 3: sites 3,6,9,... → columns 3,6,9,...
        return collect(3:3:N)
    elseif phase == 8
        # ZZ: all sites 1:N → columns 1:N
        return collect(1:N)
    end
end

"""
    _samples_per_layer(model::AnyonModel) -> Int

Return the number of sample columns needed per layer (maximum across all layer types).
"""
_samples_per_layer(model::AnyonModel{FibonacciAnyon}) = model.N ÷ 2
_samples_per_layer(model::AnyonModel{IsingAnyon}) = model.N
_samples_per_layer(model::AnyonModel{OBFAnyon}) = model.N  # Max of all layer types


struct Measurement_outcome_bulk{ET}
    states::Vector{ET}
    samples::BitMatrix
    free_energys::Vector{Float32}
end

struct Measurement_outcome_boundary{T}
    state::Vector{T}
    sample::BitVector
    free_energy::Float32
end

"""
Configuration struct for measurement evolution parameters.

# Fields
- `τ::Float64`: Measurement strength parameter
- `t₂::Int`: 2*Number of measurement layers (time steps)
- `rng::MersenneTwister`: Random number generator (default: `MersenneTwister()`)
- `mode::Symbol`: Sampling mode, one of `:sample`, `:Born` (default: `:sample`)
- `t₁::Int`: Starting layer index for evolution (default: 1)
- `verbose::Bool`: Verbosity flag for detailed output (default: false)
- `enable_τ_eff::Bool`: Whether to enable half-strength measurement for the last layer (default: true)
- `λ::Float64`: O'Brien-Fendley coupling strength (default: 0.0, pure Ising when λ=0)
"""
Base.@kwdef struct MeasureConfig
    τ::Float64
    t₂::Int
    rng::MersenneTwister  = MersenneTwister()
    mode::Symbol = :sample
    t₁::Int = 1
    verbose::Bool = false
    enable_τ_eff::Bool = true
    x₂::Int = 1
    x₁::Int = 1
    cutoff::Float64 = 1e-12
    maxdim::Int = 1000
end

measurement_num(::FibonacciAnyon) = 1
measurement_num(::IsingAnyon) = 2
measurement_num(::OBFAnyon) = 2

"""
    layers_per_period(anyon_type) -> Int

Return the number of measurement layers per evolution period.
- Fibonacci: 2 layers
- Ising: 2 layers (X, ZZ)
- OBF: 8 layers (√X, √OBF₁, √OBF₂, OBF₃, √OBF₂, √OBF₁, √X, ZZ), here OBF represents XZZ + ZZX. At the end plus a final √ZZ layer.
"""
layers_per_period(::FibonacciAnyon) = 2
layers_per_period(::IsingAnyon) = 2
layers_per_period(::OBFAnyon) = 8

"""
    boundary_evolution(model::AnyonModel, state::Vector{T}, measure_config::MeasureConfig, 
                       sample::Union{Nothing, BitVector}=nothing; layer_idx::Int=1)
    boundary_evolution(model::AnyonModel, sites, state::MPS, sample::BitVector, 
                       measure_config::MeasureConfig; 
                       cutoff::Float64=1e-12, maxdim::Int=1000, layer_idx::Int=1)

Evolve a quantum state under a single layer of boundary measurements with a given trajectory.

# Methods

## Vector version (from `Measurement.jl`)
Evolve a state vector under boundary measurements.

### Arguments
- `model::AnyonModel`: The anyon model specifying system parameters.
- `state::Vector{T}`: The initial state vector.
- `measure_config::MeasureConfig`: Configuration containing measurement strength `τ` and mode.
- `sample::Union{Nothing, BitVector}=nothing`: The measurement outcomes for the layer.
- `layer_idx::Int=1`: The layer index for measurement configuration.

### Returns
- `Measurement_outcome_boundary`: A struct containing `state`, `sample`, and `free_energy`.

## MPS version (from `MPSMeasurement.jl`)
Evolve an MPS state under boundary measurements.

### Arguments
- `model::AnyonModel`: The anyon model specifying system parameters.
- `sites`: The ITensor site indices.
- `state::MPS`: The initial MPS state.
- `sample::BitVector`: The measurement outcomes for the layer.
- `measure_config::MeasureConfig`: Configuration containing measurement strength `τ` and mode (`:sample` or `:Born`).
- `cutoff::Float64=1e-12`: MPS truncation cutoff.
- `maxdim::Int=1000`: Maximum bond dimension for MPS operations.
- `layer_idx::Int=1`: The layer index for measurement configuration.

### Returns
- `Measurement_outcome_mps_boundary`: A struct containing:
  - `state::MPS`: The evolved MPS state.
  - `samples::BitVector`: The measurement outcomes for the layer.
  - `free_energy::Float64`: The free energy associated with the measurement layer.
"""
function boundary_evolution(anyon_model::AnyonModel{AT}, state::Vector{T}, measure_config::MeasureConfig, 
    sample::Union{Nothing, BitVector}=nothing; layer_idx::Int=1) where{T, AT<:AbstractAnyonType}
    
    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    τ_eff = measure_config.enable_τ_eff ? measure_config.τ / 2 : measure_config.τ
    if measure_config.mode == :sample
        N = anyon_model.N
        size(sample, 1) == measurement_num(anyon_model.anyon_type)*(N ÷ 2) || error("sample size mismatch with anyon_model $(N)")
        return _apply_measurement_layer(anyon_model, τ_eff, state, sample; layer_idx=layer_idx)
    elseif measure_config.mode == :Born
        return _sample_layer(anyon_model, τ_eff, state; layer_idx=layer_idx, rng=measure_config.rng, verbose=measure_config.verbose)
    end
end

"""
    _apply_measurement_layer(model::AnyonModel, τ::Float64, state::Vector{T}, 
                              layer_sample::BitVector; layer_idx::Int=1) where {T}

Apply deterministic measurements to a layer with given measurement outcomes.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `state::Vector{T}`: Quantum state vector
- `layer_sample::BitVector`: Measurement outcomes for the layer
- `layer_idx::Int=1`: Layer index (1-based) to determine measurement pattern

# Returns
- `Measurement_outcome_boundary`: A struct containing the post-measurement state, sample, and total free energy.
"""
function _apply_measurement_layer(anyon_model::AnyonModel{AT}, τ::Float64, state::Vector{T}, 
    layer_sample::BitVector; layer_idx::Int64=1) where {T, AT<:AbstractAnyonType}
    # Helper function to apply deterministic measurements to a layer, connect measure on each site together.

    total_free_energy = zero(real(T))

    measurement_sites, measure_anyon_model, measurement_strength = _obtain_measurement_config(anyon_model, layer_idx, τ)  

    for (idx, sign) in enumerate(layer_sample)
        # Apply measurement at site measurement_sites[idx] with outcome sign
        state = measuremap(measure_anyon_model, measurement_strength, state, measurement_sites[idx], sign)
        prob = real(dot(state, state))
        total_free_energy += -log(prob)
        state = state ./ sqrt(prob)
    end

    return Measurement_outcome_boundary(state, layer_sample, Float32(total_free_energy))
end

"""
    _sample_layer(model::AnyonModel, τ::Float64, state::Vector{T};
                   layer_idx::Int=1, rng::MersenneTwister=MersenneTwister(), 
                   verbose::Bool=false) where {T}

Perform random measurement on a layer using Born rule sampling.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `τ::Float64`: Measurement strength parameter
- `state::Vector{T}`: Quantum state vector
- `layer_idx::Int=1`: Layer index (1-based) to determine measurement pattern
- `rng::MersenneTwister`: Random number generator
- `verbose::Bool=false`: Whether to print debug information

# Returns
- `Measurement_outcome_boundary`: A struct containing the post-measurement state, sample outcomes, and free energy.
"""
function _sample_layer(anyon_model::AnyonModel{AT}, τ::Float64, state::Vector{T};
    layer_idx::Int64=1,
    rng::MersenneTwister = MersenneTwister(),
    verbose::Bool=false) where {T, AT<:AbstractAnyonType}

    measurement_sites, measure_anyon_model, measurement_strength = _obtain_measurement_config(anyon_model, layer_idx, τ)  
    n = length(measurement_sites)
    sample = BitVector(zeros(Bool, n))
    F_layer = 0.0


    for (i, site) in enumerate(measurement_sites)
        # first 0 branch
        ψ0 = measuremap(measure_anyon_model, measurement_strength, state, site, false)
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
            ψ1 = measuremap(measure_anyon_model, measurement_strength, state, site, true)
            sample[i] = 1
            state = ψ1 ./ sqrt(p1)
            F_layer += -log(p1)
            verbose && @show -log(p1)
        end
    end
    return Measurement_outcome_boundary(state, sample, Float32(F_layer))
end

"""
    bulk_evolution(model::AnyonModel, state::Vector{ET}, measure_config::MeasureConfig,
                   samples::Union{Nothing,BitMatrix}=nothing)
    bulk_evolution(model::AnyonModel, sites, state::MPS, measure_config::MeasureConfig,
                   samples::Union{Nothing,BitMatrix}=nothing;
                   cutoff::Float64=1e-10, maxdim::Int=100)

Perform bulk measurement evolution from t₁ to t₂ on quantum state.

# Methods

## Vector version (from `Measurement.jl`)
Evolve a state vector under bulk measurements.

### Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `state::Vector{ET}`: Initial quantum state vector
- `measure_config::MeasureConfig`: Configuration struct containing τ, t₁, t₂, mode, rng, etc.
- `samples::Union{Nothing,BitMatrix}=nothing`: Predefined measurement sequences for `:sample` mode

### Returns
- `Measurement_outcome_bulk`: A struct containing:
  - `states::Vector{Vector{ET}}`: Intermediate states at each time step
  - `samples::BitMatrix`: Measurement outcome sequences
  - `free_energys::Vector{Float64}`: Free energy for each layer

## MPS version (from `MPSMeasurement.jl`)
Evolve an MPS state under bulk measurements.

### Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `sites`: ITensor site indices
- `state::MPS`: Initial MPS quantum state
- `measure_config::MeasureConfig`: Configuration struct containing τ, t₁, t₂, mode, rng, etc.
- `samples::Union{Nothing,BitMatrix}=nothing`: Predefined measurement sequences for `:sample` mode
- `cutoff::Float64=1e-10`: MPS truncation cutoff
- `maxdim::Int=100`: Maximum bond dimension

### Returns
- `Measurement_outcome_bulk`: A struct containing:
  - `states::Vector{MPS}`: Intermediate states at each time step
  - `samples::BitMatrix`: Measurement outcome sequences
  - `free_energy::Vector{Float64}`: Free energy for each layer

# Notes
- In `:Born` mode, samples are generated probabilistically via Born rule
- In `:sample` mode, `samples` must be provided as input
- (2N+1) layers of measurements correspond to N time steps of evolution
"""
function bulk_evolution(anyon_model::AnyonModel{AT},   # DRY: don't repeat yourself.
                                    state::Vector{ET},
                                    measure_config::MeasureConfig,
                                    samples::Union{Nothing,BitMatrix}=nothing) where {ET, AT<:AbstractAnyonType}
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

    n_cols = _samples_per_layer(model)  # Use max samples per layer
    Δt = measure_config.t₂ - measure_config.t₁ + 1
    τ = measure_config.τ
    enable_τ_eff = measure_config.enable_τ_eff
    rng = measure_config.rng
    verbose = measure_config.verbose
    Δt >= 0 || error("t₂ must be >= t₁")
    
    n_layers = layers_per_period(model.anyon_type)
    D = Δt * n_layers  # total number of layers

    # 1. Initialize sample matrix with max columns per layer
    samples = BitMatrix(zeros(Bool, D, n_cols))
    sample_free_energy = zeros(Float32, D)
    states = Vector{Vector{ET}}(undef, Δt)

    for period in 1:Δt
        # Apply all layers in this period
        for layer in 1:n_layers
            global_layer_idx = (period - 1) * n_layers + layer
            # Apply τ_eff only on the last layer of the last period
            τ_current = (period == Δt && layer == n_layers && enable_τ_eff) ? τ/2 : τ
            
            outcome = _sample_layer(model, τ_current, current_state; 
                                    layer_idx=global_layer_idx, rng=rng, verbose=verbose)
            current_state = outcome.state
            
            # Write samples to correct column indices for this layer
            col_indices = _get_sample_column_indices(model, global_layer_idx)
            samples[global_layer_idx, col_indices] = outcome.sample
            sample_free_energy[global_layer_idx] = outcome.free_energy
        end
        states[period] = current_state
    end

    return Measurement_outcome_bulk(states, samples, sample_free_energy)
end

function _sample_measure(model::AnyonModel{AT}, current_state::Vector{ET}, samples::BitMatrix, measure_config::MeasureConfig) where {AT, ET}

        n_cols = _samples_per_layer(model)  # Use max samples per layer
        Δt = measure_config.t₂ - measure_config.t₁ + 1
        τ = measure_config.τ
        enable_τ_eff = measure_config.enable_τ_eff
        Δt >= 0 || error("t₂ must be >= t₁")
        
        n_layers = layers_per_period(model.anyon_type)
        D = Δt * n_layers  # total number of layers

        sample_free_energy = zeros(Float32, D)
        states = Vector{Vector{ET}}(undef, Δt)
        
        # 2. Validate sample matrix dimensions
        size(samples) == (D, n_cols) || error("sample size should be ($D, $n_cols), got $(size(samples))")

        # 3. Deterministic trajectory for modes :sample
        #  Fibonacci: 2-layer period (even sites, odd sites)
        #  Ising (λ=0): 2-layer period (ZZ, X)
        #  OBF (λ≠0):   8-layer period (ZZ, X, OBF, X)

         #  If measure_operator is :Fibo, the measurement sites are half of N, circuits belike:
        #   1   1   1   1   1   1   1   1   1
        #     1   1   1   1   1   1   1   1    
        #   1   1   1   1   1   1   1   1   1
        #   -τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ-τ- (head tail concatenation)

        # If measure_operator is :X or :ZZ, the measurement sites are N, circuits belike:
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
            # √M₁ᵒ M₁ᵉ √M₁ᵒ √M₁ᵒ M₁ᵉ √M₁ᵒ ⋯ √M₁ᵒ M₁ᵉ √M₁ᵒ→ M₁ᵉ M₁ᵒ M₁ᵉ M₁ᵒ ⋯ M₁ᵉ √M₁ᵒ. 
            # √X ZZ √X √X ZZ √X ⋯ √X ZZ √X→ √X ZZ X ZZ ⋯ X ZZ √X. To ensure each layer is hermitian, first layer doesn't matter.
            # Or √ZZ X √ZZ √ZZ X √ZZ ⋯ X √ZZ X √ZZ→ X ZZ X ZZ ⋯ X √ZZ, also works (we choose this one here).
            for layer in 1:n_layers
                global_layer_idx = (period - 1) * n_layers + layer
                # Apply τ_eff only on the last layer of the last period
                τ_current = (period == Δt && layer == n_layers && enable_τ_eff) ? τ/2 : τ
                
                # Read samples from correct column indices for this layer
                col_indices = _get_sample_column_indices(model, global_layer_idx)
                layer_sample = BitVector(samples[global_layer_idx, col_indices])
                
                outcome = _apply_measurement_layer(
                                model, τ_current, current_state,
                                layer_sample; layer_idx=global_layer_idx)
                current_state = outcome.state
                sample_free_energy[global_layer_idx] = outcome.free_energy
            end
            states[period] = current_state
        end
    return Measurement_outcome_bulk(states, samples, sample_free_energy)
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