"""
    measure_basismap(::Type{T}, τ::Float64, state::T, i::Int, sign::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}

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
    
    if anyon_type == :Fibo
        
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
                return state, cstτ+coef
            elseif state & mask == str100
                return state, cstτ+coef
            elseif state & mask == str101
                return state, cstτ-coef
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
                    return state, cstτ+coef
                elseif state & mask == str100
                    return state, cstτ+coef
                elseif state & mask == str101
                    return state, cstτ-coef
                end
            elseif i == N #count from the left
            mask=bmask(T, N, 2, 1)
            str100, str101, str010, str001, str000 = bmask(T,2), bmask(T, N, 2), bmask(T, 1), bmask(T, N), T(0)
                if state & mask == str000
                    return state, X(state,i), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2)
                elseif state & mask == str010
                    return state, X(state,i), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2)
                elseif state & mask == str001
                    return state, cstτ+coef
                elseif state & mask == str100
                    return state, cstτ+coef
                elseif state & mask == str101
                    return state, cstτ-coef
                end
            end
        end
    elseif anyon_type == :IsingX
        if τ >= 1e2
            cstτ = 0.5
            coef = sign == 0 ? 0.5 : -0.5
        else
            cstτ = cosh(τ/2) / √(2cosh(τ))
            coef = sign == 0 ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
        end

        return state, X(state,i), cstτ, coef

    elseif anyon_type == :IsingZZ
        if τ >= 1e2
            cstτ = 0.5
            coef = sign == 0 ? 0.5 : -0.5
        else
            cstτ = cosh(τ/2) / √(2cosh(τ))
            coef = sign == 0 ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
        end

        if 1<= i <= N-1
            if ((state >> (N - i)) & 1) == ((state >> (N - i -1)) & 1)
                return state, cstτ+coef
            else
                return state, cstτ-coef
            end
        end

        if pbc && i == N
            if (state & 1) == (state >> (N-1) & 1)
                return state, cstτ+coef
            else
                return state, cstτ-coef
            end
        end
    elseif (anyon_type ∈ (:reset, :resetFibo) && τ >= 1e2)|| anyon_type == :IsingZ
        if τ >= 1e2
            cstτ = 0.5
            coef = sign == 0 ? 0.5 : -0.5
        else
            cstτ = cosh(τ/2) / √(2cosh(τ))
            coef = sign == 0 ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
        end

        return state, (state[N - i + 1] == 0) ? cstτ + coef : cstτ - coef
    else
        error("Unknown measure class: $anyon_type")
    end
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
        outcome = measure_basismap(T, τ, basis[i], idx, sign, pbc; anyon_type = anyon_type)

        if length(outcome) == 4
            s1, s2, w1, w2 = outcome
            j2 = searchsortedfirst(basis, s2)
            Bmatrix[i, i] += w1
            Bmatrix[i, j2] += w2
        else
            s, w = outcome
            Bmatrix[i, i] += w
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
        output = measure_basismap(T, τ, basis[i], idx, sign, pbc, anyon_type=anyon_type)
        if length(output) == 4
            outputstate1, outputstate2, output1, output2=output
            j2=searchsortedfirst(basis, outputstate2)
            mapped_state[i]+=output1*state[i] # outputstate1 is the same as basis[i]
            mapped_state[j2]+=output2*state[i]
        else
            outputstate, output1=output # outputstate is the same as basis[i]
            mapped_state[i]+=output1*state[i]
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
function measurement_enumeration(::Type{T}, τ::Float64, initial_state::Vector{ET}, measurement_sites::Vector{Int}, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}
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

"""
    apply_measurement_layer!(N::Int64, state::Vector{T}, τ::Float64, layer_sample::Vector{Int64}, layer_idx::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo)

Generate measurement samples at 1 layer with post_selection outcomes, i.e., with time step 1.

# Arguments
- `N::Int`: Chain length N
- `state::Vector{ET}`: Initial quantum state vector
- `τ::Float64`: Measurement strength parameter
- `layer_sample::Vector{Int64}`: Measurement outcomes for the layer
- `layer_idx::Int64`: Layer index (1-based) to determine measurement pattern
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
  (post-measurement states, measurement sequences, free energies)

Samples measurement outcomes with given measurement outcomes.
"""
# Helper function to apply measurements to a layer
function apply_measurement_layer!(N::Int64, state::Vector{T}, τ::Float64, layer_sample::Vector{Int64}, layer_idx::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo, return_free_energy::Bool=false) where {T}
    total_free_energy = zero(real(T))
    if anyon_type == :Fibo
        if layer_idx % 2 == 1
            measurement_sites = collect(2:2:N)  # odd sites anyons, even sites qubits
        else
            measurement_sites = collect(1:2:N)  # even sites anyons, odd sites qubits
        end
        for (idx, sign) in enumerate(layer_sample)
            state = measuremap(N, τ, state, measurement_sites[idx], sign, pbc, anyon_type = anyon_type)
            prob_p = real(dot(state, state))
            total_free_energy += -log(prob_p)

            state ./= sqrt(prob_p)
        end

        return return_free_energy ? (total_free_energy,  state) :  state

    elseif anyon_type ∈ (:IsingX, :IsingZZ, :IsingZ)
        # measure at all sites!!!
        measurement_sites = collect(1:N)
        if layer_idx % 2 == 1
            # odd layers: measure X
            for (idx, sign) in enumerate(layer_sample)
                state = measuremap(N, τ, state, measurement_sites[idx], sign, pbc, anyon_type = :IsingX)
                prob_p = real(dot(state, state))
                total_free_energy += -log(prob_p)
                state ./= sqrt(prob_p)
            end
        
        else
            # even layers: measure ZZ
            for (idx, sign) in enumerate(layer_sample)
                state = measuremap(N, τ, state, measurement_sites[idx], sign, pbc, anyon_type = :IsingZZ)
                prob_p = real(dot(state, state))
                total_free_energy += -log(prob_p)
                state ./= sqrt(prob_p)
            end
        end

        return return_free_energy ? (total_free_energy, state) : state
    else
        error("Unknown measure class: $anyon_type")
    end
end

function generate_state(τ::Float64, state::Vector{T}, sample::ET, pbc::Bool=true; temp::Bool=false, anyon_type::Symbol=:Fibo, return_free_energy::Bool=false) where{T, ET}

    if ET == Vector{Int} # single layer
        if anyon_type == :Fibo
            N = 2 * length(sample)
        else
            N = length(sample)
        end
        return apply_measurement_layer!(N, state, τ, sample, 1, pbc, anyon_type=anyon_type, return_free_energy=return_free_energy)

    elseif ET == Matrix{Int} # multiple layers
        D = size(sample, 1)
        N = (anyon_type == :Fibo) ? 2 * size(sample, 2) : size(sample, 2)

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
        statelis = temp ? Vector{Vector{T}}(undef, D) : nothing
        layer_free_energy = zeros(real(T), D)  # free energy of each layer

        for layer in 1:D
            # √M₁ᵉ √M₁ᵒ √M₁ᵉ √M₁ᵉ √M₁ᵒ √M₁ᵉ ⋯ √M₁ᵉ √M₁ᵒ √M₁ᵉ→ M₁ᵉ M₁ᵒ M₁ᵉ M₁ᵒ ⋯ M₁ᵉ M₁ᵒ √M₁ᵉ. 
            # √X √ZZ √X √X √ZZ √X ⋯ √X √ZZ √X→ X ZZ X ZZ ⋯ X ZZ √X. To ensure each layer is hermitian, first layer is doesn't matter.
            τ_eff = (layer == D) ? τ/2 : τ

            if return_free_energy
                state_measured, ΔF = apply_measurement_layer!(
                                N, state, τ_eff,
                                sample[layer, :], layer, pbc;
                                anyon_type=anyon_type,
                                return_free_energy=true)
                layer_free_energy[layer] += ΔF
            else
                state_measured = apply_measurement_layer!(
                            N, state, τ_eff,
                            sample[layer, :], layer, pbc;
                            anyon_type=anyon_type,
                            return_free_energy=false)
            end

            temp && (statelis[layer] = state_measured)
        end

        # 统一返回
        if temp
            return return_free_energy ? (statelis, layer_free_energy) : statelis
        else
            return return_free_energy ? (state_measured, layer_free_energy) : state_measured
        end
    else
        error("Unsupported sample type: $ET")
    end
end

# return (state_measured, sample, free_energy) of one layer
function _sample_layer!(state, τ_eff, sites, rng, N, pbc, anyon_type)
    n = length(sites)
    sample = zeros(Int, n)
    F_layer = 0.0

    for (i, site) in enumerate(sites)
        # first 0 branch
        ψ0 = measuremap(N, τ_eff, state, site, 0, pbc, anyon_type = anyon_type)
        p0 = real(dot(ψ0, ψ0))
        p1 = 1 - p0

        if rand(rng) < p0
            sample[i] = 0
            state_measured = ψ0 ./ sqrt(p0)
            F_layer += -log(p0)
        else
            # else 1 branch
            ψ1 = measuremap(N, τ_eff, state, site, 1, pbc, anyon_type = anyon_type)
            sample[i] = 1
            state_measured = ψ1 ./ sqrt(p1)
            F_layer += -log(p1)
        end
    end
    return state_measured, sample, F_layer
end

"""
    measure_evolution!(state, τ, sites, D=1;
             rng = Random.default_rng(),
             pbc = true,
             anyon_type = :Fibo,
             mode = :prob,      # :prob random sampling, :post0 all zeros, :post1 all ones, :sample given sample
             temp = false)     # return all intermediate states

Perform measurement evolution on the initial state with specified parameters, returning final states, measurement sequences, and free energies, with options for different sampling modes and temporal settings.
# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `D::Int=1`: Number of measurement layers (time steps)
- `rng`: Random number generator (default: `Random.default_rng()`)
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type
- `mode::Symbol=:prob`: Sampling mode, one of `:prob`, `:post0`, `:post1`, `:sample`, `:Born`
- `temp::Bool=false`: If true, return all intermediate states; if false, return only final state
- `sample_mat::Union{Nothing,Matrix{Int}}=nothing`: Predefined measurement sequences for `:sample` mode

# Returns
- `Tuple{Union{Vector{ET}, Vector{Vector{ET}}}, Vector{Vector{Int}}, Vector{Float64}}`: 
  (final state or all intermediate states, measurement sequences, free energies)
- states        : final state or all intermediate states (depends on temp)
- samples       : measurement sequences (Vector{Vector{Int}} or Matrix{Int})
- free_energies : cumulative free energies for each layer (or each sample)
"""
function measure_evolution!(N::Int,
                  τ::Float64,
                  state::Vector{ET},
                  D::Int = 1;
                  rng = Random.default_rng(),
                  pbc::Bool = true,
                  anyon_type::Symbol = :Fibo,
                  mode::Symbol = :prob,
                  temp::Bool = false,
                  sample_mat::Union{Nothing,Matrix{Int}}=nothing) where {ET}

    n_measure = anyon_type == :Fibo ? N÷2 : N

    # ---------- Sample_mat decided according to mode ----------
    if mode == :Born
        sample_mat = zeros(Int, D, n_measure)   # filled in during sampling
    elseif mode == :post0
        sample_mat = zeros(Int, D, n_measure)
    elseif mode == :post1
        sample_mat = ones(Int, D, n_measure)
    elseif mode == :sample
        isnothing(sample_mat) && error("When mode=:sample sample_mat must be ::Matrix{Int}")
        size(sample_mat) == (D, n_measure) ||
            error("sample_mat size should be ($D, $n_measure)")
    else
        error("mode must be ∈ [:Born, :post0, :post1, :sample]")
    end

    # 2. Born trajectory (only for :Born)
    if mode == :Born
        current_state = copy(state)
        free_energies = zeros(D)
        states = temp ? Vector{Vector{ET}}(undef, D) : nothing

        for layer in 1:D
            τ_eff = (layer == D) ? τ / 2 : τ

            # Random sampling for this layer
            sites = anyon_type == :Fibo ? collect(2:2:N) : collect(1:N) # measurement sites for this layer
            seq = _sample_layer!(current_state, τ_eff, sites, rng,
                                 N, pbc, anyon_type)

            sample_mat[layer, :] .= seq[1]         
            free_energies[layer] = seq[2]      
            if temp
                states[layer] = copy(current_state)
            end
        end
        final_state = temp ? states : current_state
    else
        # :post0, :post1, :sample, directly deterministic trajectory
        final_state, free_energy = generate_state(
            τ, state, sample_mat, pbc;
            anyon_type = anyon_type,
            return_free_energy = true,
            temp = temp)
        free_energies = temp ? free_energy : [free_energy]
    end

    # 3. Return in same format
    samples = [sample_mat[layer, :] for layer in 1:D]
    return final_state, samples, free_energies
end

"""
    Boundary_measure(::Type{T}, τ::Float64, state::Vector{ET}, measurement_sites::Vector{Int}, num_samples::Int=1000, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}

Generate measurement samples at boundary sites with probabilistic outcomes, i.e., without time axis evolution.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `measurement_sites::Vector{Int}`: Sites to perform measurements
- `num_samples::Int=1000`: Number of measurement samples to generate
- `rng::MersenneTwister`: Random number generator
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
  (post-measurement states, measurement sequences, free energies)

Samples measurement outcomes probabilistically based on Born rule.
"""
function Boundary_measure(::Type{T}, τ::Float64, state::Vector{ET}, measurement_sites::Vector{Int}, num_samples::Int=1000, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}, ET}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    
    num_sites = length(measurement_sites)
    
    sample_measured_states = Vector{Vector{Float64}}(undef, num_samples)
    samples = Vector{Vector{Int64}}(undef, num_samples)
    sample_free_energy = Vector{Float64}(undef, num_samples)

    for sample_idx in 1:num_samples
        current_sequence = Vector{Int64}(undef, num_sites)
        current_state = copy(state)  
        total_free_energy = 0.0
        
        # meaure from the left to the right
        for (site_idx, measurement_site) in enumerate(measurement_sites)
           
            state_after_p = measuremap(T, τ, current_state, measurement_site, 0, pbc, anyon_type=anyon_type)


            prob_p = real(dot(state_after_p, state_after_p))
            prob_m = 1 - prob_p
            
            random_number = rand(rng)
            if random_number < prob_p
                current_sequence[site_idx] = 0
                current_state = state_after_p ./ sqrt(prob_p)
                total_free_energy += -log(prob_p)
            else
                state_after_m = measuremap(T, τ, current_state, measurement_site, 1, pbc, anyon_type = anyon_type)
                current_sequence[site_idx] = 1
                current_state = state_after_m ./ sqrt(prob_m)
                total_free_energy += -log(prob_m)
            end
        end
        
        sample_measured_states[sample_idx] = current_state
        samples[sample_idx] = current_sequence
        sample_free_energy[sample_idx] = total_free_energy
    end
    
    return sample_measured_states, samples, sample_free_energy
end

Boundary_measure(N::Int, τ::Float64, state::Vector{ET}, measurement_sites::Vector{Int},num_samples::Int=1000, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET} = Boundary_measure(BitStr{N, Int}, τ, state, measurement_sites, num_samples, rng, pbc, anyon_type = anyon_type)

"""
    Boundarypost_selection(N::Int64, τ::Float64, state::Vector{ET}, measurement_sites::Vector{Int}, sign::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo)

Generate measurement samples at boundary sites with post_selection outcomes, i.e., given spatial evolution without time axis evolution.

# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `measurement_sites::Vector{Int}`: Sites to perform measurements
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
  (post-measurement states, measurement sequences, free energies)

Samples measurement outcomes with given measurement outcomes.
"""
function Boundarypost_selection(N::Int64, τ::Float64, state::Vector{ET}, measurement_sites::Vector{Int}, sign::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    num_sites = length(measurement_sites)

    current_sequence = Vector{Int64}(undef, num_sites)
    current_state = copy(state)  
    total_free_energy = 0.0
    
    # meaure from the left to the right
    for (site_idx, measurement_site) in enumerate(measurement_sites)
       
        state_after_measure = measuremap(N, τ, current_state, measurement_site, sign, pbc, anyon_type=anyon_type)

        prob = real(dot(state_after_measure, state_after_measure))

        current_sequence[site_idx] = sign
        current_state = state_after_measure ./ sqrt(prob)
        total_free_energy += -log(prob)
    end
    
    return current_state, current_sequence, total_free_energy
end

"""
    Bulkmeasure(N::Int64, τ::Float64, state::Vector{ET}, D::Int64, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true; anyon_type::Symbol=:Fibo)

Generate measurement samples at bulk sites with probabilistic outcomes, i.e., with time axis evolution, together with spatial evolution axis as bulk.

# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `D::Int64`: Number of measurement layers (depth), or time step (over 2)
- `rng::MersenneTwister`: Random number generator
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
  (post-measurement states, measurement sequences, free energies)

Samples measurement outcomes probabilistically based on Born rule.
"""
function Bulkmeasure(N::Int64, τ::Float64, state::Vector{ET}, D::Int64, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true; anyon_type::Symbol=:Fibo) where {ET}
    
    sample_free_energy = Vector{Float64}(undef, D)
    sample_measured_states = Vector{Vector{ET}}(undef, D)
    
    current_state = copy(state)  
    
    if anyon_type == :Fibo
        sample = zeros(Int, D, div(N,2))
        for layer in 1:D
            current_sequence = zeros(Int, div(N,2))
            total_free_energy = 0.0
            
            # Alternating measurement pattern
            if layer % 2 == 1
                measurement_sites = collect(2:2:N)  # odd layers: even sites
            else
                measurement_sites = collect(1:2:N)  # even layers: odd sites
            end
            
            measurement_τ = (layer == D) ? τ/2 : τ
            # measure :sqrtp is Pi 0/tau, measure :sqrtm is Pi 1, measure 0 is Pi 0/tau, measure 1 is Pi 1
                
            for (site_idx, measurement_site) in enumerate(measurement_sites)
                state_after_p = measuremap(N, measurement_τ, current_state, measurement_site, 0, pbc, anyon_type = anyon_type)

                prob_sqrtp = real(dot(state_after_p, state_after_p))
                prob_sqrtm = 1 - prob_sqrtp
                random_number = rand(rng)  
                if random_number < prob_sqrtp
                    current_sequence[site_idx] = 0
                    current_state = state_after_p ./ sqrt(prob_sqrtp)
                    total_free_energy += -log(prob_sqrtp)
                else
                    state_after_m = measuremap(N, measurement_τ, current_state, measurement_site, 1, pbc, anyon_type = anyon_type)
                    current_sequence[site_idx] = 1
                    current_state = state_after_m ./ sqrt(prob_sqrtm)
                    total_free_energy += -log(prob_sqrtm)
                end
            end
    
            sample_measured_states[layer] = current_state
            sample[layer, :] = current_sequence
            sample_free_energy[layer] = total_free_energy
        end

    elseif anyon_type == :IsingX || anyon_type == :IsingZZ
        sample = zeros(Int, D, N)
        for layer in 1:D
            current_sequence = zeros(Int, N)
            total_free_energy = 0.0
            measurement_sites = collect(1:N)            
            measurement_τ = (layer == D) ? τ/2 : τ
            # Alternating measurement pattern
            if layer % 2 == 1
               # odd layers: even sites, measure X
               for (site_idx, measurement_site) in enumerate(measurement_sites)
                   state_after_p = measuremap(N, measurement_τ, current_state, measurement_site, 0, pbc, anyon_type = :IsingX)

                   prob_sqrtp = real(dot(state_after_p, state_after_p))
                   prob_sqrtm = 1 - prob_sqrtp
                   random_number = rand(rng)  
                   if random_number < prob_sqrtp
                       current_sequence[site_idx] = 0
                       current_state = state_after_p ./ sqrt(prob_sqrtp)
                       total_free_energy += -log(prob_sqrtp)
                   else
                       state_after_m = measuremap(N, measurement_τ, current_state, measurement_site, 1, pbc, anyon_type = :IsingX)
                       current_sequence[site_idx] = 1
                       current_state = state_after_m ./ sqrt(prob_sqrtm)
                       total_free_energy += -log(prob_sqrtm)
                   end
               end
       
               sample_measured_states[layer] = current_state
               sample[layer, :] = current_sequence
               sample_free_energy[layer] = total_free_energy
            else
               # even layers: odd sites, measure ZZ
               for (site_idx, measurement_site) in enumerate(measurement_sites)
                   state_after_p = measuremap(N, measurement_τ, current_state, measurement_site, 0, pbc, anyon_type = :IsingZZ)

                   prob_sqrtp = real(dot(state_after_p, state_after_p))
                   prob_sqrtm = 1 - prob_sqrtp
                   random_number = rand(rng)  
                   if random_number < prob_sqrtp
                       current_sequence[site_idx] = 0
                       current_state = state_after_p ./ sqrt(prob_sqrtp)
                       total_free_energy += -log(prob_sqrtp)
                   else
                       state_after_m = measuremap(N, measurement_τ, current_state, measurement_site, 1, pbc, anyon_type = :IsingZZ)
                       current_sequence[site_idx] = 1
                       current_state = state_after_m ./ sqrt(prob_sqrtm)
                       total_free_energy += -log(prob_sqrtm)
                   end
               end
       
               sample_measured_states[layer] = current_state
               sample[layer, :] = current_sequence
               sample_free_energy[layer] = total_free_energy
            end
        
        end
    else
        error("Unknown measure class: $anyon_type")
    end

    return sample_measured_states, sample, sample_free_energy
end

"""
    Bulkpost_selection(N::Int64, τ::Float64, state::Vector{ET}, D::Int64, sign::Int64, pbc::Bool=true; anyon_type::Symbol=:Fibo)

Generate measurement samples at bulk sites with post_selection outcomes, i.e., with time axis evolution, together with spatial evolution axis as bulk.

# Arguments
- `N::Int`: Chain length N
- `τ::Float64`: Measurement strength parameter
- `state::Vector{ET}`: Initial quantum state vector
- `D::Int64`: Number of measurement layers (depth), or time step (over 2)
- `sign::Int64`: Measurement sign, 0 for positive, 1 for negative
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: anyon type

# Returns
- `Tuple{Vector{Vector{Float64}}, Vector{Vector{Int64}}, Vector{Float64}}`: 
  (post-measurement states, measurement sequences, free energies)

Samples measurement outcomes with given measurement outcomes.
"""
function Bulkpost_selection(
        N::Int64,
        τ::Float64,
        state::Vector{ET},
        D::Int64,
        sign::Int64,
        pbc::Bool=true;
        anyon_type::Symbol=:Fibo
    ) where {ET}

    # 1. Build the sample, all 0 or all 1, each layer N/2 or N measurements, depending on anyon_type, total D layers.
    n_measure_per_layer = (anyon_type == :Fibo) ? N ÷ 2 : N

    sample = sign == 1 ? ones(Int, D, n_measure_per_layer) :
                              zeros(Int, D, n_measure_per_layer)

    # 2. generate_state to run all the layers
    sample_measured_states, sample_free_energy = generate_state(
        τ, state, sample, pbc,
        anyon_type=anyon_type,
        return_free_energy=true,
        temp=true)          # need all the intermediate states, generate_state return each layer free energy, when temp=true return Vector{Vector{ET}}, each layer state    

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