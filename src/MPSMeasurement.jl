"""
MPS-based implementation for Fibonacci chain measurements using ITensor
"""

"""
    initial_mps(N::Int; pbc::Bool=true) -> MPS

Generate the initial state of Fibonacci chain as an MPS.
"""
function initial_mps(N::Int; pbc::Bool=true)
    # Create sites for Fibonacci anyons (using S=1/2 fermions to approximate)
    sites = siteinds("S=1/2", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]
    
    # Create MPS from product state
    ψ0 = randomMPS(sites, state)
    
    return ψ0
end


"""
    measurement_operator_mps(sites, i::Int, τ::Float64, sign::Symbol; pbc::Bool=true) -> MPO

Create measurement operator at site i with parameter τ as an MPO.
"""
function measurement_operator_mps(sites, i::Int, τ::Float64, sign::Symbol; pbc::Bool=true)
    N = length(sites)
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (:p, :m) "sign must be either :p or :m"
    
    # Golden ratio
    ϕ = (1 + √5) / 2
    
    # Calculate coefficients based on τ
    if τ >= 1e2
        cstτ = 0.5
        coef = sign == :p ? 0.5 : -0.5
    else
        cstτ = (exp(τ) + 1) / (2 * √(exp(2τ) + 1))
        coef = sign == :p ? (exp(τ) - 1) / (2 * √(exp(2τ) + 1)) : (1 - exp(τ)) / (2 * √(exp(2τ) + 1))
    end
    
    os = OpSum()
    
    # Local measurement terms based on neighboring configurations
    if 2 <= i <= N-1
        # Identity term
        os += cstτ, "I", i
        
        # Terms depending on neighboring sites
        # This is a simplified version - you may need to adjust based on exact Fibonacci rules
        os += coef, "n", i-1, "X", i, "n", i+1
        os += coef * (1 - 2 * ϕ^(-1)), "I", i-1, "X", i, "I", i+1
        
    elseif pbc
        if i == 1
            os += cstτ, "I", i
            os += coef, "n", N, "X", i, "n", 2
        elseif i == N
            os += cstτ, "I", i
            os += coef, "n", N-1, "X", i, "n", 1
        end
    end
    
    return MPO(os, sites)
end

"""
    apply_measurement_mps(ψ::MPS, sites, i::Int, τ::Float64, sign::Symbol; pbc::Bool=true) -> MPS, Float64

Apply measurement operator to MPS state and return the resulting state and probability.
"""
function apply_measurement_mps(ψ::MPS, sites, i::Int, τ::Float64, sign::Symbol; pbc::Bool=true)
    # Create measurement operator
    M = measurement_operator_mps(sites, i, τ, sign; pbc=pbc)
    
    # Apply measurement operator
    ψ_measured = apply(M, ψ; cutoff=1e-12)
    
    # Calculate probability (norm squared)
    prob = real(inner(ψ_measured, ψ_measured))
    
    # Normalize the state
    if prob > 1e-12
        ψ_normalized = (1/√prob) * ψ_measured
    else
        ψ_normalized = ψ_measured
    end
    
    return ψ_normalized, prob
end

"""
    mps_sampling(ψ::MPS, sites, measurement_sites::Vector{Int}, τ::Float64; 
                num_samples::Int=1000, pbc::Bool=true) -> Vector{Vector{Symbol}}, Vector{Float64}

Perform sampling measurements on MPS state.
"""
function mps_sampling(ψ::MPS, sites, measurement_sites::Vector{Int}, τ::Float64; 
                     num_samples::Int=1000, pbc::Bool=true)
    num_sites = length(measurement_sites)
    samples = Vector{Vector{Symbol}}(undef, num_samples)
    sample_weights = Vector{Float64}(undef, num_samples)
    
    for sample_idx in 1:num_samples
        current_sequence = Vector{Symbol}(undef, num_sites)
        current_state = copy(ψ)
        total_weight = 1.0
        
        # Measure from left to right
        for (site_idx, measurement_site) in enumerate(measurement_sites)
            # Apply both measurement outcomes
            ψ_p, prob_p = apply_measurement_mps(current_state, sites, measurement_site, τ, :p, pbc)
            prob_m = 1 - prob_p
            
            # Sample based on probabilities
            random_number = rand()
            if random_number < prob_p
                current_sequence[site_idx] = :p
                current_state = ψ_p
                total_weight *= prob_p
            else
                ψ_m, _ = apply_measurement_mps(current_state, sites, measurement_site, τ, :m, pbc)
                current_sequence[site_idx] = :m
                current_state = ψ_m
                total_weight *= prob_m
            end
        end
        
        samples[sample_idx] = current_sequence
        sample_weights[sample_idx] = total_weight
    end
    
    return samples, sample_weights
end

"""
    mps_measurement_enumeration(ψ::MPS, sites, measurement_sites::Vector{Int}, τ::Float64; 
                               pbc::Bool=true) -> Vector{MPS}, Vector{Vector{Symbol}}, Vector{Float64}

Enumerate all possible measurement trajectories on MPS state.
"""
function mps_measurement_enumeration(ψ::MPS, sites, measurement_sites::Vector{Int}, τ::Float64; 
                                   pbc::Bool=true)
    # Initialize with single initial state
    current_level_states = [copy(ψ)]
    current_level_trajectories = [Symbol[]]
    current_level_probabilities = [1.0]
    
    for (measurement_idx, site) in enumerate(measurement_sites)
        next_level_states = Vector{MPS}()
        next_level_trajectories = Vector{Vector{Symbol}}()
        next_level_probabilities = Vector{Float64}()
        
        # Branch for each current state
        for (state_idx, state) in enumerate(current_level_states)
            current_trajectory = current_level_trajectories[state_idx]
            current_prob = current_level_probabilities[state_idx]
            
            # Apply :p measurement
            ψ_p, prob_p = apply_measurement_mps(state, sites, site, τ, :p, pbc)
            if prob_p > 1e-12
                new_trajectory_p = [current_trajectory; :p]
                new_prob_p = current_prob * prob_p
                push!(next_level_states, ψ_p)
                push!(next_level_trajectories, new_trajectory_p)
                push!(next_level_probabilities, new_prob_p)
            end
            
            # Apply :m measurement
            ψ_m, prob_m = apply_measurement_mps(state, sites, site, τ, :m, pbc)
            if prob_m > 1e-12
                new_trajectory_m = [current_trajectory; :m]
                new_prob_m = current_prob * prob_m
                push!(next_level_states, ψ_m)
                push!(next_level_trajectories, new_trajectory_m)
                push!(next_level_probabilities, new_prob_m)
            end
        end
        
        current_level_states = next_level_states
        current_level_trajectories = next_level_trajectories
        current_level_probabilities = next_level_probabilities
    end
    
    return current_level_states, current_level_trajectories, current_level_probabilities
end

"""
    mps_bulk_measurement(ψ::MPS, sites, N::Int, τ::Float64, D::Int; pbc::Bool=true)

Perform bulk measurements on MPS with D layers.
"""
function mps_bulk_measurement(ψ::MPS, sites, N::Int, τ::Float64, D::Int; pbc::Bool=true)
    samples = Vector{Vector{Symbol}}(undef, D)
    sample_weights = Vector{Float64}(undef, D)
    sample_measured_states = Vector{MPS}(undef, D)
    
    current_state = copy(ψ)
    
    for layer in 1:D
        current_sequence = Vector{Symbol}(undef, div(N, 2))
        total_weight = 0.0
        
        # Alternating measurement pattern
        if layer % 2 == 1
            measurement_sites = collect(2:2:N)  # odd layers: even sites
        else
            measurement_sites = collect(1:2:N)  # even layers: odd sites
        end
        
        measurement_τ = (layer == D) ? τ/2 : τ
        
        for (site_idx, measurement_site) in enumerate(measurement_sites)
            ψ_p, prob_p = apply_measurement_mps(current_state, sites, measurement_site, measurement_τ, :p, pbc)
            prob_m = 1 - prob_p
            
            # Sample measurement outcome
            random_number = rand()
            if random_number < prob_p
                current_sequence[site_idx] = :p
                current_state = ψ_p
                total_weight += -log(prob_p)
            else
                ψ_m, _ = apply_measurement_mps(current_state, sites, measurement_site, measurement_τ, :m, pbc)
                current_sequence[site_idx] = :m
                current_state = ψ_m
                total_weight += -log(prob_m)
            end
        end
        
        sample_measured_states[layer] = current_state
        samples[layer] = current_sequence
        sample_weights[layer] = total_weight
    end
    
    return sample_measured_states, samples, sample_weights
end

"""
    calculate_entanglement_entropy_mps(ψ::MPS, b::Int) -> Float64

Calculate entanglement entropy of MPS state with bipartition at bond b.
"""
function calculate_entanglement_entropy_mps(ψ::MPS, b::Int)
    # Perform SVD at bond b
    orthogonalize!(ψ, b)
    U, S, V = svd(ψ[b], (linkind(ψ, b-1), siteind(ψ, b)))
    
    # Calculate entanglement entropy from singular values
    SvN = 0.0
    for n in 1:dim(S, 1)
        p = S[n, n]^2
        if p > 1e-12
            SvN -= p * log(p)
        end
    end
    
    return SvN
end
