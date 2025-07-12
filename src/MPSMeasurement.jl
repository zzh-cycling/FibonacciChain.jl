"""
MPS-based implementation for Fibonacci chain measurements using ITensor
"""

"""
    fibonacci_mps_ground_state(N::Int; pbc::Bool=true) -> MPS

Generate the ground state of Fibonacci chain as an MPS.
"""
function initial_mps(N::Int)
    # Create sites for Fibonacci anyons
    sites = siteinds("Qubit", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]
    
    # Create MPS from product state
    ψ0 = randomMPS(sites, state)
    
    return ψ0
end

function fibonacci_mps_ground_state(N::Int; pbc::Bool=true)
    # Create sites for Fibonacci anyons (using S=1/2 fermions to approximate)
    sites = siteinds("Qubit", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]
    
    # Create MPS from product state
    ψ0 = randomMPS(sites, state)
    
    # Create Fibonacci Hamiltonian
    H = fibonacci_hamiltonian_mps(sites; pbc=pbc)
    
    # Find ground state using DMRG
    sweeps = Sweeps(200)
    setmaxdim!(sweeps, 5000)
    setcutoff!(sweeps, 1E-10)
    
    energy, ψ = dmrg(H, ψ0, sweeps)
    
    return ψ, energy
end

"""
    fibonacci_hamiltonian_mps(sites; pbc::Bool=true) -> MPO

Create Fibonacci chain Hamiltonian as an MPO using ITensor.
"""
function fibonacci_hamiltonian_mps(sites; pbc::Bool=true)
    N = length(sites)
    os = OpSum()
    
    # Golden ratio
    ϕ = (1 + √5) / 2
    coef = 1/2
    # Three-body interactions for Fibonacci chain
    for i in 2:(N-1)
        # Add three-body terms based on Fibonacci fusion rules
        os += coef, "Proj0", i-1, "Z", i, "Proj1", i+1
        os += coef, "Proj1", i-1, "Z", i, "Proj0", i+1
        os += -coef, "Proj1", i-1, "Z", i, "Proj1", i+1
        os += coef * (1 - 2 * ϕ^(-1)), "Proj0", i-1, "Z", i, "Proj0", i+1
        os += coef * (-2 * ϕ^(-3/2)), "Proj0", i-1, "X", i, "Proj0", i+1
    end
    
    # Periodic boundary conditions
    if pbc && N > 2
        # H1 term
        os += coef, "Proj0", N, "Z", 1, "Proj1", 2
        os += coef, "Proj1", N, "Z", 1, "Proj0", 2
        os += -coef, "Proj1", N, "Z", 1, "Proj1", 2
        os += coef * (1 - 2 * ϕ^(-1)), "Proj0", N, "Z", 1, "Proj0", 2
        os += coef * (-2 * ϕ^(-3/2)), "Proj0", N, "X", 1, "Proj0", 2
        # HN term (wrap around)
        os += coef, "Proj0", N-1, "Z", N, "Proj1", 1
        os += coef, "Proj1", N-1, "Z", N, "Proj0", 1
        os += -coef, "Proj1", N-1, "Z", N, "Proj1", 1
        os += coef * (1 - 2 * ϕ^(-1)), "Proj0", N-1, "Z", N, "Proj0", 1
        os += coef * (-2 * ϕ^(-3/2)), "Proj0", N-1, "X", N, "Proj0", 1
    end
    
    return MPO(os, sites)
end

"""
    measurement_operator_mps(sites, i::Int, τ::Float64, sign::Int64; pbc::Bool=true) -> MPO

Create measurement operator at site i with parameter τ as an MPO.
"""
function measurement_operator_mps(sites, i::Int, τ::Float64, sign::Int64; pbc::Bool=true)
    N = length(sites)
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (0, 1) "sign must be either 0 or 1"
    
    # Golden ratio
    ϕ = (1 + √5) / 2
    
    # Calculate coefficients based on τ
    if τ >= 1e2
        cstτ = 0.5
        coef = sign == 0 ? 0.5 : -0.5
    else
        cstτ = (exp(τ) + 1) / (2 * √(exp(2τ) + 1))
        coef = sign == 0 ? (exp(τ) - 1) / (2 * √(exp(2τ) + 1)) : (1 - exp(τ)) / (2 * √(exp(2τ) + 1))
    end
    
    os = OpSum()
    
    # Local measurement terms based on neighboring configurations
    if 2 <= i <= N-1
        # Identity term
        os += cstτ, "I", i
        # PσP terms
        os += coef, "Proj0", i-1, "Z", i, "Proj1", i+1
        os += coef, "Proj1", i-1, "Z", i, "Proj0", i+1
        os += -coef, "Proj1", i-1, "Z", i, "Proj1", i+1
        os += coef * (1 - 2 * ϕ^(-1)), "Proj0", i-1, "Z", i, "Proj0", i+1
        os += coef * (-2 * ϕ^(-3/2)), "Proj0", i-1, "X", i, "Proj0", i+1
    elseif pbc
        if i == 1
            os += cstτ, "I", i
            os += coef, "Proj0", N, "Z", i, "Proj1", i+1
            os += coef, "Proj1", N, "Z", i, "Proj0", i+1
            os += -coef, "Proj1", N, "Z", i, "Proj1", i+1
            os += coef * (1 - 2 * ϕ^(-1)), "Proj0", N, "Z", i, "Proj0", i+1
            os += coef * (-2 * ϕ^(-3/2)), "Proj0", N, "X", i, "Proj0", i+1
        elseif i == N
            os += cstτ, "I", i
            os += coef, "Proj0", i-1, "Z", i, "Proj1", 1
            os += coef, "Proj1", i-1, "Z", i, "Proj0", 1
            os += -coef, "Proj1", i-1, "Z", i, "Proj1", 1
            os += coef * (1 - 2 * ϕ^(-1)), "Proj0", i-1, "Z", i, "Proj0", 1
            os += coef * (-2 * ϕ^(-3/2)), "Proj0", i-1, "X", i, "Proj0", 1
        end
    end
    
    return MPO(os, sites)
end

"""
    apply_measurement_mps(ψ::MPS, sites, i::Int, τ::Float64, sign::Int64; pbc::Bool=true) -> MPS, Float64

Apply measurement operator to MPS state and return the resulting state and probability.
"""
function apply_measurement_mps(ψ::MPS, sites, i::Int, τ::Float64, sign::Int64; pbc::Bool=true)
    # Create measurement operator
    M = measurement_operator_mps(sites, i, τ, sign; pbc=pbc)
    
    # Apply measurement operator
    ψ_measured = apply(M, ψ; cutoff=1e-12)
    
    # Calculate probability (norm squared)
    prob = real(inner(ψ_measured, ψ_measured))
    
    # Normalize the state
    ψ_normalized = normalize(ψ_measured)
    
    return ψ_normalized, prob
end

"""
    mps_measurement_enumeration(ψ::MPS, sites, measurement_sites::Vector{Int}, τ::Float64; 
                               pbc::Bool=true) -> Vector{MPS}, Vector{Vector{Int64}}, Vector{Float64}

Enumerate all possible measurement trajectories on MPS state.
"""
function mps_measurement_enumeration(ψ::MPS, sites, measurement_sites::Vector{Int}, τ::Float64; 
pbc::Bool=true)
    # Initialize with single initial state
    current_level_trajectories = [Int64[]]
    current_level_probabilities = [1.0]
    
    for (measurement_idx, site) in enumerate(measurement_sites)
        next_level_states = Vector{MPS}()
        next_level_trajectories = Vector{Vector{Int64}}()
        next_level_probabilities = Vector{Float64}()
        
        # Branch for each current state
        for (state_idx, state) in enumerate(current_level_states)
            current_trajectory = current_level_trajectories[state_idx]
            current_prob = current_level_probabilities[state_idx]
            
            # Apply 0 measurement
            ψ_p, prob_p = apply_measurement_mps(state, sites, site, τ, 0, pbc)
            if prob_p > 1e-12
                new_trajectory_p = [current_trajectory; 0]
                new_prob_p = current_prob * prob_p
                push!(next_level_states, ψ_p)
                push!(next_level_trajectories, new_trajectory_p)
                push!(next_level_probabilities, new_prob_p)
            end
            
            # Apply 1 measurement
            ψ_m, prob_m = apply_measurement_mps(state, sites, site, τ, 1, pbc)
            if prob_m > 1e-12
                new_trajectory_m = [current_trajectory; 1]
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
    mps_boundary_measure(ψ::MPS, sites, measurement_sites::Vector{Int}, τ::Float64; 
                num_samples::Int=1000, pbc::Bool=true) -> Vector{Vector{Int64}}, Vector{Float64}

Perform boundary measurements on MPS state.
"""
function mps_boundary_measure(ψ::MPS, sites, measurement_sites::Vector{Int}, τ::Float64; num_samples::Int=1000, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true)
    num_sites = length(measurement_sites)
    samples = Vector{Vector{Int64}}(undef, num_samples)
    sample_free_energy = Vector{Float64}(undef, num_samples)
    
    for sample_idx in 1:num_samples
        current_sequence = Vector{Int64}(undef, num_sites)
        current_state = copy(ψ)
        total_free_energy = 0.0
        
        # Measure from left to right
        for (site_idx, measurement_site) in enumerate(measurement_sites)
            # Apply both measurement outcomes
            ψ_p, prob_p = apply_measurement_mps(current_state, sites, measurement_site, τ, 0, pbc)
            prob_m = 1 - prob_p
            
            # Sample based on probabilities
            random_number = rand(rng)
            if random_number < prob_p
                current_sequence[site_idx] = 0
                current_state = ψ_p
                total_free_energy += -log(prob_p)
            else
                ψ_m, _ = apply_measurement_mps(current_state, sites, measurement_site, τ, 1, pbc)
                current_sequence[site_idx] = 1
                current_state = ψ_m
                total_free_energy += -log(prob_m)
            end
        end
        
        samples[sample_idx] = current_sequence
        sample_free_energy[sample_idx] = total_free_energy
    end
    
    return samples, sample_free_energy
end

"""
    mps_bulk_measurement(ψ::MPS, sites, N::Int, τ::Float64, D::Int; pbc::Bool=true)

Perform bulk measurements on MPS with D layers.
"""
function mps_bulk_measurement(ψ::MPS, sites, N::Int, τ::Float64, D::Int; rng::MersenneTwister=MersenneTwister(),pbc::Bool=true)
    samples = Vector{Vector{Int64}}(undef, D)
    sample_free_energy = Vector{Float64}(undef, D)
    sample_measured_states = Vector{MPS}(undef, D)
    
    current_state = copy(ψ)
    
    for layer in 1:D
        current_sequence = Vector{Int64}(undef, div(N, 2))
        total_free_energy = 0.0
        
        # Alternating measurement pattern
        if layer % 2 == 1
            measurement_sites = collect(2:2:N)  # odd layers: even sites
        else
            measurement_sites = collect(1:2:N)  # even layers: odd sites
        end
        
        measurement_τ = (layer == D) ? τ/2 : τ

        for (site_idx, measurement_site) in enumerate(measurement_sites)
            ψ_p, prob_p = apply_measurement_mps(current_state, sites, measurement_site, measurement_τ, 0, pbc)
            prob_m = 1 - prob_p
            
            # Sample measurement outcome
            random_number = rand(rng)
            if random_number < prob_p
                current_sequence[site_idx] = 0
                current_state = ψ_p
                total_free_energy += -log(prob_p)
            else
                ψ_m, _ = apply_measurement_mps(current_state, sites, measurement_site, measurement_τ, 1, pbc)
                current_sequence[site_idx] = 1
                current_state = ψ_m
                total_free_energy += -log(prob_m)
            end
        end
        
        sample_measured_states[layer] = current_state
        samples[layer] = current_sequence
        sample_free_energy[layer] = total_free_energy
    end
    
    return sample_measured_states, samples, sample_free_energy
end

# Helper function to apply measurements to a layer
function apply_measurement_layer_mps!(state::MPS, N::Int64, τ::Float64, layer_sample::Vector{Int64}, layer_idx::Int64, pbc::Bool=true) 
    if layer_idx % 2 == 1
        measurement_sites = collect(2:2:N)  # odd sites anyons, even sites qubits
    else
        measurement_sites = collect(1:2:N)  # even sites anyons, odd sites qubits
    end
    for (idx, measurement_type) in enumerate(layer_sample)
        state = apply_measurement_mps(state, measurement_sites[idx], τ, measurement_type, pbc)
        normalize!(state)
    end
    return state
end

function Generate_state_mps(τ::Float64, state::MPS, sample::ET, temp::Bool=false, pbc::Bool=true) where{ET}

    if ET == Vector{Int}
        N = 2 * length(sample)
        return apply_measurement_layer!(state, N, τ, sample, 1, pbc)
        
    elseif ET == Matrix{Int}
        D, N = size(sample, 1), 2 * size(sample, 2)
        statelis = temp ? Vector{Vector{T}}(undef, D) : nothing
        # if ET is Vector{Int64} and temp is true, we return temporary states.
        for layer in 1:D
            τ_eff = (layer == D) ? τ/2 : τ
            state = apply_measurement_layer_mps!(state, N, τ_eff, sample[layer, :], layer, pbc)
            
            if temp
                statelis[layer] = copy(state)
            end
        end
        
        return temp ? statelis : state
    end
end

"""
    calculate_entanglement_entropy_mps(ψ::MPS, b::Int) -> Float64

Calculate entanglement entropy of MPS state with bipartition at bond b.
"""
function ee_mps(ψ::MPS, b::Int)
    # Perform SVD at bond b
    ψ = orthogonalize(ψ, b)
    # U, S, V = svd(ψ[b], linkind(ψ, b-1), siteind(ψ, b))
    U, S, V = svd(ψ[b], (linkinds(ψ, b-1)..., siteinds(ψ, b)...))
    # 
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

function eelis_Fibo_mps(N::Int64, ψ::MPS)
    splitlis=Vector(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        EE_lis[m]=ee_mps(ψ, splitlis[m])
    end
    return EE_lis
end