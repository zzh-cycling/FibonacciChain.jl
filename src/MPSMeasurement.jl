"""
MPS-based implementation for Fibonacci chain measurements using ITensor.

This module provides Matrix Product State (MPS) implementations for efficient 
simulation of large Fibonacci anyon chains with measurement protocols.
"""

"""
    fibonacci_mps_ground_state(N::Int; pbc::Bool=true, sweep_times=20, maxdim=50, cutoff=1e-10)

Find ground state of Fibonacci chain Hamiltonian using DMRG.

# Arguments
- `N::Int`: System size
- `pbc::Bool=true`: Periodic boundary conditions
- `sweep_times=20`: Number of DMRG sweeps
- `maxdim=50`: Maximum bond dimension
- `cutoff=1e-10`: Truncation cutoff

# Returns
- `MPS`: Ground state as Matrix Product State
- `Float64`: Ground state energy

# Examples
```jldoctest
julia> using FibonacciChain, ITensorMPS, ITensors

julia> # Find ground state for small system
       N = 8;

julia> ψ_gs, E0 = fibonacci_mps_ground_state(N, pbc=true, maxdim=10);

julia> # Check that we got an MPS
       ψ_gs isa MPS
true

julia> # Energy should be real and negative
       E0 < 0 && imag(E0) ≈ 0
true
```
"""
function initial_mps(N::Int)
    # Create sites for Fibonacci anyons
    sites = siteinds("Qubit", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]
    
    # Create MPS from product state
    ψ0 = randomMPS(sites, state)
    
    return ψ0, sites
end

function fibonacci_mps_ground_state(N::Int; pbc::Bool=true, sweep_times=20, maxdim=50, cutoff=1e-10)
    # Create sites for Fibonacci anyons (using S=1/2 fermions to approximate)
    sites = siteinds("Qubit", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]
    
    # Create MPS from product state
    ψ0 = randomMPS(sites, state)
    
    # Create Fibonacci Hamiltonian
    H = fibonacci_hamiltonian_mps(sites; pbc=pbc)
    
    # Find ground state using DMRG
    sweeps = Sweeps(sweep_times)
    setmaxdim!(sweeps, maxdim)
    setcutoff!(sweeps, cutoff)
    
    energy, ψ = dmrg(H, ψ0, sweeps)
    
    return ψ, energy
end

"""
    fibonacci_hamiltonian_mps(sites; pbc::Bool=true)

Construct Fibonacci chain Hamiltonian as Matrix Product Operator (MPO).

# Arguments
- `sites`: ITensor site indices
- `pbc::Bool=true`: Periodic boundary conditions

# Returns
- `MPO`: Hamiltonian with three-body interactions based on Fibonacci fusion rules
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
    measurement_operator_mps(sites, i::Int, τ::Float64, sign::Int64; pbc::Bool=true, anyon_type::Symbol=:Fibo)

Create local measurement operator at site i as Matrix Product Operator.

# Arguments
- `sites`: ITensor site indices
- `i::Int`: Measurement site
- `τ::Float64`: Evolution time parameter
- `sign::Int64`: Measurement outcome (0 or 1)
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `ITensor`: Local measurement operator incorporating neighboring site correlations
"""
function measurement_operator_mps(sites, i::Int, τ::Float64, sign::Int64; pbc::Bool=true, anyon_type::Symbol=:Fibo)
    N = length(sites)
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (0, 1) "sign must be either 0 or 1"
    
    if anyon_type == :Fibo
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
    
        
        s_im1_idx = i == 1 && pbc ? N : i - 1 #i=1 and pbc, return N, otherwise i-1
        s_i_idx = i
        s_ip1_idx = i == N && pbc ? 1 : i + 1 #i=N and pbc, return 1, otherwise i+1
    
        s_im1 = sites[s_im1_idx]
        s_i = sites[s_i_idx]
        s_ip1 = sites[s_ip1_idx]
    
        id_i = op("I", s_i)
        id_im1 = op("I", s_im1)
        id_ip1 = op("I", s_ip1)
    
        M_local = cstτ * (id_im1 * id_i * id_ip1)
    
        P0_im1 = op("Proj0", s_im1)
        P1_im1 = op("Proj1", s_im1)
        P0_ip1 = op("Proj0", s_ip1)
        P1_ip1 = op("Proj1", s_ip1)
        Z_i = op("Z", s_i)
        X_i = op("X", s_i)
    
        # Local measurement terms based on neighboring configurations
        M_local += coef * (P0_im1 * Z_i * P1_ip1)
        M_local += coef * (P1_im1 * Z_i * P0_ip1)
        M_local += -coef * (P1_im1 * Z_i * P1_ip1)
        M_local += coef * (1 - 2 * ϕ^(-1)) * (P0_im1 * Z_i * P0_ip1)
        M_local += coef * (-2 * ϕ^(-3/2)) * (P0_im1 * X_i * P0_ip1)
    
        
        return M_local

    end
end

"""
    apply_measurement_mps(ψ::MPS, sites, i::Int, τ::Float64, sign::Int64; pbc::Bool=true, cutoff::Float64=1e-10, maxdim::Int=100, anyon_type::Symbol=:Fibo)

Apply measurement operator to MPS state and return post-measurement state.

# Arguments
- `ψ::MPS`: Input quantum state
- `sites`: ITensor site indices
- `i::Int`: Measurement site
- `τ::Float64`: Evolution time parameter
- `sign::Int64`: Measurement outcome
- `pbc::Bool=true`: Periodic boundary conditions
- `cutoff::Float64=1e-10`: MPS truncation cutoff
- `maxdim::Int=100`: Maximum bond dimension
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `MPS`: Post-measurement quantum state
- `Float64`: Measurement probability
"""
function apply_measurement_mps(ψ::MPS, sites, i::Int, τ::Float64, sign::Int64; pbc::Bool=true, cutoff::Float64=1e-10, maxdim::Int=100, anyon_type::Symbol=:Fibo)
    # Create measurement operator
    M = measurement_operator_mps(sites, i, τ, sign; pbc=pbc, anyon_type=anyon_type)
    
    # Apply measurement operator, initial state \psi should be normalized
    ψ_measured = apply(M, ψ; cutoff=cutoff, maxdim=maxdim)
    
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
pbc::Bool=true, cutoff::Float64=1e-10, maxdim::Int=100, anyon_type::Symbol=:Fibo)
    # Initialize with single initial state
    current_level_trajectories = [Int64[]]
    current_level_probabilities = [1.0]
    current_level_states = [copy(ψ)]

    for (measurement_idx, site) in enumerate(measurement_sites)
        next_level_states = Vector{MPS}()
        next_level_trajectories = Vector{Vector{Int64}}()
        next_level_probabilities = Vector{Float64}()
        
        # Branch for each current state
        for (state_idx, state) in enumerate(current_level_states)
            current_trajectory = current_level_trajectories[state_idx]
            current_prob = current_level_probabilities[state_idx]
            
            # Apply 0 measurement
            ψ_p, prob_p = apply_measurement_mps(state, sites, site, τ, 0; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=anyon_type)
            # if prob_p > 1e-12
                new_trajectory_p = [current_trajectory; 0]
                new_prob_p = current_prob * prob_p
                push!(next_level_states, ψ_p)
                push!(next_level_trajectories, new_trajectory_p)
                push!(next_level_probabilities, new_prob_p)
            # end
            
            # Apply 1 measurement
            ψ_m, prob_m = apply_measurement_mps(state, sites, site, τ, 1; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=anyon_type)
            # if prob_m > 1e-12
                new_trajectory_m = [current_trajectory; 1]
                new_prob_m = current_prob * prob_m
                push!(next_level_states, ψ_m)
                push!(next_level_trajectories, new_trajectory_m)
                push!(next_level_probabilities, new_prob_m)
            # end
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
function mps_boundary_measure(ψ::MPS, sites, measurement_sites::Vector{Int}, τ::Float64; num_samples::Int=1000, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true, cutoff::Float64=1e-10, maxdim::Int=100, anyon_type::Symbol=:Fibo)
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
            ψ_p, prob_p = apply_measurement_mps(current_state, sites, measurement_site, τ, 0; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=anyon_type)
            prob_m = 1 - prob_p
            
            # Sample based on probabilities
            random_number = rand(rng)
            if random_number < prob_p
                current_sequence[site_idx] = 0
                current_state = ψ_p
                total_free_energy += -log(prob_p)
            else
                ψ_m, _ = apply_measurement_mps(current_state, sites, measurement_site, τ, 1; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=anyon_type)
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
function mps_bulk_measurement(ψ::MPS, sites, N::Int, τ::Float64, D::Int64; rng::MersenneTwister=MersenneTwister(),pbc::Bool=true, cutoff::Float64=1e-10, maxdim::Int=100, anyon_type::Symbol=:Fibo)
    sample = zeros(Int, D, div(N,2))
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
            ψ_p, prob_p = apply_measurement_mps(current_state, sites, measurement_site, measurement_τ, 0; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=anyon_type)
            prob_m = 1 - prob_p
            
            # Sample measurement outcome
            random_number = rand(rng)
            if random_number < prob_p
                current_sequence[site_idx] = 0
                current_state = ψ_p
                total_free_energy += -log(prob_p)
            else
                ψ_m, _ = apply_measurement_mps(current_state, sites, measurement_site, measurement_τ, 1; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=anyon_type)
                current_sequence[site_idx] = 1
                current_state = ψ_m
                total_free_energy += -log(prob_m)
            end
        end
        
        sample_measured_states[layer] = current_state
        sample[layer, :] = current_sequence
        sample_free_energy[layer] = total_free_energy
    end
    
    return sample_measured_states, sample, sample_free_energy
end

# Helper function to apply measurements to a layer
function apply_measurement_layer_mps!(N::Int64, sites, ψ::MPS, τ::Float64, layer_sample::Vector{Int64}, layer_idx::Int64; pbc::Bool=true, cutoff::Float64=1e-10, maxdim::Int=100, anyon_type::Symbol=:Fibo) 
    if layer_idx % 2 == 1
        measurement_sites = collect(2:2:N)  # odd sites anyons, even sites qubits
    else
        measurement_sites = collect(1:2:N)  # even sites anyons, odd sites qubits
    end
    for (idx, measurement_type) in enumerate(layer_sample)
        ψ, prob = apply_measurement_mps(ψ, sites, measurement_sites[idx], τ, measurement_type; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=anyon_type)
    end
    return ψ
end

function generate_state_mps(τ::Float64, sites, state::MPS, sample::ET, temp::Bool=false; pbc::Bool=true, cutoff::Float64=1e-12, maxdim::Int=1000, anyon_type::Symbol=:Fibo) where{ET}

    if ET == Vector{Int}
        N = 2 * length(sample)
        return apply_measurement_layer_mps!(N, sites, state, τ, sample, 1; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=anyon_type)
        
    elseif ET == Matrix{Int}
        D, N = size(sample, 1), 2 * size(sample, 2)
        statelis = temp ? Vector{MPS}(undef, D) : nothing
        # if ET is Vector{Int64} and temp is true, we return temporary states.
        for layer in 1:D
            τ_eff = (layer == D) ? τ/2 : τ
            state = apply_measurement_layer_mps!(N, sites, state, τ_eff, sample[layer, :], layer; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=anyon_type)

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
    
    if abs(SvN) <1e-14
        SvN = 0.0
    end
    return SvN
end

function anyon_eelis_mps(N::Int64, ψ::MPS)
    splitlis=Vector(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        EE_lis[m]=ee_mps(ψ, splitlis[m])
    end
    return EE_lis
end