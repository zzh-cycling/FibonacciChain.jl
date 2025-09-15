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

julia> N = 8;

julia> ψ_gs, E0 = fibonacci_mps_ground_state(N, pbc=true, maxdim=10, outputlevel=0);

julia> ψ_gs isa MPS
true

julia> E0 < 0 && imag(E0) ≈ 0
true
```
"""
function initial_mps(N::Int)
    # Create sites for Fibonacci anyons
    sites = siteinds("Qubit", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]
    
    # Create MPS from product state
    ψ0 = random_mps(sites, state)
    
    return ψ0, sites
end

function fibonacci_mps_ground_state(N::Int; pbc::Bool=true, sweep_times=5, maxdim=5, cutoff=1e-10, outputlevel=2)
    # Create sites for Fibonacci anyons (using S=1/2 fermions to approximate)
    sites = siteinds("Qubit", N)
    
    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]
    
    # Create MPS from product state
    ψ0 = random_mps(sites, state)
    
    # Create Fibonacci Hamiltonian
    H = fibonacci_hamiltonian_mps(sites; pbc=pbc)
    
    # Find ground state using DMRG
    sweeps = Sweeps(sweep_times)
    setmaxdim!(sweeps, maxdim)
    setcutoff!(sweeps, cutoff)
    
    energy, ψ = dmrg(H, ψ0, sweeps, outputlevel = outputlevel)
    
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

function _apply_measurement_layer!(N::Int64, τ::Float64, sites, ψ::MPS, layer_sample::Vector{Int64}, layer_idx::Int64, pbc::Bool=true; cutoff::Float64=1e-10, maxdim::Int=100, anyon_type::Symbol=:Fibo, return_free_energy::Bool=false) 
    # Helper function to apply measurements to a layer
    measurement_sites, measure_type = _obtain_measurement_config(N, layer_idx, anyon_type)
    F_layer = 0.0

    if return_free_energy
        for (idx, sign) in enumerate(layer_sample)
            ψ, prob = apply_measurement_mps(ψ, sites, measurement_sites[idx], τ, sign; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=measure_type)
            F_layer += -log(prob)
        end
        return ψ, F_layer
    else
        for (idx, sign) in enumerate(layer_sample)
            ψ, prob = apply_measurement_mps(ψ, sites, measurement_sites[idx], τ, sign; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=measure_type)
        end
        return ψ
    end
end

function _sample_layer!(N::Int64, τ_eff::Float64, sites, ψ::MPS,
    rng::MersenneTwister = MersenneTwister(), 
    layer_idx::Int64=1, 
    pbc::Bool=true; 
    anyon_type::Symbol=:Fibo, cutoff::Float64=1e-10, maxdim::Int=100,
    return_free_energy::Bool=true)

    measurement_sites, measure_type = _obtain_measurement_config(N, layer_idx, anyon_type)  
    n = length(measurement_sites)
    sample_layer = zeros(Int, n)
    F_layer = 0.0

    final_state = copy(ψ)
    for (i, site) in enumerate(measurement_sites)
        # first 0 branch
        ψ0, p0 = apply_measurement_mps(ψ, sites, site, τ_eff, 0; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=measure_type)
        p1 = 1 - p0

        if rand(rng) < p0
            sample_layer[i] = 0
            final_state = ψ0
            F_layer += -log(p0)
        else
            # else 1 branch
            ψ1, p1 = apply_measurement_mps(ψ, sites, site, τ_eff, 1; pbc=pbc, cutoff=cutoff, maxdim=maxdim, anyon_type=measure_type)
            sample_layer[i] = 1
            final_state = ψ1
            F_layer += -log(p1)
        end
    end
    return return_free_energy ? (final_state, sample_layer, F_layer) : (final_state, sample_layer)
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

function measure_evolution!(N::Int,
                  τ::Float64,
                  sites,
                  state::MPS,
                  t₂::Int = 1;
                  rng = Random.default_rng(),
                  pbc::Bool = true,
                  anyon_type::Symbol = :Fibo,
                  mode::Symbol = :prob,
                  temp::Bool = false,
                  t₁::Int = 1,
                  return_free_energy::Bool = true, cutoff::Float64=1e-10, maxdim::Int=100,
                  sample::Union{Nothing,Matrix{Int}}=nothing)

    n_measure = anyon_type == :Fibo ? N÷2 : N

    # ---------- Sample decided according to mode ----------
    Δt = t₂ - t₁ + 1 # number of layers to evolve
    Δt > 0 || error("t₂ must be >= t₁")
    sample_free_energy = zeros(Δt) # free energy of each layer
    states = temp ? Vector{MPS}(undef, Δt) : nothing # states of each layer
    current_state = copy(state)

    if mode == :Born
         # 1. Initialize sample matrix
        sample = zeros(Int, Δt, n_measure)   # to be filled during sampling

        # 2. Born trajectory (only for :Born)
        for layer in t₁:t₂
            τ_eff = (layer == t₂) ? τ / 2 : τ

            # Random sampling for this layer
            if return_free_energy
                current_state, sample_layer, F_layer = _sample_layer!(N, τ_eff, sites, current_state, rng, layer, pbc, anyon_type = anyon_type, return_free_energy = return_free_energy, cutoff=cutoff, maxdim=maxdim)
                sample_free_energy[layer-t₁+1] = F_layer
            else
                current_state, sample_layer = _sample_layer!(N, τ_eff, sites, current_state, rng, layer, pbc, anyon_type = anyon_type, return_free_energy = return_free_energy, cutoff=cutoff, maxdim=maxdim)
            end

            sample[layer-t₁+1, :] .= sample_layer
            temp && (states[layer-t₁+1] = current_state)

        end
        
    elseif mode == :sample
        isnothing(sample) && error("When mode=:sample sample must be ::Matrix{Int}")
        size(sample) == (Δt, n_measure) ||
            error("sample size should be ($Δt, $n_measure)")

        for layer in t₁:t₂
            τ_eff = (layer == t₂) ? τ/2 : τ

            if return_free_energy
                current_state, ΔF = _apply_measurement_layer!(
                                N, τ_eff, sites, current_state, 
                                sample[layer, :], layer, pbc; cutoff=cutoff, maxdim=maxdim,
                                anyon_type=anyon_type, return_free_energy=true)
                sample_free_energy[layer] += ΔF
            else
                current_state = _apply_measurement_layer!(
                            N, τ_eff, sites, current_state, 
                            sample[layer, :], layer, pbc; cutoff=cutoff, maxdim=maxdim,
                            anyon_type=anyon_type, return_free_energy=false)
            end

            temp && (states[layer-t₁+1] = current_state)
        end

    else
        error("mode must be ∈ [:Born, :sample]")
    end
    final_state = temp ? states : current_state
    
    return return_free_energy ? (final_state, sample, sample_free_energy) : (final_state, sample)
    
end

"""
    mps_boundary_measure(τ::Float64, ψ::MPS, sites, layer_idx::Int=1; 
                num_samples::Int=1000, pbc::Bool=true, rng::MersenneTwister=MersenneTwister(),
                cutoff::Float64=1e-10, maxdim::Int=100, anyon_type::Symbol=:Fibo) 
                

Perform boundary measurements on MPS state.
# Arguments
- `τ::Float64`: Evolution time parameter
- `ψ::MPS`: Input quantum state
- `sites`: ITensor site indices
- `layer_idx`: the layer index
- `num_samples::Int=1000`: Number of measurement samples
- `pbc::Bool=true`: Periodic boundary conditions
- `rng::MersenneTwister=MersenneTwister()`: Random number generator
- `cutoff::Float64=1e-10`: MPS truncation cutoff
- `maxdim::Int=100`: Maximum bond dimension
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `Vector{Vector{Int64}}`: List of measurement outcomes for each sample
- `Vector{Float64}`: Free energy associated with each sample
"""
function mps_boundary_measure(τ::Float64, ψ::MPS, sites, layer_idx::Int=1; num_samples::Int=1000, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true, cutoff::Float64=1e-10, maxdim::Int=100, anyon_type::Symbol=:Fibo)
    measurement_sites, _ = _obtain_measurement_config(length(sites), layer_idx, anyon_type)
    samples = zeros(Int, num_samples, length(measurement_sites))
    sample_free_energy = Vector{Float64}(undef, num_samples)
    
    for sample_idx in 1:num_samples
        final_state, sample, free_energy  =  _sample_layer!(length(sites), τ, sites, ψ, rng, layer_idx, pbc; anyon_type=anyon_type, return_free_energy=true, cutoff=cutoff, maxdim=maxdim)

        samples[sample_idx, :] = sample
        sample_free_energy[sample_idx] = free_energy
    end
    
    return samples, sample_free_energy
end

"""
    mps_bulk_measure(ψ::MPS, sites, N::Int, τ::Float64, D::Int; pbc::Bool=true)

Perform bulk measurements on MPS with D layers.
"""
function mps_bulk_measure(N::Int, τ::Float64, ψ::MPS, sites, D::Int64; rng::MersenneTwister=MersenneTwister(),pbc::Bool=true, cutoff::Float64=1e-10, maxdim::Int=100, anyon_type::Symbol=:Fibo, return_free_energy::Bool=true)
    sample = zeros(Int, D, div(N,2))
    sample_free_energy = Vector{Float64}(undef, D)
    sample_measured_states = Vector{MPS}(undef, D)
    
    if return_free_energy
        final_state, sample, sample_free_energy = measure_evolution!(N, τ, sites, ψ, D; rng=rng, pbc=pbc, anyon_type=anyon_type, mode=:sample, temp=true, t₁=1, return_free_energy=true, sample=sample, cutoff=cutoff, maxdim=maxdim)
    else
        final_state, sample = measure_evolution!(N, τ, sites, ψ, D; rng=rng, pbc=pbc, anyon_type=anyon_type, mode=:sample, temp=true, t₁=1,return_free_energy=false, sample=sample, cutoff=cutoff, maxdim=maxdim)
    end
    
    
    return return_free_energy ? (final_state, sample, sample_free_energy) : (final_state, sample)
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


"""
    anyon_eelis_mps(N::Int64, ψ::MPS)

Calculate entanglement entropy profile along the chain for MPS state.

# Arguments
- `N::Int64`: System size
- `ψ::MPS`: MPS state

# Returns
- `Vector{Float64}`: Entanglement entropy at each bipartition from left to right

# Examples
```jldoctest
julia> using FibonacciChain, ITensorMPS, ITensors

julia> N = 6;

julia> ψ_gs, E0 = fibonacci_mps_ground_state(N, pbc=true, maxdim=10, outputlevel=0);

julia> # Calculate entanglement entropy profile
       ee_profile = anyon_eelis_mps(N, ψ_gs);

julia> length(ee_profile) == N - 1  # Profile has N-1 points
true

julia> all(x -> x ≥ 0, ee_profile)  # All entropies are non-negative
true
```
"""
function anyon_eelis_mps(N::Int64, ψ::MPS)
    splitlis=Vector(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        EE_lis[m]=ee_mps(ψ, splitlis[m])
    end
    return EE_lis
end