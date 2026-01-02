"""
MPS-based implementation for anyon chain measurements using ITensor.

This module provides Matrix Product State (MPS) implementations for efficient
simulation of large anyon chains with measurement protocols.
"""

"""
    anyon_mps_gst(model::AnyonModel; sweep_times=5, maxdim=5, cutoff=1e-10, outputlevel=1)

Find ground state of anyon chain Hamiltonian using DMRG.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters (N, pbc, anyon_type)
- `sweep_times=5`: Number of DMRG sweeps
- `maxdim=5`: Maximum bond dimension
- `cutoff=1e-10`: Truncation cutoff
- `outputlevel=1`: Verbosity level

# Returns
- `MPS`: Ground state as Matrix Product State
- `Float64`: Ground state energy

# Examples
```jldoctest
julia> using FibonacciChain, ITensorMPS, ITensors

julia> model = AnyonModel(FibonacciAnyon(), 8; pbc=true);

julia> ψ_gs, E0 = anyon_mps_gst(model, maxdim=10, outputlevel=0);

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

function anyon_mps_gst(model::AnyonModel{AT}; sweep_times=5, maxdim=5, cutoff=1e-10, outputlevel=1, kwargs...) where AT <: AbstractAnyonType
    # Create sites for anyons (using S=1/2 fermions to approximate)
    N = model.N
    sites = siteinds("Qubit", N)

    # Create initial product state (vacuum state)
    state = ["0" for _ in 1:N]
    
    # Create MPS from product state
    ψ0 = random_mps(sites, state)
    
    # Create anyon Hamiltonian
    H = anyon_ham(model, sites, kwargs...)

    # Find ground state using DMRG
    sweeps = Sweeps(sweep_times)
    setmaxdim!(sweeps, maxdim)
    setcutoff!(sweeps, cutoff)
    
    energy, ψ = dmrg(H, ψ0, sweeps, outputlevel = outputlevel)
    
    return ψ, energy
end

"""
    anyon_ham(model::AnyonModel, sites; kwargs...)

Construct anyon chain Hamiltonian as Matrix Product Operator (MPO).

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `sites`: ITensor site indices
- `kwargs...`: Additional parameters (e.g., `J`, `h` for Ising model)

# Returns
- `MPO`: Hamiltonian as Matrix Product Operator
  - Fibonacci: three-body interactions based on Fibonacci fusion rules
  - Ising: transverse field Ising model with ZZ and X terms
"""
function anyon_ham(model::AnyonModel{AT}, sites::Vector{<:Index}; kwargs...) where AT <: FibonacciAnyon
    N = length(sites)
    os = OpSum()
    pbc = model.pbc
    
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

function anyon_ham(model::AnyonModel{AT}, sites::Vector{<:Index}; kwargs...) where AT <: IsingAnyon
    J, h = get(kwargs, :J, 1.0), get(kwargs, :h, 1.0)
    N = length(sites)
    os = OpSum()
    pbc = model.pbc

    for i in 1:N
        os += -h, "X", i
    end

    for i in 1:N-1
        os += -J, "Z", i, "Z", i+1
    end
    if pbc && N > 2
        os += -J, "Z", N, "Z", 1
    end
    
    return MPO(os, sites)
end

"""
    measurement_operator_mps(model::AnyonModel, sites, i::Int, τ::Float64, sign::Bool)

Create local measurement operator at site i as ITensor.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `sites`: ITensor site indices
- `i::Int`: Measurement site
- `τ::Float64`: Measurement strength parameter
- `sign::Bool`: Measurement outcome (false for +, true for -)

# Returns
- `ITensor`: Local measurement operator incorporating neighboring site correlations
"""
function measurement_operator_mps(model::AnyonModel{AT}, sites, i::Int, τ::Float64, sign::Bool;) where AT <: FibonacciAnyon
    @assert model.measure_operator ∈ [:Ferro, :Antiferro] "measure_operator must be :Ferro or :Antiferro"
    N = length(sites)
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (true, false) "sign must be either true or false"
    pbc = model.pbc

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

    
    s_im1_idx = (i == 1 && pbc) ? N : i - 1 #i=1 and pbc, return N, otherwise i-1
    s_i_idx = i
    s_ip1_idx = (i == N && pbc) ? 1 : i + 1 #i=N and pbc, return 1, otherwise i+1

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

function measurement_operator_mps(model::AnyonModel{AT}, sites, i::Int, τ::Float64, sign::Bool;) where AT <: IsingAnyon
    @assert model.measure_operator in [:X, :ZZ] "measure_operator must be either :X or :ZZ"
    pbc = model.pbc
    N = length(sites)
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (true, false) "sign must be either true or false"

    if model.measure_operator == :X
        if τ >= 1e2
            cstτ = 0.5
            coef = sign == 0 ? 0.5 : -0.5
        else
            cstτ = cosh(τ/2) / √(2cosh(τ))
            coef = sign == 0 ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
        end
        
        M_local = cstτ * op("I", sites[i]) + coef * op("X", sites[i])

    elseif model.measure_operator == :ZZ
        if τ >= 1e2
            cstτ = 0.5
            coef = sign == 0 ? 0.5 : -0.5
        else
            cstτ = (exp(τ) + 1) / (2 * √(exp(2τ) + 1))
            coef = sign == 0 ? (exp(τ) - 1) / (2 * √(exp(2τ) + 1)) : (1 - exp(τ)) / (2 * √(exp(2τ) + 1))
        end
        
        idx_p1 = (i == N && pbc) ? 1 : i + 1 #i=N and pbc, return 1, otherwise i+1
        Z_i = op("Z", sites[i])
        Z_ip1 = op("Z", sites[idx_p1])

        M_local = cstτ * (op("I", sites[i]) * op("I", sites[idx_p1])) + coef * (Z_i * Z_ip1)

    end
end

#   elseif (anyon_type ∈ (:reset, :resetFibo) && τ >= 1e2)|| anyon_type == :IsingZ
#         if τ >= 1e2
#             cstτ = 0.5
#             coef = sign == 0 ? 0.5 : -0.5
#         else
#             cstτ = cosh(τ/2) / √(2cosh(τ))
#             coef = sign == 0 ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
#         end
        
#         M_local = cstτ * op("I", sites[i]) + coef * op("Z", sites[i])

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
function measuremap(model::AnyonModel{AT}, ψ::MPS, sites, i::Int, τ::Float64, sign::Bool; cutoff::Float64=1e-10, maxdim::Int=100) where AT <: AbstractAnyonType
    # Create measurement operator
    M = measurement_operator_mps(model, sites, i, τ, sign)

    # Apply measurement operator, initial state \psi should be normalized
    ψ_measured = apply(M, ψ; cutoff=cutoff, maxdim=maxdim)
    
    # Calculate probability (norm squared)
    prob = real(inner(ψ_measured, ψ_measured))
    
    # Normalize the state
    ψ_normalized = normalize(ψ_measured)
    
    return ψ_normalized, prob
end

"""
    boundary_evolution(model::AnyonModel, state::Vector{T}, measure_config::MeasureConfig, 
                       sample::Union{Nothing, BitVector}=nothing; layer_idx::Int=1)
    boundary_evolution(model::AnyonModel, sites::Vector{<:Index}, state::MPS, measure_config::MeasureConfig, sample::Union{Nothing, BitVector} =nothing; cutoff::Float64=1e-12, maxdim::Int=1000, layer_idx::Int=1)

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
function boundary_evolution(anyon_model::AnyonModel{AT}, sites::Vector{<:Index}, state::MPS, measure_config::MeasureConfig, sample::Union{Nothing, BitVector} =nothing; cutoff::Float64=1e-12, maxdim::Int=1000, layer_idx::Int=1) where AT <: AbstractAnyonType
    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    if measure_config.mode == :sample
        N = anyon_model.N
        size(sample, 1) == measurement_num(anyon_model.anyon_type)*(N ÷ 2) || error("sample size mismatch with anyon_model $(N)")
        return _apply_measurement_layer_mps(anyon_model, measure_config.τ, sites, state, sample, layer_idx; cutoff=cutoff, maxdim=maxdim)
    elseif measure_config.mode == :Born
        return _sample_layer_mps(anyon_model, measure_config.τ, sites, state, measure_config.rng, layer_idx, verbose=measure_config.verbose, cutoff=cutoff, maxdim=maxdim)
    end
end

function _apply_measurement_layer_mps(model::AnyonModel{AT}, τ::Float64, sites::Vector{<:Index}, ψ::MPS, layer_sample::BitVector, layer_idx::Int64; cutoff::Float64=1e-10, maxdim::Int=100) where AT <: AbstractAnyonType
    # Helper function to apply measurements to a layer
    measurement_sites, measure_anyon_model = _obtain_measurement_config(model, layer_idx)
    F_layer = 0.0

    for (idx, sign) in enumerate(layer_sample)
        ψ, prob = measuremap(measure_anyon_model, ψ, sites, measurement_sites[idx], τ, sign; cutoff=cutoff, maxdim=maxdim)
        F_layer += -log(prob)
    end
    return Measurement_outcome_mps_boundary(ψ, layer_sample, F_layer)
end

function _sample_layer_mps(model::AnyonModel{AT}, τ_eff::Float64, sites, ψ::MPS,
    rng::MersenneTwister = MersenneTwister(), 
    layer_idx::Int64=1;
    cutoff::Float64=1e-10, maxdim::Int=100,
    verbose::Bool=false) where AT <: AbstractAnyonType

    measurement_sites, measure_anyon_model = _obtain_measurement_config(model, layer_idx)  
    n = length(measurement_sites)
    sample_layer = BitVector(zeros(Bool, n))
    F_layer = 0.0

    final_state = copy(ψ)
    for (i, site) in enumerate(measurement_sites)
        # first 0 branch
        ψ0, p0 = measuremap(measure_anyon_model, ψ, sites, site, τ_eff, false; cutoff=cutoff, maxdim=maxdim)
        p1 = 1 - p0

        randomNumber = rand(rng)
        verbose && @show randomNumber
        if randomNumber < p0
            sample_layer[i] = 0
            ψ = ψ0
            F_layer += -log(p0)
            verbose && @show -log(p0)
        else
            # else 1 branch
            ψ1, p1_re = measuremap(measure_anyon_model, ψ, sites, site, τ_eff, true; cutoff=cutoff, maxdim=maxdim)
            sample_layer[i] = 1
            ψ = ψ1
            F_layer += -log(p1)
            verbose && @show -log(p1)
        end
    end
    return Measurement_outcome_mps_boundary(ψ, sample_layer, F_layer)
end

"""
    mps_measurement_enumeration(model::AnyonModel, ψ::MPS, sites, measurement_sites::Vector{Int}, 
                                τ::Float64; cutoff::Float64=1e-10, maxdim::Int=100)

Enumerate all possible measurement trajectories on MPS state.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `ψ::MPS`: Input quantum state
- `sites`: ITensor site indices
- `measurement_sites::Vector{Int}`: Sites to measure
- `τ::Float64`: Measurement strength parameter
- `cutoff::Float64=1e-10`: MPS truncation cutoff
- `maxdim::Int=100`: Maximum bond dimension

# Returns
- `Vector{MPS}`: Final states for each trajectory
- `Vector{Vector{Bool}}`: Measurement trajectories
- `Vector{Float64}`: Probabilities for each trajectory
"""
function mps_measurement_enumeration(model::AnyonModel{AT}, ψ::MPS, sites, measurement_sites::Vector{Int}, τ::Float64; cutoff::Float64=1e-10, maxdim::Int=100) where AT <: AbstractAnyonType
    # Initialize with single initial state
    current_level_trajectories = [Bool[]]
    current_level_probabilities = [1.0]
    current_level_states = [copy(ψ)]

    for (measurement_idx, site) in enumerate(measurement_sites)
        next_level_states = Vector{MPS}()
        next_level_trajectories = Vector{Vector{Bool}}()
        next_level_probabilities = Vector{Float64}()
        
        # Branch for each current state
        for (state_idx, state) in enumerate(current_level_states)
            current_trajectory = current_level_trajectories[state_idx]
            current_prob = current_level_probabilities[state_idx]
            
            # Apply 0 measurement
            ψ_p, prob_p = measuremap(model, state, sites, site, τ, false; cutoff=cutoff, maxdim=maxdim)
            # if prob_p > 1e-12
                new_trajectory_p = [current_trajectory; false]
                new_prob_p = current_prob * prob_p
                push!(next_level_states, ψ_p)
                push!(next_level_trajectories, new_trajectory_p)
                push!(next_level_probabilities, new_prob_p)
            # end
            
            # Apply 1 measurement
            ψ_m, prob_m = measuremap(model, state, sites, site, τ, true; cutoff=cutoff, maxdim=maxdim)
            # if prob_m > 1e-12
                new_trajectory_m = [current_trajectory; true]
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
    add_reference_qubits!(model::AnyonModel, ψ::MPS, sites, site_idx::Int=1;
                          k_new::Int=1, verbose::Bool=false, entangle_way::Symbol=:copy)

Add reference qubit(s) to MPS at specified site index.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `ψ::MPS`: Input quantum state
- `sites`: ITensor site indices
- `site_idx::Int=1`: Site index to entangle with reference qubit
- `k_new::Int=1`: Number of reference qubits to add (0 or 1)
- `verbose::Bool=false`: Verbosity flag
- `entangle_way::Symbol=:copy`: Method to entangle reference qubit (`:copy` or `:reset`)

# Returns
For `:copy` mode:
- `MPS`: New MPS with reference qubit added
- `Int`: Index of added reference qubit

For `:reset` mode:
- `MPS`: New MPS with reference qubit added
- `Int`: Index of added reference qubit  
- `Float64`: Measurement probability
"""
function add_reference_qubits!(model::AnyonModel{AT}, ψ::MPS, sites, site_idx::Int=1;
                               k_new::Int=1, 
                               verbose::Bool=false,
                               entangle_way::Symbol=:copy) where AT <: AbstractAnyonType
    1 ≤ site_idx ≤ length(ψ) || error("site_idx out of range")

    N = length(sites)         
    d = dim(sites[site_idx])          # local physical dimension

    to_add_block = copy(ψ[site_idx])
    to_add_block = replaceind!(to_add_block, siteind(ψ, site_idx), Index(d, "Qubit, Site, n=$(N+1)"))

    ψ_new = typeof(ψ)(undef, length(ψ)+1)
    ψ_new[1] = to_add_block
    ψ_new[2:end] = ψ[1:end]

    # ---------- 2. 根据 entangle_way 做 Bell 对 ----------
    if entangle_way == :copy
        # 把原 site_idx 的态拷贝到参考比特，然后做 |0⟩+|1⟩  Bell 对
        # 步骤：先逐位做 CNOT（Ref→site_idx），再把 Ref 转到 |+⟩
        for n in 1:(site_idx)            # 把参考比特“搬”到 target 旁边
            n < site_idx+1 && move_site!(ψ_new, n, n+1)   # 右移 1 格
        end
        r = 1                            # 参考比特现在与 site_idx 相邻
        t = r + 1                        # target 物理比特

        # |+⟩_Ref
        apply!(ψ_new, r, op("H", sRef))
        # CNOT_{Ref,target}
        apply!(ψ_new, [r, t], op("CNOT", sRef, s[t]))

    elseif entangle_way == :reset
        # 先对原比特做 Z 测量，再根据结果把参考比特设为 |0⟩ 或 |1⟩
        # 这里用 ITensors 的“测量+坍缩”接口
        result = sample!(ψ_new, site_idx+1, "Z")  # +1 因为插了 Ref
        verbose && @info "measurement result = $result"

        # 把参考比特设成 |result⟩
        apply!(ψ_new, 1, op(result==1 ? "X" : "I", sRef))

        # 再做 H 和 CNOT 形成 Bell 对
        apply!(ψ_new, 1, op("H", sRef))
        apply!(ψ_new, [1, 2], op("CNOT", sRef, s[2]))

        # 返回测量概率（简化：假设本征值 ±1 各占 50%，可再精确算）
        prob = 0.5
        verbose && @show prob
        return ψ_new, 1, prob   # 多返回一个概率
    else
        error("unknown entangle_way = $entangle_way")
    end

    # ---------- 3. 可选：把参考比特移回最左端 ----------
    for n in site_idx:-1:2
        move_site!(ψ_new, n, n-1)
    end

    verbose && @show maxlinkdim(ψ_new)
    return ψ_new, 1
end

"""
    move_site!(ψ::MPS, i::Int, j::Int)

Move MPS tensor from position i to position j by sequential SWAP operations.
Maintains canonical form throughout. Only for adjacent moves.
"""
function move_site!(ψ::MPS, i::Int, j::Int)
    step = j > i ? 1 : -1
    cur = i
    while cur != j
        nxt = cur + step
        # 做 SWAP 门
        apply!(ψ, [min(cur,nxt), max(cur,nxt)], op("SWAP", site_type(ψ,cur), site_type(ψ,nxt)))
        cur = nxt
    end
end

"""
    bulk_evolution(model::AnyonModel, state::Vector{ET}, measure_config::MeasureConfig,
                   samples::Union{Nothing,BitMatrix}=nothing)
    bulk_evolution(model::AnyonModel, sites::Vector{<:Index}, state::MPS, measure_config::MeasureConfig,
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
"""
function bulk_evolution(model::AnyonModel{AT},
                  sites::Vector{<:Index},
                  state::MPS,
                  measure_config::MeasureConfig,
                  samples::Union{Nothing,BitMatrix}=nothing; cutoff::Float64=1e-10, maxdim::Int=100) where AT <: AbstractAnyonType

    # ---------- Sample decided according to mode ----------
    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    current_state = copy(state)
    if mode == :Born
        return _born_measure_mps(model, sites, current_state, measure_config; cutoff=cutoff, maxdim=maxdim)
    elseif mode == :sample
        return _sample_measure_mps(model, sites, current_state, samples, measure_config; cutoff=cutoff, maxdim=maxdim)
    end
end

"""
    _born_measure_mps(model, sites, current_state, measure_config; ...)

Evolve an MPS state using probabilistic Born-rule sampling for a specified number of time steps.

This internal helper function is called by `bulk_evolution` when `mode` is `:Born`.

# Arguments
- `model::AnyonModel`: The anyon model.
- `sites`: ITensor site indices.
- `current_state::MPS`: The initial MPS state.
- `measure_config::MeasureConfig`: Configuration containing `τ`, `t₁`, `t₂`, `rng`, etc.
- `cutoff::Float64=1e-10`: MPS truncation cutoff.
- `maxdim::Int=100`: Maximum bond dimension.

# Returns
- `Measurement_outcome_mps_bulk`: A struct containing:
  - `states::Vector{MPS}`: Intermediate states at each full time step.
  - `samples::BitMatrix`: The generated measurement outcome sequences.
  - `free_energy::Vector{Float64}`: The free energy for each measurement layer.
"""
function _born_measure_mps(model::AnyonModel{AT}, sites::Vector{<:Index}, current_state::MPS, measure_config::MeasureConfig; cutoff::Float64=1e-10, maxdim::Int=100) where AT <: AbstractAnyonType

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

    # 1. Initialize sample matrix
    samples = BitMatrix(undef, (D, n_measure))   # to be filled during sampling
    sample_free_energy = zeros(D)
    states = Vector{MPS}(undef, Δt)

    for period in 1:Δt
        # Random sampling for this period
        τ_eff = (period == Δt && enable_τ_eff) ? τ/2 : τ

        outcome1 = _sample_layer_mps(
            model, τ, sites, current_state, rng, 2*period-1; cutoff=cutoff, maxdim=maxdim, verbose=verbose)
        current_state = outcome1.state
        samples[2*period-1, :] = outcome1.sample
        sample_free_energy[2*period-1] = outcome1.free_energy

        outcome2 = _sample_layer_mps(
            model, τ_eff, sites, current_state, rng, 2*period; cutoff=cutoff, maxdim=maxdim, verbose=verbose)
        current_state = outcome2.state
        samples[2*period, :] = outcome2.sample
        sample_free_energy[2*period] = outcome2.free_energy

        states[period] = current_state
    end

    return Measurement_outcome_mps_bulk(states, samples, sample_free_energy)
end

"""
    _sample_measure_mps(model, sites, current_state, samples, measure_config; ...)

Evolve an MPS state using a predefined measurement trajectory.

This internal helper function is called by `bulk_evolution` when `mode` is `:sample`.

# Arguments
- `model::AnyonModel`: The anyon model.
- `sites`: ITensor site indices.
- `current_state::MPS`: The initial MPS state.
- `samples::BitMatrix`: The predefined measurement outcomes.
- `measure_config::MeasureConfig`: Configuration containing `τ`, `t₁`, `t₂`, etc.
- `cutoff::Float64=1e-10`: MPS truncation cutoff.
- `maxdim::Int=100`: Maximum bond dimension.

# Returns
- `Measurement_outcome_mps_bulk`: A struct containing:
  - `states::Vector{MPS}`: Intermediate states at each full time step.
  - `samples::BitMatrix`: The input measurement outcome sequences.
  - `free_energy::Vector{Float64}`: The free energy for each measurement layer.
"""
function _sample_measure_mps(model::AnyonModel{AT}, sites::Vector{<:Index}, current_state::MPS, samples::BitMatrix, measure_config::MeasureConfig; cutoff::Float64=1e-10, maxdim::Int=100) where AT <: AbstractAnyonType

    n_measure = measurement_num(model.anyon_type)*(model.N÷2)
    τ = measure_config.τ
    t₁ = measure_config.t₁
    t₂ = measure_config.t₂
    enable_τ_eff = measure_config.enable_τ_eff

    Δt = t₂ - t₁ + 1
    Δt >= 0 || error("t₂ must be >= t₁")
    D = Δt * 2 # number of layers to evolve

    sample_free_energy = zeros(D)
    states = Vector{MPS}(undef, Δt)
    
    # Validate sample matrix
    isnothing(samples) && error("When mode=:sample samples must be ::BitMatrix")
    size(samples) == (D, n_measure) || error("sample size should be ($D, $n_measure)")

    for period in 1:Δt
        τ_eff = (period == Δt && enable_τ_eff) ? τ/2 : τ

        outcome1 = _apply_measurement_layer_mps(
            model, τ, sites, current_state, samples[2*period-1, :], 2*period-1; cutoff=cutoff, maxdim=maxdim)
        current_state = outcome1.state
        sample_free_energy[2*period-1] = outcome1.free_energy

        outcome2 = _apply_measurement_layer_mps(
            model, τ_eff, sites, current_state, samples[2*period, :], 2*period; cutoff=cutoff, maxdim=maxdim)
        current_state = outcome2.state
        sample_free_energy[2*period] = outcome2.free_energy

        states[period] = current_state
    end

    return Measurement_outcome_mps_bulk(states, samples, sample_free_energy)
end

struct Measurement_outcome_mps_bulk
    states::Vector{MPS}
    samples::BitMatrix
    free_energys::Vector{Float64}
end

struct Measurement_outcome_mps_boundary
    state::MPS
    sample::BitVector
    free_energy::Float64
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
function reference_evolution(model::AnyonModel{AT}, sites, forward::Vector{MPS}, sample::Matrix{Bool}, 
    x₂::Int, t₁, t₂, measure_config::MeasureConfig; x₁::Int=1, verbose::Bool=false) where AT <: AbstractAnyonType
    
    N = model.N
    τ = measure_config.τ
    rng = measure_config.rng
    mode = measure_config.mode
    
    n_measure = measurement_num(model.anyon_type)*(N÷2)
    Δt = size(sample, 1) ÷ 2
    D = size(sample, 1)   # D is the number of layers, while Δt is the true time(# period)

    @assert size(sample, 2) == n_measure "Sample size spatial dimension must be $n_measure, but got $(size(sample, 2))"
    @assert 1 <= t₁ <= t₂ <= D÷2 "Time slice t₁ must before time slice t₂, both must be in the range [1, $(D÷2)]"
    @assert 1 <= x₁ <= x₂ <= N "Site index x₁ must be smaller than site index x₂, both must be in the range [1, $(N)]"
    @assert mode ∈ [:sample, :Born] "mode must be either :sample or :Born, but got $mode"

    δt = t₂ - t₁ 
    δx = abs(x₂ - x₁) 
    state = forward[t₁]
    statelis = Vector{MPS}(undef, Δt) 
    view(statelis, 1:t₁) .= view(forward, 1:t₁)
    sample_layer = BitMatrix(undef, (size(sample, 1), n_measure))
    view(sample_layer, 1:t₁, :) .= view(sample, 1:t₁, :)
    sample_free_energy = zeros(Float64, D)


    if δt > 0 && δx > 0 # 3 ref qubits, both spatial and temporal correlation, actually 3-point correlation.
        verbose && @info "t₁ = $(t₁), t₂ = $(t₂), x₁ = $(x₁), x₂ = $(x₂), 3 refs"

        state1 = add_reference_qubits!(model, state, sites, x₁; verbose=verbose)
        state2 = add_reference_qubits!(model, state1, sites, x₂; verbose=verbose)
    
        config1 = MeasureConfig(τ=τ, t₂=(t₂-t₁), rng=rng, mode=mode, t₁=1, verbose=verbose, enable_τ_eff=false)
        final_stlis1, samples1, free_energy1 = bulk_evolution(model, sites, state2, sample[2*t₁+1:2*t₂, :], config1)

        state3 = add_reference_qubits!(model, final_stlis1[end], sites, x₂; verbose=verbose)

        config2 = MeasureConfig(τ=τ, t₂=(Δt-t₂), rng=rng, mode=mode, t₁=1, verbose=verbose, enable_τ_eff=true)
        final_stlis2, samples2, free_energy2 = bulk_evolution(model, sites, state3, sample[2*t₂+1:end, :], config2)

        view(statelis, t₁+1:t₂) .= view(final_stlis1, :)
        view(statelis, t₂+1:Δt) .= view(final_stlis2, :)

        sample_layer[2*t₁+1:2*t₂, :] .= samples1
        sample_layer[2*t₂+1:end, :] .= samples2
        sample_free_energy[2*t₁+1:2*t₂] .= free_energy1
        sample_free_energy[2*t₂+1:end] .= free_energy2

    elseif δt == 0 # 2 ref qubits, pure 2-point spatial correlation
        verbose && @info "x₁ = $(x₁), x₂ = $(x₂), δx = $(δx), at time slice t₁ = t₂ = $(t₁), 2 refs"
    
        state1 = add_reference_qubits!(model, state, sites, x₁; verbose=verbose)
        state2 = add_reference_qubits!(model, state1, sites, x₂; verbose=verbose)
        
        config2 = MeasureConfig(τ=τ, t₂=(Δt-t₁), rng=rng, mode=mode, t₁=1, verbose=verbose, enable_τ_eff=true)
        final_stlis2, samples2, free_energy2 = bulk_evolution(model, sites, state2, sample[2*t₁+1:end, :], config2)

        view(statelis, t₁+1:Δt) .= view(final_stlis2, :)
        sample_layer[2*t₁+1:end, :] .= samples2
        sample_free_energy[2*t₁+1:end] .= free_energy2

    elseif δx == 0 # 2 ref qubits, pure 2-point temporal correlation
        verbose && @info "t₁ = $(t₁), t₂ = $(t₂), δt = $(δt), at site x₁ = x₂ = $(x₂), 2 refs"

        state1 = add_reference_qubits!(model, state, sites, x₂; verbose=verbose)

        config1 = MeasureConfig(τ=τ, t₂=(t₂-t₁), rng=rng, mode=mode, t₁=1, verbose=verbose, enable_τ_eff=false)
        final_stlis1, samples1, free_energy1 = bulk_evolution(model, sites, state1, sample[2*t₁+1:2*t₂, :], config1)

        state2 = add_reference_qubits!(model, final_stlis1[end], sites, x₂; verbose=verbose)
    
        config2 = MeasureConfig(τ=τ, t₂=(Δt-t₂), rng=rng, mode=mode, t₁=1, verbose=verbose, enable_τ_eff=true)
        final_stlis2, samples2, free_energy2 = bulk_evolution(model, sites, state2, sample[2*t₂+1:end, :], config2)

        view(statelis, t₁+1:t₂) .= view(final_stlis1, :)
        view(statelis, t₂+1:Δt) .= view(final_stlis2, :)
        sample_layer[2*t₁+1:2*t₂, :] .= samples1
        sample_layer[2*t₂+1:end, :] .= sample[2*t₂+1:end, :]
        sample_free_energy[2*t₁+1:2*t₂] .= free_energy1
        sample_free_energy[2*t₂+1:end] .= free_energy2
    end

    return statelis, sample_layer, sample_free_energy
end

"""
    ee_mps(ψ::MPS, b::Int)

Calculate entanglement entropy of MPS state with bipartition at bond b.

# Arguments
- `ψ::MPS`: MPS state
- `b::Int`: Bond index for bipartition (1 ≤ b < N)

# Returns
- `Float64`: Entanglement entropy (von Neumann entropy) at bond b
"""
function ee_mps(ψ::MPS, b::Int)
    @assert 1 ≤ b < length(ψ) "Bond index b must be in the range [1, N-1]"
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
    anyon_eelis(model::AnyonModel, ψ::MPS)

Calculate entanglement entropy profile along the chain for MPS state.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `ψ::MPS`: MPS state

# Returns
- `Vector{Float64}`: Entanglement entropy at each bipartition from left to right

# Examples
```jldoctest
julia> using FibonacciChain, ITensorMPS, ITensors

julia> model = AnyonModel(FibonacciAnyon(), 6; pbc=true);

julia> ψ_gs, E0 = anyon_mps_gst(model, maxdim=10, outputlevel=0);

julia> ee_profile = anyon_eelis(model, ψ_gs);

julia> length(ee_profile) == model.N - 1  # Profile has N-1 points
true

julia> all(x -> x ≥ 0, ee_profile)  # All entropies are non-negative
true
```
"""
function anyon_eelis(model::AnyonModel, ψ::MPS)
    N = model.N
    splitlis=Vector(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        EE_lis[m]=ee_mps(ψ, splitlis[m])
    end
    return EE_lis
end