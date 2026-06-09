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
    state = ["0" for _ = 1:N]

    # Create MPS from product state
    ψ0 = productMPS(sites, state)

    return ψ0, sites
end

"""
    evenparity_mps(N::Int)

Create an MPS for the even-parity state in the Z basis, i.e. the equal-weight
superposition of all bitstrings with an even number of |1⟩.

# Returns
- `MPS`: normalized even-parity MPS
- `Vector{Index}`: site indices
"""
function evenparity_mps(N::Int)
    # Create sites for Fibonacci anyons
    sites = siteinds("Qubit", N)

    # Create initial product state in X eigenbasis |+>
    state = fill("+", N)

    # Create MPS from product state
    ψ0 = productMPS(sites, state)

    return ψ0, sites
end

function anyon_mps_gst(
    model::AnyonModel{AT};
    sweep_times = 5,
    maxdim = 5,
    cutoff = 1e-10,
    outputlevel = 1,
) where {AT<:AbstractAnyonType}
    # Create sites for anyons (using S=1/2 fermions to approximate)
    N = model.N
    sites = siteinds("Qubit", N)

    # Create initial product state (vacuum state)
    state = ["0" for _ = 1:N]

    # Create MPS from product state
    ψ0 = random_mps(sites, state)

    # Create anyon Hamiltonian
    H = anyon_ham(model, sites)

    # Find ground state using DMRG
    sweeps = Sweeps(sweep_times)
    setmaxdim!(sweeps, maxdim)
    setcutoff!(sweeps, cutoff)

    energy, ψ = dmrg(H, ψ0, sweeps, outputlevel = outputlevel)

    return ψ, energy
end

"""
    anyon_ham(model::AnyonModel, sites::Vector{<:Index}; kwargs...)

Construct anyon chain Hamiltonian as Matrix Product Operator (MPO).

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `sites::Vector{<:Index}`: ITensor site indices
- `kwargs...`: Additional parameters (e.g., `J`, `h` for Ising model)

# Returns
- `MPO`: Hamiltonian as Matrix Product Operator
  - Fibonacci: three-body interactions based on Fibonacci fusion rules
  - Ising: transverse field Ising model with ZZ and X terms
"""
function anyon_ham(model::AnyonModel{FibonacciAnyon}, sites::Vector{<:Index}; kwargs...)
    N = length(sites)
    os = OpSum()
    pbc = model.pbc

    measure_operator = model.measure_operator # Default measurement operator for Fibonacci anyons
    # Golden ratio
    ϕ = (1 + √5) / 2

    if measure_operator == :Antiferro # Default measurement 
        coef = 1/2
        # Three-body interactions for Fibonacci chain
        for i = 2:(N-1)
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
    elseif measure_operator == :Ferro # Easily found that Ferro GS basically is the most excited state of Antiferro, so H_{Ferro}= -H\_{Antiferro}
        coef = 1/2
        # Three-body interactions for Fibonacci chain (Ferro version)
        # Ferro differs from Antiferro in sign conventions
        for i = 2:(N-1)
            # Add three-body terms based on Fibonacci fusion rules (Ferro)
            os += -coef, "Proj0", i-1, "Z", i, "Proj1", i+1
            os += -coef, "Proj1", i-1, "Z", i, "Proj0", i+1
            os += coef, "Proj1", i-1, "Z", i, "Proj1", i+1  # +coef instead of -coef
            os += coef * (2 * ϕ^(-1) - 1), "Proj0", i-1, "Z", i, "Proj0", i+1  # different coefficient
            os += coef * (2 * ϕ^(-3/2)), "Proj0", i-1, "X", i, "Proj0", i+1  # positive instead of negative
        end

        # Periodic boundary conditions
        if pbc && N > 2
            # H1 term
            os += -coef, "Proj0", N, "Z", 1, "Proj1", 2
            os += -coef, "Proj1", N, "Z", 1, "Proj0", 2
            os += coef, "Proj1", N, "Z", 1, "Proj1", 2
            os += coef * (2 * ϕ^(-1) - 1), "Proj0", N, "Z", 1, "Proj0", 2
            os += coef * (2 * ϕ^(-3/2)), "Proj0", N, "X", 1, "Proj0", 2
            # HN term (wrap around)
            os += -coef, "Proj0", N-1, "Z", N, "Proj1", 1
            os += -coef, "Proj1", N-1, "Z", N, "Proj0", 1
            os += coef, "Proj1", N-1, "Z", N, "Proj1", 1
            os += coef * (2 * ϕ^(-1) - 1), "Proj0", N-1, "Z", N, "Proj0", 1
            os += coef * (2 * ϕ^(-3/2)), "Proj0", N-1, "X", N, "Proj0", 1
        end
    end

    return MPO(os, sites)

end

function anyon_ham(model::AnyonModel{IsingAnyon}, sites::Vector{<:Index})
    J = get_interaction_param(model, :J, 1.0)
    h = get_interaction_param(model, :h, 1.0)
    N = length(sites)
    os = OpSum()
    pbc = model.pbc

    for i = 1:N
        os += -h, "X", i
    end

    for i = 1:(N-1)
        os += -J, "Z", i, "Z", i+1
    end
    if pbc && N > 2
        os += -J, "Z", N, "Z", 1
    end

    return MPO(os, sites)
end

function anyon_ham(model::AnyonModel{OBFAnyon}, sites::Vector{<:Index})
    λ = get_interaction_param(model, :λ, 1.0)
    λI = get_interaction_param(model, :λI, 1.0)  # Ising coupling strength

    N = length(sites)
    os = OpSum()
    pbc = model.pbc

    for i = 1:N
        os -= λI, "X", i
    end

    for i = 1:(N-1)
        os -= λI, "Z", i, "Z", i+1
    end
    if pbc && N > 2
        os -= λI, "Z", N, "Z", 1
    end

    for i = 1:(N-2)
        os += λ/2, "X", i, "Z", i+1, "Z", i+2
        os += λ/2, "Z", i, "Z", i+1, "X", i+2
    end

    if pbc && N > 2
        # Wrap around terms
        os += λ/2, "X", N-1, "Z", N, "Z", 1
        os += λ/2, "Z", N-1, "Z", N, "X", 1
        os += λ/2, "X", N, "Z", 1, "Z", 2
        os += λ/2, "Z", N, "Z", 1, "X", 2
    end

    return MPO(os, sites)
end

"""
    measurement_operator_mps(model::AnyonModel, sites::Vector{<:Index}, i::Int, τ::Float64, sign::Bool)

Create local measurement operator at site i as ITensor.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `sites::Vector{<:Index}`: ITensor site indices
- `i::Int`: Measurement site
- `τ::Float64`: Measurement strength parameter
- `sign::Bool`: Measurement outcome (false for +, true for -)

# Returns
- `ITensor`: Local measurement operator incorporating neighboring site correlations
"""
function measurement_operator_mps(
    model::AnyonModel{FibonacciAnyon},
    sites::Vector{<:Index},
    i::Int,
    τ::Float64,
    sign::Bool;
)
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
        coef = sign ? -0.5 : 0.5
    else
        cstτ = (exp(τ) + 1) / (2 * √(exp(2τ) + 1))
        coef =
            sign ? (1 - exp(τ)) / (2 * √(exp(2τ) + 1)) : (exp(τ) - 1) / (2 * √(exp(2τ) + 1))
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
    if model.measure_operator == :Antiferro
        M_local += coef * (P0_im1 * Z_i * P1_ip1)
        M_local += coef * (P1_im1 * Z_i * P0_ip1)
        M_local += -coef * (P1_im1 * Z_i * P1_ip1)
        M_local += coef * (1 - 2 * ϕ^(-1)) * (P0_im1 * Z_i * P0_ip1)
        M_local += coef * (-2 * ϕ^(-3/2)) * (P0_im1 * X_i * P0_ip1)
    elseif model.measure_operator == :Ferro
        M_local += -coef * (P0_im1 * Z_i * P1_ip1)
        M_local += -coef * (P1_im1 * Z_i * P0_ip1)
        M_local += coef * (P1_im1 * Z_i * P1_ip1)  # +coef instead of -coef
        M_local += coef * (2 * ϕ^(-1) - 1) * (P0_im1 * Z_i * P0_ip1)  # different coefficient
        M_local += coef * (2 * ϕ^(-3/2)) * (P0_im1 * X_i * P0_ip1)  # positive instead of negative
    end

    return M_local
end

function measurement_operator_mps(
    model::AnyonModel{IsingAnyon},
    sites::Vector{<:Index},
    i::Int,
    τ::Float64,
    sign::Bool;
)
    @assert model.measure_operator in [:X, :ZZ] "measure_operator must be either :X or :ZZ"
    pbc = model.pbc
    N = length(sites)
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (true, false) "sign must be either true or false"

    if τ >= 1e2
        cstτ = 0.5
        coef = sign ? -0.5 : 0.5
    else
        cstτ = cosh(τ/2) / √(2cosh(τ))
        coef = sign ? - sinh(τ/2) / √(2cosh(τ)) : sinh(τ/2) / √(2cosh(τ))
    end

    if model.measure_operator == :X

        M_local = cstτ * op("I", sites[i]) + coef * op("X", sites[i])

    elseif model.measure_operator == :ZZ
        idx_p1 = (i == N && pbc) ? 1 : i + 1 #i=N and pbc, return 1, otherwise i+1
        Z_i = op("Z", sites[i])
        Z_ip1 = op("Z", sites[idx_p1])

        M_local = cstτ * (op("I", sites[i]) * op("I", sites[idx_p1])) + coef * (Z_i * Z_ip1)

    end

    return M_local
end

function measurement_operator_mps(
    model::AnyonModel{OBFAnyon},
    sites::Vector{<:Index},
    i::Int,
    τ::Float64,
    sign::Bool;
)
    @assert model.measure_operator in [:XZZ, :ZZX, :ZZ, :X] "measure_operator must be :XZZ, :ZZX, :ZZ, :X"
    pbc = model.pbc
    N = length(sites)
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (true, false) "sign must be either true or false"

    # Common coefficients for all operators, here in contrast to Ising case the sign convention is inverse.
    if τ >= 1e2
        cstτ = 0.5
        coef = sign ? 0.5 : -0.5
    else
        cstτ = cosh(τ/2) / √(2cosh(τ))
        coef = sign ? sinh(τ/2) / √(2cosh(τ)) : -sinh(τ/2) / √(2cosh(τ))
    end

    # Calculate indices for neighboring sites
    if pbc
        i1, i2 = mod1(i + 1, N), mod1(i + 2, N)
    else
        @assert 1 <= i <= N - 2 "Index i must be in [1, N-2] for OBC (OBF)"
        i1, i2 = i + 1, i + 2
    end

    if model.measure_operator ∈ [:ZZ, :X]
        # Inline Ising implementation to avoid model allocation
        if τ >= 1e2
            cstτ_ising = 0.5
            coef_ising = sign ? -0.5 : 0.5
        else
            cstτ_ising = cosh(τ/2) / √(2cosh(τ))
            coef_ising = sign ? -sinh(τ/2) / √(2cosh(τ)) : sinh(τ/2) / √(2cosh(τ))
        end

        if model.measure_operator == :X
            M_local = cstτ_ising * op("I", sites[i]) + coef_ising * op("X", sites[i])
        else  # :ZZ
            idx_p1 = (i == N && pbc) ? 1 : i + 1
            Z_i = op("Z", sites[i])
            Z_ip1 = op("Z", sites[idx_p1])
            M_local =
                cstτ_ising * (op("I", sites[i]) * op("I", sites[idx_p1])) +
                coef_ising * (Z_i * Z_ip1)
        end
        return M_local
    elseif model.measure_operator == :XZZ
        # XZZ: X on site i, ZZ on sites i+1 and i+2
        X_i = op("X", sites[i])
        I_i = op("I", sites[i])
        Z_i1 = op("Z", sites[i1])
        Z_i2 = op("Z", sites[i2])
        I_i1 = op("I", sites[i1])
        I_i2 = op("I", sites[i2])

        M_local = cstτ * (I_i * I_i1 * I_i2) + coef * (X_i * Z_i1 * Z_i2)
    elseif model.measure_operator == :ZZX
        # ZZX: ZZ on sites i and i+1, X on site i+2
        Z_i = op("Z", sites[i])
        Z_i1 = op("Z", sites[i1])
        X_i2 = op("X", sites[i2])
        I_i = op("I", sites[i])
        I_i1 = op("I", sites[i1])
        I_i2 = op("I", sites[i2])

        M_local = cstτ * (I_i * I_i1 * I_i2) + coef * (Z_i * Z_i1 * X_i2)
    end

    return M_local
end


function measuremap(
    model::AnyonModel{AT},
    ψ::MPS,
    sites::Vector{<:Index},
    i::Int,
    τ::Float64,
    sign::Bool;
    cutoff::Float64 = 1e-10,
    maxdim::Int = 100,
) where {AT<:AbstractAnyonType}
    # Create measurement operator
    M = measurement_operator_mps(model, sites, i, τ, sign)
    return _measuremap_with_operator(ψ, M; cutoff = cutoff, maxdim = maxdim)
end

function _measuremap_with_operator(
    ψ::MPS,
    M::ITensor;
    cutoff::Float64 = 1e-10,
    maxdim::Int = 100,
    truncate_per_event::Bool = true,
)
    # Apply measurement operator, initial state ψ should be normalized
    ψ_measured =
        truncate_per_event ? apply(M, ψ; cutoff = cutoff, maxdim = maxdim) : apply(M, ψ)

    # Calculate probability (norm squared)
    prob = real(inner(ψ_measured, ψ_measured))

    # Normalize the state in-place to avoid extra MPS allocation
    normalize!(ψ_measured)

    return ψ_measured, prob
end

function boundary_evolution(
    anyon_model::AnyonModel{AT},
    sites::Vector{<:Index},
    state::MPS,
    measure_config::MeasureConfig,
    sample::Union{Nothing,BitVector} = nothing;
    layer_idx::Int = 1,
) where {AT<:AbstractAnyonType}
    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    cutoff = measure_config.cutoff
    maxdim = measure_config.maxdim
    truncate_every_events = measure_config.truncate_every_events
    truncate_every_events >= 1 || error("truncate_every_events must be >= 1")
    if measure_config.mode == :sample
        N = anyon_model.N
        size(sample, 1) == measurement_num(anyon_model.anyon_type)*(N ÷ 2) ||
            error("sample size mismatch with anyon_model $(N)")
        return _apply_measurement_layer_mps(
            anyon_model,
            measure_config.τ,
            sites,
            state,
            sample,
            layer_idx;
            cutoff = cutoff,
            maxdim = maxdim,
            truncate_every_events = truncate_every_events,
        )
    elseif measure_config.mode == :Born
        return _stochastic_measurement_layer_mps_mps(
            anyon_model,
            measure_config.τ,
            sites,
            state,
            measure_config.rng,
            layer_idx,
            verbose = measure_config.verbose,
            cutoff = cutoff,
            maxdim = maxdim,
            truncate_every_events = truncate_every_events,
        )
    end
end

function _apply_measurement_layer_mps(
    model::AnyonModel{AT},
    τ::Float64,
    sites::Vector{<:Index},
    ψ::MPS,
    layer_sample::BitVector,
    layer_idx::Int64;
    cutoff::Float64 = 1e-10,
    maxdim::Int = 100,
    truncate_every_events::Int = 1,
    normalized::Bool = true,
) where {AT<:AbstractAnyonType}
    # Helper function to apply measurements to a layer
    measurement_sites, measure_anyon_model, measurement_strength =
        _obtain_measurement_config(model, layer_idx, τ)
    F_layer = 0.0
    n = length(measurement_sites)
    operators = Vector{ITensor}(undef, n)

    # Cache layer operators to avoid rebuilding them in each apply.
    @inbounds for k = 1:n
        operators[k] = measurement_operator_mps(
            measure_anyon_model,
            sites,
            measurement_sites[k],
            measurement_strength,
            layer_sample[k],
        )
    end

    do_per_event_truncate = (truncate_every_events == 1)

    if normalized
        if do_per_event_truncate
            # Preserve legacy behavior for exact RNG trajectory compatibility.
            @inbounds for idx = 1:n
                ψ, prob = measuremap(
                    measure_anyon_model,
                    ψ,
                    sites,
                    measurement_sites[idx],
                    measurement_strength,
                    layer_sample[idx];
                    cutoff = cutoff,
                    maxdim = maxdim,
                )
                F_layer += -log(prob)
            end
        else
            # If truncating less frequently, we apply all operators first and then truncate at the end    
            @inbounds for idx = 1:n
                if idx % truncate_every_events == 0
                    truncate_signal=true
                else
                    truncate_signal=false
                end
                ψ, prob = _measuremap_with_operator(
                    ψ,
                    operators[idx];
                    cutoff = cutoff,
                    maxdim = maxdim,
                    truncate_per_event = truncate_signal,
                )
                F_layer += -log(prob)
            end
        end
    else
        # Unnormalized evolution: apply operators without normalizing
        if do_per_event_truncate
            @inbounds for idx = 1:n
                ψ = apply(operators[idx], ψ; cutoff = cutoff, maxdim = maxdim)
            end
        else
            @inbounds for idx = 1:n
                if idx % truncate_every_events == 0
                    ψ = apply(operators[idx], ψ; cutoff = cutoff, maxdim = maxdim)
                else
                    ψ = apply(operators[idx], ψ)
                end
            end
        end
    end
    return Measurement_outcome_mps_boundary(ψ, layer_sample, Float32(F_layer))
end

function _stochastic_measurement_layer_mps_mps(
    model::AnyonModel{AT},
    τ::Float64,
    sites::Vector{<:Index},
    ψ::MPS,
    rng::MersenneTwister = MersenneTwister(),
    layer_idx::Int64 = 1;
    cutoff::Float64 = 1e-10,
    maxdim::Int = 100,
    verbose::Bool = false,
    truncate_every_events::Int = 1,
) where {AT<:AbstractAnyonType}

    measurement_sites, measure_anyon_model, measurement_strength =
        _obtain_measurement_config(model, layer_idx, τ)
    n = length(measurement_sites)
    sample_layer = BitVector(zeros(Bool, n))
    F_layer = 0.0
    operators_false = Vector{ITensor}(undef, n)
    operators_true = Vector{ITensor}(undef, n)

    # Build local operators once per layer and reuse for branch evaluations.
    @inbounds for i = 1:n
        site = measurement_sites[i]
        operators_false[i] = measurement_operator_mps(
            measure_anyon_model,
            sites,
            site,
            measurement_strength,
            false,
        )
        operators_true[i] = measurement_operator_mps(
            measure_anyon_model,
            sites,
            site,
            measurement_strength,
            true,
        )
    end

    do_per_event_truncate = (truncate_every_events == 1)

    if do_per_event_truncate
        # Preserve legacy behavior for exact RNG trajectory compatibility.
        @inbounds for idx = 1:n
            site = measurement_sites[idx]
            ψ0, p0 = measuremap(
                measure_anyon_model,
                ψ,
                sites,
                site,
                measurement_strength,
                false;
                cutoff = cutoff,
                maxdim = maxdim,
            )
            p1 = 1 - p0

            randomNumber = rand(rng)
            verbose && @show randomNumber
            if randomNumber < p0
                sample_layer[idx] = 0
                ψ = ψ0
                F_layer += -log(p0)
            else
                ψ, _ = measuremap(
                    measure_anyon_model,
                    ψ,
                    sites,
                    site,
                    measurement_strength,
                    true;
                    cutoff = cutoff,
                    maxdim = maxdim,
                )
                sample_layer[idx] = 1
                F_layer += -log(p1)
            end
        end
    else
        @inbounds for i = 1:n
            if i % truncate_every_events == 0
                truncate_signal=true
            else
                truncate_signal=false
            end
            # Compute probability of outcome 0 via measuremap
            ψ0, p0 = _measuremap_with_operator(
                ψ,
                operators_false[i];
                cutoff = cutoff,
                maxdim = maxdim,
                truncate_per_event = truncate_signal,
            )
            p1 = 1 - p0

            randomNumber = rand(rng)
            verbose && @show randomNumber
            if randomNumber < p0
                sample_layer[i] = 0
                ψ = ψ0  # reuse the already-computed branch
                F_layer += -log(p0)
                verbose && @show -log(p0)
            else
                # Discard ψ0 (goes out of scope), compute only the needed branch
                ψ0 = nothing  # release ψ0 memory before allocating ψ1
                ψ1, _ = _measuremap_with_operator(
                    ψ,
                    operators_true[i];
                    cutoff = cutoff,
                    maxdim = maxdim,
                    truncate_per_event = truncate_signal,
                )
                sample_layer[i] = 1
                ψ = ψ1
                F_layer += -log(p1)
                verbose && @show -log(p1)
            end
        end
    end

    return Measurement_outcome_mps_boundary(ψ, sample_layer, Float32(F_layer))
end

"""
    transfer_matrix_subspace_mps(model::AnyonModel, sites::Vector{<:Index}, τ::Float64, sample::BitMatrix; 
                                 n_states::Int=10, cutoff::Float64=1e-10, maxdim::Int=100)

Compute the dominant spectrum of the transfer matrix via subspace iteration, using MPS states.

This is the MPS analogue of [`transfer_matrix_subspace`](@ref). The algorithm initializes
`n_states` product states (basis vectors), then iteratively applies the transfer matrix
for each time slice's measurement outcome. At each step, states are orthogonalized via a
Gram-matrix Cholesky procedure (the MPS equivalent of QR), and the diagonal Cholesky
factors are recorded as the spectrum.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `sites::Vector{<:Index}`: ITensor site indices
- `τ::Float64`: Measurement strength parameter
- `sample::BitMatrix`: Measurement outcome sequences (rows = layers, cols = sites).
  The number of rows must be divisible by `layers_per_period(model.anyon_type)`.
- `n_states::Int=10`: Number of initial basis vectors to propagate
- `cutoff::Float64=1e-10`: MPS truncation cutoff
- `maxdim::Int=100`: Maximum bond dimension for MPS operations

# Returns
- `Matrix{Float64}`: Matrix of size `(k, t)` where `k = min(n_states, length(anyon_basis(model)))`
  and `t = size(sample,1) ÷ layers_per_period`. Each column contains `-log.(abs.(diag(L)))`
  for that time step.

# Examples
```jldoctest
julia> using FibonacciChain, ITensorMPS, ITensors

julia> L = 8; τ = atanh(0.95);

julia> model = AnyonModel(FibonacciAnyon(), L; pbc = true);

julia> sites = siteinds("Qubit", L);

julia> sample = BitMatrix(ones(Int8, 2, div(L, 2)));

julia> spectrum = transfer_matrix_subspace_mps(model, sites, τ, sample; n_states = 5);

julia> size(spectrum, 1) == 5
true
```
"""
function transfer_matrix_subspace_mps(
    model::AnyonModel{AT},
    sites::Vector{<:Index},
    τ::Float64,
    sample::BitMatrix;
    n_states::Int = 10,
    cutoff::Float64 = 1e-10,
    maxdim::Int = 100,
) where {AT<:AbstractAnyonType}
    # Here the transfer matrix is not hermitian, thus the Schur vector is not eigenvectors.
    # We need to do a QR-like projection via Gram-matrix Cholesky. When the non-hermitian
    # matrix is too ill-conditioned, the eigen fallback is used.
    n_layers = layers_per_period(model.anyon_type)
    D_layers, n_cols = size(sample)
    @assert D_layers % n_layers == 0 "Number of layers $D_layers must be divisible by $n_layers"
    t = D_layers ÷ n_layers
    n_cols == _samples_per_layer(model) ||
        error("sample size spatial dimension must be $(_samples_per_layer(model)), got $n_cols")

    basis = anyon_basis(model)
    l = length(basis)
    k = min(n_states, l)
    N = length(sites)

    # Initialize k product states (basis vectors)
    states = Vector{MPS}(undef, k)
    for i in 1:k
        buf = basis[i].buf
        state_str = [bit ? "1" : "0" for bit in reverse(digits(Bool, buf; base=2, pad=N))]
        states[i] = productMPS(sites, state_str)
    end
    spectrum_tlis = zeros(k, t)

    for step in 1:t
        sample_layer = sample[(step - 1) * n_layers + 1 : step * n_layers, :]
        config = MeasureConfig(τ = τ, mode = :sample, t₂ = 1, enable_τ_eff = false)
        for i in 1:k
            outcome = _sample_measure_mps(
                model,
                sites,
                states[i],
                sample_layer,
                config;
                cutoff = cutoff,
                maxdim = maxdim,
                normalized = false,
            )
            states[i] = outcome.state
        end

        # Compute Gram matrix
        G = zeros(Float64, k, k)
        for i in 1:k
            for j in i:k
                val = real(inner(states[i], states[j]))
                G[i, j] = val
                G[j, i] = val
            end
        end

        # Gram-matrix Cholesky (MPS analogue of QR)
        try
            F = cholesky(Hermitian(G))
            L = F.L
            Linv = inv(L)

            # Form orthonormal states: Q_i = sum_j (Linv)_{ij} ψ_j
            new_states = Vector{MPS}(undef, k)
            for i in 1:k
                ψ_new = Linv[i, 1] * states[1]
                for j in 2:k
                    ψ_new = ψ_new + Linv[i, j] * states[j]
                end
                new_states[i] = truncate(ψ_new; cutoff = cutoff, maxdim = maxdim)
            end
            states = new_states

            # Note here do not sort, will distort the spectrum
            spectrum_tlis[:, step] = -log.(abs.(diag(L)))
        catch e
            if e isa PosDefException
                # Fallback to eigen if Cholesky fails (subspace collapse)
                @show "collapse at step $step, falling back to eigen decomposition"
                F = eigen(Hermitian(G))
                vals = F.values
                vecs = F.vectors

                new_states = Vector{MPS}(undef, k)
                for i in 1:k
                    if vals[i] > 1e-14
                        coef = vecs[:, i] / sqrt(vals[i])
                        ψ_new = coef[1] * states[1]
                        for j in 2:k
                            ψ_new = ψ_new + coef[j] * states[j]
                        end
                        new_states[i] = truncate(ψ_new; cutoff = cutoff, maxdim = maxdim)
                    else
                        new_states[i] = productMPS(sites, ["0" for _ in 1:N])
                    end
                end
                states = new_states
                spectrum_tlis[:, step] = -log.(sqrt.(abs.(vals[1:k])))
            else
                rethrow(e)
            end
        end
    end

    return spectrum_tlis
end

"""
    mps_measurement_enumeration(model::AnyonModel, ψ::MPS, sites::Vector{<:Index}, measurement_sites::Vector{Int}, 
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
function mps_measurement_enumeration(
    model::AnyonModel{AT},
    ψ::MPS,
    sites::Vector{<:Index},
    measurement_sites::Vector{Int},
    τ::Float64;
    cutoff::Float64 = 1e-10,
    maxdim::Int = 100,
) where {AT<:AbstractAnyonType}
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
            ψ_p, prob_p = measuremap(
                model,
                state,
                sites,
                site,
                τ,
                false;
                cutoff = cutoff,
                maxdim = maxdim,
            )
            # if prob_p > 1e-12
            new_trajectory_p = [current_trajectory; false]
            new_prob_p = current_prob * prob_p
            push!(next_level_states, ψ_p)
            push!(next_level_trajectories, new_trajectory_p)
            push!(next_level_probabilities, new_prob_p)
            # end

            # Apply 1 measurement
            ψ_m, prob_m = measuremap(
                model,
                state,
                sites,
                site,
                τ,
                true;
                cutoff = cutoff,
                maxdim = maxdim,
            )
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
    add_reference_qubits(model::AnyonModel, ψ::MPS, sites::Vector{<:Index}, site_idx::Int=1;
                         k_new::Int=1, verbose::Bool=false) -> Tuple{MPS, Vector{<:Index}}

Add reference qubits to an MPS state at specified site using copy method for correlation measurements.

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
    |ψ⟩ → α|00⟩ + β|11⟩  (reference ⊗ system)
```

This preserves the original state's structure while creating classical correlation.

For reset method (measure system qubit first, then create Bell pair), use 
[`add_reference_qubits_reset`](@ref) instead.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `ψ::MPS`: Input MPS quantum state
- `sites`: ITensor site indices
- `site_idx::Int=1`: Site index to entangle with reference qubit (1 ≤ site_idx ≤ N)
- `k_new::Int=1`: Number of reference qubits to add (must be 0 or 1)
- `verbose::Bool=false`: Enable verbose output for debugging

# Returns
- `MPS`: New MPS with reference qubit added at the leftmost position
- `Vector{Index}`: Updated site indices including the new reference qubit

# Notes
- Each system qubit can only be entangled with one reference qubit
- To add multiple reference qubits at different sites, call this function multiple times
- The reference qubit is inserted at position 1 (leftmost) of the MPS
- Original MPS is not modified (returns a new copy)

# Example
```julia
model = AnyonModel(FibonacciAnyon(), 6)
ψ, sites = anyon_mps_gst(model)
ψ_ref, sites_ref = add_reference_qubits(model, ψ, sites, 3)  # Add ref qubit at site 3
```

See also: [`add_reference_qubits_reset`](@ref)
"""
function add_reference_qubits(
    model::AnyonModel{AT},
    ψ::MPS,
    sites::Vector{<:Index},
    site_idx::Int = 1;
    k_new::Int = 1,
    verbose::Bool = false,
) where {AT<:AbstractAnyonType}
    1 ≤ site_idx ≤ length(ψ) ||
        error("site_idx must be in range [1, $(length(ψ))], got $site_idx")
    0 ≤ k_new ≤ 1 || error("k_new must be 0 or 1, got $k_new")

    k_new == 0 && return copy(ψ), copy(sites)  # No-op if k_new == 0

    N = length(sites)

    # Create new reference qubit site index
    ref_site = Index(2, "Qubit,Site,n=ref")
    new_sites = vcat([ref_site], sites)

    # Build new MPS with reference qubit at position 1
    ψ_new = MPS(N + 1)

    # Create reference qubit tensor initialized to |+⟩ = (|0⟩ + |1⟩)/√2
    # with a link index connecting to the rest of MPS
    link_idx = Index(1, "Link,l=0")
    ref_tensor = ITensor(ref_site, link_idx)
    ref_tensor[ref_site=>1, link_idx=>1] = 1.0  # Start with |0⟩
    ψ_new[1] = ref_tensor

    # Copy rest of MPS, adjusting first tensor to have the link
    first_tensor = ψ[1]
    first_inds = inds(first_tensor)
    new_first_tensor = first_tensor * ITensor([1.0], link_idx)
    ψ_new[2] = new_first_tensor

    for i = 2:N
        ψ_new[i+1] = ψ[i]
    end

    # Move reference qubit next to target site for applying copy gate
    target_pos = site_idx + 1  # +1 because we added reference at position 1

    # Swap reference qubit (at pos 1) towards target position
    for i = 1:(target_pos-2)
        # Swap sites i and i+1
        orthogonalize!(ψ_new, i)
        s_i = new_sites[i]
        s_ip1 = new_sites[i+1]

        SWAP = ITensor(s_i, s_i', s_ip1, s_ip1')
        for a = 1:dim(s_i), b = 1:dim(s_ip1)
            SWAP[s_i=>a, s_i'=>b, s_ip1=>b, s_ip1'=>a] = 1.0
        end
        ψ_new = apply(SWAP, ψ_new; cutoff = 1e-15)
        new_sites[i], new_sites[i+1] = new_sites[i+1], new_sites[i]
    end

    # Now reference is adjacent to target (at position target_pos - 1)
    ref_pos = target_pos - 1
    s_ref = new_sites[ref_pos]
    s_tgt = new_sites[target_pos]

    # Create copy gate: target controls reference
    # |t=0,r=0⟩ → |t=0,r=0⟩
    # |t=1,r=0⟩ → |t=1,r=1⟩  (copy 1 to reference)
    # |t=0,r=1⟩ → |t=0,r=1⟩
    # |t=1,r=1⟩ → |t=1,r=0⟩  (XOR)
    CopyGate = ITensor(s_tgt, s_tgt', s_ref, s_ref')
    CopyGate[s_tgt=>1, s_tgt'=>1, s_ref=>1, s_ref'=>1] = 1.0  # |0,0⟩ → |0,0⟩
    CopyGate[s_tgt=>2, s_tgt'=>2, s_ref=>1, s_ref'=>2] = 1.0  # |1,0⟩ → |1,1⟩
    CopyGate[s_tgt=>1, s_tgt'=>1, s_ref=>2, s_ref'=>2] = 1.0  # |0,1⟩ → |0,1⟩
    CopyGate[s_tgt=>2, s_tgt'=>2, s_ref=>2, s_ref'=>1] = 1.0  # |1,1⟩ → |1,0⟩

    orthogonalize!(ψ_new, ref_pos)
    ψ_new = apply(CopyGate, ψ_new; cutoff = 1e-15)

    # Swap reference back to position 1
    for i = (ref_pos-1):-1:1
        orthogonalize!(ψ_new, i+1)
        s_i = new_sites[i]
        s_ip1 = new_sites[i+1]

        SWAP = ITensor(s_i, s_i', s_ip1, s_ip1')
        for a = 1:dim(s_i), b = 1:dim(s_ip1)
            SWAP[s_i=>a, s_i'=>b, s_ip1=>b, s_ip1'=>a] = 1.0
        end
        ψ_new = apply(SWAP, ψ_new; cutoff = 1e-15)
        new_sites[i], new_sites[i+1] = new_sites[i+1], new_sites[i]
    end

    verbose && @info "Added reference qubit at position 1, entangled with site $site_idx"
    verbose && @show maxlinkdim(ψ_new)

    return ψ_new, new_sites
end

"""
    add_reference_qubits_reset(model::AnyonModel, ψ::MPS, sites::Vector{<:Index}, site_idx::Int=1;
                               k_new::Int=1, verbose::Bool=false) -> Tuple{Float64, Float64, MPS, MPS, Vector{Index}}

Add reference qubits using reset method: first measures/resets the system qubit, then creates a Bell pair.

# Principle

1. Project the system qubit onto |0⟩ or |1⟩ basis (strong measurement)
2. Create a Bell pair between the reset qubit and reference qubit
3. Returns both branches with their probabilities

This method is useful when you want to prepare a known initial state before 
adding the reference qubit.

# Arguments
- `model::AnyonModel`: Anyon model containing system parameters
- `ψ::MPS`: Input MPS quantum state
- `sites`: ITensor site indices
- `site_idx::Int=1`: Site index for reference qubit insertion (1 ≤ site_idx ≤ N)
- `k_new::Int=1`: Number of new reference qubits to add (must be 0 or 1)
- `verbose::Bool=false`: Enable verbose output for debugging

# Returns
- `Float64`: Probability of measuring |0⟩ at site_idx
- `Float64`: Probability of measuring |1⟩ at site_idx
- `MPS`: State after measuring |0⟩ and adding Bell pair
- `MPS`: State after measuring |1⟩ and adding Bell pair
- `Vector{Index}`: Updated site indices including the new reference qubit

# Notes
- Unlike `add_reference_qubits`, this method returns both measurement branches
- The probabilities satisfy prob₀ + prob₁ = 1
- Each returned MPS is normalized

# Example
```julia
model = AnyonModel(FibonacciAnyon(), 6)
ψ, sites = anyon_mps_gst(model)
prob0, prob1, ψ0, ψ1, sites_ref = add_reference_qubits_reset(model, ψ, sites, 3)
# Use ψ0 with probability prob0, or ψ1 with probability prob1
```

See also: [`add_reference_qubits`](@ref)
"""
function add_reference_qubits_reset(
    model::AnyonModel{AT},
    ψ::MPS,
    sites::Vector{<:Index},
    site_idx::Int = 1;
    k_new::Int = 1,
    verbose::Bool = false,
) where {AT<:AbstractAnyonType}
    1 ≤ site_idx ≤ length(ψ) ||
        error("site_idx must be in range [1, $(length(ψ))], got $site_idx")
    0 ≤ k_new ≤ 1 || error("k_new must be 0 or 1, got $k_new")

    if k_new == 0
        # No reference qubit to add, just return identity
        return 1.0, 0.0, copy(ψ), copy(ψ), copy(sites)
    end

    N = length(sites)

    # Step 1: Measure the system qubit at site_idx in Z basis
    # Calculate probabilities for |0⟩ and |1⟩ outcomes
    ψ_orth = orthogonalize(ψ, site_idx)
    site_tensor = ψ_orth[site_idx]
    s = sites[site_idx]

    # Project onto |0⟩
    proj0 = ITensor(s)
    proj0[s=>1] = 1.0

    # Project onto |1⟩  
    proj1 = ITensor(s)
    proj1[s=>2] = 1.0

    # Calculate probabilities using tensor contraction
    ψ0_proj = copy(ψ_orth)
    ψ1_proj = copy(ψ_orth)

    # Apply projectors
    ψ0_proj[site_idx] = noprime(site_tensor * proj0 * dag(prime(proj0, s)))
    ψ1_proj[site_idx] = noprime(site_tensor * proj1 * dag(prime(proj1, s)))

    prob0 = real(inner(ψ0_proj, ψ0_proj))
    prob1 = real(inner(ψ1_proj, ψ1_proj))

    # Normalize probabilities (should sum to 1)
    total_prob = prob0 + prob1
    prob0 /= total_prob
    prob1 /= total_prob

    verbose && @info "Measurement probabilities: P(0) = $prob0, P(1) = $prob1"

    # Normalize the projected states
    normalize!(ψ0_proj)
    normalize!(ψ1_proj)

    # Step 2: Add reference qubit and create Bell pair for each branch
    # Create reference site
    ref_site = Index(2, "Qubit,Site,n=ref")
    new_sites = vcat([ref_site], sites)

    # Helper function to add ref and create Bell pair
    function create_bell_pair_state(ψ_collapsed::MPS, measured_val::Int)
        ψ_bell = MPS(N + 1)

        # Create Bell pair: (|00⟩ + |11⟩)/√2 for measured 0
        #                   (|01⟩ + |10⟩)/√2 for measured 1
        link_idx = Index(2, "Link,l=0")

        ref_tensor = ITensor(ref_site, link_idx)
        if measured_val == 0
            # |Φ+⟩ = (|0⟩|0⟩ + |1⟩|1⟩)/√2
            ref_tensor[ref_site=>1, link_idx=>1] = 1/sqrt(2)  # |0⟩ paired with sys |0⟩
            ref_tensor[ref_site=>2, link_idx=>2] = 1/sqrt(2)  # |1⟩ paired with sys |1⟩
        else
            # |Ψ+⟩ = (|0⟩|1⟩ + |1⟩|0⟩)/√2
            ref_tensor[ref_site=>1, link_idx=>2] = 1/sqrt(2)  # |0⟩ paired with sys |1⟩
            ref_tensor[ref_site=>2, link_idx=>1] = 1/sqrt(2)  # |1⟩ paired with sys |0⟩
        end
        ψ_bell[1] = ref_tensor

        # Connect to rest of MPS
        first_tensor = ψ_collapsed[1]
        # Add link index to first tensor
        entangle_tensor = ITensor(link_idx, sites[1], sites[1]')
        # link=1 -> |0⟩, link=2 -> |1⟩
        entangle_tensor[link_idx=>1, sites[1]=>1, sites[1]'=>1] = 1.0
        entangle_tensor[link_idx=>2, sites[1]=>2, sites[1]'=>2] = 1.0

        new_first = noprime(first_tensor * entangle_tensor)
        ψ_bell[2] = new_first

        for i = 2:N
            ψ_bell[i+1] = ψ_collapsed[i]
        end

        normalize!(ψ_bell)
        return ψ_bell
    end

    ψ0_bell = create_bell_pair_state(ψ0_proj, 0)
    ψ1_bell = create_bell_pair_state(ψ1_proj, 1)

    verbose && @info "Created Bell pairs for both measurement branches"
    verbose && @show maxlinkdim(ψ0_bell), maxlinkdim(ψ1_bell)

    return prob0, prob1, ψ0_bell, ψ1_bell, new_sites
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
        # SWAP gate
        apply!(
            ψ,
            [min(cur, nxt), max(cur, nxt)],
            op("SWAP", site_type(ψ, cur), site_type(ψ, nxt)),
        )
        cur = nxt
    end
end

# MPS method for bulk_evolution - see Measurement.jl for full docstring
function bulk_evolution(
    model::AnyonModel{AT},
    sites::Vector{<:Index},
    state::MPS,
    measure_config::MeasureConfig,
    samples::Union{Nothing,BitMatrix} = nothing;
    normalized::Bool = true,
) where {AT<:AbstractAnyonType}

    # ---------- Sample decided according to mode ----------
    mode = measure_config.mode
    mode ∈ (:sample, :Born) || error("mode must be one of :sample, :Born")

    cutoff = measure_config.cutoff
    maxdim = measure_config.maxdim
    current_state = copy(state)
    if mode == :Born
        return _born_measure_mps(
            model,
            sites,
            current_state,
            measure_config;
            cutoff = cutoff,
            maxdim = maxdim,
        )
    elseif mode == :sample
        return _sample_measure_mps(
            model,
            sites,
            current_state,
            samples,
            measure_config;
            cutoff = cutoff,
            maxdim = maxdim,
            normalized = normalized,
        )
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
  - `state::MPS`: The final MPS state.
  - `samples::BitMatrix`: The generated measurement outcome sequences.
  - `free_energys::Vector{Float32}`: The free energy for each measurement layer.
  - `entanglement_entropys::Vector{Float32}`: Half-chain entanglement entropy at each period.
"""
function _born_measure_mps(
    model::AnyonModel{AT},
    sites::Vector{<:Index},
    current_state::MPS,
    measure_config::MeasureConfig;
    cutoff::Float64 = 1e-10,
    maxdim::Int = 100,
) where {AT<:AbstractAnyonType}

    n_cols = _samples_per_layer(model)  # Use max samples per layer
    τ = measure_config.τ
    t₁ = measure_config.t₁
    t₂ = measure_config.t₂
    rng = measure_config.rng
    enable_τ_eff = measure_config.enable_τ_eff
    verbose = measure_config.verbose

    Δt = t₂ - t₁ + 1
    Δt >= 0 || error("t₂ must be >= t₁")

    n_layers = layers_per_period(model.anyon_type)
    D = Δt * n_layers  # total number of layers
    N = length(sites)

    # 1. Initialize sample matrix with max columns per layer
    samples = BitMatrix(zeros(Bool, D, n_cols))
    sample_free_energy = zeros(Float32, D)
    entanglement_entropys = zeros(Float32, Δt)

    for period = 1:Δt
        # Apply all layers in this period
        for layer = 1:n_layers
            global_layer_idx = (period - 1) * n_layers + layer
            # Apply τ_eff only on the last layer of the last period
            τ_current = (period == Δt && layer == n_layers && enable_τ_eff) ? τ/2 : τ

            outcome = _stochastic_measurement_layer_mps_mps(
                model,
                τ_current,
                sites,
                current_state,
                rng,
                global_layer_idx;
                cutoff = cutoff,
                maxdim = maxdim,
                verbose = verbose,
                truncate_every_events = measure_config.truncate_every_events,
            )
            current_state = outcome.state

            # Write samples to correct column indices for this layer
            col_indices = _get_sample_column_indices(model, global_layer_idx)
            samples[global_layer_idx, col_indices] = outcome.sample
            sample_free_energy[global_layer_idx] = outcome.free_energy
        end
        # Compute half-chain EE on-the-fly
        entanglement_entropys[period] = Float32(ee_mps(current_state, div(N, 2)))
    end

    return Measurement_outcome_mps_bulk(
        current_state,
        samples,
        sample_free_energy,
        entanglement_entropys,
    )
end

"""
    _sample_measure_mps(model, sites::Vector{<:Index}, current_state, samples, measure_config; ...)

Evolve an MPS state using a predefined measurement trajectory.

This internal helper function is called by `bulk_evolution` when `mode` is `:sample`.

# Arguments
- `model::AnyonModel`: The anyon model.
- `sites::Vector{<:Index}`: ITensor site indices.
- `current_state::MPS`: The initial MPS state.
- `samples::BitMatrix`: The predefined measurement outcomes.
- `measure_config::MeasureConfig`: Configuration containing `τ`, `t₁`, `t₂`, etc.
- `cutoff::Float64=1e-10`: MPS truncation cutoff.
- `maxdim::Int=100`: Maximum bond dimension.

# Returns
- `Measurement_outcome_mps_bulk`: A struct containing:
  - `state::MPS`: The final MPS state.
  - `samples::BitMatrix`: The input measurement outcome sequences.
  - `free_energys::Vector{Float32}`: The free energy for each measurement layer.
  - `entanglement_entropys::Vector{Float32}`: Half-chain entanglement entropy at each period.
"""
function _sample_measure_mps(
    model::AnyonModel{AT},
    sites::Vector{<:Index},
    current_state::MPS,
    samples::BitMatrix,
    measure_config::MeasureConfig;
    cutoff::Float64 = 1e-10,
    maxdim::Int = 100,
    normalized::Bool = true,
) where {AT<:AbstractAnyonType}

    n_cols = _samples_per_layer(model)  # Use max samples per layer
    τ = measure_config.τ
    t₁ = measure_config.t₁
    t₂ = measure_config.t₂
    enable_τ_eff = measure_config.enable_τ_eff

    Δt = t₂ - t₁ + 1
    Δt >= 0 || error("t₂ must be >= t₁")

    n_layers = layers_per_period(model.anyon_type)
    D = Δt * n_layers  # total number of layers

    sample_free_energy = zeros(Float32, D)
    N = length(sites)
    entanglement_entropys = zeros(Float32, Δt)

    # Validate sample matrix dimensions
    isnothing(samples) && error("When mode=:sample samples must be ::BitMatrix")
    size(samples) == (D, n_cols) ||
        error("sample size should be ($D, $n_cols), got $(size(samples))")

    for period = 1:Δt
        # Apply all layers in this period
        for layer = 1:n_layers
            global_layer_idx = (period - 1) * n_layers + layer
            # Apply τ_eff only on the last layer of the last period
            τ_current = (period == Δt && layer == n_layers && enable_τ_eff) ? τ/2 : τ

            # Read samples from correct column indices for this layer
            col_indices = _get_sample_column_indices(model, global_layer_idx)
            layer_sample = BitVector(samples[global_layer_idx, col_indices])

            outcome = _apply_measurement_layer_mps(
                model,
                τ_current,
                sites,
                current_state,
                layer_sample,
                global_layer_idx;
                cutoff = cutoff,
                maxdim = maxdim,
                truncate_every_events = measure_config.truncate_every_events,
                normalized = normalized,
            )
            current_state = outcome.state
            sample_free_energy[global_layer_idx] = outcome.free_energy
        end
        # Compute half-chain EE on-the-fly
        entanglement_entropys[period] = Float32(ee_mps(current_state, div(N, 2)))
    end

    return Measurement_outcome_mps_bulk(
        current_state,
        samples,
        sample_free_energy,
        entanglement_entropys,
    )
end

struct Measurement_outcome_mps_bulk
    state::MPS
    samples::BitMatrix
    free_energys::Vector{Float32}
    entanglement_entropys::Vector{Float32}
end

struct Measurement_outcome_mps_boundary
    state::MPS
    sample::BitVector
    free_energy::Float32
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
function reference_evolution(
    model::AnyonModel{AT},
    sites::Vector{<:Index},
    forward::MPS,
    measure_config::MeasureConfig,
    sample::BitMatrix,
) where {AT<:AbstractAnyonType}

    N = model.N
    τ = measure_config.τ
    t₁ = measure_config.t₁
    t₂ = measure_config.t₂
    x₁ = measure_config.x₁
    x₂ = measure_config.x₂
    rng = measure_config.rng
    verbose = measure_config.verbose
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
    state = forward
    sample_layer = BitMatrix(undef, (size(sample, 1), n_measure))
    sample_free_energy = zeros(Float64, D)
    final_state = state  # will be overwritten


    if δt > 0 && δx > 0 # 3 ref qubits, both spatial and temporal correlation, actually 3-point correlation.
        verbose && @info "t₁ = $(t₁), t₂ = $(t₂), x₁ = $(x₁), x₂ = $(x₂), 3 refs"

        state1, sites1 = add_reference_qubits(model, state, sites, x₁; verbose = verbose)
        state2, sites2 = add_reference_qubits(model, state1, sites1, x₂; verbose = verbose)

        config1 = MeasureConfig(
            τ = τ,
            t₂ = (t₂-t₁),
            rng = rng,
            mode = mode,
            t₁ = 1,
            verbose = verbose,
            enable_τ_eff = false,
        )
        mo1 = bulk_evolution(model, sites2, state2, config1, sample[(2*t₁+1):(2*t₂), :])

        state3, sites3 =
            add_reference_qubits(model, mo1.state, sites2, x₂; verbose = verbose)

        config2 = MeasureConfig(
            τ = τ,
            t₂ = (Δt-t₂),
            rng = rng,
            mode = mode,
            t₁ = 1,
            verbose = verbose,
            enable_τ_eff = true,
        )
        mo2 = bulk_evolution(model, sites3, state3, config2, sample[(2*t₂+1):end, :])

        final_state = mo2.state

        sample_layer[(2*t₁+1):(2*t₂), :] .= mo1.samples
        sample_layer[(2*t₂+1):end, :] .= mo2.samples
        sample_free_energy[(2*t₁+1):(2*t₂)] .= mo1.free_energys
        sample_free_energy[(2*t₂+1):end] .= mo2.free_energys

    elseif δt == 0 # 2 ref qubits, pure 2-point spatial correlation
        verbose &&
            @info "x₁ = $(x₁), x₂ = $(x₂), δx = $(δx), at time slice t₁ = t₂ = $(t₁), 2 refs"

        state1, sites1 = add_reference_qubits(model, state, sites, x₁; verbose = verbose)
        state2, sites2 = add_reference_qubits(model, state1, sites1, x₂; verbose = verbose)

        config2 = MeasureConfig(
            τ = τ,
            t₂ = (Δt-t₁),
            rng = rng,
            mode = mode,
            t₁ = 1,
            verbose = verbose,
            enable_τ_eff = true,
        )
        mo2 = bulk_evolution(model, sites2, state2, config2, sample[(2*t₁+1):end, :])

        final_state = mo2.state
        sample_layer[(2*t₁+1):end, :] .= mo2.samples
        sample_free_energy[(2*t₁+1):end] .= mo2.free_energys

    elseif δx == 0 # 2 ref qubits, pure 2-point temporal correlation
        verbose &&
            @info "t₁ = $(t₁), t₂ = $(t₂), δt = $(δt), at site x₁ = x₂ = $(x₂), 2 refs"

        state1, sites1 = add_reference_qubits(model, state, sites, x₂; verbose = verbose)

        config1 = MeasureConfig(
            τ = τ,
            t₂ = (t₂-t₁),
            rng = rng,
            mode = mode,
            t₁ = 1,
            verbose = verbose,
            enable_τ_eff = false,
        )
        mo1 = bulk_evolution(model, sites1, state1, config1, sample[(2*t₁+1):(2*t₂), :])

        state2, sites2 =
            add_reference_qubits(model, mo1.state, sites1, x₂; verbose = verbose)

        config2 = MeasureConfig(
            τ = τ,
            t₂ = (Δt-t₂),
            rng = rng,
            mode = mode,
            t₁ = 1,
            verbose = verbose,
            enable_τ_eff = true,
        )
        mo2 = bulk_evolution(model, sites2, state2, config2, sample[(2*t₂+1):end, :])

        final_state = mo2.state
        sample_layer[(2*t₁+1):(2*t₂), :] .= mo1.samples
        sample_layer[(2*t₂+1):end, :] .= mo2.samples
        sample_free_energy[(2*t₁+1):(2*t₂)] .= mo1.free_energys
        sample_free_energy[(2*t₂+1):end] .= mo2.free_energys
    end

    return final_state, sample_layer, sample_free_energy
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
    for n = 1:dim(S, 1)
        p = S[n, n]^2
        if p > 1e-12
            SvN -= p * log(p)
        end
    end

    if abs(SvN) < 1e-14
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
    splitlis=Vector(1:(N-1))
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        EE_lis[m]=ee_mps(ψ, splitlis[m])
    end
    return EE_lis
end
