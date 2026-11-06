using FibonacciChain
using JLD2
using LinearAlgebra
using Random

"""
    simulate_staggered_ssb(L, λ_even, λ_odd, t; mode=:Born, seed=1,
                                        sign=true, samples=nothing)

Prepare |0> Π1|00…0> + |1> Πτ|00…0> with the original sector weights.
Each period applies odd sites 1:2:L with λ_odd, then even sites 2:2:L with λ_even.
Each event uses K_s = exp(-λ Π_i^s) / sqrt(1+exp(-2λ)) at its layer strength.
Here sign=true labels the suppressed vacuum channel, and false the τ channel.
Use mode=:Born to sample K_s, or :sample for a fixed record (default: all sign).
Time counts complete odd-then-even periods; no endpoint half-strength layer is added.
Fixed samples have size (2t,L÷2), with odd/even rows for odd/even bonds respectively.
Reference entropy is in nats and includes t=0. λ=Inf gives the exact projective limit.
The return value additionally carries final-state observables from
[`_staggered_final_observables`](@ref): per-reference-block subsystem entropy
scaling `final_subsystem_entropy`, Y expectations (`y_expectation`,
`y_expectation_by_reference`, `reference_weights`) and the τ-sector projected
occupation correlation `occupation_correlation_tau`.
"""
function simulate_staggered_ssb(
    L::Integer,
    λ_even::Real,
    λ_odd::Real,
    t::Integer;
    mode::Symbol = :Born,
    seed::Integer = 1,
    sign::Bool = true,
    samples::Union{Nothing,BitMatrix} = nothing,
)
    L >= 4 && iseven(L) || throw(ArgumentError("L must be even and at least 4"))
    all(λ -> !isnan(λ) && λ >= 0, (λ_even, λ_odd)) ||
        throw(ArgumentError("λ_even and λ_odd must be nonnegative (Inf is the projective limit)"))
    t >= 0 || throw(ArgumentError("t must be nonnegative"))
    mode ∈ (:Born, :sample) || throw(ArgumentError("mode must be :Born or :sample"))
    mode == :Born && samples !== nothing &&
        throw(ArgumentError("A fixed record requires mode=:sample"))
    if mode == :sample
        samples = isnothing(samples) ? BitMatrix(fill(sign, 2t, L÷2)) : copy(samples)
        size(samples) == (2t, L÷2) || throw(DimensionMismatch("samples must have size (2t,L÷2)"))
    end

    model = AnyonModel(FibonacciAnyon(), Int(L); pbc = true)
    basis = anyon_basis(model)
    initial_state = zeros(length(basis))
    initial_state[1] = 1
    ϕ = (1 + √5) / 2
    Y = topological_charge_operator(model)
    Ystate = Y * initial_state
    state_y₁ = (Ystate .+ inv(ϕ) .* initial_state) ./ (ϕ + inv(ϕ))
    state_yτ = initial_state .- state_y₁
    joint_state = vcat(state_y₁, state_yτ)
    joint_state ./= norm(joint_state)

    rng = MersenneTwister(seed)
    configs = map((λ_odd, λ_even)) do λ
        MeasureConfig(
            τ = isinf(λ) ? 1000.0 : Float64(λ),
            t₂ = Int(t),
            mode = mode,
            rng = rng,
            enable_τ_eff = false,
        )
    end
    record = falses(2t, L÷2)
    sample_free_energy = zeros(2t)
    reference_entropy = zeros(t+1)
    sector_probabilities = zeros(t+1, 2)
    ρ = reference_rdm(model, [1], joint_state)
    ρ ./= real(tr(ρ))
    reference_entropy[1] = ee(ρ)
    sector_probabilities[1, :] = real.(diag(ρ))

    for step in 1:t
        for layer in 1:2
            index = 2(step-1) + layer
            # The existing Kraus map favors its sign; complement it to suppress s.
            layer_sample = mode == :sample ? .!samples[index, :] : nothing
            # Library layer 2 is odd; layer 3 is even.
            outcome = reference_boundary_evolution(
                model, joint_state, configs[layer], layer_sample; layer_idx = index + 1,
            )
            joint_state = outcome.state
            all(isfinite, joint_state) || error("Invalid state at layer $index")
            record[index, :] = .!outcome.sample
            sample_free_energy[index] = outcome.free_energy
        end
        ρ = reference_rdm(model, [1], joint_state)
        ρ ./= real(tr(ρ))
        reference_entropy[step+1] = ee(ρ)
        sector_probabilities[step+1, :] = real.(diag(ρ))
    end

    final_observables = _staggered_final_observables(model, joint_state, Y)

    return (;
        L=Int(L), λ_even=Float64(λ_even), λ_odd=Float64(λ_odd),
        γ_even=tanh(Float64(λ_even)), γ_odd=tanh(Float64(λ_odd)),
        t=Int(t), times=collect(0:t), mode, seed=Int(seed), sign,
        measurement_sites_even=collect(2:2:L), measurement_sites_odd=collect(1:2:L),
        samples=record, sample_free_energy,
        initial_reference_entropy=first(reference_entropy), reference_entropy,
        sector_probabilities, reference_density_matrix=ρ, state=joint_state,
        initial_state="all_zero_fusion_path", reference_assignment="0:y=1; 1:y=tau",
        kraus_operator="exp(-lambda*Pi_i^s)/sqrt(1+exp(-2lambda))",
        sample_convention="true: suppressed vacuum channel; false: suppressed tau channel",
        time_unit="complete_odd_then_even_period", layer_order="odd_then_even", enable_τ_eff=false,
        final_observables...,
    )
end

"""
    _staggered_final_observables(model, joint_state, Y)

Final-state observables of a `simulate_staggered_ssb` joint state (system ⊗ one
reference qubit, blocks ordered 0:y=1, 1:y=τ):

- `final_subsystem_entropy[x, r]`: von Neumann entropy of the normalized
  reference block `r` restricted to sites `1:x`, for `x = 1:L-1` (subsystem
  scaling of each sector's conditional pure state; NaN for zero-weight blocks);
- `y_expectation_by_reference[r]`: `<Y>` within block `r`;
- `reference_weights[r]`: squared norms of the blocks (match the final
  `sector_probabilities` row);
- `y_expectation`: weight-averaged `<Y>` of the physical state;
- `occupation_correlation_tau`, `occupation_weight_tau`: the
  `staggered_projected_occupation_correlation` result in the τ sector.
"""
function _staggered_final_observables(model, joint_state, Y)
    L = model.N
    blocks = reshape(joint_state, :, 2)
    subsystem_sizes = collect(1:L-1)
    final_subsystem_entropy = fill(NaN, L - 1, 2)
    y_expectation_by_reference = fill(NaN, 2)
    reference_weights = zeros(2)
    for reference in 1:2
        ψ = blocks[:, reference]
        w = sum(abs2, ψ)
        reference_weights[reference] = w
        w > 0 || continue
        y_expectation_by_reference[reference] = real(dot(ψ, Y * ψ)) / w
        ψ_normalized = ψ ./ sqrt(w)
        for (index, x) in enumerate(subsystem_sizes)
            final_subsystem_entropy[index, reference] =
                ee(anyon_rdm(model, collect(1:x), ψ_normalized))
        end
    end
    y_expectation = sum(
        w * y for (w, y) in zip(reference_weights, y_expectation_by_reference) if !isnan(y)
    ) / sum(reference_weights)
    occupation = staggered_projected_occupation_correlation(L, joint_state; Y = Y)
    return (;
        subsystem_sizes, final_subsystem_entropy,
        y_expectation_by_reference, reference_weights, y_expectation,
        occupation_correlation_tau = occupation.correlation,
        occupation_weight_tau = occupation.weight_tau,
        occupation_operator = occupation.operator,
        occupation_ordering = occupation.ordering,
        final_observables_note = "blocks of the final joint state, ordered 0:y=1, 1:y=tau; entropies in nats",
    )
end

"""
    staggered_projected_occupation_correlation(L, state; sites=collect(1:L), sector=:tau, Y=nothing)
    staggered_projected_occupation_correlation(result; kwargs...)

Equal-time C[i,j] = <O_i O_j†> with O_i = Πτ n_i Π1, where n_i is the
fusion-path bit-1 occupation. Set sector=:one for O_i = Π1 n_i Πτ.
Accepts a physical state or the joint state with one reference qubit returned
by simulate_staggered_ssb. The reference is traced out; sector weights are
retained without conditioning. A supplied Y can be reused across trajectories.

Returns the complex correlation matrix, τ-sector weight, sites and operator
metadata. For an ensemble, average abs.(correlation) per trajectory, not the
absolute value of the ensemble-averaged complex matrix. The result wrapper
evaluates the saved final state without replaying its measurement record.
"""
function staggered_projected_occupation_correlation(
    L::Integer,
    state::AbstractVector;
    sites = collect(1:L),
    sector::Symbol = :tau,
    Y = nothing,
)
    L >= 4 && iseven(L) || throw(ArgumentError("L must be even and at least 4"))
    sector ∈ (:tau, :one) || throw(ArgumentError("sector must be :tau or :one"))
    !isempty(sites) && allunique(sites) &&
        all(i -> i isa Integer && 1 <= i <= L, sites) ||
        throw(ArgumentError("sites must be distinct integers in 1:L"))
    model = AnyonModel(FibonacciAnyon(), Int(L); pbc = true)
    basis = anyon_basis(model)
    dimension = length(basis)
    length(state) ∈ (dimension, 2dimension) ||
        throw(DimensionMismatch("Expected a physical state or a state with one reference qubit"))
    norm_squared = sum(abs2, state)
    all(isfinite, state) && isfinite(norm_squared) && norm_squared > 0 ||
        throw(ArgumentError("Invalid trajectory state"))
    Y = isnothing(Y) ? topological_charge_operator(model) : Y
    size(Y) == (dimension, dimension) || throw(DimensionMismatch("Invalid Y matrix"))
    occupations = [Float64(path[L-site+1]) for path in basis, site in sites]
    states = reshape(state, dimension, :)
    ϕ = (1 + √5) / 2
    denominator = ϕ + inv(ϕ)
    tau_states = (ϕ .* states .- Y * states) ./ denominator
    sector_states = sector == :tau ? tau_states : states .- tau_states
    frame = zeros(promote_type(eltype(state), eltype(Y), Float64), length(state), length(sites))
    for reference in axes(states, 2)
        # Each block contains O_i†|ψ_reference>; the Gram matrix traces the reference.
        block = occupations .* sector_states[:, reference]
        block = sector == :tau ? (Y * block .+ inv(ϕ) .* block) ./ denominator :
            (ϕ .* block .- Y * block) ./ denominator
        frame[(reference-1)*dimension+1:reference*dimension, :] = block
    end
    return (;
        correlation = frame' * frame / norm_squared,
        weight_tau = sum(abs2, tau_states) / norm_squared,
        sites = collect(sites), sector,
        operator = sector == :tau ? "O_i = Pi_tau n_i Pi_1; n_i = fusion-path bit 1 occupation" :
            "O_i = Pi_1 n_i Pi_tau; n_i = fusion-path bit 1 occupation",
        ordering = "C[i,j] = <O_sites[i] O_sites[j]^dagger>; reference traced without sector conditioning",
    )
end

staggered_projected_occupation_correlation(result::NamedTuple; kwargs...) =
    staggered_projected_occupation_correlation(result.L, result.state; kwargs...)

"""
    staggered_lyapunov_spectrum(L, λ_even, λ_odd, samples;
                               n_states=10, sector=:all, seed=0, rank_rtol=100eps(Float64))

Running QR estimates along a measurement record from [`simulate_staggered_ssb`](@ref).
Each period applies the odd layer, then the even layer, with the physical Kraus
normalization exp(-λ Π_i^s)/sqrt(1+exp(-2λ)). The record labels the suppressed
channel, so its bits are complemented when calling the native measurement map.

QR is performed after each local map, without sorting the frame. Finite strengths
are split into factors of strength at most one, with the normalization corrected
in the logarithmic stretches. Thus small nonzero directions are retained even
when a whole period is very ill-conditioned. Time still counts complete periods.
These are running QR estimates, not finite-time singular-value exponents.

Only `λ=Inf` permits rank loss. Its numerical rank is determined from the singular
values of the projected orthonormal frame, using `rank_rtol` relative to the unit
operator norm of an orthogonal projector. The leading independent QR columns are
retained; annihilated trailing directions have absorbing `-Inf` stretches and
exponents. `numerical_rank` records the retained frame size after each period,
capped by `n_states`. Rank decisions near the tolerance require a precision or
`rank_rtol` convergence check; they are not exact algebraic rank certificates.
An unresolved leading QR flag raises an error instead of propagating a spurious
frame. A zero product returns an empty frame and `-Inf` exponents.

`sector=:trivial` selects y=1 (Y eigenvalue ϕ); `:tau` selects y=τ (-1/ϕ).
The frame is projected into that sector initially and after each period. The
latter QR correction is included in the stretches; `sector_leakage` measures
leakage before correction. Every exact Fibonacci Kraus map commutes with Y.

`free_energy_spectrum` is `-lyapunov_exponents`. In the even-only projective,
fixed-record limits (λ_even=Inf, λ_odd=0), surviving rates approach -L*log(2)/4,
including the scalar odd-layer Kraus normalization. Random-record product ranks
are not zero-mode counts of a fixed projector Hamiltonian. Ensemble averages
must retain `-Inf`; finite-only conditional means describe a different ensemble.
"""
function staggered_lyapunov_spectrum(
    L::Integer,
    λ_even::Real,
    λ_odd::Real,
    samples::BitMatrix;
    n_states::Int = 10,
    sector::Symbol = :all,
    seed::Integer = 0,
    rank_rtol::Real = 100eps(Float64),
)
    L >= 4 && iseven(L) || throw(ArgumentError("L must be even and at least 4"))
    all(λ -> !isnan(λ) && λ >= 0, (λ_even, λ_odd)) ||
        throw(ArgumentError("λ_even and λ_odd must be nonnegative (Inf is the projective limit)"))
    isfinite(rank_rtol) && 0 < rank_rtol < 1 ||
        throw(ArgumentError("rank_rtol must lie strictly between zero and one"))
    sector ∈ (:all, :trivial, :tau) ||
        throw(ArgumentError("sector must be :all, :trivial (y=1) or :tau (y=τ)"))
    D, n_cols = size(samples)
    iseven(D) || throw(ArgumentError("samples must have an even number of rows (odd+even layers)"))
    periods = D ÷ 2
    n_cols == L ÷ 2 || throw(DimensionMismatch("samples must have $(L÷2) columns, got $n_cols"))

    model = AnyonModel(FibonacciAnyon(), Int(L); pbc = true)
    basis = anyon_basis(model)
    l = length(basis)

    P = nothing
    Y = nothing
    y_eigenvalue = NaN
    sector_dimension = l
    if sector != :all
        ϕ = (1 + √5) / 2
        y_eigenvalue = sector == :trivial ? ϕ : -inv(ϕ)
        y_other = sector == :trivial ? -inv(ϕ) : ϕ
        Y = topological_charge_operator(model)
        P = Symmetric((Y - y_other * I(l)) / (y_eigenvalue - y_other))
        sector_dimension = round(Int, tr(P))
    end
    1 <= n_states <= sector_dimension || throw(ArgumentError(
        "n_states must be between 1 and the $sector sector dimension $sector_dimension",
    ))
    k = n_states

    # A generic frame avoids the missing support of coordinate basis vectors.
    rng = MersenneTwister(seed)
    states = randn(rng, Float64, l, k)
    factor = qr(P === nothing ? states : P * states)
    minimum(abs.(diag(factor.R))) > 100eps(Float64) || error(
        "the random frame is rank deficient after $sector projection; try another seed",
    )
    states = Matrix(factor.Q)[:, 1:k]
    local_log_stretches = fill(-Inf, k, periods)
    numerical_rank = zeros(Int, periods)
    sector_leakage = sector == :all ? Float64[] : fill(NaN, periods)

    for step in 1:periods
        k_now = size(states, 2)
        k_now == 0 && break
        local_log_stretches[1:k_now, step] .= 0
        for (layer, λ) in enumerate((λ_odd, λ_even))
            if λ == 0
                local_log_stretches[1:k_now, step] .-= L * log(2) / 4
                continue
            end
            # Subdivision keeps every finite local factor well-conditioned and
            # avoids the native τ >= 100 projective sentinel for finite λ.
            n_substeps = isinf(λ) ? 1 : max(1, ceil(Int, λ))
            δ = isinf(λ) ? 1000.0 : Float64(λ / n_substeps)
            correction = isinf(λ) ? 0.0 :
                (log1p(exp(-2δ)) - log1p(exp(-2λ)) / n_substeps) / 2
            sites = layer == 1 ? (1:2:L) : (2:2:L)
            for (column, site) in enumerate(sites), _ in 1:n_substeps
                propagated = similar(states)
                FibonacciChain._measuremap_frame_impl!(
                    propagated, basis, model, δ, states, site,
                    !samples[2(step-1)+layer, column],
                )
                if isinf(λ)
                    # SVD determines rank only; rotating/sorting the QR flag
                    # would invalidate accumulated directional stretches.
                    retained = count(>(rank_rtol), svdvals(propagated))
                    local_log_stretches[retained+1:k_now, step] .= -Inf
                    k_now = retained
                    states = propagated[:, 1:k_now]
                    k_now == 0 && break
                end
                factor = qr(propagated[:, 1:k_now])
                stretches = abs.(diag(factor.R))
                all(isfinite, stretches) && minimum(stretches) > 0 ||
                    error("Numerical QR breakdown at period $step, site $site")
                if isinf(λ) && minimum(stretches) <= rank_rtol
                    error("Unresolved leading QR flag at period $step, site $site; increase precision or check rank_rtol")
                end
                local_log_stretches[1:k_now, step] .+= log.(stretches) .+ correction
                states = Matrix(factor.Q)[:, 1:k_now]
            end
            k_now == 0 && break
        end
        numerical_rank[step] = k_now
        if P !== nothing && k_now > 0
            sector_leakage[step] = norm(Y * states - y_eigenvalue * states) / norm(states)
            factor = qr(P * states)
            stretches = abs.(diag(factor.R))
            minimum(stretches) > 100eps(Float64) ||
                error("Numerical sector projection breakdown at period $step")
            local_log_stretches[1:k_now, step] .+= log.(stretches)
            states = Matrix(factor.Q)[:, 1:k_now]
        end
    end

    cumulative = cumsum(local_log_stretches; dims = 2)
    lyapunov_exponents = cumulative ./ reshape(collect(1:periods), 1, :)

    return (;
        local_log_stretches, lyapunov_exponents,
        free_energy_spectrum = -lyapunov_exponents,
        numerical_rank, rank_rtol = Float64(rank_rtol),
        frame_seed = Int(seed),
        method = "running_QR_estimates_per_local_factor",
        missing_value_convention = "-Inf: projective numerical rank loss; NaN is not a collapsed direction",
        sector, sector_dimension, sector_leakage, final_frame = states,
    )
end

"""
    staggered_lyapunov_spectrum(result; kwargs...)

Convenience wrapper taking the named tuple returned by
[`simulate_staggered_ssb`](@ref).
"""
staggered_lyapunov_spectrum(result::NamedTuple; kwargs...) =
    staggered_lyapunov_spectrum(
        result.L, result.λ_even, result.λ_odd, result.samples; kwargs...,
    )

"""
    simulate_staggered_sector_lyapunov(L, λ_even, λ_odd, t;
        sector=:trivial, n_states=10, trajectory_seed=1, seed=0,
        rank_rtol=100eps(Float64))

Generate a Born trajectory from the normalized projection of the all-τ fusion
path into one topological charge sector, then compute its sector-restricted
running QR spectrum. `sector=:trivial` selects y=1 (Y=ϕ), while `:tau` selects
y=τ (Y=-1/ϕ). Each period applies the odd layer followed by the even layer.

Returns `(trajectory, spectrum)`. The trajectory contains the final physical
state, samples in the same suppressed-channel convention as
[`simulate_staggered_ssb`](@ref), layer free energies, and the normalized Y
expectation after every period. `seed` controls the initial QR frame and
`trajectory_seed` controls Born sampling independently.
"""
function simulate_staggered_sector_lyapunov(
    L::Integer,
    λ_even::Real,
    λ_odd::Real,
    t::Integer;
    sector::Symbol = :trivial,
    n_states::Int = 10,
    trajectory_seed::Integer = 1,
    seed::Integer = 0,
    rank_rtol::Real = 100eps(Float64),
)
    L >= 4 && iseven(L) || throw(ArgumentError("L must be even and at least 4"))
    all(λ -> !isnan(λ) && λ >= 0, (λ_even, λ_odd)) ||
        throw(ArgumentError("λ_even and λ_odd must be nonnegative (Inf is the projective limit)"))
    t >= 1 || throw(ArgumentError("t must be positive"))
    sector ∈ (:trivial, :tau) ||
        throw(ArgumentError("sector must be :trivial (y=1) or :tau (y=τ)"))

    model = AnyonModel(FibonacciAnyon(), Int(L); pbc = true)
    basis = anyon_basis(model)
    initial_state = zeros(Float64, length(basis))
    initial_state[1] = 1
    ϕ = (1 + √5) / 2
    y_eigenvalue = sector == :trivial ? ϕ : -inv(ϕ)
    y_other = sector == :trivial ? -inv(ϕ) : ϕ
    Y = topological_charge_operator(model)
    state = (Y * initial_state .- y_other .* initial_state) ./ (y_eigenvalue - y_other)
    state_norm = norm(state)
    state_norm > 1e-14 || error("The all-τ path has zero weight in the $sector sector")
    state ./= state_norm

    rng = MersenneTwister(trajectory_seed)
    configs = map((λ_odd, λ_even)) do λ
        MeasureConfig(
            τ = isinf(λ) ? 1000.0 : Float64(λ),
            t₂ = Int(t),
            mode = :Born,
            rng = rng,
            enable_τ_eff = false,
        )
    end
    samples = falses(2t, L÷2)
    sample_free_energy = zeros(2t)
    y_expectation_values = zeros(t)

    for step in 1:t
        for layer in 1:2
            index = 2(step-1) + layer
            outcome = reference_boundary_evolution(
                model, state, configs[layer]; layer_idx = index + 1,
            )
            state = outcome.state
            all(isfinite, state) || error("Invalid state at layer $index")
            samples[index, :] = .!outcome.sample
            sample_free_energy[index] = outcome.free_energy
        end
        y_expectation_values[step] = real(dot(state, Y * state)) / sum(abs2, state)
    end

    spectrum = staggered_lyapunov_spectrum(
        L, λ_even, λ_odd, samples;
        n_states, sector, seed, rank_rtol,
    )
    maximum(abs.(y_expectation_values .- y_eigenvalue)) < 1e-5 ||
        error("the Born trajectory left the $sector sector")
    all(leakage -> isnan(leakage) || leakage < 1e-9, spectrum.sector_leakage) ||
        error("the Lyapunov frame left the $sector sector")

    trajectory = (;
        L = Int(L), λ_even = Float64(λ_even), λ_odd = Float64(λ_odd),
        γ_even = tanh(Float64(λ_even)), γ_odd = tanh(Float64(λ_odd)),
        t = Int(t), sector, y_eigenvalue, trajectory_seed = Int(trajectory_seed),
        state, samples, sample_free_energy, y_expectation_values,
        sample_convention = "true: suppressed vacuum channel; false: suppressed tau channel",
        time_unit = "complete_odd_then_even_period", layer_order = "odd_then_even",
    )
    return trajectory, spectrum
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 5 || error(
        "Usage: julia --project exm/Bulk_measure/staggered_ssb.jl L lambda_even lambda_odd t Born|sample [seed] [sign] [output.jld2]",
    )
    L = parse(Int, ARGS[1])
    λ_even = parse(Float64, ARGS[2])
    λ_odd = parse(Float64, ARGS[3])
    t = parse(Int, ARGS[4])
    mode = Symbol(ARGS[5])
    seed = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 1
    sign = length(ARGS) >= 7 ? parse(Int, ARGS[7]) == 1 : true
    branch = mode == :Born ? "Born" : sign ? "sample_s1" : "sample_stau"
    output = length(ARGS) >= 8 ? ARGS[8] : joinpath(
        @__DIR__, "../data/Bulk_measure/staggered_ssb", branch,
        "L$L", "lambda_even$(λ_even)_lambda_odd$(λ_odd)", "t$t", "trajectory_seed$seed.jld2",
    )
    ispath(output) && error("Output already exists: $output")
    result = simulate_staggered_ssb(L, λ_even, λ_odd, t; mode, seed, sign)
    n_states = min(10, length(anyon_basis(AnyonModel(FibonacciAnyon(), L; pbc = true))))
    spectrum = staggered_lyapunov_spectrum(result; n_states = n_states)
    result = (; result...,
        n_states,
        local_log_stretches = spectrum.local_log_stretches,
        lyapunov_exponents = spectrum.lyapunov_exponents,
        free_energy_spectrum = spectrum.free_energy_spectrum,
        numerical_rank = spectrum.numerical_rank,
        rank_rtol = spectrum.rank_rtol,
        frame_seed = spectrum.frame_seed,
        lyapunov_method = spectrum.method,
        missing_value_convention = spectrum.missing_value_convention,
    )
    mkpath(dirname(abspath(output)))
    JLD2.jldsave(output; result...)
    println("Saved $output")
    println("S_R(0)=$(first(result.reference_entropy)), S_R(t)=$(last(result.reference_entropy))")
    t > 0 && println("final running QR estimates: $(spectrum.lyapunov_exponents[:, end])")
end
