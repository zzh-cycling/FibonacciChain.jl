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
    Ystate = topological_charge_operator(model) * initial_state
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
    )
end

"""
    staggered_lyapunov_spectrum(L, λ_even, λ_odd, samples; n_states=10, sector=:all)

Finite-time Lyapunov spectrum of the staggered transfer-matrix product along a
measurement record `samples` produced by [`simulate_staggered_ssb`](@ref). The
frame is propagated with the unnormalized Kraus map: each period first applies
the odd-site layer with `λ_odd`, then the even-site layer with `λ_even`, using
the complemented record bits `.!samples[row, :]`, exactly reproducing the
Kraus operators of the trajectory. After every period the frame is
QR-orthogonalized, `T_t Q_{t-1} = Q_t R_t`, the local logarithmic stretches are
`log(abs(diag(R_t)))` (no per-step sorting, to avoid mixing Oseledec
directions), and the finite-time exponents are their running averages.

With `sector=:trivial` (`y=1`, eigenvalue `ϕ`) or `:tau` (`y=τ`, eigenvalue
`-1/ϕ`) the frame is projected into that topological charge sector initially
and after every period; this is exact because every Fibonacci Kraus layer
commutes with `Y`. `λ=Inf` gives the projective limit (τ=1000).

Returns a named tuple with `local_log_stretches`, `lyapunov_exponents`,
`free_energy_spectrum` (= `-lyapunov_exponents`, the `ln Z` convention of
`lyapunov_spectrum`), `sector_dimension`, `sector_leakage` and `final_frame`.

In the projective postselected limits the leading flat exponents count the
zero modes of the corresponding projector Hamiltonian: suppressing the vacuum
channel on even bonds only (λ_even=Inf, λ_odd=0, all-`true` record) gives
`Lucas(L/2)` flat exponents (triplet manifold), while suppressing the τ
channel (all-`false` record) gives 2 flat exponents (singlet SSB doublet).
"""
function staggered_lyapunov_spectrum(
    L::Integer,
    λ_even::Real,
    λ_odd::Real,
    samples::BitMatrix;
    n_states::Int = 10,
    sector::Symbol = :all,
)
    L >= 4 && iseven(L) || throw(ArgumentError("L must be even and at least 4"))
    sector ∈ (:all, :trivial, :tau) ||
        throw(ArgumentError("sector must be :all, :trivial (y=1) or :tau (y=τ)"))
    D, n_cols = size(samples)
    iseven(D) || throw(ArgumentError("samples must have an even number of rows (odd+even layers)"))
    periods = D ÷ 2
    n_cols == L ÷ 2 || throw(DimensionMismatch("samples must have $(L÷2) columns, got $n_cols"))
    τ_even = isinf(λ_even) ? 1000.0 : Float64(λ_even)
    τ_odd = isinf(λ_odd) ? 1000.0 : Float64(λ_odd)

    model = AnyonModel(FibonacciAnyon(), Int(L); pbc = true)
    basis = anyon_basis(model)
    l = length(basis)

    # Optional topological sector projector P = (Y - ȳI)/(y - ȳ)
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
    n_states >= 1 || throw(ArgumentError("n_states must be positive"))
    n_states <= sector_dimension || error(
        "requested n_states=$n_states, but the $sector sector dimension is $sector_dimension",
    )
    k = n_states

    # Initial frame: first k basis vectors, projected and orthonormalized
    states = zeros(Float64, l, k)
    for i in 1:k
        states[i, i] = 1.0
    end
    factor = qr(P === nothing ? states : P * states)
    minimum(abs.(diag(factor.R)[1:k])) > 1e-12 || error(
        "the first $k basis vectors are linearly dependent after $sector projection",
    )
    states = Matrix(factor.Q)[:, 1:k]

    # In the projective limit (λ=Inf) the transfer matrix is rank-deficient:
    # directions outside the zero-mode manifold are killed exactly. The frame
    # is truncated adaptively; collapsed directions keep NaN stretches (their
    # Lyapunov exponents are -∞), so the surviving flat exponents directly
    # count the zero modes of the corresponding projector Hamiltonian.
    local_log_stretches = fill(NaN, k, periods)
    active = collect(1:k)
    sector_leakage = sector == :all ? Float64[] : zeros(Float64, periods)

    for step in 1:periods
        k_now = length(active)
        propagated = similar(states)
        for i in 1:k_now
            # Odd-site layer with λ_odd: record row 2step-1; even library
            # layer_idx addresses odd sites (cf. simulate_staggered_ssb, which
            # uses layer_idx = index + 1 and complemented record bits).
            outcome = FibonacciChain._apply_measurement_layer(
                model, τ_odd, states[:, i],
                BitVector(.!samples[2step-1, :]); layer_idx = 2, normalized = false,
            )
            # Even-site layer with λ_even: record row 2step; odd layer_idx.
            outcome = FibonacciChain._apply_measurement_layer(
                model, τ_even, outcome.state,
                BitVector(.!samples[2step, :]); layer_idx = 3, normalized = false,
            )
            propagated[:, i] = outcome.state
        end

        factor = qr(propagated)
        stretches = abs.(diag(factor.R)[1:k_now])
        alive = stretches .> 1e-12 * max(maximum(stretches), 1.0)
        local_log_stretches[active[alive], step] = log.(stretches[alive])
        active = active[alive]
        isempty(active) && error("Lyapunov frame collapsed entirely at period $step")
        states = Matrix(factor.Q)[:, alive]
        if P !== nothing
            states = Matrix(qr(P * states).Q)[:, 1:length(active)]
            sector_leakage[step] = norm(Y * states - y_eigenvalue * states) / norm(states)
        end
    end

    cumulative = cumsum(local_log_stretches; dims = 2)
    lyapunov_exponents = cumulative ./ reshape(collect(1:periods), 1, :)

    return (;
        local_log_stretches, lyapunov_exponents,
        free_energy_spectrum = -lyapunov_exponents,
        sector, sector_dimension, sector_leakage, final_frame = states,
    )
end

"""
    staggered_lyapunov_spectrum(result; n_states=10, sector=:all)

Convenience wrapper taking the named tuple returned by
[`simulate_staggered_ssb`](@ref).
"""
staggered_lyapunov_spectrum(result::NamedTuple; n_states::Int = 10, sector::Symbol = :all) =
    staggered_lyapunov_spectrum(
        result.L, result.λ_even, result.λ_odd, result.samples; n_states, sector,
    )

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
    )
    mkpath(dirname(abspath(output)))
    JLD2.jldsave(output; result...)
    println("Saved $output")
    println("S_R(0)=$(first(result.reference_entropy)), S_R(t)=$(last(result.reference_entropy))")
    println("final Lyapunov exponents: $(spectrum.lyapunov_exponents[:, end])")
end
