using FibonacciChain
using JLD2
using LinearAlgebra
using Random
using Statistics

const HYBRID_DATA_ROOT =
    joinpath("exm", "data", "Bulk_measure", "hybrid_dynamics")

function hybrid_sector_projector(model::AnyonModel{FibonacciAnyon}, sector::Symbol)
    model.pbc || error("Topological sectors require periodic boundaries")
    sector ∈ (:trivial, :tau) || throw(
        ArgumentError("sector must be :trivial or :tau, got $sector"),
    )
    ϕ = (1 + √5) / 2
    y = sector == :trivial ? ϕ : -inv(ϕ)
    y_other = sector == :trivial ? -inv(ϕ) : ϕ
    Y = topological_charge_operator(model)
    dim = size(Y, 1)
    projector = Matrix(Symmetric((Y - y_other * I(dim)) / (y - y_other)))
    return projector, Y, y
end

"""
    hybrid_sector_state(model, sector) -> Vector

Project the all-`τ` fusion path into one topological sector. `sector=:trivial`
denotes `y=1` with `Y` eigenvalue `ϕ`; `sector=:tau` denotes `y=τ` with
eigenvalue `-1/ϕ`.
"""
function hybrid_sector_state(
    model::AnyonModel{FibonacciAnyon},
    sector::Symbol,
)
    projector, _, _ = hybrid_sector_projector(model, sector)
    dim = size(projector, 1)
    state = zeros(Float64, dim)
    state[1] = 1.0
    projected = projector * state
    norm(projected) > 1e-14 || error("all-τ state has zero weight in $sector")
    normalize!(projected)
    return projected
end

"""Convert measurement sharpness `γ ∈ [0,1]` to the repository parameter `τ`."""
function hybrid_measurement_strength(γ::Real)
    0 <= γ <= 1 || throw(ArgumentError("γ must satisfy 0 ≤ γ ≤ 1, got $γ"))
    return γ == 1 ? Inf : atanh(Float64(γ))
end

"""
    exclusive_hybrid_config(p, periods, seed; γ=1, mode=:Born,
                            random_angles=true, θ=π)

Configuration for the exclusive brick protocol. Every brick is independently
a measurement with probability `p`, or a unitary with probability `1-p`.
`γ=1` gives projectors exactly. With `random_angles=true`, every unitary has an
independent `θ[x,t] ~ Uniform(0,2π)`.
"""
function exclusive_hybrid_config(
    p::Real,
    periods::Integer,
    seed::Integer;
    γ::Real = 1.0,
    mode::Symbol = :Born,
    random_angles::Bool = true,
    θ::Real = π,
)
    periods >= 1 || throw(ArgumentError("periods must be positive"))
    return HybridConfig(
        τ = hybrid_measurement_strength(γ),
        t₂ = Int(periods),
        mode = mode,
        rng = MersenneTwister(seed),
        enable_τ_eff = false,
        p = p,
        θ = θ,
        random_angles = random_angles,
    )
end

function hybrid_quarter_tmi(model::AnyonModel, state::Vector)
    L = model.N
    L % 4 == 0 || error("quarter-chain I₃ requires L divisible by 4, got L=$L")
    q = L ÷ 4
    A = collect(1:q)
    B = collect(q + 1:2q)
    C = collect(2q + 1:3q)
    return tri_mutual_information(model, (A, B, C), state)
end

"""
    simulate_hybrid_sector_exact(L, p, periods, seed; sector=:trivial, γ=1)

Generate one exact Born trajectory in a fixed topological sector. The returned
data contain the full replayable circuit, half-chain entropy versus time,
final entropy profile, final quarter-chain `I₃`, and the final `Y` expectation.
"""
function simulate_hybrid_sector_exact(
    L::Integer,
    p::Real,
    periods::Integer,
    seed::Integer;
    sector::Symbol = :trivial,
    γ::Real = 1.0,
    random_angles::Bool = true,
    θ::Real = π,
)
    iseven(L) || error("The staggered PBC circuit requires even L, got $L")
    model = AnyonModel(FibonacciAnyon(), Int(L); pbc = true)
    initial_state = hybrid_sector_state(model, sector)
    config = exclusive_hybrid_config(
        p,
        periods,
        seed;
        γ = γ,
        random_angles = random_angles,
        θ = θ,
    )
    outcome = bulk_evolution(model, initial_state, config)
    Y = topological_charge_operator(model)
    final_y = real(dot(outcome.state, Y * outcome.state))
    expected_y = sector == :trivial ? (1 + √5) / 2 : -2 / (1 + √5)
    isapprox(final_y, expected_y; atol = 1e-8, rtol = 1e-8) ||
        error("hybrid trajectory left the $sector sector: ⟨Y⟩=$final_y")

    return (
        state = outcome.state,
        schedule = outcome.schedule,
        sample_free_energy = outcome.free_energys,
        halfchain_entropy = outcome.entanglement_entropys,
        final_entropy_profile = anyon_eelis(model, outcome.state),
        final_quarter_tmi = L % 4 == 0 ? hybrid_quarter_tmi(model, outcome.state) : NaN,
        final_y_expectation = final_y,
        sector = sector,
    )
end

"""
    simulate_hybrid_learnability_exact(L, p, periods, seed; γ=1,
                                       sampled_sector=:mixture)

Run the common-record two-hypothesis Bayesian evolution. The circuit and
record are sampled from the equal-prior sector mixture by default. All returned
time series are single-trajectory estimators and must be averaged over both
circuit randomness and records.
"""
function simulate_hybrid_learnability_exact(
    L::Integer,
    p::Real,
    periods::Integer,
    seed::Integer;
    γ::Real = 1.0,
    sampled_sector::Symbol = :mixture,
    random_angles::Bool = true,
    θ::Real = π,
)
    iseven(L) || error("The staggered PBC circuit requires even L, got $L")
    model = AnyonModel(FibonacciAnyon(), Int(L); pbc = true)
    state_y1 = hybrid_sector_state(model, :trivial)
    state_ytau = hybrid_sector_state(model, :tau)
    config = exclusive_hybrid_config(
        p,
        periods,
        seed;
        γ = γ,
        random_angles = random_angles,
        θ = θ,
    )
    # Keep the equal-prior draw independent of all circuit/Born random numbers.
    sector_rng = MersenneTwister(UInt(seed) ⊻ UInt(0x6a09e667))
    outcome = hybrid_bayesian_evolution(
        model,
        state_y1,
        state_ytau,
        config;
        sampled_sector = sampled_sector,
        sector_rng = sector_rng,
    )
    truth_state = outcome.sampled_sector == :tau ? outcome.state_ytau : outcome.state_y1

    return (
        outcome = outcome,
        truth_halfchain_entropy = ee(
            anyon_rdm(model, collect(1:(L ÷ 2)), truth_state),
        ),
        truth_quarter_tmi = L % 4 == 0 ? hybrid_quarter_tmi(model, truth_state) : NaN,
    )
end

"""
    simulate_hybrid_lyapunov_exact(L, p, periods, seed;
                                   sector=:trivial, n_states=10, γ=1)

Generate a Born circuit/record in one topological sector and compute the exact
Lyapunov spectrum of that same hybrid transfer-matrix product, including the
stored random unitary angle at every unitary brick.
"""
function simulate_hybrid_lyapunov_exact(
    L::Integer,
    p::Real,
    periods::Integer,
    seed::Integer;
    sector::Symbol = :trivial,
    n_states::Int = 10,
    γ::Real = 1.0,
    random_angles::Bool = true,
    θ::Real = π,
)
    trajectory = simulate_hybrid_sector_exact(
        L,
        p,
        periods,
        seed;
        sector = sector,
        γ = γ,
        random_angles = random_angles,
        θ = θ,
    )
    model = AnyonModel(FibonacciAnyon(), Int(L); pbc = true)
    config = exclusive_hybrid_config(
        p,
        periods,
        seed;
        γ = γ,
        mode = :sample,
        random_angles = random_angles,
        θ = θ,
    )
    spectrum = hybrid_lyapunov_spectrum(
        model,
        config,
        trajectory.schedule;
        n_states = n_states,
        sector = sector,
    )
    maximum(spectrum.sector_leakage) < 1e-9 ||
        error("hybrid Lyapunov frame left the $sector sector")
    return trajectory, spectrum
end

function _density_matrix_entropy(ρ::AbstractMatrix)
    eigenvalues = eigvals(Hermitian((ρ + ρ') / 2))
    return -sum(x -> x > 1e-14 ? x * log(x) : 0.0, real.(eigenvalues))
end

"""
    simulate_hybrid_purification_exact(L, p, periods, seed;
                                       sector=:trivial, γ=1)

Start from the maximally mixed state `ρ₀=P_y/d_y` in one topological sector and
sample the exclusive hybrid circuit directly with the mixed-state Born rule.
Returns the von Neumann entropy and purity after every period. This is an exact
dense-density-matrix diagnostic and is therefore intended for small `L`.
"""
function simulate_hybrid_purification_exact(
    L::Integer,
    p::Real,
    periods::Integer,
    seed::Integer;
    sector::Symbol = :trivial,
    γ::Real = 1.0,
    random_angles::Bool = true,
    θ::Real = π,
)
    iseven(L) || error("The staggered PBC circuit requires even L, got $L")
    model = AnyonModel(FibonacciAnyon(), Int(L); pbc = true)
    P, Y, expected_y = hybrid_sector_projector(model, sector)
    sector_dimension = round(Int, real(tr(P)))
    ρ = ComplexF64.(P / sector_dimension)
    config = exclusive_hybrid_config(
        p,
        periods,
        seed;
        γ = γ,
        random_angles = random_angles,
        θ = θ,
    )
    mc = config.measurement
    n_layers = FibonacciChain.layers_per_period(model)
    D = periods * n_layers
    ncols = L ÷ 2
    masks = falses(D, ncols)
    outcomes = falses(D, ncols)
    angles = fill(NaN, D, ncols)
    entropies = zeros(Float64, periods)
    purities = zeros(Float64, periods)
    y_expectations = zeros(Float64, periods)

    for period in 1:periods
        for layer in 1:n_layers
            glayer = (period - 1) * n_layers + layer
            event_sites, measure_model, strength =
                FibonacciChain._obtain_measurement_config(model, glayer, mc.τ)
            cols = FibonacciChain._get_sample_column_indices(model, glayer)
            for (k, site) in enumerate(event_sites)
                col = cols[k]
                is_measurement = config.p == 1 ? true :
                                 config.p == 0 ? false : rand(mc.rng) < config.p
                masks[glayer, col] = is_measurement
                if is_measurement
                    M0 = ComplexF64.(
                        FibonacciChain.measure_matrix(
                            measure_model,
                            strength,
                            site,
                            false,
                        ),
                    )
                    ρ0 = M0 * ρ * M0'
                    p0 = clamp(real(tr(ρ0)), 0.0, 1.0)
                    outcome = rand(mc.rng) >= p0
                    outcomes[glayer, col] = outcome
                    if outcome
                        M1 = ComplexF64.(
                            FibonacciChain.measure_matrix(
                                measure_model,
                                strength,
                                site,
                                true,
                            ),
                        )
                        ρnext = M1 * ρ * M1'
                        probability = real(tr(ρnext))
                    else
                        ρnext = ρ0
                        probability = p0
                    end
                    probability > 0 || error("sampled a zero-probability mixed-state outcome")
                    ρ = ρnext / probability
                else
                    θ_event = config.random_angles ? 2π * rand(mc.rng) : config.θ
                    angles[glayer, col] = θ_event
                    Ptau = ComplexF64.(
                        FibonacciChain.measure_matrix(
                            measure_model,
                            Inf,
                            site,
                            false,
                        ),
                    )
                    P1 = ComplexF64.(
                        FibonacciChain.measure_matrix(
                            measure_model,
                            Inf,
                            site,
                            true,
                        ),
                    )
                    U = P1 + cis(θ_event) * Ptau
                    ρ = U * ρ * U'
                end
            end
        end
        ρ = (ρ + ρ') / 2
        ρ ./= real(tr(ρ))
        entropies[period] = _density_matrix_entropy(ρ)
        purities[period] = real(tr(ρ * ρ))
        y_expectations[period] = real(tr(ρ * Y))
    end
    maximum(abs.(y_expectations .- expected_y)) < 1e-8 ||
        error("mixed-state trajectory left the $sector sector")

    return (
        state = ρ,
        schedule = HybridGateSchedule(masks, outcomes, angles),
        entropy = entropies,
        purity = purities,
        y_expectation = y_expectations,
        sector = sector,
        sector_dimension = sector_dimension,
    )
end

function _mean_and_stderr(values::AbstractMatrix)
    n = size(values, 2)
    average = vec(mean(values; dims = 2))
    stderr = n == 1 ? zeros(size(values, 1)) : vec(std(values; dims = 2)) / sqrt(n)
    return average, stderr
end

"""
    hybrid_sector_ensemble(L, p, periods, seed_start, seed_end;
                           sector=:trivial, γ=1)

Average fixed-sector entanglement observables over circuit and Born-record
randomness. This is the companion dataset for locating `p_c` independently of
the record-fidelity estimate of `p_#`.
"""
function hybrid_sector_ensemble(
    L::Integer,
    p::Real,
    periods::Integer,
    seed_start::Integer,
    seed_end::Integer;
    sector::Symbol = :trivial,
    γ::Real = 1.0,
    random_angles::Bool = true,
    θ::Real = π,
)
    seed_end >= seed_start || error("seed_end must be at least seed_start")
    seeds = collect(Int(seed_start):Int(seed_end))
    n = length(seeds)
    halfchain_entropy = zeros(Float64, periods, n)
    entropy_profiles = zeros(Float64, L - 1, n)
    quarter_tmi = zeros(Float64, n)
    y_expectation = zeros(Float64, n)
    measurement_fraction = zeros(Float64, n)

    for (column, seed) in enumerate(seeds)
        result = simulate_hybrid_sector_exact(
            L,
            p,
            periods,
            seed;
            sector = sector,
            γ = γ,
            random_angles = random_angles,
            θ = θ,
        )
        halfchain_entropy[:, column] = result.halfchain_entropy
        entropy_profiles[:, column] = result.final_entropy_profile
        quarter_tmi[column] = result.final_quarter_tmi
        y_expectation[column] = result.final_y_expectation
        measurement_fraction[column] =
            count(result.schedule.measurement_mask) / length(result.schedule.measurement_mask)
    end

    mean_halfchain, stderr_halfchain = _mean_and_stderr(halfchain_entropy)
    mean_profile, stderr_profile = _mean_and_stderr(entropy_profiles)
    return (
        seeds = seeds,
        halfchain_entropies = halfchain_entropy,
        mean_halfchain_entropy = mean_halfchain,
        stderr_halfchain_entropy = stderr_halfchain,
        final_entropy_profiles = entropy_profiles,
        mean_final_entropy_profile = mean_profile,
        stderr_final_entropy_profile = stderr_profile,
        final_quarter_tmis = quarter_tmi,
        mean_final_quarter_tmi = mean(quarter_tmi),
        stderr_final_quarter_tmi = n == 1 ? 0.0 : std(quarter_tmi) / sqrt(n),
        final_y_expectations = y_expectation,
        measurement_fractions = measurement_fraction,
        mean_measurement_fraction = mean(measurement_fraction),
    )
end

"""
    hybrid_learnability_ensemble(L, p, periods, seed_start, seed_end; γ=1)

Average the Ma–Vasseur record diagnostics over equal-prior mixture records.
`F_t` is the mean `record_fidelity_estimator`, while `H(Y|M_t)` is the mean
`conditional_entropy_estimator`. The function also returns the ensemble Bayes
error, mutual information, final truth-state half-chain entropy and final
`I₃`. `learning_time` is the first period satisfying `-log(F_t) ≥ 1`, or `0`
if the threshold is not reached.
"""
function hybrid_learnability_ensemble(
    L::Integer,
    p::Real,
    periods::Integer,
    seed_start::Integer,
    seed_end::Integer;
    γ::Real = 1.0,
    random_angles::Bool = true,
    θ::Real = π,
)
    seed_end >= seed_start || error("seed_end must be at least seed_start")
    seeds = collect(Int(seed_start):Int(seed_end))
    n = length(seeds)
    fidelity = zeros(Float64, periods, n)
    conditional_entropy = similar(fidelity)
    mutual_information = similar(fidelity)
    bayes_error = similar(fidelity)
    log_likelihood = similar(fidelity)
    halfchain_entropy = zeros(Float64, n)
    quarter_tmi = zeros(Float64, n)
    sampled_sectors = Vector{String}(undef, n)

    for (column, seed) in enumerate(seeds)
        result = simulate_hybrid_learnability_exact(
            L,
            p,
            periods,
            seed;
            γ = γ,
            random_angles = random_angles,
            θ = θ,
        )
        outcome = result.outcome
        fidelity[:, column] = outcome.record_fidelity_estimator
        conditional_entropy[:, column] = outcome.conditional_entropy_estimator
        mutual_information[:, column] = outcome.mutual_information_estimator
        bayes_error[:, column] = outcome.bayes_error_estimator
        log_likelihood[:, column] = outcome.log_likelihood_ratio
        halfchain_entropy[column] = result.truth_halfchain_entropy
        quarter_tmi[column] = result.truth_quarter_tmi
        sampled_sectors[column] = String(outcome.sampled_sector)
    end

    mean_fidelity, stderr_fidelity = _mean_and_stderr(fidelity)
    mean_conditional_entropy, stderr_conditional_entropy =
        _mean_and_stderr(conditional_entropy)
    mean_mutual_information, stderr_mutual_information =
        _mean_and_stderr(mutual_information)
    mean_bayes_error, stderr_bayes_error = _mean_and_stderr(bayes_error)
    crossing = findfirst(x -> x >= 1, -log.(mean_fidelity))

    return (
        seeds = seeds,
        sampled_sectors = sampled_sectors,
        fidelity_estimators = fidelity,
        conditional_entropy_estimators = conditional_entropy,
        mutual_information_estimators = mutual_information,
        bayes_error_estimators = bayes_error,
        log_likelihood_ratios = log_likelihood,
        mean_record_fidelity = mean_fidelity,
        stderr_record_fidelity = stderr_fidelity,
        mean_conditional_entropy = mean_conditional_entropy,
        stderr_conditional_entropy = stderr_conditional_entropy,
        mean_mutual_information = mean_mutual_information,
        stderr_mutual_information = stderr_mutual_information,
        mean_bayes_error = mean_bayes_error,
        stderr_bayes_error = stderr_bayes_error,
        final_halfchain_entropies = halfchain_entropy,
        final_quarter_tmis = quarter_tmi,
        mean_final_halfchain_entropy = mean(halfchain_entropy),
        stderr_final_halfchain_entropy =
            n == 1 ? 0.0 : std(halfchain_entropy) / sqrt(n),
        mean_final_quarter_tmi = mean(quarter_tmi),
        stderr_final_quarter_tmi = n == 1 ? 0.0 : std(quarter_tmi) / sqrt(n),
        learning_time = isnothing(crossing) ? 0 : crossing,
    )
end

function hybrid_output_dir(L::Integer, p::Real; γ::Real = 1.0)
    return joinpath(HYBRID_DATA_ROOT, "L$(L)", "p$(Float64(p))", "gamma$(Float64(γ))")
end

function save_hybrid_sector_trajectory(
    L::Integer,
    p::Real,
    periods::Integer,
    seed::Integer;
    sector::Symbol = :trivial,
    γ::Real = 1.0,
)
    result = simulate_hybrid_sector_exact(
        L,
        p,
        periods,
        seed;
        sector = sector,
        γ = γ,
    )
    output_dir = hybrid_output_dir(L, p; γ = γ)
    mkpath(output_dir)
    output_path = joinpath(
        output_dir,
        "sector_$(sector)_periods$(periods)_seed$(seed).jld2",
    )
    jldsave(
        output_path;
        L = Int(L),
        p = Float64(p),
        γ = Float64(γ),
        periods = Int(periods),
        seed = Int(seed),
        sector = String(sector),
        measurement_mask = result.schedule.measurement_mask,
        measurement_outcomes = result.schedule.outcomes,
        unitary_angles = result.schedule.unitary_angles,
        sample_free_energy = result.sample_free_energy,
        halfchain_entropy = result.halfchain_entropy,
        final_entropy_profile = result.final_entropy_profile,
        final_quarter_tmi = result.final_quarter_tmi,
        final_y_expectation = result.final_y_expectation,
    )
    return output_path
end

function save_hybrid_sector_ensemble(
    L::Integer,
    p::Real,
    periods::Integer,
    seed_start::Integer,
    seed_end::Integer;
    sector::Symbol = :trivial,
    γ::Real = 1.0,
)
    result = hybrid_sector_ensemble(
        L,
        p,
        periods,
        seed_start,
        seed_end;
        sector = sector,
        γ = γ,
    )
    output_dir = hybrid_output_dir(L, p; γ = γ)
    mkpath(output_dir)
    output_path = joinpath(
        output_dir,
        "entanglement_sector$(sector)_periods$(periods)_seeds$(seed_start)-$(seed_end).jld2",
    )
    jldsave(
        output_path;
        L = Int(L),
        p = Float64(p),
        γ = Float64(γ),
        periods = Int(periods),
        seed_start = Int(seed_start),
        seed_end = Int(seed_end),
        sector = String(sector),
        result...,
    )
    return output_path
end

function save_hybrid_learnability_ensemble(
    L::Integer,
    p::Real,
    periods::Integer,
    seed_start::Integer,
    seed_end::Integer;
    γ::Real = 1.0,
)
    result = hybrid_learnability_ensemble(
        L,
        p,
        periods,
        seed_start,
        seed_end;
        γ = γ,
    )
    output_dir = hybrid_output_dir(L, p; γ = γ)
    mkpath(output_dir)
    output_path = joinpath(
        output_dir,
        "learnability_periods$(periods)_seeds$(seed_start)-$(seed_end).jld2",
    )
    jldsave(
        output_path;
        L = Int(L),
        p = Float64(p),
        γ = Float64(γ),
        periods = Int(periods),
        seed_start = Int(seed_start),
        seed_end = Int(seed_end),
        result...,
    )
    return output_path
end

function save_hybrid_lyapunov_trajectory(
    L::Integer,
    p::Real,
    periods::Integer,
    n_states::Integer,
    seed::Integer;
    sector::Symbol = :trivial,
    γ::Real = 1.0,
)
    trajectory, spectrum = simulate_hybrid_lyapunov_exact(
        L,
        p,
        periods,
        seed;
        sector = sector,
        n_states = Int(n_states),
        γ = γ,
    )
    output_dir = hybrid_output_dir(L, p; γ = γ)
    mkpath(output_dir)
    output_path = joinpath(
        output_dir,
        "lyapunov_sector$(sector)_periods$(periods)_nstates$(n_states)_seed$(seed).jld2",
    )
    jldsave(
        output_path;
        L = Int(L),
        p = Float64(p),
        γ = Float64(γ),
        periods = Int(periods),
        n_states = Int(n_states),
        seed = Int(seed),
        sector = String(sector),
        measurement_mask = trajectory.schedule.measurement_mask,
        measurement_outcomes = trajectory.schedule.outcomes,
        unitary_angles = trajectory.schedule.unitary_angles,
        local_log_stretches = spectrum.local_log_stretches,
        lyapunov_exponents = spectrum.lyapunov_exponents,
        free_energy_spectrum = spectrum.free_energy_spectrum,
        sector_dimension = spectrum.sector_dimension,
        sector_leakage = spectrum.sector_leakage,
    )
    return output_path
end

function save_hybrid_purification_trajectory(
    L::Integer,
    p::Real,
    periods::Integer,
    seed::Integer;
    sector::Symbol = :trivial,
    γ::Real = 1.0,
)
    result = simulate_hybrid_purification_exact(
        L,
        p,
        periods,
        seed;
        sector = sector,
        γ = γ,
    )
    output_dir = hybrid_output_dir(L, p; γ = γ)
    mkpath(output_dir)
    output_path = joinpath(
        output_dir,
        "purification_sector$(sector)_periods$(periods)_seed$(seed).jld2",
    )
    jldsave(
        output_path;
        L = Int(L),
        p = Float64(p),
        γ = Float64(γ),
        periods = Int(periods),
        seed = Int(seed),
        sector = String(sector),
        sector_dimension = result.sector_dimension,
        measurement_mask = result.schedule.measurement_mask,
        measurement_outcomes = result.schedule.outcomes,
        unitary_angles = result.schedule.unitary_angles,
        entropy = result.entropy,
        purity = result.purity,
        y_expectation = result.y_expectation,
    )
    return output_path
end

function print_usage()
    println("Exclusive hybrid Fibonacci circuit (exact backend)")
    println()
    println("Sector trajectory:")
    println("  julia --project=. exm/Bulk_measure/hybrid_dynamics.jl sector L p periods seed [trivial|tau] [gamma]")
    println()
    println("Learnability ensemble:")
    println("  julia --project=. exm/Bulk_measure/hybrid_dynamics.jl learnability L p periods seed_start seed_end [gamma]")
    println()
    println("Fixed-sector entanglement ensemble:")
    println("  julia --project=. exm/Bulk_measure/hybrid_dynamics.jl entanglement L p periods seed_start seed_end [trivial|tau] [gamma]")
    println()
    println("Sector-resolved hybrid Lyapunov spectrum:")
    println("  julia --project=. exm/Bulk_measure/hybrid_dynamics.jl lyapunov L p periods n_states seed [trivial|tau] [gamma]")
    println()
    println("Mixed-state purification in one sector:")
    println("  julia --project=. exm/Bulk_measure/hybrid_dynamics.jl purification L p periods seed [trivial|tau] [gamma]")
end

function main(args = ARGS)
    isempty(args) && return print_usage()
    action = Symbol(args[1])
    if action == :sector
        length(args) in 5:7 || return print_usage()
        L = parse(Int, args[2])
        p = parse(Float64, args[3])
        periods = parse(Int, args[4])
        seed = parse(Int, args[5])
        sector = length(args) >= 6 ? Symbol(args[6]) : :trivial
        γ = length(args) >= 7 ? parse(Float64, args[7]) : 1.0
        output = save_hybrid_sector_trajectory(
            L,
            p,
            periods,
            seed;
            sector = sector,
            γ = γ,
        )
        println("saved: $output")
    elseif action == :learnability
        length(args) in 6:7 || return print_usage()
        L = parse(Int, args[2])
        p = parse(Float64, args[3])
        periods = parse(Int, args[4])
        seed_start = parse(Int, args[5])
        seed_end = parse(Int, args[6])
        γ = length(args) == 7 ? parse(Float64, args[7]) : 1.0
        output = save_hybrid_learnability_ensemble(
            L,
            p,
            periods,
            seed_start,
            seed_end;
            γ = γ,
        )
        println("saved: $output")
    elseif action == :entanglement
        length(args) in 6:8 || return print_usage()
        L = parse(Int, args[2])
        p = parse(Float64, args[3])
        periods = parse(Int, args[4])
        seed_start = parse(Int, args[5])
        seed_end = parse(Int, args[6])
        sector = length(args) >= 7 ? Symbol(args[7]) : :trivial
        γ = length(args) >= 8 ? parse(Float64, args[8]) : 1.0
        output = save_hybrid_sector_ensemble(
            L,
            p,
            periods,
            seed_start,
            seed_end;
            sector = sector,
            γ = γ,
        )
        println("saved: $output")
    elseif action == :lyapunov
        length(args) in 6:8 || return print_usage()
        L = parse(Int, args[2])
        p = parse(Float64, args[3])
        periods = parse(Int, args[4])
        n_states = parse(Int, args[5])
        seed = parse(Int, args[6])
        sector = length(args) >= 7 ? Symbol(args[7]) : :trivial
        γ = length(args) >= 8 ? parse(Float64, args[8]) : 1.0
        output = save_hybrid_lyapunov_trajectory(
            L,
            p,
            periods,
            n_states,
            seed;
            sector = sector,
            γ = γ,
        )
        println("saved: $output")
    elseif action == :purification
        length(args) in 5:7 || return print_usage()
        L = parse(Int, args[2])
        p = parse(Float64, args[3])
        periods = parse(Int, args[4])
        seed = parse(Int, args[5])
        sector = length(args) >= 6 ? Symbol(args[6]) : :trivial
        γ = length(args) >= 7 ? parse(Float64, args[7]) : 1.0
        output = save_hybrid_purification_trajectory(
            L,
            p,
            periods,
            seed;
            sector = sector,
            γ = γ,
        )
        println("saved: $output")
    else
        print_usage()
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
