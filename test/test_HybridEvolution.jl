using FibonacciChain
using Random
using LinearAlgebra
using Test
using ITensors
using ITensorMPS

@testset "Hybrid monitored Fibonacci evolution" begin
    N = 6
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    state = zeros(length(anyon_basis(model)))
    state[1] = 1.0

    @testset "validation and exact Born replay" begin
        @test_throws ArgumentError HybridConfig(τ = 1.0, t₂ = 1, p = -0.1)
        @test_throws ArgumentError HybridConfig(τ = 1.0, t₂ = 1, p = 1.1)

        born_config = HybridConfig(
            τ = 1.0,
            t₂ = 2,
            mode = :Born,
            rng = MersenneTwister(7),
            enable_τ_eff = false,
            p = 0.5,
            θ = 0.7,
        )
        born = bulk_evolution(model, state, born_config)
        @test size(born.schedule.measurement_mask) == (4, N ÷ 2)
        @test any(born.schedule.measurement_mask)
        @test any(.!born.schedule.measurement_mask)
        @test all(isnan, born.schedule.unitary_angles[born.schedule.measurement_mask])
        @test all(isfinite, born.schedule.unitary_angles[.!born.schedule.measurement_mask])
        @test norm(born.state) ≈ 1 atol = 1e-12

        replay_config = HybridConfig(
            τ = 1.0,
            t₂ = 2,
            mode = :sample,
            enable_τ_eff = false,
            p = 0.5,
            θ = 0.7,
        )
        replay = bulk_evolution(model, state, replay_config, born.schedule)
        @test replay.state ≈ born.state atol = 1e-13
        @test replay.free_energys ≈ born.free_energys atol = 1e-7
    end

    @testset "independent random angles and replay" begin
        config = HybridConfig(
            τ = Inf,
            t₂ = 2,
            mode = :Born,
            rng = MersenneTwister(19),
            enable_τ_eff = false,
            p = 0.0,
            random_angles = true,
        )
        generated = bulk_evolution(model, state, config)
        @test !any(generated.schedule.measurement_mask)
        @test all(isfinite, generated.schedule.unitary_angles)
        @test all(x -> 0 <= x < 2π, generated.schedule.unitary_angles)
        @test length(unique(generated.schedule.unitary_angles)) > 1

        replay_config = HybridConfig(
            τ = Inf,
            t₂ = 2,
            mode = :sample,
            enable_τ_eff = false,
            p = 0.0,
            random_angles = true,
        )
        replay = bulk_evolution(model, state, replay_config, generated.schedule)
        @test replay.state ≈ generated.state atol = 1e-13
    end

    @testset "limits p=1, p=0, and θ=0" begin
        seed = 42
        measurement_config = MeasureConfig(
            τ = 1.0,
            t₂ = 1,
            mode = :Born,
            rng = MersenneTwister(seed),
            enable_τ_eff = false,
        )
        legacy = bulk_evolution(model, state, measurement_config)
        hybrid = bulk_evolution(
            model,
            state,
            HybridConfig(
                MeasureConfig(
                    τ = 1.0,
                    t₂ = 1,
                    mode = :Born,
                    rng = MersenneTwister(seed),
                    enable_τ_eff = false,
                );
                p = 1.0,
                θ = 0.9,
            ),
        )
        @test hybrid.schedule.outcomes == legacy.samples
        @test all(hybrid.schedule.measurement_mask)
        @test hybrid.state ≈ legacy.state atol = 1e-13
        @test hybrid.free_energys ≈ legacy.free_energys atol = 1e-7

        no_gates = HybridGateSchedule(falses(2, N ÷ 2), falses(2, N ÷ 2))
        identity_run = bulk_evolution(
            model,
            state,
            HybridConfig(
                τ = 1.0,
                t₂ = 1,
                mode = :sample,
                enable_τ_eff = false,
                p = 0.0,
                θ = 0.0,
            ),
            no_gates,
        )
        @test identity_run.state ≈ state atol = 1e-14
        @test iszero(identity_run.free_energys)
    end

    @testset "exact/MPS schedule cross-check" begin
        schedule = HybridGateSchedule(
            BitMatrix([0 1 0; 1 0 1]),
            BitMatrix([0 0 0; 1 0 0]),
            [0.2 NaN 1.1; NaN 2.3 NaN],
        )
        config = HybridConfig(
            τ = 1.0,
            t₂ = 1,
            mode = :sample,
            enable_τ_eff = false,
            p = 0.5,
            θ = 0.7,
            cutoff = 1e-14,
            maxdim = 100,
        )
        exact = bulk_evolution(model, state, config, schedule)
        ψ, sites = initial_mps(N)
        mps = bulk_evolution(model, sites, ψ, config, schedule)

        dense = vec(Array(reduce(*, mps.state).tensor))
        reverse_N_bits(x) = sum(((x >> (j - 1)) & 1) << (N - j) for j = 1:N)
        legal_indices = [reverse_N_bits(Int(b.buf)) + 1 for b in anyon_basis(model)]
        dense_legal = dense[legal_indices]
        overlap = dot(exact.state, dense_legal)
        dense_legal .*= conj(overlap) / abs(overlap)

        @test dense_legal ≈ exact.state atol = 1e-10
        @test mps.free_energys ≈ exact.free_energys atol = 1e-6
        @test real(inner(mps.state, mps.state)) ≈ 1 atol = 1e-11
    end


    @testset "Bayesian common-record diagnostics" begin
        Y = topological_charge_operator(model)
        ϕ = (1 + √5) / 2
        eig = eigen(Symmetric(Y))
        state_y1 = eig.vectors[:, findfirst(x -> isapprox(x, ϕ), eig.values)]
        state_ytau = eig.vectors[:, findfirst(x -> isapprox(x, -inv(ϕ)), eig.values)]

        # With no measurements the record contains no sector information.
        unitary_only = hybrid_bayesian_evolution(
            model,
            state_y1,
            state_ytau,
            HybridConfig(
                τ = Inf,
                t₂ = 3,
                mode = :Born,
                rng = MersenneTwister(31),
                enable_τ_eff = false,
                p = 0.0,
                random_angles = true,
            );
            sampled_sector = :mixture,
            sector_rng = MersenneTwister(32),
        )
        @test iszero(unitary_only.log_likelihood_ratio)
        @test unitary_only.posterior_y1 == fill(0.5, 3)
        @test unitary_only.record_fidelity_estimator == ones(3)
        @test unitary_only.conditional_entropy_estimator ≈ fill(log(2), 3)
        @test iszero(unitary_only.mutual_information_estimator)
        @test unitary_only.bayes_error_estimator == fill(0.5, 3)
        @test unitary_only.y_expectation_y1 ≈ fill(ϕ, 3) atol = 1e-10
        @test unitary_only.y_expectation_ytau ≈ fill(-inv(ϕ), 3) atol = 1e-10

        replay_config = HybridConfig(
            τ = Inf,
            t₂ = 3,
            mode = :sample,
            enable_τ_eff = false,
            p = 0.0,
            random_angles = true,
        )
        replay = hybrid_bayesian_evolution(
            model,
            state_y1,
            state_ytau,
            replay_config,
            unitary_only.schedule,
        )
        @test replay.state_y1 ≈ unitary_only.state_y1 atol = 1e-13
        @test replay.state_ytau ≈ unitary_only.state_ytau atol = 1e-13
        @test replay.log_likelihood_ratio == unitary_only.log_likelihood_ratio

        # Enumerate both outcomes of a schedule containing exactly one
        # measurement. This directly verifies
        # E_Pmix[sech(ℓ/2)] = Σ_m sqrt(P₁(m)Pτ(m)).
        mask = falses(2, N ÷ 2)
        mask[1, 1] = true
        angles = zeros(2, N ÷ 2)
        angles[mask] .= NaN
        enum_config = HybridConfig(
            τ = Inf,
            t₂ = 1,
            mode = :sample,
            enable_τ_eff = false,
            p = 1.0,
        )
        weighted_estimator = 0.0
        bhattacharyya = 0.0
        total_mixture_probability = 0.0
        for observed_outcome in (false, true)
            outcomes = falses(2, N ÷ 2)
            outcomes[1, 1] = observed_outcome
            result = hybrid_bayesian_evolution(
                model,
                state_y1,
                state_ytau,
                enum_config,
                HybridGateSchedule(mask, outcomes, angles),
            )
            p_y1 = exp(last(result.log_record_probability_y1))
            p_ytau = exp(last(result.log_record_probability_ytau))
            p_mix = (p_y1 + p_ytau) / 2
            total_mixture_probability += p_mix
            weighted_estimator += p_mix * last(result.record_fidelity_estimator)
            bhattacharyya += sqrt(p_y1 * p_ytau)
        end
        @test total_mixture_probability ≈ 1 atol = 1e-12
        @test weighted_estimator ≈ bhattacharyya atol = 1e-12
    end


    @testset "hybrid Lyapunov spectrum" begin
        periods = 4
        unitary_config = HybridConfig(
            τ = Inf,
            t₂ = periods,
            mode = :Born,
            rng = MersenneTwister(51),
            enable_τ_eff = false,
            p = 0.0,
            random_angles = true,
        )
        generated = bulk_evolution(model, state, unitary_config)
        unitary_spectrum = hybrid_lyapunov_spectrum(
            model,
            unitary_config,
            generated.schedule;
            n_states = 3,
            sector = :trivial,
        )
        @test size(unitary_spectrum.lyapunov_exponents) == (3, periods)
        @test maximum(abs, unitary_spectrum.lyapunov_exponents) < 1e-12
        @test maximum(unitary_spectrum.sector_leakage) < 1e-10

        # At p=1 the hybrid routine reduces exactly to the existing
        # measurement-only sector-resolved Lyapunov calculation.
        τ = 0.7
        samples = BitMatrix([
            0 1 0
            1 0 1
            1 1 0
            0 0 1
            1 0 0
            0 1 1
            0 0 0
            1 1 1
        ])
        mask = trues(size(samples))
        measurement_schedule = HybridGateSchedule(mask, samples)
        measurement_config = HybridConfig(
            τ = τ,
            t₂ = periods,
            mode = :sample,
            enable_τ_eff = false,
            p = 1.0,
        )
        hybrid = hybrid_lyapunov_spectrum(
            model,
            measurement_config,
            measurement_schedule;
            n_states = 3,
            sector = :trivial,
        )
        measurement_only = lyapunov_spectrum_topological_sector(
            model,
            τ,
            samples;
            n_states = 3,
            sector = :trivial,
        )
        @test hybrid.local_log_stretches ≈
              measurement_only.local_log_stretches atol = 1e-12
        @test hybrid.lyapunov_exponents ≈
              measurement_only.lyapunov_exponents atol = 1e-12
    end
end
