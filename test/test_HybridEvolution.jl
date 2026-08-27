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
end
