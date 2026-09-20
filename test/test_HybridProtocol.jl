using FibonacciChain, LinearAlgebra, Test
using JLD2
include("../exm/HybridEvolution/protocol.jl")

@testset "Hybrid JLD2 output round trip" begin
    results = [samples_generate(6, 0.0, 1, seed; save_schedule = true) for seed in 1:2]
    mktempdir() do directory
        save_ensemble(directory, results)
        data = load(joinpath(directory, "trajectories.jld2"))
        @test size(data["S_half"]) == (2, 2)
        @test data["trajectory_seed"] == [1, 2]
        @test data["time"] == [0, 1]
        @test data["S_half"][2, :] == results[2].entropy
        @test data["Y_expectation"][1, :] == results[1].y_expectation
        summary = load(joinpath(directory, "summary.jld2"))
        @test summary["S_mean"] ≈ vec(mean(data["S_half"]; dims = 1))
        @test load(joinpath(directory, "sharpening.jld2"), "censored")
        @test isnan(load(joinpath(directory, "sharpening.jld2"), "t_sharp"))
        record = load(joinpath(directory, "schedule_seed1.jld2"))
        schedule = HybridGateSchedule(record["measurement_mask"], record["outcomes"],
            record["unitary_angles"])
        @test isequal(schedule.unitary_angles, results[1].schedule.unitary_angles)
        model = AnyonModel(FibonacciAnyon(), 6; pbc = true)
        state = zeros(length(anyon_basis(model)))
        state[1] = 1
        replay = bulk_evolution(model, state,
            HybridConfig(τ = Inf, t₂ = 1, mode = :sample, enable_τ_eff = false), schedule)
        @test last(replay.entanglement_entropys) == last(results[1].entropy)
        @test all(endswith(".jld2"), readdir(directory))
    end
    mktempdir() do directory
        mps = samples_generate_mps(6, 0.0, 1, 1)
        save_ensemble(directory, [mps])
        @test load(joinpath(directory, "mps_diagnostics.jld2"), "final_bond_dimension") ==
            [mps.final_bond_dimension]
        @test all(isnan, load(joinpath(directory, "summary.jld2"), "S_sem"))
        @test !isfile(joinpath(directory, "schedule_seed1.jld2"))
    end
end

@testset "Coherent hybrid protocol" begin
    model = AnyonModel(FibonacciAnyon(), 6; pbc = true)
    Y = topological_charge_operator(model)
    initial = zeros(ComplexF64, size(Y, 1))
    initial[1] = 1
    @test Y * Y ≈ I + Y
    for site in 1:6
        P = FibonacciChain.measure_matrix(model, Inf, site, true)
        Q = FibonacciChain.measure_matrix(model, Inf, site, false)
        U = P + cis(0.37) * Q
        @test P + Q ≈ I
        @test P' * P + Q' * Q ≈ I
        @test U' * U ≈ I
        @test norm(Y * P - P * Y) < 1e-12
        @test norm(Y * U - U * Y) < 1e-12
    end
    for p in (0.0, 0.5, 1.0)
        r = samples_generate(6, p, 3, 17; save_schedule = true)
        @test 0 < r.initial_weight < 1
        @test r.times == 0:3
        config_born = HybridConfig(τ = Inf, t₂ = 3, mode = :Born,
            rng = MersenneTwister(17), enable_τ_eff = false, p = p,
            random_angles = true, track_y_expectation = true)
        direct = bulk_evolution(model, initial, config_born)
        @test direct.schedule.measurement_mask == r.schedule.measurement_mask
        @test direct.schedule.outcomes == r.schedule.outcomes
        @test isequal(direct.schedule.unitary_angles, r.schedule.unitary_angles)
        @test Float64.(direct.y_expectation_values) == r.y_expectation[2:end]
        @test Float64.(direct.entanglement_entropys) == r.entropy[2:end]
        current = copy(initial)
        for t in 1:3
            rows = (2t - 1):2t
            schedule = HybridGateSchedule(r.schedule.measurement_mask[rows, :],
                r.schedule.outcomes[rows, :], r.schedule.unitary_angles[rows, :])
            config = HybridConfig(τ = Inf, t₂ = 1, mode = :sample, enable_τ_eff = false)
            evolved = bulk_evolution(model, current, config, schedule)
            @test isempty(evolved.y_expectation_values)
            current = evolved.state
            @test r.y_expectation[t + 1] ≈ real(dot(current, Y * current)) atol = 1e-6
            @test r.entropy[t + 1] ≈ evolved.entanglement_entropys[1] atol = 1e-6
        end
        p == 0 && @test all(y -> isapprox(y, r.y_expectation[1]; atol = 1e-6), r.y_expectation)
        repeated = samples_generate(6, p, 3, 17)
        @test repeated.entropy == r.entropy
        @test repeated.y_expectation == r.y_expectation
    end
    sparse = samples_generate(6, 0.5, 3, 17; stride = 2)
    @test sparse.times == [0, 2, 3]
    @test_throws ArgumentError samples_generate(5, 0.5, 2, 1)
    legacy = HybridMeasurementOutcome(initial,
        HybridGateSchedule(falses(2, 3), falses(2, 3)), zeros(Float32, 2), zeros(Float32, 1))
    @test isempty(legacy.y_expectation_values)
end

@testset "MPS hybrid protocol" begin
    for p in (0.0, 0.5, 1.0)
        exact = samples_generate(6, p, 2, 17; save_schedule = true)
        mps = samples_generate_mps(6, p, 2, 17;
            cutoff = 1e-14, maxdim = 64, save_schedule = true)
        @test mps.schedule.measurement_mask == exact.schedule.measurement_mask
        @test mps.schedule.outcomes == exact.schedule.outcomes
        @test isequal(mps.schedule.unitary_angles, exact.schedule.unitary_angles)
        @test mps.entropy ≈ exact.entropy atol = 1e-6
        @test mps.y_expectation ≈ exact.y_expectation atol = 1e-6
        @test mps.initial_weight ≈ exact.initial_weight atol = 1e-12
        @test mps.final_bond_dimension <= 64
    end
    model = AnyonModel(FibonacciAnyon(), 6; pbc = true)
    record = samples_generate(6, 0.5, 2, 19; save_schedule = true)
    for interval in (1, 2), enforce in (false, true)
        state, sites = initial_mps(6)
        config = HybridConfig(τ = Inf, t₂ = 2, mode = :sample,
            enable_τ_eff = false, track_y_expectation = true,
            cutoff = 1e-14, maxdim = 64, truncate_every_events = interval,
            enforce_fibonacci_constraint = enforce)
        replay = bulk_evolution(model, sites, state, config, record.schedule)
        @test replay.entanglement_entropys ≈ record.entropy[2:end] atol = 1e-6
        @test replay.y_expectation_values ≈ record.y_expectation[2:end] atol = 1e-6
        Y = topological_charge_mpo(sites)
        @test last(replay.y_expectation_values) ≈
            real(inner(prime(replay.state), Y, replay.state)) atol = 1e-6
        @test real(inner(replay.state, replay.state)) ≈ 1 atol = 1e-10
    end
    sparse = samples_generate_mps(6, 0.5, 3, 17; stride = 2)
    @test sparse.times == [0, 2, 3]
    @test_throws ArgumentError samples_generate_mps(6, 0.5, 1, 1; maxdim = 0)
end
