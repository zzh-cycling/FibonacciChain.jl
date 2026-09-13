using FibonacciChain, ITensorMPS, ITensors, Test, Random, LinearAlgebra

# Reference the pre-optimization Born algorithm: reconstruct gates per event,
# contract the whole MPS for probabilities, and retain the old RNG/branch order.
function legacy_mps_born(model, sites, initial, config)
    ψ = deepcopy(initial)
    rng = copy(config.rng)
    periods = config.t₂ - config.t₁ + 1
    nlayers = FibonacciChain.layers_per_period(model)
    samples = falses(periods*nlayers, FibonacciChain._samples_per_layer(model))
    energies = zeros(Float64, periods*nlayers)
    entropies = Float64[]
    stride = config.truncate_every_events
    function branch(ψ, operator)
        result = stride == 1 ?
            apply(operator, ψ; cutoff=config.cutoff, mindim=config.mindim, maxdim=config.maxdim) :
            apply(operator, ψ; cutoff=0.0)
        probability = real(inner(result, result))
        normalize!(result)
        return result, probability
    end
    for period in 1:periods, layer in 1:nlayers
        row = (period-1)*nlayers + layer
        τ = row == periods*nlayers && config.enable_τ_eff ? config.τ/2 : config.τ
        positions, local_model, strength = FibonacciChain._obtain_measurement_config(model, row, τ)
        cols = FibonacciChain._get_sample_column_indices(model, row)
        for (k, site) in enumerate(positions)
            M0 = FibonacciChain._measurement_operator_mps_application(local_model, sites, site, strength, false)
            ψ0, p0 = branch(ψ, M0)
            if rand(rng) < p0
                ψ = ψ0
                energies[row] -= log(p0)
            else
                M1 = FibonacciChain._measurement_operator_mps_application(local_model, sites, site, strength, true)
                ψ, _ = branch(ψ, M1)
                samples[row, cols[k]] = true
                energies[row] -= log(1-p0)
            end
            if stride > 1 && (k % stride == 0 || k == length(positions))
                ψ = truncate(ψ; cutoff=config.cutoff, mindim=config.mindim, maxdim=config.maxdim)
                normalize!(ψ)
            end
        end
        layer == nlayers && push!(entropies, ee_mps(ψ, length(sites)÷2))
    end
    return (; state=ψ, samples, energies, entropies, rng)
end

@testset "MPS optimization preserves seeded Born trajectories" begin
    models = (
        AnyonModel(FibonacciAnyon(), 8; pbc=true),
        AnyonModel(SpinHalf(), 6; model_type=:Ising, pbc=true),
        AnyonModel(SpinHalf(), 6; model_type=:OBF, λ=0.3, pbc=true),
    )
    for model in models, stride in (1, 2), half_layer in (false, true)
        @testset "$(typeof(model)), stride=$stride, half_layer=$half_layer" begin
            ψ, sites = initial_mps(model.N)
            original = deepcopy(ψ)
            config = MeasureConfig(τ=0.8, t₂=3, mode=:Born, rng=MersenneTwister(36),
                cutoff=1e-12, mindim=1, maxdim=4, truncate_every_events=stride,
                enable_τ_eff=half_layer)
            expected = legacy_mps_born(model, sites, ψ, config)
            actual = bulk_evolution(model, sites, ψ, config)
            @test actual.samples == expected.samples
            @test any(actual.samples) && !all(actual.samples)
            @test actual.free_energys ≈ expected.energies atol=2e-6 rtol=1e-6
            @test actual.entanglement_entropys ≈ expected.entropies atol=2e-6 rtol=1e-6
            @test abs(inner(actual.state, expected.state))^2 ≈ 1 atol=1e-9
            @test norm(actual.state) ≈ 1 atol=1e-12
            @test maxlinkdim(actual.state) <= config.maxdim
            @test rand(config.rng) == rand(expected.rng)
            @test abs(inner(ψ, original))^2 ≈ 1 atol=1e-12

            replay = MeasureConfig(τ=config.τ, t₂=config.t₂, mode=:sample,
                cutoff=config.cutoff, maxdim=config.maxdim, truncate_every_events=stride,
                enable_τ_eff=half_layer)
            replayed = bulk_evolution(model, sites, ψ, replay, actual.samples)
            @test abs(inner(actual.state, replayed.state))^2 ≈ 1 atol=1e-9
            # Finite-bond compression makes the true-branch norm differ from
            # 1-p0. Compare free energies to the legacy Born path above, not to
            # the replay convention, which uses the selected branch's norm.
        end
    end
end

@testset "Canonical MPS probability agrees with full contraction" begin
    rng = MersenneTwister(36)
    sites = siteinds("Qubit", 6)
    tensor = ITensor(randn(rng, ComplexF64, ntuple(_ -> 2, 6)), sites...)
    ψ = MPS(tensor, sites; cutoff=0.0)
    normalize!(ψ)
    model = AnyonModel(FibonacciAnyon(), 6; pbc=true)
    for i in (1, 3, 6), sign in (false, true), deferred in (false, true)
        M = FibonacciChain._measurement_operator_mps_application(model, sites, i, 0.7, sign)
        expected = deferred ? apply(M, ψ; cutoff=0.0) : apply(M, ψ; cutoff=1e-12, maxdim=3)
        p_expected = real(inner(expected, expected))
        normalize!(expected)
        actual, p = FibonacciChain._measuremap_with_operator(ψ, M;
            cutoff=1e-12, maxdim=3, truncate_per_event=!deferred)
        @test p ≈ p_expected atol=1e-12 rtol=1e-12
        @test abs(inner(actual, expected))^2 ≈ 1 atol=1e-12
    end
end
