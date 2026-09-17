using FibonacciChain, ITensorMPS, ITensors, Test, Random, LinearAlgebra

@testset "Spin boundary MPOs equal the dense Pauli measurement" begin
    for N in (3, 4, 6), kind in (:Ising, :OBF)
        sites = siteinds("Qubit", N)
        operators = kind == :Ising ? (:ZZ,) : (:ZZ, :XZZ, :ZZX)
        for operator in operators, i in (N-1, N), sign in (false, true), τ in (0.0, 0.8, Inf)
            model = AnyonModel(SpinHalf(), N; model_type=kind, pbc=true, measure_operator=operator)
            gate = measurement_operator_mps(model, sites, i, τ, sign)
            optimized = FibonacciChain._measurement_operator_mps_application(model, sites, i, τ, sign)
            actual = optimized isa MPO ? prod(optimized) : optimized
            expected = gate
            for site in sites
                hasind(expected, site) || (expected *= op("I", site))
                hasind(actual, site) || (actual *= op("I", site))
            end
            @test norm(actual - expected) <= 1e-13 * norm(expected)
            if optimized isa MPO
                @test maxlinkdim(optimized) == 2
            end
        end
    end
end

@testset "Compact Fibonacci boundary MPO is the exact local operator" begin
    for N in (3, 4, 6), convention in (:Antiferro, :Ferro)
        model = AnyonModel(FibonacciAnyon(), N; pbc=true, measure_operator=convention)
        sites = siteinds("Qubit", N)
        for i in (1, N), sign in (false, true), τ in (0.0, 1e-8, 0.8, Inf)
            mpo = FibonacciChain.measurement_operator_mpo(model, sites, i, τ, sign)
            expected = measurement_operator_mps(model, sites, i, τ, sign)
            for j in setdiff(1:N, (mod1(i-1, N), i, mod1(i+1, N)))
                expected *= op("I", sites[j])
            end
            # Compare the full operator, including states outside the fusion constraint.
            @test norm(prod(mpo) - expected) <= 1e-13 * norm(expected)
            @test linkdims(mpo) == (i == 1 ? [3; fill(2, N-2)] : [fill(2, N-2); 3])
        end
    end
end

# Reference the pre-optimization Born algorithm: reconstruct gates per event,
# contract the whole MPS for probabilities, and retain the old RNG/branch order.
function legacy_measurement_operator(model, sites, i, τ, sign)
    if model isa Union{AnyonModel{SpinHalf,:Ising},AnyonModel{SpinHalf,:OBF}}
        return measurement_operator_mps(model, sites, i, τ, sign)
    end
    if !(model isa AnyonModel{FibonacciAnyon} && model.pbc && i in (1, length(sites)))
        return FibonacciChain._measurement_operator_mps_application(model, sites, i, τ, sign)
    end
    ϕ = (1 + √5)/2
    cst = τ >= 100 ? 0.5 : (exp(τ)+1)/(2sqrt(exp(2τ)+1))
    coef = τ >= 100 ? 0.5 : (exp(τ)-1)/(2sqrt(exp(2τ)+1))
    coef *= sign ? -1 : 1
    coef *= model.measure_operator == :Antiferro ? 1 : -1
    im1, ip1 = mod1(i-1, length(sites)), mod1(i+1, length(sites))
    os = OpSum()
    os += cst, "I", 1
    os += coef, "Proj0", im1, "Z", i, "Proj1", ip1
    os += coef, "Proj1", im1, "Z", i, "Proj0", ip1
    os += -coef, "Proj1", im1, "Z", i, "Proj1", ip1
    os += coef*(1-2/ϕ), "Proj0", im1, "Z", i, "Proj0", ip1
    os += coef*(-2*ϕ^(-3/2)), "Proj0", im1, "X", i, "Proj0", ip1
    return MPO(os, sites)
end

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
            M0 = legacy_measurement_operator(local_model, sites, site, strength, false)
            ψ0, p0 = branch(ψ, M0)
            if rand(rng) < p0
                ψ = ψ0
                energies[row] -= log(p0)
            else
                M1 = legacy_measurement_operator(local_model, sites, site, strength, true)
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
            # Spin boundary MPOs avoid the old intermediate SWAP truncations.
            # Compare their trajectories in the exact six-site limit; compressed
            # spin trajectories are checked against dense evolution separately.
            maxdim = model isa AnyonModel{FibonacciAnyon} ? 4 : 8
            config = MeasureConfig(τ=0.8, t₂=3, mode=:Born, rng=MersenneTwister(36),
                cutoff=1e-12, mindim=1, maxdim=maxdim, truncate_every_events=stride,
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

function reference_spin_replay(model, sites, initial, samples; maxdim=nothing)
    state = isnothing(maxdim) ? prod(initial) : deepcopy(initial)
    for row in axes(samples, 1)
        τ = row == size(samples, 1) ? 0.4 : 0.8
        positions, local_model, strength = FibonacciChain._obtain_measurement_config(model, row, τ)
        cols = FibonacciChain._get_sample_column_indices(model, row)
        for (k, i) in enumerate(positions)
            gate = measurement_operator_mps(local_model, sites, i, strength, samples[row, cols[k]])
            state = isnothing(maxdim) ? apply(gate, state) :
                apply(gate, state; cutoff=1e-12, maxdim=maxdim)
            state /= norm(state)
        end
    end
    return state
end

@testset "Compressed spin trajectories converge to dense evolution" begin
    for kind in (:Ising, :OBF)
        model = AnyonModel(SpinHalf(), 6; model_type=kind, pbc=true, λ=0.3)
        initial, sites = initial_mps(6)
        config = MeasureConfig(τ=0.8, t₂=3, mode=:Born, rng=MersenneTwister(36),
            maxdim=8, cutoff=1e-12)
        samples = bulk_evolution(model, sites, initial, config).samples
        exact = reference_spin_replay(model, sites, initial, samples)
        errors = Float64[]
        for chi in (2, 4, 8)
            replay = MeasureConfig(τ=0.8, t₂=3, mode=:sample, maxdim=chi, cutoff=1e-12)
            actual = bulk_evolution(model, sites, initial, replay, samples).state
            legacy = reference_spin_replay(model, sites, initial, samples; maxdim=chi)
            error = 1 - abs(inner(exact, prod(actual)))^2
            old_error = 1 - abs(inner(exact, prod(legacy)))^2
            push!(errors, error)
            @test maxlinkdim(actual) <= chi
            @test norm(actual) ≈ 1 atol=1e-12
            # For these fixed trajectories, avoiding SWAP truncations improves
            # accuracy. This is a regression example, not a universal error bound.
            @test error <= old_error + 1e-10
        end
        @test errors[2] < errors[1]
        @test errors[3] < 1e-9
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
        old_M = legacy_measurement_operator(model, sites, i, 0.7, sign)
        expected = deferred ? apply(old_M, ψ; cutoff=0.0) : apply(old_M, ψ; cutoff=1e-12, maxdim=3)
        p_expected = real(inner(expected, expected))
        normalize!(expected)
        actual, p = FibonacciChain._measuremap_with_operator(ψ, M;
            cutoff=1e-12, maxdim=3, truncate_per_event=!deferred)
        @test p ≈ p_expected atol=1e-12 rtol=1e-12
        @test abs(inner(actual, expected))^2 ≈ 1 atol=1e-12
    end
end
