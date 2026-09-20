using FibonacciChain, LinearAlgebra, Random, Test

@testset "Ising chain with categorical KW reference" begin
    rng = MersenneTwister(420)
    X = [0.0 1.0; 1.0 0.0]
    Z = [1.0 0.0; 0.0 -1.0]
    id2 = Matrix{Float64}(I, 2, 2)
    for L in 2:4
        model = AnyonModel(SpinHalf(), L; model_type = :Ising, pbc = true)
        dim = 2^L
        D = kramers_wannier_operator(model)
        Y = [zeros(dim, dim) D'; D zeros(dim, dim)]
        pauli(P, j) = reduce(kron, [i == j ? P : id2 for i in 1:L])
        eta = reduce(kron, fill(X, L))
        @test Y' ≈ Y
        @test Y^2 ≈ kron(id2, Matrix{Float64}(I, dim, dim) + eta) atol = 1e-12
        state = normalize(randn(rng, ComplexF64, 2dim))
        @test categorical_kw_expectation(model, 3state) ≈ real(dot(state, Y * state)) atol = 1e-12
        chain = normalize(randn(rng, ComplexF64, dim))
        @test ising_reference_state(chain) ≈ kron([1, 1] / sqrt(2), chain)
        @test categorical_kw_expectation(model, [chain; zeros(dim)]) ≈ 0 atol = 1e-12

        # Independent dense Pauli calculation checks all physical gates,
        # the shifted B site, shared outcomes, and joint normalization.
        basis = anyon_basis(model)
        mx = AnyonModel(SpinHalf(), L; model_type = :Ising, pbc = true, measure_operator = :X)
        mz = AnyonModel(SpinHalf(), L; model_type = :Ising, pbc = true, measure_operator = :ZZ)
        for layer in 1:2, site in 1:L
            next = mod1(site + 1, L)
            A = isodd(layer) ? pauli(X, site) : pauli(Z, site) * pauli(Z, next)
            B = isodd(layer) ? pauli(Z, site) * pauli(Z, next) : pauli(X, next)
            O = [A zeros(dim, dim); zeros(dim, dim) B]
            @test Y * O ≈ O * Y atol = 1e-12
            completeness = zeros(ComplexF64, 2dim, 2dim)
            for bit in (false, true)
                M = exp((bit ? -1 : 1) * 0.7 / 2 * O) / sqrt(2cosh(0.7))
                buffer = similar(state)
                FibonacciChain._ising_reference_measure!(buffer, state, basis, mx, mz, 0.7, layer, site, bit)
                @test buffer ≈ M * state atol = 1e-12
                completeness += M' * M
            end
            @test completeness ≈ I atol = 1e-12
        end

        for half_boundary in (false, true)
            config = MeasureConfig(τ = 0.7, t₂ = 3, mode = :Born,
                rng = MersenneTwister(17), enable_τ_eff = half_boundary)
            result = ising_reference_evolution(model, state, config)
            replay_config = MeasureConfig(τ = 0.7, t₂ = 3, mode = :sample,
                enable_τ_eff = half_boundary)
            replay = ising_reference_evolution(model, state, replay_config, result.samples)
            @test result.state ≈ replay.state atol = 1e-12
            @test result.kw_expectation_values ≈ replay.kw_expectation_values atol = 1e-12
            dense = copy(state)
            dense_kw = [real(dot(dense, Y * dense))]
            # Reproduce the RNG draws with dense joint Born probabilities.
            dense_rng = MersenneTwister(17)
            for layer in 1:6
                tau = layer == 6 && half_boundary ? 0.35 : 0.7
                for site in 1:L
                    next = mod1(site + 1, L)
                    A = isodd(layer) ? pauli(X, site) : pauli(Z, site) * pauli(Z, next)
                    B = isodd(layer) ? pauli(Z, site) * pauli(Z, next) : pauli(X, next)
                    O = [A zeros(dim, dim); zeros(dim, dim) B]
                    M0 = exp(tau / 2 * O) / sqrt(2cosh(tau))
                    bit = rand(dense_rng) >= norm(M0 * dense)^2
                    @test bit == result.samples[layer, site]
                    M = bit ? exp(-tau / 2 * O) / sqrt(2cosh(tau)) : M0
                    dense = normalize(M * dense)
                end
                iseven(layer) && push!(dense_kw, real(dot(dense, Y * dense)))
            end
            @test result.state ≈ dense atol = 1e-11
            @test result.kw_expectation_values ≈ dense_kw atol = 1e-11
        end

        # Exact charge eigenstates stay sharp; an A-only state has zero mean.
        even = normalize(chain + eta * chain)
        for sign in (-1, 1)
            eigenstate = [even / sqrt(2); sign * D * even / 2]
            result = ising_reference_evolution(model, eigenstate,
                MeasureConfig(τ = 0.9, t₂ = 4, mode = :Born, rng = MersenneTwister(8)))
            @test result.kw_expectation_values ≈ fill(sign * sqrt(2), 5) atol = 1e-11
        end
        result = ising_reference_evolution(model, [chain; zeros(dim)],
            MeasureConfig(τ = 1000.0, t₂ = 2, mode = :Born, rng = MersenneTwister(3)))
        @test all(iszero, result.kw_expectation_values)
        @test norm(result.state) ≈ 1 atol = 1e-12
        @test_throws DimensionMismatch categorical_kw_expectation(model, chain)
        @test_throws ArgumentError categorical_kw_expectation(model, zeros(2dim))
        @test_throws DimensionMismatch ising_reference_evolution(model, state,
            MeasureConfig(τ = 0.7, t₂ = 2, mode = :sample), falses(1, L))
    end
end
