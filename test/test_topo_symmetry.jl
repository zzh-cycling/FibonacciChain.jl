using FibonacciChain
using LinearAlgebra
using Random
using Test

function Fibonacci_number(L)
    return (1/sqrt(5)) * (((1+sqrt(5))/2)^L - ((1-sqrt(5))/2)^L)
end

function Lucas_number(L)
    return Fibonacci_number(L+1) + Fibonacci_number(L-1)
end

@testset "Topological symmetry of Fibonacci-chain dynamics" begin
    L = 6
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    Y = topological_charge_operator(model)
    projector_strength = 1000.0
    weak_strength = atanh(0.8)
    atol = 1e-10

    @testset "Local projectors commute with Y" begin
        for site in 1:L
            Πp = FibonacciChain.measure_matrix(model, projector_strength, site, false)
            Πm = FibonacciChain.measure_matrix(model, projector_strength, site, true)

            @test Πp * Y ≈ Y * Πp atol = atol
            @test Πm * Y ≈ Y * Πm atol = atol
        end
    end

    @testset "Hamiltonian commutes with Y" begin
        H = anyon_ham(model)
        @test H * Y ≈ Y * H atol = atol
    end

    @testset "Weak measurement operators commute with Y" begin
        for site in 1:L
            Mplus = FibonacciChain.measure_matrix(model, weak_strength, site, false)
            Mminus = FibonacciChain.measure_matrix(model, weak_strength, site, true)

            @test Mplus * Y ≈ Y * Mplus atol = atol
            @test Mminus * Y ≈ Y * Mminus atol = atol
        end
    end

    @testset "post selection fusion-tree transfer matrix commutes with Y" begin
        n_periods = 4
        samples = BitMatrix(ones(2n_periods, L ÷ 2))

        dim = length(anyon_basis(model))
        dynamics_transfer_matrix = Matrix{Float64}(I, dim, dim)

        for period in 1:n_periods
            rows = (2period-1):(2period)
            period_transfer_matrix =
                transfer_matrix(model, weak_strength, BitMatrix(samples[rows, :]))
            dynamics_transfer_matrix =
                period_transfer_matrix * dynamics_transfer_matrix
        end

        @test dynamics_transfer_matrix * Y ≈
              Y * dynamics_transfer_matrix atol = atol
    end

    @testset "Random-outcome fusion-tree transfer matrix commutes with Y" begin
        rng = MersenneTwister(1234)
        n_periods = 4
        samples = BitMatrix(rand(rng, Bool, 2n_periods, L ÷ 2))

        dim = length(anyon_basis(model))
        dynamics_transfer_matrix = Matrix{Float64}(I, dim, dim)

        for period in 1:n_periods
            rows = (2period-1):(2period)
            period_transfer_matrix =
                transfer_matrix(model, weak_strength, BitMatrix(samples[rows, :]))
            dynamics_transfer_matrix =
                period_transfer_matrix * dynamics_transfer_matrix
        end

        @test dynamics_transfer_matrix * Y ≈
              Y * dynamics_transfer_matrix atol = atol
    end
end


@testset "Temperley Lieb algebra" begin
    N = 8
    τ = 1000.0
    ϕ = (1 + √5) / 2
    model = AnyonModel(FibonacciAnyon(), N)
    Xlis = ϕ .* [FibonacciChain.measure_matrix(model, τ, idx, true) for idx = 1:N] # s=1

    # X_i ^2 = d X_i
    @test all(Xlis[i] * Xlis[i] ≈ ϕ .* Xlis[i] for i = 1:N)
    # X_i * X_{i+1} * X_i = X_i
    @test all(Xlis[i] * Xlis[i+1] * Xlis[i] ≈ Xlis[i] for i = 1:(N-1))
    # X_i * X_{i-1} * X_i = X_i
    @test all(Xlis[i] * Xlis[i-1] * Xlis[i] ≈ Xlis[i] for i = 2:N)
    # [X_i, X_{j}] = 0, |i-j|>=2
    @test all(Xlis[i] * Xlis[j] ≈ Xlis[j] * Xlis[i] for i = 1:N for j = (i+2):(N-1))
end

@testset "Y and projector trace" begin
    L = 12  # finite-size deviations from the trace formulas are ~4e-5 at L=8, ~1e-7 at L=12
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    Y = topological_charge_operator(model)
    τ = 1000.0
    ϕ = (1 + √5) / 2

    H = anyon_ham(model)
    @test norm(H*Y - Y*H) ≈ 0.0 atol=1e-10 # Check that the Hamiltonian commutes with the topological charge operator

    Πplis = [FibonacciChain.measure_matrix(model, τ, idx, true) for idx = 1:1] # s=1
    Πmlis = [FibonacciChain.measure_matrix(model, τ, idx, false) for idx = 1:1] # s=τ
    Πp = Πplis[1]
    Πm = Πmlis[1]
    
    @test norm(Πp*Y - Y*Πp) ≈ 0.0 atol=1e-10 # local measurement operator commutes with the topological charge operator
    @test norm(Πm*Y - Y*Πm) ≈ 0.0 atol=1e-10
    
    energy, states = eigen((Y+Y')/2)
    y1_index = findall(x -> isapprox(x, ϕ), energy)
    yτ_index = findall(x -> isapprox(x, 1-ϕ), energy)
    
    projector_y1 = sum([states[:, i] * states[:, i]' for i in y1_index])
    projector_yτ = sum([states[:, i] * states[:, i]' for i in yτ_index])

    y = (-1/ϕ)^L
    
    @test tr(Πp) ≈ Lucas_number(L-2)
    @test tr(Πm) ≈ Lucas_number(L-1)
    @test tr(Πp*Y) ≈ y*ϕ^2
    @test tr(Πm*Y) ≈ -y*ϕ

    theoretical_Πm = (y+1/ϕ)/(ϕ+1/ϕ)^2*(Lucas_number(L-1)/ϕ -y*ϕ)/Fibonacci_number(L-1) + 
    (ϕ-y)/(ϕ+1/ϕ)^2*(Lucas_number(L-1)*ϕ + y*ϕ)/Fibonacci_number(L+1)
    theoretical_Πp = (y+1/ϕ)/(ϕ+1/ϕ)^2*(Lucas_number(L-2)/ϕ -y*ϕ)/Fibonacci_number(L-1) + 
    (ϕ-y)/(ϕ+1/ϕ)^2*(Lucas_number(L-2)*ϕ + y*ϕ)/Fibonacci_number(L+1)
    @test theoretical_Πm ≈ tr(Πm*projector_y1)/ length(y1_index) *(y+1/ϕ)/(ϕ+1/ϕ) + tr(Πm*projector_yτ)/ length(yτ_index) *(ϕ - y)/(ϕ+1/ϕ)
    @test theoretical_Πp ≈ tr(Πp*projector_y1)/ length(y1_index) *(y+1/ϕ)/(ϕ+1/ϕ) + tr(Πp*projector_yτ)/ length(yτ_index) *(ϕ - y)/(ϕ+1/ϕ) atol=1e-6
end

@testset "Topological charge sharpening" begin
    L = 6
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    Y = topological_charge_operator(model)
    ϕ = (1 + √5) / 2

    values, vectors = eigen(Symmetric(Y))
    state_y₁ = vectors[:, findfirst(x -> isapprox(x, ϕ), values)] # A y=1 eigenstate of the topological charge operator.
    state_yτ = vectors[:, findfirst(x -> isapprox(x, -inv(ϕ)), values)] # A y=τ eigenstate of the topological charge operator.

    @testset "bulk evolution Y expectation" begin
        born_config = MeasureConfig(
            τ = 0.7,
            t₂ = 3,
            rng = MersenneTwister(2026),
            mode = :Born,
            enable_τ_eff = false,
            track_y_expectation = true,
        )
        born_outcome = bulk_evolution(model, state_y₁, born_config)

        @test length(born_outcome.y_expectation_values) == born_config.t₂
        @test born_outcome.y_expectation_values ≈ fill(ϕ, born_config.t₂) atol = 1e-6

        sample_config = MeasureConfig(
            τ = born_config.τ,
            t₂ = born_config.t₂,
            mode = :sample,
            enable_τ_eff = false,
            track_y_expectation = true,
        )
        sample_outcome =
            bulk_evolution(model, state_y₁, sample_config, born_outcome.samples)
        unnormalized_outcome =
            bulk_evolution(model, state_y₁, sample_config, born_outcome.samples, false)

        @test sample_outcome.y_expectation_values ≈
              born_outcome.y_expectation_values atol = 1e-6
        @test unnormalized_outcome.y_expectation_values ≈
              born_outcome.y_expectation_values atol = 1e-6

        tracking_off_config = MeasureConfig(
            τ = born_config.τ,
            t₂ = born_config.t₂,
            mode = :sample,
            enable_τ_eff = false,
        )
        tracking_off_outcome =
            bulk_evolution(model, state_y₁, tracking_off_config, born_outcome.samples)
        @test isempty(tracking_off_outcome.y_expectation_values)
    end

    config = MeasureConfig(
        # At τ = 0 every local Kraus operator is proportional to the identity.
        # A Born record therefore has identical likelihood in the y = 1 and
        # y = τ sectors and cannot update their posterior probabilities.
        τ = 0.0,
        t₂ = 3,
        rng = MersenneTwister(1234),
        mode = :Born,
        enable_τ_eff = false,
    )
    # Here p₁ = pτ = 1/2. Since τ = 0 supplies no charge information, this
    # maximally mixed ancilla remains unchanged after every measurement period:
    # Sₐ(t) = log(2), rather than merely Sₐ(0) = log(2).
    balanced_state = (state_y₁ + state_yτ) / √2
    balanced_outcome = topological_charge_sharpening(model, balanced_state, config)
    balanced_ancilla = reference_rdm(model, [1], balanced_outcome.state)

    @test balanced_ancilla ≈ Matrix{Float64}(I, 2, 2) / 2 atol = 1e-10
    @test balanced_outcome.entanglement_entropys ≈ fill(log(2), 3) atol = 1e-6
    @test size(balanced_outcome.samples) == (6, L ÷ 2)

    # The dimension-weighted prior is p_q = dim(H_q)/dim(H). For L = 6 the
    # two Y sectors have dimensions 5 and 13, giving Sₐ(0) ≈ 0.5908. Again,
    # τ = 0 keeps this value fixed, indirectly testing the t = 0 entropy;
    # `entanglement_entropys` itself starts after the first full period.
    n_y₁ = count(x -> isapprox(x, ϕ), values)
    n_yτ = count(x -> isapprox(x, -inv(ϕ)), values)
    p_y₁ = n_y₁ / (n_y₁ + n_yτ)
    p_yτ = n_yτ / (n_y₁ + n_yτ)
    dimension_weighted_state = √(p_y₁) .* state_y₁ + √(p_yτ) .* state_yτ
    dimension_weighted_outcome =
        topological_charge_sharpening(model, dimension_weighted_state, config)
    dimension_weighted_ancilla =
        reference_rdm(model, [1], dimension_weighted_outcome.state)
    initial_charge_entropy = -p_y₁ * log(p_y₁) - p_yτ * log(p_yτ)

    @test (n_y₁, n_yτ) == (5, 13)
    @test dimension_weighted_ancilla ≈ [p_yτ 0.0; 0.0 p_y₁] atol = 1e-10
    @test initial_charge_entropy ≈ 0.5908 atol = 1e-4
    @test dimension_weighted_outcome.entanglement_entropys ≈
          fill(initial_charge_entropy, 3) atol = 1e-6

    sharpening_config = MeasureConfig(
        τ = 1.0,
        t₂ = 4,
        rng = MersenneTwister(1234),
        mode = :Born,
        enable_τ_eff = false,
    )
    # A single Born trajectory performs Bayesian updates of the sector weights.
    # Its entropy need not be monotone for every outcome, but this fixed-seed
    # trajectory sharpens the charge to Sₐ(t_final) < Sₐ(0)/2.
    sharpening_outcome =
        topological_charge_sharpening(model, dimension_weighted_state, sharpening_config)
    @test all(x -> 0 <= x <= log(2), sharpening_outcome.entanglement_entropys)
    @test sharpening_outcome.entanglement_entropys[end] < initial_charge_entropy / 2

    # Replay the Born-generated record as a fixed post-selected trajectory.
    # For an identical sequence of Kraus outcomes, :sample and :Born must give
    # the same normalized state, free energies, and reference-entropy history.
    sample_config = MeasureConfig(
        τ = 1.0,
        t₂ = 4,
        mode = :sample,
        enable_τ_eff = false,
    )
    sample_outcome = topological_charge_sharpening(
        model,
        dimension_weighted_state,
        sample_config,
        sharpening_outcome.samples,
    )

    @test sample_outcome.samples == sharpening_outcome.samples
    @test sample_outcome.state ≈ sharpening_outcome.state atol = 1e-10
    @test sample_outcome.free_energys ≈ sharpening_outcome.free_energys atol = 1e-6
    @test sample_outcome.entanglement_entropys ≈
          sharpening_outcome.entanglement_entropys atol = 1e-6
    @test_throws ErrorException topological_charge_sharpening(
        model,
        dimension_weighted_state,
        sample_config,
    )

    # A definite y = 1 charge maps to the pure ancilla state |1⟩, so its
    # reference entropy is zero before and after the τ = 0 evolution.
    sector_outcome = topological_charge_sharpening(model, state_y₁, config)
    sector_ancilla = reference_rdm(model, [1], sector_outcome.state)

    @test sector_ancilla ≈ [0.0 0.0; 0.0 1.0] atol = 1e-10
    @test sector_outcome.entanglement_entropys ≈ zeros(3) atol = 1e-7

    # When the initial state is a "product" state. The initial charge entropy is not 0.5908.
    initial_st = zeros(length(anyon_basis(model))); initial_st[1] = 1.0
    product_outcome = topological_charge_sharpening(model, initial_st, config)
    product_ancilla = reference_rdm(model, [1], product_outcome.state)
    @test ee(product_ancilla) ≈ 0.611974857421813 atol = 1e-10
end
