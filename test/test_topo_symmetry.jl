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
