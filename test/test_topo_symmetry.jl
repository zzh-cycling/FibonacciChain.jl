using FibonacciChain
using LinearAlgebra
using JLD
using Test

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[7] = log(1 + √2)  # atanh(1/√2) = log(1 + √2)
τlis[end] = 1000.0     # Last value is for γ=1

function get_cfg_params_Born(ind, L)
    cfg = Dict(
        1 => (2500L, 1000, 750L),
        2 => (500L, 100, 120L),
        3 => (200L, 40, 80L),
        4 => (100L, 40, 40L),
        5 => (80L, 32, 20L),
        6 => (45L, 20, 15L),
        7 => (35L, 14, 10L),
        8 => (25L, 10, 5L),
        9 => (8L, 4, 2L),
        10 => (8L, 4, 2L),
        11 => (5L, 2, 1L),
    )
    D, step, start = get(cfg, ind, (5L, 2, L))
    inds = collect(1:step:div(D, 2))
    avg_range = start:(div(D, 2)-5)
    return D, inds, avg_range
end


function check_commute_tf_Born(
    L::Int,
    τ_idx::Int,
    index::Int;
    )   
    τ = τlis[τ_idx]
    D, inds, avg_range = get_cfg_params_Born(τ_idx, L)
    t = div(D, 2L)
    path = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τ_idx)/t$(t)_samples$(index).jld"
    data = load(path)
    sample = data["sample"]
    sample = similar(sample)
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    
    l = length(anyon_basis(model))
    M = zeros(l, l)

    for i in 1:l
        st = zeros(length(anyon_basis(model)))
        st[i] = 1.0
        config = MeasureConfig(τ = τ, mode = :sample, t₂ = t*L, enable_τ_eff=false)
        outcome = bulk_evolution(model, st, config, sample, false)
        final_state = outcome.state
        M[i, :] = final_state
    end
    
    Y = topological_charge_operator(model)
    
    M_normalized = M ./ minimum(abs.(M))
    return norm(M_normalized*Y - Y*M_normalized)
end

function check_commute_tf_ps(
    L::Int,
    τ_idx::Int,
    )   
    τ = τlis[τ_idx]
    D, inds, avg_range = get_cfg_params_Born(τ_idx, L)
    t = div(D, 2L)
    
    sample = BitMatrix(ones(D, div(L,2)))
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    
    l = length(anyon_basis(model))
    M = zeros(l, l)

    for i in 1:l
        st = zeros(length(anyon_basis(model)))
        st[i] = 1.0
        config = MeasureConfig(τ = τ, mode = :sample, t₂ = t*L, enable_τ_eff=false)
        outcome = bulk_evolution(model, st, config, sample, false)
        final_state = outcome.state
        M[i, :] = final_state
    end
    
    Y = topological_charge_operator(model)
    
    M_normalized = M ./ minimum(abs.(M))
    return norm(M_normalized*Y - Y*M_normalized)
end

function compute_state_components(L)
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0
    Y = topological_charge_operator(model)
    stp  = (Y*st .+ st ./ ϕ) ./ sqrt(5)   # P₊|0000⟩
    stm  = (ϕ .* st .- Y*st) ./ sqrt(5)   # P₊|0000⟩
    w₊, w₋ = norm(stp)^2, norm(stm)^2
    return w₊, w₋
    
    # vals= eigvals(Y)
    # ϕ = (1 + √5)/2
    # ϕ_count= count(x->isapprox(x, ϕ), vals)
    # mϕ_count= count(x->isapprox(x, 1-ϕ), vals)
    # @show ϕ_count/mϕ_count
end


function measure_eigenvalues(L)
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    Y = topological_charge_operator(model)
    τ = 1000.0
    ϕ = (1 + √5) / 2
    Πplis = [FibonacciChain.measure_matrix(model, τ, idx, true) for idx = 1:1] # s=1

    Πmlis = [FibonacciChain.measure_matrix(model, τ, idx, false) for idx = 1:1] # s=τ
    
    energy, states = eigen((Y+Y')/2)
    y1_index = findall(x -> isapprox(x, ϕ), energy)
    yτ_index = findall(x -> isapprox(x, 1-ϕ), energy)
    Πp = Πplis[1]
    Πm = Πmlis[1]
    GS = states[:, 1]
    Πp_GS = GS' * (Πp * GS)
    Πm_GS = GS' * (Πm * GS)

    projector_y1 = sum([states[:, i] * states[:, i]' for i in y1_index])
    projector_yτ = sum([states[:, i] * states[:, i]' for i in yτ_index])
    # return Πp_GS, Πm_GS

    return tr(Πp * projector_y1)/ length(y1_index), tr(Πp*projector_yτ)/ length(yτ_index), tr(Πm * projector_y1)/ length(y1_index), tr(Πm*projector_yτ)/ length(yτ_index)

    y = (-1/ϕ)^L
    theoretical_Πm = (y+1/ϕ)/(ϕ+1/ϕ)^2*((Fibonacci_number(L)+Fibonacci_number(L-2))/ϕ -y*ϕ)/Fibonacci_number(L-1) + 
    (ϕ-y)/(ϕ+1/ϕ)^2*((Fibonacci_number(L)+Fibonacci_number(L-2))*ϕ + y*ϕ)/Fibonacci_number(L+1)
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

@testset begin
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

function Fibonacci_number(L)
    return (1/sqrt(5)) * (((1+sqrt(5))/2)^L - ((1-sqrt(5))/2)^L)
end

function Lucas_number(L)
    return Fibonacci_number(L+1) + Fibonacci_number(L-1)
end

function check_dynamics(L::Int, τ_idx::Int, index::Int)
    τ = τlis[τ_idx]
    D, inds, avg_range = get_cfg_params_Born(τ_idx, L)
    t = div(D, 2L)
    path = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τ_idx)/t$(t)_samples$(index).jld"
    data = load(path)
    sample = data["sample"]
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0
    config = MeasureConfig(τ = τ, mode = :sample, t₂ = t*L, enable_τ_eff=false)
    outcome = bulk_evolution(model, st, config, sample, false)
    final_state = outcome.state

    normalized_final_state = final_state ./ norm(final_state)
    Y = topological_charge_operator(model)
    initial_y = st' * (Y * st)
    final_y = normalized_final_state' * (Y * normalized_final_state)

    ϕ = (1 + √5) / 2
    return (initial_y+1/ϕ)/(ϕ+1/ϕ), (final_y+1/ϕ)/(ϕ+1/ϕ)
end