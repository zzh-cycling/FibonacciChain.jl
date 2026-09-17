using FibonacciChain
using Test
using LinearAlgebra
using BitBasis
using Random
using Arpack

@testset "Heisenberg basis" begin
    N = 4
    basis_obc = anyon_basis(AnyonModel(SpinHalf(), N; model_type=:Heisenberg, pbc = false))
    @test length(basis_obc) == 2^N
    @test basis_obc[1] == BitStr{N}(0b0000)
    @test basis_obc[end] == BitStr{N}(0b1111)
    # spin-1/2 basis does not depend on boundary condition
    @test anyon_basis(AnyonModel(SpinHalf(), N; model_type=:Heisenberg, pbc = true)) == basis_obc

    # Only the :SWAP measurement operator is allowed (and is the default)
    @test AnyonModel(SpinHalf(), N; model_type=:Heisenberg).measure_operator == :SWAP
    @test_throws AssertionError AnyonModel(SpinHalf(), N; model_type=:Heisenberg, measure_operator = :X)
end

@testset "Heisenbergmap" begin
    N = 3
    T = BitStr{N,Int}

    # Aligned pair on bond (1,2): diagonal +JΔ only
    out = FibonacciChain.Heisenbergmap(T(bit"000"), 1, false)
    @test out == (T(bit"000"), T(bit"000"), 1.0, 0.0)

    # Anti-aligned pair on bond (1,2): diagonal -JΔ, exchange with 2J
    out = FibonacciChain.Heisenbergmap(T(bit"010"), 1, false)
    @test out == (T(bit"010"), T(bit"100"), -1.0, 2.0)

    # No field term: the diagonal weight does not depend on the spectator bit
    out = FibonacciChain.Heisenbergmap(T(bit"011"), 1, false)
    @test out == (T(bit"011"), T(bit"101"), -1.0, 2.0)

    # OBC last site: no bond, zero weights
    out = FibonacciChain.Heisenbergmap(T(bit"001"), N, false)
    @test out == (T(bit"001"), T(bit"001"), 0.0, 0.0)

    # PBC wraps the bond (N,1)
    out = FibonacciChain.Heisenbergmap(T(bit"001"), N, true)
    @test out == (T(bit"001"), T(bit"100"), -1.0, 2.0)

    # Anisotropic coupling
    out = FibonacciChain.Heisenbergmap(T(bit"010"), 1, false; J = 1.2, Δ = 0.7)
    @test out == (T(bit"010"), T(bit"100"), -1.2 * 0.7, 2.4)

    # Ferromagnetic coupling flips all signs
    out = FibonacciChain.Heisenbergmap(T(bit"010"), 1, false; J = -1.0)
    @test out == (T(bit"010"), T(bit"100"), 1.0, -2.0)
end

@testset "actingHam Heisenberg" begin
    T = BitStr{2,Int}

    # N=2 OBC, isotropic AFM: H = XX + YY + ZZ
    model = AnyonModel(SpinHalf(), 2; model_type=:Heisenberg, pbc = false)
    output = FibonacciChain.actingHam(model, bit"00")
    @test output == Dict(T(bit"00") => 1.0)
    output = FibonacciChain.actingHam(model, bit"01")
    @test output == Dict(T(bit"01") => -1.0, T(bit"10") => 2.0)

    # Anisotropic: H = J(XX + YY + Δ ZZ)
    model = AnyonModel(SpinHalf(), 2; model_type=:Heisenberg, pbc = false, Δ = 0.5)
    output = FibonacciChain.actingHam(model, bit"00")
    @test output == Dict(T(bit"00") => 0.5)
    output = FibonacciChain.actingHam(model, bit"01")
    @test output == Dict(T(bit"01") => -0.5, T(bit"10") => 2.0)

    # Ferromagnetic: J = -1
    model = AnyonModel(SpinHalf(), 2; model_type=:Heisenberg, pbc = false, J = -1.0)
    output = FibonacciChain.actingHam(model, bit"01")
    @test output == Dict(T(bit"01") => 1.0, T(bit"10") => -2.0)

    # N=2 PBC: the single bond is counted twice (i=1 and i=2)
    model = AnyonModel(SpinHalf(), 2; model_type=:Heisenberg, pbc = true)
    output = FibonacciChain.actingHam(model, bit"01")
    @test output == Dict(T(bit"01") => -2.0, T(bit"10") => 4.0)
end

@testset "anyon_ham Heisenberg" begin
    X = Float64[0 1; 1 0]
    Y = ComplexF64[0 -im; im 0]
    Z = Float64[1 0; 0 -1]
    Id = Float64[1 0; 0 1]
    ⊗(a, b) = kron(a, b)

    # Reference kron construction for N = 3, anisotropic (no field term)
    J, Δ = 1.2, 0.7
    H_obc_ref =
        J * (X ⊗ X ⊗ Id + Id ⊗ X ⊗ X) + J * (Y ⊗ Y ⊗ Id + Id ⊗ Y ⊗ Y) +
        J * Δ * (Z ⊗ Z ⊗ Id + Id ⊗ Z ⊗ Z)
    model = AnyonModel(SpinHalf(), 3; model_type=:Heisenberg, pbc = false, J = J, Δ = Δ)
    @test anyon_ham(model) ≈ H_obc_ref

    # PBC adds the (3,1) bond
    H_pbc_ref = H_obc_ref + J * (X ⊗ Id ⊗ X) + J * (Y ⊗ Id ⊗ Y) + J * Δ * (Z ⊗ Id ⊗ Z)
    model_pbc = AnyonModel(SpinHalf(), 3; model_type=:Heisenberg, pbc = true, J = J, Δ = Δ)
    @test anyon_ham(model_pbc) ≈ H_pbc_ref

    # Δ defaults to 1 (isotropic XXX)
    @test anyon_ham(AnyonModel(SpinHalf(), 3; model_type=:Heisenberg, pbc = false, J = J)) ≈
          anyon_ham(AnyonModel(SpinHalf(), 3; model_type=:Heisenberg, pbc = false, J = J, Δ = 1.0))

    # Singlet/triplet spectrum of the isotropic two-site chain
    @test eigvals(anyon_ham(AnyonModel(SpinHalf(), 2; model_type=:Heisenberg, pbc = false))) ≈
          [-3.0, 1.0, 1.0, 1.0]

    # Hermiticity and dense/sparse consistency
    H = anyon_ham(AnyonModel(SpinHalf(), 4; model_type=:Heisenberg, pbc = true, Δ = 0.3))
    @test ishermitian(H)
    H_sparse = anyon_ham_sparse(AnyonModel(SpinHalf(), 4; model_type=:Heisenberg, pbc = true, Δ = 0.3))
    @test norm(Matrix(H_sparse) - H) < 1e-10
end

@testset "anyon_ham Heisenberg momentum sectors" begin
    # The union of all momentum-sector spectra reproduces the full spectrum
    model = AnyonModel(SpinHalf(), 4; model_type=:Heisenberg, pbc = true, J = 1.0, Δ = 0.3)
    E_full = eigvals(anyon_ham(model))
    E_sectors = vcat([eigvals(anyon_ham(model, k)) for k = 0:3]...)
    @test sort(E_sectors) ≈ E_full
end

@testset "Heisenberg SWAP measurement" begin
    # The :SWAP bond gate at strength τ applies exp(s·τ·h_bond) with
    # h_bond = J (XX + YY + Δ ZZ), s = +1 for sign = false and -1 for sign = true.
    X = Float64[0 1; 1 0]
    Y = ComplexF64[0 -im; im 0]
    Z = Float64[1 0; 0 -1]
    ⊗(a, b) = kron(a, b)

    τ = 0.7
    J, Δ = 1.3, 0.6
    h4 = J * (X ⊗ X + Y ⊗ Y + Δ * (Z ⊗ Z))
    model = AnyonModel(SpinHalf(), 2; model_type = :Heisenberg, pbc = false, J = J, Δ = Δ)
    @test FibonacciChain.measure_matrix(model, τ, 1, false) ≈ exp(τ * h4)
    @test FibonacciChain.measure_matrix(model, τ, 1, true) ≈ exp(-τ * h4)

    # Hand-computed weights on the anti-aligned state |01⟩ at τ = 0.5, J = Δ = 1:
    # h_bond = -I + 2σ_x on the {|01⟩, |10⟩} block
    model = AnyonModel(SpinHalf(), 2; model_type = :Heisenberg, pbc = false)
    T = BitStr{2,Int}
    s1, s2, w1, w2 = measure_basismap(model, 0.5, T(bit"01"), 1, false)
    @test (s1, s2) == (T(bit"01"), T(bit"10"))
    @test w1 ≈ exp(-0.5) * cosh(1.0)
    @test w2 ≈ exp(-0.5) * sinh(1.0)
    # sign = true flips the exponent and the off-diagonal branch sign
    s1, s2, w1, w2 = measure_basismap(model, 0.5, T(bit"01"), 1, true)
    @test w1 ≈ exp(0.5) * cosh(1.0)
    @test w2 ≈ -exp(0.5) * sinh(1.0)
    # Aligned states are eigenstates of the bond term
    s1, s2, w1, w2 = measure_basismap(model, 0.5, T(bit"00"), 1, false)
    @test (s1, s2) == (T(bit"00"), T(bit"00"))
    @test w1 ≈ exp(0.5) && w2 == 0.0

    # OBC: no bond starts at site N
    @test_throws AssertionError measure_basismap(model, 0.5, T(bit"01"), 2, false)
end

@testset "Heisenberg SWAP measurement circuit" begin
    N = 4
    model = AnyonModel(SpinHalf(), N; model_type = :Heisenberg, pbc = true)

    # Staggered SWAP-bond layers: odd layers on bonds (1,2), (3,4), ...; even on (2,3), (4,1)
    sites1, m1, τ1 = FibonacciChain._obtain_measurement_config(model, 1, 0.5)
    @test sites1 == [1, 3]
    @test m1.measure_operator == :SWAP && τ1 == 0.5
    sites2, m2, _ = FibonacciChain._obtain_measurement_config(model, 2, 0.5)
    @test sites2 == [2, 4]
    @test m2.params[:J] == 1.0 && m2.params[:Δ] == 1.0

    # OBC: the even layer loses the dangling bond start i = N
    model_obc = AnyonModel(SpinHalf(), N; model_type = :Heisenberg, pbc = false, Δ = 0.5)
    sites1, _, _ = FibonacciChain._obtain_measurement_config(model_obc, 1, 0.5)
    @test sites1 == [1, 3]
    sites2, m2o, _ = FibonacciChain._obtain_measurement_config(model_obc, 2, 0.5)
    @test sites2 == [2]
    @test m2o.params[:Δ] == 0.5
    @test FibonacciChain._get_sample_column_indices(model_obc, 2) == [1]

    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0
    config = MeasureConfig(τ = 0.5, t₂ = 2, rng = MersenneTwister(100), mode = :Born)
    outcome = bulk_evolution(model, st, config)
    @test size(outcome.samples) == (4, N ÷ 2)  # 2 layers per period × 2 periods, one sample per bond
    @test length(outcome.state) == 2^N
    @test norm(outcome.state) ≈ 1.0
    outcome_obc = bulk_evolution(model_obc, st, config)
    @test size(outcome_obc.samples) == (4, N ÷ 2)
    @test norm(outcome_obc.state) ≈ 1.0

    # Sample-mode (postselected) evolution consumes the same per-layer layout
    outcome_s = bulk_evolution(model, st, MeasureConfig(τ = 0.5, t₂ = 2, mode = :sample), falses(4, N ÷ 2))
    @test norm(outcome_s.state) ≈ 1.0

    # The transfer-matrix path dispatches generically through measure_basismap
    tm = transfer_matrix(model, 0.5, falses(4, N ÷ 2))
    @test size(tm) == (2^N, 2^N)

    # AFM Born run from a Néel state generates entanglement
    N6 = 6
    model6 = AnyonModel(SpinHalf(), N6; model_type = :Heisenberg, pbc = true)
    basis6 = anyon_basis(model6)
    st6 = zeros(length(basis6))
    st6[findfirst(==(BitStr{N6}(0b010101)), basis6)] = 1.0
    out6 = bulk_evolution(model6, st6, MeasureConfig(τ = 0.5, t₂ = 3, rng = MersenneTwister(42), mode = :Born))
    @test norm(out6.state) ≈ 1.0
    ee6 = anyon_eelis(model6, out6.state)
    @test all(>=(0), ee6) && maximum(ee6) > 0.5
end

@testset "Heisenberg AFM central charge" begin
    # Ground-state entanglement scaling of the critical AFM chain: plain linear
    # regression of S(x) against x_c = (1/3) log[(L/π) sin(πx/L)] on interior
    # cuts; the slope is the central charge c.
    L = 16
    cuts = collect(3:(L-3))
    xc = (1 / 3) .* log.(L / π .* sin.(π .* cuts ./ L))
    Xmat = hcat(ones(length(cuts)), xc)
    # Measured at L = 16: c ≈ 1.012 for Δ = 0 and c ≈ 1.021 for Δ = 1
    # (marginal-operator log corrections); tolerances are set around these values.
    for (Δ, tol) in ((0.0, 0.02), (1.0, 0.05))
        model = AnyonModel(SpinHalf(), L; model_type = :Heisenberg, pbc = true, J = 1.0, Δ = Δ)
        H = anyon_ham_sparse(model)
        E, vecs = Arpack.eigs(H, nev = 1, which = :SR)
        S = anyon_eelis(model, vecs[:, 1])
        c = (Xmat \ S[cuts])[2]
        @test abs(c - 1) < tol
    end
end

@testset "Heisenberg FM ground state degeneracy" begin
    # FM isotropic chain: the fully symmetric spin-L/2 multiplet gives L+1
    # exactly degenerate ground states, separated by a finite magnon gap.
    L = 8
    model = AnyonModel(SpinHalf(), L; model_type = :Heisenberg, pbc = true, J = -1.0, Δ = 1.0)
    H = anyon_ham_sparse(model)
    # ncv is raised above the default so Arpack resolves the full degenerate multiplet
    E = sort(real.(Arpack.eigs(H, nev = L + 2, ncv = 40, which = :SR)[1]))
    @test E[1] ≈ -Float64(L)  # every bond in the triplet state contributes JΔ = -1
    @test maximum(abs.(E[1:(L+1)] .- E[1])) < 1e-8
    @test E[L+2] - E[1] > 0.1  # measured ≈ 1.17 at L = 8
end

@testset "Heisenberg space-time isotropic bond gate" begin
    X = [0 1; 1 0]
    Y = [0 -im; im 0]
    Z = [1 0; 0 -1]
    h4(Δ) = kron(X, X) + kron(Y, Y) + Δ * kron(Z, Z)
    I4 = Matrix{ComplexF64}(I, 4, 4)

    # Space-direction dual: view U as the tensor T[o1, o2, i1, i2] (outputs o1, o2;
    # inputs i1, i2) and re-pair the legs of site 1 against those of site 2,
    # Ũ[(o1, i1), (o2, i2)] = T[o1, o2, i1, i2].
    dual(U) = reshape(permutedims(reshape(U, 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)

    # At the full-swap angle θ = π/4 the gate is dual-unitary for every Δ: the
    # two-qubit dual-unitary classification fixes only the XX and YY angles to
    # π/4 (the iSWAP family), leaving the ZZ angle — here πΔ/4 — unconstrained.
    # At Δ = 0 the gate is iSWAP†, at Δ = 1 it is SWAP up to a phase, and the
    # duals take the same form.
    U0 = exp(-im * π / 4 * h4(0.0))
    @test U0 * U0' ≈ I4
    @test dual(U0) * dual(U0)' ≈ I4 atol = 1e-10
    @test dual(U0) ≈ ComplexF64[1 0 0 0; 0 0 -im 0; 0 -im 0 0; 0 0 0 1] atol = 1e-12

    U1 = exp(-im * π / 4 * h4(1.0))
    @test U1 * U1' ≈ I4
    @test dual(U1) * dual(U1)' ≈ I4 atol = 1e-10
    @test dual(U1) ≈ exp(-im * π / 4) * ComplexF64[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1] atol = 1e-12

    # Δ = 0.5 at θ = π/4 is dual-unitary as well (same iSWAP-family point) ...
    U5 = exp(-im * π / 4 * h4(0.5))
    @test dual(U5) * dual(U5)' ≈ I4 atol = 1e-10

    # ... whereas a partial-swap gate (θ = π/8) is genuinely not dual-unitary
    U8 = exp(-im * π / 8 * h4(1.0))
    @test U8 * U8' ≈ I4
    @test maximum(abs.(dual(U8) * dual(U8)' - I4)) > 1e-2
end
