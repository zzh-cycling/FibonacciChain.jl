using FibonacciChain
using Test
using LinearAlgebra
using BitBasis
using Random

@testset "Heisenberg basis" begin
    N = 4
    basis_obc = anyon_basis(AnyonModel(SpinHalf(), N; model_type=:Heisenberg, pbc = false))
    @test length(basis_obc) == 2^N
    @test basis_obc[1] == BitStr{N}(0b0000)
    @test basis_obc[end] == BitStr{N}(0b1111)
    # spin-1/2 basis does not depend on boundary condition
    @test anyon_basis(AnyonModel(SpinHalf(), N; model_type=:Heisenberg, pbc = true)) == basis_obc
end

@testset "Heisenbergmap" begin
    N = 3
    T = BitStr{N,Int}

    # Aligned pair on bond (1,2): diagonal +Jz only
    out = FibonacciChain.Heisenbergmap(T(bit"000"), 1, false)
    @test out == (T(bit"000"), T(bit"000"), 1.0, 0.0)

    # Anti-aligned pair on bond (1,2): diagonal -Jz, exchange with 2J
    out = FibonacciChain.Heisenbergmap(T(bit"010"), 1, false)
    @test out == (T(bit"010"), T(bit"100"), -1.0, 2.0)

    # Field enters the diagonal weight
    out = FibonacciChain.Heisenbergmap(T(bit"010"), 1, false; h = 0.5)
    @test out == (T(bit"010"), T(bit"100"), -1.5, 2.0)  # z_1 = +1

    # OBC last site: no bond, field term only
    out = FibonacciChain.Heisenbergmap(T(bit"001"), N, false; h = 0.5)
    @test out == (T(bit"001"), T(bit"001"), 0.5, 0.0)  # z_3 = -1

    # PBC wraps the bond (N,1)
    out = FibonacciChain.Heisenbergmap(T(bit"001"), N, true)
    @test out == (T(bit"001"), T(bit"100"), -1.0, 2.0)

    # Anisotropic coupling
    out = FibonacciChain.Heisenbergmap(T(bit"010"), 1, false; J = 1.2, Jz = 0.7)
    @test out == (T(bit"010"), T(bit"100"), -0.7, 2.4)
end

@testset "actingHam Heisenberg" begin
    T = BitStr{2,Int}

    # N=2 OBC, isotropic: H = XX + YY + ZZ
    model = AnyonModel(SpinHalf(), 2; model_type=:Heisenberg, pbc = false)
    output = FibonacciChain.actingHam(model, bit"00")
    @test output == Dict(T(bit"00") => 1.0)
    output = FibonacciChain.actingHam(model, bit"01")
    @test output == Dict(T(bit"01") => -1.0, T(bit"10") => 2.0)

    # N=2 OBC with field: H = XX + YY + ZZ - h(Z1 + Z2)
    model = AnyonModel(SpinHalf(), 2; model_type=:Heisenberg, pbc = false, h = 0.5)
    output = FibonacciChain.actingHam(model, bit"00")
    @test output == Dict(T(bit"00") => 0.0)  # 1 - 0.5 - 0.5
    output = FibonacciChain.actingHam(model, bit"01")
    @test output == Dict(T(bit"01") => -1.0, T(bit"10") => 2.0)  # -1 - 0.5 + 0.5

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

    # Reference kron construction for N = 3, anisotropic with field
    J, Jz, h = 1.2, 0.7, 0.4
    H_obc_ref =
        J * (X ⊗ X ⊗ Id + Id ⊗ X ⊗ X) + J * (Y ⊗ Y ⊗ Id + Id ⊗ Y ⊗ Y) +
        Jz * (Z ⊗ Z ⊗ Id + Id ⊗ Z ⊗ Z) -
        h * (Z ⊗ Id ⊗ Id + Id ⊗ Z ⊗ Id + Id ⊗ Id ⊗ Z)
    model = AnyonModel(SpinHalf(), 3; model_type=:Heisenberg, pbc = false, J = J, Jz = Jz, h = h)
    @test anyon_ham(model) ≈ H_obc_ref

    # PBC adds the (3,1) bond
    H_pbc_ref = H_obc_ref + J * (X ⊗ Id ⊗ X) + J * (Y ⊗ Id ⊗ Y) + Jz * (Z ⊗ Id ⊗ Z)
    model_pbc = AnyonModel(SpinHalf(), 3; model_type=:Heisenberg, pbc = true, J = J, Jz = Jz, h = h)
    @test anyon_ham(model_pbc) ≈ H_pbc_ref

    # Jz defaults to J (isotropic XXX)
    @test anyon_ham(AnyonModel(SpinHalf(), 3; model_type=:Heisenberg, pbc = false, J = J)) ≈
          anyon_ham(AnyonModel(SpinHalf(), 3; model_type=:Heisenberg, pbc = false, J = J, Jz = J))

    # Singlet/triplet spectrum of the isotropic two-site chain
    @test eigvals(anyon_ham(AnyonModel(SpinHalf(), 2; model_type=:Heisenberg, pbc = false))) ≈
          [-3.0, 1.0, 1.0, 1.0]

    # Hermiticity and dense/sparse consistency
    H = anyon_ham(AnyonModel(SpinHalf(), 4; model_type=:Heisenberg, pbc = true, Jz = 0.3, h = 0.2))
    @test ishermitian(H)
    H_sparse = anyon_ham_sparse(AnyonModel(SpinHalf(), 4; model_type=:Heisenberg, pbc = true, Jz = 0.3, h = 0.2))
    @test norm(Matrix(H_sparse) - H) < 1e-10
end

@testset "anyon_ham Heisenberg momentum sectors" begin
    # The union of all momentum-sector spectra reproduces the full spectrum
    model = AnyonModel(SpinHalf(), 4; model_type=:Heisenberg, pbc = true, J = 1.0, h = 0.3)
    E_full = eigvals(anyon_ham(model))
    E_sectors = vcat([eigvals(anyon_ham(model, k)) for k = 0:3]...)
    @test sort(E_sectors) ≈ E_full
end

@testset "Heisenberg measurement operators" begin
    N = 3
    τ = 1.0
    # :X and :ZZ measurement matrices coincide with the Ising ones (shared spin-1/2 basis)
    for mop in (:X, :ZZ)
        model_h = AnyonModel(SpinHalf(), N; model_type=:Heisenberg, pbc = true, measure_operator = mop)
        model_i = AnyonModel(SpinHalf(), N; model_type=:Ising, pbc = true, measure_operator = mop)
        @test FibonacciChain.measure_matrix(model_h, τ, 2, false) ==
              FibonacciChain.measure_matrix(model_i, τ, 2, false)
        @test FibonacciChain.measure_matrix(model_h, τ, 2, true) ==
              FibonacciChain.measure_matrix(model_i, τ, 2, true)
    end

    # Weak-measurement circuit runs end-to-end (Ising-style X/ZZ alternating layers)
    model = AnyonModel(SpinHalf(), N; model_type=:Heisenberg, pbc = true)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0
    config = MeasureConfig(τ = 1e3, t₂ = 2, rng = MersenneTwister(100), mode = :Born)
    outcome = bulk_evolution(model, st, config)
    @test size(outcome.samples) == (4, N)  # 2 layers per period × 2 periods
    @test length(outcome.state) == 2^N
    @test norm(outcome.state) ≈ 1.0
end
