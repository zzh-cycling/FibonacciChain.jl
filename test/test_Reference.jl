using FibonacciChain
using Test
using BitBasis
using LinearAlgebra 
using Random

@testset "build_extended_basis" begin
    N = 3
    k_total = 2
    T = BitStr{N+k_total, Int}
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    basis = anyon_basis(model)
    extended_basis = FibonacciChain.build_extended_basis(k_total,basis)

    @test length(extended_basis) == 4 * length(basis)
    @test extended_basis == T.([0, 1, 2, 4, 8, 9, 10, 12, 16, 17, 18, 20, 24, 25, 26, 28])

    k_total = 1
    T = BitStr{N+k_total, Int}
    extended_basis = FibonacciChain.build_extended_basis(k_total, basis)

    @test length(extended_basis) == 2 * length(basis)
    @test extended_basis == T.([0, 1, 2, 4, 8, 9, 10, 12])

    k_total = 0
    extended_basis = FibonacciChain.build_extended_basis(k_total, basis)
    @test length(extended_basis) == length(basis)
    @test extended_basis == basis

    # Test for Ising basis
    model_Ising = AnyonModel(IsingAnyon(), N, pbc=true)
    basis_ising = anyon_basis(model_Ising)
    extended_basis_ising = FibonacciChain.build_extended_basis(0, basis_ising)
    @testset extended_basis_ising == basis_ising

    extended_basis_ising1 = FibonacciChain.build_extended_basis(1, basis_ising)
    @test extended_basis_ising1 == anyon_basis(AnyonModel(IsingAnyon(), N+1, pbc=true))

    extended_basis_ising2 = FibonacciChain.build_extended_basis(2, basis_ising)
    @test extended_basis_ising2 == anyon_basis(AnyonModel(IsingAnyon(), N+2, pbc=true))
end

@testset "reference_measure_basismap" begin
    # Now main system is Fibonacci anyon
    N = 3
    τ = 1.0
    sign = false
    k_old = 1
    T = BitStr{N, Int}
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    basislis = anyon_basis(model)
    l = length(basislis)
    output11 = FibonacciChain.reference_measure_basismap.(Ref(model), T, T, τ, basislis, 1, sign, k_old=0)
    output12 = FibonacciChain.measure_basismap.(Ref(model), τ, basislis, 1, sign)
    output21 = FibonacciChain.reference_measure_basismap.(Ref(model), T, T, τ, basislis, 2, sign, k_old=0)
    output22 = FibonacciChain.measure_basismap.(Ref(model), τ, basislis, 2, sign)
    output31 = FibonacciChain.reference_measure_basismap.(Ref(model), T, T, τ, basislis, 3, sign, k_old=0)
    output32 = FibonacciChain.measure_basismap.(Ref(model), τ, basislis, 3, sign)
    @test all([all(output11[i] .≈ output12[i]) for i in 1:l])
    @test all([all(output21[i] .≈ output22[i]) for i in 1:l])
    @test all([all(output31[i] .≈ output32[i]) for i in 1:l])

    extended_basis = FibonacciChain.build_extended_basis(1, basislis)
    output13 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+1, Int}, τ, extended_basis, 1, sign, k_old=1)
    output23 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+1, Int}, τ, extended_basis, 2, sign, k_old=1)
    output33 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+1, Int}, τ, extended_basis, 3, sign, k_old=1)

    @test all([all([all(output13[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output23[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output33[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:1])

    extended_basis2 = FibonacciChain.build_extended_basis(2, basislis)
    output14 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+2, Int}, τ, extended_basis2, 1, sign, k_old=2)
    output24 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+2, Int}, τ, extended_basis2, 2, sign, k_old=2)
    output34 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+2, Int}, τ, extended_basis2, 3, sign, k_old=2)

    @test all([all([all(output14[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output24[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output34[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:3])
end

@testset "reference_measure_basismap_Ising" begin
    # Now main system is Ising anyon
    N = 3
    τ = 1000.0
    sign = false
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)
    k_old = 1
    T = BitStr{N, Int}
    basislis = anyon_basis(model)
    l = length(basislis)
    
    output11 = FibonacciChain.reference_measure_basismap.(Ref(model), T, T, τ, basislis, 1, sign, k_old=0)
    output12 = FibonacciChain.measure_basismap.(Ref(model), τ, basislis, 1, sign)
    output21 = FibonacciChain.reference_measure_basismap.(Ref(model), T, T, τ, basislis, 2, sign, k_old=0)
    output22 = FibonacciChain.measure_basismap.(Ref(model), τ, basislis, 2, sign)
    output31 = FibonacciChain.reference_measure_basismap.(Ref(model), T, T, τ, basislis, 3, sign, k_old=0)
    output32 = FibonacciChain.measure_basismap.(Ref(model), τ, basislis, 3, sign)
    @test all([all(output11[i] .≈ output12[i]) for i in 1:l])
    @test all([all(output21[i] .≈ output22[i]) for i in 1:l])
    @test all([all(output31[i] .≈ output32[i]) for i in 1:l])

    extended_basis = FibonacciChain.build_extended_basis(1, basislis)
    output13 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+1}, τ, extended_basis, 1, sign, k_old=1)
    output23 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+1}, τ, extended_basis, 2, sign, k_old=1)
    output33 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+1}, τ, extended_basis, 3, sign, k_old=1)

    @test all([all([all(output13[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output23[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output33[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:1])

    extended_basis2 = FibonacciChain.build_extended_basis(2, basislis)
    output14 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+2, Int}, τ, extended_basis2, 1, sign, k_old=2)
    output24 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+2, Int}, τ, extended_basis2, 2, sign, k_old=2)
    output34 = FibonacciChain.reference_measure_basismap.(Ref(model), T, BitStr{N+2, Int}, τ, extended_basis2, 3, sign, k_old=2)

    @test all([all([all(output14[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output24[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output34[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:3])

    # Test the IsingZZ
    model_ZZ = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:ZZ)
    output12zz = FibonacciChain.measure_basismap.(Ref(model_ZZ), τ, basislis, 1, sign)
    output22zz = FibonacciChain.measure_basismap.(Ref(model_ZZ), τ, basislis, 2, sign)
    output32zz = FibonacciChain.measure_basismap.(Ref(model_ZZ), τ, basislis, 3, sign)

    extended_basis = FibonacciChain.build_extended_basis(1, basislis)
    output13zz = FibonacciChain.reference_measure_basismap.(Ref(model_ZZ), T, BitStr{N+1}, τ, extended_basis, 1, sign, k_old=1)
    output23zz = FibonacciChain.reference_measure_basismap.(Ref(model_ZZ), T, BitStr{N+1}, τ, extended_basis, 2, sign, k_old=1)
    output33zz = FibonacciChain.reference_measure_basismap.(Ref(model_ZZ), T, BitStr{N+1}, τ, extended_basis, 3, sign, k_old=1)

    @test all([all([all(output13zz[i+j*l] .≈ output12zz[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output23zz[i+j*l] .≈ output22zz[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output33zz[i+j*l] .≈ output32zz[i]) for i in 1:l]) for j in 0:1])
end

@testset "reference_measuremap" begin
    N = 3
    τ = 1000.0
    sign = false
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    k_old = 1
    st = ones(4)/2;
    ϕ = (1 + √5) / 2  

    add_st = FibonacciChain.add_reference_qubits_reset(model, st, 1)[3]

    ext_basis = FibonacciChain.build_extended_basis(1, anyon_basis(model))
    output13 = FibonacciChain.reference_measuremap(model, τ, add_st, 1, sign, k_old=1, extended_basis=ext_basis)
    output23 = FibonacciChain.reference_measuremap(model, τ, add_st, 2, sign, k_old=1, extended_basis=ext_basis)
    output33 = FibonacciChain.reference_measuremap(model, τ, add_st, 3, sign, k_old=1, extended_basis=ext_basis)
    @test output13 == 0.5*[(1-ϕ^(-1)), 1, 1, -ϕ^(-3/2), -ϕ^(-3/2), 0, 0, ϕ^(-1)]
    @test output23 == 0.5*[(1-ϕ^(-1)- ϕ^(-3/2)), 1, ϕ^(-1)-ϕ^(-3/2), 0, 0 , 0, 0, 1]
    @test output33 == 0.5*[(1-ϕ^(-1)- ϕ^(-3/2)), ϕ^(-1)-ϕ^(-3/2), 1, 0, 0, 0, 0, 1]

    output13 = FibonacciChain.reference_measuremap(model, τ, st, 1, false, k_old=0, extended_basis=anyon_basis(model))
    output23 = FibonacciChain.reference_measuremap(model, τ, st, 2, false, k_old=0, extended_basis=anyon_basis(model))
    output33 = FibonacciChain.reference_measuremap(model, τ, st, 3, true, k_old=0, extended_basis=anyon_basis(model))
    @test output13 == measuremap(model, τ, st, 1, false)
    @test output23 == measuremap(model, τ, st, 2, false)
    @test output33 == measuremap(model, τ, st, 3, true)
end

@testset "reference_measuremap_Ising" begin
    N = 3
    τ = 1000.0
    sign = false
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)
    k_old = 1
    T = BitStr{N, Int}
    st = zeros(2^N); st[1] = 1 
    add_st = FibonacciChain.add_reference_qubits_reset(model, st, 1)[3]

    ext_basis = FibonacciChain.build_extended_basis(1, anyon_basis(model))
    output13 = FibonacciChain.reference_measuremap(model, τ, add_st, 1, sign, k_old=1, extended_basis=ext_basis)
    output23 = FibonacciChain.reference_measuremap(model, τ, add_st, 2, sign, k_old=1, extended_basis=ext_basis)

    model2 = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:ZZ)
    output33 = FibonacciChain.reference_measuremap(model2, τ, add_st, 3, sign, k_old=1, extended_basis=ext_basis)
    @test output13[[1, 5, 9, 13]] ≈ 1/2√2*ones(4)
    @test output23[[1, 3, 13, 15]] ≈ 1/2√2*ones(4)
    @test output33[[1]] ≈ [1/√2]
end

@testset "concat_bell_pair" begin
    N = 3   
    st = ones(4)/2; 
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    basis_F = anyon_basis(model)
    ext_basis = FibonacciChain.build_extended_basis(1, basis_F)
    bell_st = FibonacciChain.concat_bell_pair(N, st, basis_F, ext_basis, site_idx = 1, k_old = 0)
    @test bell_st[[1,2,3,8]] == 1/2*ones(4)
end

@testset "add_reference_qubits" begin
    N = 3
    st = ones(4)/2; 
    # Test add_reference_qubit reset to 0
    add_model = AnyonModel(FibonacciAnyon(), N+1, pbc=true)
    add_st = FibonacciChain.add_reference_qubits_reset(add_model, ones(7)./√7, 1)[3]
    @test add_st[[collect(1:5)..., 13,14]] ≈ 0.37796447300922725*ones(7)


    add_stp = FibonacciChain.add_reference_qubits_reset(add_model, add_st, 1)[3]
    @test add_stp[[collect(1:5)..., 20,21]] ≈ 0.37796447300922725*ones(7)

    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    add_st1 = FibonacciChain.add_reference_qubits_reset(model, st, 1)[3]
    add_st2 = FibonacciChain.add_reference_qubits_reset(model, st, 2)[3]
    add_st3 = FibonacciChain.add_reference_qubits_reset(model, st, 3)[3]
    @test add_st1 ≈ [0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5] 
    @test add_st2 ≈ [0.5, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0]
    @test add_st3 ≈ [0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.0]

    model2 = AnyonModel(FibonacciAnyon(), N, pbc=true)
    add2_st1 = FibonacciChain.add_reference_qubits_reset(model2, add_st1, 1)[3]
    add2_st2 = FibonacciChain.add_reference_qubits_reset(model2, add_st2, 2)[3]
    add2_st3 = FibonacciChain.add_reference_qubits_reset(model2, add_st3, 3)[3]

    @test add2_st1[[1,2,3,12]] == [0.5, 0.5, 0.5, 0.5]
    @test add2_st2[[1,2,4,11]] == [0.5, 0.5, 0.5, 0.5]
    @test add2_st3[[1,3,4,10]] == [0.5, 0.5, 0.5, 0.5]
end

@testset "add_reference_qubits_Ising" begin
    # Test for Ising basis
    N = 3
    st_ising = zeros(2^N); st_ising[1] = 1/√2; st_ising[end] = 1/√2; # set the GHZ state.
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:Z)
    add_st_ising1 = FibonacciChain.add_reference_qubits_reset(model, st_ising, 1)[3]
    add_st_ising2 = FibonacciChain.add_reference_qubits_reset(model, st_ising, 2)[3]
    add_st_ising3 = FibonacciChain.add_reference_qubits_reset(model, st_ising, 3)[3]
    @test add_st_ising1[[1,13]] == [1/√2, 1/√2]
    @test add_st_ising2[[1,11]] == [1/√2, 1/√2]
    @test add_st_ising3[[1,10]] == [1/√2, 1/√2]

    st_ising = zeros(2^N); st_ising[1] = 1/√2; st_ising[end] = 1/√2; # set the last two qubits to be in the Bell state
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)
    add_st_ising1 = FibonacciChain.add_reference_qubits_reset(model, st_ising, 1)[3]
    add_st_ising2 = FibonacciChain.add_reference_qubits_reset(model, st_ising, 2)[3]
    add_st_ising3 = FibonacciChain.add_reference_qubits_reset(model, st_ising, 3)[3]
    @test add_st_ising1[[1, 4, 13, 16]] == 1/2*ones(4)
    @test add_st_ising2[[1, 6, 11, 16]] == 1/2*ones(4)
    @test add_st_ising3[[1, 7, 10, 16]] == 1/2*ones(4)
end

@testset "add_reference_qubits_Ising_copy" begin
    # Test for Ising basis
    N = 3
    st_ising = zeros(2^N); st_ising[1] = 1/√2; st_ising[end] = 1/√2; # set the GHZ state.
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:Z)
    add_st_ising1 = FibonacciChain.add_reference_qubits(model, st_ising, 1) # For copy way to entangle, the result is irrelevant to measurement outcome.
    add_st_ising2 = FibonacciChain.add_reference_qubits(model, st_ising, 2)
    add_st_ising3 = FibonacciChain.add_reference_qubits(model, st_ising, 3)
    @test add_st_ising1 == add_st_ising2 == add_st_ising3
    @test add_st_ising1[[1,16]] ≈ [1/√2, 1/√2]

    st_ising = ones(2^N); st_ising /= norm(st_ising); # set to be plus state
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)
    add_st_ising1 = FibonacciChain.add_reference_qubits(model, st_ising, 1)
    add_st_ising2 = FibonacciChain.add_reference_qubits(model, st_ising, 2)
    add_st_ising3 = FibonacciChain.add_reference_qubits(model, st_ising, 3)
    @test add_st_ising1[[collect(1:4)..., collect(13:16)...]] ≈ 1/(2√2)*ones(8)
    @test add_st_ising2[[1, 2, 5, 6, 11, 12, 15, 16]] ≈ 1/(2√2)*ones(8)
    @test add_st_ising3[[collect(1:2:7)..., collect(10:2:16)...]] ≈ 1/(2√2)*ones(8)
end

@testset "reference_boundary_evolution" begin
    model = AnyonModel(FibonacciAnyon(), 4, pbc=true)
    τ = 1000.0
    st = zeros(length(anyon_basis(model))); st[1] = 1

    # First post_selection all zero
    add_st = FibonacciChain.add_reference_qubits(model, st, 1) # still all 00000
    config = MeasureConfig(τ = τ, t₂ = 1, mode=:sample)
    sample = BitVector(zeros(Int, length(2:2:4)))
    boundary_evo = reference_boundary_evolution(model, add_st, config, sample)
    output1, free_energy = boundary_evo.state, boundary_evo.free_energy
    @test length(output1) == 2*length(anyon_basis(model))
    @test free_energy ≈ 1.9248473002384139 atol = 1e-10

    model_Ising = AnyonModel(IsingAnyon(), 4, pbc=true, measure_operator=:X)
    st = zeros(length(anyon_basis(model_Ising))); st[1] = 1

    add_st = FibonacciChain.add_reference_qubits(model_Ising, st, 1)
    sample = BitVector(zeros(Int, length(1:4)))
    boundary_evo_Ising = reference_boundary_evolution(model_Ising, add_st, config, sample)
    output1, free_energy = boundary_evo_Ising.state, boundary_evo_Ising.free_energy
    @test output1[1:16] == 1/4 .* ones(16) # |0++++⟩
    @test free_energy ≈  2.7725887222397816 atol = 1e-10 # log(16) = 4log(2)

    model_Ising = AnyonModel(IsingAnyon(), 2, pbc=true, measure_operator=:ZZ)
    st = zeros(length(anyon_basis(model_Ising))); st[1] = 1

    add_st = FibonacciChain.add_reference_qubits(model_Ising, st, 1)
    sample = BitVector(zeros(Int, length(1:2)))
    boundary_evo_Ising = reference_boundary_evolution(model_Ising, add_st, config, sample, layer_idx=2) # Need to note here must be layer_idx=2 for ZZ measurement, otherwise always measure X
    output1, free_energy = boundary_evo_Ising.state, boundary_evo_Ising.free_energy
    @test output1 == add_st # |000⟩
end

@testset "reference_apply_measurement_layer" begin
    N = 6
    τ = 1000.0
    sign = false
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    k_old = 1
    st = zeros(18); st[1] = 1
    add_st = FibonacciChain.add_reference_qubits(model, st, 1)
    sample = BitVector(zeros(Int, length(2:2:N)))

    basis_F = anyon_basis(model)
    ext_basis = FibonacciChain.build_extended_basis(1, basis_F)

    measure_outcome = FibonacciChain._reference_apply_measurement_layer(model, τ, add_st, sample, 1, extended_basis=ext_basis)
    output1, free_energy = measure_outcome.state, measure_outcome.free_energy
    @test length(output1) == 2*length(basis_F)
    @test free_energy ≈ 2.8872709503576206
end

@testset "reference_sample_layer" begin
    N = 6
    τ = 1000.0
    sign = false
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)
    k_old = 1
    st = zeros(2^N); st[1] = 1 
    add_st = FibonacciChain.add_reference_qubits(model, st, 1)

    basis_F = anyon_basis(model)
    ext_basis = FibonacciChain.build_extended_basis(1, basis_F)
    measure_outcome = FibonacciChain._reference_sample_layer(model, τ, add_st, MersenneTwister(100), 1, extended_basis=ext_basis)
    output1, sample_layer = measure_outcome.state, measure_outcome.sample
    @test length(output1) == 2*length(basis_F)
    @test sample_layer == [1, 1, 1, 0, 0, 1]
end


@testset "reference_bulk_evolution" begin
    N = 8
    τ = 1000.0
    sign = true
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    k_old = 1
    st = zeros(length(anyon_basis(model))); st[1] = 1
    add_st = FibonacciChain.add_reference_qubits_reset(model, st, 1)[3]
    sample = BitMatrix(zeros(Int, (4, length(2:2:N)))) # Four layer sample, 2 peroid evolution

    config = MeasureConfig(τ = τ, t₂ = 2, mode=:sample)
    mo = reference_bulk_evolution(model, add_st, config, sample)
    output1, samples, free_energy = mo.states, mo.samples, mo.free_energys
    @test length(output1) == 2
    @test samples == zeros(Int, (4, length(2:2:N)))

    config = MeasureConfig(τ = τ, t₂ = 2, rng =MersenneTwister(100), mode=:Born)
    mo_born = reference_bulk_evolution(model, add_st, config)
    output1, samples, free_energy = mo_born.states, mo_born.samples, mo_born.free_energys
    @test length(output1) == 2
    @test samples == [1 1 1 1; 0 1 0 1; 0 0 0 0; 0 0 1 0]
    @test free_energy ≈ [2.617994480798357, 2.563763819200177, 1.415477184509482, 2.8031245965479576]
end

@testset "_reference_born_measure" begin
    model = AnyonModel(FibonacciAnyon(), 6)
    t = 10
    measure_config = MeasureConfig(τ=1000.0, t₂=t, rng=MersenneTwister(42), mode=:Born)
    state = zeros(length(anyon_basis(model))); state[1] = 1.0

    add_st = add_reference_qubits(model, state, 1)
    ext_basis = FibonacciChain.build_extended_basis(1, anyon_basis(model))
    mo = FibonacciChain._reference_born_measure(model, add_st, measure_config, extended_basis=ext_basis)
    sts, samples, free_energys = mo.states, mo.samples, mo.free_energys
    @test size(samples) == (20, 3)
    @test free_energys[1] ≈  1.9248473002384137 atol=1e-6
end

@testset "_reference_sample_measure" begin
    N = 6
    model = AnyonModel(IsingAnyon(), N, measure_operator=:X)
    t = 10
    measure_config = MeasureConfig(τ=1000.0, t₂=t, rng=MersenneTwister(42), mode=:sample)
    state = zeros(length(anyon_basis(model))); state[1] = 1.0
    samples = BitMatrix(zeros(2t, N))
    add_st = add_reference_qubits(model, state, 1)

    ext_basis = FibonacciChain.build_extended_basis(1, anyon_basis(model))
    measure_outcome = FibonacciChain._reference_sample_measure(model, add_st, samples, measure_config, extended_basis=ext_basis)
    @test measure_outcome.states[end][[1,64]] ≈ 1/√2 .* ones(2)
    @test measure_outcome.free_energys[end] ≈ 5log(2) atol=1e-6
end

@testset "reference_rdm" begin
    N = 3
    st = ones(4)/2;
    site = 1

    # Add a ref qubit to site 1 with reset 0
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    add_site1 = FibonacciChain.add_reference_qubits_reset(model, st, site)[3]
    
    rdm = reference_rdm(model, [1], add_site1)
    @test rdm == [0.75 0.0; 0.0 0.25]

    # Trace the first ref qubit with IsingAnyon
    model_Ising = AnyonModel(IsingAnyon(), 4, pbc=true, measure_operator=:X)
    full_st = zeros(2^(N+1));
    inds = [1, 3, 5, 10]
    full_st[inds] .= 0.5
    @test anyon_rdm(model_Ising, [1], full_st) == rdm

    # Add another ref qubit to site 1 with reset 0, the 1st rdm should be the same as above, but the 2nd rdm should be |0><0|
    model = AnyonModel(FibonacciAnyon(), N, pbc=true)
    add_st2 = FibonacciChain.add_reference_qubits_reset(model, add_site1, site)[3]
    rdm2 = reference_rdm(model, [1], add_st2)
    rdm3 = reference_rdm(model, [2], add_st2)
    @test rdm ≈ rdm2
    @test rdm3 ≈ [1.0 0.0; 0.0 0.0]
    rdm4 = reference_rdm(model, [1, 2], add_st2)
    @test diag(rdm4) ≈ [0.75, 0.0, 0.25, 0.0]

    # add a ref qubit to a Bell state, the entanglement entropy should not change
    N=4
    model = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:Z)
    st = zeros(2^N); st[1] = 1; st[end] = 1; st ./= norm(st)
    @test ee(anyon_rdm(model, collect(1:div(N,2)), st)) ≈ log(2)
    add_st = FibonacciChain.add_reference_qubits_reset(model, st, 1)[3]
    rdm = reference_rdm(model, collect(1:2), add_st, traceref=false)
    @test ee(rdm) ≈ log(2)
    # Only for GHZ state, the entanglement entropy is not changed after adding reference qubits, if is W state, Haar Random state, NO.
end

@testset "reference_rdm_IsingZ" begin
    N = 3
    st = ones(2^N); st /= norm(st)  # Normalize the state
    site = 1
    model_Ising = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)
    add_site1 = FibonacciChain.add_reference_qubits_reset(model_Ising, st, site)[3]

    # Trace out ref of bell pair, the rdm should be maximally mixed
    rdm = reference_rdm(model_Ising, [1], add_site1)
    @test rdm ≈ [0.5 0.0; 0.0 0.5]
    
    # Trace to get the system, the rdm should be pure Bell state in X basis
    rdm_system = reference_rdm(model_Ising, [2, 3], add_site1, traceref=false)
    @test rdm_system ≈ ones(4,4)/4

    add_st2 = FibonacciChain.add_reference_qubits_reset(model_Ising, add_site1, site)[3]
    rdm2 = reference_rdm(model_Ising, [1], add_st2)
    rdm3 = reference_rdm(model_Ising, [2], add_st2)
    @test rdm ≈ rdm2
    @test rdm3  ≈ [1.0 0.0; 0.0 0.0]
 
    rdm4 = reference_rdm(model_Ising, [1, 2], add_st2)
    @test diag(rdm4) ≈ [0.4999999999999999, 0.0, 0.4999999999999999, 0.0]
end

@testset "reference_rdm_IsingX" begin
    N = 3
    st = zeros(2^N); st[1] = 1; st[end] = 1; st ./= norm(st)  # Bell state
    site = 1
    model_IsingX = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:X)
    add_site1 = FibonacciChain.add_reference_qubits_reset(model_IsingX, st, site)[3]

    # Trace out ref of bell pair, the rdm should be maximally mixed, no matter IsingX or IsingZ
    rdm = reference_rdm(model_IsingX, [1], add_site1)
    @test rdm ≈ [0.5 0.0; 0.0 0.5]
    
    # Trace to get the system, the rdm should be pure Bell state
    rdm_system = reference_rdm(model_IsingX, [2, 3], add_site1, traceref=false)
    @test rdm_system ≈ [0.5 0.0 0.0 0.5; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.5 0.0 0.0 0.5]

    # add_st2 = FibonacciChain.add_reference_qubits(N, add_site1, site, anyon_type = anyon_type, entangle_way = :reset)[3]
    # rdm2 = reference_rdm(N, [1], add_st2, anyon_type = anyon_type)
    # rdm3 = reference_rdm(N, [2], add_st2, anyon_type = anyon_type)
    # @test rdm ≈ rdm2 
    # @test rdm3 ≈ [0.5 0.5; 0.5 0.5]
 
    # rdm4 = reference_rdm(N, [1, 2], add_st2, anyon_type = anyon_type)
    # @test rdm4 ≈ [0.24999999999999994 0.24999999999999994 0.0 0.0; 0.24999999999999994 0.24999999999999994 0.0 0.0; 0.0 0.0 0.24999999999999994 0.24999999999999994; 0.0 0.0 0.24999999999999994 0.24999999999999994]
end

@testset "spatial_corr_matrix and reference rdm" begin
    # The spatial correlation should be the same after adding reference qubits
    N=6
    model_IsingZ = AnyonModel(IsingAnyon(), N, pbc=true, measure_operator=:Z)
    mes = zeros(length(anyon_basis(model_IsingZ)));
    mes[1] = 1/√2 
    mes[end] = 1/√2
    sclis = [spatial_correlation(model_IsingZ, mes, i, j) for i in 1:N for j in 1:N if j!=i]
    @test sclis ≈ log(2) * ones(length(sclis))  # GHZ state

    add_mes = FibonacciChain.add_reference_qubits_reset(model_IsingZ, mes, 1)[3]

    add_mes2 = FibonacciChain.add_reference_qubits_reset(model_IsingZ, add_mes, 1)[3]
    ρ1=reference_rdm(model_IsingZ, [1], add_mes2)
    ρ2=reference_rdm(model_IsingZ, [2], add_mes2)
    ρ12=reference_rdm(model_IsingZ, [1, 2], add_mes2) 
    I = ee(ρ1) + ee(ρ2) - ee(ρ12)
    # Two qubit each form a bell pair, while together is pure state
    @test I ≈ 0.0 atol=1e-12

    H=anyon_ham(model_IsingZ)
    energy, states = eigen(H)
    GS = states[:,1]
    add_GS1 = add_reference_qubits(model_IsingZ, GS, 1)
    add_GS2 = add_reference_qubits(model_IsingZ, add_GS1, 1)
    ρ1=reference_rdm(model_IsingZ, [1], add_GS2)
    ρ2=reference_rdm(model_IsingZ, [2], add_GS2)
    ρ12=reference_rdm(model_IsingZ, [1, 2], add_GS2) 
    I = ee(ρ1) + ee(ρ2) - ee(ρ12)
    @test I ≈ log(2) atol=1e-12

    model_Fibo = AnyonModel(FibonacciAnyon(), N, pbc=true)
    H_Fibo=anyon_ham(model_Fibo)
    energy_Fibo, states_Fibo = eigen(H_Fibo)
    GS_Fibo = states_Fibo[:,1]
    add_GS1_Fibo = add_reference_qubits(model_Fibo, GS_Fibo, 1)
    add_GS2_Fibo = add_reference_qubits(model_Fibo, add_GS1_Fibo, 1)
    ρ1_Fibo=reference_rdm(model_Fibo, [1], add_GS2_Fibo)
    ρ2_Fibo=reference_rdm(model_Fibo, [2], add_GS2_Fibo)
    ρ12_Fibo=reference_rdm(model_Fibo, [1, 2], add_GS2_Fibo)
    I_Fibo = ee(ρ1_Fibo) + ee(ρ2_Fibo) - ee(ρ12_Fibo)
    @test I_Fibo ≈ 0.5895144857350505 atol=1e-12

    # Test for copy way to entangle reference qubits to GHZ state, in contrast, the reset way above give I = 0, so :reset maybe fail. ref to next testset
    add_mes_copy = add_reference_qubits(model_IsingZ, mes, 1)
    # Counting for system, need to add traceref = false, add a ref qubit of Z eigenst does not change the spatial correlation
    ρ = reference_rdm(model_IsingZ, collect(1:N), add_mes_copy, traceref=false)
    sclis_ref = [spatial_correlation(model_IsingZ, ρ, i, j) for i in 1:N for j in 1:N if j!=i]

    @test sclis_ref ≈ sclis

    plus_st = ones(length(anyon_basis(model_IsingZ)));
    minus_st = zeros(length(anyon_basis(model_IsingZ)));

    for i in eachindex(minus_st)
        num_1 = count_ones(anyon_basis(model_IsingZ)[i])
        if iseven(num_1)
            minus_st[i] = 1
        else
            minus_st[i] = -1
        end
    end

    mes = plus_st + minus_st
    mes /= norm(mes)
    sc = spatial_correlation(model_IsingZ, mes, 2, 4)

    add_mes_copy = FibonacciChain.add_reference_qubits(model_IsingZ, mes, 1)
    # Counting for system, need to add traceref = false, add a ref qubit of X eigenst seems change the spatial correlation
    ρ = reference_rdm(model_IsingZ, collect(1:N), add_mes_copy, traceref=false)
    sc_ref = spatial_correlation(model_IsingZ, ρ, 2, 4)
    @test sc_ref ≈ sc
end

@testset "Reference Qubit Correlations" begin
    # ==========================================================
    # Test 4: Reference qubit protocol preserves correlations
    # ==========================================================
    # Adding a reference qubit via :copy should not change the 
    # spatial correlations in the original system
    N_ref = 4
    model_ref = AnyonModel(IsingAnyon(), N_ref, pbc=true, measure_operator=:Z)
    test_state = zeros(2^N_ref)
    test_state[1] = 1/√2
    test_state[end] = 1/√2  # GHZ state
    
    # Correlation before adding reference qubit
    I_before = spatial_correlation(model_ref, test_state, 1, 2)
    
    # Add reference qubit via :copy
    state_with_ref = add_reference_qubits(model_ref, test_state, 1)
    ρ_system = reference_rdm(model_ref, collect(1:N_ref), state_with_ref, traceref=false)
    I_after = spatial_correlation(model_ref, ρ_system, 1, 2)
    
    @test I_before ≈ I_after  atol=1e-10
    @test I_before ≈ log(2) atol=1e-10

    # ==========================================================
    # Test 5: Bell pair with reference - Temporal correlation
    # ==========================================================
    # When adding two reference qubits to a Bell state at the same site,
    # the two reference qubits form a Bell pair themselves
    # Mutual information I(ref1:ref2) = 2log(2)
    N_temp = 4
    model_temp = AnyonModel(IsingAnyon(), N_temp, pbc=true, measure_operator=:Z)
    bell_st = zeros(2^N_temp)
    bell_st[1] = 1/√2
    bell_st[end] = 1/√2  # GHZ state
    
    # Add first reference qubit at site 1
    st_ref1 = add_reference_qubits(model_temp, bell_st, 1)
    # Add second reference qubit at same site
    st_ref2 = add_reference_qubits(model_temp, st_ref1, 1)
    
    ρ1 = reference_rdm(model_temp, [1], st_ref2)
    ρ2 = reference_rdm(model_temp, [2], st_ref2)
    ρ12 = reference_rdm(model_temp, [1, 2], st_ref2)
    I_refs = ee(ρ1) + ee(ρ2) - ee(ρ12)
    
    @test I_refs ≈ log(2) atol=1e-10

    # ==========================================================
    # Test 6: Product state with reference - Zero temporal correlation  
    # ==========================================================
    # For a product state, reference qubits at the same site should
    # have zero mutual information (no temporal correlation)
    N_prod2 = 4
    model_prod2 = AnyonModel(IsingAnyon(), N_prod2, pbc=true, measure_operator=:X)
    prod_st = zeros(2^N_prod2)
    prod_st[1] = 1.0  # |0000⟩ in Z basis = product state in X basis
    
    # Add two reference qubits via :copy at site 1
    st_prod_ref1 = add_reference_qubits(model_prod2, prod_st, 1)
    st_prod_ref2 = add_reference_qubits(model_prod2, st_prod_ref1, 1)
    
    ρ1_prod = reference_rdm(model_prod2, [1], st_prod_ref2)
    ρ2_prod = reference_rdm(model_prod2, [2], st_prod_ref2)
    ρ12_prod = reference_rdm(model_prod2, [1, 2], st_prod_ref2)
    I_refs_prod = ee(ρ1_prod) + ee(ρ2_prod) - ee(ρ12_prod)
    
    # The second reference should be in |0⟩ state (pure), so I = 0
    @test abs(I_refs_prod) < 1e-10

    # ==========================================================
    # Test 8: Maximally mixed state - Zero correlations  
    # ==========================================================
    # A maximally mixed state has I(i:j) = 0 for all pairs
    N_mix = 12
    model_mix = AnyonModel(IsingAnyon(), N_mix, pbc=true, measure_operator=:Z)
    # Create maximally mixed pure state approximation via density matrix
    # We test that a random superposition has bounded correlations
    Random.seed!(42)
    mixed_approx = randn(ComplexF64, 2^N_mix)
    mixed_approx /= norm(mixed_approx)
    
    # For a Haar-random state, the expected correlation is small but not exactly 0
    # Test that the calculation completes without error
    I_rand = spatial_correlation(model_mix, mixed_approx, 1, 2)
    @test I_rand ≈ 0.001089674525010187 atol=1e-2
end

function compute_post_selection_Ising(L::Int64, τ::Float64, D::Int64=5L, δt::Int=2; sign::Bool=false)
    model = AnyonModel(IsingAnyon(), L, pbc=true, measure_operator=:X)
    sample = sign ? BitMatrix(ones(Int, 2D, L)) : BitMatrix(zeros(Int, 2D, L))

    initial_state = ones(length(anyon_basis(model)))
    initial_state /= norm(initial_state) # initial state is all |++++⟩ state

    # Here should use τ to evolve at time slice inserting reference qubit, instead of τ_eff
    config = MeasureConfig(τ = τ, t₂ = D, mode=:sample, enable_τ_eff=false)
    mo = bulk_evolution(model, initial_state, config, sample)
    statelis = mo.states

    ref_sample = sign ? BitMatrix(ones(Int, 2*(D+δt+D), L)) : BitMatrix(zeros(Int, 2*(D+δt+D), L))

    if δt == 0
        ref_config = MeasureConfig(τ = τ, t₂ = D, t₁ = D, mode=:sample, x₂=L÷2+1) # It doesn't need config set here, it will be overwrite in reference_evolution  # seems need to
        mo = reference_evolution(model, statelis, ref_config, ref_sample) # to compute temporal correlation, add ref qubit at site L/2+1
    else
        ref_config = MeasureConfig(τ = τ, t₂ = D+δt, t₁ = D, mode=:sample, x₂=L÷2+1, x₁ = L÷2+1)
        mo = reference_evolution(model, statelis, ref_config, ref_sample) 
    end


    return mo.states[end]
end

function compute_Born_Ising(L::Int64, τ::Float64, D::Int64=5L, δt::Int=2; sign::Bool=false)
    model = AnyonModel(IsingAnyon(), L, pbc=true, measure_operator=:X)
    sample = sign ?  BitMatrix(ones(Int, 2D, L)) : BitMatrix(zeros(Int, 2D, L))

    initial_state = ones(length(anyon_basis(model)))
    initial_state /= norm(initial_state) # initial state is all zero state

    config = MeasureConfig(τ = τ, t₂ = D, mode=:sample, enable_τ_eff=false)
    mo = bulk_evolution(model, initial_state, config, sample)
    statelis = mo.states

    
    if δt == 0
        ref_config = MeasureConfig(τ = τ, t₂ = D, t₁ = D, rng = MersenneTwister(100), mode=:Born, x₂=L÷2+1)
        mo = reference_evolution(model, statelis, ref_config)
        ref2stlis, sample_layer, sample_free_energy = mo.states, mo.samples, mo.free_energys  # to compute temporal correlation, add ref qubit at site L/2+1
    else
        ref_config = MeasureConfig(τ = τ, t₂ = D+δt, t₁ = D, rng = MersenneTwister(100), mode=:Born, x₂=L÷2+1, x₁ = L÷2+1)
        mo = reference_evolution(model, statelis, ref_config)
        ref2stlis, sample_layer, sample_free_energy = mo.states, mo.samples, mo.free_energys  # to compute temporal correlation, add ref qubit at site L/2+1
    end


    return ref2stlis[end], sample_layer, sample_free_energy
end

@testset "reference_evolution" begin
    L = 8
    τ = log(1+√5)/2
    D = 10L
    t = log(1+√2)/π
    model = AnyonModel(IsingAnyon(), L, pbc=true, measure_operator=:X)
    ref2st = compute_post_selection_Ising(L, τ, D, 0, sign=true)
    spatial_corr, _ = ref_correlation(model, ref2st, spatial = true)
    ref2st = compute_post_selection_Ising(L, τ, D, round(Int, t*L), sign=true)
    _, temporal_corr1 = ref_correlation(model, ref2st, temporal = true)
    ref2st = compute_post_selection_Ising(L, τ, D, round(Int, t*L)+1, sign=true)
    _, temporal_corr2 = ref_correlation(model, ref2st, temporal = true)
    @test temporal_corr1/spatial_corr ≈ 1.1397709251465786
    @test temporal_corr2/spatial_corr ≈ 1.050524715257348

    ref2st, sample_layer, sample_free_energy = compute_Born_Ising(L, τ, D, 0, sign=true)
    spatial_corr, _ = ref_correlation(model, ref2st, spatial = true)
    ref2st, sample_layer, sample_free_energy = compute_Born_Ising(L, τ, D, 4, sign=true)
    _, temporal_corr = ref_correlation(model, ref2st, temporal = true)
    @test temporal_corr/spatial_corr ≈ 9.881413517280196
end