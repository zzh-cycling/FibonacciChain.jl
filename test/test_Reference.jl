using FibonacciChain
using Test
using BitBasis
using LinearAlgebra 
using Random

@testset "build_extended_basis" begin
    N = 3
    pbc = true
    k_total = 2
    T = BitStr{N+k_total, Int}
    basis = anyon_basis(N)
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
    basis_ising = anyon_basis(N, pbc, anyon_type=:IsingX)
    extended_basis_ising = FibonacciChain.build_extended_basis(0, basis_ising)
    @testset extended_basis_ising == basis_ising

    extended_basis_ising1 = FibonacciChain.build_extended_basis(1, basis_ising)
    @test extended_basis_ising1 == anyon_basis(N+1, pbc, anyon_type=:IsingX)

    extended_basis_ising2 = FibonacciChain.build_extended_basis(2, basis_ising)
    @test extended_basis_ising2 == anyon_basis(N+2, pbc, anyon_type=:IsingX)
end

@testset "reference_measure_basismap" begin
    N = 3
    τ = 1.0
    sign = 0
    pbc = true
    k_old = 1
    T = BitStr{N, Int}
    basislis = anyon_basis(N, pbc, anyon_type=:Fibo)
    l = length(basislis)
    output11 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 1, sign, pbc, k_old=0)
    output12 = FibonacciChain.measure_basismap.(T, τ, basislis, 1, sign, pbc)
    output21 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 2, sign, pbc, k_old=0)
    output22 = FibonacciChain.measure_basismap.(T, τ, basislis, 2, sign, pbc)
    output31 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 3, sign, pbc, k_old=0)
    output32 = FibonacciChain.measure_basismap.(T, τ, basislis, 3, sign, pbc)
    @test all([all(output11[i] .≈ output12[i]) for i in 1:l])
    @test all([all(output21[i] .≈ output22[i]) for i in 1:l])
    @test all([all(output31[i] .≈ output32[i]) for i in 1:l])

    extended_basis = FibonacciChain.build_extended_basis(1, basislis)
    output13 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 1, sign, pbc, k_old=1)
    output23 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 2, sign, pbc, k_old=1)
    output33 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 3, sign, pbc, k_old=1)

    @test all([all([all(output13[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output23[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output33[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:1])

    extended_basis2 = FibonacciChain.build_extended_basis(2, basislis)
    output14 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 1, sign, pbc, k_old=2)
    output24 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 2, sign, pbc, k_old=2)
    output34 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 3, sign, pbc, k_old=2)

    @test all([all([all(output14[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output24[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output34[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:3])
end

@testset "reference_measure_basismap_Ising" begin
    N = 3
    τ = 1000.0
    sign = 0
    pbc = true
    k_old = 1
    T = BitStr{N, Int}
    anyon_type1 = :IsingX
    anyon_type2 = :IsingZZ
    basislis = anyon_basis(N, pbc, anyon_type=anyon_type1)
    l = length(basislis)
    output11 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 1, sign, pbc, k_old=0, anyon_type = anyon_type1)
    output12 = FibonacciChain.measure_basismap.(T, τ, basislis, 1, sign, pbc, anyon_type = anyon_type1)
    output21 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 2, sign, pbc, k_old=0, anyon_type = anyon_type1)
    output22 = FibonacciChain.measure_basismap.(T, τ, basislis, 2, sign, pbc, anyon_type = anyon_type1)
    output31 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 3, sign, pbc, k_old=0, anyon_type = anyon_type1)
    output32 = FibonacciChain.measure_basismap.(T, τ, basislis, 3, sign, pbc, anyon_type = anyon_type1)
    @test all([all(output11[i] .≈ output12[i]) for i in 1:l])
    @test all([all(output21[i] .≈ output22[i]) for i in 1:l])
    @test all([all(output31[i] .≈ output32[i]) for i in 1:l])

    extended_basis = FibonacciChain.build_extended_basis(1, basislis)
    output13 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 1, sign, pbc, k_old=1, anyon_type = anyon_type1)
    output23 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 2, sign, pbc, k_old=1, anyon_type = anyon_type1)
    output33 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 3, sign, pbc, k_old=1, anyon_type = anyon_type1)

    @test all([all([all(output13[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output23[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output33[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:1])

    extended_basis2 = FibonacciChain.build_extended_basis(2, basislis)
    output14 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 1, sign, pbc, k_old=2, anyon_type = anyon_type1)
    output24 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 2, sign, pbc, k_old=2, anyon_type = anyon_type1)
    output34 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 3, sign, pbc, k_old=2, anyon_type = anyon_type1)

    @test all([all([all(output14[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output24[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output34[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:3])

    # Test the IsingZZ
    output12zz = FibonacciChain.measure_basismap.(T, τ, basislis, 1, sign, pbc, anyon_type = anyon_type2)
    output22zz = FibonacciChain.measure_basismap.(T, τ, basislis, 2, sign, pbc, anyon_type = anyon_type2)
    output32zz = FibonacciChain.measure_basismap.(T, τ, basislis, 3, sign, pbc, anyon_type = anyon_type2)

    extended_basis = FibonacciChain.build_extended_basis(1, basislis)
    output13zz = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 1, sign, pbc, k_old=1, anyon_type = anyon_type2)
    output23zz = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 2, sign, pbc, k_old=1, anyon_type = anyon_type2)
    output33zz = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 3, sign, pbc, k_old=1, anyon_type = anyon_type2)

    @test all([all([all(output13zz[i+j*l] .≈ output12zz[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output23zz[i+j*l] .≈ output22zz[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output33zz[i+j*l] .≈ output32zz[i]) for i in 1:l]) for j in 0:1])
end

@testset "reference_measuremap" begin
    N = 3
    τ = 1000.0
    sign = 0
    pbc = true
    k_old = 1
    T = BitStr{N, Int}
    st = ones(4)/2;
    ϕ = (1 + √5) / 2  
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1, entangle_way = :reset)[3]

    ext_basis = FibonacciChain.build_extended_basis(1, anyon_basis(N, pbc))
    output13 = FibonacciChain.reference_measuremap(T, τ, add_st, 1, sign, pbc, k_old=1, extended_basis=ext_basis)
    output23 = FibonacciChain.reference_measuremap(T, τ, add_st, 2, sign, pbc, k_old=1, extended_basis=ext_basis)
    output33 = FibonacciChain.reference_measuremap(T, τ, add_st, 3, sign, pbc, k_old=1, extended_basis=ext_basis)
    @test output13 == 0.5*[(1-ϕ^(-1)), 1, 1, -ϕ^(-3/2), -ϕ^(-3/2), 0, 0, ϕ^(-1)]
    @test output23 == 0.5*[(1-ϕ^(-1)- ϕ^(-3/2)), 1, ϕ^(-1)-ϕ^(-3/2), 0, 0 , 0, 0, 1]
    @test output33 == 0.5*[(1-ϕ^(-1)- ϕ^(-3/2)), ϕ^(-1)-ϕ^(-3/2), 1, 0, 0, 0, 0, 1]

    output13 = FibonacciChain.reference_measuremap(T, τ, st, 1, 0, pbc, k_old=0, extended_basis=anyon_basis(N, pbc))
    output23 = FibonacciChain.reference_measuremap(T, τ, st, 2, 0, pbc, k_old=0, extended_basis=anyon_basis(N, pbc))
    output33 = FibonacciChain.reference_measuremap(T, τ, st, 3, 1, pbc, k_old=0, extended_basis=anyon_basis(N, pbc))
    @test output13 == measuremap(N, τ, st, 1, 0, pbc)
    @test output23 == measuremap(N, τ, st, 2, 0, pbc)
    @test output33 == measuremap(N, τ, st, 3, 1, pbc)
end

@testset "reference_measuremap_Ising" begin
    N = 3
    τ = 1000.0
    sign = 0
    pbc = true
    k_old = 1
    T = BitStr{N, Int}
    st = zeros(2^N); st[1] = 1 
    anyon_type1 = :IsingX
    anyon_type2 = :IsingZZ
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1, anyon_type = anyon_type1, entangle_way = :reset)[3]

    ext_basis = FibonacciChain.build_extended_basis(1, anyon_basis(N, pbc, anyon_type=anyon_type1))
    output13 = FibonacciChain.reference_measuremap(T, τ, add_st, 1, sign, pbc, k_old=1, anyon_type = anyon_type1, extended_basis=ext_basis)
    output23 = FibonacciChain.reference_measuremap(T, τ, add_st, 2, sign, pbc, k_old=1, anyon_type = anyon_type1, extended_basis=ext_basis)
    output33 = FibonacciChain.reference_measuremap(T, τ, add_st, 3, sign, pbc, k_old=1, anyon_type = anyon_type2, extended_basis=ext_basis)
    @test output13[[1, 5, 9, 13]] ≈ 1/2√2*ones(4)
    @test output23[[1, 3, 13, 15]] ≈ 1/2√2*ones(4)
    @test output33[[1]] ≈ [1/√2]
end

@testset "concat_bell_pair" begin
    N = 3
    st = ones(4)/2; 
    basis_F = anyon_basis(N)
    ext_basis = FibonacciChain.build_extended_basis(1, basis_F)
    bell_st = FibonacciChain.concat_bell_pair!(N, st, basis_F, ext_basis, site_idx = 1, k_old = 0)
    @test bell_st[[1,2,3,8]] == 1/2*ones(4)
end

@testset "add_reference_qubits" begin
    N = 3
    st = ones(4)/2; 
    entangle_way = :reset
    # Test add_reference_qubit reset to 0
    add_st = FibonacciChain.add_reference_qubits!(N+1, ones(7)./√7, 1, entangle_way = entangle_way)[3]
    @test add_st[[collect(1:5)..., 13,14]] ≈ 0.37796447300922725*ones(7)

    add_stp = FibonacciChain.add_reference_qubits!(N+1, add_st, 1, entangle_way = entangle_way)[3]
    @test add_stp[[collect(1:5)..., 20,21]] ≈ 0.37796447300922725*ones(7)

    add_st1 = FibonacciChain.add_reference_qubits!(N, st, 1, entangle_way = entangle_way)[3]
    add_st2 = FibonacciChain.add_reference_qubits!(N, st, 2, entangle_way = entangle_way)[3]
    add_st3 = FibonacciChain.add_reference_qubits!(N, st, 3, entangle_way = entangle_way)[3]
    @test add_st1 ≈ [0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5] 
    @test add_st2 ≈ [0.5, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0]
    @test add_st3 ≈ [0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.0]

    add2_st1 = FibonacciChain.add_reference_qubits!(N, add_st1, 1, entangle_way = entangle_way)[3]
    add2_st2 = FibonacciChain.add_reference_qubits!(N, add_st2, 2, entangle_way = entangle_way)[3]
    add2_st3 = FibonacciChain.add_reference_qubits!(N, add_st3, 3, entangle_way = entangle_way)[3]

    @test add2_st1[[1,2,3,12]] == [0.5, 0.5, 0.5, 0.5]
    @test add2_st2[[1,2,4,11]] == [0.5, 0.5, 0.5, 0.5]
    @test add2_st3[[1,3,4,10]] == [0.5, 0.5, 0.5, 0.5]
end

@testset "add_reference_qubits_Ising" begin
    # Test for Ising basis
    N = 3
    st_ising = zeros(2^N); st_ising[1] = 1/√2; st_ising[end] = 1/√2; # set the GHZ state.
    add_st_ising1 = FibonacciChain.add_reference_qubits!(N, st_ising, 1, anyon_type=:IsingZ, entangle_way = :reset)[3]
    add_st_ising2 = FibonacciChain.add_reference_qubits!(N, st_ising, 2, anyon_type=:IsingZ, entangle_way = :reset)[3]
    add_st_ising3 = FibonacciChain.add_reference_qubits!(N, st_ising, 3, anyon_type=:IsingZ, entangle_way = :reset)[3]
    @test add_st_ising1[[1,13]] == [1/√2, 1/√2]
    @test add_st_ising2[[1,11]] == [1/√2, 1/√2]
    @test add_st_ising3[[1,10]] == [1/√2, 1/√2]

    st_ising = zeros(2^N); st_ising[1] = 1/√2; st_ising[end] = 1/√2; # set the last two qubits to be in the Bell state
    add_st_ising1 = FibonacciChain.add_reference_qubits!(N, st_ising, 1, anyon_type=:IsingX, entangle_way = :reset)[3]
    add_st_ising2 = FibonacciChain.add_reference_qubits!(N, st_ising, 2, anyon_type=:IsingX, entangle_way = :reset)[3]
    add_st_ising3 = FibonacciChain.add_reference_qubits!(N, st_ising, 3, anyon_type=:IsingX, entangle_way = :reset)[3]
    @test add_st_ising1[[1, 4, 13, 16]] == 1/2*ones(4)
    @test add_st_ising2[[1, 6, 11, 16]] == 1/2*ones(4)
    @test add_st_ising3[[1, 7, 10, 16]] == 1/2*ones(4)
end

@testset "add_reference_qubits_Ising_copy" begin
    # Test for Ising basis
    N = 3
    st_ising = zeros(2^N); st_ising[1] = 1/√2; st_ising[end] = 1/√2; # set the GHZ state.
    add_st_ising1 = FibonacciChain.add_reference_qubits!(N, st_ising, 1, anyon_type=:IsingZ) # For copy way to entangle, the result is irrelevant to measurement outcome.
    add_st_ising2 = FibonacciChain.add_reference_qubits!(N, st_ising, 2, anyon_type=:IsingZ)
    add_st_ising3 = FibonacciChain.add_reference_qubits!(N, st_ising, 3, anyon_type=:IsingZ)
    @test add_st_ising1 == add_st_ising2 == add_st_ising3
    @test add_st_ising1[[1,16]] ≈ [1/√2, 1/√2]

    st_ising = ones(2^N); st_ising /= norm(st_ising); # set to be plus state
    add_st_ising1 = FibonacciChain.add_reference_qubits!(N, st_ising, 1, anyon_type=:IsingX)
    add_st_ising2 = FibonacciChain.add_reference_qubits!(N, st_ising, 2, anyon_type=:IsingX)
    add_st_ising3 = FibonacciChain.add_reference_qubits!(N, st_ising, 3, anyon_type=:IsingX)
    @test add_st_ising1[[collect(1:4)..., collect(13:16)...]] ≈ 1/(2√2)*ones(8)
    @test add_st_ising2[[1, 2, 5, 6, 11, 12, 15, 16]] ≈ 1/(2√2)*ones(8)
    @test add_st_ising3[[collect(1:2:7)..., collect(10:2:16)...]] ≈ 1/(2√2)*ones(8)
end

@testset "reference_apply_measurement_layer!" begin
    N = 6
    τ = 1000.0
    sign = 0
    pbc = true
    k_old = 1
    st = zeros(18); st[1] = 1
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1, entangle_way = :copy)
    sample = zeros(Int, length(2:2:N))

    basis_F = anyon_basis(N, pbc)
    ext_basis = FibonacciChain.build_extended_basis(1, basis_F)
    output1, free_energy = reference_apply_measurement_layer!(N, τ, add_st, sample, 1, pbc, extended_basis=ext_basis, return_free_energy=true)
    @test length(output1) == 2*length(basis_F)
    @test free_energy ≈ 2.8872709503576206
end

@testset "reference_sample_layer" begin
    N = 6
    τ = 1000.0
    sign = 0
    pbc = true
    k_old = 1
    st = zeros(2^N); st[1] = 1 
    anyon_type = :IsingX
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1, anyon_type = anyon_type, entangle_way = :copy)
    sample = zeros(Int, length(2:2:N))

    basis_F = anyon_basis(N, pbc, anyon_type = anyon_type)
    ext_basis = FibonacciChain.build_extended_basis(1, basis_F)
    output1, sample_layer = reference_sample_layer!(N, τ, add_st, MersenneTwister(100), 1, pbc, anyon_type = anyon_type, extended_basis=ext_basis)
    @test length(output1) == 2*length(basis_F)
    @test sample_layer == [1, 1, 1, 0, 0, 1]
end


@testset "reference_generate_state" begin
    N = 8
    τ = 1000.0
    sign = 1
    pbc = true
    k_old = 1
    st = zeros(length(anyon_basis(N, pbc))); st[1] = 1
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1, entangle_way = :reset)[3]
    sample = zeros(Int, (3, length(2:2:N)))

    output1, samples = reference_generate_state(τ, add_st, sample, pbc)
    @test length(output1) == 3
    @test samples == zeros(Int, (3, length(2:2:N)))

    output1, samples = reference_generate_state(τ, add_st, sample, pbc, rng = MersenneTwister(100), mode=:Born)
    @test length(output1) == 3
    @test samples == [1 1 1 1; 0 1 0 1; 0 0 0 0]
end

@testset "reference_rdm" begin
    N = 3
    st = ones(4)/2;
    site = 1

    # Add a ref qubit to site 1 with reset 0
    add_site1 = FibonacciChain.add_reference_qubits!(N, st, site, entangle_way = :reset)[3]
    
    rdm = reference_rdm(N, [1], add_site1)
    @test rdm == [0.75 0.0; 0.0 0.25]

    # Trace the first ref qubit with anyon_type = :IsingX
    full_st = zeros(2^(N+1));
    inds = [1, 3, 5, 10]
    full_st[inds] .= 0.5
    @test anyon_rdm(4, [1], full_st, anyon_type=:IsingX) == rdm

    # Add another ref qubit to site 1 with reset 0, the 1st rdm should be the same as above, but the 2nd rdm should be |0><0|
    add_st2 = FibonacciChain.add_reference_qubits!(N, add_site1, site, entangle_way = :reset)[3]
    rdm2 = reference_rdm(N, [1], add_st2)
    rdm3 = reference_rdm(N, [2], add_st2)
    @test rdm ≈ rdm2
    @test rdm3 ≈ [1.0 0.0; 0.0 0.0]

    rdm4 = reference_rdm(N, [1, 2], add_st2)
    @test diag(rdm4) ≈ [0.75, 0.0, 0.25, 0.0]

    # add a ref qubit to a Bell state, the entanglement entropy should not change
    N=4
    st = zeros(2^N); st[1] = 1; st[end] = 1; st ./= norm(st)
    @test ee(anyon_rdm(N, collect(1:div(N,2)), st, anyon_type=:IsingZ)) ≈ log(2)
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1, anyon_type=:IsingZ, entangle_way = :reset)[3]
    rdm = reference_rdm(N, collect(1:2), add_st, anyon_type=:IsingZ, traceref=false)
    @test ee(rdm) ≈ log(2)
    # Only for GHZ state, the entanglement entropy is not changed after adding reference qubits, if is W state, Haar Random state, NO.
end

@testset "reference_rdm_IsingZ" begin
    N = 3
    st = ones(2^N); st /= norm(st)  # Normalize the state
    site = 1
    anyon_type = :IsingZ
    add_site1 = FibonacciChain.add_reference_qubits!(N, st, site, anyon_type = anyon_type, entangle_way = :reset)[3]

    # Trace out ref of bell pair, the rdm should be maximally mixed
    rdm = reference_rdm(N, [1], add_site1, anyon_type = anyon_type)
    @test rdm ≈ [0.5 0.0; 0.0 0.5]
    
    # Trace to get the system, the rdm should be pure Bell state in X basis
    rdm_system = reference_rdm(N, [2, 3], add_site1, anyon_type = anyon_type, traceref=false)
    @test rdm_system ≈ ones(4,4)/4

    add_st2 = FibonacciChain.add_reference_qubits!(N, add_site1, site, anyon_type = anyon_type, entangle_way = :reset)[3]
    rdm2 = reference_rdm(N, [1], add_st2, anyon_type = anyon_type)
    rdm3 = reference_rdm(N, [2], add_st2, anyon_type = anyon_type)
    @test rdm == rdm2
    @test rdm3  ≈ [1.0 0.0; 0.0 0.0]
 
    rdm4 = reference_rdm(N, [1, 2], add_st2, anyon_type = anyon_type)
    @test diag(rdm4) ≈ [0.4999999999999999, 0.0, 0.4999999999999999, 0.0]
end

@testset "reference_rdm_IsingX" begin
    N = 3
    st = zeros(2^N); st[1] = 1; st[end] = 1; st ./= norm(st)  # Bell state
    site = 1
    anyon_type = :IsingX
    add_site1 = FibonacciChain.add_reference_qubits!(N, st, site, anyon_type = anyon_type, entangle_way = :reset)[3]

    # Trace out ref of bell pair, the rdm should be maximally mixed, no matter IsingX or IsingZ
    rdm = reference_rdm(N, [1], add_site1, anyon_type = anyon_type)
    @test rdm ≈ [0.5 0.0; 0.0 0.5]
    
    # Trace to get the system, the rdm should be pure Bell state
    rdm_system = reference_rdm(N, [2, 3], add_site1, anyon_type = anyon_type, traceref=false)
    @test rdm_system ≈ [0.5 0.0 0.0 0.5; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.5 0.0 0.0 0.5]

    # add_st2 = FibonacciChain.add_reference_qubits!(N, add_site1, site, anyon_type = anyon_type, entangle_way = :reset)[3]
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
    anyon_type = :IsingZ
    mes = zeros(length(anyon_basis(N, anyon_type=anyon_type)));
    mes[1] = 1/√2 
    mes[end] = 1/√2
    sclis = [spatial_correlation(N, mes, i, j, anyon_type = anyon_type) for i in 1:N for j in 1:N if j!=i]
    @test sclis ≈ log(2) * ones(length(sclis))  

    add_mes = FibonacciChain.add_reference_qubits!(N, mes, 1, anyon_type = anyon_type, entangle_way = :reset)[3]

    add_mes2 = FibonacciChain.add_reference_qubits!(N, add_mes, 1, anyon_type = anyon_type, entangle_way = :reset)[3]
    ρ1=reference_rdm(N, [1], add_mes, anyon_type = anyon_type)
    ρ2=reference_rdm(N, [2], add_mes, anyon_type = anyon_type)
    ρ12=reference_rdm(N, [1, 2], add_mes, anyon_type = anyon_type) 
    I = ee(ρ1) + ee(ρ2) - ee(ρ12)
    # Two qubit each form a bell pair, while together is pure state
    @test I ≈ 2*log(2) atol=1e-12

    add_mes_copy = FibonacciChain.add_reference_qubits!(N, mes, 1, anyon_type = anyon_type, entangle_way = :copy)
    # Counting for system, need to add traceref = false, add a ref qubit of Z eigenst does not change the spatial correlation
    ρ = reference_rdm(N, collect(1:N), add_mes_copy, anyon_type = anyon_type, traceref=false)
    sclis_ref = [spatial_correlation(N, ρ, i, j, anyon_type = anyon_type) for i in 1:N for j in 1:N if j!=i]

    @test sclis_ref ≈ sclis

    plus_st = ones(length(anyon_basis(N, anyon_type=anyon_type)));
    minus_st = zeros(length(anyon_basis(N, anyon_type=anyon_type)));

    for i in eachindex(minus_st)
        num_1 = count_ones(anyon_basis(N, anyon_type=anyon_type)[i])
        if iseven(num_1)
            minus_st[i] = 1
        else
            minus_st[i] = -1
        end
    end

    mes = plus_st + minus_st
    mes /= norm(mes)
    sc = spatial_correlation(N, mes, 2, 4, anyon_type = anyon_type) 

    add_mes_copy = FibonacciChain.add_reference_qubits!(N, mes, 1, anyon_type = anyon_type, entangle_way = :copy)
    # Counting for system, need to add traceref = false, add a ref qubit of X eigenst seems change the spatial correlation
    ρ = reference_rdm(N, collect(1:N), add_mes_copy, anyon_type = anyon_type, traceref=false)
    sc_ref = spatial_correlation(N, ρ, 2, 4, anyon_type = anyon_type)
    @test sc_ref ≈ sc
end

function compute_post_selection_Ising(L::Int64, τ::Float64, D::Int64=5L, δt::Int=2; sign::Int64=0, entangle_way::Symbol=:copy, rng = MersenneTwister(100))
    pbc = true
    anyon_type = :IsingX
    sample = (sign == 1) ? ones(Int, D, L) : zeros(Int, D, L)

    initial_state = ones(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
    initial_state /= norm(initial_state) # initial state is all zero state

    statelis = generate_state(τ, initial_state, sample, temp= true, anyon_type=anyon_type)

    ref_sample = (sign == 0) ? zeros(Int, D+δt+D+1, L) : ones(Int, D+δt+D+1, L)
    
    if entangle_way == :copy
        if δt == 0
            ref2st, sample_layer = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D, anyon_type=:IsingX, return_free_energy=false, rng=rng) # to compute temporal correlation, add ref qubit at site L/2+1

        else
            ref2st, sample_layer = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D+δt, anyon_type=:IsingX, x₁ = L÷2+1, return_free_energy=false, rng=rng,) # to compute temporal correlation, add ref qubit at site L/2+1
        end
    end

    return ref2st
end

function compute_Born_Ising(L::Int64, τ::Float64, D::Int64=5L, δt::Int=2; sign::Int64=0, entangle_way::Symbol=:copy, mode::Symbol=:Born, rng = MersenneTwister(100))
    pbc = true
    anyon_type = :IsingX
    sample = (sign == 1) ? ones(Int, D, L) : zeros(Int, D, L)

    initial_state = ones(length(anyon_basis(BitStr{L, Int}, pbc, anyon_type=anyon_type)))
    initial_state /= norm(initial_state) # initial state is all zero state

    statelis = generate_state(τ, initial_state, sample, temp= true, anyon_type=anyon_type)

    ref_sample = (sign == 0) ? zeros(Int, D+δt+D+1, L) : ones(Int, D+δt+D+1, L)
    
    if entangle_way == :copy
        if δt == 0
            ref2st, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D, anyon_type=:IsingX, return_free_energy=true, rng=rng, mode=mode) # to compute temporal correlation, add ref qubit at site L/2+1

        else
            ref2st, sample_layer, sample_free_energy = reference_evolution(L, τ, statelis, ref_sample, L÷2+1, D, D+δt, anyon_type=:IsingX, x₁ = L÷2+1, return_free_energy=true, rng=rng, mode=mode) # to compute temporal correlation, add ref qubit at site L/2+1
        end
    end

    return ref2st, sample_layer, sample_free_energy
end

@testset "reference_evolution" begin
    L = 8
    τ = log(1+√5)/2
    D = 10L
    ref2st = compute_post_selection_Ising(L, τ, D, 0, sign=1, entangle_way=:copy)
    spatial_corr, _ = ref_correlation(L, ref2st, anyon_type=:IsingX, spatial = true)
    ref2st = compute_post_selection_Ising(L, τ, D, 4, sign=1, entangle_way=:copy)
    _, temporal_corr = ref_correlation(L, ref2st, anyon_type=:IsingX, temporal = true)
    @test temporal_corr/spatial_corr ≈ 1.3473583771849356

    ref2st, sample_layer, sample_free_energy = compute_Born_Ising(L, τ, D, 0, sign=1, entangle_way=:copy, rng = MersenneTwister(100))
    spatial_corr, _ = ref_correlation(L, ref2st, anyon_type=:IsingX, spatial = true)
    ref2st, sample_layer, sample_free_energy = compute_Born_Ising(L, τ, D, 4, sign=1, entangle_way=:copy, rng = MersenneTwister(100))
    _, temporal_corr = ref_correlation(L, ref2st, anyon_type=:IsingX, temporal = true)
    @test temporal_corr/spatial_corr ≈ 6.69459188579601
end