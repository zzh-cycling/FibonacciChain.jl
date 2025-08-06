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
    basis = Fibonacci_basis(N)
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
    basis_ising = Fibonacci_basis(N, pbc, measure_class=:IsingX)
    extended_basis_ising = FibonacciChain.build_extended_basis(0, basis_ising)
    @testset extended_basis_ising == basis_ising

    extended_basis_ising1 = FibonacciChain.build_extended_basis(1, basis_ising)
    @test extended_basis_ising1 == Fibonacci_basis(N+1, pbc, measure_class=:IsingX)

    extended_basis_ising2 = FibonacciChain.build_extended_basis(2, basis_ising)
    @test extended_basis_ising2 == Fibonacci_basis(N+2, pbc, measure_class=:IsingX)
end

@testset "add_reference_qubits" begin
    N = 3
    st = ones(4)/2; 
    # Test add_reference_qubit with Proj 0
    seed = 90
    add_st1 = FibonacciChain.add_reference_qubits!(N, st, 1, MersenneTwister(seed))
    add_st2 = FibonacciChain.add_reference_qubits!(N, st, 2, MersenneTwister(seed))
    add_st3 = FibonacciChain.add_reference_qubits!(N, st, 3, MersenneTwister(seed))
    @test add_st1 ≈ [0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5] 
    @test add_st2 ≈ [0.5, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0]
    @test add_st3 ≈ [0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.0]

    add_st = FibonacciChain.add_reference_qubits!(N+1, ones(7)./√7, 1, MersenneTwister(90))
    @test add_st == [0.37796447300922725, 0.37796447300922725, 0.37796447300922725, 0.37796447300922725, 0.37796447300922725, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.37796447300922725, 0.37796447300922725]
    add_st = FibonacciChain.add_reference_qubits!(N+1, ones(7)./√7, 1, MersenneTwister(100))
    @test add_st[[1,3,13,14]] == [0.5, 0.5, 0.5, 0.5]
    add_st = FibonacciChain.add_reference_qubits!(N+1, add_st, 1, MersenneTwister(90))
    @test add_st[[1,3, 20, 21]] == [0.5, 0.5, 0.5, 0.5]

    add2_st1 = FibonacciChain.add_reference_qubits!(N, add_st1, 1, MersenneTwister(seed))
    add2_st2 = FibonacciChain.add_reference_qubits!(N, add_st2, 2, MersenneTwister(seed))
    add2_st3 = FibonacciChain.add_reference_qubits!(N, add_st3, 3, MersenneTwister(seed))

    @test add2_st1[[1,2,3,12]] == [0.5, 0.5, 0.5, 0.5]
    @test add2_st2[[1,2,4,11]] == [0.5, 0.5, 0.5, 0.5]
    @test add2_st3[[1,3,4,10]] == [0.5, 0.5, 0.5, 0.5]

    seed = 100
    add2_st1 = FibonacciChain.add_reference_qubits!(N, add_st1, 1, MersenneTwister(seed))
    add2_st2 = FibonacciChain.add_reference_qubits!(N, add_st2, 2, MersenneTwister(seed))
    add2_st3 = FibonacciChain.add_reference_qubits!(N, add_st3, 3, MersenneTwister(seed))

    @test add2_st1[[5,16]] == [1/√2, 1/√2]
    @test add2_st2[[5,15]] == [1/√2, 1/√2]
    @test add2_st3[[5,14]] == [1/√2, 1/√2]

    add_st1 = FibonacciChain.add_reference_qubits!(N, st, 1, MersenneTwister(seed))
    add_st2 = FibonacciChain.add_reference_qubits!(N, st, 2, MersenneTwister(seed))
    add_st3 = FibonacciChain.add_reference_qubits!(N, st, 3, MersenneTwister(seed))
    @test add_st1[[1,8]] ≈ [1/√2, 1/√2]
    @test add_st2[[1,7]] ≈ [1/√2, 1/√2]
    @test add_st3[[1,6]] ≈ [1/√2, 1/√2]

end

@testset "add_reference_qubits_Ising" begin
    # Test for Ising basis
    N = 3
    seed = 90
    st_ising = zeros(2^N); st_ising[1] = 1/√2; st_ising[end] = 1/√2; # set the last two qubits to be in the Bell state
    add_st_ising1 = FibonacciChain.add_reference_qubits!(N, st_ising, 1, MersenneTwister(seed), measure_class=:IsingX)
    add_st_ising2 = FibonacciChain.add_reference_qubits!(N, st_ising, 2, MersenneTwister(seed), measure_class=:IsingX)
    add_st_ising3 = FibonacciChain.add_reference_qubits!(N, st_ising, 3, MersenneTwister(seed), measure_class=:IsingX)
    @test add_st_ising1[[1,13]] == [1/√2, 1/√2]
    @test add_st_ising2[[1,11]] == [1/√2, 1/√2]
    @test add_st_ising3[[1,10]] == [1/√2, 1/√2]

    seed=100
    add_st_ising4 = FibonacciChain.add_reference_qubits!(N, st_ising, 1, MersenneTwister(seed), measure_class=:IsingX)
    add_st_ising5 = FibonacciChain.add_reference_qubits!(N, st_ising, 2, MersenneTwister(seed), measure_class=:IsingX)
    add_st_ising6 = FibonacciChain.add_reference_qubits!(N, st_ising, 3, MersenneTwister(seed), measure_class=:IsingX)
    @test add_st_ising4[[4,end]] == [1/√2, 1/√2]
    @test add_st_ising5[[6,end]] == [1/√2, 1/√2]
    @test add_st_ising6[[7,end]] == [1/√2, 1/√2]
end

@testset "reference_measure_basismap" begin
    N = 3
    τ = 1.0
    sign = 0
    pbc = true
    k_old = 1
    T = BitStr{N, Int}
    basislis = Fibonacci_basis(N, pbc, measure_class=:Fibo)
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
    measure_class1 = :IsingX
    measure_class2 = :IsingZZ
    basislis = Fibonacci_basis(N, pbc, measure_class=measure_class1)
    l = length(basislis)
    output11 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 1, sign, pbc, k_old=0, measure_class = measure_class1)
    output12 = FibonacciChain.measure_basismap.(T, τ, basislis, 1, sign, pbc, measure_class = measure_class1)
    output21 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 2, sign, pbc, k_old=0, measure_class = measure_class1)
    output22 = FibonacciChain.measure_basismap.(T, τ, basislis, 2, sign, pbc, measure_class = measure_class1)
    output31 = FibonacciChain.reference_measure_basismap.(T, τ, basislis, 3, sign, pbc, k_old=0, measure_class = measure_class1)
    output32 = FibonacciChain.measure_basismap.(T, τ, basislis, 3, sign, pbc, measure_class = measure_class1)
    @test all([all(output11[i] .≈ output12[i]) for i in 1:l])
    @test all([all(output21[i] .≈ output22[i]) for i in 1:l])
    @test all([all(output31[i] .≈ output32[i]) for i in 1:l])

    extended_basis = FibonacciChain.build_extended_basis(1, basislis)
    output13 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 1, sign, pbc, k_old=1, measure_class = measure_class1)
    output23 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 2, sign, pbc, k_old=1, measure_class = measure_class1)
    output33 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 3, sign, pbc, k_old=1, measure_class = measure_class1)

    @test all([all([all(output13[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output23[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:1])
    @test all([all([all(output33[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:1])

    extended_basis2 = FibonacciChain.build_extended_basis(2, basislis)
    output14 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 1, sign, pbc, k_old=2, measure_class = measure_class1)
    output24 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 2, sign, pbc, k_old=2, measure_class = measure_class1)
    output34 = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis2, 3, sign, pbc, k_old=2, measure_class = measure_class1)

    @test all([all([all(output14[i+j*l] .≈ output12[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output24[i+j*l] .≈ output22[i]) for i in 1:l]) for j in 0:3])
    @test all([all([all(output34[i+j*l] .≈ output32[i]) for i in 1:l]) for j in 0:3])

    # Test the IsingZZ
    output12zz = FibonacciChain.measure_basismap.(T, τ, basislis, 1, sign, pbc, measure_class = measure_class2)
    output22zz = FibonacciChain.measure_basismap.(T, τ, basislis, 2, sign, pbc, measure_class = measure_class2)
    output32zz = FibonacciChain.measure_basismap.(T, τ, basislis, 3, sign, pbc, measure_class = measure_class2)

    extended_basis = FibonacciChain.build_extended_basis(1, basislis)
    output13zz = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 1, sign, pbc, k_old=1, measure_class = measure_class2)
    output23zz = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 2, sign, pbc, k_old=1, measure_class = measure_class2)
    output33zz = FibonacciChain.reference_measure_basismap.(T, τ, extended_basis, 3, sign, pbc, k_old=1, measure_class = measure_class2)

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
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1, MersenneTwister(90))

    output13 = FibonacciChain.reference_measuremap(T, τ, add_st, 1, sign, pbc, k_old=1)
    output23 = FibonacciChain.reference_measuremap(T, τ, add_st, 2, sign, pbc, k_old=1)
    output33 = FibonacciChain.reference_measuremap(T, τ, add_st, 3, sign, pbc, k_old=1)
    @test output13 == 0.5*[(1-ϕ^(-1)), 1, 1, -ϕ^(-3/2), -ϕ^(-3/2), 0, 0, ϕ^(-1)]
    @test output23 == 0.5*[(1-ϕ^(-1)- ϕ^(-3/2)), 1, ϕ^(-1)-ϕ^(-3/2), 0, 0 , 0, 0, 1]
    @test output33 == 0.5*[(1-ϕ^(-1)- ϕ^(-3/2)), ϕ^(-1)-ϕ^(-3/2), 1, 0, 0, 0, 0, 1]

    output13 = FibonacciChain.reference_measuremap(T, τ, st, 1, 0, pbc, k_old=0)
    output23 = FibonacciChain.reference_measuremap(T, τ, st, 2, 0, pbc, k_old=0)
    output33 = FibonacciChain.reference_measuremap(T, τ, st, 3, 1, pbc, k_old=0)
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
    measure_class1 = :IsingX
    measure_class2 = :IsingZZ
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1, measure_class = measure_class1)

    output13 = FibonacciChain.reference_measuremap(T, τ, add_st, 1, sign, pbc, k_old=1, measure_class = measure_class1)
    output23 = FibonacciChain.reference_measuremap(T, τ, add_st, 2, sign, pbc, k_old=1, measure_class = measure_class1)
    output33 = FibonacciChain.reference_measuremap(T, τ, add_st, 3, sign, pbc, k_old=1, measure_class = measure_class2)
    @test output13[[1, 5, 9, 13]] == 1/2√2*ones(4)
    @test output23[[1, 3, 13, 15]] == 1/2√2*ones(4)
    @test output33[[1]] == [1/√2]
end

@testset "reference_generate_state" begin
    N = 8
    τ = 1000.0
    sign = 1
    pbc = true
    k_old = 1
    st = zeros(length(Fibonacci_basis(N, pbc))); st[1] = 1
    add_st = FibonacciChain.add_reference_qubits!(N, st, 1)
    sample = zeros(Int, (3, length(2:2:N)))

    output1 = reference_generate_state(τ, add_st, sample, pbc, k_old=1)
    @test length(output1) == 3
end

@testset "reference_rdm" begin
    N = 3
    st = ones(4)/2;
    site = 1
    add_site2 = FibonacciChain.add_reference_qubits!(N, st, site, MersenneTwister(100))

    rdm = reference_rdm(N, [1], add_site2)
    @test rdm ≈ [0.5 0.0; 0.0 0.5]

    add_site1 = FibonacciChain.add_reference_qubits!(N, st, site, MersenneTwister(90))
    
    rdm = reference_rdm(N, [1], add_site1)
    @test rdm == [0.75 0.0; 0.0 0.25]


    full_st = zeros(2^(N+1));
    inds = [1, 3, 5, 10]
    full_st[inds] .= 0.5
    @test rdm_Fibo(4, [1], full_st, measure_class=:IsingX) == rdm

    add_st2 = FibonacciChain.add_reference_qubits!(N, add_site1, site, MersenneTwister(90))
    rdm2 = reference_rdm(N, [1], add_st2)
    rdm3 = reference_rdm(N, [2], add_st2)
    @test rdm == rdm2
    @test rdm2 + 0.25 * [0 1;1 0] == rdm3

    rdm4 = reference_rdm(N, [1, 2], add_st2)
    @test rdm4 == [0.5 0.25 0.0; 0.25 0.25 0.0; 0.0 0.0 0.25]
end

@testset "reference_rdm_Ising" begin
    N = 3
    st = ones(2^N); st /= norm(st)  # Normalize the state
    site = 1
    measure_class = :IsingX
    add_site1 = FibonacciChain.add_reference_qubits!(N, st, site, measure_class = measure_class)
    
    rdm = reference_rdm(N, [1], add_site1, measure_class = measure_class)
    rdm_system = reference_rdm(N, [2, 3], add_site1, measure_class = measure_class)
    @test rdm ≈ [0.5 0.0; 0.0 0.5]
  

    add_st2 = FibonacciChain.add_reference_qubits!(N, add_site1, site, measure_class = measure_class)
    rdm2 = reference_rdm(N, [1], add_st2, measure_class = measure_class)
    rdm3 = reference_rdm(N, [2], add_st2, measure_class = measure_class)
    @test rdm == rdm2
    @test rdm2 == rdm3

    rdm4 = reference_rdm(N, [1, 2], add_st2, measure_class = measure_class)
    @test rdm4 ≈ [0.4999999999999999 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.4999999999999999]
end

@testset "spatial_corr_matrix and reference rdm" begin
    # The spatial correlation should be the same after adding reference qubits
    N=6
    mes = zeros(length(Fibonacci_basis(N, true)))
    mes[13] = 1/√2 
    mes[end] = 1/√2
    sclis = [spatial_correlation(N, mes, i, j) for i in 1:N for j in 1:N if j!=i]

    add_mes = FibonacciChain.add_reference_qubits!(N, mes, 1)
    # Counting for system
    ρ = reference_rdm(N, collect(1:N), add_mes, traceref=false)
    sclis_ref = [spatial_correlation(N, ρ, i, j) for i in 1:N for j in 1:N if j!=i]

    @test sclis ≈ sclis_ref
end