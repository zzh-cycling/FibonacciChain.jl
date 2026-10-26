using FibonacciChain, Test, Random, LinearAlgebra
import MPSKit, TensorKit
using ITensors, ITensorMPS

# Independent representations, common BitBasis ordering for comparisons.
function kit_test_itvector(ψ, sites)
    a = Array(prod(ψ), sites...)
    return vec(permutedims(a, length(sites):-1:1))
end

function kit_test_itstate(v, N)
    sites = siteinds("Qubit", N)
    a = permutedims(reshape(v, ntuple(_ -> 2, N)), N:-1:1)
    return MPS(ITensor(a, sites...), sites; cutoff = 0.0), sites
end

function kit_test_matrix(op, N)
    a = convert(Array, convert(TensorKit.TensorMap, op))
    axes = (reverse(1:N)..., reverse(N+1:2N)...)
    return reshape(permutedims(reshape(a, ntuple(_ -> 2, 2N)), axes), 2^N, 2^N)
end

@testset "MPSKit encoded representation and validation" begin
    N = 4
    rng = MersenneTwister(153)
    v = normalize!(randn(rng, ComplexF64, 2^N))
    ψ = mpskit_state(v, N)
    @test mps_statevector(ψ) ≈ v atol = 2e-14
    it, isites = kit_test_itstate(v, N)
    @test kit_test_itvector(it, isites) ≈ v atol = 2e-14
    @test mps_statevector(mpskit_state(it)) ≈ v atol = 2e-14
    @test first(initial_mps(N)) isa MPSKit.FiniteMPS
    projected = FibonacciChain._kit_compress(
        FibonacciChain.fibonacci_constraint_projector_mpo(mpskit_sites(N)) * ψ;
        cutoff=1e-14, maxdim=32)
    pv = mps_statevector(projected)
    for coefficient in (-1.7, 0.2im, 1e-30)
        scaled = FibonacciChain._kit_scale(projected, coefficient)
        combined = FibonacciChain._kit_add(ψ, scaled)
        @test mps_statevector(combined) ≈ v + coefficient*pv atol=3e-13
        compressed = FibonacciChain._kit_compress(scaled; cutoff=1e-14, maxdim=32)
        @test mps_statevector(compressed)/coefficient ≈ pv atol=3e-13
    end
    for b in 1:N-1
        @test ee_mps(ψ, b) ≈ ee_mps(it, b) atol = 2e-12
    end
    for constructor in (initial_mps, evenparity_mps)
        k, _ = constructor(N; backend = :mpskit)
        i, s = constructor(N; backend = :itensor)
        @test mps_statevector(k) ≈ kit_test_itvector(i, s) atol = 1e-14
    end
    @test_throws ArgumentError initial_mps(N; backend = :unknown)
    @test_throws DimensionMismatch mpskit_state(v[1:end-1], N)
    model = AnyonModel(FibonacciAnyon(), N)
    k, s = initial_mps(N; backend = :mpskit)
    @test_throws DimensionMismatch bulk_evolution(model, s, k,
        MeasureConfig(; τ = 0.8, t₂ = 1), falses(1, 2))
    @test_throws ArgumentError bulk_evolution(model, s, k,
        MeasureConfig(; τ = 0.8, t₂ = 1, mode = :Born); normalized = false)
end

@testset "MPSKit operators against ITensor and exact basis" begin
    N = 4
    ks = mpskit_sites(N)
    its = siteinds("Qubit", N)
    v = normalize!(randn(MersenneTwister(52), ComplexF64, 2^N))
    k = mpskit_state(v, N)
    it, its = kit_test_itstate(v, N)
    models = AnyonModel[AnyonModel(FibonacciAnyon(), N; measure_operator = convention, pbc)
        for convention in (:Antiferro, :Ferro) for pbc in (false, true)]
    append!(models, [AnyonModel(SpinHalf(), N; model_type = kind, pbc, λ = 0.3)
        for kind in (:Ising, :OBF) for pbc in (false, true)])
    for model in models
        h = anyon_ham(model, ks)
        hi = anyon_ham(model, its)
        @test MPSKit.expectation_value(k, h) ≈ inner(prime(it), hi, it) atol = 2e-12
        raw = Array(prod(hi), prime.(its)..., its...)
        expected = reshape(permutedims(raw, (reverse(1:N)..., reverse(N+1:2N)...)), 2^N, 2^N)
        @test kit_test_matrix(h, N) ≈ expected atol = 2e-12
        basis = anyon_basis(model)
        indices = [Int(b.buf)+1 for b in basis]
        # Legacy dense OBF OBC construction calls OBFmap beyond N-2. The full
        # independent ITensor MPO matrix above is the oracle for that case.
        if !(model isa AnyonModel{SpinHalf,:OBF} && !model.pbc)
            @test kit_test_matrix(h, N)[indices, indices] ≈ anyon_ham(model) atol = 2e-12
        end
    end
    for convention in (:Antiferro, :Ferro), i in 1:N, τ in (0.0, 0.8, Inf), bit in (false, true)
        model = AnyonModel(FibonacciAnyon(), N; measure_operator = convention)
        q = measurement_operator_mps(model, ks, i, τ, bit)
        qi = FibonacciChain._measurement_operator_mps_application(model, its, i, τ, bit)
        expected = apply(qi, it; cutoff = 0.0)
        @test mps_statevector(q*k) ≈ kit_test_itvector(expected, its) atol = 3e-12
    end
    for pbc in (false, true)
        p = FibonacciChain.fibonacci_constraint_projector_mpo(ks; pbc)
        pi = FibonacciChain.fibonacci_constraint_projector_mpo(its; pbc)
        @test mps_statevector(p*k) ≈ kit_test_itvector(apply(pi, it; cutoff=0.0), its) atol=3e-12
        y = topological_charge_mpo(ks; pbc)
        yi = topological_charge_mpo(its; pbc)
        @test mps_statevector(y*k) ≈ kit_test_itvector(apply(yi, it; cutoff=0.0), its) atol=3e-12
        @test fibonacci_constraint_violation(k; pbc) ≈ fibonacci_constraint_violation(it; pbc) atol=3e-12
        @test fibonacci_constraint_violation(0.2*k; pbc) ≈ fibonacci_constraint_violation(k; pbc) atol=3e-12
    end
end

@testset "MPSKit Born, replay and unnormalized trajectories" begin
    for kind in (:Fibonacci, :Ising, :OBF), stride in (1, 2)
        N = 6
        model = kind == :Fibonacci ? AnyonModel(FibonacciAnyon(), N) :
            AnyonModel(SpinHalf(), N; model_type = kind, λ = 0.3)
        k, ks = initial_mps(N; backend = :mpskit)
        it, its = initial_mps(N; backend = :itensor)
        config(mode) = MeasureConfig(; τ = 0.8, t₂ = 2, mode, rng = MersenneTwister(36),
            cutoff = 1e-14, maxdim = 32, truncate_every_events = stride)
        a = bulk_evolution(model, ks, k, config(:Born))
        b = bulk_evolution(model, its, it, config(:Born))
        @test a.samples == b.samples
        @test a.free_energys ≈ b.free_energys atol = 3e-5
        @test a.entanglement_entropys ≈ b.entanglement_entropys atol = 2e-6
        @test abs(dot(mps_statevector(a.state), kit_test_itvector(b.state, its))) ≈ 1 atol=2e-11
        replay = bulk_evolution(model, ks, k, config(:sample), a.samples)
        @test mps_statevector(replay.state) ≈ mps_statevector(a.state) atol=2e-11
        @test replay.free_energys ≈ a.free_energys atol=3e-5
        ar = bulk_evolution(model, ks, k, config(:sample), a.samples; normalized=false)
        br = bulk_evolution(model, its, it, config(:sample), a.samples; normalized=false)
        @test mps_statevector(ar.state) ≈ kit_test_itvector(br.state, its) atol=2e-11
        scaled = bulk_evolution(model, ks, 1e-30*k, config(:sample), a.samples;
            normalized=false, compute_entropy=false)
        @test mps_statevector(scaled.state) / 1e-30 ≈ mps_statevector(ar.state) rtol=2e-10 atol=1e-15
        @test mps_statevector(k) ≈ [1; zeros(2^N-1)] atol=1e-14
    end
end

@testset "MPSKit hybrid schedules and charge tracking" begin
    N = 6
    model = AnyonModel(FibonacciAnyon(), N)
    k, ks = initial_mps(N; backend=:mpskit)
    it, its = initial_mps(N; backend=:itensor)
    config(mode) = HybridConfig(; τ=0.8, t₂=2, mode, p=0.5, random_angles=true,
        rng=MersenneTwister(7), maxdim=32, cutoff=1e-14, track_y_expectation=true,
        enforce_fibonacci_constraint=true)
    a = bulk_evolution(model, ks, k, config(:Born))
    b = bulk_evolution(model, its, it, config(:Born))
    @test a.schedule.measurement_mask == b.schedule.measurement_mask
    @test a.schedule.outcomes == b.schedule.outcomes
    @test isequal(a.schedule.unitary_angles, b.schedule.unitary_angles)
    @test a.free_energys ≈ b.free_energys atol=3e-6
    @test a.entanglement_entropys ≈ b.entanglement_entropys atol=3e-6
    @test a.y_expectation_values ≈ b.y_expectation_values atol=3e-6
    @test abs(dot(mps_statevector(a.state), kit_test_itvector(b.state, its))) ≈ 1 atol=3e-11
    replay = bulk_evolution(model, ks, k, config(:sample), a.schedule)
    @test mps_statevector(replay.state) ≈ mps_statevector(a.state) atol=3e-11
    @test fibonacci_constraint_violation(a.state) < 1e-12
end

@testset "MPSKit DMRG against ED and ITensor" begin
    for kind in (:Fibonacci, :Ising)
        model = kind == :Fibonacci ? AnyonModel(FibonacciAnyon(), 6) :
            AnyonModel(SpinHalf(), 6; model_type=:Ising)
        k, ek = anyon_mps_gst(model; backend=:mpskit, sweep_times=12, maxdim=32, cutoff=1e-12, outputlevel=0)
        it, ei = anyon_mps_gst(model; backend=:itensor, sweep_times=12, maxdim=32, cutoff=1e-12, outputlevel=0)
        ed = eigmin(Hermitian(anyon_ham(model)))
        @test ek ≈ ed atol=2e-8
        @test ek ≈ ei atol=2e-8
        @test ee_mps(k, 3) ≈ ee_mps(it, 3) atol=2e-6
        kind == :Fibonacci && (@test fibonacci_constraint_violation(k) < 1e-8)
    end
end

@testset "Native Fibonacci category foundations" begin
    p, q = fibonacci_category_projectors()
    φ = (1+sqrt(5.0))/2
    @test norm(p*p-p) < 1e-14
    @test norm(q*q-q) < 1e-14
    @test norm(p*q) < 1e-14
    @test norm(p+q-TensorKit.id(TensorKit.domain(p))) < 1e-14
    @test tr(p) ≈ 1
    @test tr(q) ≈ φ
    τsector = TensorKit.FibonacciAnyon(:τ)
    sectors = (TensorKit.FibonacciAnyon(:I), τsector)
    F = [TensorKit.Fsymbol(τsector,τsector,τsector,τsector,a,b) for a in sectors, b in sectors]
    @test F'*F ≈ I atol=1e-14
    swap = [0.0 1.0; 1.0 0.0]
    vacuum = swap*F*Diagonal([1.0,0.0])*F'*swap
    Q = (1-2/φ)*Diagonal([1.0,-1.0])-2φ^(-3/2)*swap
    @test vacuum ≈ (I-Q)/2 atol=1e-14
    for θ in (0.0, 0.7, π)
        u = p + cis(θ)*q
        @test norm(u'*u-p-q) < 1e-14
    end
    for τ in (0.0, 0.8, Inf)
        a,b = FibonacciChain._kit_measurement_coefficients(τ,false)
        k0 = (a-b)*p+(a+b)*q
        k1 = (a+b)*p+(a-b)*q
        @test norm(k0'*k0+k1'*k1-p-q) < 1e-14
    end
end

@testset "MPSKit Lyapunov frame and charge ancilla" begin
    N = 4
    model = AnyonModel(FibonacciAnyon(), N)
    k, ks = initial_mps(N; backend=:mpskit)
    it, its = initial_mps(N; backend=:itensor)
    record = BitMatrix([0 1; 1 0; 0 0; 1 1])
    v = normalize!(randn(MersenneTwister(72), ComplexF64, 2^N))
    random_k, random_it = mpskit_state(v, N), first(kit_test_itstate(v,N))
    for i in (1, N)
        ak, _ = add_reference_qubits(model, random_k, ks, i)
        ai, asi = add_reference_qubits(model, random_it, siteinds(random_it), i)
        @test mps_statevector(ak) ≈ kit_test_itvector(ai, asi) atol=3e-12
    end
    lk = lyapunov_spectrum_mps(model, ks, 0.4, record; n_states=3, maxdim=32, cutoff=1e-14)
    li = lyapunov_spectrum_mps(model, its, 0.4, record; n_states=3, maxdim=32, cutoff=1e-14)
    @test lk ≈ li atol=2e-10
    for sector in (:trivial, :tau)
        lk = lyapunov_spectrum_mps(model, ks, 0.4, record; n_states=1, sector, maxdim=32, cutoff=1e-14)
        li = lyapunov_spectrum_mps(model, its, 0.4, record; n_states=1, sector, maxdim=32, cutoff=1e-14)
        @test lk ≈ li atol=2e-10
    end
    config(mode) = MeasureConfig(; τ=0.8, t₂=2, mode, maxdim=32, cutoff=1e-14,
        rng=MersenneTwister(32), enforce_fibonacci_constraint=true)
    a = topological_charge_sharpening(model, ks, k, config(:Born))
    b = topological_charge_sharpening(model, its, it, config(:Born))
    @test a.samples == b.samples
    @test a.free_energys ≈ b.free_energys atol=2e-6
    @test a.entanglement_entropys ≈ b.entanglement_entropys atol=2e-6
    @test reference_rdm(model, a.state) ≈ reference_rdm(model, b.state) atol=2e-11
    r = topological_charge_sharpening(model, ks, k, config(:sample), a.samples)
    @test reference_rdm(model, r.state) ≈ reference_rdm(model, a.state) atol=2e-11
    sk, rk, pk = mps_measurement_enumeration(model, k, ks, [1, 3], 0.8; cutoff=1e-14, maxdim=32)
    si, ri, pi = mps_measurement_enumeration(model, it, its, [1, 3], 0.8; cutoff=1e-14, maxdim=32)
    @test rk == ri
    @test pk ≈ pi atol=2e-12
    @test sum(pk) ≈ 1 atol=2e-12
end

@testset "MPSKit finite bond convergence" begin
    N = 6
    model = AnyonModel(FibonacciAnyon(), N)
    k, ks = initial_mps(N; backend=:mpskit)
    record = BitMatrix([0 1 0; 1 0 1; 0 0 1; 1 1 0])
    results = [bulk_evolution(model, ks, k, MeasureConfig(; τ=0.8, t₂=2,
        maxdim=χ, cutoff=0.0), record).state for χ in (2, 4, 8)]
    v = mps_statevector(results[end])
    errors = [max(0, 1-abs(dot(v,mps_statevector(ψ)))^2) for ψ in results]
    @test errors[2] <= errors[1] + 1e-12
    @test errors[end] < 1e-12
    for (χ, ψ) in zip((2,4,8), results)
        @test maximum(TensorKit.dim(MPSKit.right_virtualspace(ψ,i)) for i in 1:N-1) <= χ
        @test norm(ψ) ≈ 1 atol=1e-12
    end
end
