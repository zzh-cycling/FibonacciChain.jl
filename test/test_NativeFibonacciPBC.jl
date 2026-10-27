using FibonacciChain, Test, Random, LinearAlgebra
import TensorKit, MPSKit
using TensorKit: ⊗

function native_test_matrix(op,sites)
    t = convert(TensorKit.TensorMap,op)
    c = TensorKit.FibonacciAnyon(sites.flux == :trivial ? :I : :τ)
    return Matrix(TensorKit.block(t,c))
end

@testset "Canonical flux PBC Hamiltonians" begin
    expected4 = Dict(:trivial=>[-3.236067977499789,-0.7639320225002105],
        :tau=>[-3.0,-2.642684597778774,-1.0,-1.0,-0.35731540222122404])
    φ = (1+sqrt(5.0))/2
    for N in (4,6)
        model = AnyonModel(FibonacciAnyon(),N)
        encodedH = anyon_ham(model)
        yeig = eigen(Hermitian(topological_charge_operator(model)))
        allenergies = Float64[]
        for flux in (:trivial,:tau)
            sites = native_pbc_sites(N;flux)
            H = native_pbc_dense(model;flux,construction=:mpo)
            reference = native_pbc_dense(model;flux,construction=:braid)
            @test H ≈ H' atol=2e-11
            @test eigvals(Hermitian(H)) ≈ eigvals(Hermitian(reference)) atol=2e-11
            # The two braiding prescriptions in the supplied comment are
            # opposite complex conventions in τ flux, not identical matrices.
            @test H ≈ (flux == :trivial ? reference : conj(reference)) atol=2e-11
            N == 4 && (@test eigvals(Hermitian(H)) ≈ expected4[flux] atol=2e-11)
            y = flux == :trivial ? φ : -1/φ
            U = yeig.vectors[:,findall(x->isapprox(x,y;atol=1e-10),yeig.values)]
            @test size(H,1) == size(U,2)
            @test eigvals(Hermitian(H)) ≈ eigvals(Hermitian(U'*encodedH*U)) atol=2e-11
            append!(allenergies,eigvals(Hermitian(H)))
            projectors = [native_test_matrix(measurement_operator_mps(model,sites,i,Inf,true),sites) for i in 1:N]
            @test -sum(projectors) ≈ H atol=2e-11
            for i in 1:N
                p = projectors[i]
                @test p*p ≈ p atol=2e-11
                @test p*projectors[mod1(i+1,N)]*p ≈ p/φ^2 atol=2e-11
            end
            ferro = AnyonModel(FibonacciAnyon(),N;measure_operator=:Ferro)
            @test native_pbc_dense(ferro;flux) ≈ -H atol=2e-11
        end
        @test sort(allenergies) ≈ eigvals(Hermitian(encodedH)) atol=2e-11
    end
    # Dropping the winding around the τ flux leg changes the Hamiltonian.
    sites = native_pbc_sites(4;flux=:tau)
    p,_ = fibonacci_category_projectors()
    h = -p
    wrong = MPSKit.FiniteMPOHamiltonian(sites.lattice,
        [(i,i+1)=>h for i in 1:3]...,(1,4)=>h)
    @test norm(eigvals(Hermitian(native_test_matrix(wrong,sites))) - expected4[:tau]) > 1e-3
    @test_throws ArgumentError native_pbc_sites(4;flux=:invalid)
    @test_throws ArgumentError anyon_ham(AnyonModel(FibonacciAnyon(),4;pbc=false),sites)
    @test_throws ArgumentError native_pbc_dense(AnyonModel(FibonacciAnyon(),14))
end

@testset "Canonical native states and DMRG" begin
    model = AnyonModel(FibonacciAnyon(),4)
    for flux in (:trivial,:tau)
        sites = native_pbc_sites(4;flux)
        ψ = native_pbc_state(4;flux,seed=31)
        @test length(ψ.mps) == 4+Int(flux == :tau)
        @test MPSKit.left_virtualspace(ψ.mps,1) == TensorKit.unitspace(sites[1])
        @test MPSKit.right_virtualspace(ψ.mps,length(ψ.mps)) ==
            (flux == :trivial ? TensorKit.unitspace(sites[1]) : sites[1])
        v = native_pbc_amplitudes(ψ)
        @test norm(v) ≈ 1 atol=2e-12
        @test norm(ψ) ≈ 1 atol=2e-12
        @test native_pbc_amplitudes(native_pbc_state(4;flux,seed=31)) ≈ v atol=2e-12
        @test native_pbc_amplitudes(native_pbc_state(v,4;flux)) ≈ v atol=3e-12
        @test native_pbc_amplitudes(native_pbc_state(1e-30*v,4;flux))/1e-30 ≈ v atol=3e-12
        @test dot(ψ,(0.2+0.3im)*ψ) ≈ 0.2+0.3im atol=3e-12
        H = native_pbc_dense(model;flux)
        @test native_pbc_energy(model,ψ) ≈ real(dot(v,H*v)) atol=3e-12
        ground,energy = anyon_mps_gst(model;backend=:native_pbc,flux,
            sweep_times=12,maxdim=32,cutoff=1e-12,outputlevel=0,seed=51)
        @test energy ≈ eigmin(Hermitian(H)) atol=2e-9
        @test norm(H*native_pbc_amplitudes(ground)-energy*native_pbc_amplitudes(ground)) < 1e-7
        @test native_pbc_energy(model,ground) ≈ energy atol=2e-11
        @test isfinite(native_schmidt_entropy(ground,2))
        @test_throws ArgumentError ee_mps(ground,2)
        @test_throws ArgumentError mps_statevector(ground)
        @test_throws DimensionMismatch native_pbc_state(ones(3),4;flux)
        initial,returned_sites = initial_mps(4;backend=:native_pbc,flux)
        @test initial.sites.flux == returned_sites.flux == flux
    end
    @test dot(native_pbc_state(4;flux=:trivial),native_pbc_state(4;flux=:tau)) == 0
    @test_throws ArgumentError initial_mps(4;backend=:mpskit,flux=:tau)
    @test_throws ArgumentError native_pbc_state(zeros(2),4)
end

@testset "Native flux PBC measurement trajectories" begin
    model = AnyonModel(FibonacciAnyon(),4)
    for flux in (:trivial,:tau)
        ψ = native_pbc_state(4;flux,seed=92)
        sites = ψ.sites
        v = native_pbc_amplitudes(ψ)
        for i in 1:4, τ in (0.0,0.8,Inf)
            probabilities = Float64[]
            matrices = Matrix{ComplexF64}[]
            for bit in (false,true)
                op = measurement_operator_mps(model,sites,i,τ,bit)
                matrix = native_test_matrix(op,sites)
                push!(matrices,matrix)
                expected = matrix*v
                out,prob = measuremap(model,ψ,sites,i,τ,bit;cutoff=1e-14,maxdim=32)
                push!(probabilities,prob)
                @test prob ≈ norm(expected)^2 atol=3e-12
                @test native_pbc_amplitudes(out)*sqrt(prob) ≈ expected atol=4e-12
            end
            @test sum(probabilities) ≈ 1 atol=4e-12
            @test sum(m'*m for m in matrices) ≈ I atol=4e-12
        end
        for stride in (1,2)
            config(mode) = MeasureConfig(;τ=0.8,t₂=2,mode,maxdim=32,cutoff=1e-14,
                truncate_every_events=stride,track_y_expectation=true,rng=MersenneTwister(32))
            out = bulk_evolution(model,sites,ψ,config(:Born))
            expected = copy(v)
            energies = Float64[]
            for row in 1:4
                strength = row == 4 ? 0.4 : 0.8
                positions = isodd(row) ? (2,4) : (1,3)
                energy = 0.0
                for (k,i) in enumerate(positions)
                    matrix = native_test_matrix(measurement_operator_mps(model,sites,i,strength,out.samples[row,k]),sites)
                    expected = matrix*expected
                    energy -= log(norm(expected)^2)
                    normalize!(expected)
                end
                push!(energies,energy)
            end
            @test native_pbc_amplitudes(out.state) ≈ expected atol=1e-11
            @test out.free_energys ≈ energies atol=2e-6
            @test out.y_expectation_values ≈ fill(flux == :trivial ? (1+sqrt(5))/2 : -2/(1+sqrt(5)),2) atol=2e-7
            replay = bulk_evolution(model,sites,ψ,config(:sample),out.samples)
            @test native_pbc_amplitudes(replay.state) ≈ expected atol=1e-11
            raw = bulk_evolution(model,sites,ψ,config(:sample),out.samples;normalized=false)
            @test native_pbc_amplitudes(raw.state) ≈ exp(-sum(energies)/2)*expected atol=1e-11
            tiny = bulk_evolution(model,sites,1e-30*ψ,config(:sample),out.samples;normalized=false)
            @test native_pbc_amplitudes(tiny.state)/1e-30 ≈ native_pbc_amplitudes(raw.state) atol=1e-12
        end
        @test native_pbc_amplitudes(ψ) ≈ v atol=2e-12
        @test_throws ArgumentError bulk_evolution(model,sites,ψ,MeasureConfig(τ=0.8,t₂=1,mode=:Born);compute_entropy=true)
        out = bulk_evolution(model,sites,ψ,MeasureConfig(τ=0.8,t₂=1,mode=:Born,rng=MersenneTwister(82));
            compute_entropy=true,entropy=:native_schmidt)
        @test all(isfinite,out.entanglement_entropys)
    end
end
