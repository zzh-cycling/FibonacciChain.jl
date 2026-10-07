using Test, Random, LinearAlgebra
include(joinpath(@__DIR__, "..", "trajectory_equivalence.jl"))

@testset "Deformed PEPS and trajectory equivalence" begin
    psi = fibonacci_stringnet_peps(1,1)
    ref = exact_stringnet_reference(psi)
    ne = length(psi.links)
    edges = collect(1:ne)
    vac,tau = Tuple(zeros(Int,ne)),Tuple(ones(Int,ne))
    phi = (1+sqrt(5))/2
    base = stringnet_target_distribution(ref)
    @test base.probabilities[vac] ≈ 1/(1+phi^2)
    @test base.probabilities[tau] ≈ phi^2/(1+phi^2)
    @test length(base.probabilities) == 2
    @test_throws ArgumentError exact_stringnet_reference(psi;max_elements=2)

    @testset "Local filters versus exact reweighting" begin
        for bc in (:open,:cylinder_x)
            small = bc == :open ? psi : fibonacci_stringnet_peps(2,1;bc)
            reference = bc == :open ? ref : exact_stringnet_reference(small)
            for j in (-Inf,-0.43,0.0,0.7)
                es = [1,2]
                filtered = deform_stringnet_peps(small;edges=es,J=j)
                direct = exact_stringnet_reference(filtered.peps)
                target = stringnet_target_distribution(reference;edges=es,J=j)
                @test direct.input_norm2*exp(2filtered.log_scale)/reference.input_norm2 ≈ target.deformation_norm_ratio
                for (key,a) in target.amplitudes
                    idx = 1+sum(Int(key[e]) << (e-1) for e in eachindex(key))
                    @test direct.amplitudes[idx] ≈ a atol=1e-12
                end
            end
        end
        projected = stringnet_target_distribution(ref;edges,J=-Inf)
        @test projected.p_zigzag_tau ≈ base.probabilities[tau]
        @test !isapprox(projected.p_zigzag_tau,base.probabilities[tau]^ne)
        @test projected.probabilities[tau] ≈ 1
        @test length(projected.probabilities) == 1
        rawref = exact_stringnet_reference(fibonacci_stringnet_peps(1,1;normalize=false))
        @test rawref.input_norm2 ≈ 1+phi^2
        @test stringnet_target_distribution(rawref;edges,J=-Inf).p_zigzag_tau ≈ projected.p_zigzag_tau
        huge = stringnet_target_distribution(ref;edges=[1],J=-1000)
        @test isfinite(huge.log_probabilities[vac])
        @test huge.probabilities[vac] == 0
        underflow = compare_trajectories(ref,[collect(vac)];edges=[1],J=-1000,bootstrap_replicates=0)
        @test underflow.summary.empirical.kl_empirical_to_target ≈ -huge.log_probabilities[vac]
    end

    @testset "Marginals, conditioning, and encoding" begin
        marginal = TrajectoryEdgeMap(psi,reshape([1,2],1,2);outcome_labels=(1,0))
        target = stringnet_target_distribution(ref;mapping=marginal)
        @test target.probabilities[(0,0)] ≈ base.probabilities[vac]
        @test target.amplitudes === nothing
        @test length(target.unresolved_edges) == ne-2
        stats = compare_trajectories(ref,[trues(1,2)];mapping=marginal,bootstrap_replicates=0)
        @test only(filter(r -> r.count==1,stats.rows)).edge_labels == (0,0)
        conditioned = TrajectoryEdgeMap(psi,[1];fixed_edges=Dict(e=>1 for e in 2:ne))
        ct = stringnet_target_distribution(ref;mapping=conditioned)
        @test ct.probabilities[(1,)] ≈ 1
        @test ct.fixed_event_probability ≈ base.probabilities[tau]
        @test ct.amplitudes !== nothing
        @test_throws ArgumentError compare_exact_amplitudes(ref,Dict((0,0)=>1);mapping=marginal)
        @test_throws ArgumentError compare_trajectories(ref,[trues(1,2)];mapping=marginal,
            log_probabilities=[log(base.probabilities[vac])],phases=[1],bootstrap_replicates=0)
        @test_throws DimensionMismatch compare_trajectories(ref,[ones(2)];mapping=marginal)
        @test_throws ArgumentError TrajectoryEdgeMap(psi,[1,1])
        @test_throws ArgumentError TrajectoryEdgeMap(psi,[1];fixed_edges=Dict(1=>0))
        @test_throws ArgumentError deform_stringnet_peps(psi;edges=[ne+1],J=0)
        @test_throws ArgumentError deform_stringnet_peps(psi;edges=[1],J=Inf)
    end

    @testset "Sampling statistics and phases" begin
        # This tension makes the two single-hexagon configurations equiprobable.
        j = log(phi)/ne
        records = vcat([collect(vac) for _ in 1:100],[collect(tau) for _ in 1:100])
        result = compare_trajectories(ref,records;edges,J=j,log_probabilities=fill(-log(2),200),
            phases=ones(200),rng=MersenneTwister(42),bootstrap_replicates=39)
        @test result.summary.empirical.total_variation < 1e-12
        @test result.summary.tv_null_pvalue == 1
        @test result.summary.importance.bhattacharyya ≈ 1
        @test result.summary.importance.classical_fidelity_unbiased ≈ 1
        @test result.summary.importance.quantum.fidelity_unbiased ≈ 1
        @test result.summary.importance.importance_ess ≈ 200
        @test result.summary.importance.sampled_target_mass ≈ 1
        unequal = compare_trajectories(ref,records;log_probabilities=fill(-log(2),200),bootstrap_replicates=0)
        expected_b = (sqrt(base.probabilities[vac])+sqrt(base.probabilities[tau]))/sqrt(2)
        @test unequal.summary.importance.bhattacharyya ≈ expected_b
        @test unequal.summary.importance.sampled_target_mass ≈ 1
        @test unequal.summary.importance.kl_p_to_q ≈ -log(2)-log(base.probabilities[vac]*base.probabilities[tau])/2
        restricted = compare_trajectories(ref,records[1:100];log_probabilities=zeros(100),bootstrap_replicates=0)
        @test restricted.summary.importance.sampled_target_mass ≈ base.probabilities[vac]
        no_phase = compare_trajectories(ref,records;edges,J=j,bootstrap_replicates=0)
        @test no_phase.summary.importance === nothing
        correlated = compare_trajectories(ref,records;edges,J=j,iid=false,log_probabilities=fill(-log(2),200))
        @test correlated.summary.tv_null_pvalue === nothing
        @test correlated.summary.importance.bhattacharyya_standard_error === nothing
        @test correlated.summary.importance.classical_fidelity_unbiased === nothing
        wrong = compare_trajectories(ref,records[1:100];edges,J=j,rng=MersenneTwister(42),bootstrap_replicates=99)
        @test wrong.summary.tv_null_pvalue <= 0.01
        @test wrong.summary.unseen_target_mass ≈ 0.5
        unsupported = compare_trajectories(ref,[[1,0,0,0,0,0]];bootstrap_replicates=1)
        @test unsupported.summary.forbidden_target_records == 1
        @test unsupported.summary.branching_violations == 1
        @test unsupported.summary.tv_null_pvalue == 0
        @test unsupported.summary.empirical.kl_empirical_to_target == Inf
        badlogp = fill(-log(2),200); badlogp[2] = -1
        @test_throws ArgumentError compare_trajectories(ref,records;log_probabilities=badlogp,bootstrap_replicates=0)
        target = stringnet_target_distribution(ref;edges,J=j)
        phased = Dict(k=>a*exp(0.31im) for (k,a) in target.amplitudes)
        exact = compare_exact_amplitudes(ref,phased;edges,J=j)
        @test exact.fidelity ≈ 1
        @test exact.max_amplitude_error < 1e-12
        flipped = copy(target.amplitudes); flipped[vac] *= -1
        exact_flipped = compare_exact_amplitudes(ref,flipped;edges,J=j)
        @test exact_flipped.classical_fidelity ≈ 1
        @test exact_flipped.fidelity < 1e-25
        estimated_flipped = compare_trajectories(ref,records;edges,J=j,
            log_probabilities=fill(-log(2),200),phases=vcat(fill(-1,100),ones(100)),bootstrap_replicates=0)
        @test estimated_flipped.summary.importance.quantum.fidelity_plugin < 1e-25
        fit = fit_string_tension(ref,records;edges,J_grid=[-0.5,0,j,0.5])
        @test fit.best_J == j
        @test fit.identifiable
        pinned = TrajectoryEdgeMap(psi,[1];fixed_edges=Dict(e=>1 for e in 2:ne))
        flat = fit_string_tension(ref,[[1]];mapping=pinned,edges=[2],J_grid=[-Inf,-1,0])
        @test !flat.identifiable
        @test flat.best_J === nothing
        vacuum_boundary = TrajectoryEdgeMap(psi,[1];fixed_edges=Dict(e=>0 for e in 2:ne))
        excluded = fit_string_tension(ref,[[0]];mapping=vacuum_boundary,edges=[2],J_grid=[-Inf,0])
        @test excluded.log_likelihood[1] == -Inf
        @test excluded.best_J == 0
    end

    @testset "CFT trajectory schema" begin
        data = Dict("sample"=>Bool[0 1;1 0;0 0;1 1],"sample_free_energy"=>Float32[1,2,3,4],
            "L"=>4,"periods"=>2,"tau"=>2.0,"seed"=>15,"backend"=>"exact","initial_state"=>"TCI_GS")
        record = cft_trajectory_record(data;terminal_half_strength=true)
        @test record.log_probability == -10
        @test record.outcome_labels == (1,0)
        @test record.measurement_sites == [2 4;1 3;2 4;1 3]
        @test record.measurement_strengths == [2,2,2,1]
        @test record.metadata.seed == 15
        @test record.phases === nothing
        @test cft_trajectory_record(data).measurement_strengths === nothing
        @test cft_trajectory_record(data;layers=1:2).log_probability == -3
        @test cft_trajectory_record(data;layers=2:3).log_probability === nothing
        @test cft_trajectory_record(data;columns=[1]).log_probability === nothing
        @test cft_trajectory_record(data;layers=[2,1]).log_probability === nothing
        @test_throws ArgumentError cft_trajectory_record(data;layers=[1,1])
        bad = copy(data); bad["sample_free_energy"] = [1]
        @test_throws DimensionMismatch cft_trajectory_record(bad)
        firstevent = cft_trajectory_record(data;layers=[1],columns=[1])
        mapping = TrajectoryEdgeMap(psi,reshape([1],1,1);outcome_labels=record.outcome_labels)
        compared = compare_trajectories(ref,[firstevent];mapping,bootstrap_replicates=0)
        @test only(filter(r -> r.count==1,compared.rows)).edge_labels == (1,)
    end
end

include(joinpath(@__DIR__, "..", "cylinder_transfer.jl"))
@testset "Periodic F network, boundary sectors, and coherent records" begin
    for L in (4,6,8)
        ref = fibonacci_cylinder_reference(L;tau=atanh(0.95))
        model = FibonacciChain.AnyonModel(FibonacciChain.FibonacciAnyon(),L;pbc=true)
        @test ref.Y ≈ FibonacciChain.topological_charge_operator(model) atol=1e-12
        @test ref.Y*ref.Y ≈ I+ref.Y atol=1e-12
        ps = fibonacci_flux_projectors(ref.Y)
        @test ps.vacuum_flux+ps.tau_flux ≈ I atol=1e-12
        @test ps.vacuum_flux*ps.tau_flux ≈ zero(ref.Y) atol=1e-12
        @test ps.vacuum_flux^2 ≈ ps.vacuum_flux atol=1e-12
        @test cylinder_sector_weights(ref).vacuum_flux ≈ 1 atol=1e-12
        for i in 1:L
            P = Matrix(ref.vacuum_projectors[i])
            @test P ≈ FibonacciChain.measure_matrix(model,1000.0,i,true) atol=1e-12
            @test P^2 ≈ P atol=1e-12
            @test P*ref.Y ≈ ref.Y*P atol=1e-12
            for t in (0.0,0.3,2.0,Inf), sign in (false,true)
                f = fibonacci_outcome_filter(t)
                hi,lo = f[1,1],f[1,2]
                K = sign ? lo*I+(hi-lo)*P : hi*I+(lo-hi)*P
                @test K ≈ FibonacciChain.measure_matrix(model,t,i,sign) atol=1e-12
            end
        end
    end
    @test_throws ArgumentError fibonacci_cylinder_reference(40;tau=1)
    @test_throws ArgumentError fibonacci_cylinder_reference(5;tau=1)
    @test_throws ArgumentError fibonacci_outcome_filter(-1)
    ref = fibonacci_cylinder_reference(4;tau=0.7)
    other = fibonacci_cylinder_reference(4;tau=0.7,sector=:tau_flux)
    @test cylinder_sector_weights(other).tau_flux ≈ 1 atol=1e-12
    for input in (ref,other)
        state = cylinder_record_state(input,2)
        @test state.norm2 ≈ 1 atol=1e-12
        @test 0 < state.record_purity < 1-1e-5
        for k in eachindex(state.probabilities)
            physical = reshape([_sn_bit(k-1,e) for e in 1:4],2,2)
            raw = physical .== 0
            path = cylinder_record_probability(input,raw)
            @test exp(path.log_probability) ≈ state.probabilities[k] atol=1e-12
            @test getproperty(path.sector_weights,input.sector == :tau_flux ? :tau_flux : :vacuum_flux) ≈ 1 atol=1e-12
            @test state.joint_amplitudes[:,k] ≈ exp(path.log_probability/2)*path.final_state atol=1e-12
        end
        selected = cylinder_record_state(input,2;final_state=input.initial_state)
        @test sum(abs2,values(selected.amplitudes)) ≈ 1
        @test 0 < selected.postselection_probability <= 1
    end
    # Independent matrix action: weak record state = product of non-diagonal
    # physical filters acting on the projective coherent record state.
    strong = fibonacci_cylinder_reference(4;tau=Inf,initial_state=ref.initial_state)
    hard = cylinder_record_state(strong,2)
    soft = cylinder_record_state(ref,2)
    filters = [fibonacci_outcome_filter(r==2 ? ref.tau/2 : ref.tau) for c in 1:2 for r in 1:2]
    R = reduce(kron,reverse(filters))
    @test soft.joint_amplitudes ≈ hard.joint_amplitudes*transpose(R) atol=1e-12
    @test fibonacci_outcome_filter(Inf) ≈ I
    @test fibonacci_outcome_filter(0.0) ≈ ones(2,2)/sqrt(2)
    for t in (0.0,0.2,atanh(0.95),Inf)
        d = fibonacci_outcome_deformation(t)
        @test d.filter ≈ d.scale*d.rotation*Diagonal([exp(d.J),1])*d.rotation' atol=1e-12
    end
    @test_throws ArgumentError cylinder_record_state(ref,20;max_elements=100)
    prefix = cylinder_record_state(ref,1;total_layers=2)
    for k in 1:4
        raw = reshape([_sn_bit(k-1,e)==0 for e in 1:2],1,2)
        path = cylinder_record_probability(ref,raw;total_layers=2)
        @test exp(path.log_probability) ≈ prefix.probabilities[k]
    end
    raw = Bool[1 0;0 1]
    replay = cylinder_record_probability(ref,raw)
    data = Dict("L"=>4,"tau"=>0.7,"periods"=>1,"initial_state"=>"TCI_GS",
        "sample"=>raw,"sample_free_energy"=>Float32.(-replay.layer_log_probabilities),"seed"=>7)
    @test abs(only(audit_cft_trajectories(ref,[data])).log_probability_difference) < 1e-6
    # Actual small prefix inspected on hpc2ust, no full dataset copied.
    remote = fibonacci_cylinder_reference(8;tau=atanh(0.95))
    @test remote.initial_energy ≈ -6.196346083200089 atol=1e-12
    observed = cylinder_record_probability(remote,Bool[1 0 0 0;1 0 0 1];total_layers=64)
    @test observed.layer_log_probabilities ≈ [-4.3262553215026855,-3.1630616188049316] atol=3e-7 rtol=0
    @test observed.sector_weights.vacuum_flux ≈ 1 atol=1e-12
end
