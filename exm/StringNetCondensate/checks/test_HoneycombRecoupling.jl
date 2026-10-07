using Test, LinearAlgebra, Random
isdefined(@__MODULE__, :RecoupledHoneycombPEPS) || include(joinpath(@__DIR__, "..", "honeycomb_recoupling.jl"))

@testset "Honeycomb recoupling and cylinder boundaries" begin
    @testset "Local tetrahedral identity" begin
        for a in 0:1,b in 0:1,d in 0:1,s in 0:1
            @test fib_gsymbol(1,1,s,d,b,a) ≈ fib_fsymbol(a,1,1,d,b,s)/sqrt(_quantum_dimension(b)*_quantum_dimension(s)) atol=1e-14
            for bp in 0:1
                lhs = _honeycomb_vertex(s,a,b,d)*_honeycomb_vertex(s,a,bp,d)
                rhs = φ/sqrt(_quantum_dimension(s)) *
                      (_quantum_dimension(a)*_quantum_dimension(d))^(1/3) /
                      (_quantum_dimension(b)*_quantum_dimension(bp))^(1/3) *
                      fib_fsymbol(a,1,1,d,b,s)*fib_fsymbol(a,1,1,d,bp,s)
                @test lhs ≈ rhs atol=1e-14
            end
        end
    end

    @testset "Actual G-tensor contraction versus independent circuit" begin
        rng = MersenneTwister(310)
        for sector in (:vacuum_flux,:tau_flux), t in (Inf,0.8), nr in (2,3)
            ref = fibonacci_cylinder_reference(4;tau=t,sector)
            psi = recoupled_honeycomb_peps(4,nr;tau=t)
            actual = contract_recoupled_honeycomb(psi,ref)
            expected = cylinder_record_state(ref,nr)
            @test actual.joint_amplitudes ≈ expected.joint_amplitudes atol=5e-13
            @test actual.norm2 ≈ 1 atol=5e-13
            cache = HoneycombTransferCache(ref)
            for _ in 1:3
                record = bitrand(rng,nr,2)
                bits = sum(Int(!record[e])<<(e-1) for e in eachindex(record))
                direct = recoupled_honeycomb_record(psi,ref,record).final_vector
                @test direct ≈ expected.joint_amplitudes[:,bits+1] atol=5e-13
                transfer = honeycomb_record_probability(cache,record)
                @test direct ≈ exp(transfer.log_weight/2)*transfer.final_state atol=5e-13
                scalar = recoupled_honeycomb_record(psi,ref,record;final_state=ref.initial_state).amplitude
                @test scalar ≈ dot(ref.initial_state,direct) atol=5e-13
            end
            @test length(cache.operators)<=4*ref.L
        end
        # Complex, sector-mixed input tests phase and boundary conjugation.
        template = fibonacci_cylinder_reference(4;tau=0.4)
        v = randn(rng,ComplexF64,length(template.initial_state))
        ref = fibonacci_cylinder_reference(4;tau=0.4,initial_state=v)
        final = normalize(randn(rng,ComplexF64,length(v)))
        psi = recoupled_honeycomb_peps(4,2;tau=0.4,total_layers=64)
        got = contract_recoupled_honeycomb(psi,ref;final_state=final)
        expected = cylinder_record_state(ref,2;total_layers=64)
        @test got.joint_amplitudes ≈ expected.joint_amplitudes atol=5e-13
        @test got.postselected_amplitudes ≈ vec(final'*expected.joint_amplitudes) atol=5e-13
        # Tau=0, where all outcomes have equal probability.
        ref0 = fibonacci_cylinder_reference(4;tau=0)
        uniform = contract_recoupled_honeycomb(recoupled_honeycomb_peps(4,2;tau=0),ref0)
        @test uniform.probabilities ≈ fill(1/16,16) atol=1e-13
    end

    @testset "Boundary metric and virtual Y" begin
        ref = fibonacci_cylinder_reference(6;tau=0.7)
        q1,q2 = honeycomb_boundary_factors(ref,1),honeycomb_boundary_factors(ref,2)
        @test q1.*q2 ≈ ones(length(q1))
        for r in 1:3,side in (:input,:output)
            b = honeycomb_boundary_symmetry(ref,r;side)
            @test b.loop^2 ≈ I+b.loop atol=1e-12
            @test b.loop'*b.metric ≈ b.metric*b.loop atol=1e-12
            @test b.vacuum_flux*b.tau_flux ≈ zero(ref.Y) atol=1e-12
            q = honeycomb_boundary_factors(ref,r)
            v = ref.initial_state .* (side == :input ? inv.(q) : q)
            @test b.loop*v ≈ φ*v atol=1e-12
        end
        @test norm(honeycomb_vacuum_cap(ref,1)) ≈ 1
        # The loop can be moved through the actual paired-G transfer, with
        # the different input and output boundary coordinate conventions.
        cache = HoneycombTransferCache(ref)
        G = Matrix{Float64}(I,length(ref.initial_state),length(ref.initial_state))
        for r in 1:3,c in 1:ref.L÷2
            i = isodd(r) ? 2c : 2c-1
            G = _honeycomb_transfer(cache,i,ref.tau,isodd(r+c))*G
        end
        yin = honeycomb_boundary_symmetry(ref,1;side=:input).loop
        yout = honeycomb_boundary_symmetry(ref,3;side=:output).loop
        @test yout*G ≈ G*yin atol=1e-12
    end

    @testset "Existing vacuum-capped honeycomb PEPS: every amplitude" begin
        for (L,nr) in ((4,3),(4,4),(6,3))
            template = fibonacci_cylinder_reference(L;tau=Inf)
            ref = fibonacci_cylinder_reference(L;tau=Inf,initial_state=honeycomb_vacuum_cap(template,1))
            final = honeycomb_vacuum_cap(ref,nr)
            geometry = honeycomb_record_geometry(L,nr)
            old = fibonacci_stringnet_peps(nr-2,L÷2;bc=:cylinder_y,normalize=false)
            @test geometry.lattice.edges == old.lattice.edges
            @test sort(vcat(vec(geometry.record_edges[2:end-1,:]),geometry.zigzag_edges)) == collect(eachindex(old.links))
            @test all(v -> old.lattice.vertex_edges[v][3]==0,geometry.input_vertices)
            raw = HoneycombTransferCache(ref;representation=:stringnet)
            circuit = HoneycombTransferCache(ref)
            n = (nr-2)*(L÷2)
            sn,circ = ComplexF64[],ComplexF64[]
            for bits in 0:(1<<n)-1
                cfg = ones(Int,length(old.links))
                record = trues(nr,L÷2) # physical boundary rungs fixed vacuum
                for c in 1:L÷2,r in 2:nr-1
                    s = (bits>>((c-1)*(nr-2)+r-2))&1
                    cfg[geometry.record_edges[r,c]] = s
                    record[r,c] = s==0
                end
                amplitude = peps_amplitude(old,cfg)
                g = honeycomb_record_probability(raw,record;final_state=final)
                a = exp(g.log_weight/2)*g.boundary_overlap
                @test a/φ^(L/2) ≈ amplitude atol=5e-12
                k = honeycomb_record_probability(circuit,record;final_state=final)
                b = exp(k.log_weight/2)*k.boundary_overlap
                fugacity = prod(sqrt(_quantum_dimension(!x)) for x in record)/φ^length(record)
                @test b ≈ fugacity*a atol=5e-12
                push!(sn,amplitude); push!(circ,b)
            end
            normalize!(sn); normalize!(circ)
            # The rung dimension factor is physically necessary, not a scalar.
            @test abs2(dot(sn,circ)) < 1-1e-4
        end
    end

    @testset "Selected wider record, batch interface, and guards" begin
        ref = fibonacci_cylinder_reference(8;tau=atanh(0.95))
        record = Bool[1 0 0 0;1 0 0 1]
        psi = recoupled_honeycomb_peps(8,2;tau=ref.tau,total_layers=64)
        cache = HoneycombTransferCache(ref)
        g = honeycomb_record_probability(cache,record;total_layers=64)
        f = cylinder_record_probability(ref,record;total_layers=64)
        @test g.log_weight ≈ f.log_probability atol=1e-12
        @test g.final_state ≈ f.final_state atol=1e-12
        a = recoupled_honeycomb_record(psi,ref,record).final_vector
        @test a ≈ exp(g.log_weight/2)*g.final_state atol=1e-12
        saved = Dict("L"=>8,"tau"=>ref.tau,"periods"=>1,"sample"=>record,
                     "sample_free_energy"=>Float32.(-f.layer_log_probabilities),"seed"=>4695)
        # This saved record is now a complete two-layer run without half strength.
        audit = compare_honeycomb_trajectories(cache,[saved];terminal_half_strength=false)
        @test abs(only(audit).difference) < 3e-7
        @test_throws ArgumentError recoupled_honeycomb_peps(5,2)
        @test_throws ArgumentError honeycomb_record_geometry(4,2)
        @test_throws DimensionMismatch recoupled_honeycomb_record(psi,ref,trues(2,2))
        @test_throws ArgumentError contract_recoupled_honeycomb(psi,ref;max_elements=10)
        @test_throws ArgumentError compare_honeycomb_trajectories(HoneycombTransferCache(ref;representation=:stringnet),[saved])
    end
end
