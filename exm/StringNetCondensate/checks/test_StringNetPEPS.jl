using Test, ITensors, LinearAlgebra, Random
isdefined(@__MODULE__, :StringNetPEPS) || include(joinpath(@__DIR__, "..", "peps.jl"))

# Independent exact application of plaquette loop-insertion operators to |0>.
# Levin-Wen matrix element: product_v F(m_v, j_v, s, k_prev, j_prev, k_v),
# where j/k are the old/new boundary labels and m is the unchanged outer leg.
function apply_tau_loop(lattice, state, p)
    es = lattice.plaquettes[p]
    stems = map(1:6) do n
        prev, cur = es[mod1(n-1,6)], es[n]
        v = only(intersect(lattice.edges[prev], lattice.edges[cur]))
        only(filter(e -> e != prev && e != cur, lattice.vertex_edges[v]))
    end
    out = Dict{Int,Float64}()
    for (bits, amplitude) in state, newlabels in 0:63
        factor = amplitude
        newbits = bits
        for n in 1:6
            prev = mod1(n-1,6)
            j = (bits >> (es[n]-1)) & 1
            jp = (bits >> (es[prev]-1)) & 1
            k = (newlabels >> (n-1)) & 1
            kp = (newlabels >> (prev-1)) & 1
            m = stems[n] == 0 ? 0 : (bits >> (stems[n]-1)) & 1
            factor *= fib_fsymbol(m,j,1,kp,jp,k)
            iszero(factor) && break
            mask = 1 << (es[n]-1)
            newbits = (newbits & ~mask) | (k << (es[n]-1))
        end
        iszero(factor) || (out[newbits] = get(out,newbits,0.0) + factor)
    end
    return out
end

function projector_reference(lattice; order=eachindex(lattice.plaquettes))
    state = Dict(0 => 1.0)
    for p in order
        loop = apply_tau_loop(lattice,state,p)
        for (bits, amplitude) in loop
            state[bits] = get(state,bits,0.0) + φ * amplitude
        end
    end
    filter!(pair -> abs(last(pair)) > 1e-13, state)
    return state
end

# Close missing vacuum edges by exterior vacuum tadpoles. This supplies a
# trivalent planar graph for the independently implemented diagram evaluator.
function planar_reference(lattice, cfg)
    lattice.bc == :open || error("spherical evaluation does not fix cylinder sector")
    edges = [(u,v,Int(cfg[e])) for (e,(u,v)) in enumerate(lattice.edges)]
    vertex_edges = [collect(es) for es in lattice.vertex_edges]
    for v in eachindex(lattice.vertices)
        if vertex_edges[v][3] == 0
            far = length(vertex_edges) + 1
            stem = length(edges) + 1
            push!(edges,(v,far,0),(far,far,0))
            vertex_edges[v][3] = stem
            push!(vertex_edges,[stem,stem+1])
        end
    end
    return eval_diagram(build_diagram(length(vertex_edges),edges,vertex_edges))
end

@testset "Fibonacci PEPS" begin
    @testset "G-symbol and geometry" begin
        for i in 0:1, j in 0:1, k in 0:1, a in 0:1, b in 0:1, c in 0:1
            @test fib_gsymbol(i,j,k,a,b,c) ≈ fib_gsymbol(j,k,i,b,c,a) atol=1e-14
        end
        for bc in (:open,:cylinder_x,:cylinder_y), Lx in 2:4, Ly in 2:4
            lat = HoneycombPatch(Lx,Ly;bc)
            @test length(lat.plaquettes) == Lx*Ly
            @test length(lat.vertices)-length(lat.edges)+Lx*Ly == (bc == :open ? 1 : 0)
            @test all(length(unique(p)) == 6 for p in lat.plaquettes)
            @test all(mod(lat.vertices[u][1],3) != mod(lat.vertices[v][1],3) for (u,v) in lat.edges)
        end
        @test_throws ArgumentError fibonacci_stringnet_peps(0,1)
        @test_throws ArgumentError fibonacci_stringnet_peps(1,1;bc=:torus)
        @test_throws ArgumentError fibonacci_stringnet_peps(1,2;bc=:cylinder_x)
        @test_throws ArgumentError fibonacci_stringnet_peps(2,1;bc=:cylinder_y)
    end

    @testset "Dense state versus independent plaquette projectors" begin
        for (Lx,Ly,bc) in ((1,1,:open),(2,1,:open),(1,2,:open),
                           (2,1,:cylinder_x),(1,2,:cylinder_y),
                           (2,2,:cylinder_x),(2,2,:cylinder_y))
            psi = fibonacci_stringnet_peps(Lx,Ly;bc)
            lat = psi.lattice
            ne = length(lat.edges)
            ref = projector_reference(lat)
            nrm = sqrt(sum(abs2, values(ref)))
            @test nrm ≈ (1 + φ^2)^(Lx*Ly/2) atol=1e-11
            expected = zeros(2^ne)
            for (bits,a) in ref
                expected[bits+1] = a/nrm
            end
            dense = ITensors.@set_warn_order 40 contract_peps(psi)
            actual = vec(Array(dense,psi.physical_sites...))
            @test actual ≈ expected atol=1e-11
            @test norm(actual) ≈ 1.0 atol=1e-12
            for p in eachindex(lat.plaquettes)
                loop = apply_tau_loop(lat,ref,p)
                # B^tau |Psi> = phi |Psi>; hence B_p |Psi> = |Psi>.
                @test all(isapprox(get(loop,k,0.0),φ*get(ref,k,0.0);atol=1e-10)
                          for k in union(keys(loop),keys(ref)))
            end
            @test all(count(T -> hasind(T,s),psi.tensors) == 1 for s in psi.physical_sites)
            @test all(count(T -> hasind(T,s),psi.tensors) == 2 for s in psi.links)
            @test sort(vcat(psi.physical_edges...)) == collect(1:ne)
            @test all(dim(psi.links[e]) == (0 in lat.edge_faces[e] ? 2 : 5) for e in 1:ne)
            @test peps_amplitude(psi,zeros(Int,ne)) ≈ 1/nrm atol=1e-12
            single = zeros(Int,ne); single[1] = 1
            @test abs(peps_amplitude(psi,single)) < 1e-12
            if bc == :open
                for (bits,a) in ref
                    cfg = [(bits >> (e-1)) & 1 for e in 1:ne]
                    @test planar_reference(lat,cfg) ≈ a atol=1e-11
                end
            end
        end
    end

    @testset "Two-dimensional amplitudes and input guards" begin
        rng = MersenneTwister(103)
        for bc in (:open,:cylinder_x,:cylinder_y)
            psi = fibonacci_stringnet_peps(3,2;bc,normalize=false)
            lat = psi.lattice
            ref = projector_reference(lat)
            reverse_ref = projector_reference(lat;order=reverse(eachindex(lat.plaquettes)))
            @test all(isapprox(get(ref,k,0.0),get(reverse_ref,k,0.0);atol=1e-10)
                      for k in union(keys(ref),keys(reverse_ref)))
            bits_to_test = shuffle(rng, collect(keys(ref)))[1:12]
            for bits in bits_to_test
                cfg = [(bits >> (e-1)) & 1 for e in eachindex(lat.edges)]
                @test peps_amplitude(psi,cfg) ≈ ref[bits] atol=1e-10
                if bc == :open
                    @test planar_reference(lat,cfg) ≈ ref[bits] atol=1e-10
                end
            end
            @test peps_amplitude(psi,zeros(Int,length(lat.edges))) ≈ 1.0 atol=1e-12
            @test_throws ArgumentError contract_peps(psi;max_elements=64)
            @test_throws ArgumentError peps_amplitude(psi,zeros(Int,length(lat.edges));max_elements=1)
            @test_throws DimensionMismatch peps_amplitude(psi,[0])
            @test_throws ArgumentError peps_amplitude(psi,fill(2,length(lat.edges)))
        end
    end
end
