using Test, Random, LinearAlgebra
isdefined(@__MODULE__, :Diagram) || include(joinpath(@__DIR__, "..", "stringnet.jl"))

# Independent deletion-contraction of the dual graph. Vacuum primal edges
# contract dual faces; occupied primal edges prohibit equal face colors.
function chromatic(n, edges, q, memo)
    key = (n, Tuple(sort!(collect(edges))))
    return get!(memo, key) do
        any(a == b for (a, b) in edges) && return 0.0
        isempty(edges) && return q^n
        a, b = first(edges)
        deleted = setdiff(edges, Set([(a, b)]))
        relabel(v) = v == b ? a : v - (v > b)
        contracted = Set(minmax(relabel(u), relabel(v)) for (u, v) in deleted)
        chromatic(n, deleted, q, memo) - chromatic(n - 1, contracted, q, memo)
    end
end

function chromatic_weight(lad, cfg, memo)
    D = ladder_diagram(lad, cfg)
    fs = faces(D)
    face_id = zeros(Int, length(D.alive))
    for (i, face) in enumerate(fs), h in face
        face_id[h] = i
    end
    # Separate reference union implementation, independent of entropy helpers.
    colors = collect(1:length(fs))
    for e in eachindex(cfg)
        cfg[e] == 0 || continue
        a, b = colors[face_id[2e - 1]], colors[face_id[2e]]
        replace!(colors, b => a)
    end
    names = sort!(unique(colors))
    ids = Dict(v => i for (i, v) in enumerate(names))
    edges = Set(minmax(ids[colors[face_id[2e - 1]]], ids[colors[face_id[2e]]])
                for e in eachindex(cfg) if cfg[e] == 1)
    return chromatic(length(names), edges, φ + 2, memo) / (φ + 2)
end

function dense_entropy(cfgs, amps, L, ell)
    C = zeros(eltype(amps), 2^(3ell), 2^(3(L - ell)))
    for (cfg, amp) in zip(cfgs, amps)
        a, b = 0, 0
        for rail in 0:2, i in 1:L
            if i <= ell
                a = 2a + cfg[rail * L + i]
            else
                b = 2b + cfg[rail * L + i]
            end
        end
        C[a + 1, b + 1] += amp
    end
    return -sum(p > 0 ? p * log(p) : 0.0 for p in abs2.(svdvals(C)))
end

@testset "Fibonacci string-net condensate" begin
    @testset "Fusion and F-symbol coherence" begin
        @test !fusion_allowed(0, -1, 3)
        @test fib_fsymbol(1, 1, 1, 0, 1, 1) == 1.0
        for a in 0:1, b in 0:1, c in 0:1, d in 0:1
            left = [e for e in 0:1 if fusion_allowed(a,b,e) && fusion_allowed(e,c,d)]
            right = [f for f in 0:1 if fusion_allowed(b,c,f) && fusion_allowed(a,f,d)]
            F = [fib_fsymbol(a,b,c,d,e,f) for e in left, f in right]
            @test F * F' ≈ Matrix{Float64}(I, length(left), length(left)) atol=1e-14
        end
        for a in 0:1, b in 0:1, c in 0:1, d in 0:1,
            e in 0:1, f in 0:1, g in 0:1, h in 0:1, t in 0:1
            lhs = sum(fib_fsymbol(a,b,c,f,e,n) * fib_fsymbol(a,n,d,t,f,g) *
                      fib_fsymbol(b,c,d,g,n,h) for n in 0:1)
            rhs = fib_fsymbol(e,c,d,t,f,h) * fib_fsymbol(a,b,h,t,e,g)
            @test lhs ≈ rhs atol=1e-14
        end
    end

    @testset "Enumeration and sphere evaluation" begin
        @test_throws ArgumentError Ladder(1)
        @test !config_valid(Ladder(2), [0])
        for L in 2:5
            lad = Ladder(L)
            cfgs = ladder_configs(lad)
            # Brute-force all 2^(3L) bit strings, independent of transfer paths.
            brute = Set(Tuple((s >> (e - 1)) & 1 for e in 1:3L)
                        for s in 0:(2^(3L)-1)
                        if config_valid(lad, [(s >> (e - 1)) & 1 for e in 1:3L]))
            @test Set(Tuple.(cfgs)) == brute
            @test length(cfgs) == length(brute)
            memo_a, memo_b = Dict{String,Float64}(), Dict{String,Float64}()
            chromemo = Dict{Any,Float64}()
            raw = Float64[]
            for cfg in cfgs
                D = ladder_diagram(lad, cfg)
                a = eval_diagram(D, memo_a)
                b = eval_diagram(D, memo_b; chooser=choose_fmove_edge_alt)
                @test a ≈ b atol=1e-11
                @test a^2 ≈ chromatic_weight(lad, cfg, chromemo) atol=1e-10
                push!(raw, a)
                if L == 3
                    # Check every local move, including vacuum labels that the
                    # recursive evaluator normally absorbs before recoupling.
                    for edge in 1:3L
                        he = 2edge - 1
                        _, i, j = vertex_halfedges(D, he)
                        _, m, k = vertex_halfedges(D, D.he_twin[he])
                        recoupled = sum(fib_fsymbol(D.he_label[i], D.he_label[j],
                            D.he_label[m], D.he_label[k], D.he_label[he], n) *
                            eval_diagram(fmove_diagram(D, he, i, j, m, k, n), memo_a)
                            for n in 0:1 if fusion_allowed(D.he_label[j], D.he_label[m], n) &&
                                           fusion_allowed(D.he_label[k], D.he_label[i], n))
                        @test a ≈ recoupled atol=1e-11
                    end
                end
            end
            @test sum(abs2, raw) ≈ (1 + φ^2)^(L + 1) rtol=1e-12
            _, amps = stringnet_ground_state(lad)
            @test amps ≈ normalize(raw) atol=1e-12
            @test eval_diagram(ladder_diagram(lad, zeros(Int, 3L))) ≈ 1.0
            winding = zeros(Int, 3L)
            winding[1:L] .= 1
            @test eval_diagram(ladder_diagram(lad, winding)) ≈ φ
        end
        # Vacuum self-loop connected to a theta graph by a vacuum stem.
        # Suppression must retain the tau loop on the remote side.
        D = build_diagram(4, [(1,1,0), (1,2,0), (2,3,1), (2,4,1),
                              (3,4,1), (3,4,0)], [[1,2], [2,3,4], [3,5,6], [4,6,5]])
        @test eval_diagram(D) ≈ φ
        @test D.alive == trues(12) # public evaluation must not mutate its input
        absorbed = copy(D)
        @test absorb_vacuum!(absorbed, 1)
        @test all(all(absorbed.alive[k] && absorbed.he_vert[k] == absorbed.he_vert[h]
                      for k in vertex_halfedges(absorbed, h))
                  for h in eachindex(absorbed.alive)
                  if absorbed.alive[h] && absorbed.he_vert[h] != 0)
        @test eval_diagram(absorbed) ≈ φ
        # Two disconnected plain loops also exercise subdiagram's null next.
        loops = build_diagram(0, [(0,0,1), (0,0,1)], Vector{Int}[])
        @test eval_diagram(loops) ≈ φ^2
        @test eval_diagram(subdiagram(loops, [1,2])) ≈ φ
        # Canonicalization must ignore half-edge/vertex numbering, including
        # sparse vertex IDs, while retaining the rotation system and labels.
        D = ladder_diagram(Ladder(6), ones(Int, 18))
        rng = MersenneTwister(71)
        permutation = randperm(rng, length(D.alive))
        renamed = copy(D)
        for h in eachindex(permutation)
            p = permutation[h]
            renamed.he_twin[p] = permutation[D.he_twin[h]]
            renamed.he_next[p] = permutation[D.he_next[h]]
            renamed.he_vert[p] = 100 + 3D.he_vert[h]
            renamed.he_label[p] = D.he_label[h]
        end
        renamed.nv = maximum(renamed.he_vert)
        @test canonical_key(D) == canonical_key(renamed)
        @test eval_diagram(D) ≈ eval_diagram(renamed) atol=1e-12
    end

    @testset "Schmidt entropy" begin
        rng = MersenneTwister(831)
        for L in 2:4
            lad = Ladder(L)
            cfgs, amps = stringnet_ground_state(lad)
            # Complex random states test conjugation and avoid assuming any
            # translation symmetry or fixed-point structure in arc_entropy.
            random_amps = normalize(randn(rng, ComplexF64, length(cfgs)))
            for ell in 0:L, a in (amps, random_amps)
                @test arc_entropy(cfgs, a, lad, ell) ≈
                      dense_entropy(cfgs, a, L, ell) atol=1e-11
            end
        end
        lad = Ladder(6)
        cfgs, amps = stringnet_ground_state(lad)
        ss = [arc_entropy(cfgs, amps, lad, ell) for ell in 1:5]
        @test ss ≈ reverse(ss) atol=1e-11
        @test maximum(ss[2:end-1]) - minimum(ss[2:end-1]) < 1e-11
        @test ss[1] < ss[2] # overlapping boundaries of a single column
        @test_throws ArgumentError arc_entropy(cfgs, amps, lad, -1)
        @test_throws ArgumentError arc_entropy(cfgs, amps, lad, 7)
        @test_throws DimensionMismatch arc_entropy(cfgs, amps[2:end], lad, 1)
        @test_throws ArgumentError arc_entropy(cfgs, 2amps, lad, 1)
        # Bits beyond UInt64 must not silently alias. One qubit on either
        # side changes, producing a Bell pair even for a very unequal cut.
        c0, c1 = zeros(Int, 69), zeros(Int, 69)
        c1[1] = c1[69] = 1
        @test arc_entropy([c0,c1], fill(inv(sqrt(2)),2), Ladder(23),1) ≈ log(2)
        # Repeated configurations add coherently before normalization.
        @test arc_entropy([c0,c0,c1], [0.5,0.5,1] ./ sqrt(2), Ladder(23),1) ≈ log(2)
    end
end

# Larger-system coverage. Smaller sizes are already covered above by
# dense-entropy comparisons and the L=6 plateau checks.
@testset "L=8 ladder entropy" begin
    lad = Ladder(8)
    cfgs, amps = stringnet_ground_state(lad)
    ss = [arc_entropy(cfgs, amps, lad, ell) for ell in 1:7]
    @test sum(abs2, amps) ≈ 1.0 atol=1e-12
    @test ss ≈ reverse(ss) atol=1e-10
    @test maximum(ss[2:end-1]) - minimum(ss[2:end-1]) < 1e-10
end
