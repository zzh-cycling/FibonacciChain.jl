# Fibonacci string-net condensate (Fib-SNC) on a quasi-1D strip of plaquettes.
#
# The condensate (Levin-Wen fixed-point wavefunction, Levin & Wen 2005;
# Minev et al. 2025, arXiv:2406.12820) is constructed on a ring of L square
# plaquettes: two horizontal rails (top t_i, bottom b_i) joined by L vertical
# rungs v_i, i = 1..L (indices mod L). Every vertex is trivalent:
#   T_i = (v_i, t_i, t_{i-1}),  B_i = (v_i, b_{i-1}, b_i).
# Edges carry labels 0 (vacuum) or 1 (tau string). Allowed vertex triples obey
# the Fibonacci fusion rules {0,0,0}, {0,τ,τ}, {τ,τ,τ}.
#
# The wavefunction is built by direct diagram evaluation:
#   Ψ(G) = Eval(ladder diagram with edge labels G),
# using the Fibonacci diagrammatic rules:
#   vacuum absorption, plain τ loop = d_τ = φ, τ tadpole = 0,
#   τ digon = √φ (from theta(τ,τ,τ) = φ^{3/2}), and F-moves with
#   F^{τττ}_τ = [[φ^{-1}, φ^{-1/2}], [φ^{-1/2}, -φ^{-1}]] (rows/cols {0,1});
#   all F-symbols with a vacuum upper index equal 1.
# The annulus is capped by vacuum disks when evaluating the diagram on the
# sphere: winding loops also evaluate to φ. Other annular sectors/boundary
# conditions require additional data and are not represented here.

using LinearAlgebra
using SparseArrays

const φ = (1 + sqrt(5)) / 2
const SQRTφ = sqrt(φ)

# --------------------------------------------------------------------------
# Fibonacci F-symbol: F^{ijm}_{k;ln}, labels in {0 (vacuum), 1 (τ)}.
# Recoupling ((i×j→l)×m→k) = Σ_n F^{ijm}_{k;ln} (i×(j×m→n)→k).
# --------------------------------------------------------------------------

fusion_allowed(a::Integer, b::Integer, c::Integer) =
    (0 <= a <= 1 && 0 <= b <= 1 && 0 <= c <= 1) && a + b + c != 1

function fib_fsymbol(i::Integer, j::Integer, m::Integer, k::Integer, l::Integer, n::Integer)
    fusion_allowed(i, j, l) || return 0.0
    fusion_allowed(l, m, k) || return 0.0
    fusion_allowed(j, m, n) || return 0.0
    fusion_allowed(i, n, k) || return 0.0
    (i == 0 || j == 0 || m == 0) && return 1.0
    # i = j = m = τ
    if k == 0
        return 1.0                  # l = n = τ; the block is one-dimensional
    else
        l == 0 && n == 0 && return φ^(-1)
        l == 0 && n == 1 && return φ^(-1 / 2)
        l == 1 && n == 0 && return φ^(-1 / 2)
        return -φ^(-1)
    end
end

# --------------------------------------------------------------------------
# Diagram: combinatorial map (rotation system) for planar trivalent
# string-net diagrams, with vacuum (0) edge labels and vertex-free plain
# loop components. Half-edges are stored explicitly with a twin map.
#   he_twin[h]   : the other half-edge of the same edge
#   he_vert[h]   : vertex id the half-edge is attached to (0 = none)
#   he_next[h]   : next half-edge in ccw order around he_vert[h]
#   he_label[h]  : 0 (vacuum) or 1 (τ); identical for twins
#   alive[h]     : false once the half-edge has been removed
# An edge with both half-edges vertex-less is a plain loop component.
# --------------------------------------------------------------------------

mutable struct Diagram
    he_twin::Vector{Int}
    he_vert::Vector{Int}
    he_next::Vector{Int}
    he_label::Vector{Int8}
    alive::Vector{Bool}
    nv::Int                      # number of vertex ids issued (may be stale)
end

Diagram(nhe::Int) = Diagram(zeros(Int, nhe), zeros(Int, nhe), zeros(Int, nhe),
                            zeros(Int8, nhe), falses(nhe), 0)

# Build a diagram from an edge list (va, vb, label), vertex ids 1..nv,
# 0 = no vertex on that side, plus the ccw ordered edge ids at each vertex.
function build_diagram(nv::Int, edges::Vector{Tuple{Int,Int,Int}},
                       vertex_edges::Vector{Vector{Int}})
    ne = length(edges)
    D = Diagram(2ne)
    for (e, (va, vb, lab)) in enumerate(edges)
        h1, h2 = 2e - 1, 2e
        D.he_twin[h1] = h2
        D.he_twin[h2] = h1
        D.he_vert[h1] = va
        D.he_vert[h2] = vb
        D.he_label[h1] = lab
        D.he_label[h2] = lab
        D.alive[h1] = D.alive[h2] = true
    end
    D.nv = nv
    for (v, elist) in enumerate(vertex_edges)
        isempty(elist) && continue
        hs = Int[]
        for e in elist
            va, vb, _ = edges[e]
            va == v && push!(hs, 2e - 1)
            vb == v && push!(hs, 2e)
        end
        length(hs) == 3 || error("vertex $v has $(length(hs)) half-edges")
        for k in 1:3
            D.he_next[hs[k]] = hs[mod1(k + 1, 3)]
        end
    end
    return D
end

Base.copy(D::Diagram) = Diagram(copy(D.he_twin), copy(D.he_vert),
                                copy(D.he_next), copy(D.he_label),
                                copy(D.alive), D.nv)

# half-edges of the vertex containing h, in ccw order starting from h
function vertex_halfedges(D::Diagram, h::Int)
    h2 = D.he_next[h]
    return (h, h2, D.he_next[h2])
end

# --------------------------------------------------------------------------
# Local simplification rules.
# --------------------------------------------------------------------------

function delete_edge!(D::Diagram, h::Int)
    D.alive[h] = false
    D.alive[D.he_twin[h]] = false
    return D
end

# Splice: join half-edges o1, o2 into a single edge.
function splice!(D::Diagram, o1::Int, o2::Int)
    D.he_twin[o1] = o2
    D.he_twin[o2] = o1
    return D
end

# Absorb a vacuum edge (half-edge h). Returns false if a vertex is
# inconsistent (diagram evaluates to zero).
function absorb_vacuum!(D::Diagram, h::Int)
    ht = D.he_twin[h]
    va, vb = D.he_vert[h], D.he_vert[ht]
    if va == 0 && vb == 0
        delete_edge!(D, h)               # plain vacuum loop: factor 1
        return true
    end
    if va == vb
        # Absorb the vacuum stem first. This suppresses BOTH endpoints;
        # detaching its far half-edge would corrupt the far vertex.
        g = only(filter(x -> x != h && x != ht, vertex_halfedges(D, h)))
        D.he_label[g] == 0 || return false
        return absorb_vacuum!(D, g)
    end
    for hv in (h, ht)
        x, y = filter(!=(hv), vertex_halfedges(D, hv))
        D.he_label[x] == D.he_label[y] || return false
        xt, yt = D.he_twin[x], D.he_twin[y]
        if xt == y
            # x and y are halves of one edge -> plain loop: keep it alive,
            # just detach it from the vertex
            D.he_vert[x] = 0
            D.he_vert[y] = 0
        else
            splice!(D, xt, yt)
            D.alive[x] = false
            D.alive[y] = false
        end
    end
    delete_edge!(D, h)
    return true
end

# Connected components of live half-edges (adjacency via twin and vertex).
function components(D::Diagram)
    seen = falses(length(D.alive))
    comps = Vector{Int}[]
    for h in eachindex(D.alive)
        (D.alive[h] && !seen[h]) || continue
        stack = [h]
        seen[h] = true
        comp = Int[]
        while !isempty(stack)
            x = pop!(stack)
            push!(comp, x)
            for y in (D.he_twin[x], D.he_next[x])
                if y != 0 && D.alive[y] && !seen[y]
                    seen[y] = true
                    push!(stack, y)
                end
            end
        end
        push!(comps, comp)
    end
    return comps
end

# Copy of the diagram restricted to one component.
function subdiagram(D::Diagram, comp::Vector{Int})
    sub = Diagram(length(comp))
    ren = Dict{Int,Int}(h => i for (i, h) in enumerate(comp))
    for (i, h) in enumerate(comp)
        sub.he_twin[i] = ren[D.he_twin[h]]
        sub.he_vert[i] = D.he_vert[h]
        sub.he_next[i] = D.he_vert[h] == 0 ? 0 : ren[D.he_next[h]]
        sub.he_label[i] = D.he_label[h]
        sub.alive[i] = true
    end
    sub.nv = D.nv
    return sub
end

# Remove all plain loops; returns accumulated factor (φ per τ loop).
function pop_plain_loops!(D::Diagram)
    factor = 1.0
    changed = false
    for h in eachindex(D.alive)
        D.alive[h] || continue
        if D.he_vert[h] == 0 && D.he_vert[D.he_twin[h]] == 0
            D.he_label[h] == 1 && (factor *= φ)
            delete_edge!(D, h)
            changed = true
        end
    end
    return factor, changed
end

# One pass of cheap reductions; returns (factor, ok, changed).
function simplify_pass!(D::Diagram)
    factor = 1.0
    changed = false
    # vacuum absorption
    for h in eachindex(D.alive)
        D.alive[h] || continue
        if D.he_label[h] == 0
            absorb_vacuum!(D, h) || return (1.0, false, true)
            changed = true
        end
    end
    f, c = pop_plain_loops!(D)
    factor *= f
    changed |= c
    # tadpoles: τ self-loop at a vertex -> zero (stem is τ post-absorption)
    for h in eachindex(D.alive)
        D.alive[h] || continue
        ht = D.he_twin[h]
        v, vt = D.he_vert[h], D.he_vert[ht]
        if v != 0 && v == vt && D.he_label[h] == 1
            return (1.0, false, true)
        end
    end
    # digons: two parallel τ edges between distinct vertices
    for h in eachindex(D.alive)
        D.alive[h] || continue
        v1, v2 = D.he_vert[h], D.he_vert[D.he_twin[h]]
        (v1 != 0 && v2 != 0 && v1 != v2) || continue
        hs1 = vertex_halfedges(D, h)
        for c in hs1
            c == h && continue
            ct = D.he_twin[c]
            D.he_vert[ct] == v2 || continue
            # digon (h, c); legs: e (third edge at v1), d (third edge at v2)
            e = only(filter(x -> x != h && x != c, hs1))
            hs2 = vertex_halfedges(D, D.he_twin[h])
            d = only(filter(x -> x != D.he_twin[h] && x != ct, hs2))
            # digon value √φ; join the legs e--d
            et, dt = D.he_twin[e], D.he_twin[d]
            factor *= SQRTφ
            if et == d
                # e and d are halves of one edge -> plain loop, keep alive
                D.he_vert[e] = 0
                D.he_vert[d] = 0
            else
                splice!(D, et, dt)
                D.alive[e] = false
                D.alive[d] = false
            end
            for x in (h, D.he_twin[h], c, ct)
                D.alive[x] = false
            end
            changed = true
            break
        end
        changed && break
    end
    return (factor, true, changed)
end

# Simplify until stable. Returns (factor, ok); ok=false means zero diagram.
function simplify!(D::Diagram)
    factor = 1.0
    while true
        f, ok, changed = simplify_pass!(D)
        factor *= f
        ok || return (factor, false)
        changed || break
    end
    return (factor, true)
end

# --------------------------------------------------------------------------
# Canonical form for memoization: minimal BFS serialization over all
# starting half-edges. Assumes a connected diagram with at least one vertex.
# --------------------------------------------------------------------------

# Keys are opaque binary strings (three UInt32 fields per live half-edge).
# Reuse dense scratch maps across roots instead of allocating dictionaries
# and formatting decimal integers at every BFS step.
function canonical_key(D::Diagram)
    nhe = length(D.alive)
    vmap, emap, queue = zeros(Int, nhe), zeros(Int, nhe), zeros(Int, nhe)
    buffer = Vector{UInt32}(undef, 3nhe)
    best = nothing
    for h0 in eachindex(D.alive)
        D.alive[h0] || continue
        key = _bfs_key!(D, h0, vmap, emap, queue, buffer)
        if best === nothing || isless(key, best)
            best = key
        end
    end
    return best === nothing ? "" : best
end

function bfs_key(D::Diagram, h0::Int)
    nhe = length(D.alive)
    return _bfs_key!(D, h0, zeros(Int, nhe), zeros(Int, nhe),
                     zeros(Int, nhe), Vector{UInt32}(undef, 3nhe))
end

function _bfs_key!(D, h0, vmap, emap, queue, buffer)
    fill!(vmap, 0)
    fill!(emap, 0)
    for h in vertex_halfedges(D, h0)
        vmap[h] = 1
    end
    queue[1] = h0
    head, tail, nextv, nexte, pos = 1, 1, 1, 0, 0
    while head <= tail
        hstart = queue[head]
        head += 1
        for h in vertex_halfedges(D, hstart)
            t = D.he_twin[h]
            if emap[h] == 0
                nexte += 1
                emap[h] = emap[t] = nexte
            end
            if vmap[t] == 0
                nextv += 1
                for ht in vertex_halfedges(D, t)
                    vmap[ht] = nextv
                end
                tail += 1
                queue[tail] = t
            end
            buffer[pos + 1] = emap[h]
            buffer[pos + 2] = D.he_label[h]
            buffer[pos + 3] = vmap[t]
            pos += 3
        end
    end
    return String(copy(reinterpret(UInt8, @view buffer[1:pos])))
end

# --------------------------------------------------------------------------
# Diagram evaluation by memoized recursion.
# --------------------------------------------------------------------------

function eval_diagram(D::Diagram, memo::Dict{String,Float64}=Dict{String,Float64}();
                      chooser::Function=choose_fmove_edge)
    return _eval_diagram!(copy(D), memo, chooser)
end

# Internal ownership-taking evaluator: recursive branches already own a copy.
function _eval_diagram!(D::Diagram, memo::Dict{String,Float64}, chooser::F) where {F}
    factor, ok = simplify!(D)
    ok || return 0.0
    comps = components(D)
    isempty(comps) && return factor
    if length(comps) > 1
        val = factor
        for c in comps
            val *= _eval_diagram!(subdiagram(D, c), memo, chooser)
        end
        return val
    end
    key = canonical_key(D)
    if haskey(memo, key)
        return factor * memo[key]
    end
    val = _eval_connected(D, memo, chooser)
    memo[key] = val
    return factor * val
end

# Connected, fully simplified diagram: apply an F-move and recurse.
function _eval_connected(D::Diagram, memo::Dict{String,Float64}, chooser::F) where {F}
    h = chooser(D)
    h == 0 && error("no F-move available for diagram")
    ht = D.he_twin[h]
    hs1 = vertex_halfedges(D, h)           # (h, i, j) ccw at v1
    i, j = hs1[2], hs1[3]
    hs2 = vertex_halfedges(D, ht)          # (ht, m, k) ccw at v2; pairs j-m, i-k
    m, k = hs2[2], hs2[3]
    il, jl = D.he_label[i], D.he_label[j]
    ml, kl = D.he_label[m], D.he_label[k]
    ll = D.he_label[h]
    total = 0.0
    for n in (0, 1)
        f = fib_fsymbol(il, jl, ml, kl, ll, n)
        f == 0.0 && continue
        Dn = fmove_diagram(D, h, i, j, m, k, n)
        total += f * _eval_diagram!(Dn, memo, chooser)
    end
    return total
end

# Recoupled diagram after an F-move on edge h:
# v1=(h,i,j), v2=(ht,m,k) are replaced by v1'=(h,j,m), v2'=(ht,k,i).
function fmove_diagram(D::Diagram, h::Int, i::Int, j::Int, m::Int, k::Int, nlab::Int)
    Dn = copy(D)
    ht = Dn.he_twin[h]
    w1, w2 = Dn.he_vert[h], Dn.he_vert[ht]
    Dn.he_label[h] = nlab
    Dn.he_label[ht] = nlab
    Dn.he_vert[j] = w1
    Dn.he_vert[m] = w1
    Dn.he_vert[k] = w2
    Dn.he_vert[i] = w2
    Dn.he_vert[h] = w1
    Dn.he_vert[ht] = w2
    Dn.he_next[h] = j
    Dn.he_next[j] = m
    Dn.he_next[m] = h
    Dn.he_next[ht] = k
    Dn.he_next[k] = i
    Dn.he_next[i] = ht
    return Dn
end

# Face tracing: returns the list of half-edges of every face (each face is
# the cycle on the left of a directed half-edge).
function faces(D::Diagram)
    seen = falses(length(D.alive))
    fs = Vector{Int}[]
    for h0 in eachindex(D.alive)
        (D.alive[h0] && !seen[h0]) || continue
        h = h0
        face = Int[]
        while true
            seen[h] = true
            push!(face, h)
            t = D.he_twin[h]
            D.he_vert[t] == 0 && break
            h = D.he_next[t]
            h == h0 && break
        end
        push!(fs, face)
    end
    return fs
end

# Choose an edge for the F-move: the first edge of a smallest face of
# size >= 3 (this strategy shrinks small faces and terminates).
function choose_fmove_edge(D::Diagram)
    best, best_sz = 0, typemax(Int)
    for face in faces(D)
        if length(face) >= 3 && length(face) < best_sz
            best_sz = length(face)
            best = face[1]
        end
    end
    return best
end

# Alternative chooser (for path-independence checks): same smallest-face
# strategy, but pick the LAST edge of that face. Different reduction path,
# still terminating.
function choose_fmove_edge_alt(D::Diagram)
    best, best_sz = 0, typemax(Int)
    for face in faces(D)
        if length(face) >= 3 && length(face) < best_sz
            best_sz = length(face)
            best = face[end]
        end
    end
    return best
end

# --------------------------------------------------------------------------
# Ladder geometry: ring of L plaquettes.
# Edges: t_i = i, b_i = L + i, v_i = 2L + i  (i = 1..L, mod L).
# Vertices: T_i = (v_i, t_i, t_{i-1}) ccw, B_i = (v_i, b_{i-1}, b_i) ccw.
# A configuration is a length-3L 0/1 vector.
# --------------------------------------------------------------------------

struct Ladder
    L::Int
    function Ladder(L::Int)
        L >= 2 || throw(ArgumentError("need L >= 2"))
        return new(L)
    end
end

t_edge(lad::Ladder, i::Int) = mod1(i, lad.L)
b_edge(lad::Ladder, i::Int) = lad.L + mod1(i, lad.L)
v_edge(lad::Ladder, i::Int) = 2lad.L + mod1(i, lad.L)

function config_valid(lad::Ladder, cfg)::Bool
    length(cfg) == 3lad.L || return false
    for i in 1:lad.L
        fusion_allowed(cfg[v_edge(lad, i)], cfg[t_edge(lad, i)],
                       cfg[t_edge(lad, i - 1)]) || return false
        fusion_allowed(cfg[v_edge(lad, i)], cfg[b_edge(lad, i - 1)],
                       cfg[b_edge(lad, i)]) || return false
    end
    return true
end

# Four transfer states (t_i, b_i). Keep one work buffer and copy only
# complete, closed configurations; no partial-configuration dictionaries.
function ladder_configs(lad::Ladder)
    L = lad.L
    transitions = [Tuple{Int,Int}[] for _ in 1:4]
    for prev in 0:3, next in 0:3, v in 0:1
        fusion_allowed(v, prev >> 1, next >> 1) || continue
        fusion_allowed(v, prev & 1, next & 1) || continue
        push!(transitions[prev + 1], (next, v))
    end
    out = Vector{Int}[]
    cfg = zeros(Int, 3L)
    for start in 0:3
        # can_close[r+1, s+1]: a path of r columns from s to start exists.
        can_close = falses(L + 1, 4)
        can_close[1, start + 1] = true
        for r in 1:L, prev in 0:3
            can_close[r + 1, prev + 1] =
                any(can_close[r, next + 1] for (next, _) in transitions[prev + 1])
        end
        _extend_configs!(out, cfg, transitions, can_close, L, 1, start)
    end
    return out
end

function _extend_configs!(out, cfg, transitions, can_close, L, i, prev)
    for (next, v) in transitions[prev + 1]
        can_close[L - i + 1, next + 1] || continue
        cfg[i] = next >> 1
        cfg[L + i] = next & 1
        cfg[2L + i] = v
        if i == L
            push!(out, copy(cfg))
        else
            _extend_configs!(out, cfg, transitions, can_close, L, i + 1, next)
        end
    end
    return out
end

# Labeled ladder diagram (rotation system) for a configuration.
function ladder_diagram(lad::Ladder, cfg)
    L = lad.L
    edges = Tuple{Int,Int,Int}[]
    for i in 1:L
        push!(edges, (i, mod1(i + 1, L), cfg[t_edge(lad, i)]))          # t_i: T_i--T_{i+1}
    end
    for i in 1:L
        push!(edges, (L + i, L + mod1(i + 1, L), cfg[b_edge(lad, i)]))  # b_i: B_i--B_{i+1}
    end
    for i in 1:L
        push!(edges, (i, L + i, cfg[v_edge(lad, i)]))                   # v_i: T_i--B_i
    end
    vertex_edges = Vector{Vector{Int}}(undef, 2L)
    for i in 1:L
        vertex_edges[i] = [v_edge(lad, i), t_edge(lad, i), t_edge(lad, i - 1)]
        vertex_edges[L+i] = [v_edge(lad, i), b_edge(lad, i - 1), b_edge(lad, i)]
    end
    build_diagram(2L, edges, vertex_edges)
end

# --------------------------------------------------------------------------
# The condensate wavefunction (vacuum sector).
# Returns (configs, amplitudes) with <Ψ|Ψ> = 1.
# --------------------------------------------------------------------------

function stringnet_ground_state(lad::Ladder; memo::Dict{String,Float64}=Dict{String,Float64}(),
                                chooser::Function=choose_fmove_edge)
    cfgs = ladder_configs(lad)
    raw = Vector{Float64}(undef, length(cfgs))
    template = ladder_diagram(lad, first(cfgs))
    for (idx, cfg) in enumerate(cfgs)
        D = copy(template)
        for e in eachindex(cfg)
            D.he_label[2e - 1] = D.he_label[2e] = cfg[e]
        end
        raw[idx] = _eval_diagram!(D, memo, chooser)
    end
    normalize!(raw)
    return cfgs, raw
end

# --------------------------------------------------------------------------
# Entanglement in the tensor product of edge-label Hilbert spaces (natural
# logarithm): A = {t_i, b_i, v_i : 1 <= i <= ell}. This is an edge partition,
# not a union of complete plaquettes. At ell=1 or L-1 its boundary vertices
# overlap, so an exact plateau across ALL ell is not an area-law criterion.
# --------------------------------------------------------------------------

function arc_entropy(cfgs, amps, lad::Ladder, ell::Int)
    0 <= ell <= lad.L || throw(ArgumentError("need 0 <= ell <= L"))
    length(cfgs) == length(amps) || throw(DimensionMismatch("configs and amplitudes"))
    isempty(amps) && throw(ArgumentError("empty state"))
    # Compact masks for each subsystem; BigInt avoids collisions past 64 bits.
    K = 3max(ell, lad.L - ell) <= 64 ? UInt64 : BigInt
    return _arc_entropy(cfgs, amps, lad, ell, K)
end

function _arc_entropy(cfgs, amps, lad, ell, ::Type{K}) where {K}
    L = lad.L
    aindex, bindex = Dict{K,Int}(), Dict{K,Int}()
    ri, ci = Int[], Int[]
    xs = typeof(float(zero(eltype(amps))))[]
    for (cfg, amp) in zip(cfgs, amps)
        length(cfg) == 3L || throw(DimensionMismatch("configuration needs 3L edges"))
        all(x -> x == 0 || x == 1, cfg) || throw(ArgumentError("labels must be 0 or 1"))
        iszero(amp) && continue
        a, b = zero(K), zero(K)
        for rail in 0:2
            for i in 1:ell
                a = (a << 1) | K(cfg[rail * L + i])
            end
            for i in (ell + 1):L
                b = (b << 1) | K(cfg[rail * L + i])
            end
        end
        push!(ri, get!(aindex, a, length(aindex) + 1))
        push!(ci, get!(bindex, b, length(bindex) + 1))
        push!(xs, amp)
    end
    na, nb = length(aindex), length(bindex)
    # Work on the smaller subsystem; transpose preserves singular values.
    if na > nb
        ri, ci = ci, ri
        na, nb = nb, na
    end
    C = sparse(ri, ci, xs, na, nb)
    dropzeros!(C)
    isapprox(sum(abs2, C), 1; atol=1e-10, rtol=1e-8) ||
        throw(ArgumentError("amplitudes must describe a normalized state"))

    # Connected components of the bipartite support give independent Schmidt
    # blocks. Join rows sharing a column, avoiding a global dense rho_A.
    parent = collect(1:na)
    for col in 1:nb
        entries = nzrange(C, col)
        isempty(entries) && continue
        root = _support_root!(parent, C.rowval[first(entries)])
        for k in entries
            parent[_support_root!(parent, C.rowval[k])] = root
        end
    end
    block_rows = [Int[] for _ in 1:na]
    block_cols = [Int[] for _ in 1:na]
    for row in 1:na
        push!(block_rows[_support_root!(parent, row)], row)
    end
    for col in 1:nb
        entries = nzrange(C, col)
        isempty(entries) && continue
        push!(block_cols[_support_root!(parent, C.rowval[first(entries)])], col)
    end
    entropy = 0.0
    for root in 1:na
        isempty(block_cols[root]) && continue
        block = Matrix(C[block_rows[root], block_cols[root]])
        for sigma in svdvals!(block)
            p = abs2(sigma)
            p > 0 && (entropy -= p * log(p))
        end
    end
    return max(0.0, entropy)
end

function _support_root!(parent, row)
    while parent[row] != row
        parent[row] = parent[parent[row]]
        row = parent[row]
    end
    return row
end
