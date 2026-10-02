using ITensors
using LinearAlgebra

# Keep the small-system diagram evaluator as an independent reference.
isdefined(@__MODULE__, :Diagram) || include(joinpath(@__DIR__, "stringnet.jl"))

"""
A patch of `Lx * Ly` hexagonal plaquettes in axial coordinates.
`bc` is `:open`, `:cylinder_x`, or `:cylinder_y`. Periodic lengths must be >= 2.
`edges[e]` gives its two endpoint vertex IDs; `edge_faces[e]` gives the faces
on its left/right when directed from its first endpoint to its second.
Face 0 is outside the patch and its auxiliary string label is fixed to vacuum.
`vertex_edges[v]` is counterclockwise, padded by a missing vacuum edge (0) at
boundary vertices. `vertex_faces[v] = (a,b,c)` denotes wedges (k,i), (i,j), (j,k).
Integer coordinates `(u,v)` embed as `(u/2, sqrt(3)*v/2)` before wrapping.
"""
struct HoneycombPatch
    Lx::Int
    Ly::Int
    bc::Symbol
    vertices::Vector{NTuple{2,Int}}
    edges::Vector{NTuple{2,Int}}
    edge_vectors::Vector{NTuple{2,Int}}
    edge_faces::Vector{NTuple{2,Int}}
    plaquettes::Vector{NTuple{6,Int}}
    vertex_edges::Vector{NTuple{3,Int}}
    vertex_faces::Vector{NTuple{3,Int}}
end

function _wrap_vertex(u, v, Lx, Ly, bc)
    if bc == :cylinder_x
        winding = fld(u, 3Lx)
        return (u - winding * 3Lx, v - winding * Lx)
    elseif bc == :cylinder_y
        return (u, mod(v, 2Ly))
    end
    return (u, v)
end

function _outgoing_faces(edges, edge_faces, e, v)
    left, right = edge_faces[e]
    return edges[e][1] == v ? (left, right) : (right, left)
end

function HoneycombPatch(Lx::Int, Ly::Int; bc::Symbol=:open)
    Lx >= 1 && Ly >= 1 || throw(ArgumentError("Lx and Ly must be positive"))
    bc in (:open, :cylinder_x, :cylinder_y) ||
        throw(ArgumentError("bc must be :open, :cylinder_x, or :cylinder_y"))
    (bc == :cylinder_x && Lx < 2 || bc == :cylinder_y && Ly < 2) &&
        throw(ArgumentError("the periodic direction needs at least two plaquettes"))
    corners = ((2,0), (1,1), (-1,1), (-2,0), (-1,-1), (1,-1))
    vertices = NTuple{2,Int}[]
    vertex_ids = Dict{NTuple{2,Int},Int}()
    edges, vectors, edge_faces = NTuple{2,Int}[], NTuple{2,Int}[], NTuple{2,Int}[]
    edge_ids = Dict{NTuple{4,Int},Int}()
    plaquettes = NTuple{6,Int}[]
    for y in 0:(Ly-1), x in 0:(Lx-1)
        face = length(plaquettes) + 1
        vs = ntuple(6) do k
            du, dv = corners[k]
            position = _wrap_vertex(3x + du, x + 2y + dv, Lx, Ly, bc)
            get!(vertex_ids, position) do
                push!(vertices, position)
                length(vertices)
            end
        end
        es = ntuple(6) do k
            next = mod1(k + 1, 6)
            u, v = vs[k], vs[next]
            du, dv = corners[next] .- corners[k]
            forward = (u, v, du, dv)
            backward = (v, u, -du, -dv)
            key = min(forward, backward)
            e = get!(edge_ids, key) do
                push!(edges, (key[1], key[2]))
                push!(vectors, (key[3], key[4]))
                push!(edge_faces, (0, 0))
                length(edges)
            end
            left, right = edge_faces[e]
            if key == forward
                left == 0 || error("overlapping plaquettes")
                edge_faces[e] = (face, right)
            else
                right == 0 || error("overlapping plaquettes")
                edge_faces[e] = (left, face)
            end
            e
        end
        push!(plaquettes, es)
    end
    incident = [Int[] for _ in vertices]
    for (e, (u,v)) in enumerate(edges)
        push!(incident[u], e)
        push!(incident[v], e)
    end
    vertex_edges = NTuple{3,Int}[]
    vertex_faces = NTuple{3,Int}[]
    for (v, es) in enumerate(incident)
        sort!(es; by=e -> begin
            du, dv = vectors[e]
            sign = edges[e][1] == v ? 1 : -1
            atan(sign * sqrt(3) * dv, sign * du)
        end)
        if length(es) == 2
            # Put the interior wedge between edges i and j; the absent k
            # points into the exterior and has fixed label 0.
            _outgoing_faces(edges, edge_faces, es[1], v)[1] == 0 && reverse!(es)
        end
        length(es) in (2,3) || error("patch is not trivalent with degree-two boundary")
        i, j = es[1:2]
        b, a = _outgoing_faces(edges, edge_faces, i, v)
        c, bcheck = _outgoing_faces(edges, edge_faces, j, v)
        b == bcheck || error("inconsistent rotation system")
        k = length(es) == 3 ? es[3] : 0
        if k == 0
            a == c == 0 || error("invalid boundary wedge")
        else
            _outgoing_faces(edges, edge_faces, k, v) == (a,c) || error("invalid third wedge")
        end
        push!(vertex_edges, (i,j,k))
        push!(vertex_faces, (a,b,c))
    end
    return HoneycombPatch(Lx, Ly, bc, vertices, edges, vectors, edge_faces,
                          plaquettes, vertex_edges, vertex_faces)
end

_quantum_dimension(i) = i == 0 ? 1.0 : φ

"""
Symmetric Fibonacci 6j-symbol in the convention with admissible triples
`(i,j,k)`, `(i,a,b)`, `(j,b,c)`, `(k,c,a)`:
`G(i,j,k,a,b,c) = F(i,j,c,a,k,b) / sqrt(d_k*d_b)`.
"""
fib_gsymbol(i, j, k, a, b, c) =
    fib_fsymbol(i,j,c,a,k,b) / sqrt(_quantum_dimension(k) * _quantum_dimension(b))

"""
Finite honeycomb PEPS with one physical qubit per original edge.
`tensors[v]` is an ITensor at honeycomb vertex v; `links[e]` connects exactly
its two endpoint tensors. Physical qubits are owned by the A sublattice
(u mod 3 == 2); B tensors have no physical leg (equivalently physical dimension 1).
Thus no physical qubit is duplicated. `physical_edges[v]` records the ownership;
`physical_sites[e]` is the corresponding dimension-two Index (0 -> 1, tau -> 2).
A bulk link has five basis states `(physical_label, left_face_label, right_face_label)`.
Boundary links have two after fixing the exterior face to vacuum.
"""
struct StringNetPEPS
    lattice::HoneycombPatch
    tensors::Vector{ITensor}
    physical_sites::Vector{Index{Int}}
    links::Vector{Index{Int}}
    bond_basis::Vector{Vector{NTuple{3,Int}}}
    physical_edges::Vector{Vector{Int}}
    normalized::Bool
end

"""
    fibonacci_stringnet_peps(Lx, Ly; bc=:open, normalize=true)

Construct an exact doubled-Fibonacci fixed-point PEPS on Lx by Ly hexagonal
plaquettes in O(Lx*Ly) time/storage, without enumerating physical configurations.
`bc=:cylinder_x` or `:cylinder_y` wraps only the named direction.

The state is proportional to `prod_p (B_p^0 + phi*B_p^tau) |vacuum>` with
exterior virtual face labels fixed to 0. On a cylinder this specifies one
closure, without a noncontractible MPO insertion; it is not a complete basis
of topological sectors. The squared norm before normalization is
`(1 + phi^2)^(Lx*Ly)`. Normalization is distributed among local tensors and
requires no network contraction.

The local weight is `(d_i*d_j*d_k)^(1/4) * (d_a*d_b*d_c)^(1/6) * G(i,j,k,a,b,c)`.
See Schotte et al., arXiv:1909.06284, Appendix A, Eqs. (20)-(21).
"""
function fibonacci_stringnet_peps(Lx::Int, Ly::Int; bc::Symbol=:open, normalize::Bool=true)
    lattice = HoneycombPatch(Lx, Ly; bc)
    basis = [NTuple{3,Int}[(i,a,b) for i in 0:1 for a in 0:1 for b in 0:1
             if fusion_allowed(i,a,b) && (left != 0 || a == 0) && (right != 0 || b == 0)]
             for (left,right) in lattice.edge_faces]
    links = [Index(length(b), "Link,StringNet,e=$e") for (e,b) in enumerate(basis)]
    physical = [Index(2, "Site,StringNet,e=$e") for e in eachindex(links)]
    owned = [mod(position[1],3) == 2 ? filter(!=(0), collect(es)) : Int[]
             for (position,es) in zip(lattice.vertices,lattice.vertex_edges)]
    tensors = Vector{ITensor}(undef, length(lattice.vertices))
    q = 1 + φ^2
    for v in eachindex(tensors)
        es = lattice.vertex_edges[v]
        fs = lattice.vertex_faces[v]
        actual_edges = filter(!=(0), collect(es))
        inds = vcat(links[actual_edges], physical[owned[v]])
        data = zeros(Float64, dim.(inds)...)
        nface_corners = count(!=(0), fs)
        scale = normalize ? q^(-nface_corners/12) : 1.0
        for i in 0:1, j in 0:1, k in 0:1, a in 0:1, b in 0:1, c in 0:1
            es[3] == 0 && k != 0 && continue
            ((fs[1] == 0 && a != 0) || (fs[2] == 0 && b != 0) ||
             (fs[3] == 0 && c != 0)) && continue
            g = fib_gsymbol(i,j,k,a,b,c)
            iszero(g) && continue
            labs = (i,j,k)
            pairs = ((a,b), (b,c), (c,a)) # right/left of outgoing half-edge
            positions = Int[]
            for (slot,e) in enumerate(actual_edges)
                right, left = pairs[slot]
                state = lattice.edges[e][1] == v ?
                    (labs[slot],left,right) : (labs[slot],right,left)
                idx = findfirst(==(state), basis[e])
                idx === nothing && error("inconsistent bond basis")
                push!(positions, idx)
            end
            isempty(owned[v]) || append!(positions, (labs[s]+1 for s in eachindex(actual_edges)))
            weight = (_quantum_dimension(i)*_quantum_dimension(j)*_quantum_dimension(k))^0.25 *
                     (_quantum_dimension(a)*_quantum_dimension(b)*_quantum_dimension(c))^(1/6)
            data[positions...] += scale * weight * g
        end
        tensors[v] = ITensor(data, inds...)
    end
    return StringNetPEPS(lattice,tensors,physical,links,basis,owned,normalize)
end

# A small, dependency-free greedy contraction for validation; construction
# never calls it. The explicit intermediate-size guard also applies to amplitudes.
function _contract_small(tensors; max_elements::Int)
    max_elements >= 1 || throw(ArgumentError("max_elements must be positive"))
    work = copy(tensors)
    while length(work) > 1
        best = nothing
        best_size = big(max_elements) + 1
        for i in 1:length(work)-1, j in i+1:length(work)
            common = commoninds(work[i],work[j])
            isempty(common) && continue
            output = uniqueinds(work[i],work[j])
            other = uniqueinds(work[j],work[i])
            size = prod(big(dim(k)) for k in output; init=big(1)) *
                   prod(big(dim(k)) for k in other; init=big(1))
            if size < best_size
                best_size, best = size, (i,j)
            end
        end
        best === nothing && throw(ArgumentError("contraction would exceed max_elements=$max_elements"))
        i,j = best
        work[i] = work[i] * work[j]
        deleteat!(work,j)
    end
    return only(work)
end

"""Contract a small PEPS to an ITensor with one open physical index per edge.
Exponential in patch size; use a PEPS contraction algorithm for larger systems.
"""
function contract_peps(psi::StringNetPEPS; max_elements::Int=1<<22)
    big(2)^length(psi.physical_sites) <= max_elements ||
        throw(ArgumentError("full state exceeds max_elements=$max_elements; contract selected amplitudes instead"))
    return _contract_small(psi.tensors; max_elements)
end

"""Contract a selected 0/1 edge configuration (ordered as `psi.lattice.edges`)."""
function peps_amplitude(psi::StringNetPEPS, config; max_elements::Int=1<<22)
    length(config) == length(psi.physical_sites) || throw(DimensionMismatch("one label per edge required"))
    all(x -> x == 0 || x == 1, config) || throw(ArgumentError("edge labels must be 0 or 1"))
    tensors = copy(psi.tensors)
    for v in eachindex(tensors), e in psi.physical_edges[v]
        tensors[v] *= onehot(psi.physical_sites[e] => config[e]+1)
    end
    return _contract_small(tensors; max_elements)[]
end
