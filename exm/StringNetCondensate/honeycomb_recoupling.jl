using ITensors, LinearAlgebra
isdefined(@__MODULE__, :FibonacciCylinderReference) || include(joinpath(@__DIR__, "cylinder_transfer.jl"))

"""Honeycomb strip obtained by resolving each brick-wall event into two
trivalent G tensors. Space (L tau worldlines) is periodic; time is open.
All zigzag physical edges are pinned tau. The physical_sites matrix labels
rung/outcome qubits, with shape (num_layers,L/2), 0=vacuum and 1=tau.
Input/output links are OPEN virtual bonds with basis (left_face,right_face).
They must be contracted with the boundary maps below, not vacuum-capped.
"""
struct RecoupledHoneycombPEPS
    L::Int
    num_layers::Int
    tensors::Vector{ITensor}
    physical_sites::Matrix{Index{Int}}
    input_links::Vector{Index{Int}}
    output_links::Vector{Index{Int}}
    zigzag_basis::Vector{NTuple{2,Int}}
    rung_basis::Vector{NTuple{3,Int}}
    representation::Symbol
    tau::Float64
    total_layers::Int
    terminal_half_strength::Bool
end

"""Explicit embedding into HoneycombPatch(num_layers-2,L/2;bc=:cylinder_y).
Interior record (r,c) is the horizontal edge centered at u=3(r-2), v=i mod L,
where i=2c for odd r and 2c-1 for even r. The first/last rows are dangling
boundary rungs, marked 0 in record_edges, and MUST NOT be treated as bulk edges.
All slanted physical edges are the tau-pinned zigzag worldlines.
"""
function honeycomb_record_geometry(L::Int,num_layers::Int)
    L>=4 && iseven(L) && num_layers>=3 || throw(ArgumentError("even L>=4 and at least 3 layers required"))
    lat = HoneycombPatch(num_layers-2,L÷2;bc=:cylinder_y)
    lookup = Dict{Tuple{Int,Int},Int}()
    zigzag = Int[]
    for (e,(a,b)) in enumerate(lat.edges)
        if lat.edge_vectors[e][2]==0
            u = (lat.vertices[a][1]+lat.vertices[b][1])÷2
            lookup[(u,lat.vertices[a][2])] = e
        else
            push!(zigzag,e)
        end
    end
    records = zeros(Int,num_layers,L÷2)
    for r in 2:num_layers-1,c in 1:L÷2
        i = isodd(r) ? 2c : 2c-1
        records[r,c] = lookup[(3(r-2),mod(i,L))]
    end
    boundary_vertices = [findfirst(==((u,mod(isodd(r) ? 2c : 2c-1,L))),lat.vertices)
                         for (r,u) in ((1,-2),(num_layers,3(num_layers-2)-1)),c in 1:L÷2]
    all(!isnothing,boundary_vertices) || error("invalid boundary embedding")
    return (lattice=lat,record_edges=records,zigzag_edges=zigzag,
            input_vertices=Int.(boundary_vertices[1,:]),output_vertices=Int.(boundary_vertices[2,:]))
end

"""Single fusion-path vector that reproduces the original vacuum-capped
PEPS at this end after fixing the adjacent boundary rung outcomes to vacuum.
This is generally NOT the TCI boundary state and need not be a Y eigenstate.
"""
function honeycomb_vacuum_cap(reference::FibonacciCylinderReference,layer::Int)
    target = Int8[isodd(i+layer) ? 1 : 0 for i in 1:reference.L]
    k = findfirst(k -> reference.labels[k,:]==target,axes(reference.labels,1))
    k === nothing && error("missing alternating fusion path")
    state = zeros(ComplexF64,size(reference.labels,1)); state[k] = 1
    return state
end

# G(τ,τ,s;d,b,a) = F(a,τ,τ;d,b,s)/sqrt(d_b*d_s).
# Use the actual symmetric-G vertex, not a product of measurement matrices.
function _honeycomb_vertex(s,a,b,d)
    ds,da,db,dd = _quantum_dimension.((s,a,b,d))
    return (φ^2*ds)^0.25*(da*db*dd)^(1/6)*fib_gsymbol(1,1,s,d,b,a)
end

"""
    recoupled_honeycomb_peps(L, num_layers; representation=:monitored,
                             tau=Inf, total_layers=num_layers,
                             terminal_half_strength=true)

Build O(L*num_layers) local honeycomb PEPS tensors, without a dense state.
:stringnet uses the unmodified standard G vertices with zigzag=tau.
:monitored applies R(tau)*diag(1,sqrt(phi))/phi on every physical rung.
At tau=Inf this is a diagonal rung tension J=-log(phi)/2, up to a scalar.
These are physical filters, not virtual gauge transformations.
The first and last rung rows lie on the cut boundaries of the usual
HoneycombPatch(num_layers-2,L/2;bc=:cylinder_y), for num_layers>=3.
"""
function recoupled_honeycomb_peps(L::Int,num_layers::Int;
        representation::Symbol=:monitored,tau::Real=Inf,
        total_layers::Int=num_layers,terminal_half_strength::Bool=true)
    L>=4 && iseven(L) || throw(ArgumentError("even L>=4 required"))
    1<=num_layers<=total_layers || throw(ArgumentError("invalid layer count"))
    representation in (:stringnet,:monitored) || throw(ArgumentError("unknown representation"))
    fibonacci_outcome_filter(tau)
    zb = [(a,b) for a in 0:1 for b in 0:1 if fusion_allowed(1,a,b)]
    rb = [(s,a,d) for s in 0:1 for a in 0:1 for d in 0:1 if fusion_allowed(s,a,d)]
    inputs = [Index(length(zb),"Link,Honeycomb,input=$i") for i in 1:L]
    wires = copy(inputs)
    sites = [Index(2,"Site,Honeycomb,row=$r,col=$c") for r in 1:num_layers,c in 1:L÷2]
    tensors = ITensor[]
    for r in 1:num_layers
        outputs = [Index(length(zb),"Link,Honeycomb,row=$r,wire=$i") for i in 1:L]
        t = terminal_half_strength && r==total_layers ? tau/2 : tau
        filter = representation == :stringnet ? Matrix{Float64}(I,2,2) :
                 fibonacci_outcome_filter(t)*Diagonal([1.0,sqrt(φ)])/φ
        for c in 1:L÷2
            i = isodd(r) ? 2c : 2c-1
            j = mod1(i+1,L)
            rung = Index(length(rb),"Link,Honeycomb,rung=$r,$c")
            lower = zeros(length(zb),length(zb),length(rb),2)
            upper = zeros(length(zb),length(zb),length(rb))
            for (k,(s,a,d)) in enumerate(rb), b in 0:1
                left,right = findfirst(==((a,b)),zb),findfirst(==((b,d)),zb)
                (left === nothing || right === nothing) && continue
                w = _honeycomb_vertex(s,a,b,d)
                upper[left,right,k] = w
                for outcome in 0:1
                    lower[left,right,k,outcome+1] = w*filter[outcome+1,s+1]
                end
            end
            push!(tensors,ITensor(lower,wires[i],wires[j],rung,sites[r,c]))
            push!(tensors,ITensor(upper,outputs[i],outputs[j],rung))
        end
        wires = outputs
    end
    return RecoupledHoneycombPEPS(L,num_layers,tensors,sites,inputs,wires,zb,rb,
                                 representation,Float64(tau),total_layers,terminal_half_strength)
end

"""q_r(x)=prod_unmeasured d(x)^(1/3) prod_measured d(x)^(-1/3).
The paired G transfer is C_s Q_r P_layer(s) Q_r; adjacent Q cancel because
odd/even layers alternate. Q_1^{-1} and Q_T^{-1} are the boundary maps.
"""
function honeycomb_boundary_factors(reference::FibonacciCylinderReference,layer::Int)
    layer>=1 || throw(ArgumentError("layer must be positive"))
    return [prod(_quantum_dimension(reference.labels[k,i])^(iseven(i+layer) ? 1/3 : -1/3)
                 for i in 1:reference.L) for k in axes(reference.labels,1)]
end

"""Y acting in the boundary coordinates of the G-tensor PEPS. The input
encoding is Q_1^{-1}, while the undecoded output is Q_T. The similarity-
transformed loop is Hermitian in the returned metric, not generally in the
Euclidean metric. This realizes the chain's two flux projectors on the actual
PEPS virtual cut; it does not resolve all four bulk tube-algebra sectors.
"""
function honeycomb_boundary_symmetry(reference::FibonacciCylinderReference,layer::Int;
        side::Symbol=:output)
    side in (:input,:output) || throw(ArgumentError("invalid boundary side"))
    q = honeycomb_boundary_factors(reference,layer)
    scale = side == :input ? inv.(q) : q
    loop = Diagonal(scale)*reference.Y*Diagonal(inv.(scale))
    return (loop=loop,metric=Diagonal(inv.(scale).^2),
            vacuum_flux=(loop+I/φ)/(φ+1/φ),tau_flux=(φ*I-loop)/(φ+1/φ))
end

function _honeycomb_boundary_positions(psi,reference,k)
    return [something(findfirst(==((Int(reference.labels[k,mod1(i-1,psi.L)]),
                                    Int(reference.labels[k,i]))),psi.zigzag_basis)) for i in 1:psi.L]
end

"""Encode a fusion-path boundary vector as an open-virtual-boundary ITensor.
The returned tensor includes Q_layer^{-1}; output bras are conjugated.
No normalization is applied to the supplied vector, preserving its amplitude.
"""
function honeycomb_boundary_tensor(psi::RecoupledHoneycombPEPS,reference::FibonacciCylinderReference,
        state=reference.initial_state;side::Symbol=:input,max_elements::Int=1<<22)
    psi.L==reference.L || throw(DimensionMismatch("cylinder circumference"))
    side in (:input,:output) || throw(ArgumentError("side must be :input or :output"))
    length(state)==size(reference.labels,1) || throw(DimensionMismatch("fusion boundary vector"))
    all(isfinite,state) || throw(ArgumentError("nonfinite boundary amplitude"))
    big(length(psi.zigzag_basis))^psi.L<=max_elements || throw(ArgumentError("boundary tensor exceeds max_elements"))
    links = side == :input ? psi.input_links : psi.output_links
    q = honeycomb_boundary_factors(reference,side == :input ? 1 : psi.num_layers)
    data = zeros(ComplexF64,dim.(links)...)
    for k in eachindex(state)
        data[_honeycomb_boundary_positions(psi,reference,k)...] = (side == :input ? state[k] : conj(state[k]))/q[k]
    end
    return ITensor(data,links...)
end

"""Contract a small open-boundary honeycomb PEPS to a matrix of coherent
amplitudes (final fusion path, physical record). The output decoding includes
Q_T^{-1}; therefore its Euclidean norm is the physical Born norm. With a
final_state, also return its unnormalized boundary-postselected amplitudes.
The :stringnet result keeps its physical rung weights; :monitored must match
cylinder_record_state before normalization, including phase.
"""
function contract_recoupled_honeycomb(psi::RecoupledHoneycombPEPS,
        reference::FibonacciCylinderReference;final_state=nothing,max_elements::Int=1<<22)
    big(length(psi.zigzag_basis))^psi.L*big(2)^length(psi.physical_sites)<=max_elements ||
        throw(ArgumentError("open-boundary output exceeds max_elements"))
    tensors = vcat(psi.tensors,[honeycomb_boundary_tensor(psi,reference;max_elements)])
    tensor = ITensors.@set_warn_order 64 _contract_small(tensors;max_elements)
    data = Array(tensor,psi.output_links...,vec(psi.physical_sites)...)
    flat = reshape(data,3^psi.L,:)
    linear = LinearIndices(Tuple(fill(3,psi.L)))
    q = honeycomb_boundary_factors(reference,psi.num_layers)
    W = Matrix{ComplexF64}(undef,size(reference.labels,1),size(flat,2))
    for k in axes(W,1)
        index = linear[_honeycomb_boundary_positions(psi,reference,k)...]
        W[k,:] = flat[index,:]/q[k]
    end
    amplitudes = final_state === nothing ? nothing : vec(final_state'*W)
    return (joint_amplitudes=W,probabilities=vec(sum(abs2,W;dims=1)),
            norm2=sum(abs2,W),postselected_amplitudes=amplitudes,
            record_shape=size(psi.physical_sites))
end

"""Contract one raw false=tau/true=vacuum record. This avoids enumerating
all outcome strings. The result retains the final chain, or contracts a
specified final boundary. It uses actual local G tensors and exact contraction.
"""
function recoupled_honeycomb_record(psi::RecoupledHoneycombPEPS,
        reference::FibonacciCylinderReference,record::AbstractMatrix;
        final_state=nothing,max_elements::Int=1<<22)
    size(record)==size(psi.physical_sites) || throw(DimensionMismatch("record shape"))
    all(x -> x in (0,1),record) || throw(ArgumentError("nonbinary record"))
    tensors = copy(psi.tensors)
    for r in 1:psi.num_layers,c in 1:psi.L÷2
        k = 2*((r-1)*(psi.L÷2)+c)-1
        tensors[k] *= onehot(psi.physical_sites[r,c] => 2-Int(record[r,c]))
    end
    push!(tensors,honeycomb_boundary_tensor(psi,reference;max_elements))
    if final_state !== nothing
        push!(tensors,honeycomb_boundary_tensor(psi,reference,final_state;side=:output,max_elements))
        return (amplitude=_contract_small(tensors;max_elements)[],final_vector=nothing)
    end
    tensor = _contract_small(tensors;max_elements)
    q = honeycomb_boundary_factors(reference,psi.num_layers)
    state = [tensor[(psi.output_links .=> _honeycomb_boundary_positions(psi,reference,k))...]/q[k]
             for k in axes(reference.labels,1)]
    return (amplitude=nothing,final_vector=state)
end

"""Cached sparse paired-G transfers for many trajectories at a fixed width.
Matrices are built from the two local G vertices, independently of the circuit
projectors. The cache is local to one reference and not thread-safe during its
first fill. State-vector memory scales with the fusion basis, not 2^(L*T/2).
"""
mutable struct HoneycombTransferCache
    reference::FibonacciCylinderReference
    representation::Symbol
    operators::Dict{Tuple{Int,Float64,Bool},SparseMatrixCSC{Float64,Int}}
    positions::Dict{Tuple,Int}
end

function HoneycombTransferCache(reference::FibonacciCylinderReference;representation::Symbol=:monitored)
    representation in (:stringnet,:monitored) || throw(ArgumentError("unknown representation"))
    return HoneycombTransferCache(reference,representation,
        Dict{Tuple{Int,Float64,Bool},SparseMatrixCSC{Float64,Int}}(),
        Dict(Tuple(reference.labels[k,:])=>k for k in axes(reference.labels,1)))
end

function _honeycomb_transfer(cache::HoneycombTransferCache,i,t,raw)
    return get!(cache.operators,(i,Float64(t),Bool(raw))) do
        ref = cache.reference
        filter = cache.representation == :stringnet ? Matrix{Float64}(I,2,2) :
                 fibonacci_outcome_filter(t)*Diagonal([1.0,sqrt(φ)])/φ
        rows,cols,values = Int[],Int[],Float64[]
        for k in axes(ref.labels,1)
            cfg = collect(ref.labels[k,:])
            a,b,d = cfg[mod1(i-1,ref.L)],cfg[i],cfg[mod1(i+1,ref.L)]
            for bp in 0:1
                value = sum(filter[2-Int(raw),s+1]*_honeycomb_vertex(s,a,b,d)*
                            _honeycomb_vertex(s,a,bp,d) for s in 0:1)
                iszero(value) && continue
                cfg[i] = bp
                push!(rows,cache.positions[Tuple(cfg)]); push!(cols,k); push!(values,value)
            end
        end
        n = size(ref.labels,1)
        sparse(rows,cols,values,n,n)
    end
end

"""Log Born weight from a row contraction of the local honeycomb G tensors.
For representation=:monitored it is normalized p(s). For :stringnet it is an
UNNORMALIZED target weight; do not use it as log q(s) without its partition sum.
Optional final_state returns the normalized-final-vector boundary overlap;
the scalar amplitude is exp(log_weight/2)*boundary_overlap.
The initial/final Q factors and the last-layer strength convention are explicit.
"""
function honeycomb_record_probability(cache::HoneycombTransferCache,record::AbstractMatrix;
        total_layers::Int=size(record,1),terminal_half_strength::Bool=true,final_state=nothing)
    ref = cache.reference
    nr,nc = size(record)
    nc==ref.L÷2 && 1<=nr<=total_layers || throw(DimensionMismatch("complete chronological prefix required"))
    all(x -> x in (0,1),record) || throw(ArgumentError("nonbinary outcome"))
    state = ref.initial_state ./ honeycomb_boundary_factors(ref,1)
    nrm = norm(state)
    logweight = 2log(nrm)
    state ./= nrm
    buffer = similar(state)
    for r in 1:nr
        t = terminal_half_strength && r==total_layers ? ref.tau/2 : ref.tau
        for c in 1:nc
            i = isodd(r) ? 2c : 2c-1
            mul!(buffer,_honeycomb_transfer(cache,i,t,record[r,c]),state)
            nrm = norm(buffer)
            iszero(nrm) && return (log_weight=-Inf,normalized_probability=cache.representation==:monitored,
                                   final_state=nothing,boundary_overlap=nothing)
            logweight += 2log(nrm)
            buffer ./= nrm
            state,buffer = buffer,state
        end
    end
    state ./= honeycomb_boundary_factors(ref,nr)
    nrm = norm(state)
    logweight += 2log(nrm)
    state ./= nrm
    overlap = nothing
    if final_state !== nothing
        length(final_state)==length(state) || throw(DimensionMismatch("final boundary"))
        all(isfinite,final_state) || throw(ArgumentError("nonfinite final boundary"))
        overlap = dot(final_state,state)
    end
    return (log_weight=logweight,normalized_probability=cache.representation==:monitored,
            final_state=state,boundary_overlap=overlap)
end

"""Compare a batch of saved CFT trajectories with the boundary-matched
honeycomb PEPS. The cache is reused across records. Float32 saved free energies
and MPS boundary/truncation errors limit comparison to stored log p(s).
"""
function compare_honeycomb_trajectories(cache::HoneycombTransferCache,datasets;
        terminal_half_strength::Bool=true)
    cache.representation==:monitored || throw(ArgumentError("batch probability audit requires the normalized monitored representation"))
    rows = NamedTuple[]
    for data in datasets
        data["L"]==cache.reference.L && isapprox(data["tau"],cache.reference.tau) ||
            throw(ArgumentError("mixed L/tau ensemble"))
        record = cft_trajectory_record(data;terminal_half_strength)
        replay = honeycomb_record_probability(cache,record.samples;terminal_half_strength)
        saved = record.log_probability
        delta = saved === nothing ? nothing : replay.log_weight-saved
        push!(rows,(seed=get(data,"seed",nothing),stored_log_probability=saved,
                    peps_log_probability=replay.log_weight,difference=delta))
    end
    return rows
end
