using Random
using LinearAlgebra
isdefined(@__MODULE__, :StringNetPEPS) || include(joinpath(@__DIR__, "peps.jl"))

"""
    cft_trajectory_record(data; layers=:, columns=:, terminal_half_strength=nothing)

Adapt a dictionary returned by `JLD2.load` for monitored_dynamics_cft.jl.
Preserve raw Bool outcomes: false=tau, true=vacuum, so construct the PEPS map
with `outcome_labels=(1,0)`. `measurement_sites` locates events on the chain;
it does NOT assert a honeycomb embedding or identify an outcome with an edge.

The stored layer free energies give a marginal log probability ONLY for all
columns of a chronological prefix starting at layer 1. Crops/suffixes return
`log_probability=nothing`. Float32 storage and MPS truncation limit accuracy.
The producer does not save enable_tau_eff; pass terminal_half_strength=true
only after verifying that generating convention (the current default).
No coherent amplitude or final boundary state can be recovered from this data.
"""
function cft_trajectory_record(data::AbstractDict; layers=Colon(),columns=Colon(),
                               terminal_half_strength::Union{Nothing,Bool}=nothing)
    raw = data["sample"]
    raw isa AbstractMatrix || throw(ArgumentError("sample must be a layer-by-column matrix"))
    all(x -> x in (0,1),raw) || throw(ArgumentError("sample contains nonbinary outcomes"))
    L = Int(data["L"])
    L > 0 && iseven(L) && size(raw,2) == L÷2 || throw(DimensionMismatch("Fibonacci CFT record must have L/2 columns"))
    nr,nc = size(raw)
    if haskey(data,"periods")
        nr == 2Int(data["periods"]) || throw(DimensionMismatch("two layers per Fibonacci period required"))
    end
    rows = layers isa Colon ? collect(1:nr) : collect(Int,layers)
    cols = columns isa Colon ? collect(1:nc) : collect(Int,columns)
    for (ids,n) in ((rows,nr),(cols,nc))
        !isempty(ids) && all(i -> 1 <= i <= n,ids) && length(unique(ids)) == length(ids) ||
            throw(ArgumentError("layer/column selections must be nonempty, unique, and in bounds"))
    end
    prefix = rows == collect(1:length(rows)) && cols == collect(1:nc)
    logp = nothing
    if haskey(data,"sample_free_energy")
        fe = data["sample_free_energy"]
        fe isa AbstractVector && length(fe) == nr || throw(DimensionMismatch("one free energy per full layer required"))
        all(x -> isfinite(x) && x >= 0,fe) || throw(ArgumentError("layer free energies must be finite and nonnegative"))
        prefix && (logp = -sum(Float64,fe[rows]))
    end
    sites = [isodd(r) ? 2c : 2c-1 for r in rows,c in cols]
    tau = haskey(data,"tau") ? Float64(data["tau"]) : nothing
    strengths = tau === nothing || terminal_half_strength === nothing ? nothing :
        [terminal_half_strength && r == nr ? tau/2 : tau for r in rows]
    metadata = (; (Symbol(k) => data[k] for k in
        ("L","periods","seed","backend","initial_state","initial_state_method","tau","gamma","chi")
        if haskey(data,k))...)
    return (samples=BitMatrix(raw[rows,cols]),outcome_labels=(1,0),
            layers=rows,columns=cols,measurement_sites=sites,measurement_strengths=strengths,
            log_probability=logp,probability_scope=logp === nothing ? :unavailable : :prefix_marginal,
            phases=nothing,metadata=metadata)
end

# This workflow deliberately requires an explicit record-to-edge map. A circuit
# BitMatrix does not determine the geometric embedding or the outcome convention.
struct TrajectoryEdgeMap
    num_edges::Int
    edge_ids::Vector{Int}
    record_shape::Tuple{Vararg{Int}}
    fixed_labels::Vector{Int8} # -1 means unspecified, hence marginalized
    outcome_labels::NTuple{2,Int8}
end

"""
    TrajectoryEdgeMap(psi, record_edges; fixed_edges=Dict(), outcome_labels=(0,1))

`record_edges` has the same shape as one trajectory's record; its entries are
unique PEPS edge IDs. `outcome_labels` maps raw false/true (0/1) to physical
vacuum/tau labels. `fixed_edges` conditions the target on additional known
labels. All remaining edges are marginalized, not silently assigned vacuum.
Records may be arrays or objects with a `.samples` array. Julia's column-major
`vec` order is used for both the edge map and records.
"""
function TrajectoryEdgeMap(psi::StringNetPEPS, record_edges::AbstractArray{<:Integer};
                           fixed_edges=Dict{Int,Int}(), outcome_labels=(0,1))
    ne = length(psi.links)
    es = collect(Int, vec(record_edges))
    isempty(es) && throw(ArgumentError("record_edges must not be empty"))
    all(e -> 1 <= e <= ne, es) || throw(ArgumentError("record edge outside lattice"))
    length(unique(es)) == length(es) || throw(ArgumentError("record edge IDs must be unique"))
    Tuple(outcome_labels) in ((0,1),(1,0)) || throw(ArgumentError("outcome_labels must be (0,1) or (1,0)"))
    fixed = fill(Int8(-1), ne)
    for (e,s) in pairs(fixed_edges)
        e isa Integer && 1 <= e <= ne || throw(ArgumentError("fixed edge outside lattice"))
        s in (0,1) || throw(ArgumentError("fixed labels must be 0 or 1"))
        e in es && throw(ArgumentError("edge $e is both recorded and fixed"))
        fixed[e] = s
    end
    return TrajectoryEdgeMap(ne,es,size(record_edges),fixed,Int8.(Tuple(outcome_labels)))
end

function _sn_couplings(psi, edges, J)
    es = collect(Int, edges)
    all(e -> 1 <= e <= length(psi.links), es) || throw(ArgumentError("deformation edge outside lattice"))
    length(unique(es)) == length(es) || throw(ArgumentError("duplicate deformation edge"))
    js = J isa Real ? fill(Float64(J),length(es)) : Float64.(collect(J))
    length(js) == length(es) || throw(DimensionMismatch("one J per deformation edge required"))
    all(j -> isfinite(j) || j == -Inf, js) || throw(ArgumentError("J must be finite or -Inf"))
    couplings = zeros(length(psi.links))
    couplings[es] = js
    return couplings
end

"""Local filtered PEPS; `exp(log_scale) * peps` represents the literal D(J)|base>.
Finite positive J is stabilized by removing a common scalar from each filter.
`peps` is unnormalized. J=-Inf implements the exact tau projector.
"""
struct DeformedStringNetPEPS
    base::StringNetPEPS
    peps::StringNetPEPS
    couplings::Vector{Float64}
    log_scale::Float64
end

"""Apply diagonal physical filters diag(exp(J_e),1) on selected edges.
Accepts a scalar J or one J per edge; no global contraction or normalization.
"""
function deform_stringnet_peps(psi::StringNetPEPS; edges, J=-Inf)
    js = _sn_couplings(psi,edges,J)
    tensors = copy(psi.tensors)
    log_scale = 0.0
    for v in eachindex(tensors), e in psi.physical_edges[v]
        j = js[e]
        iszero(j) && continue
        shift = max(0.0,j)
        log_scale += shift
        s = psi.physical_sites[e]
        gate = ITensor([exp(j-shift) 0.0; 0.0 exp(-shift)],prime(s),s)
        tensors[v] = replaceind(tensors[v]*gate,prime(s),s)
    end
    filtered = StringNetPEPS(psi.lattice,tensors,psi.physical_sites,psi.links,
                             psi.bond_basis,psi.physical_edges,false)
    return DeformedStringNetPEPS(psi,filtered,js,log_scale)
end

"""Exact base amplitudes for small patches, in edge-bit order (edge 1 is LSB).
Amplitudes are normalized explicitly; `input_norm2` records the input norm.
Building this table is exponential and guarded by `max_elements`.
"""
struct ExactStringNetReference
    peps::StringNetPEPS
    amplitudes::Vector{ComplexF64}
    input_norm2::Float64
end

function exact_stringnet_reference(psi::StringNetPEPS; max_elements::Int=1<<20)
    dense = ITensors.@set_warn_order 64 contract_peps(psi;max_elements)
    a = ComplexF64.(vec(Array(dense,psi.physical_sites...)))
    nrm = norm(a)
    isfinite(nrm) && nrm > 0 || throw(ArgumentError("base PEPS has zero or nonfinite norm"))
    return ExactStringNetReference(psi,a/nrm,nrm^2)
end

_sn_logadd(a,b) = a == -Inf ? b : b == -Inf ? a : max(a,b)+log1p(exp(-abs(a-b)))
_sn_logsum(xs) = foldl(_sn_logadd,xs;init=-Inf)
_sn_bit(bits,e) = Int8((bits >> (e-1)) & 1)

function _sn_check_map(reference, mapping)
    mapping.num_edges == length(reference.peps.links) || throw(DimensionMismatch("map and reference lattice"))
end

"""
    stringnet_target_distribution(reference; mapping, edges=Int[], J=0)

Exact normalized distribution of recorded edges in D(J)|SN>, conditioned on
`mapping.fixed_labels`, tracing all other edges. No normalization over the
observed sample set is used. `deformation_norm_ratio` is <SN|D^2|SN>;
`joint_norm_ratio` includes additional fixed-edge conditioning. In the pure
projector case `p_zigzag_tau` equals the former, not a product of marginals.
Normalized complex amplitudes are returned only when no edge is unobserved.
"""
function stringnet_target_distribution(reference::ExactStringNetReference;
        mapping=TrajectoryEdgeMap(reference.peps,collect(eachindex(reference.peps.links))),
        edges=Int[], J=0.0)
    _sn_check_map(reference,mapping)
    js = _sn_couplings(reference.peps,edges,J)
    active = findall(!iszero,js)
    fixed = findall(>=(0),mapping.fixed_labels)
    unresolved = setdiff(collect(1:mapping.num_edges),union(mapping.edge_ids,fixed))
    logweights = fill(-Inf,length(reference.amplitudes))
    for idx in eachindex(logweights)
        a = reference.amplitudes[idx]
        iszero(a) && continue
        bits = idx-1
        value = 2log(abs(a))
        for e in active
            _sn_bit(bits,e) == 0 && (value += 2js[e])
        end
        logweights[idx] = value
    end
    logZ = _sn_logsum(logweights)
    logZ == -Inf && throw(DomainError(J,"deformation has zero norm"))
    isfinite(logZ) || throw(ArgumentError("deformation has zero norm or unrepresentable log norm"))
    for idx in eachindex(logweights)
        all(_sn_bit(idx-1,e) == mapping.fixed_labels[e] for e in fixed) || (logweights[idx] = -Inf)
    end
    logjoint = _sn_logsum(logweights)
    logjoint == -Inf && throw(DomainError(mapping.fixed_labels,"fixed-edge event has zero target probability"))
    isfinite(logjoint) || throw(ArgumentError("unrepresentable conditional log norm"))
    logbins = Dict{Tuple,Float64}()
    phases = Dict{Tuple,ComplexF64}()
    for idx in eachindex(logweights)
        lw = logweights[idx]
        isfinite(lw) || continue
        key = Tuple(_sn_bit(idx-1,e) for e in mapping.edge_ids)
        logbins[key] = _sn_logadd(get(logbins,key,-Inf),lw-logjoint)
        if isempty(unresolved)
            a = reference.amplitudes[idx]
            phases[key] = a/abs(a)
        end
    end
    probs = Dict(key => exp(value) for (key,value) in logbins)
    # Preserve log probabilities even when tiny probabilities underflow.
    amps = isempty(unresolved) ? Dict(key => sqrt(probs[key])*phases[key] for key in keys(probs)) : nothing
    strict = any(==(-Inf),js) && all(j -> j == 0 || j == -Inf,js)
    return (probabilities=probs, log_probabilities=logbins, amplitudes=amps,
            phases=isempty(unresolved) ? phases : nothing, mapping=mapping,
            couplings=js, unresolved_edges=unresolved,
            deformation_log_norm_ratio=logZ, deformation_norm_ratio=exp(logZ),
            joint_log_norm_ratio=logjoint, joint_norm_ratio=exp(logjoint),
            fixed_event_probability=exp(logjoint-logZ),
            p_zigzag_tau=strict ? exp(logZ) : nothing)
end

function _sn_record_key(trajectory,mapping)
    record = trajectory isa AbstractArray ? trajectory :
        hasproperty(trajectory,:samples) ? trajectory.samples :
        throw(ArgumentError("each trajectory must be an array or have a .samples array"))
    size(record) == mapping.record_shape || throw(DimensionMismatch("record and edge map shapes differ"))
    all(x -> x == 0 || x == 1,record) || throw(ArgumentError("records must contain 0/1 or Bool outcomes"))
    return Tuple(mapping.outcome_labels[Int(x)+1] for x in vec(record))
end

function _sn_branching_violation(key,mapping,lattice)
    cfg = copy(mapping.fixed_labels)
    cfg[mapping.edge_ids] = collect(key)
    for es in lattice.vertex_edges
        labs = ntuple(i -> es[i] == 0 ? Int8(0) : cfg[es[i]],3)
        all(>=(0),labs) && !fusion_allowed(labs...) && return true
    end
    return false
end

function _sn_metrics(counts,probs,n; logprobs=nothing)
    tv,js,kl,bc = 0.0,0.0,0.0,0.0
    for key in union(keys(counts),keys(probs))
        p,q = get(counts,key,0)/n,get(probs,key,0.0)
        tv += abs(p-q)/2
        m = (p+q)/2
        if p > 0
            js += p*log(p/m)/2
            logq = logprobs === nothing ? log(q) : get(logprobs,key,-Inf)
            kl += p*(log(p)-logq)
        end
        q > 0 && (js += q*log(q/m)/2)
        bc += sqrt(p*q)
    end
    return (total_variation=tv,jensen_shannon=js,kl_empirical_to_target=kl,
            classical_fidelity=bc^2)
end

# Finite-sample errors assume independent complete trajectories, not independent
# measurement events within a trajectory. Never reweight Born samples by p(s).
function _sn_mean_se(xs)
    mu = sum(xs)/length(xs)
    se = length(xs)>1 ? sqrt(sum(abs2(x-mu) for x in xs)/(length(xs)-1)/length(xs)) : NaN
    return (mean=mu,standard_error=se)
end

function _sn_importance(keys_sample,target,logp,phases; iid=true)
    n = length(keys_sample)
    length(logp) == n || throw(DimensionMismatch("one log p(s) per trajectory required"))
    all(x -> isfinite(x) && x <= 1e-10,logp) || throw(ArgumentError("log p(s) must be finite and <= 0"))
    phases !== nothing && target.amplitudes === nothing &&
        throw(ArgumentError("pure-state overlap requires all non-fixed edges to be recorded"))
    if phases !== nothing
        length(phases) == n || throw(DimensionMismatch("one phase per trajectory required"))
        all(x -> isfinite(x) && isapprox(abs(x),1;atol=1e-10),phases) ||
            throw(ArgumentError("phases must be unit-modulus amplitude factors, not phase angles"))
    end
    seen = Dict{Tuple,Int}()
    logr = Float64[]
    for (i,key) in enumerate(keys_sample)
        if haskey(seen,key)
            j = seen[key]
            isapprox(logp[i],logp[j];atol=1e-8,rtol=1e-8) ||
                throw(ArgumentError("duplicate records have inconsistent marginal log probabilities"))
            phases !== nothing && !isapprox(phases[i],phases[j];atol=1e-10) &&
                throw(ArgumentError("duplicate records have inconsistent amplitude phases"))
        end
        seen[key] = i
        push!(logr,get(target.log_probabilities,key,-Inf)-logp[i])
    end
    r = exp.(logr/2)
    all(isfinite,r) || throw(ArgumentError("importance weights overflow; the sampling distribution poorly covers the target"))
    overlap = _sn_mean_se(r)
    # U-statistic removes the square-of-sample-mean bias. It is not clipped.
    unbiased = n>1 ? (sum(r)^2-sum(abs2,r))/(n*(n-1)) : NaN
    w = abs2.(r)
    all(isfinite,w) || throw(ArgumentError("squared importance weights overflow"))
    mass = _sn_mean_se(w)
    scaled_w = iszero(maximum(w)) ? w : w/maximum(w)
    ess = iszero(sum(scaled_w)) ? 0.0 : sum(scaled_w)^2/sum(abs2,scaled_w)
    logp_over_q = -logr
    kl = all(isfinite,logp_over_q) ? _sn_mean_se(logp_over_q) : (mean=Inf,standard_error=NaN)
    quantum = nothing
    if phases !== nothing
        z = [iszero(r[i]) ? 0.0+0im : r[i]*conj(phases[i])*target.phases[keys_sample[i]] for i in 1:n]
        m = _sn_mean_se(z)
        quantum = (overlap=m.mean,overlap_standard_error=iid ? m.standard_error : nothing,
                   fidelity_plugin=abs2(m.mean),
                   fidelity_unbiased=iid && n>1 ? (abs2(sum(z))-sum(abs2,z))/(n*(n-1)) : nothing)
    end
    return (bhattacharyya=overlap.mean,bhattacharyya_standard_error=iid ? overlap.standard_error : nothing,
            classical_fidelity_plugin=overlap.mean^2,classical_fidelity_unbiased=iid ? unbiased : nothing,
            sampled_target_mass=mass.mean,sampled_target_mass_standard_error=iid ? mass.standard_error : nothing,
            importance_ess=ess,kl_p_to_q=kl.mean,kl_standard_error=iid ? kl.standard_error : nothing,
            quantum=quantum)
end

"""
    compare_trajectories(reference, trajectories; mapping, edges=Int[], J=0,
                        log_probabilities=nothing, phases=nothing,
                        bootstrap_replicates=199, iid=true, rng=Random.default_rng())

Compare equally weighted samples drawn from p(s) with the exact deformed PEPS
marginal q(s). Report empirical TV/JS/classical fidelity, unseen target mass,
forbidden records, and an IID Monte Carlo null calibration of TV. Finite-sample
histogram differences are not an exact-state fidelity or a proof of inequivalence.

Optional *normalized marginal* log p(s) gives importance estimates of
sum sqrt(p*q), KL(p||q), and q's mass on supp(p). Optional unit phases additionally
give a quantum overlap ONLY for complete records with scalar coherent amplitudes.
Do not supply conditional path probabilities for a many-to-one marginal mapping.
Do not infer phases from sqrt(p); a retained final chain generally leaves a mixed
record state. No quantum fidelity is reported from samples alone.
"""
function compare_trajectories(reference::ExactStringNetReference,trajectories;
        mapping=TrajectoryEdgeMap(reference.peps,collect(eachindex(reference.peps.links))),
        edges=Int[],J=0.0,log_probabilities=nothing,phases=nothing,
        bootstrap_replicates::Int=199,iid::Bool=true,rng=Random.default_rng())
    bootstrap_replicates >= 0 || throw(ArgumentError("bootstrap_replicates must be nonnegative"))
    phases !== nothing && log_probabilities === nothing &&
        throw(ArgumentError("phases require actual log p(s), not histogram estimates"))
    keys_sample = [_sn_record_key(s,mapping) for s in trajectories]
    n = length(keys_sample)
    n > 0 || throw(ArgumentError("at least one trajectory required"))
    target = stringnet_target_distribution(reference;mapping,edges,J)
    counts = Dict{Tuple,Int}()
    for key in keys_sample
        counts[key] = get(counts,key,0)+1
    end
    metrics = _sn_metrics(counts,target.probabilities,n;logprobs=target.log_probabilities)
    impossible = sum((count for (key,count) in counts if !haskey(target.log_probabilities,key));init=0)
    branching = sum((count for (key,count) in counts if _sn_branching_violation(key,mapping,reference.peps.lattice));init=0)
    unseen = sum((p for (key,p) in target.probabilities if !haskey(counts,key));init=0.0)
    null_p = nothing
    if iid && bootstrap_replicates > 0
        basis = collect(keys(target.probabilities))
        cdf = cumsum([target.probabilities[k] for k in basis]); cdf[end] = 1.0
        exceed = 0
        for _ in 1:bootstrap_replicates
            nullcounts = Dict{Tuple,Int}()
            for _ in 1:n
                key = basis[searchsortedfirst(cdf,rand(rng))]
                nullcounts[key] = get(nullcounts,key,0)+1
            end
            exceed += _sn_metrics(nullcounts,target.probabilities,n).total_variation >= metrics.total_variation
        end
        null_p = impossible>0 ? 0.0 : (exceed+1)/(bootstrap_replicates+1)
    end
    importance = log_probabilities === nothing ? nothing :
        _sn_importance(keys_sample,target,Float64.(log_probabilities),phases;iid)
    rows = [(edge_labels=key,count=get(counts,key,0),empirical_probability=get(counts,key,0)/n,
             target_probability=get(target.probabilities,key,0.0),
             log_target_probability=get(target.log_probabilities,key,-Inf))
            for key in sort!(collect(union(keys(counts),keys(target.probabilities))))]
    summary = (num_trajectories=n,num_distinct_records=length(counts),
               num_target_records=length(target.probabilities),branching_violations=branching,
               forbidden_target_records=impossible,unseen_target_mass=unseen,
               empirical=metrics,tv_null_pvalue=null_p,iid_assumed=iid,
               p_zigzag_tau=target.p_zigzag_tau,
               deformation_log_norm_ratio=target.deformation_log_norm_ratio,
               fixed_event_probability=target.fixed_event_probability,
               unresolved_edges=target.unresolved_edges,importance=importance)
    return (summary=summary,rows=rows,target=target)
end

"""Grid maximum-likelihood fit using exact, normalized deformed PEPS marginals.
Pass held-out trajectories to compare_trajectories after fitting; a null p-value
computed on the fitting sample without refitting is not a calibrated fit test.
A flat conditional likelihood (e.g. all tension edges fixed tau) returns
`identifiable=false, best_J=nothing` rather than claiming a finite tension.
"""
function fit_string_tension(reference::ExactStringNetReference,trajectories;
        mapping=TrajectoryEdgeMap(reference.peps,collect(eachindex(reference.peps.links))),
        edges,J_grid=range(-3.0,0.0;length=31))
    grid = Float64.(collect(J_grid))
    isempty(grid) && throw(ArgumentError("J_grid is empty"))
    ks = [_sn_record_key(s,mapping) for s in trajectories]
    isempty(ks) && throw(ArgumentError("at least one trajectory required"))
    counts = Dict{Tuple,Int}()
    for key in ks
        counts[key] = get(counts,key,0)+1
    end
    likelihoods = Float64[]
    for j in grid
        likelihood = try
            target = stringnet_target_distribution(reference;mapping,edges,J=j)
            sum(count*get(target.log_probabilities,k,-Inf) for (k,count) in counts)
        catch err
            err isa DomainError || rethrow()
            -Inf # zero-norm/zero-conditioning candidate, not a valid fit
        end
        push!(likelihoods,likelihood)
    end
    finite = filter(isfinite,likelihoods)
    compatible = !isempty(finite)
    identifiable = compatible && (length(finite)!=length(grid) || maximum(finite)-minimum(finite)>1e-10*length(ks))
    best = identifiable ? grid[argmax(likelihoods)] : nothing
    return (J_grid=grid,log_likelihood=likelihoods,best_J=best,
            identifiable=identifiable,compatible=compatible)
end

"""Exact pure-state comparison with a dictionary of complex amplitudes keyed by
*physical* record-label tuples. Omitted keys mean zero amplitude, not unobserved
support. This is a full-vector assertion by the caller, not a sample histogram.
"""
function compare_exact_amplitudes(reference::ExactStringNetReference,our_amplitudes::AbstractDict;
        mapping=TrajectoryEdgeMap(reference.peps,collect(eachindex(reference.peps.links))),
        edges=Int[],J=0.0)
    target = stringnet_target_distribution(reference;mapping,edges,J)
    target.amplitudes === nothing && throw(ArgumentError("unobserved edges prevent a pure record-state comparison"))
    for (key,a) in our_amplitudes
        length(key)==length(mapping.edge_ids) && all(x -> x in (0,1),key) ||
            throw(ArgumentError("amplitude keys must be physical record-label tuples"))
        isfinite(a) || throw(ArgumentError("nonfinite amplitude"))
    end
    norm_our = sum(abs2,values(our_amplitudes))
    isfinite(norm_our) && norm_our>0 || throw(ArgumentError("zero or nonfinite source norm"))
    overlap = sum((conj(a)*get(target.amplitudes,key,0.0) for (key,a) in our_amplitudes);init=0.0+0im)/sqrt(norm_our)
    classical = sum((abs(a)*sqrt(get(target.probabilities,key,0.0)) for (key,a) in our_amplitudes);init=0.0)/sqrt(norm_our)
    phase = iszero(overlap) ? 1.0+0im : conj(overlap)/abs(overlap)
    maxerror = maximum(abs(get(our_amplitudes,key,0.0)/sqrt(norm_our)-phase*get(target.amplitudes,key,0.0))
                       for key in union(keys(our_amplitudes),keys(target.amplitudes)))
    return (norm_our=norm_our,fidelity=abs2(overlap),classical_fidelity=classical^2,
            global_phase=phase,max_amplitude_error=maxerror,p_zigzag_tau=target.p_zigzag_tau)
end
