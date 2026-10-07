using LinearAlgebra, SparseArrays
import FibonacciChain
isdefined(@__MODULE__, :ExactStringNetReference) || include(joinpath(@__DIR__, "trajectory_equivalence.jl"))

"""Small-width cylindrical F-symbol transfer reference. Space is periodic;
time is open. The initial vector is in FibonacciChain.anyon_basis order.
The final fusion-path index is retained unless a final vector is supplied.
The physical-edge embedding and the open-boundary triple-line PEPS realization
are implemented independently in honeycomb_recoupling.jl.
"""
struct FibonacciCylinderReference
    L::Int
    labels::Matrix{Int8}
    vacuum_projectors::Vector{SparseMatrixCSC{Float64,Int}}
    Y::Matrix{Float64}
    initial_state::Vector{ComplexF64}
    tau::Float64
    initial_energy::Float64
    sector::Symbol
end

"""The weak-measurement filter on the projective OUTCOME index (vacuum,tau).
R = [1 exp(-tau); exp(-tau) 1]/sqrt(1+exp(-2tau)). This filter is not diagonal
string tension in that basis. It becomes the identity in the projective limit.
"""
function fibonacci_outcome_filter(tau::Real)
    tau >= 0 && !isnan(tau) || throw(ArgumentError("tau must be nonnegative"))
    r = exp(-Float64(tau))
    return [1.0 r; r 1.0]/sqrt(1+r*r)
end

"""Diagonalization of the outcome filter, not an edge-basis identification.
In the (antisymmetric,symmetric) outcome basis U, R=scale*U*diag(exp(J),1)*U'.
For tau>0, J=log(tanh(tau/2)); at tau=0, J=-Inf. The rotated projective state
must also be transformed; raw Born records cannot simply be relabeled into it.
"""
function fibonacci_outcome_deformation(tau::Real)
    R = fibonacci_outcome_filter(tau)
    U = [1.0 1.0; -1.0 1.0]/sqrt(2)
    scale = R[1,1]+R[1,2]
    J = log(tanh(Float64(tau)/2))
    return (rotation=U,J=J,scale=scale,filter=R)
end

"""Spectral projectors of the boundary Fibonacci loop, Y^2=I+Y.
These resolve two fusion-algebra eigenvalues, not all doubled-Fibonacci tube
sectors. They are not assumed to be logical Pauli projectors.
"""
function fibonacci_flux_projectors(Y::AbstractMatrix;atol=1e-10)
    size(Y,1)==size(Y,2) || throw(DimensionMismatch("Y must be square"))
    isapprox(Y,Y';atol) && isapprox(Y*Y,I+Y;atol) || throw(ArgumentError("Y must be a Hermitian Fibonacci loop operator"))
    return (vacuum_flux=(Y+I/φ)/(φ+1/φ),tau_flux=(φ*I-Y)/(φ+1/φ))
end

function cylinder_sector_weights(reference::FibonacciCylinderReference,state=reference.initial_state)
    length(state)==size(reference.Y,1) || throw(DimensionMismatch("state and fusion basis"))
    n = real(dot(state,state))
    n>0 && isfinite(n) || throw(ArgumentError("state must have finite positive norm"))
    y = real(dot(state,reference.Y*state))/n
    return (y_expectation=y,vacuum_flux=(y+1/φ)/(φ+1/φ),tau_flux=(φ-y)/(φ+1/φ))
end

"""
    fibonacci_cylinder_reference(L; tau, initial_state=nothing,
                                 sector=:unrestricted, max_basis=512)

Construct each local vacuum-fusion projector independently from fib_fsymbol:
P_i(b',b) = F(a,tau,tau,d,b',vac) F(a,tau,tau,d,b,vac).
Build the periodic boundary Y from a closed product of the same F symbols.
Default input is the exact TCI ground state of the repository Hamiltonian;
sector=:vacuum_flux or :tau_flux selects its lowest state within that sector.
A supplied initial_state is projected and normalized in the requested sector.
Dense ED is limited by max_basis; no large-width or HPC calculation is launched.
"""
function fibonacci_cylinder_reference(L::Int;tau::Real,initial_state=nothing,
                                      sector::Symbol=:unrestricted,max_basis::Int=512)
    L>=4 && iseven(L) || throw(ArgumentError("even periodic L>=4 required"))
    sector in (:unrestricted,:vacuum_flux,:tau_flux) || throw(ArgumentError("unknown boundary sector"))
    fibonacci_outcome_filter(tau)
    # Lucas-number dimension, checked BEFORE enumerating the basis.
    a,b = big(2),big(1)
    for _ in 2:L
        a,b = b,a+b
        b<=max_basis || throw(ArgumentError("fusion basis exceeds max_basis=$max_basis"))
    end
    model = FibonacciChain.AnyonModel(FibonacciChain.FibonacciAnyon(),L;pbc=true)
    basis = FibonacciChain.anyon_basis(model)
    n = length(basis)
    labels = Int8[1-Int(basis[k][L-i+1]) for k in 1:n,i in 1:L]
    indices = Dict(Tuple(labels[k,:])=>k for k in 1:n)
    ps = SparseMatrixCSC{Float64,Int}[]
    for i in 1:L
        rows,cols,values = Int[],Int[],Float64[]
        for k in 1:n
            config = collect(labels[k,:])
            a,d,b = config[mod1(i-1,L)],config[mod1(i+1,L)],config[i]
            f = fib_fsymbol(a,1,1,d,b,0)
            iszero(f) && continue
            for bp in 0:1
                weight = f*fib_fsymbol(a,1,1,d,bp,0)
                iszero(weight) && continue
                config[i] = bp
                push!(rows,indices[Tuple(config)]); push!(cols,k); push!(values,weight)
            end
        end
        push!(ps,sparse(rows,cols,values,n,n))
    end
    Y = [prod(fib_fsymbol(1,labels[r,i],1,labels[c,mod1(i+1,L)],
                         labels[r,mod1(i+1,L)],labels[c,i]) for i in 1:L) for r in 1:n,c in 1:n]
    projectors = fibonacci_flux_projectors(Y)
    H = FibonacciChain.anyon_ham(model)
    if initial_state === nothing
        if sector == :unrestricted
            eig = eigen(Hermitian(H))
            state = ComplexF64.(eig.vectors[:,1])
        else
            p = getproperty(projectors,sector)
            ep = eigen(Hermitian(p))
            U = ep.vectors[:,findall(>(0.5),ep.values)]
            state = ComplexF64.(U*eigen(Hermitian(U'*H*U)).vectors[:,1])
        end
    else
        length(initial_state)==n || throw(DimensionMismatch("initial state and fusion basis"))
        state = ComplexF64.(initial_state)
        sector == :unrestricted || (state = getproperty(projectors,sector)*state)
    end
    norm(state)>1e-13 && isfinite(norm(state)) || throw(ArgumentError("initial state has zero sector weight or invalid norm"))
    normalize!(state)
    return FibonacciCylinderReference(L,labels,ps,Y,state,Float64(tau),real(dot(state,H*state)),sector)
end

"""Replay a complete record or chronological prefix with independent F tensors.
Set total_layers to the ORIGINAL record length when replaying a prefix:
the tau/2 correction is applied only on that original terminal layer.
Return log Born probability, normalized final chain, and its Y-sector weights.
Optional final_state gives the scalar boundary amplitude divided by sqrt(p);
its unnormalized log magnitude is log_probability/2 + log(abs(boundary_overlap)).
"""
function cylinder_record_probability(reference::FibonacciCylinderReference,samples::AbstractMatrix;
        total_layers::Int=size(samples,1),terminal_half_strength::Bool=true,final_state=nothing)
    nr,nc = size(samples)
    nc==reference.L÷2 && 1<=nr<=total_layers || throw(DimensionMismatch("complete layers of a chronological prefix required"))
    all(x -> x in (0,1),samples) || throw(ArgumentError("nonbinary outcome"))
    state = copy(reference.initial_state)
    logp = 0.0
    layers = Float64[]
    for r in 1:nr
        t = terminal_half_strength && r==total_layers ? reference.tau/2 : reference.tau
        filter = fibonacci_outcome_filter(t)
        high,low = filter[1,1],filter[1,2]
        layer_logp = 0.0
        for c in 1:nc
            i = isodd(r) ? 2c : 2c-1
            pstate = reference.vacuum_projectors[i]*state
            state = samples[r,c]==1 ? low*state+(high-low)*pstate : high*state+(low-high)*pstate
            nrm = norm(state)
            if iszero(nrm)
                return (log_probability=-Inf,layer_log_probabilities=vcat(layers,-Inf),
                        final_state=nothing,sector_weights=nothing,boundary_overlap=nothing)
            end
            layer_logp += 2log(nrm)
            state ./= nrm
        end
        logp += layer_logp
        push!(layers,layer_logp)
    end
    boundary_overlap = nothing
    if final_state !== nothing
        length(final_state)==length(state) || throw(DimensionMismatch("final boundary vector"))
        isfinite(norm(final_state)) && norm(final_state)>0 || throw(ArgumentError("invalid final boundary"))
        boundary_overlap = dot(final_state,state)/norm(final_state)
    end
    return (log_probability=logp,layer_log_probabilities=layers,final_state=state,
            sector_weights=cylinder_sector_weights(reference,state),boundary_overlap=boundary_overlap)
end

"""Audit JLD2.load dictionaries from one CFT ensemble against the F network.
Checks L, tau, and default TCI boundary metadata, and reuses the same reference.
No Born probabilities are inferred for arbitrary cropped records.
"""
function audit_cft_trajectories(reference::FibonacciCylinderReference,datasets;
        terminal_half_strength::Bool=true)
    rows = NamedTuple[]
    for data in datasets
        data["L"]==reference.L && isapprox(data["tau"],reference.tau) || throw(ArgumentError("mixed L/tau ensembles"))
        get(data,"initial_state",nothing)=="TCI_GS" || throw(ArgumentError("expected TCI_GS data"))
        record = cft_trajectory_record(data;terminal_half_strength)
        result = cylinder_record_probability(reference,record.samples;terminal_half_strength)
        stored = record.log_probability
        push!(rows,(seed=get(data,"seed",nothing),stored_log_probability=stored,
            f_network_log_probability=result.log_probability,
            log_probability_difference=stored === nothing ? nothing : result.log_probability-stored,
            sector_weights=result.sector_weights))
    end
    return rows
end

"""Enumerate the coherent periodic-circuit state for a small number of layers.
Columns use physical outcome labels 0=vacuum, 1=tau, with edge bits in Julia
vec(record) order. Rows retain the final fusion path. transpose(W'W) is
the reduced record density matrix; its diagonal is Born p(s), not its square
root wavefunction. `record_purity` uses the smaller boundary density matrix.
With final_state, additionally return normalized scalar record amplitudes and
the total postselection probability. Enumeration is guarded by max_elements.
"""
function cylinder_record_state(reference::FibonacciCylinderReference,num_layers::Int;
        total_layers::Int=num_layers,terminal_half_strength::Bool=true,
        final_state=nothing,max_elements::Int=1<<20)
    1<=num_layers<=total_layers || throw(ArgumentError("invalid layer count"))
    nc = reference.L÷2
    events = num_layers*nc
    n = length(reference.initial_state)
    events < 8sizeof(Int)-2 && (big(n)<<events)<=max_elements ||
        throw(ArgumentError("coherent record table exceeds max_elements=$max_elements"))
    W = zeros(ComplexF64,n,1<<events)
    function visit(state,event,bits)
        if event>events
            W[:,bits+1] = state
            return
        end
        r0,c0 = divrem(event-1,nc)
        r,c = r0+1,c0+1
        i = isodd(r) ? 2c : 2c-1
        t = terminal_half_strength && r==total_layers ? reference.tau/2 : reference.tau
        f = fibonacci_outcome_filter(t)
        high,low = f[1,1],f[1,2]
        pstate = reference.vacuum_projectors[i]*state
        visit(low*state+(high-low)*pstate,event+1,bits) # physical vacuum, raw true
        bit = r+(c-1)*num_layers-1
        visit(high*state+(low-high)*pstate,event+1,bits | (1<<bit))
    end
    visit(reference.initial_state,1,0)
    probabilities = vec(sum(abs2,W;dims=1))
    boundary_density = W*W'
    norm2 = sum(probabilities)
    purity = real(sum(abs2,boundary_density))/norm2^2
    amplitudes,postselection = nothing,nothing
    if final_state !== nothing
        length(final_state)==n || throw(DimensionMismatch("final boundary vector"))
        f = ComplexF64.(final_state)
        norm(f)>0 && isfinite(norm(f)) || throw(ArgumentError("invalid final boundary"))
        normalize!(f)
        projected = vec(f'*W)
        postselection = sum(abs2,projected)
        postselection>0 || throw(ArgumentError("zero-probability final boundary"))
        amplitudes = Dict(Tuple(_sn_bit(k-1,e) for e in 1:events)=>projected[k]/sqrt(postselection)
                          for k in eachindex(projected))
    end
    return (joint_amplitudes=W,probabilities=probabilities,norm2=norm2,
            record_purity=purity,amplitudes=amplitudes,
            postselection_probability=postselection,record_shape=(num_layers,nc))
end
