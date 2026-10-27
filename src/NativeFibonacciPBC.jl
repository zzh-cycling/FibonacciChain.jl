# Native PBC via canonical fusion trees and the flux-leg construction supplied
# by the TensorKit/MPSKit maintainer. No closed-loop tensor basis is assumed.

"""Canonical physical lattice and fixed flux (`:trivial` or `:tau`) of a Fibonacci ring."""
struct NativePBCSites{S}
    N::Int
    flux::Symbol
    lattice::Vector{S}
end

_native_tau() = TensorKit.Vect[TensorKit.FibonacciAnyon](:τ => 1)
_native_sector(sites::NativePBCSites) = TensorKit.FibonacciAnyon(sites.flux == :trivial ? :I : :τ)
_native_right(sites::NativePBCSites) = sites.flux == :trivial ? TensorKit.unitspace(first(sites.lattice)) : _native_tau()

"""
    native_pbc_sites(N; flux=:trivial)

Use N physical τ sites in trivial flux. In τ flux append a τ auxiliary site
and fix the right MPS boundary to τ. The auxiliary leg is not a physical site.
The stable encoded backend is unaffected.
"""
function native_pbc_sites(N::Int; flux::Symbol=:trivial)
    N >= 3 || throw(ArgumentError("native Fibonacci PBC requires N >= 3"))
    flux in (:trivial,:tau) || throw(ArgumentError("flux must be :trivial or :tau"))
    return NativePBCSites(N,flux,fill(_native_tau(),N+Int(flux == :tau)))
end
Base.length(sites::NativePBCSites) = length(sites.lattice)
Base.getindex(sites::NativePBCSites,i::Int) = sites.lattice[i]

"""
    NativeFibonacciPBC

A native MPSKit state in a fixed ring flux sector, represented on canonical
fusion trees. `mps` contains N sites for trivial flux or N+1 sites for τ flux.
Its amplitudes and Schmidt entropy are NOT the encoded fusion-label amplitudes
and encoded interval entropy of the stable backend.
"""
struct NativeFibonacciPBC{M,S}
    mps::M
    sites::S
end
Base.copy(ψ::NativeFibonacciPBC) = NativeFibonacciPBC(copy(ψ.mps),ψ.sites)
LinearAlgebra.norm(ψ::NativeFibonacciPBC) = norm(ψ.mps)
LinearAlgebra.normalize!(ψ::NativeFibonacciPBC) = (normalize!(ψ.mps); ψ)
LinearAlgebra.normalize(ψ::NativeFibonacciPBC) = normalize!(copy(ψ))
Base.:*(a::Number,ψ::NativeFibonacciPBC) = NativeFibonacciPBC(_kit_scale(ψ.mps,a),ψ.sites)
Base.:*(ψ::NativeFibonacciPBC,a::Number) = a*ψ
function LinearAlgebra.dot(ψ::NativeFibonacciPBC,χ::NativeFibonacciPBC)
    ψ.sites.N == χ.sites.N || throw(DimensionMismatch("ring lengths differ"))
    ψ.sites.flux == χ.sites.flux || return 0.0+0.0im
    return dot(ψ.mps,χ.mps)
end

function _native_validate(model,sites::NativePBCSites)
    model.pbc || throw(ArgumentError("native flux backend requires PBC"))
    model.N == sites.N || throw(DimensionMismatch("model and ring lengths differ"))
    model.measure_operator in (:Antiferro,:Ferro) || throw(ArgumentError("unsupported Fibonacci interaction"))
end

# The two braids form a monodromy around the auxiliary flux leg, exactly as
# in the supplied three-site MPSKit construction. Omitting it changes the model.
function _native_wrap(gate,sites::NativePBCSites)
    sites.flux == :trivial && return gate
    V = _native_tau()
    d = (4,5,6)
    M = TensorKit.braid(TensorKit.id(ComplexF64,V ⊗ V ⊗ V),((1,3,2),d),(1,-100,3,d...))
    M = TensorKit.braid(M,((1,3,2),d),(1,2,100,d...))
    return M' * (gate ⊗ TensorKit.id(ComplexF64,V)) * M
end

function _native_term(gate,sites::NativePBCSites,i::Int)
    1 <= i <= sites.N || throw(BoundsError(1:sites.N,i))
    i < sites.N && return (i,i+1) => gate
    indices = sites.flux == :trivial ? (1,sites.N) : (1,sites.N,sites.N+1)
    return indices => _native_wrap(gate,sites)
end

function anyon_ham(model::AnyonModel{FibonacciAnyon},sites::NativePBCSites)
    _native_validate(model,sites)
    vacuum,_ = fibonacci_category_projectors()
    h = (model.measure_operator == :Antiferro ? -1 : 1)*vacuum
    return MPSKit.FiniteMPOHamiltonian(sites.lattice,[_native_term(h,sites,i) for i in 1:sites.N]...)
end

function _native_random_mps(sites,bond_dim,seed)
    bond_dim >= 1 || throw(ArgumentError("bond_dim must be positive"))
    V = _native_tau()
    D = TensorKit.Vect[TensorKit.FibonacciAnyon](:I=>bond_dim,:τ=>bond_dim)
    L = length(sites)
    spaces = fill(D,L+1)
    spaces[1],spaces[end] = TensorKit.unitspace(V),_native_right(sites)
    for i in 2:L
        spaces[i] = TensorKit.infimum(TensorKit.fuse(spaces[i-1],sites[i-1]),D)
    end
    for i in L:-1:2
        spaces[i] = TensorKit.infimum(spaces[i],TensorKit.fuse(spaces[i+1],TensorKit.dual(sites[i])))
    end
    rng = MersenneTwister(seed)
    tensors = map(1:L) do i
        t = zeros(ComplexF64,spaces[i] ⊗ sites[i] ← spaces[i+1])
        for (_,block) in TensorKit.blocks(t)
            randn!(rng,block)
        end
        t
    end
    return MPSKit.FiniteMPS(tensors;normalize=true)
end

"""
    native_pbc_state(N; flux=:trivial, bond_dim=4, seed=1234)
    native_pbc_state(amplitudes, N; flux=:trivial, reference_limit=16)

Construct a seeded random native MPS, or convert a small vector in TensorKit's
canonical fusion-tree block order. `bond_dim` is the multiplicity of each
virtual sector. Vector coefficients are divided by sqrt(d_flux) to make their
Euclidean norm equal the categorical state norm. No encoded-basis identification
is implied; vector conversion is exponential and size guarded.
"""
function native_pbc_state(N::Int;flux::Symbol=:trivial,bond_dim::Int=4,seed::Int=1234)
    sites = native_pbc_sites(N;flux)
    return NativeFibonacciPBC(_native_random_mps(sites,bond_dim,seed),sites)
end
function native_pbc_state(v::AbstractVector,N::Int;flux::Symbol=:trivial,reference_limit::Int=16)
    N <= reference_limit || throw(ArgumentError("canonical amplitude conversion exceeds reference_limit"))
    sites = native_pbc_sites(N;flux)
    all(isfinite,v) || throw(ArgumentError("amplitudes must be finite"))
    norm(v) > 0 || throw(ArgumentError("cannot construct a zero MPS from amplitudes"))
    t = zeros(ComplexF64,prod(sites.lattice) ← _native_right(sites))
    block = TensorKit.block(t,_native_sector(sites))
    length(v) == size(block,1) || throw(DimensionMismatch("expected $(size(block,1)) canonical amplitudes"))
    block[:,1] .= v / sqrt(TensorKit.dim(_native_sector(sites)))
    scale = norm(t)
    with_left = TensorKit.insertleftunit(t/scale,1)
    p = MPSKit.FiniteMPS(MPSKit.decompose_localmps(with_left,TensorKit.notrunc()))
    return NativeFibonacciPBC(_kit_scale(p,scale),sites)
end

function _initial_native_pbc(N::Int;flux::Symbol=:trivial)
    ψ = native_pbc_state(N;flux)
    return ψ,ψ.sites
end

"""Return amplitudes in the canonical fusion-tree block order (small-system reference only)."""
function native_pbc_amplitudes(ψ::NativeFibonacciPBC;reference_limit::Int=16)
    ψ.sites.N <= reference_limit || throw(ArgumentError("canonical amplitude extraction exceeds reference_limit"))
    t = ψ.mps.AC[1]
    for i in 2:length(ψ.mps)
        tail = TensorKit.repartition(ψ.mps.AR[i],1,TensorKit.numind(ψ.mps.AR[i])-1)
        t = t*tail
        t = TensorKit.repartition(t,TensorKit.numind(t)-1,1)
    end
    t = TensorKit.removeunit(t,1)
    c = _native_sector(ψ.sites)
    return vec(copy(TensorKit.block(t,c))) * sqrt(TensorKit.dim(c))
end

function _native_groundstate(model;flux=:trivial,sweep_times=20,maxdim=100,cutoff=1e-10,
    outputlevel=1,seed=1234,initial_linkdim=min(8,maxdim))
    sweep_times >= 1 || throw(ArgumentError("sweep_times must be positive"))
    1 <= initial_linkdim <= maxdim || throw(ArgumentError("invalid initial_linkdim"))
    sites = native_pbc_sites(model.N;flux)
    H = anyon_ham(model,sites)
    ψ = _native_random_mps(sites,initial_linkdim,seed)
    ψ,_,_ = MPSKit.find_groundstate(ψ,H,MPSKit.DMRG2(;maxiter=sweep_times,
        verbosity=outputlevel,tol=max(cutoff,1e-12),trscheme=_kit_truncation(cutoff,1,maxdim)))
    normalize!(ψ)
    return NativeFibonacciPBC(ψ,sites),real(MPSKit.expectation_value(ψ,H))
end

"""Normalized energy in the fixed flux sector, computed directly with the native MPO."""
function native_pbc_energy(model::AnyonModel{FibonacciAnyon},ψ::NativeFibonacciPBC)
    return real(MPSKit.expectation_value(ψ.mps,anyon_ham(model,ψ.sites)))
end

"""
    native_pbc_dense(model; flux=:trivial, construction=:mpo, reference_limit=12)

Small-system Hamiltonian block. `:mpo` contracts the production local-term MPO;
`:braid` independently implements the supplied full-ring ED braiding. In τ
flux these conventions give complex-conjugate matrices, with the same spectrum.
Use `:mpo` with `native_pbc_amplitudes`; do not mix the two amplitude conventions.
"""
function native_pbc_dense(model::AnyonModel{FibonacciAnyon};flux::Symbol=:trivial,
    construction::Symbol=:mpo,reference_limit::Int=12)
    model.N <= reference_limit || throw(ArgumentError("dense construction exceeds reference_limit"))
    sites = native_pbc_sites(model.N;flux)
    _native_validate(model,sites)
    if construction == :mpo
        H = convert(TensorKit.TensorMap,anyon_ham(model,sites))
    elseif construction == :braid
        N = model.N
        V = _native_tau()
        one = TensorKit.id(ComplexF64,V)
        vacuum,_ = fibonacci_category_projectors()
        h = (model.measure_operator == :Antiferro ? -1 : 1)*vacuum
        term(i) = reduce(⊗,[fill(one,i-1);h;fill(one,N-1-i)])
        if flux == :trivial
            dom = ntuple(i->N+i,N)
            T = TensorKit.braid(TensorKit.id(ComplexF64,V^N),((N,1:N-1...),dom),(1:N-1...,100,dom...))
            H = sum(term,1:N-1) + T'*term(1)*T
        else
            dom = ntuple(i->N+1+i,N+1)
            P = V^N ⊗ V
            T = TensorKit.braid(TensorKit.id(ComplexF64,P),((1:N-1...,N+1,N),dom),(1:N-1...,100,N,dom...))
            T = TensorKit.braid(T,((N+1,1:N-1...,N),dom),(1:N...,-100,dom...))
            H = sum(i->term(i) ⊗ one,1:N-1) + T'*(term(1) ⊗ one)*T
        end
    else
        throw(ArgumentError("construction must be :mpo or :braid"))
    end
    return Matrix(TensorKit.block(H,_native_sector(sites)))
end

"""
    native_schmidt_entropy(state, bond)

MPSKit's categorical Schmidt entropy at a canonical MPS bond. In τ flux the
auxiliary leg is part of that representation. This is not the encoded PBC
interval entropy, nor a claim of equivalence to a physical annular interval RDM.
"""
function native_schmidt_entropy(ψ::NativeFibonacciPBC,b::Int)
    1 <= b < length(ψ.mps) || throw(BoundsError(1:length(ψ.mps)-1,b))
    return MPSKit.entropy(normalize(ψ.mps),b)
end
function ee_mps(::NativeFibonacciPBC,::Int)
    throw(ArgumentError("native flux Schmidt entropy differs from encoded interval entropy; call native_schmidt_entropy explicitly"))
end
function mps_statevector(::NativeFibonacciPBC)
    throw(ArgumentError("no encoded-basis isometry is assumed; use native_pbc_amplitudes for canonical fusion-tree coefficients"))
end

function measurement_operator_mps(model::AnyonModel{FibonacciAnyon},sites::NativePBCSites,
    i::Int,τ::Real,outcome::Bool)
    _native_validate(model,sites)
    p,q = fibonacci_category_projectors()
    a,b = _kit_measurement_coefficients(τ,outcome)
    model.measure_operator == :Ferro && (b = -b)
    gate = (a-b)*p+(a+b)*q
    return MPSKit.FiniteMPOHamiltonian(sites.lattice,_native_term(gate,sites,i))
end

function measuremap(model::AnyonModel{FibonacciAnyon},ψ::NativeFibonacciPBC,
    sites::NativePBCSites,i::Int,τ::Real,outcome::Bool;cutoff::Float64=1e-10,
    mindim::Int=1,maxdim::Int=100)
    ψ.sites.N == sites.N && ψ.sites.flux == sites.flux || throw(DimensionMismatch("state and lattice differ"))
    out = _kit_compress(_kit_apply(measurement_operator_mps(model,sites,i,τ,outcome),ψ.mps);
        cutoff,mindim,maxdim)
    probability = norm(out)^2
    probability > 0 || throw(ArgumentError("zero-probability outcome"))
    normalize!(out)
    return NativeFibonacciPBC(out,sites),probability
end

function bulk_evolution(model::AnyonModel{FibonacciAnyon},sites::NativePBCSites,
    ψ::NativeFibonacciPBC,mc::MeasureConfig,samples::Union{Nothing,BitMatrix}=nothing;
    normalized::Bool=true,compute_entropy::Bool=false,entropy::Symbol=:none)
    _native_validate(model,sites)
    ψ.sites.N == sites.N && ψ.sites.flux == sites.flux || throw(DimensionMismatch("state and lattice differ"))
    iseven(model.N) || throw(ArgumentError("two-layer protocol requires even N"))
    compute_entropy && entropy != :native_schmidt && throw(ArgumentError("specify entropy=:native_schmidt explicitly"))
    mc.mode in (:Born,:sample) || throw(ArgumentError("mode must be :Born or :sample"))
    mc.mode == :Born && !normalized && throw(ArgumentError("Born mode requires normalization"))
    mc.truncate_every_events >= 1 || throw(ArgumentError("truncate_every_events must be positive"))
    normalized && !isapprox(norm(ψ),1;atol=1e-10,rtol=1e-10) && throw(ArgumentError("initial state must be normalized"))
    _kit_truncation(mc.cutoff,mc.mindim,mc.maxdim)
    periods = mc.t₂-mc.t₁+1
    periods >= 0 || throw(ArgumentError("negative number of periods"))
    dims = (2periods,_samples_per_layer(model))
    born = mc.mode == :Born
    !born && (samples === nothing || size(samples) != dims) && throw(DimensionMismatch("sample dimensions"))
    record = born ? falses(dims) : copy(samples)
    energies,entropies = zeros(Float32,dims[1]),zeros(Float32,periods)
    φ = (1+sqrt(5.0))/2
    ys = mc.track_y_expectation ? fill(Float32(sites.flux == :trivial ? φ : -1/φ),periods) : Float32[]
    current,cache = copy(ψ.mps),Dict()
    compress(p) = _kit_compress(p;cutoff=mc.cutoff,mindim=mc.mindim,maxdim=mc.maxdim)
    for period in 1:periods
        for layer in 1:2
            row = 2(period-1)+layer
            strength = row == dims[1] && mc.enable_τ_eff ? mc.τ/2 : mc.τ
            positions,localmodel,τ = _obtain_measurement_config(model,row,strength)
            cols = _get_sample_column_indices(model,row)
            energy = 0.0
            for (k,i) in enumerate(positions)
                ops = get!(cache,(i,τ)) do
                    (measurement_operator_mps(localmodel,sites,i,τ,false),measurement_operator_mps(localmodel,sites,i,τ,true))
                end
                applybranch(bit) = begin
                    out = _kit_apply(ops[Int(bit)+1],current)
                    mc.truncate_every_events == 1 ? compress(out) : out
                end
                if born
                    branch0 = applybranch(false)
                    p0 = clamp(norm(branch0)^2,0.0,1.0)
                    bit = rand(mc.rng) >= p0
                    record[row,cols[k]] = bit
                    current = bit ? applybranch(true) : branch0
                    energy -= bit ? log1p(-p0) : log(p0)
                    normalize!(current)
                else
                    current = applybranch(record[row,cols[k]])
                    if normalized
                        energy -= log(norm(current)^2)
                        normalize!(current)
                    end
                end
                if mc.truncate_every_events > 1 && (k % mc.truncate_every_events == 0 || k == length(positions))
                    current = compress(current)
                    normalized && normalize!(current)
                end
            end
            energies[row] = energy
        end
        compute_entropy && (entropies[period] = native_schmidt_entropy(NativeFibonacciPBC(current,sites),model.N÷2))
    end
    return Measurement_outcome_mps_bulk(NativeFibonacciPBC(current,sites),record,energies,entropies,ys)
end
