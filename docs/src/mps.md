# Matrix Product State Methods

The default MPS backend is MPSKit/TensorKit. The ITensor implementation remains
available with `backend=:itensor`, including its original regression tests.

## Selecting a backend

```julia
using FibonacciChain, Random

model = AnyonModel(FibonacciAnyon(), 8; pbc=true)
ψ, sites = initial_mps(8) # MPSKit.FiniteMPS and TensorKit.ComplexSpace sites
config = MeasureConfig(τ=0.8, t₂=2, mode=:Born, rng=MersenneTwister(36))
result = bulk_evolution(model, sites, ψ, config)

ψ_old, sites_old = initial_mps(8; backend=:itensor)
replay = MeasureConfig(τ=0.8, t₂=2, mode=:sample)
old_result = bulk_evolution(model, sites_old, ψ_old, replay, result.samples)

ψ_gs, E0 = anyon_mps_gst(model; backend=:mpskit, maxdim=64)
H = anyon_ham(model, mpskit_sites(model.N))
```

This migration deliberately retains the **fusion-label/qubit encoding**, the
closed-path PBC Hilbert space, and ordinary encoded entropy in nats. It does not
replace the physical sites by categorical τ objects. The Fibonacci constraint
penalty/projection is consequently still required. `fibonacci_category_projectors()`
provides separate native TensorKit projectors on τ ⊗ τ for category foundations;
these are not interchangeable with the encoded PBC operators.

MPSKit implements Fibonacci, Ising and OBF Hamiltonians; measurements and
Born/fixed-record bulk and boundary evolution; Fibonacci hybrid evolution;
encoded entropy; constraint and topological-charge MPOs; Lyapunov frames;
measurement enumeration; controlled-copy reference insertion; and Fibonacci
charge sharpening with an ancilla. The older reset-reference and spacetime
reference-correlation workflows still require `backend=:itensor`. Historical
research scripts explicitly select ITensor to preserve their diagnostics and
stored state types. A separate experimental native PBC backend is described below.

All encoded MPSKit evolution, compression, DMRG and overlaps use MPSKit/TensorKit. No
trajectory is converted through ITensor or an exponentially large dense vector.
`mpskit_state(old_mps)` converts an existing qubit ITensor MPS tensor by tensor.
For small-system validation only, `mpskit_state(vector, N)` and
`mps_statevector(state)` use full 2^N amplitudes in BitBasis order (site 1 is the
most significant bit).

The truncation cutoff is interpreted as relative discarded squared weight;
the TensorKit discarded-norm tolerance is therefore `sqrt(cutoff)`. A hard
`maxdim` and a `mindim` floor are also applied. Compression order differs
between libraries, so finite-bond states and seeded trajectories need not be
bitwise identical. Compare fixed records and converged observables. The legacy
Born policy (compressed outcome-0 norm, complemented for outcome 1) is retained;
fixed-record replay instead records the selected branch norm. Unnormalized
replay retains the Kraus prefactors and does not renormalize after compression.
Rank-deficient Lyapunov frames raise an explicit error instead of inserting an
arbitrary replacement direction.

Run the independent cross-backend tests and a compilation-warmed benchmark:

```sh
julia --project=. test/test_MPSKitBackend.jl
julia --project=. exm/benchmark_mps_backends.jl 8 32 4
```

The benchmark uses one fixed record, reports wall time and cumulative allocated
bytes (not peak memory), and computes the final overlap without dense conversion.
The generic MPO-per-event implementation is a correctness baseline; it does not
yet exploit categorical block sparsity or promise a speedup over local ITensor
gates. The OBF OBC dense Hamiltonian currently has an independent boundary-index
bug; that case is checked against the full ITensor MPO matrix instead of ED.

A local warmed run on 2026-10-26 (Julia 1.12.2, MPSKit 0.13.13,
TensorKit 0.17.0, one BLAS thread; L=8, maxdim=32, four periods) gave:

| Backend | Time (s) | Allocated bytes | Final bond dimension |
|:--|--:|--:|--:|
| ITensor | 0.0243 | 37,433,024 | 8 |
| MPSKit | 0.162 | 224,007,760 | 8 |

The final fidelity was 0.9999999999999973 and the recorded Float32 free
energies agreed exactly. This small-system result establishes equivalence,
not a performance gain: the generic MPSKit MPO application was about 6.7 times
slower in this run. Larger-system scaling and native category block sparsity
require separate benchmarks.

## Native PBC with canonical fusion trees and flux legs

The opt-in `backend=:native_pbc` follows the TensorKit/MPSKit maintainer's
canonical-fusion-tree construction. It does not treat closed-loop diagrams as
an independent TensorKit basis, and it does not change the default encoded backend.

```julia
model = AnyonModel(FibonacciAnyon(), 4; pbc=true)
ψ1, E1 = anyon_mps_gst(model; backend=:native_pbc, flux=:trivial,
    maxdim=32, sweep_times=12, outputlevel=0)
ψτ, Eτ = anyon_mps_gst(model; backend=:native_pbc, flux=:tau,
    maxdim=32, sweep_times=12, outputlevel=0)
# E1 ≈ -3.2360679775; Eτ ≈ -3.0
Hτ = anyon_ham(model, native_pbc_sites(4; flux=:tau))
```

In trivial flux there are N physical τ sites and vacuum MPS boundaries.
The boundary interaction is the long-range local term on `(1, N)`.
In τ flux the canonical lattice contains N physical τ sites plus a final
auxiliary τ leg W; the left boundary is vacuum and the right boundary is τ.
The auxiliary leg is not counted in `model.N`. The seam term acts on
`(1, N, N+1)` and is `M' * (h ⊗ id(W)) * M`, where the two explicit braids
in M wind the interaction around W. Simply omitting this winding gives a
different spectrum. All production Hamiltonian, DMRG, and measurement operations
use native category tensors and MPOs, without full-state expansion.

`native_pbc_dense(model; flux, construction=:braid)` reproduces the supplied
full-ring ED braiding prescription. `construction=:mpo` contracts the production
local-term MPO. At N=4 and N=6 the two prescriptions are identical matrices in
trivial flux and complex-conjugate matrices in τ flux, with identical spectra.
Use `:mpo` when comparing against `native_pbc_amplitudes`; do not mix the amplitude
conventions. Their combined spectra and individual sector dimensions/spectra
are checked against the encoded ED Hamiltonian resolved by Y.

`native_pbc_state(N; flux, bond_dim=4, seed=1234)` constructs a normalized random
MPS using a local RNG. `bond_dim` specifies the multiplicity of each virtual
sector, not an encoded bond dimension. `initial_mps(N; backend=:native_pbc, flux)`
uses this seeded random initializer; it is not the encoded all-zero fusion path.
`native_pbc_state(vector, N; flux)` and `native_pbc_amplitudes` are small-system
reference conversions in TensorKit's canonical block order. Their sqrt(d_flux)
normalization makes the vector's Euclidean norm equal the categorical MPS norm.
The conversions are exponentially sized and default to N ≤ 16; the production
algorithms have no such reference-size cutoff.

Local measurements and Born/fixed-record bulk trajectories include the same
braided seam term. The two-layer protocol requires even N. The dense Fibonacci
measurement routine fixes AFM outcome labels even for `:Ferro`; native measurements
follow the existing MPS convention, which flips those labels for `:Ferro`.

```julia
ψ, sites = initial_mps(6; backend=:native_pbc, flux=:tau)
config = MeasureConfig(τ=0.8, t₂=2, mode=:Born, rng=MersenneTwister(36),
    track_y_expectation=true)
model = AnyonModel(FibonacciAnyon(), 6)
result = bulk_evolution(model, sites, ψ, config)
```

For a fixed flux, Y tracking returns its eigenvalue (φ or -1/φ); it is not a
numerical contraction of an independently constructed native Y MPO.
`enforce_fibonacci_constraint` needs no extra projector because fusion admissibility
is built into the tensor spaces. Deferred compression and unnormalized fixed-record
replay are supported.

**Entropy and basis boundary:** native canonical amplitudes are not the encoded
fusion-label amplitudes. An explicit isometry to the encoded basis is not yet
provided. Consequently `mps_statevector` and `ee_mps` reject native states rather
than silently changing their physical meaning. `native_schmidt_entropy(ψ, bond)`
exposes MPSKit's categorical Schmidt entropy of the canonical state, including
the auxiliary-leg representation; this does not establish equivalence to the
encoded PBC interval entropy or a physical annular interval RDM. Bulk evolution
defaults to `compute_entropy=false`; request `compute_entropy=true,
entropy=:native_schmidt` explicitly to record this observable.

Native hybrid, Lyapunov-frame, and reference-qubit workflows, coherent superpositions
of different flux sectors, and an encoded/native basis isometry remain future work.
Use the encoded backend for those workflows. Native category block sparsity is
active, but no speedup over the encoded backend is claimed without a benchmark.

```sh
julia --project test/test_NativeFibonacciPBC.jl
```

## Ground State Calculations


MPS methods allow simulation of much larger systems (N ~ 50-100) compared to exact diagonalization:

```julia
using FibonacciChain

# Find ground state for N=50 Fibonacci anyon chain
N = 50
model = AnyonModel(FibonacciAnyon(), N; pbc=true)
ψ, E0 = anyon_mps_gst(model; maxdim=100, cutoff=1e-10)
```

## Hamiltonian Construction


The MPS Hamiltonian is constructed as a Matrix Product Operator (MPO) with three-body interactions:

```julia
# Create sites and Hamiltonian MPO
sites = mpskit_sites(N)
H = anyon_ham(model, sites)
```

## MPS Measurements

### Measurement Operators


### MPS-based Protocols



## Entanglement Entropy



MPS naturally provides access to entanglement structure:

```julia
# Calculate entanglement entropy between first L and last N-L sites
L = N ÷ 2
S = ee_mps(ψ, L)

# Get full entanglement profile
ee_profile = anyon_eelis(model, ψ)
```

## State Generation and Evolution



## Advantages of MPS Methods

### Computational Efficiency
- **Memory**: Scales as O(χ²N) instead of O(2ᴺ)
- **Time**: DMRG scales as O(χ³N) for ground states
- **Bond dimension**: χ ~ 50-200 typically sufficient

### Physical Insight
- Direct access to entanglement structure
- Natural handling of gapped phases
- Efficient for 1D systems with area law

## Usage Example: Large System Simulation

```julia
using FibonacciChain

N = 80  # Much larger than exact diagonalization allows
maxdim = 200
cutoff = 1e-12

# Find ground state
model = AnyonModel(FibonacciAnyon(), N; pbc=true)
ψ_gs, E0 = anyon_mps_gst(model,
                         sweep_times=30,
                         maxdim=maxdim,
                         cutoff=cutoff)

# Calculate observables
ee_profile = anyon_eelis(model, ψ_gs)
central_ee = ee_profile[N÷2]

println("Ground state energy: $E0")
println("Central entanglement entropy: $central_ee")

# Apply measurement and evolve
τ = 0.5
measurement_site = N÷2
sites = mpskit_sites(N)
ψ_measured, prob = measuremap(model, ψ_gs, sites, measurement_site, τ, false)

println("Measurement probability: $prob")
```

## Performance Tips

1. **Bond dimension**: Start with small χ and increase gradually
2. **Convergence**: Monitor energy convergence during DMRG sweeps  
3. **Memory**: Use appropriate number precision (Float64 vs ComplexF64)
4. **Parallelization**: ITensors supports threading for large calculations

### Born evolution performance (issue #36)

Born MPS evolution caches the two measurement operators for each layer phase
within a trajectory. The final half-strength layer has a separate cache entry.
Probabilities use `norm(ψ)^2`, which contracts only the orthogonality center when
the resulting MPS is canonical, and falls back to a full contraction otherwise.
This avoids rebuilding fixed gates and repeatedly contracting the entire chain.
The public API, RNG draw order, and truncation settings are unchanged.

Keep `truncate_every_events=1` as the starting point. Deferring truncation can
increase intermediate bond dimensions substantially, making SVDs more expensive
and increasing memory use; a larger interval is not necessarily faster.

Run the reproducible benchmark from the repository root:

```sh
julia --project=. exm/benchmark_mps_issue36.jl 3 1
```

The arguments are the number of seeded repetitions and `truncate_every_events`.
The script warms up compilation, fixes BLAS to one thread, and reports elapsed
time, cumulative allocated bytes, and final maximum bond dimension. Allocated
bytes measure allocation traffic, not peak resident memory. Use the same script
with `--project=/path/to/baseline` to compare revisions, running them sequentially.

Local results on Apple Silicon with Julia 1.12.5, one Julia thread and one BLAS
thread (2026-09-13), comparing baseline `d2efae7` with commit `da00223`:

| L | maxdim | Periods | Before (s) | After (s) | Before allocated (GB) | After allocated (GB) |
|---|--------|---------|------------|-----------|-----------------------|----------------------|
| 8 | 64 | 80 | 0.571 | 0.346 | 1.706 | 0.918 |
| 16 | 32 | 32 | 0.932 | 0.651 | 3.249 | 2.118 |
| 16 | 64 | 32 | 0.996 | 0.671 | 3.418 | 2.273 |
| 32 | 64 | 64 | 26.256 | 21.307 | 60.690 | 45.861 |

Each entry is the median of seeds 1–3, with periodic Fibonacci chains,
`cutoff=1e-12`, `truncate_every_events=1`, and the default final half-strength
layer. The strength is `atanh(1/sqrt(2))` for L=8 and `atanh(0.95)` otherwise.
GB denotes decimal gigabytes of cumulative allocations. The L=32 workload takes
about 19% less time and allocates 24% fewer bytes. These are local measurements,
not an asymptotic scaling fit or a peak-memory measurement; gate application and
SVD compression still dominate larger runs.

An additional before/after comparison of all 12 seeded trajectories found
identical samples and stored Float32 free energies and entropies; final-state
fidelities differed from one by less than `1e-8`.

Regression tests compare against the original Born algorithm for Fibonacci,
Ising, and OBF models, including both outcomes, deferred truncation, and the final
half-strength layer. They check samples, free energies, entropies, final states,
RNG consumption, and replay of recorded outcomes.

### Compact periodic boundary MPOs

A subsequent profile of the L=32 benchmark put about 78% of sampled evolution
time in MPO application. The two periodic boundary measurement MPOs built by
`OpSum` carried four bond channels through the entire chain. The distant
neighbor only contributes `Proj0` or `Proj1`, so two channels suffice in the
bulk. Constructing those channels explicitly gives bond dimensions
`[3, 2, …, 2]` for a measurement centered at site 1, and the reverse for site N.
The three-channel bond carries the measured site's `I`, `Z`, and `X` terms.
This is an exact representation of the same operator, with no operator
truncation or change to the MPS compression settings.

Sequential before/after measurements against `da00223`, using the same machine,
Julia/BLAS thread counts, workloads, seeds, and warmup as above:

| L | maxdim | Before (s) | Compact MPO (s) | Before allocated (GB) | Compact MPO allocated (GB) |
|---|--------|------------|-----------------|-----------------------|---------------------------|
| 8 | 64 | 0.299 | 0.314 | 0.918 | 0.867 |
| 16 | 32 | 0.604 | 0.499 | 2.118 | 1.442 |
| 16 | 64 | 0.664 | 0.523 | 2.273 | 1.581 |
| 32 | 64 | 20.460 | 13.896 | 45.861 | 26.161 |

The L=32 case takes **32% less time** (1.47× speedup) and allocates **43% fewer
bytes** than the previously optimized version. The L=8 case shows no timing
improvement in this run; the benefit grows as boundary contractions become more
expensive. These figures measure cumulative allocation traffic, not peak memory.
The benchmark workload itself is unchanged. For the before/after validation,
the returned states were serialized outside the timed region for comparison.

All 12 trajectories retained identical sampled outcomes and stored free energies.
The largest stored entropy difference was `1.2e-7`, and all final-state
fidelities agreed with one within `1e-8`. The tests also compare the entire
boundary operator against the dense local gate, including states outside the
fusion constraint, both boundary sites, both measurement conventions and
outcomes, and strengths from zero through the projective limit. The trajectory
reference independently reconstructs the original four-channel `OpSum` MPO.

### Ising and OBF periodic measurements

Periodic spin measurements now use an exact bond-dimension-2 MPO for wrapping
`ZZ`, `XZZ`, and `ZZX` terms. Each gate has the form `c I + a P`, where `P` is a
Pauli string, so two product-operator channels suffice. This avoids moving sites
across the chain and back with SWAPs. Contiguous gates retain their local ITensor
application. Both Born sampling and recorded-outcome replay use the new path.

Run the spin benchmarks with the third argument:

```sh
julia --project=. exm/benchmark_mps_issue36.jl 3 1 Ising
julia --project=. exm/benchmark_mps_issue36.jl 3 1 OBF
```

The default remains the Fibonacci benchmark. The spin workloads use periodic
chains, `τ=0.8`, 12 periods, `cutoff=1e-12`, per-event truncation, the final
half-strength layer, and `λ=0.3` for OBF. Sequential comparisons against
`4295184` on Apple Silicon, Julia 1.12.5, one Julia thread and one BLAS thread
(2026-09-14), with compilation warmed up and medians of seeds 1–3:

| Model | L | maxdim | Before (s) | After (s) | Speedup | Before allocated (GB) | After allocated (GB) |
|-------|---|--------|------------|-----------|---------|-----------------------|----------------------|
| Ising | 6 | 16 | 0.0440 | 0.0274 | 1.61× | 0.119 | 0.075 |
| Ising | 24 | 32 | 0.475 | 0.288 | 1.65× | 1.249 | 0.786 |
| Ising | 24 | 64 | 0.724 | 0.443 | 1.63× | 1.896 | 1.228 |
| OBF | 6 | 16 | 0.339 | 0.173 | 1.96× | 0.861 | 0.443 |
| OBF | 24 | 32 | 5.915 | 2.919 | 2.03× | 13.230 | 6.677 |
| OBF | 24 | 64 | 15.619 | 8.073 | 1.93× | 29.992 | 15.962 |

Allocated GB measures cumulative traffic, not peak memory. Returned states were
saved outside the timed region to compare trajectories. All 18 runs retained
identical sampled outcomes, but finite-bond observables are not identical:
the old method also truncated during each intermediate SWAP, whereas the MPO
method compresses the resulting state in its original site order. Consequently,
matching a seed does not guarantee identical outcomes for arbitrary compressed
runs. The `maxdim`, `cutoff`, and truncation-interval settings still apply.

For L=24, the largest old/new final-state infidelity among the three seeds was
`5.1e-7` (Ising) and `4.9e-5` (OBF) at maxdim 32, falling to `7.7e-10` and
`9.7e-7` at maxdim 64. These are differences between approximations, not errors
against an exact solution. For an independent accuracy check, six-site,
three-period trajectories at seed 36 were replayed against dense evolution:

| Model | maxdim | Old SWAP infidelity | MPO infidelity |
|-------|--------|---------------------|----------------|
| Ising | 2 | 0.0171 | 0.0110 |
| Ising | 4 | 6.61e-5 | 5.73e-5 |
| OBF | 2 | 0.219 | 0.0310 |
| OBF | 4 | 0.00610 | 0.00241 |

Both methods agree with dense evolution to infidelity below `1e-9` at maxdim 8.
The improvement in these fixed trajectories is regression evidence, not a
universal ordering of truncation errors. Tests also verify complete gate
matrices, both outcomes, the OBF sign convention, zero and projective strengths,
seeded exact-limit trajectories, and the bond-dimension bound during replay.
