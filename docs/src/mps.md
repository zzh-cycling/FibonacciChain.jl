# Matrix Product State Methods

This section covers the efficient Matrix Product State (MPS) implementations for simulating larger anyon systems using the ITensors.jl library.

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
sites = siteinds("Qubit", N)
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
using FibonacciChain, ITensors

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
sites = siteinds("Qubit", N)
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
