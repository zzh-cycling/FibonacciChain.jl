# MPSKit/TensorKit backend migration report

Date: 2026-10-26. Branch: `catsymMPS`.

## Scope and physical conventions

The default MPS constructors and ground-state solver now select MPSKit/TensorKit.
The ITensor implementation remains available through `backend=:itensor`.
This is an equivalent encoded-backend migration: it retains BitBasis ordering,
the fusion-label/qubit representation, closed-path PBC Hilbert space, and ordinary
encoded von Neumann entropy in nats. It does not reinterpret encoded entropy as
categorical quantum-trace entropy.

## Implementation

- Added `src/MPSKitBackend.jl`, MPSKit/TensorKit dependencies, public conversion
  helpers, and generic state fields in measurement outcomes.
- Implemented Fibonacci, Ising, and OBF Hamiltonians; DMRG; measurement MPOs;
  Born and fixed-record trajectories; unnormalized replay; hybrid evolution;
  encoded entropy; constraint and topological-charge MPOs; Lyapunov frames;
  measurement enumeration; controlled-copy reference insertion; and Fibonacci
  charge sharpening with an ancilla.
- Kept the existing Born sampling convention and matched cutoff units through
  the square root of the discarded squared-weight tolerance. Finite-bond
  compression paths can differ between libraries.
- Protected small unnormalized amplitudes against absolute truncation thresholds
  by scaling around MPO application and compression. Added direct-sum MPS
  addition to preserve relative phases in sector projections and frame updates.
  These operations use local tensors, not dense full-state conversion.
- Added native TensorKit Fibonacci projectors and tests of fusion-channel quantum
  traces, F-matrix recoupling, unitarity, and Kraus completeness. These category
  objects are separate from the encoded dynamics.
- Preserved the existing ITensor regression tests with explicit backend selection.
  Historical research scripts also explicitly select ITensor to retain their
  state types and diagnostics.
- Updated documentation and added `exm/benchmark_mps_backends.jl`.

## Validation

All 4,841 assertions passed across staged runs, including 203 new backend/category
checks and 1,251 existing MPS assertions. Coverage includes dense operator and ED
comparisons, DMRG energies and entropy, same-record trajectories, tiny-amplitude
replay, hybrid schedules, both topological sectors, ancillas, and finite-bond
convergence. Sector projections and short evolutions agree with ED to about
1e-15 in the diagnostic small system.

`Pkg.test()` exposed two pre-existing HybridProtocol test-environment/fixture
issues: missing `Statistics` in test dependencies and expectations for the old
schedule storage/default. Added the dependency; updated the fixture to read the
record from the trajectory file and explicitly request `save_schedule=false`
when testing omitted records. The remaining test files were run with bounds
checking in a temporary environment containing the package and declared test
dependencies. After the fixture fixes, HybridProtocol passed all 150 assertions.
The final result is staged suite coverage, not a single uninterrupted successful
`Pkg.test()` invocation after the last fixture edit.

Documenter build and doctests passed locally without deployment.
`git diff --check` passed.

Reproduction commands:

```sh
julia --project -e 'using Pkg; Pkg.test()'
julia --project test/test_MPSKitBackend.jl
julia --project exm/benchmark_mps_backends.jl 8 32 4
```

## Performance evidence

One warmed local run used Julia 1.12.2, MPSKit 0.13.13, TensorKit 0.17.0,
one BLAS thread, L=8, maxdim=32, four periods, and an identical fixed record.

| Backend | Time (s) | Cumulative allocated bytes | Final bond dimension |
|:--|--:|--:|--:|
| ITensor | 0.024273125 | 37,433,024 | 8 |
| MPSKit | 0.161948958 | 224,007,760 | 8 |

Final fidelity was 0.9999999999999973; stored Float32 free energies agreed
exactly and both final entropies were 0.40336457. MPSKit was about 6.7 times
slower in this small-system run. Allocated bytes are not peak memory. This is a
correctness baseline, not evidence of a speedup or an asymptotic scaling result.

## Remaining limitations

- Native anyonic PBC dynamics and categorical block sparsity are not implemented;
  potential non-Abelian/category symmetry speedups remain unmeasured.
- Legacy reset-reference and spacetime reference-correlation workflows still
  require `backend=:itensor`.
- Rank-deficient Lyapunov frames raise an error rather than inserting an arbitrary
  replacement direction.
- The existing OBF OBC dense Hamiltonian has a boundary-index issue; that operator
  comparison uses the full independent ITensor MPO matrix instead of ED.
