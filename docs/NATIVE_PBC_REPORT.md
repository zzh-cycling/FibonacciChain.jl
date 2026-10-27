# Native Fibonacci PBC backend

This implementation follows the TensorKit/MPSKit maintainer's supplied
canonical-fusion-tree construction. It adds `backend=:native_pbc` and leaves
the default encoded MPSKit backend and the legacy ITensor tests in place.

## Representation and supported operations

- Trivial flux: N physical τ sites, vacuum boundaries, and a wrap interaction
  on `(1, N)`.
- τ flux: N physical τ sites and one auxiliary τ site, with right boundary τ.
  Two explicit braids wind the wrap interaction around the auxiliary leg;
  its MPO support is `(1, N, N+1)`.
- Hamiltonian construction, DMRG2, local measurements, Born trajectories,
  fixed-record replay, and deferred truncation use native category tensors.
  Production operations do not expand the complete state vector.
- Random initialization uses a local RNG. Small-system vector conversion
  uses canonical fusion-tree coefficients with a `sqrt(d_flux)` normalization
  to match the categorical norm.

The reference ED implementation also reproduces the supplied full-ring
braiding construction independently of the production local-term MPO.
At N=4 and N=6, the resulting matrices agree in trivial flux and are complex
conjugates in τ flux. Their spectra agree. Amplitude comparisons must use
`native_pbc_dense(...; construction=:mpo)`.

## Validation

`test/test_NativeFibonacciPBC.jl` covers:

- N=4 supplied spectra, N=4 and N=6 sector dimensions and spectra against
  the existing encoded Hamiltonian resolved by topological charge Y.
- Hermiticity, projector identities, periodic Temperley–Lieb relations,
  AFM/Ferro sign reversal, and a regression showing that omitting the flux
  winding changes the spectrum.
- Native DMRG ground energies and eigenvector residuals against ED in both
  flux sectors, canonical amplitude round trips, and tiny state amplitudes.
- All local measurement positions, including the seam, against dense
  canonical operators; completeness, trajectory replay, free energy,
  deferred compression, entropy opt-in, and nonmutation of the input state.

At N=4, the DMRG targets are `-3.236067977499789` in trivial flux and `-3.0`
in τ flux. Energy tolerance is `2e-9`; eigenvector residual tolerance is `1e-7`.
The runnable example is `exm/native_fibonacci_pbc.jl`.

Validation on 2026-10-26 with Julia 1.12.2, TensorKit 0.17.0, and MPSKit
0.13.13 passed all 289 native PBC assertions and all 203 existing encoded
backend/category/ITensor comparison assertions (492 total). The standalone
example and Documenter build, including doctests, also passed. The entire
repository test suite was not rerun for this change.

Commands:

```sh
julia --project -e 'include("test/test_NativeFibonacciPBC.jl"); include("test/test_MPSKitBackend.jl")'
julia --project exm/native_fibonacci_pbc.jl
```

## Observable boundaries

This backend represents one fixed flux sector at a time. Y tracking reports
its known eigenvalue, not an independent contraction of a native Y MPO.
There is no encoded/native basis isometry yet. Consequently native coefficients
are not returned by `mps_statevector`, and native states are rejected by
`ee_mps`. `native_schmidt_entropy` explicitly exposes the categorical canonical
bond entropy; its equivalence to encoded interval entropy or an annular RDM
has not been established.

Native hybrid, reference-qubit, Lyapunov-frame, and coherent mixed-flux
workflows remain unsupported. Category block sparsity is active, but a
performance advantage over the encoded backend has not been benchmarked.
