# Copilot instructions for FibonacciChain.jl

## Build, test, and docs commands

- Precompile project dependencies:
  - `make init`
- Update dependencies and precompile:
  - `make update`
- Run full test suite:
  - `make test`
  - Equivalent: `julia --project -e 'using Pkg; Pkg.test()'`
- Run tests with coverage:
  - `make coverage`
- Run a single test file (repo pattern is include-based):
  - `julia --project -e 'using FibonacciChain, Test; include("test/test_Basis.jl")'`
  - Replace with another file from `test/` as needed (for example `test/test_Measurement.jl`).
- Build docs:
  - `make docs`
  - Equivalent: `julia --project=docs docs/make.jl`
- Serve docs locally:
  - `make serve`

There are currently no dedicated lint/formatter commands in `Makefile` or GitHub Actions workflows.

## High-level architecture

`src/FibonacciChain.jl` is the API entrypoint: it exports the public surface and composes the package by including focused modules.

- `src/Basis.jl`
  - Core anyon types (`FibonacciAnyon`, `IsingAnyon`, `OBFAnyon`) and `AnyonModel`.
  - Basis construction (`anyon_basis`), Hamiltonian construction (`anyon_ham`), reduced density matrices, and topological-symmetry operators.
- `src/FiboSparse.jl`
  - Sparse Hamiltonian path (`anyon_ham_sparse`) for larger systems.
- `src/Observable.jl`
  - Entropy/information measures, symmetry maps, braiding operators, and correlation functions.
- `src/Measurement.jl`
  - Measurement maps and protocols on basis states (`measuremap`, `measurement_enumeration`, `transfer_matrix`, `boundary_evolution`, `bulk_evolution`).
  - Configuration and outputs are structured via typed structs (for example `MeasureConfig`).
- `src/MPSMeasurement.jl`
  - ITensors/ITensorMPS-backed path (MPS state prep, measurement operators, MPS enumeration/evolution).
- `src/ReferenceProbe.jl`
  - Reference-qubit extensions layered on top of the core measurement/basis machinery.
- `src/AnyonLadder.jl`
  - Two-chain ladder-specific mappings and RDM utilities.

Typical flow across files:
1. Define `AnyonModel` + construct basis (`Basis.jl`).
2. Build Hamiltonian using dense/sparse/MPS path (`Basis.jl`, `FiboSparse.jl`, `MPSMeasurement.jl`).
3. Compute observables (`Observable.jl`) and/or run measurement evolution (`Measurement.jl`, `MPSMeasurement.jl`).
4. Optionally use reference-qubit workflow (`ReferenceProbe.jl`) or ladder workflow (`AnyonLadder.jl`).

## Key repository-specific conventions

- Multiple dispatch on `AnyonModel{AT}` is central; behavior differs by anyon type. Prefer extending existing typed methods over adding type switches.
- Basis/state encoding heavily uses `BitBasis.BitStr`; preserve bit-encoding assumptions when modifying basis or mapping logic.
- Public API is curated through explicit `export` entries in `src/FibonacciChain.jl`; keep exports aligned with user-facing additions.
- Internal helpers are commonly prefixed with `_`; mutating functions use `!`.
- Tests are organized as include-based per-file testsets from `test/runtests.jl`. Follow that pattern when adding tests.
- Documentation and verification rely on Documenter doctests (`docs/make.jl` sets `doctest=true`, and CI runs doctests explicitly). Keep doc examples executable.

## Project context from README/docs

- The package models interacting anyons with currently supported families corresponding to `SU(2)_2` (Ising), `SU(2)_3` (Fibonacci), and `SU(2)_∞` (ordinary spin-1/2 basis via `OBFAnyon`).
