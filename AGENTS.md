# Repository Guidelines

## Project Structure & Module Organization

`src/FibonacciChain.jl` is the package entrypoint and owns the public exports. Core implementations live in focused files under `src/`: basis and Hamiltonian construction in `Basis.jl`, observables in `Observable.jl`, exact measurement dynamics in `Measurement.jl`, tensor-network methods in `MPSMeasurement.jl`, and hybrid evolution in `HybridEvolution.jl`. Tests mirror these areas in `test/test_*.jl` and are included by `test/runtests.jl`. Documentation sources are in `docs/src/`; runnable research workflows and Slurm scripts belong under `exm/`. Treat `exm/data/` and `docs/build/` as generated output, not source.

## Build, Test, and Development Commands

- `make init`: precompile the root Julia environment.
- `make test`: run the complete `Test.jl` suite through `Pkg.test()`.
- `julia --project -e 'using FibonacciChain, Test; include("test/test_Basis.jl")'`: run one test file while iterating.
- `make coverage`: run tests with Julia coverage enabled.
- `make docs`: build the Documenter site and execute doctests.
- `make serve`: serve documentation locally with LiveServer.

Avoid `make update` unless dependency upgrades are intentional; it updates the active environment.

## Coding Style & Naming Conventions

Use four-space indentation and follow the formatting of adjacent Julia code; the repository has no enforced formatter. Prefer multiple dispatch on `AnyonModel` and concrete basis types over large runtime type switches. Use `CamelCase` for types, `snake_case` for functions and variables, a trailing `!` for mutating functions, and a leading `_` for internal helpers. Add user-facing APIs to the explicit export list in `src/FibonacciChain.jl`. Preserve `BitBasis.BitStr` encoding and periodic-boundary conventions when changing basis logic.

## Testing Guidelines

Add focused `@testset`s to the matching `test/test_*.jl` file and register new files in `test/runtests.jl`. Use deterministic RNG seeds for stochastic dynamics. Check analytic identities with small-system exact calculations and use tolerances explicitly for floating-point or MPS comparisons. There is no stated coverage threshold, but every bug fix and public behavior change should include a regression test. Keep Documenter examples executable.

## Commit & Pull Request Guidelines

Recent commits use short imperative summaries such as `Add KW duality tracking` or `Fix initial QR`. Keep each commit scoped and mention the affected subsystem. Pull requests should explain the physical/numerical motivation, summarize API or data-format changes, list commands run, and link relevant issues or papers. Include plots or convergence evidence when numerical behavior changes; do not commit generated datasets unless they are intentionally reviewed artifacts.


## Data Storage, Compression, and Cleanup

### Minimal physical schema

- Store physical results and only the metadata required to interpret or replay them: system size, circuit geometry and depth, physical parameters, initial state, seeds, sample counts, time axes, observable definitions, and missing-data coverage.
- **Do not persist SHA256 digests, code-commit IDs, driver hashes, wall-clock timestamps, run dates, source-file paths used only for provenance, or verbose audit payloads in simulation or repacked data files.** Keep necessary verification information in external task receipts under `tmp/` or the cleanup operation directory.
- **Minimize hash verification (尽量减少 hash 校验).** Do not repeatedly hash whole directory trees. Prefer seed coverage, dimensions, physical invariants, deterministic readback, and representative replay. When destructive replacement needs integrity checks, hash each source at most once where practical, store it only in an external receipt, and avoid a second full-tree hash pass unless a specific inconsistency requires it.
- Do not store repeated static metadata independently in every seed. Use shared metadata, stacked physical arrays, and explicit per-seed variants or presence masks only when values actually differ. Do not require all seeds to have identical optional key sets.
- Avoid duplicated arrays and alias keys. Born measurement-only records normally need `samples`, not additional identical `outcomes` and all-true `measurement_mask`. All measurement-only circuits must omit `unitary_angles`; a large all-NaN placeholder is not useful data. Replay code must remain compatible with legacy placeholders.
- Do not store final state vectors, reference density matrices, continuation RNG objects, or continuation-only fields such as `extension_seed` when the physical replay record is sufficient. Recompute derived quantities when needed. Preserve small final observables required by the analysis even when the large state is omitted.
- Hybrid reference circuits store angles as Float32 for generation, evolution, storage, and replay. Document the approximate-replay error budget (observed observable errors around `1e-6`) and test it explicitly. Do not apply this conversion to OBF, which has no physical unitary angles.
- Old `.jld` trajectories may require `JLD`, not `JLD2`; inspect actual formats rather than inferring them from the extension.

### Choosing a reduction

- Decide explicitly whether the task requires exact replay, approximate replay, or only statistics. Use collect-only means/SEMs and necessary covariance/bootstrap data when records are unnecessary; use chunked compressed ensembles when replay is required. Do not silently turn a replay dataset into a collect-only dataset.
- Replace millions of tiny per-seed files with bounded chunked ensembles rather than merely gzipping every file. Preserve all required physical fields and avoid redundant per-seed plus ensemble copies of the same payload.
- A longer trajectory containing an exact dynamics prefix is **not** automatically a complete copy of the shorter file. Final states, I3, final subsystem entropies, reference density matrices, and other final observables can differ. Check these fields before claiming coverage; delete the short dataset only with explicit user approval after the required observables have been addressed.
- Files with the same seeds and time axes are not necessarily duplicates. Compare values and protocol semantics before deleting legacy or historical data. Similarly, identical directory names on two hosts do not prove content coverage, and different directory layouts do not prove its absence.
- Existing compact archives may be better than a new ensemble. Reuse them instead of creating a larger parallel copy, but preserve the required active replay path.

### Verification and deletion

- Use bounded-memory streaming conversion and atomic temporary outputs. Before replacing or deleting source data, read back the result and validate seed coverage, physical dimensions, every required field, and at least one real replay for the applicable initial-state/circuit family.
- Deletion requires the user's authorization for that target. On verification failure, retain sources and repair only the failed tasks. Do not submit duplicate Slurm arrays; retry only failed or timed-out indices with appropriate time limits.
- Check active jobs before modifying data and avoid excessive concurrent I/O. Full `du` and recursive scans over millions of small files are expensive; use existing manifests and targeted per-directory statistics instead of repeatedly traversing the whole tree.
- Keep operational logs, deletion manifests, and verification reports in task-specific ignored directories, not as extra fields in scientific datasets or as maintained root documentation.
