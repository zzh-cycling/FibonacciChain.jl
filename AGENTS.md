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
