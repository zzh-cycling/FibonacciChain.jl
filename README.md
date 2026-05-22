# FibonacciChain

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://zzh-cycling.github.io/FibonacciChain.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://zzh-cycling.github.io/FibonacciChain.jl/dev/)
[![Build Status](https://github.com/zzh-cycling/FibonacciChain.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/zzh-cycling/FibonacciChain.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/zzh-cycling/FibonacciChain.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/zzh-cycling/FibonacciChain.jl)

This package is designed to simulate interacting anyons based on the idea that, like in the Heisenberg model where singlet states have lower energy, here we choose the trivial fusion outcome as the lower energy state. Currently we support three types of anyons: $SU(2)_2$, $SU(2)_3$, and $SU(2)_\infty$, corresponding to Ising anyons / Majorana fermions, Fibonacci anyons, and normal spin-$1/2$ basis, respectively.

Ref: Phys. Rev. Lett. 98, 160409 [DOI](https://doi.org/10.1103/PhysRevLett.98.160409)

## Installation

You can install the package via the Julia REPL:

```julia
using Pkg
Pkg.add("FibonacciChain")
```

Or in the Pkg REPL mode (press `]`):

```julia
pkg> add FibonacciChain
```

## Quick Start

```julia
using FibonacciChain

# Create a Fibonacci anyon chain with N=6 sites
N = 6
model = AnyonModel(FibonacciAnyon(), N; pbc=true)
basis = anyon_basis(model)

# Generate the Hamiltonian
H = anyon_ham(model)

# Find the ground state
eigenvals, eigenvecs = eigen(H)
ground_state = eigenvecs[:, 1]

# Calculate entanglement entropy profile
ee_profile = anyon_eelis(model, ground_state)
```

For more examples and API details, see the [documentation](https://zzh-cycling.github.io/FibonacciChain.jl/dev/).

## Features

- **Exact Diagonalization**: Full quantum many-body calculations for small systems.
- **Matrix Product States (MPS)**: Efficient simulation of larger systems using [ITensors.jl](https://itensor.github.io/ITensors.jl/stable/).
- **Measurement Protocols**: Quantum measurement and post-selection dynamics.
- **Topological Properties**: Topological symmetry sector analysis.
- **Anyon Operations**: Anyonic braiding, exchange statistics, and fusion operations.

## Citation

If you use this package in your research, please cite the original reference:

> Levin, M., & Wen, X. G. (2007). Phys. Rev. Lett. 98, 160409.
