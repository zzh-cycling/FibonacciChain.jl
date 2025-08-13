```@meta
CurrentModule = FibonacciChain
```

# FibonacciChain.jl

**FibonacciChain.jl** is a Julia package for simulating 1D interacting anyons chain, particularly focusing on Fibonacci anyons and related topological systems.

## Overview

This package provides comprehensive tools for studying anyon chains based on the principle that, $(2+1) D$ bulk TQFT has correspondence to $(1+1) D$ boundary CFT. Similar to the Heisenberg model where singlet states have lower energy, fusion outcomes with trivial topological charge are energetically preferred. Then we can write down a interacting Hamiltonian.

### Supported Anyon Types

The package currently supports three types of anyonic systems:

- **Ising Anyons** (SU(2)₂): Majorana fermions with non-Abelian statistics
- **Fibonacci Anyons** (SU(2)₃): Universal anyons for quantum computation
- **Spin-1/2 Systems** (SU(2)∞): Regular spin systems

### Key Features

- **Exact Diagonalization**: Full quantum many-body calculations for small systems
- **Matrix Product States (MPS)**: Efficient simulation of larger systems using ITensors.jl
- **Measurement Protocols**: Quantum measurement and post-selection dynamics
- **Topological Properties**: Utilizing topological symmetry sector.
- **Anyon Operations**: Implementation of anyonic braiding and exchange statistics, fusion operation.

## Installation

```julia
using Pkg
Pkg.add("FibonacciChain")
```

## Quick Start

```julia
using FibonacciChain

# Create a Fibonacci anyon chain with N=6 sites
N = 6
basis = anyon_basis(N, true)  # periodic boundary conditions

# Generate the Hamiltonian
H = anyon_ham(N, true)

# Find ground state
eigenvals, eigenvecs = eigen(H)
ground_state = eigenvecs[:, 1]

# Calculate entanglement entropy profile
ee_profile = anyon_eelis(N, ground_state)
```

## Documentation Sections

```@contents
Pages = [
    "basis.md", 
    "observables.md",
    "measurements.md", 
    "mps.md",
    "examples.md",
    "api.md"
]
Depth = 2
```

## Theoretical Background

The package implements the anyon chain Hamiltonian described in [Phys. Rev. Lett. 98, 160409 (2007)](https://doi.org/10.1103/PhysRevLett.98.160409), where the energy scale is determined by favoring trivial fusion channels in three-body interactions.

For Fibonacci anyons, the local Hamiltonian terms involve three consecutive anyons and depend on their fusion outcomes according to the golden ratio φ = (1+√5)/2, which characterizes the Fibonacci fusion algebra.

## API Reference

```@autodocs
Modules = [FibonacciChain]
```

## Index

```@index
```
