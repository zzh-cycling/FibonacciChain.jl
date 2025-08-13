# Basis Functions

This section covers the fundamental basis generation and manipulation functions for anyon chains.

## Anyon Basis Generation

```@docs
anyon_basis
```

The `anyon_basis` function generates the computational basis states for various anyon models:

### Fibonacci Anyons
For Fibonacci anyons, the basis consists of configurations that satisfy the fusion constraints. Each site can be in state `0` (vacuum) or `1` (τ particle), with the constraint that no two adjacent τ particles can fuse to vacuum (forbidden configurations are excluded).

### Ising Anyons  
For Ising anyons (Majorana fermions), the basis includes all possible σ and 1 configurations, representing the two types of Ising anyons.

### Usage Examples

```julia
using FibonacciChain

# Generate Fibonacci basis for 6 sites with PBC
N = 6
basis_pbc = anyon_basis(N, true, measure_class=:Fibo)

# Generate Ising anyon basis  
basis_ising = anyon_basis(N, true, measure_class=:IsingX)

# Generate basis in momentum sector k=0
basis_k0, rep_dict = anyon_basis(N, 0)
```

## Hamiltonian Construction

```@docs
anyon_ham
anyon_ham_sparse
```

## Reduced Density Matrices

```@docs
anyon_rdm
disjoint_rdm
```

## Topological Symmetries

```@docs
topological_symmetry_basismap
Fsymmetry_coef
```

The topological symmetry operations are crucial for understanding the anyonic nature of the system. They encode how the fusion outcomes transform under topological operations.
