# API Reference

This page provides a complete reference for all exported functions in FibonacciChain.jl.

## Basis and Hamiltonian Functions

```@docs
anyon_basis
anyon_ham
anyon_ham_sparse
```

## Reduced Density Matrices

```@docs
anyon_rdm
anyon_rdm_sec
disjoint_rdm
mapst_sec2tot
```

## Entanglement and Observables

```@docs
ee
ee_mps
anyon_eelis
```

## Symmetry Operations

```@docs
translation_matrix
inversion_matrix
braidingsqmap
topological_symmetry_basismap
Fsymmetry_coef
```

## Ladder System Functions

```@docs
ladderChoi
ladderrdm
ladderbraidingsqmap
laddertranslationmap
```

## Measurement Functions

```@docs
measure_basismap
measuremap
measurement_enumeration
measurement_tree_visualization
transfer_matrix
boundary_evolution
bulk_evolution
```

## MPS Functions

```@docs
measurement_operator_mps
initial_mps
mps_measurement_enumeration
```

## Reference Qubit Functions

```@docs
add_reference_qubits
add_reference_qubits_reset
reference_measuremap
reference_boundary_evolution
reference_bulk_evolution
reference_rdm
reference_evolution
```
