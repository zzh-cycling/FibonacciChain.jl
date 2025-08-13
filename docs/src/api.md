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
anyon_eelis
anyonladder_eelis
ee_mps
anyon_eelis_mps
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
laddermeasuremap
measurement_enumeration
measurement_tree_visualization
Boundary_measure
Boundarypost_selection
Bulkmeasure
Bulkpost_selection
generate_state
apply_measurement_layer!
```

## MPS Functions

```@docs
fibonacci_mps_ground_state
fibonacci_hamiltonian_mps
measurement_operator_mps
apply_measurement_mps
generate_state_mps
initial_mps
mps_measurement_enumeration
mps_boundary_measure
mps_bulk_measurement
```

## Reference Qubit Functions

```@docs
add_reference_qubits!
reference_measuremap
spatial_correlation
temporal_correlation
reference_generate_state
reference_apply_measurement_layer!
reference_evolution
reference_rdm
```
