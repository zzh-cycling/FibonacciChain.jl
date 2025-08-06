module FibonacciChain

using BitBasis, LinearAlgebra, SparseArrays, Arpack, Random
using ITensorMPS, ITensors

export Fibonacci_Ham, Fibonacci_basis, rdm_Fibo, disjoint_rdm, reference_rdm, topological_symmetry_basismap, mapst_sec2tot, rdm_Fibo_sec
export ee, eelis_Fibo_state, eelis_Fiboladder_state, translation_matrix, inversion_matrix, braidingsqmap
export ladderChoi, ladderrdm, ladderbraidingsqmap, laddertranslationmap
export measure_basismap, measuremap, laddermeasuremap, measurement_enumeration, measurement_tree_visualization, Boundary_measure, Boundarypost_selection, Bulkmeasure, Bulkpost_selection, generate_state
export Fibonacci_Ham_sparse
# MPS-based functions
export fibonacci_mps_ground_state, fibonacci_hamiltonian_mps, measurement_operator_mps, apply_measurement_mps, generate_state_mps
export initial_mps, mps_measurement_enumeration, mps_boundary_measure, mps_bulk_measurement, ee_mps, eelis_Fibo_mps
export add_reference_qubits!, reference_measuremap, spatial_correlation, temporal_correlation, reference_generate_state

include("Basis.jl") 
include("Observable.jl")
include("LadderFibo.jl")
include("Measurement.jl")
include("FiboSparse.jl")
include("MPSMeasurement.jl")
include("ReferenceProbe.jl")

end
