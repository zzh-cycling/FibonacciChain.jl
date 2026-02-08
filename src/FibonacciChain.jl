module FibonacciChain

using BitBasis, LinearAlgebra, SparseArrays, Arpack, Random
using ITensorMPS, ITensors

#  Exported types and structs
export FibonacciAnyon, IsingAnyon, OBFAnyon, AnyonModel, MeasureConfig, Measurement_outcome_bulk, Measurement_outcome_boundary, Measurement_outcome_mps
#  Basis and observable functions
export anyon_basis, anyon_ham, anyon_rdm, disjoint_rdm, topological_symmetry_basismap, mapst_sec2tot, anyon_rdm_sec, Fsymmetry_coef, topological_charge_operator
export ee, anyon_eelis, anyonladder_eelis, translation_matrix, inversion_matrix, braidingsqmap, spatial_correlation, temporal_correlation, ref_correlation, mutual_information, tri_mutual_information, trace_distance, fidelity
# ladder of anyons functions: two parallel anyon chains
export ladderChoi, ladderrdm, ladderbraidingsqmap, laddertranslationmap
export measure_basismap, measuremap, laddermeasuremap, measurement_enumeration, measurement_tree_visualization, transfer_matrix, boundary_evolution, bulk_evolution
export anyon_ham_sparse
# MPS-based functions
export anyon_mps_gst, measurement_operator_mps, apply_measurement_mps, initial_mps, evenparity_mps, mps_measurement_enumeration, ee_mps
# Reference qubits functions
export add_reference_qubits, add_reference_qubits_reset, reference_measuremap, reference_bulk_evolution, reference_boundary_evolution, reference_rdm, reference_evolution, build_extended_basis

include("Basis.jl") 
include("Observable.jl")
include("AnyonLadder.jl")
include("Measurement.jl")
include("FiboSparse.jl")
include("MPSMeasurement.jl")
include("ReferenceProbe.jl")

end
