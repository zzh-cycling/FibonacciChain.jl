using FibonacciChain
using JLD
using Statistics
using Random

sample = load("exm/data/Bulk_measure/Samples_monitored_dynamics/L12/τ0.8813735870195429/D35_Samples1.jld", "sample")

N = 12
pbc = true
measure_class = :Fibo
initial_state = zeros(length(Fibonacci_basis(BitStr{N, Int}, pbc, measure_class=measure_class)))
initial_state[1] = 1.0 # initial state is all zero state

statelis = generate_state(0.5, initial_state, sample, true, temp= true)

final_st= statelis[end]
spatial_corr = spatial_correlation(N, final_st, 1, 5, pbc)


S1_state = add_reference_qubit!(initial_state, 1, pbc, measure_class=measure_class)
S2_state = add_reference_qubit!(initial_state, 2, pbc, measure_class=measure_class)
S12_state = add_reference_qubit!(initial_state, [1, 2], pbc, measure_class=measure_class)
S1_rdm = disjoint_rdm(BitStr{N, Int}, BitStr{N, Int}, [1], [2], S1_state, pbc; measure_classA=measure_class, measure_classB=measure_class)
S2_rdm = disjoint_rdm(BitStr{N, Int}, BitStr{N, Int}, [1], [2], S2_state, pbc; measure_classA=measure_class, measure_classB=measure_class)
S12_rdm = disjoint_rdm(BitStr{N, Int}, BitStr{N, Int}, [1, 2], [1, 2], S1_state, pbc; measure_classA=measure_class, measure_classB=measure_class)
S1_ee = ee(S1_rdm)
S2_ee = ee(S2_rdm)
S12_ee = ee(S12_rdm)
I = S1_ee + S2_ee- S12_ee