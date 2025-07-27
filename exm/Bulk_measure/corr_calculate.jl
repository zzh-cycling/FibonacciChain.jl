using FibonacciChain
using JLD
using Statistics
using Random

ensemble_seed, sample = load("exm/data/Bulk_measure/Ising/monitored_EE_FEdynamics_L8_τ0.2027325540540822_D20.jld", "ensemble_seed", "sample")

rng = MersenneTwister(ensemble_seed)
N = 8
pbc = true
measure_class = :Fibo
initial_state = zeros(FibonacciChain.Fibonacci_basis(BitStr{N, Int}, pbc, measure_class=measure_class))
initial_state[1] = 1.0 # initial state is all zero state

statelis = generate_state(τ, st, samples, true, temp= true)


S1_state = add_reference_qubit!(initial_state, 1, pbc, measure_class=measure_class)
S2_state = add_reference_qubit!(initial_state, 2, pbc, measure_class=measure_class)
S1_rdm = disjoint_rdm(BitStr{N, Int}, BitStr{N, Int}, [1], [2], S1_state, pbc; measure_classA=measure_class, measure_classB=measure_class)
S2_rdm = disjoint_rdm(BitStr{N, Int}, BitStr{N, Int}, [1], [2], S2_state, pbc; measure_classA=measure_class, measure_classB=measure_class)
S1_ee = ee(S1_rdm)
S2_ee = ee(S2_rdm)
I = S1_ee + S2_ee- S12_ee