using FibonacciChain    
using BenchmarkTools
using Profile
using Random 

model_X = AnyonModel(OBFAnyon(), 12, pbc=true, measure_operator=:X)
@btime measuremap(model_X, 0.5, ones(2^12), 1, false)

model_XZZ = AnyonModel(OBFAnyon(), 12, pbc=true, measure_operator=:XZZ)
@btime measuremap(model_XZZ, 0.5, ones(2^12), 1, false)
@profile for _ = 1:1000 measuremap(model_XZZ, 0.5, ones(2^12), 1, false) end
@code_warntype measuremap(model_XZZ, 0.5, ones(2^12), 1, false)
@code_warntype FibonacciChain._apply_result(model_XZZ, 0.5, BitStr{12, Int}(1), 1, false)
@code_warntype FibonacciChain.measure_basismap(model_XZZ, 0.5, BitStr{12, Int}(1), 1, false)

t = 10
measure_config = MeasureConfig(τ=0.5, t₂=t, mode=:sample)
samples = BitMatrix(zeros(Int8, 14t, 10))
model = AnyonModel(OBFAnyon(), 10, λ = 0.5, pbc=true)
st = zeros(length(anyon_basis(model)))
st[1] = 1.0
@btime bulk_evolution(model, st, measure_config, samples)
@profile for _ = 1:5 bulk_evolution(model, st, measure_config, samples) end

@code_warntype bulk_evolution(model, st, measure_config, samples)


L=10
rng = MersenneTwister(1)
model = AnyonModel(FibonacciAnyon(), L; pbc=true)
ψ, sites = initial_mps(L)
config = MeasureConfig(τ=log(1+√2), mode=:Born, t₂=10*L, rng=rng, cutoff=1e-12, maxdim=20)
# Pre-allocate arrays for on-the-fly EE computation
last_state = Ref{MPS}()  # to capture the final state

# Callback: compute EE immediately and discard the MPS
function on_state_callback(period::Int, state::MPS)
    halfchain_EE_tlis[period] = ee_mps(state, div(L, 2))
    last_state[] = state  # only keeps reference to the latest state
end
halfchain_EE_tlis = zeros(10 * L)
final_EElis = zeros(L - 1)
@time mps_mo = bulk_evolution(model, sites, ψ, config; on_state=on_state_callback)