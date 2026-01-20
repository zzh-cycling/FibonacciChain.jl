using FibonacciChain    
using BenchmarkTools
using Profile

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