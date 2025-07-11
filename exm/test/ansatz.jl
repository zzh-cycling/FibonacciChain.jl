using FibonacciChain

N=10
initial_state = zeros(length(Fibonacci_basis(N)))
initial_state[end] = 1.0

τ = 1.0
measure_site = 3
sign=1 # measure to vacuum
measured_st = measuremap(N, τ, initial_state, measure_site, sign)

ϕ = (1 + sqrt(5)) / 2
eelis = eelis_Fibo_state(N, measured_st)

log(ϕ)