using FibonacciChain
using JLD
using Statistics

function post_selection(L::Int64, τ::Float64, D::Int64, sign::Int64=1)
    pbc = true
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    average_EElis=zeros(L-1)

    @time sample_measured_states, samples, sample_weights = Bulkpost_selection(L, τ, st, D, sign, pbc)
    EE_tlis = [eelis_Fibo_state(L, state_t) for state_t in sample_measured_states]
    final_state = sample_measured_states[end]
    average_EElis = eelis_Fibo_state(L, final_state)

    
    return average_EElis, EE_tlis
end


L = 10
τ = 0.1
D = 60L
average_EElis, EE_tlis = post_selection(L, τ, D, 1)

save("exm/data/post_selection1_L$(L)_τ$(τ)_D$(div(D,L)).jld", "average_EElis", average_EElis, "EE_tlis", EE_tlis)