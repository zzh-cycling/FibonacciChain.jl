using FibonacciChain
using JLD
using Statistics

function post_selection(L::Int64, τ::Float64, D::Int64, sign::Int64=1)
    pbc = true
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    average_EElis=zeros(L-1)

    @time sample_measured_states, sample, sample_free_energy = Bulkpost_selection(L, τ, st, D, sign, pbc)
    EE_tlis = [ee(rdm_Fibo(L, collect(1:div(L,2)), state_t)) for state_t in sample_measured_states]
    final_state = sample_measured_states[end]
    average_EElis = eelis_Fibo_state(L, final_state)

    
    return average_EElis, EE_tlis, sample_free_energy
end



for τ in τlis[end]
    for L in 8:2:20
        D = get_system_params(τ, L)[1]
        average_EElis, EE_tlis, sample_free_energy = post_selection(L, τ, D, 1)
        save("exm/data/post_selection1/τ$(τ)/L$(L)_D$(div(D,L)).jld", "average_EElis", average_EElis, "EE_tlis", EE_tlis, "sample_free_energy", sample_free_energy)
    end
end


L=8
average_EElis, EE_tlis, sample_free_energy = post_selection(L, 0.6, 20L, 1)

