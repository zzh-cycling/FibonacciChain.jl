using FibonacciChain
using JLD
using Statistics
include("../FitEntEntScal.jl")

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

function get_system_params(τ, L)
    if τ == log(1 + √2)
        D = 35L
        inds = collect(1:14:D)
        avg_range = 20L:D-5
    elseif τ == atanh(0.1)
        D = 2500L
        inds = collect(1:1000:D)
        avg_range = 1500L:D-5
    elseif τ == atanh(0.2)
        D = 500L
        inds = collect(1:100:D)
        avg_range = 250L:D-5
    elseif τ == atanh(0.3)
        D = 120L
        inds = collect(1:48:D)
        avg_range = 100L:D-5
    elseif τ == atanh(0.4)
        D = 100L
        inds = collect(1:40:D)
        avg_range = 80L:D-5
    elseif τ == atanh(0.5)
        D = 80L
        inds = collect(1:32:D)
        avg_range = 40L:D-5
    elseif τ == atanh(0.6)
        D = 45L
        inds = collect(1:20:D)
        avg_range = 30L:D-5
    elseif τ == atanh(0.8)
        D = 25L
        inds = collect(1:10:D)
        avg_range = 10L:D-5
    elseif τ == atanh(0.9) || τ == atanh(0.95)
        D = 8L
        inds = collect(1:4:D)
        avg_range = 4L:D-5
    elseif τ == atanh(0.999)
        D = 5L
        inds = collect(1:2:D)
        avg_range = 2L:D-5
    else
        D = 5L  # Default value for τ=1000.0
        inds = collect(1:2:D)
        avg_range = 2L:D-5
    end
    
    return D, inds, avg_range
end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 

for τ in τlis[end]
    for L in 8:2:20
        D = get_system_params(τ, L)[1]
        @show L
        average_EElis, EE_tlis, sample_free_energy = post_selection(L, τ, D, 0)
        save("exm/data/post_selection0/τ$(τ)/L$(L)_D$(div(D,L)).jld", "average_EElis", average_EElis, "EE_tlis", EE_tlis, "sample_free_energy", sample_free_energy)
    end
end


L=20
average_EElis, EE_tlis, sample_free_energy = post_selection(L, 1000.0, 5L, 0)
