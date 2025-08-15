using FibonacciChain
using JLD
using Statistics
include("../FitEntEntScal.jl")

function post_selection(L::Int64, τ::Float64, D::Int64, sign::Int64=1)
    pbc = true
    st=zeros(length(anyon_basis(L)))
    st[1] = 1.0
    average_EElis=zeros(L-1)

    @time sample_measured_states, sample, sample_free_energy = Bulkpost_selection(L, τ, st, D, sign, pbc)
    EE_tlis = [ee(anyon_rdm(L, collect(1:div(L,2)), state_t)) for state_t in sample_measured_states]
    final_state = sample_measured_states[end]
    average_EElis = anyon_eelis(L, final_state)

    
    return average_EElis, EE_tlis, sample_free_energy
end

function get_system_params(τ, L)
    table = Dict(
            atanh(0.1)  => (2500L, 1000, 1500L),
            atanh(0.2)  => (500L,  100, 250L),
            atanh(0.3)  => (120L,  48, 100L),
            atanh(0.4)  => (100L,  40, 80L),
            atanh(0.5)  => (80L,   32, 40L),
            atanh(0.6)  => (45L,   20, 30L),
            log(1 + √2) => (35L,   14, 20L),
            atanh(0.8)  => (25L,   10, 10L),
            atanh(0.9)  => (8L,    4, 4L),
            atanh(0.95) => (8L,    4, 4L),
            atanh(0.999)=> (5L,    2, 2L),
        )
    D, step, start = get(table, τ, (5L, 2, 2L))
    inds = collect(1:step:D)
    avg_range = start:D-5
    return D, inds, avg_range
end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 

for τ in τlis
    @show τ
    for L in 8:2:20
        D = get_system_params(τ, L)[1]
        @show L
        average_EElis, EE_tlis, sample_free_energy = post_selection(L, τ, D, 1)
        save("exm/data/post_selection1/τ$(τ)/L$(L)_D$(div(D,L)).jld", "average_EElis", average_EElis, "EE_tlis", EE_tlis, "sample_free_energy", sample_free_energy)
    end
end


# L=20
# average_EElis, EE_tlis, sample_free_energy = post_selection(L, 1000.0, 5L, 0)
