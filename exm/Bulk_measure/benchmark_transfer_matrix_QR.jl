using FibonacciChain
using LinearAlgebra
using KrylovKit

include(joinpath(@__DIR__, "config.jl"))

# Only for ps
# function ps_spectrum()
#     L = 8
#     τ = atanh(0.95)
#     # sample = BitMatrix(ones(Int8, 2, div(8, 2)))
#     sample = BitMatrix(vcat(ones(Int8, 1, div(L, 2)), zeros(Int8, 1, div(L, 2))))
#     model = AnyonModel(FibonacciAnyon(), L)
#     T = transfer_matrix(model, τ, sample)
#     energy = eigen(T).values
#     # For non-symmetric matrices, eigenvalues are NOT sorted.
#     # We must explicitly sort by magnitude to pick the largest ones.
#     return -log.(energy[end-9:end])
# end

# This file is used to benchmark Born case (non-hermitian matrix) spectrum
# Method 1: build full tf for ED spectrum
function transfer_matrix_ED(
        L::Int,
        τ_idx::Int,
        sample::BitMatrix
    )
    τ = τlis[τ_idx]

    D = 64
    t = div(D, 2)
    
    model = fib_model(L)
    basis = anyon_basis(model)
    l = length(basis)

    TM = zeros(Float64, l, l)

    for i = 1:l
        st = zeros(Float64, l)
        st[i] = 1.0

        TM[:, i] = sample_evolution_unnormalized(
                model, st, sample; τ = τ, enable_τ_eff = false
            )
    end

    return TM
end

function transfer_matrix_QR(
        L::Int,
        τ_idx::Int,
        sample::BitMatrix
    )
    τ = τlis[τ_idx]

    D = 64
    t = div(D, 2)
    
    model = fib_model(L)
    basis = anyon_basis(model)
    l = length(basis)
    k = min(10, l)

    # Initialize k product states (basis vectors)
    states = zeros(Float64, l, k)
    for i in 1:k
        states[i, i] = 1.0
    end

    free_energies = zeros(k, t)
    for step in 1:t
        sample_layer = sample[2step - 1:2step, :]

        # Apply transfer matrix to each state
        for i in 1:k
            states[:, i] = sample_evolution_unnormalized(
                model, states[:, i], sample_layer; τ = τ, enable_τ_eff = false
            )
        end

        # # QR orthogonalize
        
        Q, R = qr(states)
        states = Q[:, 1:k]
        free_energies[:, step] = -log.(abs.(diag(R)))
    end

    return free_energies
end

# 定义线性映射
free_energies = transfer_matrix_QR(8, 10, sample);

# sample = BitMatrix(ones(Int8, 64, div(8, 2)))
# sample = BitMatrix(vcat(ones(Int8, 1, div(L, 2)), zeros(Int8, 1, div(L, 2))))
sample = load("exm/data/Bulk_measure/monitored_dynamics/L8/gammaind10/t4_samples1.jld", "sample")
sample_free_energy = load("exm/data/Bulk_measure/monitored_dynamics/L8/gammaind10/t4_samples1.jld","sample_free_energy")

@show [sample_free_energy[2:2:end] free_energies[1,:]]


# # Method 2: subspace. Method 3: Krylov
# @show [-log.(eigvals(transfer_matrix_ED(8, 10, sample))[end-9:end]) ./32 -log.(transfer_matrix_subspace(model, τlis[10], sample)[1:10, end]) energies]
# # @show [sort(transfer_matrix_Born(8, 10, 10), rev = true) sort(ps_spectrum(sample), rev=true)]
