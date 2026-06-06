using FibonacciChain
using LinearAlgebra
using KrylovKit

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[7] = log(1 + √2)  # atanh(1/√2) = log(1 + √2)
τlis[end] = 1000.0     # Last value is for γ=1

function transfer_matrix_Born(
        L::Int,
        τ_idx::Int,
        sample::BitMatrix
    )
    τ = τlis[τ_idx]

    D = 64
    t = div(D, 2)
    
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    basis = anyon_basis(model)
    l = length(basis)
    k = min(10, l)

    # Initialize k product states (basis vectors)
    states = zeros(Float64, l, k)
    for i in 1:k
        states[i, i] = 1.0
    end

    for step in 1:t
        sample_layer = sample[2step - 1:2step, :]

        # Apply transfer matrix to each state
        for i in 1:k
            states[:, i] = sample_evolution_unnormalized(
                model, states[:, i], sample_layer; τ = τ, enable_τ_eff = false
            )
        end

        # Normalize each state
        probs = vec(sum(abs2, states, dims = 1))
        for i in 1:k
            norm_i = sqrt(probs[i])
            if norm_i > 0
                states[:, i] ./= norm_i
            end
        end

        # QR orthogonalize
        Q = qr(states).Q
        states = Q[:, 1:k]
    end

    # Rayleigh-Ritz projection onto the converged subspace
    # H[i,j] = <states[:, i] | T | states[:, j]>
    H = zeros(k, k)
    sample_layer = sample[2t - 1:2t, :]
    for j in 1:k
        w = sample_evolution_unnormalized(
            model, states[:, j], sample_layer; τ = τ, enable_τ_eff = false
        )
        for i in 1:k
            H[i, j] = dot(states[:, i], w)
        end
    end
    ritz_values = eigvals(H)
    return -log.(ritz_values)
end

function ps_spectrum(sample::BitMatrix)
    L = 8
    τ = atanh(0.95)

    model = AnyonModel(FibonacciAnyon(), L)
    T = transfer_matrix(model, τ, sample)
    energy = eigen(T).values
    # For non-symmetric matrices, eigenvalues are NOT sorted.
    # We must explicitly sort by magnitude to pick the largest ones.
    return -log.(energy[end-9:end])
end

# 定义线性映射
function Tmap(v)
    L = 8
    model = AnyonModel(FibonacciAnyon(), L)
    τ = atanh(0.95)
    # sample = BitMatrix(ones(Int8, 2, div(L, 2)))
    sample = load("exm/data/Bulk_measure/monitored_dynamics/L8/gammaind10/t4_samples1.jld", "sample")

    return sample_evolution_unnormalized(model, v, sample; τ=τ, enable_τ_eff=false)
end


sample = BitMatrix(ones(Int8, 64, div(8, 2)))
# sample = BitMatrix(vcat(ones(Int8, 1, div(L, 2)), zeros(Int8, 1, div(L, 2))))
# sample = load("exm/data/Bulk_measure/monitored_dynamics/L8/gammaind10/t4_samples1.jld", "sample")
# 求解前 k 个最大模本征值
st = zeros(47);st[1] = 1
vals, vecs, info = eigsolve(Tmap, st, 10, :LM; tol=1e-100, krylovdim=200);
# eigsolve(Tmap, st, 10, :LR; tol=1e-50,krylovdim=80, maxiter=2000);
energies = real.(-log.(vals))[1:10]/size(sample, 1)*2


@show [transfer_matrix_Born(8, 10, sample) ps_spectrum(sample) sort(energies, rev=true)]
# @show [sort(transfer_matrix_Born(8, 10, 10), rev = true) sort(ps_spectrum(sample), rev=true)]
