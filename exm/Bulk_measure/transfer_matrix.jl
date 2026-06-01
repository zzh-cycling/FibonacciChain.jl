using FibonacciChain
using LinearAlgebra
using KrylovKit
using JLD, JLD2
using Statistics
using Plots
using LaTeXStrings
using LsqFit

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[7] = log(1 + √2)  # atanh(1/√2) = log(1 + √2)
τlis[end] = 1000.0  # Last value is for γ=1

function get_cfg_params_Born(ind, L)
    cfg = Dict(
        1 => (2500L, 1000, 750L),
        2 => (500L, 100, 120L),
        3 => (200L, 40, 80L),
        4 => (100L, 40, 40L),
        5 => (80L, 32, 20L),
        6 => (45L, 20, 15L),
        7 => (35L, 14, 10L),
        8 => (25L, 10, 5L),
        9 => (8L, 4, 2L),
        10 => (8L, 4, 2L),
        11 => (5L, 2, 1L),
    )
    D, step, start = get(cfg, ind, (5L, 2, L))
    inds = collect(1:step:div(D, 2))
    avg_range = start:(div(D, 2)-5)
    return D, inds, avg_range
end

function get_mps_params_Born(τind, L)
    cfg = if L <= 32
        Dict(
            1 => (1250, 1000, 600),
            2 => (250, 100, 150),
            3 => (40, 48, 30),
            4 => (28, 40, 30),
            5 => (40, 32, 24),
            6 => (22, 20, 15),
            7 => (10, 14, 7),
            8 => (12, 10, 8),
            9 => (3, 4, 3),
            10 => (4, 4, 2.5),
            11 => (3, 2, 2),
        )
    else
        Dict(
            1 => (700, 1000, 500),
            2 => (150, 100, 100),
            3 => (40, 48, 30),
            4 => (28, 40, 22),
            5 => (20, 32, 16),
            6 => (12, 20, 10),
            7 => (10, 14, 7),
            8 => (7, 10, 5.5),
            9 => (3, 4, 2.5),
            10 => (2, 4, 1.5),
            11 => (2, 2, 1.5),
        )
    end
    t, step, start = get(cfg, τind, (2, 2, 1))
    inds = collect(1:step:(t*L))
    avg_range = Int(start*L):2:(Int(t*L)-4)
    return t, inds, avg_range
end

function get_default_chi_Born(ind, L)
    if L == 32
        chi64_table = Dict(3 => 150, 4 => 150, 7 => 150, 9 => 200)
        return get(chi64_table, ind, 80)
    elseif L == 48
        chi48_table = Dict(1 => 150)
        return get(chi48_table, ind, 200)
    elseif L == 128 && ind == 10
        return 300
    elseif L == 64
        chi64_table = Dict(
            3 => 250,
            4 => 250,
            5 => 300,
            6 => 175,
            7 => 250,
            8 => 300,
            9 => 200,
            10 => 250,
        )
        return get(chi64_table, ind, 110)
    end
end

# ---------------------------------------------------------------------------
# Unnormalized measurement-layer application (no state normalization)
# ---------------------------------------------------------------------------
function _apply_measurement_layer_unnormalized(
    anyon_model::AnyonModel{AT},
    τ::Float64,
    state::Vector{T},
    layer_sample::BitVector;
    layer_idx::Int64 = 1,
) where {T, AT}

    measurement_sites, measure_anyon_model, measurement_strength =
        FibonacciChain._obtain_measurement_config(anyon_model, layer_idx, τ)

    mop = anyon_model.measure_operator
    N = anyon_model.N
    pbc = anyon_model.pbc
    if mop ∈ (:Ferro, :Antiferro)
        @assert pbc || 2 <= measurement_sites[1] || measurement_sites[end] <= N-1 "Index idx must be in [2, N-1] for open BC (Fibonacci)"
    elseif mop == :ZZ
        @assert pbc || 1 <= measurement_sites[1] || measurement_sites[end] <= N-1 "Index idx must be in [1, N-1] for open BC (ZZ)"
    elseif mop ∈ (:X, :Z, :reset, :resetFibo)
        @assert pbc || 1 <= measurement_sites[1] || measurement_sites[end] <= N "Index idx must be in [1, N] for open BC (X)"
    elseif mop ∈ (:XZZ, :ZZX)
        @assert pbc || 1 <= measurement_sites[1] || measurement_sites[end] <= N-2 "Index idx must be in [1, N-2] for open BC (OBF)"
    end

    basis = FibonacciChain.anyon_basis(measure_anyon_model)
    l = length(basis)
    buf = Vector{T}(undef, l)
    current_state = copy(state)

    for (idx, sign) in enumerate(layer_sample)
        fill!(buf, zero(T))
        FibonacciChain._measuremap_impl!(
            buf,
            basis,
            measure_anyon_model,
            measurement_strength,
            current_state,
            measurement_sites[idx],
            sign,
        )
        current_state, buf = buf, current_state
    end

    return current_state
end

# ---------------------------------------------------------------------------
# Unnormalized bulk evolution: apply a fixed sample without normalizing.
# This gives a linear map T(s) acting on the state vector.
# ---------------------------------------------------------------------------
function bulk_evolution_unnormalized(
    model::AnyonModel{AT},
    state::Vector{ET},
    samples::BitMatrix;
    τ::Float64 = 1.0,
    enable_τ_eff::Bool = true,
) where {ET,AT}

    n_layers = FibonacciChain.layers_per_period(model.anyon_type)
    D_layers, n_cols = size(samples)

    @assert D_layers % n_layers == 0 "Number of layers $D_layers must be divisible by $n_layers"
    Δt = D_layers ÷ n_layers

    n_cols == FibonacciChain._samples_per_layer(model) ||
        error("sample size spatial dimension must be $n_cols, got $(size(samples, 2))")

    current_state = copy(state)

    for period = 1:Δt
        for layer = 1:n_layers
            global_layer_idx = (period - 1) * n_layers + layer
            τ_current = (period == Δt && layer == n_layers && enable_τ_eff) ? τ/2 : τ

            col_indices = FibonacciChain._get_sample_column_indices(model, global_layer_idx)
            layer_sample = BitVector(samples[global_layer_idx, col_indices])

            current_state = _apply_measurement_layer_unnormalized(
                model,
                τ_current,
                current_state,
                layer_sample;
                layer_idx = global_layer_idx,
            )
        end
    end

    return current_state
end

# ---------------------------------------------------------------------------
# Build the transfer matrix for a fixed sample.
# Column i = T(s) * e_i  (evolved basis vector, unnormalized).
# ---------------------------------------------------------------------------
function transfer_matrix(
    L::Int,
    τ_idx::Int,
    index::Int;
    sample_path::Union{String,Nothing} = nothing,
)
    τ = τlis[τ_idx]

    # Default path matches the original hard-coded convention.
    # Override with sample_path=... if your data lives elsewhere.
    path = something(
        sample_path,
        "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τ_idx)/t4_samples$(index).jld"
    )

    data = load(path)
    sample = data["sample"]

    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    model_basis = anyon_basis(model)
    l = length(model_basis)

    T = zeros(l, l)
    for i in 1:l
        st = zeros(l)
        st[i] = 1.0
        final_st = bulk_evolution_unnormalized(model, st, sample; τ = τ)
        T[:, i] = final_st
    end

    return T
end

# Example usage (uncomment to run):
λlis = Vector{Vector{ComplexF64}}(undef, 1000);
for i in 1:1000
    TF = transfer_matrix(12, 10, i)
    @show i
    energy, states = eigen(TF)
    λlis[i] = energy
end
save("transfer_matrix_spectrum_L12.jld", "λlis", λlis)

# Llis = collect(8:2:12)

function plot_scaling_dimension(Llis)
    λ₀Llis = zeros(length(Llis))
    λ₀_errLlis = zeros(length(Llis))
    λ₁Llis = zeros(length(Llis))
    λ₁_errLlis = zeros(length(Llis))
    for (idx, L) in enumerate(Llis)
        data = load("transfer_matrix_spectrum_L$(L).jld", "λlis")
        λ₀lis = [real(λ[end]) for λ in λlis]
        λ₁lis = [real(λ[end-1]) for λ in λlis]
        λ₀Llis[idx] = mean(-log.(λ₀lis)/4L)
        λ₀_errLlis[idx] = std(-log.(λ₀lis)/4L)/sqrt(length(λ₀lis))
        λ₁Llis[idx] = mean(-log.(λ₁lis)/4L)
        λ₁_errLlis[idx] = std(-log.(λ₁lis)/4L)/sqrt(length(λ₁lis))
    end
    
    fit_model(x, p) = 2π * p[1] .*x .+ p[2]
    
    fig = scatter(1 ./Llis, λ₀Llis, yerr = λ₀_errLlis, label = false, xticks = (vcat(0, 1 ./Llis), latexstring.(vcat("0", ["\\frac{1}{$(i)}" for i in Llis]))), markershape=:circle)
    xlabel!(L"1/L")
    ylabel!(L"ΔF/T")
    fit = curve_fit(fit_model, 1 ./Llis, λ₀Llis,  [0.5, 0.0])
    x_fit = range(0, stop = 1 ./minimum(Llis), length = 100)
    y_fit = fit_model(x_fit, fit.param)
    plot!(x_fit, y_fit, label = latexstring("\\lambda_0, αΔ=$(round(fit.param[1], digits=3))"), linewidth=2)
    
    scatter!(1 ./Llis, λ₁Llis, yerr = λ₁_errLlis, label =false, markershape=:circle)
    fit2 = curve_fit(fit_model, 1 ./Llis, λ₁Llis,  [0.5, 0.0])
    y_fit2 = fit_model(x_fit, fit2.param)
    plot!(x_fit, y_fit2, label = latexstring("\\lambda_1, αΔ=$(round(fit2.param[1], digits=3))"), linewidth=2)
    return fig
end

fig = plot_scaling_dimension(collect(8:2:12))