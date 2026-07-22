using FibonacciChain
using LinearAlgebra
using JLD

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[7] = log(1 + √2)  # atanh(1/√2) = log(1 + √2)
τlis[end] = 1000.0     # Last value is for γ=1

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

function check_commute_state(
    L::Int,
    τ_idx::Int,
    index::Int;
    )   
    τ = τlis[τ_idx]
    D, inds, avg_range = get_cfg_params_Born(τ_idx, L)
    t = div(D, 2L)
    path = "exm/data/Bulk_measure/monitored_dynamics/L$(L)/gammaind$(τ_idx)/t$(t)_samples$(index).jld"
    data = load(path)
    sample = data["sample"]
    sample = similar(sample)
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    
    l = length(anyon_basis(model))
    M = zeros(l, l)

    for i in 1:l
        st = zeros(length(anyon_basis(model)))
        st[i] = 1.0
        config = MeasureConfig(τ = τ, mode = :sample, t₂ = t*L)
        outcome = bulk_evolution(model, st, config, sample)
        final_state = outcome.state
        M[i, :] = final_state
    end
    
    Y = topological_charge_operator(model)
    
    
    return M, Y
end