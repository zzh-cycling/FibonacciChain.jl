using FibonacciChain
using LinearAlgebra
using KrylovKit
using JLD, JLD2


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

function transfer_matrix(L::Int64, τ_idx::Int64, index::Int64)
        try
            τ = τlis[τ_idx]
            data = load("exm/data/Bulk_measure/monitored_dynamics/L12/gammaind10/t4_samples1.jld")
            sample = data["sample"]
        
            D, _, _ = get_cfg_params_Born(τ_idx, L)
            t = div(D, 2)
            model = AnyonModel(FibonacciAnyon(), L; pbc = true)
            config = MeasureConfig(τ = τ, mode = :sample, t₂ = t)
            
            model_basis = anyon_basis(model)
            l = length(model_basis)

            transfer_matrix = zeros(l, l)
            for i in 1:l 
                st = zeros(l)
                st[i] = 1.0
    
                outcome = bulk_evolution(model, st, config, sample)
                final_st = outcome.state
                transfer_matrix[:, i] = final_st
            end


            
            return transfer_matrix
            # return (L, τ_idx, index, :success, nothing)
        catch e
            return (L, τ_idx, index, :failed, e)
        end
end

TF = transfer_matrix(12, 10, 1)
energy, states = eigen(TF)