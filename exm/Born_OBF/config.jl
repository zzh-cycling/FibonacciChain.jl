# Shared configuration for the Born_OBF examples.
# Include this file (also inside @everywhere blocks) instead of redefining
# these settings in each script.

using FibonacciChain

const γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
const τlis = let τ = atanh.(γlis)
    τ[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)
    τ[findfirst(γlis .== 1/√2)] = log(1 + √2)
    τ
end

const λlis = vcat(collect(0.0:0.1:1.5), [atanh(0.999), 11.0])

function obf_model(L::Int, λ::Float64)
    if λ >= 10.0
        return AnyonModel(OBFAnyon(), L; λI = 0.0, pbc = true)
    else
        return AnyonModel(OBFAnyon(), L; λ = λ, pbc = true)
    end
end

function get_born_dynamics_params(ind, L, λ)
    if L >= 18
        if ind == 1
            cfg = Dict(11.0 => (250, 14, 40))
            t, step, start = get(cfg, λ, (150, 14, 20))
        elseif ind == 7
            cfg = Dict(
                atanh(0.999) => (5, 14, 1),
                11.0 => (10, 14, 2),
            )
            t, step, start = get(cfg, λ, (8, 14, 1))
        end
        # Default parameters for MPS
    else
        if ind == 1
            cfg = Dict(11.0 => (400, 14, 50))
            t, step, start = get(cfg, λ, (200, 14, 30))
        elseif ind == 7
            cfg = Dict(11.0 => (18, 14, 2))
            t, step, start = get(cfg, λ, (18, 14, 2))
        end
    end
    inds = collect(1:step:t*L)
    return t, inds, start
end

function get_born_dynamics_params(ind, L, λ)
    if L < 18
        if ind == 1
            cfg = Dict(11.0 => (400, 14, 50))
            t, step, start = get(cfg, λ, (200, 14, 30))
        elseif ind == 7
            cfg = Dict(11.0 => (18, 14, 2))
            t, step, start = get(cfg, λ, (18, 14, 2))
        end
    elseif 18 <= L <=64
        if ind == 1
            cfg = Dict(11.0 => (250, 14, 40))
            t, step, start = get(cfg, λ, (150, 14, 20))
        elseif ind == 7
            cfg = Dict(
                atanh(0.999) => (5, 14, 1),
                11.0 => (10, 14, 2),
            )
            t, step, start = get(cfg, λ, (8, 14, 1))
        end
        # Default parameters for MPS
    elseif L >=128
        if ind == 1
            cfg = Dict(11.0 => (400, 14, 50))
            t, step, start = get(cfg, λ, (200, 14, 30))
        elseif ind == 7
            cfg = Dict(11.0 => (18, 14, 2))
            t, step, start = get(cfg, λ, (4, 14, 1))
        end
    end
    inds = collect(1:step:t*L)
    return t, inds, start
end

function chi_table(ind::Int, L::Int, λ::Float64)
    if λ == atanh(0.999)
        cfg = Dict(
            1 => Dict(
                8 => 35,
            ),
            7 => Dict(
                24 => 200,
            ),
        )
    else
        cfg = Dict(
            1 => Dict(
                8 => 35,
            ),
            7 => Dict(
                24 => 64,
                32 => 96,
                64 => 128,
                128 => 128,
            ),
        )
    end
    return get(get(cfg, ind, Dict()), L, 35)
end
