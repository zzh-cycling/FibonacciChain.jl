# Shared configuration for the Born_Heisenberg examples.
# Include this file (also inside @everywhere blocks) instead of redefining
# these settings in each script.

using FibonacciChain

const γlis = [0.01, 0.1, 0.2, 1/√2, 1]
const τlis = let τ = atanh.(γlis)
    τ[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)
    τ[findfirst(γlis .== 1/√2)] = log(1 + √2)
    τ
end

const Δlis = unique!(sort(vcat(collect(-1.0:0.2:1.0), collect(-1.04:0.02:-0.96), collect(0.96:0.02:1.04))))

heisenberg_model(L::Int, Δ::Float64; J::Float64 = 1.0) =
    AnyonModel(SpinHalf(), L; model_type = :Heisenberg, pbc = true, J = J, Δ = Δ)

function get_born_dynamics_params(ind, L, Δ)
    if L < 18
        if ind == 1
            if Δ < -0.0
                cfg = Dict(-1.0 => (80, 2, 20))
                t, step, start = get(cfg, Δ, (80, 2, 20))
            else
                cfg = Dict(2.0 => (40, 2, 10))
                t, step, start = get(cfg, Δ, (30, 2, 8))
            end 
        elseif ind == 2
            # γ = 0.1: weak measurement
            if Δ < -0.0
                cfg = Dict(-1.0 => (10, 2, 2))
                t, step, start = get(cfg, Δ, (10, 2, 2))
            else
                cfg = Dict(2.0 => (5, 2, 1)) 
                t, step, start = get(cfg, Δ, (5, 2, 1))
            end
        else
            cfg = Dict(2.0 => (5, 2, 1)) 
            t, step, start = get(cfg, Δ, (5, 2, 1))
        end
    elseif L >= 18
        if ind == 1
            cfg = Dict(2.0 => (250, 2, 40))
            t, step, start = get(cfg, Δ, (150, 2, 20))
        elseif ind == 7
            cfg = Dict(2.0 => (10, 2, 2))
            t, step, start = get(cfg, Δ, (8, 2, 1))
        elseif ind == 10
            cfg = Dict(0.9 => (6, 2, 3))
            t, step, start = get(cfg, Δ, (4, 2, 2))
        end
    end
    inds = collect(1:step:t*L)
    return t, inds, start
end
