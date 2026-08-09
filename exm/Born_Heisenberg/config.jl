# Shared configuration for the Born_Heisenberg examples.
# Include this file (also inside @everywhere blocks) instead of redefining
# these settings in each script.

using FibonacciChain

const γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
const τlis = let τ = atanh.(γlis)
    τ[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)
    τ[findfirst(γlis .== 1/√2)] = log(1 + √2)
    τ
end

const Δlis = collect(-1.0:0.1:2.0)

heisenberg_model(L::Int, Δ::Float64; J::Float64 = 1.0) =
    AnyonModel(SpinHalf(), L; model_type = :Heisenberg, pbc = true, J = J, Δ = Δ)

function get_born_dynamics_params(ind, L, Δ)
    if L < 18
        if ind == 1
            cfg = Dict(2.0 => (400, 2, 50))
            t, step, start = get(cfg, Δ, (200, 2, 30))
        elseif ind == 7
            cfg = Dict(2.0 => (18, 2, 2))
            t, step, start = get(cfg, Δ, (18, 2, 2))
        end
    elseif L >= 18
        if ind == 1
            cfg = Dict(2.0 => (250, 2, 40))
            t, step, start = get(cfg, Δ, (150, 2, 20))
        elseif ind == 7
            cfg = Dict(2.0 => (10, 2, 2))
            t, step, start = get(cfg, Δ, (8, 2, 1))
        end
    end
    inds = collect(1:step:t*L)
    return t, inds, start
end
