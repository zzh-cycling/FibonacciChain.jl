using FibonacciChain
using ITensorMPS, ITensors
using Random
using BenchmarkTools

function get_system_params(τind, L)
    cfg = Dict(
        1 => (1250, 1000, 1000),
        2 => (250, 100, 200),
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
    t, step, start = get(cfg, τind, (2, 2, 1))
    inds = collect(1:step:(t*step))
    avg_range = Int(start*L):2:(Int(t*L)-4)
    # avg_range = 20*L:2:25*L
    return t, inds, avg_range
end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0
L = 32
τind = 10
τ = τlis[τind]
χ = 64
index=1
t, _, _ = get_system_params(τind, L)
rng = MersenneTwister(index)
model = AnyonModel(FibonacciAnyon(), L; pbc = true)
ψ, sites = initial_mps(L)
# config = MeasureConfig(τ=τ, mode=:Born, t₂=t*L, rng=rng, cutoff=1e-12, maxdim=χ)

if length(ARGS) == 0
    println("No arguments provided.")
    println(
        "Usage: julia -p N monitored_dynamics_mps.jl mode L τ_idx index_start index_end",
    )
    println("Example: julia -p 16 monitored_dynamics_mps.jl 2 10 7 1 100")
else
    n = parse(Int64, ARGS[1])
    config = MeasureConfig(
        τ = τ,
        mode = :Born,
        t₂ = t*L,
        rng = rng,
        cutoff = 1e-12,
        maxdim = χ,
        truncate_every_events = n,
    )
    @time mps_mo = bulk_evolution(model, sites, ψ, config)
    # n = 1, 27s, n = 2, 82s, n = 4, 223s, n = 8, 1144s
end
