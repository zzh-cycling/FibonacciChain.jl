# Shared configuration for the Born_Ising examples.
# Include this file instead of redefining these settings in each script.

using FibonacciChain

const γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
const τlis = let τ = atanh.(γlis)
    τ[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)
    τ[findfirst(γlis .== 1/√2)] = log(1 + √2)
    τ
end

ising_model(L::Int) = AnyonModel(SpinHalf(), L; model_type = :Ising, pbc = true, measure_operator = :X)

# Shared `(t, inds, avg_range)` table keyed by τ index. Here `t` is the
# number of complete measurement periods; `inds` and `avg_range` are also
# expressed in period indices.
function get_cfg_params_Born(ind, L)
    cfg = Dict(
        7 => (100L, 2, 20L),
    )
    t, step, start = get(cfg, ind, (24L, 10, 5L))
    inds = collect(1:step:t)
    avg_range = start:t-5
    return t, inds, avg_range
end

# Shared seed sanity check
function check_duplicates(seeds)
    if length(seeds) != length(unique(seeds))
        duplicates = findall(x -> count(==(x), seeds) > 1, unique(seeds))
        duplicate_values = unique(seeds)[duplicates]
        println("WARNING: Found duplicate seeds: $duplicate_values")
        return true
    else
        println("No duplicate seeds found in $(length(seeds)) seeds.")
        return false
    end
end
