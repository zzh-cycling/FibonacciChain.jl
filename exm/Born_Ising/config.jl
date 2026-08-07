# Shared configuration for the Born_Ising examples.
# Include this file instead of redefining these settings in each script.

using FibonacciChain

const γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
const τlis = let τ = atanh.(γlis)
    τ[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)
    τ[findfirst(γlis .== 1/√2)] = log(1 + √2)
    τ
end

ising_model(L::Int) = AnyonModel(IsingAnyon(), L; pbc = true, measure_operator = :X)

# Shared (D, inds, avg_range) table keyed by τ index
function get_cfg_params_Born(ind, L)
    cfg = Dict(
        1 => (2500L, 1000, 750L),
        2 => (500L,  100, 120L),
        3 => (200L,  40, 80L),
        4 => (100L,  40, 40L),
        5 => (80L,   32, 20L),
        6 => (45L,   20, 15L),
        7 => (35L,   14, 10L),
        8 => (25L,   10, 5L),
        9 => (8L,    4, 2L),
        10 => (8L,    4, 2L),
        11 => (5L,    2, 1L),
        12 => (5L, 2, L),
    )
    D, step, start = get(cfg, ind, (24L, 10, 5L))
    inds = collect(1:step:div(D,2))
    avg_range = start:2:div(D,2)-5
    return D, inds, avg_range
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
