# Shared configuration for the Bulk_measure examples.
# Include this file (also inside @everywhere blocks) instead of redefining
# these settings in each script.

using FibonacciChain

const γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
const τlis = let τ = atanh.(γlis)
    τ[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)
    τ[findfirst(γlis .== 1/√2)] = log(1 + √2)
    τ
end
# Extended γ/τ grid (finer around γ ≈ 0.8); superset of config.jl's γlis/τlis
γlis_ext = vcat([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1], collect(0.77:0.01:0.79), collect(0.81:0.01:0.82), [0.825], collect(0.83:0.01:0.84))
τlis_ext = atanh.(γlis_ext)
τlis_ext[7] = log(1 + √2)
τlis_ext[12] = 1000.0

fib_model(L::Int) = AnyonModel(FibonacciAnyon(), L; pbc = true)

ising_model(L::Int) = AnyonModel(SpinHalf(), L; model_type = :Ising, pbc = true, measure_operator = :X)

# Shared δt lists for the spatial-temporal correlation runs
# (identical in parrellel.jl and parrellel_collect.jl).
function get_δtL_Born(τ, L)
    if L == 6
        table = Dict(
            atanh(0.1) => vcat(collect(1:150), collect(152:2:620)),
            atanh(0.2) => (collect(50:2:120)),
            atanh(0.3) => (collect(1:45)),
            atanh(0.4) => (collect(1:35)),
            atanh(0.5) => (collect(1:25)),
            atanh(0.6) => (collect(1:10)),
        )
        δtlis = get(table, τ, collect(1:10))
    elseif L == 8
        table = Dict(
            atanh(0.1) => collect(50:25:500),
            atanh(0.2) => (collect(65:2:145)),
            atanh(0.3) => (collect(1:54)),
            atanh(0.4) => (collect(1:45)),
            atanh(0.5) => (collect(1:30)),
            atanh(0.6) => (collect(1:16)),
            atanh(1/√2) => (collect(1:12)),
        )
        δtlis = get(table, τ, collect(1:10))
    elseif L == 10
        table = Dict(
            atanh(0.1) => collect(100:25:500),
            atanh(0.2) => (collect(80:2:130)),
            atanh(0.3) => (collect(30:60)),
            atanh(0.4) => (collect(1:35)),
            atanh(0.5) => (collect(1:24)),
            atanh(0.6) => (collect(1:22)),
            atanh(1/√2) => (collect(1:16)),
        )
        δtlis = get(table, τ, collect(1:10))
    elseif L == 12
        table = Dict(
            atanh(0.1) => collect(300:25:600),
            atanh(0.2) => (collect(100:2:136)),
            atanh(0.3) => vcat([10, 20, 30], collect(40:2:64), collect(70:10:90)),
            atanh(0.4) => vcat([4, 8, 12], collect(15:35), [40, 45]),
            atanh(0.5) => (collect(1:20)),
            atanh(0.6) => (collect(1:25)),
            atanh(1/√2) => (collect(1:18)),
        )
        δtlis = get(table, τ, collect(1:15))
    elseif L == 14
        table = Dict(
            atanh(0.1) => (collect(50:25:550)),
            atanh(0.2) => (collect(120:2:146)),
            atanh(0.3) => vcat([10, 20, 30], collect(40:2:70), collect(80:10:90)),
            atanh(0.4) => vcat([5, 10, 15, 20], collect(25:39), [45, 50]),
            atanh(0.5) => (collect(5:24)),
            atanh(0.6) => (collect(1:28)),
            atanh(1/√2) => (collect(1:8)),
            atanh(0.8) => (collect(1:15)),
        )
        δtlis = get(table, τ, collect(1:8))
    elseif L == 16
        table = Dict(
            atanh(0.1) => sort(vcat(collect(80:100:780), [640, 740])),
            atanh(0.2) => vcat([100, 110, 120], collect(135:2:160), [170, 180, 190, 200]),
            atanh(0.3) => sort(vcat(collect(10:10:50), collect(55:4:65), collect(67:2:75), collect(80:10:90))),
            atanh(0.4) => vcat([5, 10, 15, 20, 25, 30], (collect(32:42)), [45, 50]),
            atanh(0.5) => (collect(4:22)),
            atanh(0.6) => vcat(collect(1:16), [20, 25]),
            atanh(1/√2) => (collect(1:8)),
            atanh(0.8) => (collect(1:16)),
        )
        δtlis = get(table, τ, collect(1:8))
    elseif L == 18
        table = Dict(
            atanh(0.1) => vcat(collect(500:50:650), [660, 670, 730, 800, 850]),
            atanh(0.2) => sort(vcat([130, 140, 158, 189, 190], collect(150:5:180))),
            atanh(0.3) => sort(vcat([10, 20, 30, 40, 50], collect(71:74), collect(60:5:90))),
            atanh(0.4) => sort(vcat([5, 15], [38, 39, 41, 42], collect(25:5:50))),
            atanh(0.5) => vcat([5, 10], collect(20:28)),
            atanh(0.6) => vcat([2], collect(5:15), [20, 25]),
            atanh(1/√2) => (collect(1:10)),
            atanh(0.8) => (collect(1:10)),
            atanh(0.9) => (collect(1:10)),
        )
        δtlis = get(table, τ, collect(1:8))
    elseif L == 20
        table = Dict(
            atanh(0.1) => collect(650:10:750),
            atanh(0.2) => vcat([150, 160, 170], collect(173:2:193), [220, 210, 220]),
            atanh(0.3) => vcat([10, 30, 50], collect(68:4:88), [95, 100]),
            atanh(0.4) => vcat([10, 20], collect(38:48), [60]),
            atanh(0.5) => vcat([5, 15, 20], collect(24:34), [40]),
            atanh(0.6) => vcat([4,8], collect(12:20)),
            atanh(1/√2) => (collect(1:12)),
        )
        δtlis = get(table, τ, collect(1:10))
    elseif L == 22
        table = Dict(
            atanh(0.4) =>vcat([10, 20, 30, 35, 37], collect(41:2:49), collect(50:2:52), [60, 65, 80]),
            atanh(0.5) => sort(vcat([5, 10], collect(16:2:36), [40])),
            atanh(0.6) => sort(vcat(collect(5:2:19), [20])),
            atanh(1/√2) => (collect(1:15)),
        )
        δtlis = get(table, τ, collect(1:10))
    elseif L == 24
        table = Dict(
            atanh(0.5) => sort(vcat(collect(25:2:35), collect(24:2:28), [5, 10, 15, 20])),
            atanh(0.6) => sort(vcat(collect(10:2:20), [5, 25])),
            atanh(1/√2) => (collect(1:15)),
        )
        δtlis = get(table, τ, collect(1:10))
    elseif L == 26
        table = Dict(
            atanh(1/√2) => (collect(1:15)),
            atanh(0.999) => (collect(1:8)),
            1000.0 => (collect(1:8)),
        )
        δtlis = get(table, τ, collect(1:10))
    elseif L == 28
        table = Dict(
            atanh(0.999) => (collect(1:8)),
            1000.0 => (collect(1:8)),
        )
        δtlis = get(table, τ, collect(1:10))
    else
        δtlis = collect(1:10)
    end
    return δtlis
end

# Shared (D, inds, avg_range) table keyed by τ
# (identical in parrellel.jl, parrellel_collect.jl, parrellel_dynamics.jl and ps_prob.jl).
function get_system_params(τ, L)
    cfg = Dict(
        atanh(0.1) => (2500L, 1000, 1500L),
        atanh(0.2) => (500L, 100, 250L),
        atanh(0.3) => (120L, 48, 100L),
        atanh(0.4) => (100L, 40, 80L),
        atanh(0.5) => (80L, 32, 40L),
        atanh(0.6) => (45L, 20, 30L),
        log(1 + √2) => (35L, 14, 20L),
        atanh(0.8) => (25L, 10, 10L),
        atanh(0.9) => (8L, 4, 4L),
        atanh(0.95) => (8L, 4, 4L),
        atanh(0.999) => (5L, 2, 2L),
    )
    D, step, start = get(cfg, τ, (5L, 2, 2L))
    inds = collect(1:step:D)
    avg_range = start:(div(D, 2)-5)
    return D, inds, avg_range
end

# Shared (D, inds, avg_range) table keyed by τ index
# (identical in defect_scaling_dimension.jl and transfer_matrix.jl).
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

# Shared MPS evolution parameters keyed by τ index
# (identical in defect_scaling_dimension.jl, transfer_matrix.jl and measure_correlation.jl).
function get_mps_params_Born(τind, L)
    cfg = if L <= 32
        Dict(
            1 => (750, 1000, 300),
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
            12 => (2, 2, 1),
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
            12 => (2, 2, 1)
        )
    end
    t, step, start = get(cfg, τind, (8, 2, 2))
    inds = collect(1:step:(t*L))
    avg_range = Int(start*L):2:(Int(t*L)-4)
    return t, inds, avg_range
end

# Shared seed sanity check (identical in monitored_dynamics.jl
# and monitored_dynamics_mps.jl).
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

# Shared default bond dimensions
# (identical in measure_correlation.jl, defect_scaling_dimension.jl
# and monitored_dynamics_mps.jl).
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
