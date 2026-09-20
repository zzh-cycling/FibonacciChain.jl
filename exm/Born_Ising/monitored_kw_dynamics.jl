# Run with the package available in the active Julia environment, plus JLD2.
# Workers can be supplied with `julia -p N`; no workers are started implicitly.
using Distributed
using JLD2

@everywhere begin
    using FibonacciChain
    using LinearAlgebra
    using Random
    using JLD2

    function kw_samples_generate(task)
        L, tau, periods, seed, output_dir, half_boundary = task
        path = joinpath(output_dir, "t$(periods)_samples$(seed).jld2")
        ispath(path) && error("Refusing to overwrite $path")
        model = AnyonModel(SpinHalf(), L; model_type = :Ising, pbc = true)

        # Edit the chain/reference preparation here if needed, together with
        # initial_state below. The reference is the leading qubit.
        chain = zeros(Float64, 2^L)
        chain[1] = 1.0
        state = ising_reference_state(chain) # |+>_reference ⊗ |0...0>_chain
        initial_state = "reference_plus_chain_zero"
        config = MeasureConfig(τ = tau, t₂ = periods, mode = :Born,
            rng = MersenneTwister(seed), enable_τ_eff = half_boundary)
        result = ising_reference_evolution(model, state, config)

        # Small, explicit schema. KW includes t=0; sample has 2periods rows.
        # No final state, posterior, sector assignments, or provenance paths.
        JLD2.jldsave(path;
            sample = result.samples,
            seed,
            kw_expectation_values = result.kw_expectation_values,
            L, tau, periods, initial_state, half_boundary)
        return path
    end
end

"""Generate independent trajectories using pmap, as in monitored_dynamics.jl."""
function kw_samples_generate(L::Int, tau::Real, periods::Int, seeds, output_dir;
    half_boundary::Bool = false)
    L >= 2 && periods >= 1 || throw(ArgumentError("require L >= 2 and periods >= 1"))
    isfinite(tau) && tau >= 0 || throw(ArgumentError("tau must be finite and nonnegative"))
    seed_list = collect(Int, seeds)
    !isempty(seed_list) && all(>=(0), seed_list) && length(unique(seed_list)) == length(seed_list) ||
        throw(ArgumentError("seeds must be nonempty, nonnegative, and unique"))
    directory = abspath(output_dir)
    mkpath(directory)
    paths = [joinpath(directory, "t$(periods)_samples$(seed).jld2") for seed in seed_list]
    any(ispath, paths) && error("A requested trajectory already exists in $directory")
    # Prevent mixing parameters if this directory already contains trajectories.
    for file in readdir(directory; join = true)
        occursin(r"^t\d+_samples\d+\.jld2$", basename(file)) || continue
        data = JLD2.load(file)
        (data["L"], data["tau"], data["periods"], data["half_boundary"], data["initial_state"]) ==
            (L, Float64(tau), periods, half_boundary, "reference_plus_chain_zero") ||
            error("Incompatible trajectory parameters in $file; use another directory")
    end
    tasks = [(L, Float64(tau), periods, seed, directory, half_boundary) for seed in seed_list]
    return pmap(kw_samples_generate, tasks)
end

"""Collect raw KW trajectories (rows = seeds, columns = t=0:periods)."""
function kw_samples_collect(output_dir, output_file)
    ispath(output_file) && error("Refusing to overwrite $output_file")
    files = filter(f -> occursin(r"^t\d+_samples\d+\.jld2$", basename(f)),
        readdir(output_dir; join = true))
    isempty(files) && error("No KW trajectories in $output_dir")
    rows = [JLD2.load(file) for file in files]
    sort!(rows; by = row -> row["seed"])
    parameters = ("L", "tau", "periods", "initial_state", "half_boundary")
    first_row = first(rows)
    all(row -> all(key -> row[key] == first_row[key], parameters), rows) ||
        error("Cannot combine different trajectory parameters")
    seeds = Int[row["seed"] for row in rows]
    length(unique(seeds)) == length(seeds) || error("Duplicate trajectory seeds")
    periods = first_row["periods"]
    all(row -> length(row["kw_expectation_values"]) == periods + 1, rows) ||
        error("KW time-series length mismatch")
    kw = permutedims(hcat([row["kw_expectation_values"] for row in rows]...))
    mkpath(dirname(abspath(output_file)))
    JLD2.jldsave(output_file; seed = seeds, kw_expectation_values = kw,
        L = first_row["L"], tau = first_row["tau"], periods,
        initial_state = first_row["initial_state"], half_boundary = first_row["half_boundary"])
    return output_file
end

if abspath(PROGRAM_FILE) == @__FILE__
    usage = """
    From the checkout root:
    Generate: julia --project=. [-p N] exm/Born_Ising/monitored_kw_dynamics.jl 1 L TAU PERIODS FIRST_SEED COUNT [OUTPUT_DIR] [HALF_BOUNDARY]
    Collect:  julia --project=. exm/Born_Ising/monitored_kw_dynamics.jl 2 INPUT_DIR OUTPUT_FILE
    HALF_BOUNDARY defaults to false (all layers have full strength).
    TAU is the measurement strength, not gamma; gamma = tanh(TAU).
    """
    isempty(ARGS) && error(usage)
    if ARGS[1] == "1"
        length(ARGS) in 6:8 || error(usage)
        L = parse(Int, ARGS[2])
        tau = parse(Float64, ARGS[3])
        periods = parse(Int, ARGS[4])
        first_seed, count = parse.(Int, ARGS[5:6])
        first_seed >= 0 && count > 0 || error("require nonnegative FIRST_SEED and positive COUNT")
        output_dir = length(ARGS) >= 7 ? ARGS[7] : joinpath(@__DIR__, "..", "data",
            "Born_Ising", "categorical_kw", "L$(L)", "tau$(tau)")
        half_boundary = length(ARGS) == 8 ? parse(Bool, ARGS[8]) : false
        paths = kw_samples_generate(L, tau, periods,
            first_seed:(first_seed + count - 1), output_dir; half_boundary)
        println("Saved $(length(paths)) KW trajectories to $(abspath(output_dir))")
    elseif ARGS[1] == "2"
        length(ARGS) == 3 || error(usage)
        println("Saved KW dynamics to ", kw_samples_collect(ARGS[2], ARGS[3]))
    else
        error(usage)
    end
end
