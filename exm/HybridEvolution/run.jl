using Distributed
using LinearAlgebra
using TOML

include("protocol.jl")

function print_usage()
    println("""
Coherent Fibonacci hybrid MIPT / charge-sharpening scan.
Usage: julia --project=. exm/HybridEvolution/run.jl [options]

  --sizes 8,10,12       Even lengths; comma lists or start:step:stop
  --rates 0:0.1:1       Measurement probabilities
  --trajectories 100    Trajectories per (L,p)
  --seed 1              First seed; consecutive distinct seeds across all jobs
  --shared-seeds        Reuse the same seed range at every (L,p) point
  --time-factor 4       Run ceil(time-factor * L) complete periods
  --periods N           Override scaled time with a fixed positive duration
  --stride 1            Observation interval (always includes 0 and final time)
  --workers 0           Local worker processes; 0 uses serial map
  --backend exact       exact or mps
  --cutoff 1e-12        MPS truncation cutoff
  --mindim 1            Minimum MPS bond dimension
  --maxdim 256          Maximum MPS bond dimension
  --truncate-every 1    MPS truncation interval in events (resets each layer)
  --enforce-constraint  Reproject MPS to legal Fibonacci paths after each layer
  --epsilon 0.05        Distance to either Y eigenvalue for sharpening
  --fraction 0.9        Ensemble fraction defining t_sharp
  --output PATH         New output directory (each trajectory is saved on completion)
  --save-schedule       Store replayable measurement records (enabled by default)
  --help               Show this help

Workers receive the active Julia project and use one BLAS thread.
Existing workers started with julia -p are also supported.
Each completed trajectory is written to L{size}_p{rate}/trajectory_seed{seed}.jld2.
Trajectory files include measurement_mask, outcomes, and unitary_angles for replay.
After a full (L,p) group finishes, ensemble summaries are written in the same directory.
After the scan, averaged.jld2 contains only the point summaries and run configuration.
""")
end

function parse_grid(::Type{T}, value) where T
    if occursin(':', value)
        parts = parse.(T, split(value, ':'))
        length(parts) == 3 || error("Use start:step:stop for ranges")
        start, step, stop = parts
        step > 0 && stop >= start || error("Range must be increasing")
        return collect(start:step:stop)
    end
    return parse.(T, split(value, ','))
end

function options(args)
    defaults = Dict("sizes" => "8,10,12", "rates" => "0:0.1:1",
        "trajectories" => "100", "seed" => "1", "time-factor" => "4",
        "periods" => "0", "stride" => "1", "workers" => "0",
        "epsilon" => "0.05", "fraction" => "0.9", "output" => "",
        "backend" => "exact", "cutoff" => "1e-12", "mindim" => "1",
        "maxdim" => "256", "truncate-every" => "1")
    save_schedule = true
    enforce_constraint = false
    shared_seeds = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--save-schedule"
            save_schedule = true
        elseif arg == "--shared-seeds"
            shared_seeds = true
        elseif arg == "--enforce-constraint"
            enforce_constraint = true
        else
            startswith(arg, "--") || error("Expected --option, got $arg")
            key = arg[3:end]
            haskey(defaults, key) || error("Unknown option: $arg")
            i < length(args) || error("Missing value for $arg")
            i += 1
            defaults[key] = args[i]
        end
        i += 1
    end
    sizes = parse_grid(Int, defaults["sizes"])
    rates = parse_grid(Float64, defaults["rates"])
    trajectories = parse(Int, defaults["trajectories"])
    seed = parse(Int, defaults["seed"])
    factor = parse(Float64, defaults["time-factor"])
    periods = parse(Int, defaults["periods"])
    stride = parse(Int, defaults["stride"])
    workers = parse(Int, defaults["workers"])
    epsilon = parse(Float64, defaults["epsilon"])
    fraction = parse(Float64, defaults["fraction"])
    backend = defaults["backend"]
    backend in ("exact", "mps") || error("backend must be exact or mps")
    cutoff = parse(Float64, defaults["cutoff"])
    mindim = parse(Int, defaults["mindim"])
    maxdim = parse(Int, defaults["maxdim"])
    truncate_every = parse(Int, defaults["truncate-every"])
    isfinite(cutoff) && cutoff >= 0 && 1 <= mindim <= maxdim && truncate_every >= 1 ||
        error("Invalid MPS truncation settings")
    all(L -> L >= 4 && iseven(L), sizes) && !isempty(sizes) || error("Invalid sizes")
    all(p -> isfinite(p) && 0 <= p <= 1, rates) && !isempty(rates) || error("Invalid rates")
    length(unique(sizes)) == length(sizes) && length(unique(rates)) == length(rates) ||
        error("Duplicate sizes/rates")
    trajectories > 0 && seed >= 0 && stride > 0 && workers >= 0 && periods >= 0 ||
        error("Invalid count, seed, stride, workers or periods")
    isfinite(factor) && factor > 0 || error("time-factor must be positive and finite")
    0 < epsilon < sqrt(5.0) / 2 || error("Invalid epsilon")
    0 < fraction <= 1 || error("fraction must be in (0,1]")
    isempty(defaults["output"]) && error("Specify --output PATH")
    output = abspath(defaults["output"])
    ispath(output) && error("Output already exists: $output")
    return (; sizes, rates, trajectories, seed, factor, periods, stride,
        workers, epsilon, fraction, output, save_schedule, shared_seeds,
        backend, cutoff, mindim, maxdim, truncate_every, enforce_constraint)
end

function main(args = ARGS)
    (isempty(args) || "--help" in args) && return print_usage()
    opt = options(args)
    added_workers = Int[]
    try
        if opt.workers > 0
            project = dirname(Base.active_project())
            append!(added_workers, addprocs(opt.workers; exeflags = `--project=$project --threads=1`))
        end
        BLAS.set_num_threads(1)
        for worker in filter(!=(myid()), Distributed.workers())
            remotecall_wait(Core.eval, worker, Main, quote
                using LinearAlgebra
                BLAS.set_num_threads(1)
                include($(joinpath(@__DIR__, "protocol.jl")))
            end)
        end
        mkpath(opt.output)
        metadata = Dict(string(k) => v for (k, v) in pairs(opt))
        metadata["julia_version"] = string(VERSION)
        metadata["initial_state"] = "all-zero coherent fusion path"
        metadata["protocol"] = "HybridConfig Born; projective; random angles; even then odd"
        metadata["seed_assignment"] = opt.shared_seeds ?
            "same consecutive seed range at every (L,p) point" :
            "size-major, rate-next, trajectory-minor, starting at seed"
        open(joinpath(opt.output, "config.toml"), "w") do io
            TOML.print(io, metadata)
        end
        offset = 0
        for L in opt.sizes, p in opt.rates
            duration = opt.periods > 0 ? opt.periods : ceil(Int, opt.factor * L)
            first_seed = opt.seed + (opt.shared_seeds ? 0 : offset)
            seeds = collect(first_seed:(first_seed + opt.trajectories - 1))
            if !opt.shared_seeds
                offset += opt.trajectories
            end
            directory = joinpath(opt.output, "L$(L)_p$(p)")
            mkpath(directory)
            tasks = [(L, p, duration, seed, opt.stride, opt.save_schedule) for seed in seeds]
            processor = process_task_and_save
            if opt.backend == "mps"
                settings = (; cutoff = opt.cutoff, mindim = opt.mindim, maxdim = opt.maxdim,
                    truncate_every_events = opt.truncate_every,
                    enforce_fibonacci_constraint = opt.enforce_constraint)
                tasks = [(task..., settings) for task in tasks]
                processor = process_task_mps_and_save
            end
            jobs = [(task, directory) for task in tasks]
            results = nprocs() > 1 ? pmap(processor, jobs; batch_size = 1) :
                map(processor, jobs)
            save_ensemble(directory, results;
                epsilon = opt.epsilon, fraction = opt.fraction)
            println("Saved L=$L p=$p, $(length(results)) trajectories: $directory")
        end
        println("Saved averaged data: ",
            save_averaged(opt.output, opt.sizes, opt.rates, opt.trajectories))
    finally
        isempty(added_workers) || rmprocs(added_workers)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
