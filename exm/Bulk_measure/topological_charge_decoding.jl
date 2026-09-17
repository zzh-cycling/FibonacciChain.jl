using FibonacciChain
using JLD2
using LinearAlgebra
using Random

const PHI = (1 + sqrt(5.0)) / 2
const Y_EIGENVALUE = Dict(:y1 => PHI, :ytau => -inv(PHI))
const OTHER_Y_EIGENVALUE = Dict(:y1 => -inv(PHI), :ytau => PHI)
const TRAJECTORY_RE = r"^periods(\d+)_trajectory_seed(\d+)\.jld2$"
const RAW_SCHEMA_VERSION = 3

"""Parse comma-separated integers and integer ranges."""
function parse_int_spec(spec::AbstractString)
    values = Int[]
    for token in split(spec, ',')
        fields = split(strip(token), ':')
        if length(fields) == 1
            push!(values, parse(Int, fields[1]))
        elseif length(fields) == 2
            append!(values, parse(Int, fields[1]):parse(Int, fields[2]))
        elseif length(fields) == 3
            append!(
                values,
                parse(Int, fields[1]):parse(Int, fields[2]):parse(Int, fields[3]),
            )
        else
            error("invalid integer specification: " * token)
        end
    end
    values = sort!(unique(values))
    isempty(values) && error("integer specification must not be empty")
    first(values) >= 1 || error("all requested integers must be positive")
    return values
end

function validate_k_values(k_values::AbstractVector{<:Integer})
    isempty(k_values) && error("at least one K value is required")
    all(>(0), k_values) || error("all K values must be positive")
    issorted(k_values) || error("K values must be sorted")
    length(unique(k_values)) == length(k_values) || error("K values must be unique")
    return Int.(k_values)
end

function validate_record(path::AbstractString)
    data = JLD2.load(path)
    required = ("L", "periods", "tau", "gamma", "sample", "seed")
    for key in required
        haskey(data, key) || error("missing key " * key * " in " * path)
    end
    L = Int(data["L"])
    periods = Int(data["periods"])
    sample = BitMatrix(data["sample"])
    iseven(L) || error("L must be even, got " * string(L) * " in " * path)
    size(sample) == (2periods, L ÷ 2) || error(
        "sample in " * path * " has size " * string(size(sample)) *
        "; expected " * string((2periods, L ÷ 2)),
    )
    return (
        path = abspath(path),
        L = L,
        periods = periods,
        tau = Float64(data["tau"]),
        gamma = Float64(data["gamma"]),
        seed = Int(data["seed"]),
        initial_state = String(get(data, "initial_state", "unknown")),
        initial_state_method = String(get(data, "initial_state_method", "unknown")),
        sample = sample,
    )
end

function validate_record_set(paths::Vector{String})
    isempty(paths) && error("no trajectory paths were supplied")
    records = validate_record.(paths)
    reference = first(records)
    for record in records[2:end]
        for field in (:L, :periods, :tau, :gamma, :initial_state, :initial_state_method)
            getproperty(record, field) == getproperty(reference, field) || error(
                "trajectory metadata mismatch for field " * string(field),
            )
        end
    end
    seeds = getproperty.(records, :seed)
    length(unique(seeds)) == length(seeds) ||
        error("duplicate trajectory seeds: " * string(seeds))
    return records
end

function trajectory_paths(data_dir::AbstractString, requested_seeds::Vector{Int})
    isdir(data_dir) || error("trajectory directory does not exist: " * data_dir)
    index = Dict{Int,String}()
    for filename in readdir(data_dir)
        matched = match(TRAJECTORY_RE, filename)
        isnothing(matched) && continue
        seed = parse(Int, matched.captures[2])
        haskey(index, seed) && error("duplicate seed " * string(seed) * " in " * data_dir)
        index[seed] = joinpath(data_dir, filename)
    end
    missing = filter(seed -> !haskey(index, seed), requested_seeds)
    isempty(missing) || error("missing trajectory seeds in " * data_dir * ": " * string(missing))
    return [index[seed] for seed in requested_seeds]
end

"""
Construct a deterministic orthonormal prefix in each degenerate Y eigenspace.
Only this prefix is needed by the replay routines; no derived decoder
quantities are formed here.
"""
function sector_frames(L::Int, k_max::Int, basis_seed::Int)
    k_max >= 1 || error("k_max must be positive")
    model = AnyonModel(FibonacciAnyon(), L; pbc = true, measure_operator = :Antiferro)
    Y = topological_charge_operator(model)
    decomposition = eigen(Symmetric(Y))
    tolerance = 1e-8

    frames = Dict{Symbol,Matrix{Float64}}()
    dimensions = Dict{Symbol,Int}()
    residuals = Dict{Symbol,Float64}()
    rng = MersenneTwister(basis_seed)

    for sector in (:y1, :ytau)
        target = Y_EIGENVALUE[sector]
        indices = findall(value -> abs(value - target) <= tolerance, decomposition.values)
        dimension = length(indices)
        dimensions[sector] = dimension
        k_max <= dimension || error(
            "requested K_max=" * string(k_max) * ", but L=" * string(L) *
            " " * string(sector) * " has dimension " * string(dimension),
        )
        eigenspace = decomposition.vectors[:, indices]
        frame = if k_max == dimension
            eigenspace
        else
            gaussian = randn(rng, dimension, k_max)
            coefficients = Matrix(qr(gaussian).Q)[:, 1:k_max]
            eigenspace * coefficients
        end
        orthogonality_error = norm(frame' * frame - I(k_max))
        orthogonality_error <= 1e-9 * max(1, k_max) || error(
            string(sector) * " frame is not orthonormal: error=" *
            string(orthogonality_error),
        )
        residual = norm(Y * frame - target * frame) / norm(frame)
        residual <= 1e-9 || error(
            string(sector) * " Y-eigenstate residual is " * string(residual),
        )
        frames[sector] = frame
        residuals[sector] = residual
    end

    projector_dimension_y1 = round(
        Int,
        tr((Y - OTHER_Y_EIGENVALUE[:y1] * I(size(Y, 1))) /
           (Y_EIGENVALUE[:y1] - OTHER_Y_EIGENVALUE[:y1])),
    )
    projector_dimension_y1 == dimensions[:y1] ||
        error("inconsistent y=1 dimension")
    dimensions[:y1] + dimensions[:ytau] == size(Y, 1) ||
        error("sector dimensions do not exhaust the Hilbert space")
    return model, frames, dimensions, residuals
end

"""
Replay one recorded trajectory through the package's exact `bulk_evolution`
`:sample` path. `bulk_evolution` already normalizes after every measurement and
returns the corresponding layer free energies, so `log P` is their negative sum.
The returned `Float32` values are promoted before accumulation.
"""
function replay_log_probability_dynamics(
    model,
    tau::Float64,
    initial_state::AbstractVector,
    sample::BitMatrix;
    enable_tau_eff::Bool = true,
)
    n_layers = FibonacciChain.layers_per_period(model)
    size(sample, 1) % n_layers == 0 || error("sample has incomplete periods")
    n_periods = size(sample, 1) ÷ n_layers
    config = MeasureConfig(
        τ = tau,
        t₂ = n_periods,
        mode = :sample,
        enable_τ_eff = enable_tau_eff,
    )
    outcome = bulk_evolution(model, Vector{Float64}(initial_state), config, sample)
    logp = zeros(Float64, n_periods + 1)
    for period in 1:n_periods
        rows = ((period - 1) * n_layers + 1):(period * n_layers)
        logp[period + 1] = logp[period] - sum(
            value -> Float64(value),
            @view(outcome.free_energys[rows]),
        )
    end
    return logp
end

function replay_frame_dynamics(
    model,
    tau::Float64,
    frame::AbstractMatrix,
    sample::BitMatrix;
    enable_tau_eff::Bool = true,
)
    n_layers = FibonacciChain.layers_per_period(model)
    size(sample, 1) % n_layers == 0 || error("sample has incomplete periods")
    n_periods = size(sample, 1) ÷ n_layers
    logp = Matrix{Float64}(undef, n_periods + 1, size(frame, 2))
    Threads.@threads :dynamic for state_index in axes(frame, 2)
        logp[:, state_index] = replay_log_probability_dynamics(
            model,
            tau,
            @view(frame[:, state_index]),
            sample;
            enable_tau_eff = enable_tau_eff,
        )
    end
    return logp
end

function replay_frame(
    model,
    tau::Float64,
    frame::AbstractMatrix,
    sample::BitMatrix;
    enable_tau_eff::Bool = true,
)
    dynamics = replay_frame_dynamics(
        model,
        tau,
        frame,
        sample;
        enable_tau_eff = enable_tau_eff,
    )
    return copy(@view dynamics[end, :])
end

function validate_truth_sector(truth_sector::Symbol)
    truth_sector in (:y1, :ytau, :unknown) || error(
        "truth_sector must be :y1, :ytau, or :unknown",
    )
    return truth_sector
end

function save_static_logp(
    output_path::AbstractString,
    records,
    logp_state_y1::Array{Float64,2},
    logp_state_ytau::Array{Float64,2},
    k_values::Vector{Int},
    basis_seed::Int,
    truth_sector::Symbol,
    enable_tau_eff::Bool,
)
    first_record = first(records)
    mkpath(dirname(output_path))
    JLD2.jldsave(
        output_path;
        schema_version = RAW_SCHEMA_VERSION,
        description = "State-resolved topological-charge log probabilities",
        observable = "log P(record | phi_y,a)",
        data_layout = "logp_state_y[trajectory,state]",
        L = first_record.L,
        periods = first_record.periods,
        tau = first_record.tau,
        gamma = first_record.gamma,
        enable_tau_eff = enable_tau_eff,
        truth_sector = String(truth_sector),
        basis_seed = basis_seed,
        k_values = k_values,
        k_max = size(logp_state_y1, 2),
        trajectory_seeds = getproperty.(records, :seed),
        logp_state_y1 = logp_state_y1,
        logp_state_ytau = logp_state_ytau,
    )
    return output_path
end

function save_dynamic_logp(
    output_path::AbstractString,
    records,
    logp_state_y1::Array{Float64,3},
    logp_state_ytau::Array{Float64,3},
    k_values::Vector{Int},
    basis_seed::Int,
    truth_sector::Symbol,
    enable_tau_eff::Bool,
    layers_per_period::Int,
)
    first_record = first(records)
    mkpath(dirname(output_path))
    JLD2.jldsave(
        output_path;
        schema_version = RAW_SCHEMA_VERSION,
        description = "State-resolved topological-charge log-probability dynamics",
        observable = "log P(record prefix | phi_y,a)",
        data_layout = "logp_state_y[trajectory,time_period,state]",
        L = first_record.L,
        periods = first_record.periods,
        time_periods = collect(0:first_record.periods),
        layers_per_period = layers_per_period,
        tau = first_record.tau,
        gamma = first_record.gamma,
        enable_tau_eff = enable_tau_eff,
        truth_sector = String(truth_sector),
        basis_seed = basis_seed,
        k_values = k_values,
        k_max = size(logp_state_y1, 3),
        trajectory_seeds = getproperty.(records, :seed),
        logp_state_y1 = logp_state_y1,
        logp_state_ytau = logp_state_ytau,
    )
    return output_path
end

function decode_records(
    paths::Vector{String},
    output_path::AbstractString,
    k_values::Vector{Int};
    basis_seed::Int = 314159,
    truth_sector::Symbol = :y1,
    enable_tau_eff::Bool = true,
)
    isfile(output_path) && return (status = :skipped, output_path = output_path)
    k_values = validate_k_values(k_values)
    truth_sector = validate_truth_sector(truth_sector)
    records = validate_record_set(paths)
    L = first(records).L
    k_max = maximum(k_values)

    @info "Dense diagonalization of Y" L k_max basis_seed
    model, frames, _, _ = sector_frames(L, k_max, basis_seed)
    logp_state_y1 = Matrix{Float64}(undef, length(records), k_max)
    logp_state_ytau = similar(logp_state_y1)

    for (record_index, record) in enumerate(records)
        @info "Replaying record" record_index n_records = length(records) seed = record.seed
        logp_state_y1[record_index, :] = replay_frame(
            model,
            record.tau,
            frames[:y1],
            record.sample;
            enable_tau_eff = enable_tau_eff,
        )
        logp_state_ytau[record_index, :] = replay_frame(
            model,
            record.tau,
            frames[:ytau],
            record.sample;
            enable_tau_eff = enable_tau_eff,
        )
    end

    save_static_logp(
        output_path,
        records,
        logp_state_y1,
        logp_state_ytau,
        k_values,
        basis_seed,
        truth_sector,
        enable_tau_eff,
    )
    println("done: output=" * output_path * " records=" * string(length(records)) *
            " K_max=" * string(k_max))
    return (status = :ok, output_path = output_path)
end

function decode_records_dynamics(
    paths::Vector{String},
    output_path::AbstractString,
    k_values::Vector{Int};
    basis_seed::Int = 314159,
    truth_sector::Symbol = :y1,
    selected_seeds::Vector{Int} = Int[],
    enable_tau_eff::Bool = true,
)
    isfile(output_path) && return (status = :skipped, output_path = output_path)
    k_values = validate_k_values(k_values)
    truth_sector = validate_truth_sector(truth_sector)
    records = validate_record_set(paths)
    L = first(records).L
    periods = first(records).periods
    k_max = maximum(k_values)

    @info "Dense diagonalization of Y for dynamics" L k_max basis_seed
    model, frames, _, _ = sector_frames(L, k_max, basis_seed)
    layers_per_period = FibonacciChain.layers_per_period(model)
    logp_state_y1 = Array{Float64}(undef, length(records), periods + 1, k_max)
    logp_state_ytau = similar(logp_state_y1)

    for (record_index, record) in enumerate(records)
        @info "Replaying dynamics" record_index n_records = length(records) seed = record.seed
        current_y1 = replay_frame_dynamics(
            model,
            record.tau,
            frames[:y1],
            record.sample;
            enable_tau_eff = enable_tau_eff,
        )
        current_ytau = replay_frame_dynamics(
            model,
            record.tau,
            frames[:ytau],
            record.sample;
            enable_tau_eff = enable_tau_eff,
        )
        @views logp_state_y1[record_index, :, :] .= current_y1
        @views logp_state_ytau[record_index, :, :] .= current_ytau
    end

    isempty(selected_seeds) || all(
        seed -> seed in getproperty.(records, :seed),
        selected_seeds,
    ) || error("selected_seeds contains a seed absent from the record set")
    save_dynamic_logp(
        output_path,
        records,
        logp_state_y1,
        logp_state_ytau,
        k_values,
        basis_seed,
        truth_sector,
        enable_tau_eff,
        layers_per_period,
    )
    println("done: output=" * output_path * " records=" * string(length(records)) *
            " periods=" * string(periods) * " K_max=" * string(k_max))
    return (status = :ok, output_path = output_path)
end

function usage()
    println("Usage:")
    println("  julia --project=. --threads=N exm/Bulk_measure/topological_charge_decoding.jl single TRAJECTORY OUTPUT K_SPEC [BASIS_SEED=314159] [TRUTH=y1]")
    println("  julia --project=. --threads=N exm/Bulk_measure/topological_charge_decoding.jl ensemble DATA_DIR OUTPUT SEED_SPEC K_SPEC [BASIS_SEED=314159] [TRUTH=y1]")
    println("  julia --project=. --threads=N exm/Bulk_measure/topological_charge_decoding.jl dynamics DATA_DIR OUTPUT SEED_SPEC K_SPEC [BASIS_SEED=314159] [TRUTH=y1]")
    println("  julia --project=. --threads=N exm/Bulk_measure/topological_charge_decoding.jl scan DATA_ROOT OUTPUT_DIR L_SPEC GAMMA_DIR SEED_SPEC K_SPEC [BASIS_SEED=314159] [TRUTH=y1]")
end

function main(args::Vector{String})
    isempty(args) && return usage()
    action = Symbol(lowercase(args[1]))
    if action == :single
        length(args) in 4:6 || error("single expects 3 to 5 arguments after the action")
        input_path, output_path = args[2], args[3]
        k_values = parse_int_spec(args[4])
        basis_seed = length(args) >= 5 ? parse(Int, args[5]) : 314159
        truth_sector = length(args) >= 6 ? Symbol(lowercase(args[6])) : :y1
        return decode_records(
            [input_path],
            output_path,
            k_values;
            basis_seed = basis_seed,
            truth_sector = truth_sector,
        )
    elseif action == :ensemble
        length(args) in 5:7 || error("ensemble expects 4 to 6 arguments after the action")
        data_dir, output_path = args[2], args[3]
        seeds = parse_int_spec(args[4])
        k_values = parse_int_spec(args[5])
        basis_seed = length(args) >= 6 ? parse(Int, args[6]) : 314159
        truth_sector = length(args) >= 7 ? Symbol(lowercase(args[7])) : :y1
        return decode_records(
            trajectory_paths(data_dir, seeds),
            output_path,
            k_values;
            basis_seed = basis_seed,
            truth_sector = truth_sector,
        )
    elseif action == :dynamics
        length(args) in 5:7 || error("dynamics expects 4 to 6 arguments after the action")
        data_dir, output_path = args[2], args[3]
        seeds = parse_int_spec(args[4])
        k_values = parse_int_spec(args[5])
        basis_seed = length(args) >= 6 ? parse(Int, args[6]) : 314159
        truth_sector = length(args) >= 7 ? Symbol(lowercase(args[7])) : :y1
        return decode_records_dynamics(
            trajectory_paths(data_dir, seeds),
            output_path,
            k_values;
            basis_seed = basis_seed,
            truth_sector = truth_sector,
        )
    elseif action == :scan
        length(args) in 7:9 || error("scan expects 6 to 8 arguments after the action")
        data_root, output_dir = args[2], args[3]
        sizes = parse_int_spec(args[4])
        gamma_dir = args[5]
        seeds = parse_int_spec(args[6])
        k_values = parse_int_spec(args[7])
        basis_seed = length(args) >= 8 ? parse(Int, args[8]) : 314159
        truth_sector = length(args) >= 9 ? Symbol(lowercase(args[9])) : :y1
        for L in sizes
            data_dir = joinpath(data_root, "L" * string(L), gamma_dir)
            output_path = joinpath(
                output_dir,
                "topological_charge_decoding_L" * string(L) * "_" * gamma_dir *
                "_seeds" * string(first(seeds)) * "-" * string(last(seeds)) *
                "_bseed" * string(basis_seed) * ".jld2",
            )
            result = decode_records(
                trajectory_paths(data_dir, seeds),
                output_path,
                k_values;
                basis_seed = basis_seed,
                truth_sector = truth_sector,
            )
            println("L=" * string(L) * " status=" * string(result.status) *
                    " output=" * string(result.output_path))
        end
        return nothing
    end
    error("unknown action " * string(action))
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
