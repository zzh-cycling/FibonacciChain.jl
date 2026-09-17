using Distributed

const BORN_OBF_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BORN_OBF_CONFIG)

@everywhere begin
    using JLD2
    using Statistics

    const OBF_LAYER_LABELS = [
        "XZZ_1",
        "ZZX_1",
        "XZZ_2",
        "ZZX_2",
        "XZZ_3",
        "ZZX_3",
        "X",
        "ZZX_3",
        "XZZ_3",
        "ZZX_2",
        "XZZ_2",
        "ZZX_1",
        "XZZ_1",
        "ZZ",
    ]

    const OBF_LAYERS_PER_PERIOD = 14

    """
        obf_active_sites(L, layer)

    Return the physical sites whose outcomes are stored in an OBF measurement
    layer. Other columns in the trajectory `BitMatrix` are zero padding and must
    not be included in probability estimates.
    """
    function obf_active_sites(L::Int, layer::Int)
        1 <= layer <= OBF_LAYERS_PER_PERIOD ||
            throw(ArgumentError("OBF layer must be in 1:14, got $layer"))
        model = obf_model(L, 0.0)
        return FibonacciChain._get_sample_column_indices(model, layer)
    end

    """
        analyze_obf_trajectory(sample; start_period=1, stop_period=nothing)

    Compute layer-resolved measurement-outcome statistics for one OBF Born
    trajectory. For layer `l`, the returned histogram estimates

        P_l(K = k),   K = sum(s_i for i in A_l),

    where `A_l` contains only the active sites of that layer. The site-resolved
    array estimates `P_l(s_i = 1)` and is `NaN` on inactive, padded columns.
    The connected spatial correlator uses cyclic separation within the ordered
    active-site list,

        C_l(d) = mean(sigma_a sigma_{a+d}) - mean(sigma)^2,

    with `sigma = 2s - 1`. For the X and ZZ layers this is the physical lattice
    separation; for XZZ/ZZX it is separation on the corresponding three-site
    measurement sublattice.
    """
    function analyze_obf_trajectory(
        sample::AbstractMatrix;
        start_period::Int = 1,
        stop_period::Union{Int,Nothing} = nothing,
    )
        n_rows, L = size(sample)
        n_rows % OBF_LAYERS_PER_PERIOD == 0 || error(
            "OBF sample must have a multiple of 14 rows, got size $(size(sample))",
        )
        all(x -> x == 0 || x == 1, sample) ||
            throw(ArgumentError("measurement outcomes must be binary (0 or 1)"))

        n_periods = n_rows ÷ OBF_LAYERS_PER_PERIOD
        last_period = isnothing(stop_period) ? n_periods : stop_period
        1 <= start_period <= last_period <= n_periods || throw(
            ArgumentError(
                "period window must satisfy 1 ≤ start_period ≤ stop_period ≤ $n_periods; " *
                "got $start_period:$last_period",
            ),
        )

        period_range = start_period:last_period
        n_stationary_periods = length(period_range)
        probability_profile = fill(NaN, OBF_LAYERS_PER_PERIOD, L)
        count_distribution = zeros(Float64, OBF_LAYERS_PER_PERIOD, L + 1)
        mean_probability = zeros(Float64, OBF_LAYERS_PER_PERIOD)
        mean_magnetization = zeros(Float64, OBF_LAYERS_PER_PERIOD)
        connected_spatial_correlation =
            fill(NaN, OBF_LAYERS_PER_PERIOD, L ÷ 2)
        active_site_mask = falses(OBF_LAYERS_PER_PERIOD, L)
        active_site_count = zeros(Int, OBF_LAYERS_PER_PERIOD)

        for layer = 1:OBF_LAYERS_PER_PERIOD
            sites = obf_active_sites(L, layer)
            active_site_mask[layer, sites] .= true
            active_site_count[layer] = length(sites)

            first_row = (first(period_range) - 1) * OBF_LAYERS_PER_PERIOD + layer
            last_row = (last(period_range) - 1) * OBF_LAYERS_PER_PERIOD + layer
            rows = first_row:OBF_LAYERS_PER_PERIOD:last_row
            outcomes = Float64.(sample[rows, sites])

            site_probabilities = vec(mean(outcomes; dims = 1))
            probability_profile[layer, sites] = site_probabilities
            mean_probability[layer] = mean(outcomes)
            mean_magnetization[layer] = 2 * mean_probability[layer] - 1

            spins = 2 .* outcomes .- 1
            for distance = 1:(length(sites) ÷ 2)
                pair_mean = mean(
                    mean(
                        spins[:, site_index] .*
                        spins[:, mod1(site_index + distance, length(sites))],
                    ) for
                    site_index = eachindex(sites)
                )
                connected_spatial_correlation[layer, distance] =
                    pair_mean - mean_magnetization[layer]^2
            end

            for k in vec(sum(outcomes; dims = 2))
                count_distribution[layer, Int(k) + 1] += 1
            end
            count_distribution[layer, :] ./= n_stationary_periods
        end

        return (
            probability_profile = probability_profile,
            count_distribution = count_distribution,
            mean_probability = mean_probability,
            mean_magnetization = mean_magnetization,
            connected_spatial_correlation = connected_spatial_correlation,
            active_site_mask = active_site_mask,
            active_site_count = active_site_count,
            start_period = start_period,
            stop_period = last_period,
            n_stationary_periods = n_stationary_periods,
        )
    end

    function obf_sample_directory(
        backend::Symbol,
        L::Int,
        τ_idx::Int,
        λ::Float64;
        χ::Union{Int,Nothing} = nothing,
    )
        if backend == :exact
            τ = τlis[τ_idx]
            return "exm/data/OBF/Born_dynamics_records/L$(L)/τ$(τ)/λ$(λ)"
        elseif backend == :mps
            isnothing(χ) && throw(ArgumentError("χ is required for the MPS backend"))
            base_directory =
                "exm/data/OBF/Born_dynamics_records_mps/L$(L)/gammaind$(τ_idx)/λ$(λ)"
            current_directory = joinpath(base_directory, "chi$(χ)")
            if !isdir(current_directory) && isdir(base_directory)
                @warn "Using legacy MPS directory layout without a chi subdirectory" directory =
                    base_directory
                return base_directory
            end
            return current_directory
        else
            throw(ArgumentError("backend must be :exact or :mps, got $backend"))
        end
    end

    function _sample_index(filename::AbstractString)
        matched = match(r"_samples(\d+)", filename)
        return isnothing(matched) ? typemax(Int) : parse(Int, matched.captures[1])
    end

    function _sample_time(filename::AbstractString)
        matched = match(r"^t(\d+)_samples", filename)
        return isnothing(matched) ? nothing : parse(Int, matched.captures[1])
    end

    function obf_sample_files(
        backend::Symbol,
        L::Int,
        τ_idx::Int,
        λ::Float64;
        χ::Union{Int,Nothing} = nothing,
        sample_num::Int = 0,
    )
        t, _, _ = get_born_dynamics_params(τ_idx, L, λ)
        directory = obf_sample_directory(backend, L, τ_idx, λ; χ = χ)
        isdir(directory) || error("Sample directory not found: $directory")

        prefix = "t$(t)_samples"
        suffix = backend == :exact ? ".jld2" : "_chi$(χ).jld2"
        filenames = filter(
            name -> startswith(name, prefix) && endswith(name, suffix),
            readdir(directory),
        )
        if isempty(filenames)
            # Parameter tables change over time. Permit an old dataset only when
            # the directory contains a single, unambiguous trajectory duration.
            legacy_filenames = filter(
                name -> !isnothing(_sample_time(name)) && endswith(name, suffix),
                readdir(directory),
            )
            legacy_times = unique(_sample_time.(legacy_filenames))
            length(legacy_times) == 1 || error(
                "No files for current t=$t and no unique legacy t in $directory; " *
                "found legacy times $legacy_times",
            )
            filenames = legacy_filenames
            @warn "Using legacy OBF trajectories whose t differs from the current config" configured_t =
                t data_t = only(legacy_times)
        end
        sort!(filenames; by = _sample_index)
        isempty(filenames) && error("No sample files matching $(prefix)*$(suffix) in $directory")

        sample_num >= 0 || throw(ArgumentError("sample_num must be nonnegative"))
        if sample_num > 0
            filenames = filenames[1:min(sample_num, length(filenames))]
        end
        return joinpath.(directory, filenames)
    end

    function _process_obf_sample_file(
        path::AbstractString,
        start_period::Int,
        discard_final_periods::Int,
    )
        try
            sample = load(path, "sample")
            n_periods = size(sample, 1) ÷ OBF_LAYERS_PER_PERIOD
            stop_period = n_periods - discard_final_periods
            statistics = analyze_obf_trajectory(
                sample;
                start_period = start_period,
                stop_period = stop_period,
            )
            return (path = path, status = :success, statistics = statistics, error = nothing)
        catch error
            message = sprint(showerror, error, catch_backtrace())
            return (path = path, status = :failed, statistics = nothing, error = message)
        end
    end

    function _ensemble_stderr(values::AbstractArray, n_samples::Int)
        output_size = size(values)[1:(end - 1)]
        if n_samples == 1
            return fill(NaN, output_size)
        end
        return dropdims(std(values; dims = ndims(values)); dims = ndims(values)) ./
               sqrt(n_samples)
    end

    """
        measure_obf_distribution(backend, L, τ_idx, λ; χ=nothing,
                                 sample_num=0, discard_final_periods=2,
                                 output_root="exm/data/OBF/measurement_distribution")

    Analyze all matching independent Born trajectories and save the ensemble
    probability distribution. `sample_num=0` uses every matching file. Standard
    errors are computed across trajectory-level, time-averaged estimators, so
    correlated time slices are not incorrectly counted as independent samples.
    """
    function measure_obf_distribution(
        backend::Symbol,
        L::Int,
        τ_idx::Int,
        λ::Float64;
        χ::Union{Int,Nothing} = nothing,
        sample_num::Int = 0,
        discard_final_periods::Int = 2,
        output_root::AbstractString = "exm/data/OBF/measurement_distribution",
    )
        1 <= τ_idx <= length(τlis) ||
            throw(ArgumentError("τ_idx must be in 1:$(length(τlis))"))
        discard_final_periods >= 0 ||
            throw(ArgumentError("discard_final_periods must be nonnegative"))

        files = obf_sample_files(
            backend,
            L,
            τ_idx,
            λ;
            χ = χ,
            sample_num = sample_num,
        )
        _, _, start = get_born_dynamics_params(τ_idx, L, λ)
        start_period = max(1, start * L)

        println(
            "Processing $(length(files)) OBF trajectories for L=$L, γ=$(γlis[τ_idx]), " *
            "λ=$λ with $(nprocs()) processes...",
        )
        results = pmap(
            path -> _process_obf_sample_file(path, start_period, discard_final_periods),
            files,
        )

        failed = filter(result -> result.status == :failed, results)
        for result in failed
            @warn "Failed to process trajectory" path = result.path error = result.error
        end
        valid = filter(result -> result.status == :success, results)
        isempty(valid) && error("No valid OBF trajectories were found")

        statistics = getproperty.(valid, :statistics)
        n_valid = length(statistics)
        profiles = cat(getproperty.(statistics, :probability_profile)...; dims = 3)
        distributions = cat(getproperty.(statistics, :count_distribution)...; dims = 3)
        spatial_correlations =
            cat(getproperty.(statistics, :connected_spatial_correlation)...; dims = 3)
        layer_probabilities = hcat(getproperty.(statistics, :mean_probability)...)
        layer_magnetizations = hcat(getproperty.(statistics, :mean_magnetization)...)

        probability_profile = dropdims(mean(profiles; dims = 3); dims = 3)
        probability_profile_stderr = _ensemble_stderr(profiles, n_valid)
        count_distribution = dropdims(mean(distributions; dims = 3); dims = 3)
        count_distribution_stderr = _ensemble_stderr(distributions, n_valid)
        connected_spatial_correlation =
            dropdims(mean(spatial_correlations; dims = 3); dims = 3)
        connected_spatial_correlation_stderr =
            _ensemble_stderr(spatial_correlations, n_valid)
        mean_probability = vec(mean(layer_probabilities; dims = 2))
        mean_probability_stderr = vec(
            n_valid == 1 ? fill(NaN, OBF_LAYERS_PER_PERIOD) :
            std(layer_probabilities; dims = 2) ./ sqrt(n_valid),
        )
        mean_magnetization = vec(mean(layer_magnetizations; dims = 2))
        mean_magnetization_stderr = vec(
            n_valid == 1 ? fill(NaN, OBF_LAYERS_PER_PERIOD) :
            std(layer_magnetizations; dims = 2) ./ sqrt(n_valid),
        )

        first_statistics = first(statistics)
        active_site_mask = first_statistics.active_site_mask
        active_site_count = first_statistics.active_site_count
        n_stationary_periods = getproperty.(statistics, :n_stationary_periods)
        k_values = collect(0:L)
        magnetization_values = fill(NaN, OBF_LAYERS_PER_PERIOD, L + 1)
        valid_k_mask = falses(OBF_LAYERS_PER_PERIOD, L + 1)
        for layer = 1:OBF_LAYERS_PER_PERIOD
            n_active = active_site_count[layer]
            valid_k_mask[layer, 1:(n_active + 1)] .= true
            magnetization_values[layer, 1:(n_active + 1)] =
                (2 .* collect(0:n_active) .- n_active) ./ n_active
        end

        data_times = unique(_sample_time(basename(path)) for path in files)
        length(data_times) == 1 || error("Cannot mix trajectories with different t: $data_times")
        t = only(data_times)
        output_directory = joinpath(output_root, "L$(L)", "gammaind$(τ_idx)", "λ$(λ)")
        if backend == :mps
            output_directory = joinpath(output_directory, "chi$(χ)")
        end
        mkpath(output_directory)
        backend_suffix = backend == :exact ? "exact" : "mps_chi$(χ)"
        output_path = joinpath(
            output_directory,
            "distribution_t$(t)_$(backend_suffix)_samples$(n_valid).jld2",
        )

        save(
            output_path,
            "probability_profile",
            probability_profile,
            "probability_profile_stderr",
            probability_profile_stderr,
            "count_distribution",
            count_distribution,
            "count_distribution_stderr",
            count_distribution_stderr,
            "connected_spatial_correlation",
            connected_spatial_correlation,
            "connected_spatial_correlation_stderr",
            connected_spatial_correlation_stderr,
            "spatial_separation",
            collect(1:(L ÷ 2)),
            "mean_probability",
            mean_probability,
            "mean_probability_stderr",
            mean_probability_stderr,
            "mean_magnetization",
            mean_magnetization,
            "mean_magnetization_stderr",
            mean_magnetization_stderr,
            "active_site_mask",
            active_site_mask,
            "active_site_count",
            active_site_count,
            "valid_k_mask",
            valid_k_mask,
            "k_values",
            k_values,
            "magnetization_values",
            magnetization_values,
            "layer_labels",
            OBF_LAYER_LABELS,
            "n_valid_trajectories",
            n_valid,
            "n_failed_trajectories",
            length(failed),
            "n_stationary_periods",
            n_stationary_periods,
            "start_period",
            start_period,
            "discard_final_periods",
            discard_final_periods,
            "L",
            L,
            "tau_index",
            τ_idx,
            "tau",
            τlis[τ_idx],
            "gamma",
            γlis[τ_idx],
            "lambda",
            λ,
            "backend",
            String(backend),
            "chi",
            χ,
        )

        println("Done. Valid trajectories: $n_valid; failed: $(length(failed))")
        println("Saved to: $output_path")
        return output_path
    end
end

function print_usage()
    println("Usage:")
    println("  julia -p N measure_distribution.jl exact L tau_idx lambda [sample_num]")
    println("  julia -p N measure_distribution.jl mps   L tau_idx lambda chi [sample_num]")
    println()
    println("sample_num = 0 (default) uses all matching trajectory files.")
    println("Examples:")
    println("  julia -p 8 measure_distribution.jl exact 12 7 0.8 1000")
    println("  julia -p 8 measure_distribution.jl mps 64 7 0.8 128 1000")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS) || ARGS[1] in ("-h", "--help")
        print_usage()
        exit(isempty(ARGS) ? 1 : 0)
    end

    backend = Symbol(lowercase(ARGS[1]))
    backend in (:exact, :mps) || error("backend must be 'exact' or 'mps'")
    minimum_arguments = backend == :exact ? 4 : 5
    length(ARGS) >= minimum_arguments || begin
        print_usage()
        error("Not enough command-line arguments")
    end

    L = parse(Int, ARGS[2])
    τ_idx = parse(Int, ARGS[3])
    λ = parse(Float64, ARGS[4])
    if backend == :exact
        sample_num = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 0
        measure_obf_distribution(backend, L, τ_idx, λ; sample_num = sample_num)
    else
        χ = parse(Int, ARGS[5])
        sample_num = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 0
        measure_obf_distribution(
            backend,
            L,
            τ_idx,
            λ;
            χ = χ,
            sample_num = sample_num,
        )
    end
end
