using JLD2
using Printf
using Statistics

const TRAJECTORY_RE = r"^periods(\d+)_trajectory_seed(\d+)\.jld2$"

function trajectory_index(data_dir::AbstractString)
    isdir(data_dir) || error("Trajectory directory does not exist: $data_dir")
    indexed = Dict{Int,String}()
    periods = Set{Int}()
    for filename in readdir(data_dir)
        match_result = match(TRAJECTORY_RE, filename)
        isnothing(match_result) && continue
        period_count = parse(Int, match_result.captures[1])
        seed = parse(Int, match_result.captures[2])
        haskey(indexed, seed) && error("Duplicate trajectory seed $seed in $data_dir")
        indexed[seed] = joinpath(data_dir, filename)
        push!(periods, period_count)
    end
    isempty(indexed) && error("No trajectory files found in $data_dir")
    length(periods) == 1 || error("Mixed period counts in $data_dir: $(sort!(collect(periods)))")
    return indexed, only(periods)
end

function validate_metadata(path::AbstractString, expected_periods::Int)
    data = JLD2.load(path)
    sample = data["sample"]
    ndims(sample) == 2 || error("sample in $path is not a matrix")
    size(sample, 1) == 2expected_periods || error(
        "sample in $path has $(size(sample, 1)) layers, expected $(2expected_periods)",
    )
    return (
        L = Int(data["L"]),
        periods = Int(data["periods"]),
        tau_idx = Int(data["tau_idx"]),
        tau = Float64(data["tau"]),
        gamma = Float64(data["gamma"]),
        sample_size = size(sample),
        initial_state = String(data["initial_state"]),
    )
end

function load_records(paths::Vector{String})
    first_sample = JLD2.load(first(paths), "sample")
    n_layers, n_sites = size(first_sample)
    records = BitArray(undef, length(paths), n_layers, n_sites)
    Threads.@threads :dynamic for trajectory_index in eachindex(paths)
        sample = JLD2.load(paths[trajectory_index], "sample")
        size(sample) == (n_layers, n_sites) || error(
            "Inconsistent sample size $(size(sample)) in $(paths[trajectory_index]); " *
            "expected $((n_layers, n_sites))",
        )
        @views records[trajectory_index, :, :] .= sample
    end
    return records
end

function parity_layers(
    n_layers::Int,
    burnin_periods::Int,
    parity::Symbol;
    tail_trim_periods::Int = 0,
)
    0 <= burnin_periods < n_layers ÷ 2 || error(
        "burnin_periods must lie in 0:$(n_layers ÷ 2 - 1), got $burnin_periods",
    )
    0 <= tail_trim_periods < n_layers ÷ 2 || error(
        "tail_trim_periods must lie in 0:$(n_layers ÷ 2 - 1), got $tail_trim_periods",
    )
    burnin_periods + tail_trim_periods < n_layers ÷ 2 || error(
        "Burn-in and tail trim leave no measurement periods",
    )
    target_parity = parity === :odd ? 1 : parity === :even ? 0 : error(
        "parity must be :odd or :even",
    )
    first_layer = 2burnin_periods + 1
    last_layer = n_layers - 2tail_trim_periods
    return [layer for layer in first_layer:last_layer if layer % 2 == target_parity]
end

@inline function low_mask(nbits::Int)
    0 <= nbits <= 128 || error("UInt128 packing requires 0 <= nbits <= 128")
    nbits == 128 && return typemax(UInt128)
    return (one(UInt128) << nbits) - one(UInt128)
end

"""
Return per-trajectory magnetizations and raw autocorrelations.

For each trajectory, the raw correlation at lag `k` is averaged over every
valid time origin and all measurement channels in the chosen layer parity.
The selected time series is packed into `UInt128`, which is sufficient for the
current 70-period steady-state window and makes the XOR/Hamming calculation
substantially cheaper than scalar Boolean loops.
"""
function trajectory_statistics(
    records::BitArray{3},
    layers::Vector{Int},
    max_lag::Int,
)
    n_trajectories, _, n_sites = size(records)
    n_times = length(layers)
    n_times <= 128 || error(
        "Selected window has $n_times time points; UInt128 implementation supports at most 128",
    )
    0 <= max_lag < n_times || error("max_lag must lie in 0:$(n_times - 1)")

    magnetization = zeros(Float64, n_trajectories)
    raw_correlation = zeros(Float64, n_trajectories, max_lag + 1)
    spin_records = Matrix{Float64}(undef, n_trajectories, n_times * n_sites)
    masks = [low_mask(n_times - lag) for lag in 0:max_lag]

    Threads.@threads :dynamic for trajectory in 1:n_trajectories
        number_of_ones = 0
        @inbounds for site in 1:n_sites
            packed = zero(UInt128)
            for time_index in 1:n_times
                outcome = records[trajectory, layers[time_index], site]
                spin_records[trajectory, (time_index-1)*n_sites+site] =
                    outcome ? 1.0 : -1.0
                outcome && (packed |= one(UInt128) << (time_index - 1))
            end
            number_of_ones += count_ones(packed)
            for lag in 0:max_lag
                mismatch = count_ones(xor(packed, packed >> lag) & masks[lag+1])
                raw_correlation[trajectory, lag+1] += mismatch
            end
        end
        magnetization[trajectory] =
            2number_of_ones / (n_times * n_sites) - 1
        for lag in 0:max_lag
            pair_count = (n_times - lag) * n_sites
            raw_correlation[trajectory, lag+1] =
                1 - 2raw_correlation[trajectory, lag+1] / pair_count
        end
    end
    return magnetization, raw_correlation, spin_records
end

function summarize_profile(
    magnetization::Vector{Float64},
    raw::Matrix{Float64},
    spin_records::Matrix{Float64},
    n_times::Int,
    n_sites::Int,
)
    n_trajectories, n_lags = size(raw)
    length(magnetization) == n_trajectories || error("Trajectory count mismatch")
    size(spin_records) == (n_trajectories, n_times * n_sites) ||
        error("Spin-record matrix size mismatch")

    mean_spin_by_coordinate = vec(mean(spin_records; dims = 1))
    mean_spin = mean(mean_spin_by_coordinate)
    # Build all time-dependent centering terms at once.  For lag k,
    # centered_trajectory[r,k] equals the time/site average of
    # (s_r(t,j)-mu(t,j))*(s_r(t+k,j)-mu(t+k,j)).  The matrix product evaluates
    # the two linear-in-s terms for every trajectory and lag.
    centering_weights = zeros(Float64, n_times * n_sites, n_lags)
    mean_product = zeros(Float64, n_lags)
    @inbounds for lag in 0:n_lags-1
        pair_count = (n_times - lag) * n_sites
        for time_index in 1:n_times-lag
            for site in 1:n_sites
                first_coordinate = (time_index-1)*n_sites + site
                second_coordinate = (time_index+lag-1)*n_sites + site
                first_mean = mean_spin_by_coordinate[first_coordinate]
                second_mean = mean_spin_by_coordinate[second_coordinate]
                centering_weights[first_coordinate, lag+1] += second_mean / pair_count
                centering_weights[second_coordinate, lag+1] += first_mean / pair_count
                mean_product[lag+1] += first_mean * second_mean / pair_count
            end
        end
    end
    centered_trajectory = raw .- spin_records * centering_weights .+ mean_product'

    # Correct the O(1/N) bias from estimating the coordinate-wise means with
    # the same trajectory ensemble.  This is numerically tiny at N=30000 but
    # makes the estimator exactly equal to the usual unbiased covariance.
    covariance_correction = n_trajectories / (n_trajectories - 1)
    centered_trajectory .*= covariance_correction
    connected = vec(mean(centered_trajectory; dims = 1))
    influence_connected = centered_trajectory .- connected'
    stderr_connected = vec(std(influence_connected; dims = 1)) ./ sqrt(n_trajectories)

    connected[1] > 0 || error(
        "Zero outcome variance: C(0)=$(connected[1]); normalized correlation is undefined",
    )
    normalized = connected ./ connected[1]
    influence_normalized = Matrix{Float64}(undef, n_trajectories, n_lags)
    @inbounds for lag_index in 1:n_lags
        for trajectory in 1:n_trajectories
            influence_normalized[trajectory, lag_index] =
                (influence_connected[trajectory, lag_index] -
                 normalized[lag_index] * influence_connected[trajectory, 1]) /
                connected[1]
        end
    end
    stderr_normalized =
        vec(std(influence_normalized; dims = 1)) ./ sqrt(n_trajectories)

    mean_spin_by_time = zeros(Float64, n_times)
    @inbounds for time_index in 1:n_times
        coordinates = (time_index-1)*n_sites+1:time_index*n_sites
        mean_spin_by_time[time_index] = mean(@view mean_spin_by_coordinate[coordinates])
    end

    return (
        mean_spin = mean_spin,
        probability_one = (mean_spin + 1) / 2,
        probability_one_by_time = (mean_spin_by_time .+ 1) ./ 2,
        connected = connected,
        stderr_connected = stderr_connected,
        normalized = normalized,
        stderr_normalized = stderr_normalized,
        influence_connected = influence_connected,
        influence_normalized = influence_normalized,
    )
end

function mean_profile(odd, even)
    connected = (odd.connected .+ even.connected) ./ 2
    influence_connected =
        (odd.influence_connected .+ even.influence_connected) ./ 2
    stderr_connected =
        vec(std(influence_connected; dims = 1)) ./ sqrt(size(influence_connected, 1))
    normalized = connected ./ connected[1]
    influence_normalized = similar(influence_connected)
    @inbounds for lag_index in axes(influence_connected, 2)
        influence_normalized[:, lag_index] .=
            (influence_connected[:, lag_index] .-
             normalized[lag_index] .* influence_connected[:, 1]) ./ connected[1]
    end
    stderr_normalized =
        vec(std(influence_normalized; dims = 1)) ./ sqrt(size(influence_connected, 1))
    return (
        mean_spin = (odd.mean_spin + even.mean_spin) / 2,
        probability_one = (odd.probability_one + even.probability_one) / 2,
        probability_one_by_time =
            (odd.probability_one_by_time .+ even.probability_one_by_time) ./ 2,
        connected = connected,
        stderr_connected = stderr_connected,
        normalized = normalized,
        stderr_normalized = stderr_normalized,
        influence_connected = influence_connected,
        influence_normalized = influence_normalized,
    )
end

function analyze_sector(
    records::BitArray{3},
    burnin_periods::Int,
    max_lag::Int;
    tail_trim_periods::Int = 0,
)
    n_layers = size(records, 2)
    odd_layers = parity_layers(
        n_layers,
        burnin_periods,
        :odd;
        tail_trim_periods = tail_trim_periods,
    )
    even_layers = parity_layers(
        n_layers,
        burnin_periods,
        :even;
        tail_trim_periods = tail_trim_periods,
    )
    length(odd_layers) == length(even_layers) || error(
        "Odd/even steady-state windows have different lengths: " *
        "$(length(odd_layers)) and $(length(even_layers))",
    )
    odd_magnetization, odd_raw, odd_spin_records =
        trajectory_statistics(records, odd_layers, max_lag)
    even_magnetization, even_raw, even_spin_records =
        trajectory_statistics(records, even_layers, max_lag)
    odd = summarize_profile(
        odd_magnetization,
        odd_raw,
        odd_spin_records,
        length(odd_layers),
        size(records, 3),
    )
    even = summarize_profile(
        even_magnetization,
        even_raw,
        even_spin_records,
        length(even_layers),
        size(records, 3),
    )
    return (odd = odd, even = even, mean = mean_profile(odd, even))
end

function profile_matrix(profiles, field::Symbol)
    return hcat([getproperty(getproperty(profiles, parity), field) for parity in (:odd, :even, :mean)]...)
end

function positive_sequence_time(normalized::AbstractVector{<:Real})
    total = 0.5
    for value in @view normalized[2:end]
        value > 0 || break
        total += value
    end
    return total
end

function write_tsv(
    path::AbstractString,
    y1_profiles,
    ytau_profiles,
    paired_stderr_connected::Matrix{Float64},
    independent_stderr_connected::Matrix{Float64},
    paired_stderr_normalized::Matrix{Float64},
    independent_stderr_normalized::Matrix{Float64},
)
    open(path, "w") do io
        println(
            io,
            join(
                [
                    "parity",
                    "lag_periods",
                    "C_y1",
                    "stderr_C_y1",
                    "rho_y1",
                    "stderr_rho_y1",
                    "C_ytau",
                    "stderr_C_ytau",
                    "rho_ytau",
                    "stderr_rho_ytau",
                    "delta_C_ytau_minus_y1",
                    "paired_stderr_delta_C",
                    "independent_stderr_delta_C",
                    "delta_rho_ytau_minus_y1",
                    "paired_stderr_delta_rho",
                    "independent_stderr_delta_rho",
                ],
                '\t',
            ),
        )
        for (parity_index, parity) in enumerate((:odd, :even, :mean))
            y1 = getproperty(y1_profiles, parity)
            ytau = getproperty(ytau_profiles, parity)
            for lag_index in eachindex(y1.connected)
                values = (
                    parity,
                    lag_index - 1,
                    y1.connected[lag_index],
                    y1.stderr_connected[lag_index],
                    y1.normalized[lag_index],
                    y1.stderr_normalized[lag_index],
                    ytau.connected[lag_index],
                    ytau.stderr_connected[lag_index],
                    ytau.normalized[lag_index],
                    ytau.stderr_normalized[lag_index],
                    ytau.connected[lag_index] - y1.connected[lag_index],
                    paired_stderr_connected[lag_index, parity_index],
                    independent_stderr_connected[lag_index, parity_index],
                    ytau.normalized[lag_index] - y1.normalized[lag_index],
                    paired_stderr_normalized[lag_index, parity_index],
                    independent_stderr_normalized[lag_index, parity_index],
                )
                println(io, join(values, '\t'))
            end
        end
    end
end

function compare_sectors(
    y1_dir::AbstractString,
    ytau_dir::AbstractString,
    output_prefix::AbstractString;
    burnin_periods::Int = 10,
    max_lag::Int = 67,
    max_samples::Int = 0,
    tail_trim_periods::Int = 2,
)
    y1_index, y1_periods = trajectory_index(y1_dir)
    ytau_index, ytau_periods = trajectory_index(ytau_dir)
    y1_periods == ytau_periods || error(
        "Period mismatch: y=1 has $y1_periods, y=tau has $ytau_periods",
    )
    y1_seeds = sort!(collect(keys(y1_index)))
    ytau_seeds = sort!(collect(keys(ytau_index)))
    y1_seeds == ytau_seeds || error(
        "The two sectors do not have identical seed sets; paired comparison is invalid",
    )
    if max_samples > 0
        max_samples <= length(y1_seeds) || error(
            "Requested max_samples=$max_samples but only $(length(y1_seeds)) exist",
        )
        y1_seeds = y1_seeds[1:max_samples]
    end

    y1_metadata = validate_metadata(y1_index[first(y1_seeds)], y1_periods)
    ytau_metadata = validate_metadata(ytau_index[first(y1_seeds)], ytau_periods)
    comparable_fields = (:L, :periods, :tau_idx, :tau, :gamma, :sample_size)
    for field in comparable_fields
        getproperty(y1_metadata, field) == getproperty(ytau_metadata, field) || error(
            "Metadata mismatch for $field: $(getproperty(y1_metadata, field)) versus " *
            "$(getproperty(ytau_metadata, field))",
        )
    end

    n_steady_times = y1_metadata.periods - burnin_periods - tail_trim_periods
    max_lag < n_steady_times || error(
        "max_lag=$max_lag must be smaller than the $n_steady_times retained periods",
    )
    @info "Loading y=1 records" trajectories = length(y1_seeds) threads = Threads.nthreads()
    y1_records = load_records([y1_index[seed] for seed in y1_seeds])
    @info "Loading y=tau records" trajectories = length(y1_seeds) threads = Threads.nthreads()
    ytau_records = load_records([ytau_index[seed] for seed in y1_seeds])

    @info "Computing y=1 temporal correlations"
    y1_profiles = analyze_sector(
        y1_records,
        burnin_periods,
        max_lag;
        tail_trim_periods = tail_trim_periods,
    )
    @info "Computing y=tau temporal correlations"
    ytau_profiles = analyze_sector(
        ytau_records,
        burnin_periods,
        max_lag;
        tail_trim_periods = tail_trim_periods,
    )

    parities = (:odd, :even, :mean)
    n_lags = max_lag + 1
    paired_stderr_connected = zeros(n_lags, length(parities))
    independent_stderr_connected = zeros(n_lags, length(parities))
    paired_stderr_normalized = zeros(n_lags, length(parities))
    independent_stderr_normalized = zeros(n_lags, length(parities))
    n_trajectories = length(y1_seeds)
    for (parity_index, parity) in enumerate(parities)
        y1 = getproperty(y1_profiles, parity)
        ytau = getproperty(ytau_profiles, parity)
        paired_stderr_connected[:, parity_index] =
            vec(std(ytau.influence_connected .- y1.influence_connected; dims = 1)) ./
            sqrt(n_trajectories)
        independent_stderr_connected[:, parity_index] =
            hypot.(y1.stderr_connected, ytau.stderr_connected)
        paired_stderr_normalized[:, parity_index] =
            vec(std(ytau.influence_normalized .- y1.influence_normalized; dims = 1)) ./
            sqrt(n_trajectories)
        independent_stderr_normalized[:, parity_index] =
            hypot.(y1.stderr_normalized, ytau.stderr_normalized)
    end

    output_dir = dirname(output_prefix)
    output_dir == "." || mkpath(output_dir)
    jld2_path = output_prefix * ".jld2"
    tsv_path = output_prefix * ".tsv"
    mean_spin_y1 = [getproperty(y1_profiles, p).mean_spin for p in parities]
    mean_spin_ytau = [getproperty(ytau_profiles, p).mean_spin for p in parities]
    probability_one_y1 = [getproperty(y1_profiles, p).probability_one for p in parities]
    probability_one_ytau = [getproperty(ytau_profiles, p).probability_one for p in parities]
    probability_one_by_time_y1 =
        hcat([getproperty(y1_profiles, p).probability_one_by_time for p in parities]...)
    probability_one_by_time_ytau =
        hcat([getproperty(ytau_profiles, p).probability_one_by_time for p in parities]...)
    connected_y1 = profile_matrix(y1_profiles, :connected)
    connected_ytau = profile_matrix(ytau_profiles, :connected)
    normalized_y1 = profile_matrix(y1_profiles, :normalized)
    normalized_ytau = profile_matrix(ytau_profiles, :normalized)
    stderr_connected_y1 = profile_matrix(y1_profiles, :stderr_connected)
    stderr_connected_ytau = profile_matrix(ytau_profiles, :stderr_connected)
    stderr_normalized_y1 = profile_matrix(y1_profiles, :stderr_normalized)
    stderr_normalized_ytau = profile_matrix(ytau_profiles, :stderr_normalized)
    integrated_time_y1 = [
        positive_sequence_time(getproperty(y1_profiles, p).normalized) for p in parities
    ]
    integrated_time_ytau = [
        positive_sequence_time(getproperty(ytau_profiles, p).normalized) for p in parities
    ]

    JLD2.jldsave(
        jld2_path;
        L = y1_metadata.L,
        tau_idx = y1_metadata.tau_idx,
        tau = y1_metadata.tau,
        gamma = y1_metadata.gamma,
        periods = y1_metadata.periods,
        burnin_periods = burnin_periods,
        tail_trim_periods = tail_trim_periods,
        retained_periods = n_steady_times,
        max_lag = max_lag,
        number_of_trajectories = n_trajectories,
        seeds = y1_seeds,
        sector_y1_initial_state = y1_metadata.initial_state,
        sector_ytau_initial_state = ytau_metadata.initial_state,
        parity = String.(parities),
        lag_periods = collect(0:max_lag),
        mean_spin_y1 = mean_spin_y1,
        mean_spin_ytau = mean_spin_ytau,
        probability_one_y1 = probability_one_y1,
        probability_one_ytau = probability_one_ytau,
        probability_one_by_retained_time_y1 = probability_one_by_time_y1,
        probability_one_by_retained_time_ytau = probability_one_by_time_ytau,
        connected_y1 = connected_y1,
        connected_ytau = connected_ytau,
        stderr_connected_y1 = stderr_connected_y1,
        stderr_connected_ytau = stderr_connected_ytau,
        normalized_y1 = normalized_y1,
        normalized_ytau = normalized_ytau,
        stderr_normalized_y1 = stderr_normalized_y1,
        stderr_normalized_ytau = stderr_normalized_ytau,
        delta_connected_ytau_minus_y1 = connected_ytau - connected_y1,
        paired_stderr_delta_connected = paired_stderr_connected,
        independent_stderr_delta_connected = independent_stderr_connected,
        delta_normalized_ytau_minus_y1 = normalized_ytau - normalized_y1,
        paired_stderr_delta_normalized = paired_stderr_normalized,
        independent_stderr_delta_normalized = independent_stderr_normalized,
        positive_sequence_integrated_time_y1 = integrated_time_y1,
        positive_sequence_integrated_time_ytau = integrated_time_ytau,
    )
    write_tsv(
        tsv_path,
        y1_profiles,
        ytau_profiles,
        paired_stderr_connected,
        independent_stderr_connected,
        paired_stderr_normalized,
        independent_stderr_normalized,
    )

    println("\nTemporal measurement-outcome correlation comparison")
    println("  y=1:   $y1_dir ($(y1_metadata.initial_state))")
    println("  y=tau: $ytau_dir ($(ytau_metadata.initial_state))")
    println("  trajectories: $n_trajectories (identical seed set)")
    println(
        "  retained periods: $n_steady_times after burn-in $burnin_periods " *
        "and tail trim $tail_trim_periods",
    )
    println("  lag range: 0:$max_lag periods")
    for (parity_index, parity) in enumerate(parities)
        @printf(
            "  %-4s  P(1): %.8f -> %.8f; tau_int(+): %.6f -> %.6f\n",
            String(parity),
            probability_one_y1[parity_index],
            probability_one_ytau[parity_index],
            integrated_time_y1[parity_index],
            integrated_time_ytau[parity_index],
        )
        @printf(
            "        max P(1) drift across retained times: %.8f, %.8f\n",
            maximum(@view probability_one_by_time_y1[:, parity_index]) -
            minimum(@view probability_one_by_time_y1[:, parity_index]),
            maximum(@view probability_one_by_time_ytau[:, parity_index]) -
            minimum(@view probability_one_by_time_ytau[:, parity_index]),
        )
    end
    println("  JLD2: $jld2_path")
    println("  TSV:  $tsv_path")
    return jld2_path, tsv_path
end

function usage()
    println(
        "Usage: julia --threads=N --project=. " *
        "exm/Bulk_measure/measurement_outcome_temporal_sector_comparison.jl " *
        "Y1_DIR YTAU_DIR OUTPUT_PREFIX [BURNIN_PERIODS=10] [MAX_LAG=67] " *
        "[MAX_SAMPLES=0] [TAIL_TRIM_PERIODS=2]",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    3 <= length(ARGS) <= 7 || begin
        usage()
        exit(1)
    end
    compare_sectors(
        ARGS[1],
        ARGS[2],
        ARGS[3];
        burnin_periods = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 10,
        max_lag = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 67,
        max_samples = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 0,
        tail_trim_periods = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 2,
    )
end
