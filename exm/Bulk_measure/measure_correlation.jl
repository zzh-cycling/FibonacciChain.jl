using Distributed

# Add workers if not already added
nprocs() == 1 && addprocs(Sys.CPU_THREADS - 1)

@everywhere using JLD
@everywhere using LinearAlgebra
@everywhere using Statistics

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[7] = log(1 + √2)  # atanh(1/√2) = log(1 + √2)
τlis[end] = 1000.0  # Last value is for γ=1

# Broadcast parameters to all workers
@everywhere γlis = $γlis
@everywhere τlis = $τlis

L_list = collect(8:2:20)

@everywhere function data_fig_path(inds, L, index, τlis)
    τ = τlis[inds]
    D, _, _ = get_cfg_params_Born(τ, L)
    DATA_path = "exm/data/Bulk_measure/Samples_monitored_dynamics/L$(L)/τ$(τ)/D$(div(D,L))_Samples$(index).jld"
    return DATA_path
end

@everywhere function save_data_filename(L, inds)
    τ = τlis[inds]
    D, _, _ = get_cfg_params_Born(τ, L)
    return "exm/data/Bulk_measure/records_correlation/L$(L)/τ$(τ)/record_correlation_D$(div(D,L)).jld"
end

@everywhere function get_cfg_params_Born(τ, L)
    cfg = Dict(
        atanh(0.1)  => (2500L, 1000, 750L),
        atanh(0.2)  => (500L,  100, 120L),
        atanh(0.3)  => (120L,  48, 50L),
        atanh(0.4)  => (100L,  40, 40L),
        atanh(0.5)  => (80L,   32, 20L),
        atanh(0.6)  => (45L,   20, 15L),
        log(1 + √2) => (35L,   14, 10L),
        atanh(0.8)  => (25L,   10, 5L),
        atanh(0.9)  => (8L,    4, 2L),
        atanh(0.95) => (8L,    4, 2L),
        atanh(0.999)=> (5L,    2, 1L),
    )
    D, step, start = get(cfg, τ, (5L, 2, L))
    inds = collect(1:step:div(D,2))
    avg_range = start:div(D,2)-5
    return D, inds, avg_range
end

# """
#     layer_corr(sample::Vector{Int8})

# Compute connected spatial correlation function C(d) for a PBC chain.

# For a chain of spins s_i ∈ {-1, +1}, computes the connected correlator:
#     C(d) = ⟨s_i s_{i+d}⟩ - ⟨s_i⟩⟨s_{i+d}⟩
#          = (1/N) ∑_{i=1}^{N} s_i ⋅ s_{i+d} - m²

# where m = (1/N) ∑_i s_i is the magnetization, and the index i+d wraps 
# around due to periodic boundary conditions.

# # Arguments
# - `sample::Vector{Int8}`: Measurement outcomes (0 or 1)

# # Returns
# - `Vector{Float64}`: Connected correlation function C(d) for d = 1, 2, ..., N÷2
#   (only up to N÷2 due to PBC symmetry: C(d) = C(N-d))
# """
@everywhere function layer_corr(sample::Vector{Int8})
    # Convert 0,1 to -1,+1
    spins = 2 .* sample .- 1

    N = length(spins)
    # Compute magnetization ⟨s⟩
    m = mean(spins)
    
    # For PBC, correlation at distance d equals that at N-d, so only compute up to N÷2
    max_d = N ÷ 2
    corr_dlis = zeros(max_d)
    
    for d in 1:max_d
        # Spatial average: sum over all starting positions i
        for i in 1:N
            # PBC: index wraps around using mod1
            j = mod1(i + d, N)  # mod1 gives 1-based index
            corr_dlis[d] += spins[i] * spins[j]
        end
        corr_dlis[d] /= N  # ⟨s_i s_{i+d}⟩
        corr_dlis[d] -= m^2  # Subtract ⟨s_i⟩⟨s_{i+d}⟩ = m²
    end

    return corr_dlis, m
end


@everywhere function process_single_sample(inds, L, sample_idx, τlis)
    data_path = data_fig_path(inds, L, sample_idx, τlis)
    @info "Processing sample $sample_idx"
    if !isfile(data_path)
        @warn "File not found: $data_path"
        return nothing
    end
    
    sample_data = load(data_path, "sample")  # (time_steps, L) matrix
    _, _, avg_range = get_cfg_params_Born(τlis[inds], L)
    
    max_d = L ÷ 4
    corr_time_avg_odd = zeros(max_d)
    corr_time_avg_even = zeros(max_d)
    m_time_avg_odd = 0.0
    m_time_avg_even = 0.0
    
    n_odd = 0
    n_even = 0
    
    for t in avg_range[1]:2:avg_range[end]
        sample_t = Int8.(sample_data[t, :])
        corr_t, m_t = layer_corr(sample_t)
        corr_time_avg_odd .+= corr_t
        m_time_avg_odd += m_t
        n_odd += 1
    end

    for t in avg_range[2]:2:avg_range[end]
        sample_t = Int8.(sample_data[t, :])
        corr_t, m_t = layer_corr(sample_t)
        corr_time_avg_even .+= corr_t
        m_time_avg_even += m_t
        n_even += 1
    end

    corr_time_avg_odd ./= n_odd
    corr_time_avg_even ./= n_even
    m_time_avg_odd /= n_odd
    m_time_avg_even /= n_even
    
    GC.gc() 
    return corr_time_avg_odd, corr_time_avg_even, m_time_avg_odd, m_time_avg_even
end

# """
#     average_corr(inds, L; sample_num=10000)

# Compute average spatial correlation function over time window and sample ensemble.
# Parallelized version using Distributed.jl.

# # Arguments
# - `inds`: Index for γ/τ parameter
# - `L`: System size
# - `sample_num`: Number of samples to process (default: 10000)

# # Returns
# - `corr_avg_odd::Vector{Float64}`: Ensemble-averaged correlation for odd layers
# - `corr_std_odd::Vector{Float64}`: Standard error for odd layers
# - `corr_avg_even::Vector{Float64}`: Ensemble-averaged correlation for even layers
# - `corr_std_even::Vector{Float64}`: Standard error for even layers
# - `m_avg_odd::Float64`: Ensemble-averaged magnetization for odd layers
# - `m_std_odd::Float64`: Standard error for odd magnetization
# - `m_avg_even::Float64`: Ensemble-averaged magnetization for even layers
# - `m_std_even::Float64`: Standard error for even magnetization
# """
function average_corr(inds, L; sample_num=10000)
    max_d = L ÷ 4
    
    println("Processing $sample_num samples for L=$L, γ=$(γlis[inds]) with $(nprocs()) processes...")
    
    # Parallel map over all samples
    results = pmap(i -> process_single_sample(inds, L, i, τlis), 1:sample_num)
    
    # Filter out failed samples (nothing values)
    valid_results = filter(!isnothing, results)
    n_valid = length(valid_results)
    
    if n_valid == 0
        error("No valid samples found!")
    end
    
    if n_valid < sample_num
        @warn "Only $n_valid out of $sample_num samples were valid"
    end
    
    # Separate odd and even results (correlation and magnetization)
    corr_odd_list = [r[1] for r in valid_results]
    corr_even_list = [r[2] for r in valid_results]
    m_odd_list = [r[3] for r in valid_results]
    m_even_list = [r[4] for r in valid_results]
    
    # Stack correlation results into matrices
    corr_ensemble_odd = hcat(corr_odd_list...)'   # (n_valid, max_d)
    corr_ensemble_even = hcat(corr_even_list...)'  # (n_valid, max_d)
    
    # Average correlation over all samples
    corr_avg_odd = mean(corr_ensemble_odd, dims=1)[:]
    corr_std_odd = std(corr_ensemble_odd, dims=1)[:] ./ sqrt(n_valid)
    
    corr_avg_even = mean(corr_ensemble_even, dims=1)[:]
    corr_std_even = std(corr_ensemble_even, dims=1)[:] ./ sqrt(n_valid)
    
    # Average magnetization over all samples
    m_avg_odd = mean(m_odd_list)
    m_std_odd = std(m_odd_list) / sqrt(n_valid)
    
    m_avg_even = mean(m_even_list)
    m_std_even = std(m_even_list) / sqrt(n_valid)
    
    println("Done! Valid samples: $n_valid")
    println("Magnetization: odd = $m_avg_odd ± $m_std_odd, even = $m_avg_even ± $m_std_even")
    
    # Save results
    save_path = save_data_filename(L, inds)
    mkpath(dirname(save_path))
    save(save_path, 
         "corr_avg_odd", corr_avg_odd, 
         "corr_std_odd", corr_std_odd,
         "corr_avg_even", corr_avg_even,
         "corr_std_even", corr_std_even,
         "m_avg_odd", m_avg_odd,
         "m_std_odd", m_std_odd,
         "m_avg_even", m_avg_even,
         "m_std_even", m_std_even)
    println("Saved to: $save_path")
    
    return corr_avg_odd, corr_std_odd, corr_avg_even, corr_std_even, m_avg_odd, m_std_odd, m_avg_even, m_std_even
end

# """
#     main(inds, L_list; sample_num=10000)

# Main function to compute and save correlation functions for multiple system sizes.
# """
# function main(inds, L_list; sample_num=10000)
#     println("="^60)
#     println("Computing correlation functions")
#     println("γ = $(γlis[inds]), τ = $(τlis[inds])")
#     println("System sizes: $L_list")
#     println("Number of samples: $sample_num")
#     println("Number of workers: $(nprocs())")
#     println("="^60)
    
#     for L in L_list
#         println("\n--- L = $L ---")
#         @time corr_avg, corr_std = average_corr(inds, L; sample_num=sample_num)
        
#         # Save results
#         save_path = save_data_filename(L, inds)
#         mkpath(dirname(save_path))
#         save(save_path, 
#              "corr_avg", corr_avg, 
#              "corr_std", corr_std,
#              "L", L,
#              "gamma", γlis[inds],
#              "tau", τlis[inds])
#         println("Saved to: $save_path")
#     end
    
#     println("\n" * "="^60)
#     println("All done!")
#     println("="^60)
# end



function corr_summary(L)
    max_d = L ÷ 2 - 1
    corr_avg_odd_list = zeros(length(γlis), max_d)
    corr_std_odd_list = zeros(length(γlis), max_d)
    corr_avg_even_list = zeros(length(γlis), max_d)
    corr_std_even_list = zeros(length(γlis), max_d)
    m_avg_even_list = zeros(length(γlis), L ÷ 2)
    m_std_even_list = zeros(length(γlis), L ÷ 2)
    m_avg_odd_list = zeros(length(γlis), L ÷ 2)
    m_std_odd_list = zeros(length(γlis), L ÷ 2)
    for inds in 1:length(γlis)
            save_path = save_data_filename(L, inds)
            if isfile(save_path)
                data = load(save_path)
                corr_avg_odd = data["corr_avg_odd"]
                corr_avg_even = data["corr_avg_even"]
                corr_std_odd = data["corr_std_odd"]
                corr_std_even = data["corr_std_even"]
                m_avg_odd = data["mlis_avg_odd"]
                m_avg_even = data["mlis_avg_even"]
                m_std_odd = data["mlis_std_odd"]
                m_std_even = data["mlis_std_even"]
                corr_avg_odd_list[inds, :] = corr_avg_odd
                corr_avg_even_list[inds, :] = corr_avg_even
                corr_std_odd_list[inds, :] = corr_std_odd
                corr_std_even_list[inds, :] = corr_std_even
                m_std_odd_list[inds, :] = m_std_odd
                m_std_even_list[inds, :] = m_std_even
                m_avg_odd_list[inds, :] = m_avg_odd
                m_avg_even_list[inds, :] = m_avg_even
            end
        save_path_summary = "exm/data/Bulk_measure/records_correlation/L$(L)/summary_correlation.jld"
        mkpath(dirname(save_path_summary))
        save(save_path_summary,
            "corr_avg_odd_list", corr_avg_odd_list,
            "corr_avg_even_list", corr_avg_even_list,
            "corr_std_odd_list", corr_std_odd_list,
            "corr_std_even_list", corr_std_even_list,
            "m_std_odd_list", m_std_odd_list,
            "m_std_even_list", m_std_even_list,
            "m_avg_odd_list", m_avg_odd_list,
            "m_avg_even_list", m_avg_even_list)
    end
end

@everywhere function process_single_sample_time(inds, L, sample_idx, τlis)
    # Compute time-averaged correlation at site 1 only (no spatial average)
    # C(d) = ⟨s_1 s_{1+d}⟩_t - ⟨s_1⟩_t ⟨s_{1+d}⟩_t
    # Also compute magnetization at each site: m_i = ⟨s_i⟩_t
    data_path = data_fig_path(inds, L, sample_idx, τlis)
    @info "Processing sample $sample_idx"
    if !isfile(data_path)
        @warn "File not found: $data_path"
        return nothing
    end
    
    _, _, avg_range = get_cfg_params_Born(τlis[inds], L)
    sample_data = load(data_path, "sample")
    
    # Extract odd/even time slices and convert to ±1
    record_data_odd = 2 .* Int8.(sample_data[avg_range[1:2:end], :]) .- 1
    record_data_even = 2 .* Int8.(sample_data[avg_range[2:2:end], :]) .- 1
    
    max_d = L ÷ 2 - 1
    corr_time_odd = zeros(max_d)
    corr_time_even = zeros(max_d)
    
    # Compute magnetization at each site: m_i = ⟨s_i⟩_t (time average)
    mlis_odd = vec(mean(record_data_odd, dims=1))   # (L,) vector
    mlis_even = vec(mean(record_data_even, dims=1)) # (L,) vector
    
    # Odd layers: compute connected correlation C(d) = ⟨s_1 s_{1+d}⟩_t - ⟨s_1⟩_t ⟨s_{1+d}⟩_t
    for d in 1:max_d
        corr_time_odd[d] = mean(record_data_odd[:, 1] .* record_data_odd[:, d+1]) - mlis_odd[1] * mlis_odd[d+1]
    end
    
    # Even layers: compute connected correlation
    for d in 1:max_d
        corr_time_even[d] = mean(record_data_even[:, 1] .* record_data_even[:, d+1]) - mlis_even[1] * mlis_even[d+1]
    end
    
    GC.gc()
    return corr_time_odd, corr_time_even, mlis_odd, mlis_even
end

function average_corr_time(inds, L; sample_num=10000)
    # Compute time-averaged correlation at site 1 only (no spatial average)
    # C(d) = ⟨s_1 s_{1+d}⟩ - ⟨s_1⟩⟨s_{1+d}⟩
    # Also returns spatial magnetization profile: mlis[i] = ⟨s_i⟩
    
    println("Processing $sample_num samples for L=$L, γ=$(γlis[inds]) with $(nprocs()) processes...")
    
    # Parallel map over all samples
    results = pmap(i -> process_single_sample_time(inds, L, i, τlis), 1:sample_num)
    
    # Filter out failed samples
    valid_results = filter(!isnothing, results)
    n_valid = length(valid_results)
    
    if n_valid == 0
        error("No valid samples found!")
    end
    
    if n_valid < sample_num
        @warn "Only $n_valid out of $sample_num samples were valid"
    end
    
    # Separate results
    corr_odd_list = [r[1] for r in valid_results]
    corr_even_list = [r[2] for r in valid_results]
    mlis_odd_list = [r[3] for r in valid_results]
    mlis_even_list = [r[4] for r in valid_results]
    
    # Stack into matrices
    corr_ensemble_odd = hcat(corr_odd_list...)'    # (n_valid, max_d)
    corr_ensemble_even = hcat(corr_even_list...)'  # (n_valid, max_d)
    mlis_ensemble_odd = hcat(mlis_odd_list...)'    # (n_valid, L)
    mlis_ensemble_even = hcat(mlis_even_list...)'  # (n_valid, L)
    
    # Average correlation over samples
    corr_avg_odd = mean(corr_ensemble_odd, dims=1)[:]
    corr_std_odd = std(corr_ensemble_odd, dims=1)[:] ./ sqrt(n_valid)
    
    corr_avg_even = mean(corr_ensemble_even, dims=1)[:]
    corr_std_even = std(corr_ensemble_even, dims=1)[:] ./ sqrt(n_valid)
    
    # Average magnetization over samples (spatial profile)
    mlis_avg_odd = mean(mlis_ensemble_odd, dims=1)[:]
    mlis_std_odd = std(mlis_ensemble_odd, dims=1)[:] ./ sqrt(n_valid)
    
    mlis_avg_even = mean(mlis_ensemble_even, dims=1)[:]
    mlis_std_even = std(mlis_ensemble_even, dims=1)[:] ./ sqrt(n_valid)
    
    println("Done! Valid samples: $n_valid")
    println("Mean magnetization: odd = $(mean(mlis_avg_odd)), even = $(mean(mlis_avg_even))")
    
    # Save with consistent key names
    τ = τlis[inds]
    D, _, _ = get_cfg_params_Born(τ, L)
    save_path = "exm/data/Bulk_measure/records_correlation/L$(L)/τ$(τ)/record_correlation_D$(div(D,L)).jld"
    mkpath(dirname(save_path))
    save(save_path,
         "corr_avg_odd", corr_avg_odd,
         "corr_std_odd", corr_std_odd,
         "corr_avg_even", corr_avg_even,
         "corr_std_even", corr_std_even,
         "mlis_avg_odd", mlis_avg_odd,
         "mlis_std_odd", mlis_std_odd,
         "mlis_avg_even", mlis_avg_even,
         "mlis_std_even", mlis_std_even)
    println("Saved to: $save_path")
    
    return corr_avg_odd, corr_std_odd, corr_avg_even, corr_std_even, mlis_avg_odd, mlis_std_odd, mlis_avg_even, mlis_std_even
end

# Command line interface
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia -p <nprocs> layer_corr.jl <gamma_index> [sample_num]")
        println("  gamma_index: 1-12, corresponding to γ = $γlis")
        println("  sample_num: number of samples (default: 10000)")
        exit(1)
    end
    
    inds = parse(Int, ARGS[1])
    sample_num = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10000
    
    if !(1 <= inds <= length(γlis))
        error("gamma_index must be between 1 and $(length(γlis))")
    end
    
    # main(inds, L_list; sample_num=sample_num)
    average_corr_time(inds, 20; sample_num=sample_num) 
end

