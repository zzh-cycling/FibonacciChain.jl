using FibonacciChain
using ITensorMPS
using Random
using BenchmarkTools
using LinearAlgebra
using Printf
using Statistics

"""
Example usage of MPS-based Fibonacci chain measurements
"""
struct BenchmarkResult
    cutoff::Float64
    maxdim::Int
    execution_time::Float64
    memory_usage::Float64
    free_energy::Vector{Float64}
    final_entanglement_entropy::Vector{Float64}
    bond_dimensions::Vector{Int}
    convergence_error::Float64
    # 与精确对角化结果的比较
    free_energy_error::Float64
    ee_error::Float64
    success::Bool
    error_message::String
end

function mps_example(cutoff=1e-10, maxdim=50)

    # System parameters
    N = 20  # Number of sites
    τ = atanh(0.6)  # Measurement parameter
    pbc = true  # Periodic boundary conditions
    
    ψ, sites = initial_mps(N)
    
    # Perform bulk measurements with alternating pattern
    D = 45N  # Number of layers
    
    elapsed_time = @elapsed bulk_states, bulk_samples, bulk_weights = mps_bulk_measurement(
        ψ, sites, N, τ, D; rng=MersenneTwister(100), pbc=pbc, cutoff=cutoff, maxdim=maxdim)
    
    EElis = eelis_Fibo_mps(N, bulk_states[end])
    println("Performed $D layers of bulk measurements")
    for layer in 1:D
        println("Layer $layer:")
        println("  Measurement pattern: $([bulk_samples[layer, :]])")
        println("  Log probability: $(bulk_weights[layer])")
        println("  Final state bond dims: $(linkdims(bulk_states[layer]))")
        println("  Final state entanglement entropy: $(EElis)")
    end
    
    println("\n=== Example completed ===")

    return bulk_weights, EElis, elapsed_time
end

function get_exact_results(N, τ, D; rng_seed=100, pbc=true)
    rng = MersenneTwister(rng_seed)
    initial_st = zeros(length(Fibonacci_basis(N)))
    initial_st[1] = 1.0
    
    bulk_states, bulk_samples, bulk_weights = Bulkmeasure(N, τ, initial_st, D, rng, pbc)
    
    exact_ee = eelis_Fibo_state(N, bulk_states[end])
    
    return bulk_weights, exact_ee, bulk_states, bulk_samples
end
"""
    run_mps_benchmark(N, τ, D, cutoff, maxdim; rng_seed=100, pbc=true, exact_results=nothing)

运行单个MPS参数配置的基准测试
"""
function run_mps_benchmark(N, τ, D, cutoff, maxdim; rng_seed=100, pbc=true, exact_results=nothing)
    try
        # 设置随机种子确保可重复性
        rng = MersenneTwister(rng_seed)
        
        # 如果没有提供精确结果，计算一次
        if exact_results === nothing
            exact_free_energy, exact_ee, _, _, _ = get_exact_results(N, τ, D; rng_seed=rng_seed, pbc=pbc)
        else
            exact_free_energy, exact_ee = exact_results
        end
        
        # 生成初始MPS
        ψ, sites = initial_mps(N)
        
        # 测量执行时间和内存使用
        result = @timed begin
            bulk_states, bulk_samples, bulk_weights = mps_bulk_measurement(
                ψ, sites, N, τ, D; rng=rng, pbc=pbc, cutoff=cutoff, maxdim=maxdim)
            
            # 计算最终的纠缠熵
            final_state = bulk_states[end]
            ee_list = eelis_Fibo_mps(N, final_state)
            
            # 计算最终的自由能
            final_free_energy = bulk_weights
            
            # 获取键维度
            bond_dims = linkdims(final_state)
            
            (final_state, ee_list, final_free_energy, bond_dims)
        end
        
        final_state, ee_list, final_free_energy, bond_dims = result.value
        execution_time = result.time
        memory_usage = result.bytes / 1024^2  # MB
        
        # 计算收敛误差（基于键维度的变化）
        max_bond_dim = maximum(bond_dims)
        convergence_error = max_bond_dim >= maxdim ? 1.0 : 0.1 * max_bond_dim / maxdim
        
        # 计算与精确结果的误差
        free_energy_error = sqrt(mean(final_free_energy .- exact_free_energy).^2)
        ee_error = sqrt(mean((ee_list .- exact_ee).^2))  # RMSE
        
        return BenchmarkResult(
            cutoff, maxdim, execution_time, memory_usage,
            final_free_energy, ee_list, bond_dims, convergence_error,
            free_energy_error, ee_error,
            true, ""
        )
        
    catch e
        return BenchmarkResult(
            cutoff, maxdim, -1.0, -1.0, NaN, Float64[], Int[], 1.0,
            NaN, NaN,
            false, string(e)
        )
    end
end

"""
    parameter_sweep(N=20, τ=1.0, D=10; rng_seed=100, pbc=true)

对不同的cutoff和maxdim参数进行扫描
"""
function parameter_sweep(N=20, τ=1.0, D=10; rng_seed=100, pbc=true)
    # 定义参数范围
    cutoffs = [1e-12, 1e-10, 1e-10, 1e-8, 1e-6]
    maxdims = [50,60, 100,150,  200]
    
    println("=== MPS Parameter Sweep ===")
    println("System: N=$N, τ=$τ, D=$D")
    println("Cutoffs: $cutoffs")
    println("Max dimensions: $maxdims")
    println("Total combinations: $(length(cutoffs) * length(maxdims))")
    println()
    
    # 首先计算精确对角化结果作为基准
    println("Computing exact diagonalization benchmark...")
    exact_free_energy, exact_ee, _, _, _ = get_exact_results(N, τ, D; rng_seed=rng_seed, pbc=pbc)
    @printf("Exact results: Final Free Energy = %.6f, EE = %.4f ± %.4f\n", 
           exact_free_energy[end], mean(exact_ee), std(exact_ee))
    println()
    
    results = BenchmarkResult[]
    total_combinations = length(cutoffs) * length(maxdims)
    current_combo = 0
    
    for cutoff in cutoffs
        for maxdim in maxdims
            current_combo += 1
            print("[$current_combo/$total_combinations] Testing cutoff=$cutoff, maxdim=$maxdim... ")
            
            result = run_mps_benchmark(N, τ, D, cutoff, maxdim; 
                                     rng_seed=rng_seed, pbc=pbc, 
                                     exact_results=(exact_free_energy, exact_ee))
            push!(results, result)
            
            if result.success
                @printf("✓ Time: %.3fs, ΔF: %.2e, ΔEE: %.2e\n", 
                       result.execution_time, result.free_energy_error, result.ee_error)
            else
                println("✗ Failed: $(result.error_message)")
            end
        end
    end
    
    return results, (exact_free_energy, exact_ee)
end

function find_pareto_optimal_accuracy(results::Vector{BenchmarkResult})
    pareto_results = BenchmarkResult[]
    
    for result in results
        is_dominated = false
        for other in results
            if (other.execution_time <= result.execution_time && 
                other.free_energy_error <= result.free_energy_error &&
                (other.execution_time < result.execution_time || 
                 other.free_energy_error < result.free_energy_error))
                is_dominated = true
                break
            end
        end
        
        if !is_dominated
            push!(pareto_results, result)
        end
    end
    
    # 按执行时间排序
    sort!(pareto_results, by=r -> r.execution_time)
    return pareto_results
end

"""
    analyze_results(results::Vector{BenchmarkResult}, exact_results=nothing)

分析基准测试结果，包括与精确解的比较
"""
function analyze_results(results::Vector{BenchmarkResult}, exact_results=nothing)
    println("\n=== Benchmark Analysis ===")
    
    # 过滤成功的结果
    successful_results = filter(r -> r.success, results)
    
    if isempty(successful_results)
        println("No successful results to analyze!")
        return
    end
    
    println("Successful runs: $(length(successful_results))/$(length(results))")
    println()
    
    if exact_results !== nothing
        exact_free_energy, exact_ee = exact_results
        @printf("Exact results: Free Energy = %.6f, EE = %.4f ± %.4f\n", 
               exact_free_energy[end], mean(exact_ee), std(exact_ee))
        println()
    end
    
    # 按不同标准排序
    precision_sorted = sort(successful_results, by=r -> r.free_energy_error)  # 按自由能误差
    ee_precision_sorted = sort(successful_results, by=r -> r.ee_error)        # 按纠缠熵误差
    speed_sorted = sort(successful_results, by=r -> r.execution_time)
    memory_sorted = sort(successful_results, by=r -> r.memory_usage)
    
    println("=== Top 5 Most Accurate (Free Energy) ===")
    for (i, result) in enumerate(precision_sorted[1:min(5, end)])
        @printf("%d. cutoff=%.0e, maxdim=%d: ΔF=%.2e, time=%.3fs, memory=%.1fMB\n",
               i, result.cutoff, result.maxdim, result.free_energy_error, 
               result.execution_time, result.memory_usage)
    end
    
    println("\n=== Top 5 Most Accurate (Entanglement Entropy) ===")
    for (i, result) in enumerate(ee_precision_sorted[1:min(5, end)])
        @printf("%d. cutoff=%.0e, maxdim=%d: ΔEE=%.2e, time=%.3fs, memory=%.1fMB\n",
               i, result.cutoff, result.maxdim, result.ee_error,
               result.execution_time, result.memory_usage)
    end
    
    println("\n=== Top 5 Fastest ===")
    for (i, result) in enumerate(speed_sorted[1:min(5, end)])
        @printf("%d. cutoff=%.0e, maxdim=%d: time=%.3fs, ΔF=%.2e, ΔEE=%.2e\n",
               i, result.cutoff, result.maxdim, result.execution_time,
               result.free_energy_error, result.ee_error)
    end
    
    println("\n=== Top 5 Memory Efficient ===")
    for (i, result) in enumerate(memory_sorted[1:min(5, end)])
        @printf("%d. cutoff=%.0e, maxdim=%d: memory=%.1fMB, time=%.3fs, ΔF=%.2e\n",
               i, result.cutoff, result.maxdim, result.memory_usage,
               result.execution_time, result.free_energy_error)
    end
    
    # 计算Pareto前沿（时间vs精度）
    println("\n=== Pareto Optimal (Time vs Free Energy Accuracy) ===")
    pareto_results = find_pareto_optimal_accuracy(successful_results)
    for (i, result) in enumerate(pareto_results)
        @printf("%d. cutoff=%.0e, maxdim=%d: time=%.3fs, ΔF=%.2e, ΔEE=%.2e\n",
               i, result.cutoff, result.maxdim, result.execution_time,
               result.free_energy_error, result.ee_error)
    end
    
    # 推荐参数
    println("\n=== Parameter Recommendations ===")
    
    # 设定精度阈值
    fe_threshold = 1e-4  # 自由能误差阈值
    ee_threshold = 1e-4  # 纠缠熵误差阈值
    
    # 找到满足精度要求的结果
    accurate_results = filter(r -> r.free_energy_error < fe_threshold && r.ee_error < ee_threshold, successful_results)
    
    if !isempty(accurate_results)
        # 在满足精度要求的结果中找最快的
        fastest_accurate = argmin(r -> r.execution_time, accurate_results)
        @printf("Best for production (ΔF<%.0e, ΔEE<%.0e): cutoff=%.0e, maxdim=%d (time=%.3fs)\n",
               fe_threshold, ee_threshold, fastest_accurate.cutoff, fastest_accurate.maxdim, 
               fastest_accurate.execution_time)
    else
        println("No results meet the accuracy threshold. Consider looser requirements or higher maxdim.")
    end
    
    # 最佳平衡（时间和精度的加权几何平均）
    best_balance = argmin(r ->sqrt(r.execution_time * (r.free_energy_error + r.ee_error)), successful_results)
    @printf("Best overall balance: cutoff=%.0e, maxdim=%d (time=%.3fs, ΔF=%.2e, ΔEE=%.2e)\n",
           best_balance.cutoff, best_balance.maxdim, best_balance.execution_time,
           best_balance.free_energy_error, best_balance.ee_error)
    
    # 快速原型设计
    fast_prototype = speed_sorted[1]
    @printf("For fast prototyping: cutoff=%.0e, maxdim=%d (time=%.3fs, ΔF=%.2e)\n",
           fast_prototype.cutoff, fast_prototype.maxdim, fast_prototype.execution_time,
           fast_prototype.free_energy_error)
    
    # 高精度计算
    high_precision = precision_sorted[1]
    @printf("For high precision: cutoff=%.0e, maxdim=%d (ΔF=%.2e, ΔEE=%.2e)\n",
           high_precision.cutoff, high_precision.maxdim, 
           high_precision.free_energy_error, high_precision.ee_error)
    
    return successful_results, pareto_results, best_balance
end

"""
    multi_seed_benchmark(N=20, τ=1.0, D=10, seeds=1:10, cutoffs=[1e-12, 1e-10, 1e-8], maxdims=[50, 100, 200])

对多个随机种子进行MPS参数benchmark，对结果做平均统计
"""
function multi_seed_benchmark(N=20, τ=1.0, D=10, seeds=1:10, cutoffs=[1e-12, 1e-10, 1e-8], maxdims=[50, 100, 200])
    println("=== Multi-Seed MPS Benchmark ===")
    println("System: N=$N, τ=$τ, D=$D")
    println("Seeds: $seeds ($(length(seeds)) seeds)")
    println("Cutoffs: $cutoffs")
    println("Max dimensions: $maxdims")
    println()
    
    # 存储每个参数组合的多seed结果
    param_results = Dict()
    
    # 计算每个seed的精确结果
    exact_results_per_seed = []
    for seed in seeds
        exact_fe, exact_ee, _, _ = get_exact_results(N, τ, D; rng_seed=seed, pbc=true)
        push!(exact_results_per_seed, (exact_fe[end], exact_ee))
    end
    
    total_runs = length(cutoffs) * length(maxdims) * length(seeds)
    current_run = 0
    
    for cutoff in cutoffs
        for maxdim in maxdims
            param_key = (cutoff, maxdim)
            param_results[param_key] = []
            
            for (i, seed) in enumerate(seeds)
                current_run += 1
                print("[$current_run/$total_runs] Testing cutoff=$cutoff, maxdim=$maxdim, seed=$seed... ")
                
                exact_fe, exact_ee = exact_results_per_seed[i]
                result = run_mps_benchmark(N, τ, D, cutoff, maxdim; 
                                         rng_seed=seed, pbc=true, 
                                         exact_results=(exact_fe, exact_ee))
                
                push!(param_results[param_key], result)
                
                if result.success
                    @printf("✓ Time: %.3fs, ΔF: %.2e, ΔEE: %.2e\n", result.execution_time, result.free_energy_error, result.ee_error)
                else
                    println("✗ Failed")
                end
            end
        end
    end
    
    # 分析多seed结果
    println("\n=== Multi-Seed Analysis ===")
    println("Parameter combination results (mean ± std):")
    println("Format: cutoff, maxdim → time(s), ΔF_error, ΔEE_error, success_rate")
    println()
    
    summary_results = []
    
    for cutoff in cutoffs
        for maxdim in maxdims
            param_key = (cutoff, maxdim)
            seed_results = param_results[param_key]
            successful = filter(r -> r.success, seed_results)
            
            if !isempty(successful)
                times = [r.execution_time for r in successful]
                fe_errors = [r.free_energy_error for r in successful]
                ee_errors = [r.ee_error for r in successful]
                
                mean_time = mean(times)
                std_time = std(times)
                mean_fe_error = mean(fe_errors)
                std_fe_error = std(fe_errors)
                mean_ee_error = mean(ee_errors)
                std_ee_error = std(ee_errors)
                success_rate = length(successful) / length(seed_results)
                
                @printf("%.0e, %3d → %.3f±%.3f, %.2e±%.2e, %.2e±%.2e, %.1f%%\n",
                       cutoff, maxdim, mean_time, std_time, 
                       mean_fe_error, std_fe_error, mean_ee_error, std_ee_error, 
                       success_rate * 100)
                
                push!(summary_results, (cutoff, maxdim, mean_time, mean_fe_error, mean_ee_error, success_rate))
            else
                @printf("%.0e, %3d → All failed\n", cutoff, maxdim)
            end
        end
    end
    
    # 推荐最优参数
    println("\n=== Recommendations ===")
    if !isempty(summary_results)
        # 过滤成功率100%的结果
        reliable_results = filter(r -> r[6] == 1.0, summary_results)
        
        if !isempty(reliable_results)
            # 最佳精度
            best_accuracy = argmin(r -> r[4] + r[5], reliable_results)  # 最小误差和
            @printf("Best accuracy: cutoff=%.0e, maxdim=%d (ΔF=%.2e, ΔEE=%.2e)\n",
                   best_accuracy[1], best_accuracy[2], best_accuracy[4], best_accuracy[5])
            
            # 最佳速度
            best_speed = argmin(r -> r[3], reliable_results)  # 最短时间
            @printf("Best speed: cutoff=%.0e, maxdim=%d (time=%.3fs)\n",
                   best_speed[1], best_speed[2], best_speed[3])
            
            # 最佳平衡
            best_balance = argmin(r -> sqrt(r[3] * (r[4] + r[5])), reliable_results)
            @printf("Best balance: cutoff=%.0e, maxdim=%d (time=%.3fs, errors=%.2e)\n",
                   best_balance[1], best_balance[2], best_balance[3], best_balance[4] + best_balance[5])
        else
            println("No parameters achieved 100% success rate across all seeds")
        end
    end
    
    return param_results, summary_results
end

N=20;
τ=1.0;
D=10;
rng_seed=100;
results, (exact_free_energy, exact_ee) = parameter_sweep(N, τ, D; rng_seed=rng_seed)
successful_results, pareto_results, best_balance = analyze_results(results, (exact_free_energy, exact_ee));


param_results, summary = multi_seed_benchmark(20, 1.0, 100, 1:10, [1e-12, 1e-10, 1e-10, 1e-8, 1e-6], [50,60, 100,150, 200])