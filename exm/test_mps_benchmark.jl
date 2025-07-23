using FibonacciChain
using ITensorMPS
using Random
using LinearAlgebra
using Printf
using Statistics

"""
Simple MPS Parameter Test Script
简单的MPS参数测试脚本
"""

function test_single_configuration(N=12, τ=1.0, D=5, cutoff=1e-10, maxdim=50)
    println("Testing: N=$N, τ=$τ, D=$D, cutoff=$cutoff, maxdim=$maxdim")
    
    try
        # 生成初始MPS
        ψ, sites = initial_mps(N)
        println("✓ Initial MPS created")
        
        # 测量执行时间
        start_time = time()
        
        # 执行bulk measurement
        bulk_states, bulk_samples, bulk_weights = mps_bulk_measurement(
            ψ, sites, N, τ, D; 
            rng=MersenneTwister(100), 
            pbc=true, 
            cutoff=cutoff, 
            maxdim=maxdim
        )
        
        execution_time = time() - start_time
        println("✓ Bulk measurement completed in $(execution_time:.3f) seconds")
        
        # 计算纠缠熵
        final_state = bulk_states[end]
        ee_list = eelis_Fibo_mps(N, final_state)
        bond_dims = linkdims(final_state)
        
        # 显示结果
        println("Results:")
        @printf("  Execution time: %.3f seconds\n", execution_time)
        @printf("  Final free energy: %.6f\n", sum(bulk_weights))
        @printf("  Max entanglement entropy: %.4f\n", maximum(ee_list))
        @printf("  Mean entanglement entropy: %.4f\n", mean(ee_list))
        @printf("  Max bond dimension: %d\n", maximum(bond_dims))
        @printf("  Mean bond dimension: %.1f\n", mean(bond_dims))
        
        println("✓ Test successful")
        return true, execution_time, sum(bulk_weights), maximum(ee_list), maximum(bond_dims)
        
    catch e
        println("✗ Test failed: $e")
        return false, -1.0, NaN, NaN, -1
    end
end

function compare_parameters()
    println("=== MPS Parameter Comparison ===")
    
    # 测试参数组合
    test_configs = [
        (1e-12, 50,  "High precision, medium bond"),
        (1e-10, 50,  "Medium precision, medium bond"),
        (1e-8,  50,  "Low precision, medium bond"),
        (1e-10, 100, "Medium precision, high bond"),
        (1e-10, 20,  "Medium precision, low bond"),
    ]
    
    results = []
    
    for (i, (cutoff, maxdim, description)) in enumerate(test_configs)
        println("\n[$i/$(length(test_configs))] Testing: $description")
        println("  Parameters: cutoff=$cutoff, maxdim=$maxdim")
        
        success, time, free_energy, max_ee, max_bond = test_single_configuration(
            12, 1.0, 5, cutoff, maxdim
        )
        
        if success
            push!(results, (cutoff, maxdim, time, free_energy, max_ee, max_bond, description))
        end
    end
    
    # 分析结果
    if !isempty(results)
        println("\n=== Comparison Results ===")
        println("Cutoff\t\tMaxdim\tTime(s)\tFree Energy\tMax EE\tMax Bond\tDescription")
        println("-" ^ 80)
        
        for (cutoff, maxdim, time, fe, ee, bond, desc) in results
            @printf("%.0e\t%d\t%.3f\t%.6f\t%.4f\t%d\t\t%s\n", 
                   cutoff, maxdim, time, fe, ee, bond, desc)
        end
        
        # 找到最快的
        fastest = minimum(results, key=x -> x[3])
        println("\nFastest configuration:")
        @printf("  cutoff=%.0e, maxdim=%d, time=%.3fs (%s)\n", 
               fastest[1], fastest[2], fastest[3], fastest[7])
        
        # 找到最精确的（最高键维度）
        highest_precision = maximum(results, key=x -> x[6])
        println("\nHighest precision configuration:")
        @printf("  cutoff=%.0e, maxdim=%d, max_bond=%d (%s)\n", 
               highest_precision[1], highest_precision[2], highest_precision[6], highest_precision[7])
    end
    
    return results
end

function detailed_layer_analysis(cutoff=1e-10, maxdim=50)
    println("=== Detailed Layer-by-Layer Analysis ===")
    println("Parameters: cutoff=$cutoff, maxdim=$maxdim")
    
    N, τ, D = 12, 1.0, 8
    ψ, sites = initial_mps(N)
    
    println("\nRunning measurement with layer tracking...")
    
    bulk_states, bulk_samples, bulk_weights = mps_bulk_measurement(
        ψ, sites, N, τ, D; 
        rng=MersenneTwister(100), 
        pbc=true, 
        cutoff=cutoff, 
        maxdim=maxdim
    )
    
    println("\nLayer-by-layer results:")
    println("Layer\tFree Energy\tMax EE\t\tMean EE\t\tMax Bond\tMean Bond")
    println("-" ^ 70)
    
    cumulative_fe = 0.0
    for layer in 1:D
        state = bulk_states[layer]
        ee_list = eelis_Fibo_mps(N, state)
        bond_dims = linkdims(state)
        
        cumulative_fe += bulk_weights[layer]
        
        @printf("%d\t%.6f\t%.4f\t\t%.4f\t\t%d\t\t%.1f\n",
               layer, cumulative_fe, maximum(ee_list), mean(ee_list), 
               maximum(bond_dims), mean(bond_dims))
    end
    
    return bulk_states, bulk_samples, bulk_weights
end

function convergence_test()
    println("=== Convergence Test ===")
    println("Testing how results change with different maxdim values")
    
    N, τ, D = 12, 1.0, 5
    cutoff = 1e-10
    maxdims = [20, 30, 50, 100, 200]
    
    results = []
    
    for maxdim in maxdims
        println("Testing maxdim=$maxdim...")
        
        success, time, free_energy, max_ee, max_bond = test_single_configuration(
            N, τ, D, cutoff, maxdim
        )
        
        if success
            push!(results, (maxdim, time, free_energy, max_ee, max_bond))
        end
    end
    
    if length(results) > 1
        println("\nConvergence analysis:")
        println("Maxdim\tTime(s)\tFree Energy\tΔFE\t\tMax EE\t\tΔEE\t\tMax Bond")
        println("-" ^ 75)
        
        for (i, (maxdim, time, fe, ee, bond)) in enumerate(results)
            if i == 1
                @printf("%d\t%.3f\t%.6f\t--\t\t%.4f\t\t--\t\t%d\n", 
                       maxdim, time, fe, ee, bond)
            else
                prev_fe = results[i-1][3]
                prev_ee = results[i-1][4]
                delta_fe = abs(fe - prev_fe)
                delta_ee = abs(ee - prev_ee)
                
                @printf("%d\t%.3f\t%.6f\t%.6f\t%.4f\t\t%.6f\t%d\n", 
                       maxdim, time, fe, delta_fe, ee, delta_ee, bond)
            end
        end
    end
    
    return results
end

# 主函数
function main()
    println("MPS Parameter Optimization Test Suite")
    println("=====================================")
    
    println("\n1. Single configuration test:")
    test_single_configuration()
    
    println("\n2. Parameter comparison:")
    compare_parameters()
    
    println("\n3. Detailed layer analysis:")
    detailed_layer_analysis()
    
    println("\n4. Convergence test:")
    convergence_test()
    
    println("\nAll tests completed!")
end

# 如果直接运行此文件
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
