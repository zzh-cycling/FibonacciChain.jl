#!/bin/bash

# 测试脚本 - 验证并发控制是否正常工作

echo "=== Testing Job Control System ==="

# 设置小的参数进行测试
CONCURRENCY=3
total_jobs=0
completed_jobs=0
failed_jobs=0
task_counter=0

declare -a job_pids=()

# 调试函数：显示当前进程状态
debug_show_status() {
    echo "=== Debug Status ==="
    echo "Task counter: $task_counter"
    echo "Active PIDs (${#job_pids[@]}): ${job_pids[*]}"
    echo "Actual running processes:"
    
    local running_count=0
    for pid in "${job_pids[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
            echo "  PID $pid: RUNNING"
            ((running_count++))
        else
            echo "  PID $pid: NOT RUNNING"
        fi
    done
    
    echo "Running count: $running_count"
    echo "Total jobs: $total_jobs, Completed: $completed_jobs, Failed: $failed_jobs"
    echo "==================="
}

# 函数：等待一个作业完成并更新计数
wait_for_job() {
    if [[ ${#job_pids[@]} -gt 0 ]]; then
        # 检查是否有作业完成
        local new_job_pids=()
        local jobs_found=false
        
        for i in "${!job_pids[@]}"; do
            pid=${job_pids[i]}
            if ! kill -0 "$pid" 2>/dev/null; then
                # 进程已完成，检查退出状态（非阻塞）
                wait "$pid" 2>/dev/null
                exit_code=$?
                if [[ $exit_code -eq 0 ]]; then
                    ((completed_jobs++))
                    echo "✓ Job PID $pid completed successfully (Total completed: $completed_jobs)"
                else
                    ((failed_jobs++))
                    echo "✗ Job PID $pid failed with exit code $exit_code (Total failed: $failed_jobs)"
                fi
                ((task_counter--))
                jobs_found=true
            else
                # 进程仍在运行，保留在数组中
                new_job_pids+=("$pid")
            fi
        done
        
        # 更新PID数组
        job_pids=("${new_job_pids[@]}")
        
        if [[ "$jobs_found" == "true" ]]; then
            return 0
        else
            # 如果没有作业完成，短暂等待
            sleep 0.1
            return 1
        fi
    fi
    return 1
}

# 函数：等待直到有空闲槽位
wait_for_slot() {
    while (( task_counter >= CONCURRENCY )); do
        echo "Waiting for slot... (current: $task_counter/$CONCURRENCY)"
        if ! wait_for_job; then
            sleep 0.5
        fi
        # 额外检查：清理可能的僵尸进程
        jobs -p >/dev/null 2>&1
    done
}

# 函数：提交测试作业
submit_test_job() {
    local job_id=$1
    local duration=$2
    
    # 创建一个简单的测试脚本
    local test_script="test_job_${job_id}.sh"
    cat > "$test_script" << EOF
#!/bin/bash
echo "Test job $job_id started at \$(date)"
echo "PID: \$\$"
echo "Sleeping for $duration seconds..."
sleep $duration
echo "Test job $job_id completed at \$(date)"
exit 0
EOF
    
    chmod +x "$test_script"
    
    # 提交作业到后台
    nohup ./"$test_script" > "test_job_${job_id}.log" 2>&1 &
    local pid=$!
    
    # 验证进程是否真正启动
    sleep 0.1
    if ! kill -0 "$pid" 2>/dev/null; then
        echo "✗ Failed to start test job $job_id (PID: $pid)"
        rm -f "$test_script"
        return 1
    fi
    
    # 存储PID并更新计数
    job_pids+=("$pid")
    ((task_counter++))
    ((total_jobs++))
    
    echo "🚀 Submitted test job $job_id (PID: $pid, concurrent: $task_counter/$CONCURRENCY)"
    return 0
}

echo "Starting test with CONCURRENCY=$CONCURRENCY"
echo "Will submit 10 test jobs with varying durations"

# 提交测试作业
for i in {1..10}; do
    # 等待空闲槽位
    wait_for_slot
    
    # 随机持续时间 2-6秒
    duration=$((2 + RANDOM % 5))
    
    # 提交作业
    if submit_test_job "$i" "$duration"; then
        echo "Test job $i submitted (duration: ${duration}s)"
        debug_show_status
    else
        echo "⚠️  Failed to submit test job $i"
    fi
    
    sleep 0.5
done

echo "All test jobs submitted. Waiting for completion..."

# 等待所有剩余作业完成
while [[ ${#job_pids[@]} -gt 0 ]]; do
    echo "Waiting for ${#job_pids[@]} remaining jobs..."
    wait_for_job || sleep 1
    debug_show_status
done

# 清理测试文件
rm -f test_job_*.sh test_job_*.log

# 最终统计
echo ""
echo "========== Test Results =========="
echo "Total jobs submitted: $total_jobs"
echo "Successfully completed: $completed_jobs"
echo "Failed jobs: $failed_jobs"

if [[ $total_jobs -gt 0 ]]; then
    success_rate=$(( completed_jobs * 100 / total_jobs ))
    echo "Success rate: ${success_rate}%"
fi

if [[ $completed_jobs -eq $total_jobs && $failed_jobs -eq 0 ]]; then
    echo "✅ Test PASSED - All jobs completed successfully!"
else
    echo "❌ Test FAILED - Some jobs did not complete properly"
fi

echo "Test completed!"
