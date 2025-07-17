#!/bin/bash

# 检测操作系统
OS=$(uname -s)

# 系统参数检测（兼容macOS和Linux）
if [[ "$OS" == "Darwin" ]]; then
    # macOS
    CPU_LIMIT=$(sysctl -n hw.ncpu)
    TOTAL_MEM_BYTES=$(sysctl -n hw.memsize)
    TOTAL_MEM_GB=$(( TOTAL_MEM_BYTES / 1024 / 1024 / 1024 ))
else
    # Linux
    CPU_LIMIT=$(nproc)
    TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    TOTAL_MEM_GB=$(( TOTAL_MEM_KB / 1024 / 1024 ))
fi

# 每个任务的预计内存使用（GB）
PER_TASK_MEM_GB=2

# 基于内存计算最大并发数
MEM_LIMIT=$(( TOTAL_MEM_GB / PER_TASK_MEM_GB ))

# 设定最大并发数（取CPU和内存限制的较小值）
if [[ $MEM_LIMIT -lt $CPU_LIMIT ]]; then
    CONCURRENCY=$MEM_LIMIT
else
    CONCURRENCY=$CPU_LIMIT
fi

# 安全起见，保留一些资源供系统使用
CONCURRENCY=$(( CONCURRENCY > 4 ? CONCURRENCY - 2 : CONCURRENCY / 2 ))
CONCURRENCY=$(( CONCURRENCY > 0 ? CONCURRENCY : 1 ))

echo "Operating System: $OS"
echo "Detected CPU cores: $CPU_LIMIT, Total RAM: ${TOTAL_MEM_GB}GB"
echo "Memory-based limit: $MEM_LIMIT tasks"
echo "Setting optimal concurrency to: $CONCURRENCY"

# 作业控制变量
task_counter=0
total_jobs=0
completed_jobs=0
failed_jobs=0

# 存储后台进程PID的数组
declare -a job_pids=()

# 调试函数：显示当前进程状态
debug_show_status() {
    echo "=== Status Update ==="
    echo "Task counter: $task_counter/$CONCURRENCY"
    echo "Total: $total_jobs | Completed: $completed_jobs | Failed: $failed_jobs | Running: $task_counter"
    
    if [[ ${#job_pids[@]} -gt 0 ]]; then
        echo "Active PIDs: ${job_pids[*]}"
        
        # 验证所有PID的状态
        local running_count=0
        for pid in "${job_pids[@]}"; do
            if kill -0 "$pid" 2>/dev/null; then
                ((running_count++))
            fi
        done
        
        if [[ $running_count -ne $task_counter ]]; then
            echo "⚠️  Mismatch: task_counter=$task_counter but actual running=$running_count"
        fi
    else
        echo "No active jobs"
    fi
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
        if ! wait_for_job; then
            sleep 0.5
        fi
        # 额外检查：清理可能的僵尸进程
        jobs -p >/dev/null 2>&1
    done
}

# 函数：提交新作业
submit_job() {
    local j=$1
    local i=$2
    local seed=$3
    
    # 提交作业到后台
    nohup julia --project=. exm/Bulk_measure/monitored_dynamics_mps.jl "$j" "$i" "$seed" > "logs/job_${j}_${i}.log" 2>&1 &
    local pid=$!
    
    # 验证进程是否真正启动
    sleep 0.1
    if ! kill -0 "$pid" 2>/dev/null; then
        echo "✗ Failed to start job $i for j=$j (PID: $pid)"
        ((failed_jobs++))
        return 1
    fi
    
    # 存储PID并更新计数
    job_pids+=("$pid")
    ((task_counter++))
    ((total_jobs++))
    
    echo "🚀 Submitted job $i for j=$j (PID: $pid, concurrent: $task_counter/$CONCURRENCY)"
    return 0
}

# 创建日志目录
mkdir -p logs

J_START=32
J_END=32
J_STEP=2
I_START=101
I_END=300

echo "=== SLURM Job Controller Started at $(date) ==="
echo "Starting job submission..."
echo "Target range: j=${J_START}..${J_END} (step ${J_STEP}), i=${I_START}..${I_END}"

# 启动后台状态监控
(
    while true; do
        sleep 30  # 每30秒显示一次状态
        echo ""
        echo "=== Periodic Status Check $(date) ==="
        echo "Total: $total_jobs | Completed: $completed_jobs | Failed: $failed_jobs | Running: $task_counter/$CONCURRENCY"
        if [[ ${#job_pids[@]} -gt 0 ]]; then
            echo "Active jobs: ${#job_pids[@]} PIDs"
        fi
        echo "==============================================="
    done
) &
monitor_pid=$!

for ((j=J_START; j<=J_END; j+=J_STEP)); do
    for ((i=I_START; i<=I_END; i++)); do
        # 等待空闲槽位
        wait_for_slot
        
        # 生成随机种子
        RANDOM_SEED=$(( (j + i) * 1000 ))
        
        # 提交作业
        if submit_job "$j" "$i" "$RANDOM_SEED"; then
            # 显示进度
            echo "Progress: $total_jobs/$((I_END - I_START + 1)) submitted, $completed_jobs completed, $failed_jobs failed, $task_counter running"
            
            # 每10个作业显示一次详细状态
            if (( total_jobs % 10 == 0 )); then
                debug_show_status
            fi
        else
            echo "⚠️  Failed to submit job for j=$j, i=$i"
        fi
        
        # 短暂延迟避免过快提交
        sleep 0.1
    done
done

echo "All jobs submitted. Waiting for completion..."

# 等待所有剩余作业完成
while [[ ${#job_pids[@]} -gt 0 ]]; do
    wait_for_job || sleep 1
    if (( ${#job_pids[@]} > 0 )); then
        echo "Waiting for ${#job_pids[@]} remaining jobs..."
    fi
done

# 停止后台监控
kill $monitor_pid 2>/dev/null

# 最终统计
echo ""
echo "========== Job Execution Summary =========="
echo "Total jobs submitted: $total_jobs"
echo "Successfully completed: $completed_jobs"
echo "Failed jobs: $failed_jobs"
if [[ $total_jobs -gt 0 ]]; then
    echo "Success rate: $(( completed_jobs * 100 / total_jobs ))%"
fi

if [[ $failed_jobs -gt 0 ]]; then
    echo ""
    echo "⚠️  Some jobs failed. Check log files in logs/ directory for details."
    echo "Failed job logs can be found at: logs/job_*.log"
fi

echo "All jobs completed!"

# 可选：显示系统资源使用情况
if command -v top >/dev/null 2>&1; then
    echo ""
    echo "Current system load:"
    if [[ "$OS" == "Darwin" ]]; then
        # macOS
        echo "Load average: $(uptime | awk -F'load averages:' '{ print $2 }')"
        if command -v memory_pressure >/dev/null 2>&1; then
            echo "Memory pressure: $(memory_pressure | head -n 1)"
        fi
    else
        # Linux
        echo "Load average: $(uptime | awk -F'load average:' '{ print $2 }')"
        echo "Memory usage: $(free -h | grep ^Mem)"
    fi
fi