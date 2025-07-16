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

# 函数：等待一个作业完成并更新计数
wait_for_job() {
    if [[ ${#job_pids[@]} -gt 0 ]]; then
        # 检查是否有作业完成
        for i in "${!job_pids[@]}"; do
            pid=${job_pids[i]}
            if ! kill -0 "$pid" 2>/dev/null; then
                # 进程已完成，检查退出状态
                wait "$pid"
                exit_code=$?
                if [[ $exit_code -eq 0 ]]; then
                    ((completed_jobs++))
                    echo "✓ Job PID $pid completed successfully (Total completed: $completed_jobs)"
                else
                    ((failed_jobs++))
                    echo "✗ Job PID $pid failed with exit code $exit_code (Total failed: $failed_jobs)"
                fi
                
                # 从数组中移除这个PID
                unset job_pids[i]
                job_pids=("${job_pids[@]}")  # 重新索引数组
                ((task_counter--))
                return 0
            fi
        done
        
        # 如果没有作业完成，短暂等待
        sleep 0.1
        return 1
    fi
    return 1
}

# 函数：等待直到有空闲槽位
wait_for_slot() {
    while (( task_counter >= CONCURRENCY )); do
        wait_for_job || sleep 0.5
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
    
    # 存储PID并更新计数
    job_pids+=("$pid")
    ((task_counter++))
    ((total_jobs++))
    
    echo "🚀 Submitted job $i for j=$j (PID: $pid, concurrent: $task_counter/$CONCURRENCY)"
}

# 创建日志目录
mkdir -p logs

echo "Starting job submission..."
echo "Target jobs: j=32, i=1..50 (Total: 50 jobs)"

for ((j=32; j<=32; j+=2)); do
    for ((i=61; i<=300; i++)); do
        # 等待空闲槽位
        wait_for_slot
        
        # 生成随机种子
        RANDOM_SEED=$(( (j + i) * 1000 ))
        
        # 提交作业
        submit_job "$j" "$i" "$RANDOM_SEED"
        
        # 显示进度
        echo "Progress: $total_jobs/50 submitted, $completed_jobs completed, $failed_jobs failed"
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

# 最终统计
echo ""
echo "========== Job Execution Summary =========="
echo "Total jobs submitted: $total_jobs"
echo "Successfully completed: $completed_jobs"
echo "Failed jobs: $failed_jobs"
echo "Success rate: $(( completed_jobs * 100 / total_jobs ))%"

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