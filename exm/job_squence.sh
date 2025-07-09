#!/bin/bash

# 系统参数
CPU_LIMIT=$(nproc)    # 通常为44
TOTAL_MEM_GB=503      # 503 GB
PER_TASK_MEM_GB=2     # 每个任务的预计内存使用
MEM_LIMIT=$(( TOTAL_MEM_GB / PER_TASK_MEM_GB ))  # 251

# 设定最大并发数
if [[ MEM_LIMIT -lt CPU_LIMIT ]]; then
    CONCURRENCY=$MEM_LIMIT
else
    CONCURRENCY=$CPU_LIMIT
fi

# 安全起见，保留核心供系统使用
CONCURRENCY=$(( CONCURRENCY > 1 ? CONCURRENCY - 1 : 1 ))

echo "Detected CPU cores: $CPU_LIMIT, Total RAM: ${TOTAL_MEM_GB}GB"
echo "Setting optimal concurrency to: $CONCURRENCY"

task_counter=0

for ((i=1; i<=2000; i++)); do
    # 生成一个随机种子
    RANDOM_SEED=$(( RANDOM + i * 1000 ))  # 通过任务ID来生成种子，确保不同任务之间不重复

    # 提交惰性任务，执行当前任务
    nohup julia --project=. exm/Bulk_measure/monitored_dynamics.jl "$i" "$RANDOM_SEED" &

    ((task_counter++))
    echo "Submitted job $i (concurrent: $task_counter/$CONCURRENCY)"

    # 控制并发任务数量
    if (( task_counter >= CONCURRENCY )); then
        wait -n  # 等待任意一个完成的任务
        ((task_counter--))
        echo "Job completed, current concurrent: $task_counter"
    fi
done

wait  # 等待剩余任务
echo "All jobs completed"