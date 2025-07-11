#!/bin/bash

#SBATCH -p i64m512u # u = shared nodes; ue = exclusive cores; 32 cores / node
#SBATCH -J TMI20_index${arg}
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=100
#SBATCH --exclude=cpu1-1,cpu1-99,cpu1-78,cpu1-94,cpu1-101
## SBATCH -o job-%j-%a.txt  # 使用作业数组ID来区分输出文件
#SBATCH -o job-%j.txt

##SBATCH -e job.%j.err
##SBATCH --array=1-16
## path="/hpc2hdd/home/zzhi359/quantumErgotropy/julia_shell/erg_scaling"


> resource_usage.log
(
  while true; do
    echo "$(date) - CPU: $(top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1}')% - Memory: $(free -h | grep Mem | awk '{print $3 "/" $2}')"
    sleep 60
  done
) >> resource_usage.log &

echo "Program started at: $(date)"
echo "Received parameters: index=$index, seed=$seed"
cd ..



CPU_LIMIT=$(nproc)    # 通常为44
TOTAL_MEM_GB=503      # 503 GB
PER_TASK_MEM_GB=2    # 每个任务的预计内存使用
MEM_LIMIT=$(( TOTAL_MEM_GB / PER_TASK_MEM_GB ))  # 251

# 设定最大并发数
if [[ MEM_LIMIT -lt CPU_LIMIT ]]; then
    CONCURRENCY=$MEM_LIMIT
else
    CONCURRENCY=$CPU_LIMIT
fi

for ((i=1; i<=2000; i++)); do
    # 生成一个随机种子
    RANDOM_SEED=$(( RANDOM + i * 1000 ))  # 通过任务ID来生成种子，确保不同任务之间不重复

    export JULIA_NUM_THREADS=16
    julia --project=. exm/Bulk_measure/monitored_dynamics.jl 20 $i $RANDOM_SEED

    ((task_counter++))
    echo "Submitted job $i (concurrent: $task_counter/$CONCURRENCY)"

    # 控制并发任务数量
    if (( task_counter >= CONCURRENCY )); then
        wait -n  # 等待任意一个完成的任务
        ((task_counter--))
        echo "Job completed, current concurrent: $task_counter"
    fi
done

echo "mission accomplished at: $(date)"
kill $!
