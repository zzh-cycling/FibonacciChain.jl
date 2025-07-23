#!/bin/sh
#BSUB -J gamma0.1
#BSUB -q normal
#BSUB -o output_%J.txt
#BSUB -e err_%J.txt

# 输出节点的CPU和内存信息
echo "Gathering node resource information..."
CPU_COUNT=$(nproc)  # 获取CPU核心数
TOTAL_MEM=$(free -h | grep Mem | awk '{print $2}')  # 总内存

echo "CPU Cores Available: $CPU_COUNT"
echo "Total Memory: $TOTAL_MEM"

# 假设每个任务需要1个CPU核心和2GB内存，你可以根据实际需要调整
CPU_PER_JOB=1
MEM_PER_JOB=2  # GB

# 将总内存转换为GB单位以进行计算
TOTAL_MEM_GB=$(free -g | grep Mem | awk '{print $2}')

# 计算最多可以并行的任务数量
MAX_CONCURRENT_JOBS_BY_CPU=$((CPU_COUNT / CPU_PER_JOB))
MAX_CONCURRENT_JOBS_BY_MEM=$((TOTAL_MEM_GB / MEM_PER_JOB))

# 选择更小的可用插槽数作为限制
MAX_CONCURRENT_JOBS=$(( MAX_CONCURRENT_JOBS_BY_CPU < MAX_CONCURRENT_JOBS_BY_MEM ? MAX_CONCURRENT_JOBS_BY_CPU : MAX_CONCURRENT_JOBS_BY_MEM ))

echo "Max Concurrent Jobs by CPU: $MAX_CONCURRENT_JOBS_BY_CPU"
echo "Max Concurrent Jobs by Memory: $MAX_CONCURRENT_JOBS_BY_MEM"
echo "Recommended Max Concurrent Jobs: $MAX_CONCURRENT_JOBS"

# 小规模测试，运行少量Job以验证脚本功能
echo "Starting a small test run with max concurrent jobs: $MAX_CONCURRENT_JOBS"
USER=$(whoami)

# 获取当前用户在队列中的任务数量
get_current_jobs() {
    jqueues -u $USER --noheader | wc -l  # 你可能需要用你的集群真实命令替换
}

# 任务提交测试
for ((i=1; i<=5; i++)); do  # 小规模测试，提交5个任务
    RANDOM_SEED=$(( i * 1000 ))
    jsub -J test_job -q normal -o output_test_%J.txt -e err_test_%J.txt /hpc/home/zzhi359/FibonacciChain.jl/exm/jsub_compute.sh 8 $i $RANDOM_SEED
    
    current_jobs=$(get_current_jobs)
    echo "Submitted test job $i (current jobs: $current_jobs/$MAX_CONCURRENT_JOBS)"
    
    sleep 1  # 短时间延迟，防止过快提交
done

echo "Test jobs submitted. Check outputs to ensure everything is working."