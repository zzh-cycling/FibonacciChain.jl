#!/bin/bash

#SBATCH -p i64m512u # u = shared nodes; ue = exclusive cores; 32 cores / node
#SBATCH -J N20FibonaccidynamicsSamples
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --exclude=cpu1-1,cpu1-99,cpu1-78,cpu1-94,cpu1-101
## SBATCH -o job-%j-%a.txt  # 使用作业数组ID来区分输出文件
#SBATCH -o job-%j.txt

##SBATCH -e job.%j.err
##SBATCH --array=1-16


> resource_usage.log
(
  while true; do
    echo "$(date) - CPU: $(top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1}')% - Memory: $(free -h | grep Mem | awk '{print $3 "/" $2}')"
    sleep 60
  done
) >> resource_usage.log &

echo "Program started at: $(date)"
cd /hpc2hdd/home/zzhi359/FibonacciChain.jl


for ((j=8; j<=20; j+=2)); do
  for ((i=1; i<=2000; i++)); do
      # 生成一个随机种子
      RANDOM_SEED=($(j+i) * 1000 ) # 通过任务ID来生成种子，确保不同任务之间不重复

      sbatch /hpc2hdd/home/zzhi359/FibonacciChain.jl/exm/compute.sh $j $i $RANDOM_SEED
    
  done
done

echo "mission accomplished at: $(date)"
kill $!
