#!/bin/bash

#SBATCH -p i64m512u # u = shared nodes; ue = exclusive cores; 32 cores / node
#SBATCH -J N20FibonaccidynamicsSamples
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --exclude=cpu1-1,cpu1-99,cpu1-78,cpu1-94,cpu1-101
## SBATCH -o job-%j-%a.txt  # 使用作业数组ID来区分输出文件
#SBATCH -o job-%j.txt

##SBATCH -e job.%j.err
##SBATCH --array=1-16
## path="/hpc2hdd/home/zzhi359/quantumErgotropy/julia_shell/erg_scaling"

N=$1
index=$2
# RANDOM_SEED=$3
echo "Received parameters: N=20, index=$index"
# echo "Received parameters: N=20, index=$index, seed=$RANDOM_SEED"

echo "Program started at: $(date)"
cd /hpc2hdd/home/zzhi359/FibonacciChain.jl

export JULIA_NUM_THREADS=16
julia --project=. exm/Bulk_measure/monitored_dynamics.jl $N $index
# julia --project=. exm/Bulk_measure/monitored_dynamics.jl $N $index $RANDOM_SEED
echo "Job completed at: $(date)"
