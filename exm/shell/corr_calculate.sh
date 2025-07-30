#!/bin/bash
#SBATCH -p i64m512ue # u = shared nodes; ue = exclusive cores; 32 cores / node, long_cpu, i64m512u,(1024) emergency_cpu,(512) a128m512u,(256) i64m512r(128)
#SBATCH -J N${1}FiboSamples%a
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --exclude=cpu1-1,cpu1-99,cpu1-78,cpu1-94,cpu1-101
#SBATCH -o job-%A_%a.txt
# #SBATCH --array=1-1000

N=$1
# real_index=$(( (SLURM_ARRAY_TASK_ID - 1) * 10 + 1 ))
real_index=$2


echo "Received parameters: N=$N, index=$real_index"
echo "Program started at: $(date)"

cd /hpc2hdd/home/zzhi359/FibonacciChain.jl
export JULIA_NUM_THREADS=10
julia --project=. exm/Bulk_measure/corr_calculate.jl "$N" "$real_index"

echo "Job completed at: $(date)"