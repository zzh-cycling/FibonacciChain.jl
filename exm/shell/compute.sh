#!/bin/bash

#SBATCH -p i64m512u # u = shared nodes; ue = exclusive cores; 32 cores / node, long_cpu, i64m512u,(1024) emergency_cpu,(512) a128m512u,(256) i64m512r(128)
#SBATCH -J N20FibonaccidynamicsSamples
#SBATCH  --ntasks-per-node=64
#SBATCH --nodes=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH -o job-%j.txt

##SBATCH -e job.%j.err
##SBATCH --array=1-16
## path="/hpc2hdd/home/zzhi359/quantumErgotropy/julia_shell/erg_scaling"

N=$1
index=$2
RANDOM_SEED=$3
echo "Received parameters: N=20, index=$index, seed=$RANDOM_SEED"

echo "Program started at: $(date)"
cd /hpc2hdd/home/zzhi359/FibonacciChain.jl

echo "SLURM config: nodes=${SLURM_JOB_NUM_NODES:-8}, ntasks=${SLURM_NTASKS:-512}, cpus-per-task=${SLURM_CPUS_PER_TASK:-1}"
echo "Target mode: 1 pmap task per worker (single-thread)"

julia --project=. -e 'using Pkg; Pkg.instantiate()'

export JULIA_WORKER_TIMEOUT=300
# export JULIA_NUM_THREADS=16
julia --project=. exm/Bulk_measure/monitored_dynamics.jl $N $index $RANDOM_SEED
echo "Job completed at: $(date)"
