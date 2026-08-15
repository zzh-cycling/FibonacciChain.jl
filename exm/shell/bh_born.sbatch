#!/bin/bash

#SBATCH -p i64m512u
#SBATCH -J BH_Born_gamma0.95
#SBATCH --array=1-5
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH -t 24:00:00
#SBATCH -o job-%A_%a.txt

echo "Program started at: $(date)"
cd /hpc2hdd/home/zzhi359/FibonacciChain.jl

julia --project=. -e 'using Pkg; Pkg.instantiate()'

export JULIA_WORKER_TIMEOUT=300

Llis=(8 10 12 14 16)
L=${Llis[$((SLURM_ARRAY_TASK_ID-1))]}
# SLURM_NTASKS is not always exported in this cluster's batch environment;
# fall back to nodes * tasks-per-node, and finally to the static 2x32 layout.
NTASKS=${SLURM_NTASKS:-$((${SLURM_JOB_NUM_NODES:-2} * 32))}
NWORKERS=$((NTASKS - 1))

echo "=== Born dynamics gamma=0.95 (ind=10) L=$L, 2000 samples, workers=$NWORKERS at $(date) ==="

# mode 2: generate samples; do not attempt collection if this fails
julia --project=. -p $NWORKERS exm/Born_Heisenberg/Born_dynamics.jl 2 $L 10 1 2000 || exit 1

# mode 1: collect samples and process observables
julia --project=. -p 31 exm/Born_Heisenberg/Born_dynamics.jl 1 $L 10

echo "Job completed at: $(date)"
