#!/bin/bash

#SBATCH -p i64m512u
#SBATCH -J BH_PS_gamma0.95
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH -t 06:00:00
#SBATCH -o job-%j.txt

echo "Program started at: $(date)"
cd /hpc2hdd/home/zzhi359/FibonacciChain.jl

julia --project=. -e 'using Pkg; Pkg.instantiate()'

export JULIA_WORKER_TIMEOUT=300

for N in 8 10 12 14 16; do
    echo "=== post-selection gamma=0.95 (ind=10) N=$N at $(date) ==="
    julia --project=. exm/Born_Heisenberg/post_selection.jl $N 10
done

echo "Job completed at: $(date)"
