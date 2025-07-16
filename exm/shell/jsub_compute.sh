#!/bin/sh

N=$1
index=$2
RANDOM_SEED=$3
echo "Received parameters: N=$N, index=$index, seed=$RANDOM_SEED"

echo "Program started at: $(date)"
cd /hpc/home/zzhi359/FibonacciChain.jl

export JULIA_NUM_THREADS=16
julia --project=. exm/Bulk_measure/monitored_dynamics.jl $N $index $RANDOM_SEED
echo "Job completed at: $(date)"