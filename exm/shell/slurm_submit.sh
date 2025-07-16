#!/bin/bash

#SBATCH -p i64m512u # u = shared nodes; ue = exclusive cores; 32 cores / node # emergency_cpu 128 nodes; 
#SBATCH -J FibonacciJobController
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --exclude=cpu1-1,cpu1-99,cpu1-78,cpu1-94,cpu1-101
#SBATCH -o slurm_controller_%j.txt
#SBATCH -e slurm_controller_%j.err

# =============================================================================
# Enhanced SLURM Job Submission Script with Concurrency Control
# =============================================================================

# Configuration
MAX_CONCURRENT_JOBS=50  # Maximum number of concurrent SLURM jobs
CHECK_INTERVAL=30       # Seconds between status checks
LOG_DIR="logs/slurm"    # Directory for log files
RESOURCE_LOG="resource_usage_$(date +%Y%m%d_%H%M%S).log"

# Job parameters
J_START=8
J_END=20
J_STEP=2
I_START=1
I_END=2000

# Create log directory
mkdir -p "$LOG_DIR"

# Initialize counters and tracking
declare -a submitted_jobs=()
total_submitted=0
total_completed=0
total_failed=0
total_cancelled=0

echo "=== SLURM Job Controller Started at $(date) ==="
echo "Target range: j=${J_START}..${J_END} (step ${J_STEP}), i=${I_START}..${I_END}"
echo "Max concurrent jobs: $MAX_CONCURRENT_JOBS"
echo "Log directory: $LOG_DIR"

# Start resource monitoring
echo "Starting resource monitoring..."
(
  while true; do
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    if command -v free >/dev/null 2>&1; then
      # Linux system
      cpu_usage=$(top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1}')
      memory_usage=$(free -h | grep Mem | awk '{print $3 "/" $2}')
    else
      # Fallback for other systems
      cpu_usage="N/A"
      memory_usage="N/A"
    fi
    
    active_jobs=$(squeue -u $USER -h | wc -l)
    echo "$timestamp - CPU: ${cpu_usage}% - Memory: $memory_usage - Active SLURM jobs: $active_jobs"
    sleep 60
  done
) > "$RESOURCE_LOG" &
resource_monitor_pid=$!

echo "Resource monitoring started (PID: $resource_monitor_pid), logging to: $RESOURCE_LOG"

# =============================================================================
# Job Management Functions
# =============================================================================

# Function: Get current number of active jobs
get_active_job_count() {
    squeue -u $USER -h | wc -l
}

# Function: Get job status
get_job_status() {
    local job_id=$1
    squeue -j "$job_id" -h -o "%T" 2>/dev/null
}

# Function: Wait for job slot to become available
wait_for_slot() {
    while (( $(get_active_job_count) >= MAX_CONCURRENT_JOBS )); do
        echo "Waiting for job slot... (Current active: $(get_active_job_count)/$MAX_CONCURRENT_JOBS)"
        sleep $CHECK_INTERVAL
        
        # Check and update completed jobs
        check_completed_jobs
    done
}

# Function: Check for completed jobs and update counters
check_completed_jobs() {
    local new_submitted_jobs=()
    
    for job_id in "${submitted_jobs[@]}"; do
        status=$(get_job_status "$job_id")
        
        if [[ -z "$status" ]]; then
            # Job no longer in queue - check if it completed successfully
            if sacct -j "$job_id" -n -o State | grep -q "COMPLETED"; then
                ((total_completed++))
                echo "✓ Job $job_id completed successfully (Total completed: $total_completed)"
            elif sacct -j "$job_id" -n -o State | grep -q "FAILED\|TIMEOUT\|OUT_OF_MEMORY"; then
                ((total_failed++))
                echo "✗ Job $job_id failed (Total failed: $total_failed)"
            elif sacct -j "$job_id" -n -o State | grep -q "CANCELLED"; then
                ((total_cancelled++))
                echo "⚠ Job $job_id cancelled (Total cancelled: $total_cancelled)"
            fi
        else
            # Job still active, keep tracking
            new_submitted_jobs+=("$job_id")
        fi
    done
    
    submitted_jobs=("${new_submitted_jobs[@]}")
}

# Function: Submit a job with proper naming and error handling
submit_job() {
    local j=$1
    local i=$2
    local seed=$3
    
    # Generate standardized job name and output files
    local job_name="N${j}_Fibo_S${i}"
    local output_file="${LOG_DIR}/N${j}_sample${i}_job_%j.out"
    local error_file="${LOG_DIR}/N${j}_sample${i}_job_%j.err"
    
    # Submit job and capture job ID
    local submit_output
    submit_output=$(sbatch \
        -p i64m512re \
        -J "$job_name" \
        -o "$output_file" \
        -e "$error_file" \
        /hpc2hdd/home/zzhi359/FibonacciChain.jl/exm/compute.sh "$j" "$i" "$seed" 2>&1)
    
    if [[ $? -eq 0 ]]; then
        # Extract job ID from sbatch output
        local job_id=$(echo "$submit_output" | grep -o '[0-9]\+' | tail -1)
        
        if [[ -n "$job_id" ]]; then
            submitted_jobs+=("$job_id")
            ((total_submitted++))
            echo "🚀 Submitted job $job_id: N=$j, Sample=$i, Seed=$seed (Total: $total_submitted)"
            return 0
        else
            echo "✗ Failed to extract job ID from: $submit_output"
            return 1
        fi
    else
        echo "✗ Failed to submit job for N=$j, Sample=$i: $submit_output"
        ((total_failed++))
        return 1
    fi
}


    fi
}

# =============================================================================
# Main Execution Loop
# =============================================================================

echo ""
echo "Starting job submission..."
cd /hpc2hdd/home/zzhi359/FibonacciChain.jl

start_time=$(date +%s)

for ((j=J_START; j<=J_END; j+=J_STEP)); do
    echo ""
    echo "Processing N=$j..."
    
    for ((i=I_START; i<=I_END; i++)); do
        # Wait for available slot
        wait_for_slot
        
        # Generate random seed
        RANDOM_SEED=$(( (j + i) * 1000 ))
        
        # Submit job
        submit_job "$j" "$i" "$RANDOM_SEED"
        
        # Show progress every 100 jobs
        if (( i % 100 == 0 )); then
            active_count=$(get_active_job_count)
            echo "Progress: N=$j, Sample $i/$I_END (Active: $active_count, Submitted: $total_submitted, Completed: $total_completed)"
        fi
        
        # Brief pause to avoid overwhelming the scheduler
        sleep 0.1
    done
    
    echo "Finished submitting all samples for N=$j"
done

echo ""
echo "All jobs submitted! Waiting for completion..."

# Wait for all remaining jobs to complete
while [[ ${#submitted_jobs[@]} -gt 0 ]]; do
    echo "Waiting for ${#submitted_jobs[@]} remaining jobs to complete..."
    sleep $CHECK_INTERVAL
    check_completed_jobs
done

# =============================================================================
# Cleanup and Final Report
# =============================================================================

# Stop resource monitoring
echo "Stopping resource monitoring..."
kill $resource_monitor_pid 2>/dev/null

end_time=$(date +%s)
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))

echo ""
echo "=== Final Job Execution Summary ==="
echo "Execution time: ${hours}h ${minutes}m ${seconds}s"
echo "Total jobs submitted: $total_submitted"
echo "Successfully completed: $total_completed"
echo "Failed jobs: $total_failed"
echo "Cancelled jobs: $total_cancelled"

if [[ $total_submitted -gt 0 ]]; then
    success_rate=$(( (total_completed * 100) / total_submitted ))
    echo "Success rate: ${success_rate}%"
fi

echo ""
echo "Log files location: $LOG_DIR"
echo "Resource usage log: $RESOURCE_LOG"

if [[ $total_failed -gt 0 ]]; then
    echo ""
    echo "⚠️  Some jobs failed. Check error files in $LOG_DIR for details."
    echo "Failed job logs: ${LOG_DIR}/*job_*.err"
fi

echo ""
echo "Mission accomplished at: $(date)"

# Generate summary report
summary_file="job_summary_$(date +%Y%m%d_%H%M%S).txt"
cat > "$summary_file" << EOF
=== SLURM Job Execution Summary ===
Date: $(date)
Script: $0
Parameters: j=${J_START}..${J_END} (step ${J_STEP}), i=${I_START}..${I_END}
Max concurrent jobs: $MAX_CONCURRENT_JOBS

Results:
- Total submitted: $total_submitted
- Completed: $total_completed
- Failed: $total_failed  
- Cancelled: $total_cancelled
- Success rate: ${success_rate:-0}%
- Duration: ${hours}h ${minutes}m ${seconds}s

Log directory: $LOG_DIR
Resource log: $RESOURCE_LOG
EOF

echo "Summary report saved to: $summary_file"
