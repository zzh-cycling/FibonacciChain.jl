#!/bin/bash

# =============================================================================
# SLURM Environment Test Script
# =============================================================================

echo "=== SLURM Environment Test ==="
echo "Date: $(date)"
echo "User: $USER"
echo ""

# Test 1: Check SLURM commands availability
echo "1. Checking SLURM commands..."
commands=("sbatch" "squeue" "scancel" "sinfo" "sacct" "scontrol")
for cmd in "${commands[@]}"; do
    if command -v $cmd >/dev/null 2>&1; then
        echo "  ✓ $cmd: Available"
    else
        echo "  ✗ $cmd: Not found"
    fi
done

echo ""

# Test 2: Check partition access
echo "2. Checking partition access..."
if sinfo -p i64m512re >/dev/null 2>&1; then
    echo "  ✓ Partition i64m512re: Accessible"
    sinfo -p i64m512re -o "%.10P %.5a %.10l %.6D %.6t"
else
    echo "  ✗ Partition i64m512re: Not accessible"
fi

echo ""

# Test 3: Check current jobs
echo "3. Current job status..."
job_count=$(squeue -u $USER -h | wc -l)
echo "  Active jobs: $job_count"

if [[ $job_count -gt 0 ]]; then
    echo "  Current jobs:"
    squeue -u $USER -o "%.10i %.20j %.10T %.15M"
fi

echo ""

# Test 4: Check recent job history
echo "4. Recent job history (last 24h)..."
recent_jobs=$(sacct -u $USER -S today -n | wc -l)
echo "  Jobs in last 24h: $recent_jobs"

if [[ $recent_jobs -gt 0 ]]; then
    echo "  Recent job summary:"
    sacct -u $USER -S today -s COMPLETED,FAILED,CANCELLED -n -o "State" | sort | uniq -c
fi

echo ""

# Test 5: Check compute script
echo "5. Checking compute script..."
compute_script="/hpc2hdd/home/zzhi359/FibonacciChain.jl/exm/compute.sh"
if [[ -f "$compute_script" ]]; then
    echo "  ✓ Compute script found: $compute_script"
    if [[ -x "$compute_script" ]]; then
        echo "  ✓ Compute script is executable"
    else
        echo "  ⚠ Compute script is not executable"
    fi
else
    echo "  ✗ Compute script not found: $compute_script"
fi

echo ""

# Test 6: Check project directory
echo "6. Checking project directory..."
project_dir="/hpc2hdd/home/zzhi359/FibonacciChain.jl"
if [[ -d "$project_dir" ]]; then
    echo "  ✓ Project directory found: $project_dir"
    cd "$project_dir"
    if [[ -f "Project.toml" ]]; then
        echo "  ✓ Julia project detected"
    else
        echo "  ⚠ Project.toml not found"
    fi
else
    echo "  ✗ Project directory not found: $project_dir"
fi

echo ""

# Test 7: Check log directory
echo "7. Checking log directory..."
log_dir="logs/slurm"
if [[ -d "$log_dir" ]]; then
    echo "  ✓ Log directory exists: $log_dir"
    log_count=$(find "$log_dir" -name "*.out" -o -name "*.err" | wc -l)
    echo "  Log files: $log_count"
else
    echo "  ⚠ Log directory does not exist: $log_dir"
    echo "    Creating log directory..."
    mkdir -p "$log_dir"
    if [[ -d "$log_dir" ]]; then
        echo "  ✓ Log directory created"
    else
        echo "  ✗ Failed to create log directory"
    fi
fi

echo ""

# Test 8: Check disk space
echo "8. Checking disk space..."
df -h . | head -2

echo ""

# Test 9: Test job submission (dry run)
echo "9. Testing job submission (dry run)..."
test_script=$(mktemp)
cat > "$test_script" << 'EOF'
#!/bin/bash
echo "Test job started at: $(date)"
echo "Working directory: $(pwd)"
echo "Arguments: $@"
sleep 5
echo "Test job completed at: $(date)"
EOF

chmod +x "$test_script"

echo "  Simulating job submission..."
echo "  Command: sbatch -p i64m512re -J TestJob -o test_%j.out $test_script 1 2 3000"
echo "  (This is a dry run - no actual job submitted)"

# Clean up
rm "$test_script"

echo ""

# Summary
echo "=== Test Summary ==="
echo "✓ = Good, ⚠ = Warning, ✗ = Error"
echo ""
echo "If you see any ✗ errors, please resolve them before running the main script."
echo "If you see ⚠ warnings, check if they affect your specific use case."
echo ""
echo "To submit actual jobs, run:"
echo "  sbatch exm/slurm_submit.sh"
echo ""
echo "To monitor jobs, run:"
echo "  ./exm/slurm_monitor.sh"
echo ""
echo "Test completed at: $(date)"
