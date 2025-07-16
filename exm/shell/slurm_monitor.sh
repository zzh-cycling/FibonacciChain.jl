#!/bin/bash

# =============================================================================
# SLURM Job Monitor Script
# =============================================================================

# Configuration
REFRESH_INTERVAL=10  # Seconds between updates
LOG_DIR="logs/slurm"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function: Clear screen and show header
show_header() {
    clear
    echo "==============================================================================="
    echo "                          SLURM Job Monitor"
    echo "==============================================================================="
    echo "Last updated: $(date)"
    echo "User: $USER"
    echo "Refresh interval: ${REFRESH_INTERVAL}s (Press Ctrl+C to exit)"
    echo "==============================================================================="
}

# Function: Show job statistics
show_statistics() {
    echo ""
    echo "${BLUE}=== Job Statistics ===${NC}"
    
    # Get job counts by state
    local running=$(squeue -u $USER -h -t RUNNING | wc -l)
    local pending=$(squeue -u $USER -h -t PENDING | wc -l)
    local total_active=$((running + pending))
    
    # Get completed jobs from today
    local today=$(date +%Y-%m-%d)
    local completed=$(sacct -u $USER -S $today -s COMPLETED -n | wc -l)
    local failed=$(sacct -u $USER -S $today -s FAILED,TIMEOUT,OUT_OF_MEMORY -n | wc -l)
    local cancelled=$(sacct -u $USER -S $today -s CANCELLED -n | wc -l)
    
    printf "%-15s: ${GREEN}%d${NC}\n" "Running" $running
    printf "%-15s: ${YELLOW}%d${NC}\n" "Pending" $pending
    printf "%-15s: ${BLUE}%d${NC}\n" "Total Active" $total_active
    printf "%-15s: ${GREEN}%d${NC}\n" "Completed" $completed
    printf "%-15s: ${RED}%d${NC}\n" "Failed" $failed
    printf "%-15s: ${YELLOW}%d${NC}\n" "Cancelled" $cancelled
    
    if [[ $((completed + failed + cancelled)) -gt 0 ]]; then
        local success_rate=$(( completed * 100 / (completed + failed + cancelled) ))
        printf "%-15s: ${GREEN}%d%%${NC}\n" "Success Rate" $success_rate
    fi
}

# Function: Show active jobs
show_active_jobs() {
    echo ""
    echo "${BLUE}=== Active Jobs (Top 20) ===${NC}"
    
    if [[ $(squeue -u $USER -h | wc -l) -eq 0 ]]; then
        echo "No active jobs"
        return
    fi
    
    printf "%-10s %-20s %-10s %-15s %-10s %s\n" "JOB ID" "NAME" "STATE" "TIME" "NODES" "NODELIST"
    echo "-------------------------------------------------------------------------------"
    
    squeue -u $USER -h -o "%.10i %.20j %.10T %.15M %.10D %R" | head -20 | while read line; do
        state=$(echo $line | awk '{print $3}')
        if [[ "$state" == "RUNNING" ]]; then
            echo -e "${GREEN}$line${NC}"
        elif [[ "$state" == "PENDING" ]]; then
            echo -e "${YELLOW}$line${NC}"
        else
            echo "$line"
        fi
    done
}

# Function: Show recent completions
show_recent_completions() {
    echo ""
    echo "${BLUE}=== Recent Completions (Last 10) ===${NC}"
    
    local today=$(date +%Y-%m-%d)
    local recent_jobs=$(sacct -u $USER -S $today -s COMPLETED,FAILED,CANCELLED -n -o "JobID,JobName,State,ExitCode,End" --parsable2 | tail -10)
    
    if [[ -z "$recent_jobs" ]]; then
        echo "No recent completions"
        return
    fi
    
    printf "%-10s %-20s %-10s %-8s %s\n" "JOB ID" "NAME" "STATE" "EXIT" "END TIME"
    echo "-------------------------------------------------------------------------------"
    
    echo "$recent_jobs" | while IFS='|' read jobid name state exit_code end_time; do
        # Clean up job ID (remove .batch suffix)
        clean_jobid=$(echo $jobid | sed 's/\..*$//')
        
        if [[ "$state" == "COMPLETED" ]]; then
            echo -e "${GREEN}$clean_jobid $name $state $exit_code $end_time${NC}"
        elif [[ "$state" =~ ^(FAILED|TIMEOUT|OUT_OF_MEMORY)$ ]]; then
            echo -e "${RED}$clean_jobid $name $state $exit_code $end_time${NC}"
        else
            echo -e "${YELLOW}$clean_jobid $name $state $exit_code $end_time${NC}"
        fi
    done
}

# Function: Show system resources
show_resources() {
    echo ""
    echo "${BLUE}=== System Resources ===${NC}"
    
    # Show partition info
    echo "Partition usage:"
    sinfo -p i64m512re -o "%.10P %.5a %.10l %.6D %.6t %.6c %.8z %.8m %.8d %.8w %.10f" | head -2
    
    # Show node status
    echo ""
    echo "Node status:"
    sinfo -p i64m512re -t idle,alloc,mix,drain -o "%.8T %.6D" | sort | uniq -c
}

# Function: Show error summary
show_errors() {
    if [[ ! -d "$LOG_DIR" ]]; then
        return
    fi
    
    echo ""
    echo "${BLUE}=== Recent Error Summary ===${NC}"
    
    local error_files=$(find "$LOG_DIR" -name "*.err" -mtime -1 -size +0c 2>/dev/null | head -5)
    
    if [[ -z "$error_files" ]]; then
        echo "No recent errors"
        return
    fi
    
    echo "Recent error files (non-empty, last 24h):"
    echo "$error_files" | while read file; do
        local size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null)
        echo "  $(basename "$file") (${size} bytes)"
    done
}

# Main monitoring loop
main() {
    echo "Starting SLURM job monitor..."
    echo "Press Ctrl+C to exit"
    
    # Trap Ctrl+C to clean up
    trap 'echo -e "\n\nMonitoring stopped."; exit 0' INT
    
    while true; do
        show_header
        show_statistics
        show_active_jobs
        show_recent_completions
        show_resources
        show_errors
        
        echo ""
        echo "Press Ctrl+C to exit, or wait ${REFRESH_INTERVAL}s for refresh..."
        sleep $REFRESH_INTERVAL
    done
}

# Help function
show_help() {
    echo "SLURM Job Monitor"
    echo ""
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -h, --help     Show this help message"
    echo "  -i, --interval Set refresh interval in seconds (default: $REFRESH_INTERVAL)"
    echo "  -l, --logdir   Set log directory (default: $LOG_DIR)"
    echo ""
    echo "Examples:"
    echo "  $0                    # Start with default settings"
    echo "  $0 -i 5              # Refresh every 5 seconds"
    echo "  $0 -l logs/custom    # Use custom log directory"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -i|--interval)
            REFRESH_INTERVAL="$2"
            shift 2
            ;;
        -l|--logdir)
            LOG_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h for help"
            exit 1
            ;;
    esac
done

# Validate refresh interval
if ! [[ "$REFRESH_INTERVAL" =~ ^[0-9]+$ ]] || [[ "$REFRESH_INTERVAL" -lt 1 ]]; then
    echo "Error: Refresh interval must be a positive integer"
    exit 1
fi

# Start monitoring
main
