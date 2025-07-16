#!/bin/bash

# =============================================================================
# SLURM Job Cleanup Script
# =============================================================================

# Configuration
LOG_DIR="logs/slurm"
BACKUP_DIR="backup/$(date +%Y%m%d_%H%M%S)"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Function: Show help
show_help() {
    echo "SLURM Job Cleanup Script"
    echo ""
    echo "Usage: $0 [options] [action]"
    echo ""
    echo "Actions:"
    echo "  cancel-all     Cancel all pending/running jobs"
    echo "  cancel-pending Cancel only pending jobs"
    echo "  clean-logs     Clean old log files (keeps last 7 days)"
    echo "  backup-logs    Backup log files before cleaning"
    echo "  status         Show cleanup status"
    echo ""
    echo "Options:"
    echo "  -h, --help     Show this help message"
    echo "  -f, --force    Skip confirmation prompts"
    echo "  -d, --days     Days to keep for log cleanup (default: 7)"
    echo "  --dry-run      Show what would be done without doing it"
    echo ""
    echo "Examples:"
    echo "  $0 status                    # Show current status"
    echo "  $0 cancel-pending           # Cancel pending jobs"
    echo "  $0 clean-logs --days 3      # Clean logs older than 3 days"
    echo "  $0 cancel-all --force       # Cancel all jobs without confirmation"
}

# Function: Show current status
show_status() {
    echo "${BLUE}=== SLURM Job Status ===${NC}"
    
    local running=$(squeue -u $USER -h -t RUNNING | wc -l)
    local pending=$(squeue -u $USER -h -t PENDING | wc -l)
    local total_active=$((running + pending))
    
    printf "Running jobs: ${GREEN}%d${NC}\n" $running
    printf "Pending jobs: ${YELLOW}%d${NC}\n" $pending
    printf "Total active: ${BLUE}%d${NC}\n" $total_active
    
    if [[ $total_active -gt 0 ]]; then
        echo ""
        echo "Active jobs:"
        squeue -u $USER -o "%.10i %.20j %.10T %.15M %.10D %R"
    fi
    
    echo ""
    echo "${BLUE}=== Log Directory Status ===${NC}"
    if [[ -d "$LOG_DIR" ]]; then
        local log_count=$(find "$LOG_DIR" -name "*.out" -o -name "*.err" | wc -l)
        local log_size=$(du -sh "$LOG_DIR" 2>/dev/null | cut -f1)
        printf "Log files: %d\n" $log_count
        printf "Log directory size: %s\n" $log_size
        
        echo ""
        echo "Log file age distribution:"
        find "$LOG_DIR" -name "*.out" -o -name "*.err" | while read file; do
            echo "$file"
        done | xargs stat -f "%m %N" 2>/dev/null | while read mtime file; do
            days_old=$(( ($(date +%s) - mtime) / 86400 ))
            echo "$days_old"
        done | sort -n | uniq -c | awk '{printf "  %d days old: %d files\n", $2, $1}'
    else
        echo "Log directory does not exist: $LOG_DIR"
    fi
}

# Function: Cancel jobs
cancel_jobs() {
    local job_type="$1"  # "all", "pending", or "running"
    
    local jobs_to_cancel
    case $job_type in
        "all")
            jobs_to_cancel=$(squeue -u $USER -h -o "%i")
            echo "${YELLOW}Cancelling ALL jobs...${NC}"
            ;;
        "pending")
            jobs_to_cancel=$(squeue -u $USER -h -t PENDING -o "%i")
            echo "${YELLOW}Cancelling PENDING jobs...${NC}"
            ;;
        "running")
            jobs_to_cancel=$(squeue -u $USER -h -t RUNNING -o "%i")
            echo "${YELLOW}Cancelling RUNNING jobs...${NC}"
            ;;
        *)
            echo "${RED}Error: Invalid job type: $job_type${NC}"
            return 1
            ;;
    esac
    
    if [[ -z "$jobs_to_cancel" ]]; then
        echo "No $job_type jobs to cancel"
        return 0
    fi
    
    local job_count=$(echo "$jobs_to_cancel" | wc -l)
    echo "Found $job_count job(s) to cancel:"
    echo "$jobs_to_cancel"
    
    if [[ "$FORCE" != "true" && "$DRY_RUN" != "true" ]]; then
        echo ""
        read -p "Are you sure you want to cancel these jobs? [y/N]: " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Cancelled"
            return 0
        fi
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY RUN] Would cancel jobs: $jobs_to_cancel"
        return 0
    fi
    
    local cancelled=0
    local failed=0
    
    echo "$jobs_to_cancel" | while read job_id; do
        if scancel "$job_id" 2>/dev/null; then
            echo "${GREEN}✓${NC} Cancelled job $job_id"
            ((cancelled++))
        else
            echo "${RED}✗${NC} Failed to cancel job $job_id"
            ((failed++))
        fi
    done
    
    echo ""
    echo "${GREEN}Successfully cancelled: $cancelled${NC}"
    echo "${RED}Failed to cancel: $failed${NC}"
}

# Function: Backup logs
backup_logs() {
    if [[ ! -d "$LOG_DIR" ]]; then
        echo "Log directory does not exist: $LOG_DIR"
        return 1
    fi
    
    echo "${BLUE}Creating backup of log files...${NC}"
    
    mkdir -p "$BACKUP_DIR"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY RUN] Would backup $LOG_DIR to $BACKUP_DIR"
        return 0
    fi
    
    if cp -r "$LOG_DIR"/* "$BACKUP_DIR"/ 2>/dev/null; then
        local file_count=$(find "$BACKUP_DIR" -type f | wc -l)
        echo "${GREEN}✓${NC} Backed up $file_count files to $BACKUP_DIR"
        
        # Create backup info file
        cat > "$BACKUP_DIR/backup_info.txt" << EOF
Backup created: $(date)
Original location: $LOG_DIR
Backup script: $0
User: $USER
Files backed up: $file_count
EOF
        
        return 0
    else
        echo "${RED}✗${NC} Failed to create backup"
        return 1
    fi
}

# Function: Clean old logs
clean_logs() {
    local days_to_keep=${DAYS_TO_KEEP:-7}
    
    if [[ ! -d "$LOG_DIR" ]]; then
        echo "Log directory does not exist: $LOG_DIR"
        return 1
    fi
    
    echo "${BLUE}Cleaning log files older than $days_to_keep days...${NC}"
    
    # Find old files
    local old_files=$(find "$LOG_DIR" -name "*.out" -o -name "*.err" | xargs stat -f "%m %N" 2>/dev/null | while read mtime file; do
        days_old=$(( ($(date +%s) - mtime) / 86400 ))
        if [[ $days_old -gt $days_to_keep ]]; then
            echo "$file"
        fi
    done)
    
    if [[ -z "$old_files" ]]; then
        echo "No files older than $days_to_keep days found"
        return 0
    fi
    
    local file_count=$(echo "$old_files" | wc -l)
    echo "Found $file_count old file(s):"
    echo "$old_files" | head -10
    if [[ $file_count -gt 10 ]]; then
        echo "... and $((file_count - 10)) more"
    fi
    
    if [[ "$FORCE" != "true" && "$DRY_RUN" != "true" ]]; then
        echo ""
        read -p "Delete these files? [y/N]: " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Cancelled"
            return 0
        fi
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY RUN] Would delete $file_count files"
        return 0
    fi
    
    local deleted=0
    local failed=0
    
    echo "$old_files" | while read file; do
        if rm "$file" 2>/dev/null; then
            echo "${GREEN}✓${NC} Deleted $(basename "$file")"
            ((deleted++))
        else
            echo "${RED}✗${NC} Failed to delete $(basename "$file")"
            ((failed++))
        fi
    done
    
    echo ""
    echo "${GREEN}Successfully deleted: $deleted${NC}"
    echo "${RED}Failed to delete: $failed${NC}"
}

# Parse command line arguments
FORCE=false
DRY_RUN=false
DAYS_TO_KEEP=7
ACTION=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -f|--force)
            FORCE=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -d|--days)
            DAYS_TO_KEEP="$2"
            shift 2
            ;;
        cancel-all|cancel-pending|cancel-running|clean-logs|backup-logs|status)
            ACTION="$1"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h for help"
            exit 1
            ;;
    esac
done

# Validate days parameter
if ! [[ "$DAYS_TO_KEEP" =~ ^[0-9]+$ ]] || [[ "$DAYS_TO_KEEP" -lt 0 ]]; then
    echo "Error: Days must be a non-negative integer"
    exit 1
fi

# Execute action
case $ACTION in
    "status")
        show_status
        ;;
    "cancel-all")
        cancel_jobs "all"
        ;;
    "cancel-pending")
        cancel_jobs "pending"
        ;;
    "cancel-running")
        cancel_jobs "running"
        ;;
    "backup-logs")
        backup_logs
        ;;
    "clean-logs")
        clean_logs
        ;;
    "")
        echo "No action specified. Use -h for help."
        exit 1
        ;;
    *)
        echo "Unknown action: $ACTION"
        exit 1
        ;;
esac
