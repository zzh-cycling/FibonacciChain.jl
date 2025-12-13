#!/bin/bash

# System Resource Detection Script

echo "=== System Resource Detection Test ==="

# Detect operating system
OS=$(uname -s)
echo "Operating System: $OS"

# Detect the number of CPU cores
if [[ "$OS" == "Darwin" ]]; then
    # macOS
    CPU_LIMIT=$(sysctl -n hw.ncpu)
    LOGICAL_CPU=$(sysctl -n hw.logicalcpu)
    PHYSICAL_CPU=$(sysctl -n hw.physicalcpu)
    echo "CPU Cores (Total): $CPU_LIMIT"
    echo "CPU Cores (Logical): $LOGICAL_CPU"
    echo "CPU Cores (Physical): $PHYSICAL_CPU"
else
    # Linux
    CPU_LIMIT=$(nproc)
    PHYSICAL_CPU=$(grep "^physical id" /proc/cpuinfo | sort -u | wc -l)
    echo "CPU Cores (Total): $CPU_LIMIT"
    echo "CPU Cores (Physical): $PHYSICAL_CPU"
fi

# Detect memory
if [[ "$OS" == "Darwin" ]]; then
    # macOS
    TOTAL_MEM_BYTES=$(sysctl -n hw.memsize)
    TOTAL_MEM_GB=$(( TOTAL_MEM_BYTES / 1024 / 1024 / 1024 ))
    echo "Total Memory: ${TOTAL_MEM_GB}GB (${TOTAL_MEM_BYTES} bytes)"
    
    # Current memory usage
    if command -v vm_stat >/dev/null 2>&1; then
        echo ""
        echo "Current Memory Usage:"
        vm_stat | head -6
    fi
else
    # Linux
    TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    TOTAL_MEM_GB=$(( TOTAL_MEM_KB / 1024 / 1024 ))
    echo "Total Memory: ${TOTAL_MEM_GB}GB (${TOTAL_MEM_KB}KB)"
    
    # Current memory usage
    echo ""
    echo "Current Memory Usage:"
    free -h
fi

# Calculate recommended concurrency
PER_TASK_MEM_GB=2
MEM_LIMIT=$(( TOTAL_MEM_GB / PER_TASK_MEM_GB ))

CONCURRENCY=$(( MEM_LIMIT < CPU_LIMIT ? MEM_LIMIT : CPU_LIMIT ))
CONCURRENCY=$(( CONCURRENCY > 4 ? CONCURRENCY - 2 : CONCURRENCY / 2 ))
CONCURRENCY=$(( CONCURRENCY > 0 ? CONCURRENCY : 1 ))

echo ""
echo "=== Concurrency Calculation ==="
echo "Memory-based limit (${TOTAL_MEM_GB}GB / ${PER_TASK_MEM_GB}GB per task): $MEM_LIMIT"
echo "CPU-based limit: $CPU_LIMIT"
echo "Recommended concurrency (with safety margin): $CONCURRENCY"

# System load
echo ""
echo "=== System Load ==="
if [[ "$OS" == "Darwin" ]]; then
    echo "Load average: $(uptime | awk -F'load averages:' '{ print $2 }')"
else
    echo "Load average: $(uptime | awk -F'load average:' '{ print $2 }')"
fi

# Disk space
echo ""
echo "=== Disk Space ==="
df -h .

echo ""
echo "=== Test Complete ==="
