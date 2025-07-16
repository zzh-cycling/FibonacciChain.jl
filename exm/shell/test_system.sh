#!/bin/bash

# 测试系统资源检测脚本

echo "=== System Resource Detection Test ==="

# 检测操作系统
OS=$(uname -s)
echo "Operating System: $OS"

# 检测CPU核心数
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

# 检测内存
if [[ "$OS" == "Darwin" ]]; then
    # macOS
    TOTAL_MEM_BYTES=$(sysctl -n hw.memsize)
    TOTAL_MEM_GB=$(( TOTAL_MEM_BYTES / 1024 / 1024 / 1024 ))
    echo "Total Memory: ${TOTAL_MEM_GB}GB (${TOTAL_MEM_BYTES} bytes)"
    
    # 当前内存使用
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
    
    # 当前内存使用
    echo ""
    echo "Current Memory Usage:"
    free -h
fi

# 计算推荐的并发数
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

# 系统负载
echo ""
echo "=== System Load ==="
if [[ "$OS" == "Darwin" ]]; then
    echo "Load average: $(uptime | awk -F'load averages:' '{ print $2 }')"
else
    echo "Load average: $(uptime | awk -F'load average:' '{ print $2 }')"
fi

# 磁盘空间
echo ""
echo "=== Disk Space ==="
df -h .

echo ""
echo "=== Test Complete ==="
