#!/bin/bash

# 作业问题诊断脚本

echo "=== Job Control Diagnostics ==="
echo "Date: $(date)"
echo ""

# 1. 检查当前Julia进程
echo "1. Current Julia processes:"
julia_procs=$(pgrep -f julia)
if [[ -n "$julia_procs" ]]; then
    echo "Found Julia processes:"
    ps -p $julia_procs -o pid,ppid,state,command
else
    echo "No Julia processes found"
fi

echo ""

# 2. 检查nohup进程
echo "2. Current nohup processes:"
nohup_procs=$(pgrep -f nohup)
if [[ -n "$nohup_procs" ]]; then
    echo "Found nohup processes:"
    ps -p $nohup_procs -o pid,ppid,state,command
else
    echo "No nohup processes found"
fi

echo ""

# 3. 检查日志文件
echo "3. Recent log files:"
if [[ -d "logs" ]]; then
    echo "Log directory contents:"
    ls -la logs/ | head -10
    echo ""
    
    echo "Recent log file sizes:"
    find logs -name "job_*.log" -exec ls -lh {} \; | head -5
    
    echo ""
    echo "Checking for errors in recent logs:"
    find logs -name "job_*.log" -mtime -1 -exec grep -l "ERROR\|error\|failed\|Error" {} \; | head -5
else
    echo "No logs directory found"
fi

echo ""

# 4. 系统资源
echo "4. System resources:"
if [[ "$(uname -s)" == "Darwin" ]]; then
    # macOS
    echo "CPU cores: $(sysctl -n hw.ncpu)"
    echo "Memory: $(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))GB"
    echo "Load average: $(uptime | awk -F'load averages:' '{ print $2 }')"
else
    # Linux
    echo "CPU cores: $(nproc)"
    echo "Memory: $(free -h | grep Mem)"
    echo "Load average: $(uptime | awk -F'load average:' '{ print $2 }')"
fi

echo ""

# 5. 检查脚本权限
echo "5. Script permissions:"
ls -la exm/*job*.sh 2>/dev/null || echo "No job scripts found"

echo ""

# 6. 检查Julia环境
echo "6. Julia environment:"
if command -v julia >/dev/null 2>&1; then
    echo "Julia version: $(julia --version)"
    echo "Julia executable: $(which julia)"
    
    # 检查项目环境
    if [[ -f "Project.toml" ]]; then
        echo "Project.toml found"
        echo "Julia project test:"
        julia --project=. -e "println(\"Julia project is working\")" 2>&1
    else
        echo "No Project.toml found"
    fi
else
    echo "Julia not found in PATH"
fi

echo ""

# 7. 进程树查看
echo "7. Process tree (related to current session):"
if command -v pstree >/dev/null 2>&1; then
    pstree -p $$
elif [[ "$(uname -s)" == "Darwin" ]]; then
    ps -o pid,ppid,command -p $$
    echo "Children of current process:"
    ps -o pid,ppid,command | awk -v ppid=$$ '$2 == ppid'
else
    ps --forest -o pid,ppid,cmd -p $$
fi

echo ""

# 8. 建议
echo "=== Diagnostics Complete ==="
echo ""
echo "Common issues and solutions:"
echo "1. If Julia processes are not starting:"
echo "   - Check Julia installation: julia --version"
echo "   - Check project environment: julia --project=. -e 'println(\"OK\")'"
echo "   - Check script permissions: chmod +x exm/*.sh"
echo ""
echo "2. If processes start but don't do work:"
echo "   - Check log files for errors: tail logs/job_*.log"
echo "   - Verify input files exist"
echo "   - Check memory/CPU resources"
echo ""
echo "3. If concurrent control isn't working:"
echo "   - Run test: ./exm/test_job_control.sh"
echo "   - Check for zombie processes: ps aux | grep defunct"
echo "   - Restart terminal session"
echo ""
echo "4. Performance issues:"
echo "   - Reduce CONCURRENCY in script"
echo "   - Monitor system load: top or htop"
echo "   - Check disk space: df -h"
