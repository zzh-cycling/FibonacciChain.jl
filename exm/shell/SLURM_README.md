# SLURM作业管理系统使用说明

本目录包含了完整的SLURM作业管理系统，具有并发控制、进程管理和监控功能。

## 脚本概览

### 1. `slurm_submit.sh` - 主作业提交脚本
高级的SLURM作业提交脚本，具有：
- **智能并发控制**: 限制同时运行的作业数量
- **实时进程管理**: 跟踪作业状态和完成情况
- **规范化命名**: 标准化的作业名称和日志文件
- **资源监控**: 持续监控系统资源使用
- **错误处理**: 完整的错误检测和报告

### 2. `slurm_monitor.sh` - 实时监控脚本
提供实时的作业状态监控：
- **彩色显示**: 直观的状态颜色编码
- **统计信息**: 运行、等待、完成、失败作业统计
- **资源使用**: 集群资源使用情况
- **错误摘要**: 最近错误文件的快速查看

### 3. `slurm_cleanup.sh` - 清理管理脚本
用于管理和清理作业：
- **作业取消**: 批量取消运行或等待中的作业
- **日志管理**: 清理旧日志文件，备份重要数据
- **状态查看**: 显示当前作业和日志状态

## 主要改进特性

### 🎯 **智能并发控制**
```bash
MAX_CONCURRENT_JOBS=50  # 最大并发作业数
CHECK_INTERVAL=30       # 状态检查间隔（秒）
```

### 📁 **规范化文件命名**
- **作业名称**: `N{j}_Fibo_S{i}` (例: N8_Fibo_S1)
- **输出文件**: `logs/slurm/N{j}_sample{i}_job_{jobid}.out`
- **错误文件**: `logs/slurm/N{j}_sample{i}_job_{jobid}.err`
- **资源日志**: `resource_usage_{timestamp}.log`

### 📊 **完整的状态跟踪**
- 提交、运行、完成、失败作业计数
- 成功率统计
- 执行时间记录
- 自动生成摘要报告

### 🔧 **错误处理和恢复**
- 提交失败自动重试
- 作业状态实时检查
- 详细的错误日志
- 失败作业的诊断信息

## 使用方法

### 基本使用

#### 1. 提交作业
```bash
# 编辑配置参数（如需要）
vim exm/slurm_submit.sh

# 提交到SLURM
sbatch exm/slurm_submit.sh
```

#### 2. 监控作业
```bash
# 实时监控（每10秒刷新）
./exm/slurm_monitor.sh

# 自定义刷新间隔
./exm/slurm_monitor.sh -i 5

# 指定日志目录
./exm/slurm_monitor.sh -l logs/custom
```

#### 3. 管理作业
```bash
# 查看状态
./exm/slurm_cleanup.sh status

# 取消所有等待中的作业
./exm/slurm_cleanup.sh cancel-pending

# 取消所有作业（需确认）
./exm/slurm_cleanup.sh cancel-all

# 强制取消（无需确认）
./exm/slurm_cleanup.sh cancel-all --force
```

#### 4. 日志管理
```bash
# 备份日志文件
./exm/slurm_cleanup.sh backup-logs

# 清理7天前的日志
./exm/slurm_cleanup.sh clean-logs

# 清理3天前的日志
./exm/slurm_cleanup.sh clean-logs --days 3

# 模拟运行（不实际删除）
./exm/slurm_cleanup.sh clean-logs --dry-run
```

### 高级配置

#### 修改作业参数
在 `slurm_submit.sh` 中编辑：
```bash
# 作业范围
J_START=8      # N的起始值
J_END=20       # N的结束值  
J_STEP=2       # N的步长
I_START=1      # 样本起始编号
I_END=2000     # 样本结束编号

# 并发控制
MAX_CONCURRENT_JOBS=50  # 最大并发数
CHECK_INTERVAL=30       # 检查间隔
```

#### 自定义SLURM参数
```bash
# 在submit_job函数中修改
sbatch \
    -p i64m512re \           # 分区
    --time=24:00:00 \        # 时间限制
    --mem=4G \               # 内存限制
    --cpus-per-task=2 \      # CPU数量
    -J "$job_name" \
    -o "$output_file" \
    -e "$error_file" \
    /path/to/compute.sh "$j" "$i" "$seed"
```

## 输出文件结构

```
logs/slurm/
├── N8_sample1_job_12345.out     # 标准输出
├── N8_sample1_job_12345.err     # 错误输出
├── N8_sample2_job_12346.out
├── ...
resource_usage_20250716_143022.log  # 资源使用日志
job_summary_20250716_143022.txt     # 执行摘要
```

## 监控输出说明

### 状态颜色编码
- 🟢 **绿色**: 运行中/已完成
- 🟡 **黄色**: 等待中/已取消
- 🔴 **红色**: 失败/错误

### 监控界面示例
```
===============================================================================
                          SLURM Job Monitor
===============================================================================
Last updated: 2025-07-16 14:30:22
User: zzhi359
Refresh interval: 10s (Press Ctrl+C to exit)
===============================================================================

=== Job Statistics ===
Running        : 45
Pending        : 23
Total Active   : 68
Completed      : 1250
Failed         : 12
Cancelled      : 3
Success Rate   : 96%

=== Active Jobs (Top 20) ===
JOB ID     NAME                 STATE      TIME            NODES      NODELIST
-------------------------------------------------------------------------------
12345      N8_Fibo_S1          RUNNING    0:05:23         1          cpu1-10
12346      N8_Fibo_S2          RUNNING    0:05:20         1          cpu1-11
12347      N8_Fibo_S3          PENDING    0:00:00         1          (Resources)
...
```

## 故障排除

### 常见问题

#### 1. 作业提交失败
```bash
# 检查SLURM状态
sinfo
squeue -u $USER

# 检查分区可用性
sinfo -p i64m512re

# 查看详细错误
./exm/slurm_cleanup.sh status
```

#### 2. 作业运行失败
```bash
# 查看错误日志
ls logs/slurm/*_job_*.err

# 查看特定作业日志
cat logs/slurm/N8_sample1_job_12345.err

# 检查资源使用
tail -f resource_usage_*.log
```

#### 3. 日志文件过多
```bash
# 查看日志状态
./exm/slurm_cleanup.sh status

# 备份并清理
./exm/slurm_cleanup.sh backup-logs
./exm/slurm_cleanup.sh clean-logs --days 7
```

#### 4. 作业卡住不动
```bash
# 查看作业详情
scontrol show job JOBID

# 检查节点状态
sinfo -N

# 强制取消问题作业
scancel JOBID
```

### 调试技巧

#### 查看作业历史
```bash
# 查看今天的作业
sacct -S today

# 查看特定作业详情
sacct -j JOBID --format=JobID,JobName,State,ExitCode,Start,End,Elapsed

# 查看失败作业
sacct -S today -s FAILED,TIMEOUT,OUT_OF_MEMORY
```

#### 监控系统资源
```bash
# 查看集群整体状态
sinfo -o "%.10P %.5a %.10l %.6D %.6t %.6c %.8z %.8m"

# 查看分区详情
scontrol show partition i64m512re

# 查看节点详情
scontrol show node cpu1-10
```

## 最佳实践

### 1. 作业管理
- **分批提交**: 避免一次性提交过多作业
- **监控资源**: 定期检查集群负载
- **合理并发**: 根据集群规模调整并发数
- **及时清理**: 定期清理旧日志和完成的作业

### 2. 错误处理
- **检查日志**: 及时查看错误日志
- **重新提交**: 对失败的作业重新提交
- **资源调整**: 根据失败原因调整资源需求

### 3. 性能优化
- **节点选择**: 避免使用问题节点
- **时间估算**: 合理设置作业时间限制
- **内存管理**: 避免内存溢出
- **依赖管理**: 确保所需软件和数据可用

### 4. 安全备份
- **定期备份**: 备份重要的日志和结果
- **版本控制**: 保存脚本的不同版本
- **文档记录**: 记录重要的配置变更
