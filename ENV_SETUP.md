# DDC 环境搭建指南

## 环境概述

- **GPU**: NVIDIA GeForce RTX 3060 Laptop (CUDA 12.4)
- **系统**: WSL (Windows Subsystem for Linux)
- **Conda**: mamba/conda 混用

## 已知问题

### 1. mamba 在 WSL 下有 pthread 冲突

mamba 在 WSL 环境中执行 `mamba install` 时会触发 `pthread_mutex_lock` 断言错误导致崩溃。

**临时解决**：mamba 用来创建环境和解析依赖，实际包安装改用 `pip`。

### 2. 网络不稳定

PyTorch (CUDA 版本约 768MB) 从 `download.pytorch.org` 下载时会频繁中断。

**解决**：pip 的断点续传比 conda 更可靠，多试几次总能成功。

### 3. PyTorch index URL 会被 pip 误用到所有包

`pip install torch numpy --index-url https://download.pytorch.org/whl/cu124` 会导致 numpy、matplotlib 也从 PyTorch 的 index 查找，找不到就报错。

**解决**：先装普通包，再单独装 PyTorch。

## 推荐安装步骤

### 方式一：手动一步步安装（推荐）

```bash
# 1. 创建空环境
mamba create -n ddc python=3.10 pip -y

# 2. 激活环境
conda activate ddc

# 3. 先装普通包（从 PyPI）
pip install numpy matplotlib pandas

# 4. 再单独装 PyTorch（CUDA 12.4）
pip install torch --index-url https://download.pytorch.org/whl/cu124

# 5. 验证
python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
```

### 方式二：用 conda 安装（非 GPU 版本）

如果网络条件差，或者不需要 GPU 加速：

```bash
conda create -n ddc python=3.10 -y
conda activate ddc
conda install numpy matplotlib pandas pytorch -y
```

CPU 版本的 pytorch 体积小很多，安装更顺利。模拟本身 CPU 足够，只是大数据集生成时 GPU 会快一些。

## 环境复现

### 方式一：用 requirements.txt（推荐）

项目根目录下已有 `requirements.txt`，记录了当前环境的精确依赖版本。

```bash
# venv 方式
python -m venv ddc_env
source ddc_env/bin/activate
pip install -r requirements.txt

# 或 conda + pip
conda create -n ddc python=3.10 -y
conda activate ddc
pip install -r requirements.txt
```

注意：nvidia-* 系列包只在有 CUDA 驱动的主机上装得上。无 GPU 机器装到这些包时会自动跳过，不影响 torch 本身的使用。

### 方式二：用 ddc_env.yaml

`ddc_env.yaml` 记录了目标依赖，但直接 `mamba env create` 不可用（mamba pthread 冲突问题）。建议把它当作参考文档，实际安装按上面的方式一执行。

### 重新导出 requirements.txt

如果后续有包更新，需要重新导出：

```bash
conda activate ddc
pip freeze > requirements.txt
```

## 验证 DDC

```bash
conda activate ddc
cd /home/yizhao/Desktop/github/ddc
python src/ddc.py
```

预期输出：

```
Running DDC Phase 0 Standard Pipeline...
--- T=10 Smoke Test ---
  X shape: torch.Size([11, 50])
  P shape: torch.Size([11, 50])
  Z shape: torch.Size([11, 50])
  N shape: torch.Size([11])
Smoke test passed!
```

## 关键依赖版本

| 包 | 版本 | 说明 |
|----|------|------|
| python | 3.10 | |
| torch | 2.6.0+cu124 | GPU 版本 |
| numpy | 2.2.6 | |
| matplotlib | 3.10.9 | |
| pandas | 2.3.3 | |
