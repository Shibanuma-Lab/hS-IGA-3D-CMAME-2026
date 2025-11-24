# WSL 安装指南 / WSL Setup Guide

## 在 Windows Subsystem for Linux (WSL) 上安装

### 方法 1：使用 apt 包管理器（推荐，最快）

```bash
# 更新包列表
sudo apt update

# 安装 Python 和依赖包
sudo apt install python3 python3-numpy python3-scipy

# 验证安装
python3 -c "import numpy, scipy; print('✅ 安装成功')"
```

### 方法 2：使用 pip（标准方法）

```bash
# 安装 pip
sudo apt update
sudo apt install python3-pip

# 安装依赖包
pip3 install --user numpy scipy

# 或使用 requirements.txt
pip3 install --user -r requirements.txt

# 验证安装
python3 -c "import numpy, scipy; print('✅ 安装成功')"
```

### 方法 3：使用自动化脚本

```bash
# 克隆仓库
git clone <repository-url>
cd S-IGA-circular-crack-in-3D-solid-linear

# 运行自动安装脚本
./setup.sh
```

## 常见问题

### 问题 1: pip3: command not found

**解决方案：**
```bash
sudo apt update
sudo apt install python3-pip
```

### 问题 2: Permission denied

**解决方案：**
```bash
# 使用 --user 参数安装到用户目录
pip3 install --user numpy scipy
```

### 问题 3: 包安装成功但无法导入

**解决方案：**
```bash
# 添加用户 Python 包路径到环境变量
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### 问题 4: WSL 网络连接问题导致 pip 下载失败

**解决方案：**
```bash
# 方案 A: 使用 apt 代替 pip（推荐）
sudo apt install python3-numpy python3-scipy

# 方案 B: 配置 pip 使用国内镜像
pip3 install --user -i https://pypi.tuna.tsinghua.edu.cn/simple numpy scipy
```

## 完整安装流程（WSL Ubuntu/Debian）

```bash
# 1. 更新系统
sudo apt update && sudo apt upgrade -y

# 2. 安装基本工具
sudo apt install -y python3 python3-pip git build-essential gfortran

# 3. 克隆项目
git clone <repository-url>
cd S-IGA-circular-crack-in-3D-solid-linear

# 4. 安装 Python 依赖（选择其中一种）
# 方法 A: 使用 apt（推荐）
sudo apt install -y python3-numpy python3-scipy

# 方法 B: 使用 pip
pip3 install --user -r requirements.txt

# 5. 编译 Fortran 求解器
cd sfem_linear
make
cd ..

# 6. 验证安装
python3 -c "import numpy, scipy; print('Python packages: OK')"
ls -lh sfem_linear/bin/sfem_linear && echo "Solver: OK"

# 7. 运行测试
cd circular_crack
python3 main.py --help
```

## 性能优化建议

### WSL 1 vs WSL 2

- **WSL 2** 提供更好的性能和完整的 Linux 内核
- 检查版本：`wsl -l -v`（在 Windows PowerShell 中运行）
- 升级到 WSL 2：`wsl --set-version Ubuntu 2`

### 文件系统性能

- 将项目放在 Linux 文件系统中（`/home/user/`）而不是 Windows 文件系统（`/mnt/c/`）
- Linux 文件系统的 I/O 性能更好

### 内存和 CPU

- 编辑 `%UserProfile%\.wslconfig` 配置 WSL 资源：
  ```ini
  [wsl2]
  memory=8GB
  processors=4
  ```

## 快速诊断脚本

```bash
#!/bin/bash
echo "=== WSL 环境诊断 ==="
echo ""

echo "WSL 版本:"
grep -i microsoft /proc/version

echo ""
echo "Python 版本:"
python3 --version

echo ""
echo "pip3 状态:"
which pip3 && pip3 --version || echo "未安装"

echo ""
echo "NumPy 状态:"
python3 -c "import numpy; print(f'已安装: {numpy.__version__}')" 2>/dev/null || echo "未安装"

echo ""
echo "SciPy 状态:"
python3 -c "import scipy; print(f'已安装: {scipy.__version__}')" 2>/dev/null || echo "未安装"

echo ""
echo "Fortran 编译器:"
which gfortran && gfortran --version | head -1 || echo "未安装"
```

保存为 `diagnose.sh`，运行 `chmod +x diagnose.sh && ./diagnose.sh` 进行诊断。
