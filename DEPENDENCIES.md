# 📦 项目依赖说明

## 为什么本地能直接运行？

你的本地电脑已经安装了所有必需的 Python 包：
```bash
numpy     2.1.2   (数值计算核心库)
scipy     1.14.1  (科学计算，使用 scipy.io.loadmat 读取 MATLAB 文件)
logzero   1.7.0   (日志管理库)
```

这些包可能是之前通过以下方式之一安装的：
- `pip3 install numpy scipy logzero` (系统全局安装)
- 安装其他项目时作为依赖被自动安装
- 通过 apt 安装的系统包 (`python3-numpy`, `python3-scipy`)

**新电脑报错的原因：** 新电脑是"干净环境"，没有安装任何第三方 Python 包！

---

## 项目完整依赖列表

### 🔴 外部库（需要安装）

| 库名 | 最低版本 | 用途 | 使用位置 |
|------|----------|------|----------|
| **numpy** | 1.20.0 | 数值计算、数组操作 | 几乎所有模块 |
| **scipy** | 1.7.0 | `scipy.io.loadmat` 读取 MATLAB 文件 | `fem_data_loader.py` |
| **logzero** | 1.7.0 | 日志记录和管理 | `utils/logger.py` |

### 🟢 Python 标准库（无需安装）

以下库随 Python 一起提供，不需要额外安装：
- `os` - 操作系统接口
- `sys` - 系统相关参数和函数
- `argparse` - 命令行参数解析
- `datetime` - 日期时间处理
- `logging` - 日志记录基础设施
- `csv` - CSV 文件读写
- `pathlib` - 面向对象的文件路径
- `subprocess` - 子进程管理
- `shutil` - 文件操作工具

---

## Pipfile 是否完整？

✅ **是的，Pipfile 已经完整！**

```toml
[packages]
numpy = ">=1.20.0"    # ✅ 包含
scipy = ">=1.7.0"     # ✅ 包含
logzero = ">=1.7.0"   # ✅ 包含
```

**为什么看起来库很少？**
- 项目中大量使用的 `os`, `sys`, `argparse` 等都是 Python 标准库
- 标准库不需要列在依赖文件中
- 只有外部安装的包才需要声明

---

## 依赖验证命令

### 检查所有依赖是否满足：
```bash
python3 -c "
import numpy, scipy, logzero
from scipy.io import loadmat
print(f'✅ numpy {numpy.__version__}')
print(f'✅ scipy {scipy.__version__}')
print(f'✅ logzero {logzero.__version__}')
"
```

### 检查当前环境安装的包：
```bash
python3 -m pip list | grep -E "(numpy|scipy|logzero)"
```

### 检查包的安装位置：
```bash
python3 -c "import numpy, scipy, logzero; print(numpy.__file__); print(scipy.__file__); print(logzero.__file__)"
```

---

## 新电脑安装指南

### 方法 1: 使用虚拟环境（推荐）
```bash
./setup_venv.sh
source venv/bin/activate
```

### 方法 2: 使用 Pipenv
```bash
./setup_pipenv.sh
pipenv shell
```

### 方法 3: 系统全局安装
```bash
# WSL/Ubuntu/Debian
sudo apt install python3-numpy python3-scipy
pip3 install logzero

# 或使用 pip 安装所有包
pip3 install -r requirements.txt
```

---

## 常见问题

### Q: 为什么 scipy 需要 1.7.0+？
A: 项目使用 `scipy.io.loadmat` 读取 MATLAB `.mat` 文件，较旧版本可能存在兼容性问题。

### Q: 可以用更新的版本吗？
A: 可以！`>=` 表示最低版本要求，使用更新版本（如当前的 numpy 2.1.2, scipy 1.14.1）完全没问题。

### Q: 为什么不用 requirements.txt？
A: 两种方式都支持：
- `requirements.txt` - 简单直接，pip 标准格式
- `Pipfile` - Pipenv 使用，支持更复杂的依赖管理

两个文件内容相同，选择你喜欢的工具即可。

---

## 技术细节

### scipy.io.loadmat 的使用
在 `fem_data_loader.py` 中用于读取 Fortran 求解器输出的 MATLAB 格式文件：
```python
from scipy.io import loadmat
data = loadmat(filename)
```

### logzero 的使用
在 `utils/logger.py` 中配置日志系统：
```python
from logzero import setup_logger
logger = setup_logger(name='circular_crack', logfile='simulation.log')
```

### numpy 的核心作用
- 数组操作 (`np.array`, `np.zeros`, `np.ones`)
- 线性代数计算
- 网格生成和处理
- 数值积分和插值
