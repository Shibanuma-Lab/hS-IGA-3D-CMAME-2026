# 静态分析模式使用指南

## 🎯 功能说明

`--static_only` 标志启用特殊的静态分析模式，用于计算单一静态裂纹配置（无裂纹扩展）。该模式使用 Sneddon 解析解来施加边界条件。

## 🔄 动态模式 vs 静态模式

### 动态模式（默认）
```bash
python3 main.py --step_start 0 --step_end 10
```
- 模拟裂纹动态扩展
- 使用 FEM 插值边界条件
- 参数：c = 4.0mm, hL = 0.04mm

### 静态模式
```bash
python3 main.py --static_only --step_start 5 --step_end 6
```
- 单一静态裂纹配置
- 使用 Sneddon 解析解边界条件
- 参数：c = 1.0mm, hL = 1/48 m

## 🔢 Step Number 的特殊含义

在静态模式下，step number 有特殊含义：

### 动态模式 vs 静态模式

**动态模式**：
- step 表示时间步
- 裂纹半径 a = step * hL（例如：step=100, hL=0.04mm → a=4mm）
- 不同 step 表示裂纹扩展的不同时刻

**静态模式**：
- hL = 1/n（n为整数，如 1/48, 1/24, 1/96）
- c = 1.0（归一化裂纹半径）
- **step number = n（离散化数量）**
- 裂纹半径 a = step * hL = step * (1/n) = step/n

### 示例计算

假设 hL = 1/48（n=48）：

| step | 计算 | 裂纹半径 a | 物理意义 |
|------|------|-----------|----------|
| 12 | 12 * (1/48) | 0.25 | 1/4 裂纹 |
| 24 | 24 * (1/48) | 0.50 | 1/2 裂纹 |
| 36 | 36 * (1/48) | 0.75 | 3/4 裂纹 |
| **48** | **48 * (1/48)** | **1.00** | **完整裂纹** |
| 96 | 96 * (1/48) | 2.00 | 超出定义域 |

**推荐用法**：
- 对于 hL=1/48，使用 step=48 计算完整裂纹（a=c=1.0）
- 对于 hL=1/24，使用 step=24 计算完整裂纹（a=c=1.0）
- step < n 表示部分裂纹（用于参数研究）

## 📋 静态模式参数差异

当 `--static_only` 启用时，程序会自动修改以下参数：

### 1. 局部网格参数 (const_local_mesh.py)
```python
# 动态模式          →  静态模式
hL = 0.04e-3       →  hL = 1.0 / 48.0
aL = 51            →  aL = 12
lL = 15            →  lL = 12
HL = 11            →  HL = 12
```

### 2. 全局网格参数 (simulation_params.py)
```python
# 动态模式          →  静态模式
c = 4.0e-3         →  c = 1.0e-3
WidthG = 8.0e-3    →  WidthG = 2.0
HeightG = 4.0e-3   →  HeightG = 1.0
```

### 3. 网格生成算法 (global_mesh.py)
```python
# 动态模式：均匀单元间距
nodeGx = hG * np.arange(nPtsX)

# 静态模式：非均匀单元间距
hGX = mu_G * WidthG / (nPtsX - 1)
nodeGx = hGX * np.arange(nPtsX)
```

### 4. 边界条件施加 (boundary.py)
```python
# 动态模式：使用 FEM 数据插值
bc_values = self._get_fem_bc_values(bcGs, step)

# 静态模式：使用 Sneddon 解析解
from sneddon_solution import sneddon_displacement_cartesian
disp = sneddon_displacement_cartesian(sigma_app, a, E, nu, coords)
```

## 🚀 使用示例

### 例1：基本静态分析（自动计算 step）
```bash
cd circular_crack
python3 main.py --static_only
```
- **程序会自动计算 step = round(c/hL)**
- 对于默认参数：c=1.0, hL=1/48，自动计算 step = 48
- 计算完整裂纹配置，裂纹半径 a = step × hL = 48 × (1/48) = 1.0
- **无需手动指定 step_start 和 step_end！**

**工作原理**：
- hL = 1/48，c = 1.0
- step = round(c/hL) = round(1.0 / (1/48)) = round(48) = 48
- 裂纹半径 a = step × hL = 48 × (1/48) = 1.0 = c ✓

### 例2：不同 hL 值的静态分析
如果需要不同的网格密度，修改 `const/const_local_mesh.py` 中的 `hL_static` 参数：

```python
# 在 const_local_mesh.py 中
hL_static = 1.0 / 24.0  # 改为 1/24
```

然后运行：
```bash
python3 main.py --static_only
```
- 自动计算 step = round(1.0 / (1/24)) = 24
- 裂纹半径 a = 24 × (1/24) = 1.0

### 例3：仅生成网格（不运行求解器）
```bash
python3 main.py --static_only --meshonly
```
- step 自动计算为 48（基于 c=1.0, hL=1/48）
- 生成静态网格和输入文件
- 不运行 Fortran 求解器
- 用于检查网格和边界条件

### 例4：调试模式
```bash
python3 main.py --static_only --debugmode
```
- step 自动计算
- 启用详细日志输出
- 显示边界条件详细信息

### 例5：批量静态分析（不同网格密度）
如需不同网格密度的分析，需要修改 `const_local_mesh.py` 并运行多次：

```bash
# 手动修改 hL_static 为不同值，然后运行
# hL = 1/24: step=24, a=1.0
# hL = 1/48: step=48, a=1.0
# hL = 1/96: step=96, a=1.0

for hL_divisor in 24 48 96; do
    # 修改 const_local_mesh.py 中的 hL_static = 1.0 / $hL_divisor
    # 然后运行
    python3 main.py --static_only
done
```

**注意**：静态模式下 step 总是自动计算为使 a=c，因此不同 hL 对应不同 step，但都计算相同的裂纹半径 a=1.0

## 📐 Sneddon 解析解

静态模式使用 Sneddon (1946) 的 penny-shaped crack 解析解：

### 裂纹开口位移 (COD)
在裂纹表面 (z=0, r<a):
```
u_z = (2 * σ * (1-ν²) / E) * √(a² - r²)
```

### 远场位移
在边界上 (r, z):
```
u_r = f_r(r, z, a, σ, E, ν)  # 通过椭圆积分计算
u_z = f_z(r, z, a, σ, E, ν)
```

参见 `sneddon_solution.py` 实现细节。

## 🔬 边界条件详细说明

### 对称边界条件 (bcX0, bcY0, bcZ0)
- **bcX0**: x = 0 平面，ux = 0
- **bcY0**: y = 0 平面，uy = 0
- **bcZ0**: z = 0 平面（裂纹外），uz = 0

### 外边界条件 (bcEX)
包括：
- **bcX1**: x = WidthG 平面
- **bcY1**: y = WidthG 平面
- **bcZ1**: z = HeightG 平面

这些边界上施加 Sneddon 解析解计算的位移。

## 📄 静态模式特殊处理

### input.dat 文件差异

静态模式下，`input.dat` 文件的某些参数会被强制设置：

```
# 动态模式                  →  静态模式
solution_type: 2 (dynamic)  →  1 (static)
is_restart: 0 或 1          →  0 (always fresh start)
```

### 初始条件文件

**重要**：静态模式下**不生成初始条件文件**：

```
# 动态模式生成              静态模式跳过
init.u.g.dat (位移)        ❌ 不生成
init.v.g.dat (速度)        ❌ 不生成  
init.a.g.dat (加速度)      ❌ 不生成
init.u.l.dat (局部位移)    ❌ 不生成
init.v.l.dat (局部速度)    ❌ 不生成
init.a.l.dat (局部加速度)  ❌ 不生成
```

**原因**：
1. 静态分析不涉及时间演化，无需速度和加速度
2. 每个静态配置都是独立求解，无需前一步的结果
3. 边界条件直接由 Sneddon 解析解给出

## 📊 输出文件

静态模式生成的文件（**不包括 init.*.dat 文件**）：

```
inputfiles/
└── step00048/  # step 自动计算为 48 (c=1.0, hL=1/48)
    ├── node.g.dat      # 全局网格节点
    ├── elem.g.dat      # 全局网格单元
    ├── bc.g.dat        # 全局边界条件（Sneddon 解）
    ├── node.l.dat      # 局部网格节点
    ├── elem.l.dat      # 局部网格单元
    ├── bc.l.dat        # 局部边界条件
    ├── node.v.dat      # 可视化网格节点
    ├── elem.v.dat      # 可视化网格单元
    ├── index.g.dat     # 索引文件
    ├── weights.g.dat   # NURBS 权重
    ├── input.dat       # 求解器输入文件
    └── load.dat        # 载荷文件（静态模式下为空）
    
    # ❌ 不生成以下文件（仅动态模式需要）：
    # init.u.g.dat, init.v.g.dat, init.a.g.dat
    # init.u.l.dat, init.v.l.dat, init.a.l.dat
```

## ⚙️ 材料参数

静态模式使用 `material_property.py` 中的参数：

```python
E = 3.2e9 Pa        # Young's modulus (PMMA)
ν = 0.35            # Poisson's ratio
σ_app = 1.0e6 Pa    # Applied stress (默认 1 MPa)
```

可以在 `boundary.py` 的 `_boundary_ebc_G_static()` 方法中修改 `sigma_app`。

## 🐛 故障排除

### 问题1：ModuleNotFoundError: No module named 'sneddon_solution'
**解决方案**：确保从 `circular_crack/` 目录运行：
```bash
cd circular_crack
python3 main.py --static_only ...
```

### 问题2：边界条件节点数为零
**解决方案**：检查 step 值是否合理，确保裂纹半径 a = step * hL 在合理范围内。

### 问题3：求解器报错
**解决方案**：
1. 先使用 `--meshonly` 检查输入文件
2. 查看 `logs/` 目录中的日志文件
3. 检查 `bc.g.dat` 中的边界条件值是否合理

## 📚 理论背景

### Sneddon 解 (1946)
经典的 penny-shaped crack 解析解，适用于：
- 无限域中的圆形裂纹
- 均匀拉伸应力
- 线弹性材料

### 适用场景
- 验证数值方法精度
- 基准测试
- 参数敏感性分析
- 教学演示

### 限制
- 假设无限域（实际计算域有限）
- 不考虑裂纹扩展动力学
- 线弹性假设

## 📖 参考文献

1. Sneddon, I. N. (1946). "The distribution of stress in the neighbourhood of a crack in an elastic solid". *Proceedings of the Royal Society A*, 187(1009), 229-260.

2. Tada, H., Paris, P. C., & Irwin, G. R. (2000). *The Stress Analysis of Cracks Handbook* (3rd ed.). ASME Press.

## 💡 高级用法

### 修改应力值
编辑 `boundary.py` 第 398 行：
```python
sigma_app = 2.0e6  # 改为 2 MPa
```

### 使用不同材料
编辑 `const/material_property.py`：
```python
ee = 70e9    # Aluminum
Nu = 0.33
Rho = 2700
```

### 自定义域尺寸
编辑 `const/simulation_params.py` 静态模式部分：
```python
if static_mode:
    WidthG = 2.0   # 改为 3m x 3m x 1.5m
    HeightG = 1.5
```

## ✅ 验证步骤

1. 运行静态分析
2. 检查裂纹表面的开口位移是否符合 Sneddon 公式
3. 对比远场位移与解析解
4. 计算应力强度因子并与理论值对比：K_I = σ√(πa)
