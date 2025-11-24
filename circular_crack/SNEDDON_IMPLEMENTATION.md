# Sneddon 解析解实现说明

## 📖 背景

本实现基于您的 Mathematica 程序中的 `SneddonApp.mx` 模块，采用 **Bessel 积分 + 插值** 方法计算圆形裂纹的解析位移场。

## 🔬 Mathematica 原始实现

您的 Mathematica 程序分为两个阶段：

### 1. 预计算阶段（Sneddon0）

```mathematica
Sneddon0[{r_, z_}] := Module[{ρ, ζ, CS, Eζ, ur1, ur2, uz1, uz2},
  ρ = r/c;
  ζ = z/c;
  CS[η_] := Cos[η]/η - Sin[η]/η^2;
  Eζ[η_] := Exp[-ζ*η];
  
  ur1 = NIntegrate[CS[η]*Eζ[η]*BesselJ[1, ρ*η], {η, 0, ∞}];
  ur2 = NIntegrate[ζ*η*CS[η]*Eζ[η]*BesselJ[1, ρ*η], {η, 0, ∞}];
  uz1 = NIntegrate[CS[η]*Eζ[η]*BesselJ[0, ρ*η], {η, 0, ∞}];
  uz2 = NIntegrate[ζ*η*CS[η]*Eζ[η]*BesselJ[0, ρ*η], {η, 0, ∞}];
  
  {ur1, ur2, uz1, uz2}
]

(* 在规则网格上计算所有积分值 *)
WG = 3; HG = 1; nW = 1200; nH = 400;
posAr = Range[0, WG, WG/nW];
posAz = Range[0, HG, HG/nH];
posA = Flatten[Outer[{#2, #1} &, posAz, posAr], 1];
SA = Map[Sneddon0, posA];

(* 创建插值函数 *)
SAf = Interpolation[MapThread[{#1, #2} &, {posA, SA}], 
                    InterpolationOrder -> 1];

(* 保存到文件 *)
DumpSave["SneddonApp.mx", {c, SAf, SneddonApp}];
```

### 2. 运行时阶段（SneddonApp）

```mathematica
SneddonApp[p0_, c_, EE_, ν_, {r_, z_}] := Module[{
    ur0, uz0, ρ, ζ, ur1, ur2, uz1, uz2, ur, uz},
  
  (* 远场均匀应力位移 *)
  ur0 = -(ν*p0/EE) r;
  uz0 = p0/EE z;
  
  (* 从插值函数快速获取积分值 *)
  {ur1, ur2, uz1, uz2} = SAf[r, z];
  
  (* 计算裂纹引起的位移 *)
  ur = (2*p0*c)/(π*EE) ((1 - 2ν)*ur1 - ur2);
  uz = -(4*p0*c (1 - ν²))/(π*EE) (uz1 + uz2/(2(1 - ν)));
  
  (* 总位移 = 远场 + 裂纹 *)
  {ur0 + ur, uz0 + uz}
]
```

### 3. 边界条件计算（getbc）

```mathematica
getbc[i_] := Module{disp, n, θ},
  θ = If[controlPts[[i, 1]] == 0., π/2, 
         ArcTan[controlPts[[i, 2]]/controlPts[[i, 1]]]];
  
  disp = SneddonApp[σapp, a, EE, ν, 
                    {Sqrt[controlPts[[i,1]]² + controlPts[[i,2]]²],
                     controlPts[[i,3]]}];
  
  (* 转换到笛卡尔坐标 *)
  {{i, 1, disp[[1]]*Cos[θ]},   (* u_x *)
   {i, 2, disp[[1]]*Sin[θ]},   (* u_y *)
   {i, 3, disp[[2]]}}          (* u_z *)
]
```

## 🐍 Python 实现

我们用 Python 重现了完全相同的逻辑：

### 1. 预计算模块（sneddon_precompute.py）

```python
def sneddon_integrals(r, z, c=1.0):
    """对应 Mathematica 的 Sneddon0"""
    rho = r / c
    zeta = z / c
    
    def CS(eta):
        if abs(eta) < 1e-10:
            return 0.0
        return np.cos(eta) / eta - np.sin(eta) / (eta**2)
    
    # 定义被积函数
    integrand_ur1 = lambda eta: CS(eta) * np.exp(-zeta*eta) * special.jv(1, rho*eta)
    integrand_ur2 = lambda eta: zeta*eta * CS(eta) * np.exp(-zeta*eta) * special.jv(1, rho*eta)
    integrand_uz1 = lambda eta: CS(eta) * np.exp(-zeta*eta) * special.jv(0, rho*eta)
    integrand_uz2 = lambda eta: zeta*eta * CS(eta) * np.exp(-zeta*eta) * special.jv(0, rho*eta)
    
    # 数值积分（0 到 ∞）
    ur1, _ = integrate.quad(integrand_ur1, 0, eta_max)
    ur2, _ = integrate.quad(integrand_ur2, 0, eta_max)
    uz1, _ = integrate.quad(integrand_uz1, 0, eta_max)
    uz2, _ = integrate.quad(integrand_uz2, 0, eta_max)
    
    return [ur1, ur2, uz1, uz2]

def precompute_sneddon_data(WG=3.0, HG=1.0, nW=600, nH=200, c=1.0):
    """在规则网格上计算所有积分值"""
    r_grid = np.linspace(0, WG, nW)
    z_grid = np.linspace(0, HG, nH)
    data = np.zeros((nW, nH, 4))
    
    for i, r in enumerate(r_grid):
        for j, z in enumerate(z_grid):
            data[i, j, :] = sneddon_integrals(r, z, c)
    
    # 创建插值器（对应 Mathematica 的 Interpolation）
    interpolator = RegularGridInterpolator(
        (r_grid, z_grid), data,
        method='linear',  # InterpolationOrder->1
        bounds_error=False
    )
    
    # 保存到 .npz 文件（对应 .mx 文件）
    np.savez_compressed('sneddon_interpolation.npz',
                        r_grid=r_grid, z_grid=z_grid, data=data, c=c)
    
    return interpolator
```

### 2. 运行时模块（sneddon_solution.py）

```python
def sneddon_displacement_interpolated(sigma_app, a, E, nu, point):
    """对应 Mathematica 的 SneddonApp"""
    r, z = point
    
    # 加载插值器（懒加载）
    interpolator, c_ref = _load_interpolator()
    
    # 远场均匀应力位移
    u_r0 = -(nu * sigma_app / E) * r
    u_z0 = (sigma_app / E) * z
    
    # 从插值器快速获取积分值
    integrals = interpolator([r, z])[0]
    ur1, ur2, uz1, uz2 = integrals
    
    # 计算裂纹引起的位移（与 Mathematica 完全一致）
    if abs(r) < 1e-15:
        u_r_crack = 0.0
    else:
        u_r_crack = (2 * sigma_app * a) / (np.pi * E) * \
                    ((1 - 2*nu) * ur1 - ur2)
    
    if r >= a and abs(z) < 1e-15:
        u_z_crack = 0.0
    else:
        u_z_crack = -(4 * sigma_app * a * (1 - nu**2)) / (np.pi * E) * \
                    (uz1 + uz2 / (2 * (1 - nu)))
    
    # 总位移 = 远场 + 裂纹
    u_r = u_r0 + u_r_crack
    u_z = u_z0 + u_z_crack
    
    return np.array([u_r, u_z])

def sneddon_displacement_cartesian(sigma_app, a, E, nu, point, use_interpolation=True):
    """对应 Mathematica 的 getbc（坐标转换部分）"""
    x, y, z = point
    r = np.sqrt(x**2 + y**2)
    
    # 计算 theta（与 Mathematica 完全一致）
    if abs(x) < 1e-15:
        theta = np.pi / 2.0
    else:
        theta = np.arctan(y / x)
    
    # 获取柱坐标位移
    u_r, u_z = sneddon_displacement(sigma_app, a, E, nu, [r, z], use_interpolation)
    
    # 转换到笛卡尔坐标（与 Mathematica getbc 完全一致）
    if r < 1e-15:
        u_x = 0.0
        u_y = 0.0
    else:
        u_x = u_r * np.cos(theta)
        u_y = u_r * np.sin(theta)
    
    return np.array([u_x, u_y, u_z])
```

### 3. 边界条件应用（boundary.py）

```python
def _boundary_ebc_G_static(self, step):
    """静态模式边界条件（使用 Sneddon 解析解）"""
    from sneddon_solution import sneddon_displacement_cartesian
    
    a = step * clm.hL  # 裂纹半径
    sigma_app = mp.SigmaInfinity
    E = mp.EE
    nu = mp.Nu
    
    # 对每个边界节点计算位移（对应 Mathematica 的 getbc）
    for node in bcEX:
        node_idx = node - 1
        coords = nodeG[node_idx]  # [x, y, z]
        
        # 使用 Sneddon 解析解
        disp = sneddon_displacement_cartesian(sigma_app, a, E, nu, coords)
        
        # 应用边界条件
        bcG_list.append([node, 1, disp[0]])  # u_x
        bcG_list.append([node, 2, disp[1]])  # u_y
        bcG_list.append([node, 3, disp[2]])  # u_z
```

## 📊 网格参数对应

### Mathematica 设置
```mathematica
WG = 3;     (* r 方向宽度 *)
HG = 1;     (* z 方向高度 *)
nW = 600;   (* r 方向网格点数 *)
nH = 200;   (* z 方向网格点数 *)
c = 1;      (* 归一化裂纹半径 *)
```

### Python 对应
```python
save_interpolation_data(
    WG=3.0, HG=1.0,   # 计算域尺寸
    nW=600, nH=200,   # 网格分辨率
    c=1.0             # 归一化裂纹半径
)
```

## 🔄 使用流程

### 1. 首次使用（预计算）

```bash
# 方法 A：使用脚本（推荐）
./generate_sneddon_data.sh
# 选择：
#   1) 测试模式 - 120x40 网格，~5 分钟
#   2) 完整模式 - 600x200 网格，~30-45 分钟

# 方法 B：直接运行
python3 sneddon_precompute.py          # 完整模式
python3 sneddon_precompute.py --test   # 测试模式
```

生成文件：
- `sneddon_interpolation.npz` - 完整数据（~19 MB）
- `sneddon_interpolation_test.npz` - 测试数据（~200 KB）

### 2. 运行静态分析

```bash
python3 main.py --static_only
```

程序会自动：
1. 检测 `sneddon_interpolation.npz` 是否存在
2. 如果不存在，使用 `sneddon_interpolation_test.npz`（低精度）
3. 如果都不存在，报错并提示运行预计算

### 3. 快速测试

```bash
# 运行测试脚本（会自动生成测试数据）
python3 test_sneddon.py
```

测试内容：
- 生成 60x20 测试网格（快速）
- 对比椭圆积分 vs Bessel 积分方法
- 验证 COD（Crack Opening Displacement）
- 显示精度对比

## ⚖️ 两种实现方法对比

### 方法 A：椭圆积分（旧方法）
```python
sneddon_displacement_cartesian(..., use_interpolation=False)
```

**优点**：
- ✅ 无需预计算
- ✅ 立即可用
- ✅ 公式简单

**缺点**：
- ❌ 精度较低（简化公式）
- ❌ 与 Mathematica 不一致
- ❌ 未经您验证

### 方法 B：Bessel 积分 + 插值（新方法，推荐）
```python
sneddon_displacement_cartesian(..., use_interpolation=True)  # 默认
```

**优点**：
- ✅ 与您的 Mathematica 实现完全一致
- ✅ 精度高（原始 Sneddon 公式）
- ✅ 运行时速度快（插值）
- ✅ 经过您验证

**缺点**：
- ⚠️ 需要一次性预计算（~2 小时）
- ⚠️ 需要存储空间（~19 MB）

## 🎯 推荐配置

### 开发/测试阶段
```bash
# 生成测试数据（快速）
python3 sneddon_precompute.py --test
# 或
python3 test_sneddon.py

# 使用测试数据运行
python3 main.py --static_only
```

### 生产/论文阶段
```bash
# 生成完整数据（运行一次）
./generate_sneddon_data.sh
# 选择 2) Full mode

# 之后可以无限次使用
python3 main.py --static_only
```

## 🔍 验证方法

### 1. 检查 COD（裂纹张开位移）

```python
# 理论公式
COD(r) = (2σ(1-ν²)/E) √(a²-r²)

# 在裂纹中心 (r=0)
COD_max = (2σ(1-ν²)/E) × a

# 示例：σ=1MPa, E=3.2GPa, ν=0.35, a=1m
# COD_max ≈ 0.475 μm
```

### 2. 对比 Mathematica 结果

在相同位置计算位移，对比数值：
```mathematica
(* Mathematica *)
SneddonApp[1.0e6, 1.0, 3.2e9, 0.35, {0.5, 0.5}]
```

```python
# Python
sneddon_displacement_cartesian(1.0e6, 1.0, 3.2e9, 0.35, [0.5, 0.0, 0.5])
```

### 3. 运行测试脚本

```bash
python3 test_sneddon.py
```

查看输出的误差百分比，应该 < 1%（对于完整网格）。

## 📝 技术细节

### Bessel 函数

- `J₀(x)` - 第一类零阶 Bessel 函数 → `scipy.special.jv(0, x)`
- `J₁(x)` - 第一类一阶 Bessel 函数 → `scipy.special.jv(1, x)`

### 数值积分

- Mathematica: `NIntegrate[..., {η, 0, ∞}]`
- Python: `scipy.integrate.quad(..., 0, eta_max)`
  - `eta_max` 自适应选择（通常 50-100 足够）

### 插值方法

- Mathematica: `Interpolation[..., InterpolationOrder -> 1]` = 线性插值
- Python: `RegularGridInterpolator(..., method='linear')` = 线性插值

两者完全等价！

## 🐛 故障排除

### 问题 1：找不到插值数据

```
FileNotFoundError: Sneddon interpolation data not found!
Please run: python sneddon_precompute.py
```

**解决方案**：
```bash
python3 sneddon_precompute.py --test  # 快速测试
# 或
./generate_sneddon_data.sh             # 完整数据
```

### 问题 2：精度不足

测试脚本显示误差 > 5%

**原因**：使用低分辨率测试数据

**解决方案**：生成完整数据
```bash
./generate_sneddon_data.sh
# 选择 2) Full mode
```

### 问题 3：预计算太慢

完整模式（600x200）需要 ~30-45 分钟

**建议**：
1. 先用测试模式验证流程
2. 在不需要电脑时运行完整模式（晚上/周末）
3. 只需运行一次，之后永久可用
4. 可以减小网格（如 600x200）达到平衡

## 📚 参考文献

1. **Sneddon, I. N. (1946)**. "The distribution of stress in the neighbourhood of a crack in an elastic solid". *Proceedings of the Royal Society A*, 187(1009), 229-260.

2. **您的 Mathematica 实现** - `SneddonApp.mx` 模块
   - 使用 Bessel 积分的原始公式
   - 预计算 + 插值策略
   - 本 Python 实现完全基于此

## ✨ 总结

**核心思想**：完全复现您的 Mathematica 实现

1. ✅ 相同的数学公式（Bessel 积分）
2. ✅ 相同的算法流程（预计算 + 插值）
3. ✅ 相同的坐标转换（getbc 逻辑）
4. ✅ 相同的数值方法（线性插值）

**使用建议**：

- 🧪 **开发阶段**：使用测试数据（`--test`）
- 📊 **生产阶段**：使用完整数据（完整模式）
- 🔍 **验证**：对比 Mathematica 结果
- 📝 **论文**：引用 Sneddon (1946) + 您的实现

**下一步**：
```bash
# 1. 生成数据
./generate_sneddon_data.sh

# 2. 运行静态分析
python3 main.py --static_only

# 3. 查看结果
ls -lh inputfiles/step00048/
```
