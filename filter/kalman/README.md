# 卡尔曼滤波 (Kalman Filter) 实现解析

> 代码位置：`algorithm/filter/kalman/`

---

## 标准卡尔曼滤波五步公式

离散线性卡尔曼滤波分为 **预测（Predict）** 和 **修正（Update）** 两个阶段：

| 步 | 名称 | 公式 | 代码函数 |
|----|------|------|----------|
| Eq1 | 先验状态估计 | $\hat{x}_{k\vert k-1} = F \cdot \hat{x}_{k-1\vert k-1} + B \cdot u_k$ | `UpdateXhatMinus()` |
| Eq2 | 先验协方差估计 | $P_{k\vert k-1} = F \cdot P_{k-1\vert k-1} \cdot F^\mathsf{T} + Q$ | `UpdatePminus()` |
| Eq3 | 卡尔曼增益 | $K_k = P_{k\vert k-1} \cdot H^\mathsf{T} \cdot (H \cdot P_{k\vert k-1} \cdot H^\mathsf{T} + R)^{-1}$ | `SetGain()` |
| Eq4 | 后验状态估计 | $\hat{x}_{k\vert k} = \hat{x}_{k\vert k-1} + K_k \cdot (z_k - H \cdot \hat{x}_{k\vert k-1})$ | `UpdateXhat()` |
| Eq5 | 后验协方差估计 | $P_{k\vert k} = P_{k\vert k-1} - K_k \cdot H \cdot P_{k\vert k-1}$ | `UpdateCovariance()` |

---

## 变量 ⇔ 公式参数映射

### 核心状态向量

| 代码变量 | 公式符号 | 维度 | 描述 |
|----------|----------|------|------|
| `xhat` | $\hat{x}_{k\vert k}$ | `xhatSize × 1` | **后验**状态估计 — 修正后的最终结果 |
| `xhatminus` | $\hat{x}_{k\vert k-1}$ | `xhatSize × 1` | **先验**状态估计 — 预测值 |
| `u` | $u_k$ | `uSize × 1` | 控制输入向量（可无） |
| `z` | $z_k$ | `zSize × 1` | 观测/测量向量 |

### 协方差矩阵

| 代码变量 | 公式符号 | 维度 | 描述 |
|----------|----------|------|------|
| `P` | $P_{k\vert k}$ | `xhatSize × xhatSize` | **后验**估计协方差 |
| `Pminus` | $P_{k\vert k-1}$ | `xhatSize × xhatSize` | **先验**估计协方差 |
| `Q` | $Q$ | `xhatSize × xhatSize` | **过程噪声**协方差 — 模型不确定性 |
| `R` | $R$ | `zSize × zSize` | **观测噪声**协方差 — 传感器噪声 |
| `S` | $S$ | `zSize × zSize` | 中间矩阵 $S = H \cdot P_{k\vert k-1} \cdot H^\mathsf{T} + R$ |

### 系统矩阵

| 代码变量 | 公式符号 | 维度 | 描述 |
|----------|----------|------|------|
| `F` / `FT` | $F$ / $F^\mathsf{T}$ | `xhatSize × xhatSize` | 状态转移矩阵及其转置 |
| `B` | $B$ | `xhatSize × uSize` | 控制输入矩阵 |
| `H` / `HT` | $H$ / $H^\mathsf{T}$ | `zSize × xhatSize` | 观测矩阵及其转置 |

### 增益与中间变量

| 代码变量 | 公式符号 | 维度 | 描述 |
|----------|----------|------|------|
| `K` | $K_k$ | `xhatSize × zSize` | **卡尔曼增益** — 决定预测与观测的权重 |
| `temp_matrix` | — | `xhatSize × xhatSize` | 通用矩阵临时缓冲区 |
| `temp_matrix1` | — | `xhatSize × xhatSize` | 通用矩阵临时缓冲区 |
| `temp_vector` | — | `xhatSize × 1` | 通用向量临时缓冲区 |
| `temp_vector1` | — | `xhatSize × 1` | 通用向量临时缓冲区 |

### 自动观测调节（`UseAutoAdjustment = true` 时使用）

| 代码变量 | 描述 |
|----------|------|
| `MeasuredVector` | 原始观测输入，0 值表示该通道无效 |
| `MeasurementMap` | 每个传感器通道对应的状态索引 |
| `MeasurementDegree` | 每个通道的观测系数（填入 H 矩阵对应位置） |
| `MatR_DiagonalElements` | 每个通道的观测噪声方差（R 对角元） |
| `MeasurementValidNum` | 本轮有效观测个数 |

### 协方差下限

| 代码变量 | 描述 |
|----------|------|
| `StateMinVariance` | 每个状态维度的协方差下限，防止滤波器"过度自信" |

---

## 代码 ⇔ 公式逐句对照

### Eq1: 先验状态估计 — `UpdateXhatMinus()`

```
x̂(k|k-1) = F · x̂(k-1|k-1) + B · u(k)
```

```cpp
// 有控制输入：x̂⁻ = F·x̂ + B·u
temp_vector  = F * xhat;                    // F · x̂(k-1|k-1)
temp_vector1 = B * u;                       // B · u(k)
xhatminus    = temp_vector + temp_vector1;  // x̂(k|k-1)

// 无控制输入：x̂⁻ = F·x̂  (简化)
xhatminus = F * xhat;
```

### Eq2: 先验协方差估计 — `UpdatePminus()`

```
P(k|k-1) = F · P(k-1|k-1) · Fᵀ + Q
```

```cpp
FT          = Fᵀ;              // 先算 F 转置
Pminus      = F * P;           // F · P(k-1|k-1)
temp_matrix = Pminus * FT;     // F · P · Fᵀ
Pminus      = temp_matrix + Q; // + Q
```

### Eq3: 卡尔曼增益 — `SetGain()`

```
K = P⁻ · Hᵀ · (H · P⁻ · Hᵀ + R)⁻¹
```

```cpp
HT          = Hᵀ;
temp_matrix  = H * Pminus;               // H · P⁻
temp_matrix1 = temp_matrix * HT;         // H · P⁻ · Hᵀ
S            = temp_matrix1 + R;         // S = H·P⁻·Hᵀ + R
temp_matrix1 = S⁻¹;                       // 求逆 (Gauss-Jordan)
temp_matrix  = Pminus * HT;              // P⁻ · Hᵀ
K            = temp_matrix * temp_matrix1; // K = P⁻·Hᵀ · S⁻¹
```

### Eq4: 后验状态估计 — `UpdateXhat()`

```
x̂(k|k) = x̂(k|k-1) + K · (z - H · x̂(k|k-1))
```

```cpp
temp_vector  = H * xhatminus;       // H · x̂⁻
temp_vector1 = z - temp_vector;     // 新息/残差: z - H·x̂⁻
temp_vector  = K * temp_vector1;    // K · (新息)
xhat         = xhatminus + temp_vector;  // x̂⁺ = x̂⁻ + K·新息
```

### Eq5: 后验协方差估计 — `UpdateCovariance()`

```
P(k|k) = P(k|k-1) - K · H · P(k|k-1)
```

> 等价于 $(I - K \cdot H) \cdot P_{k\vert k-1}$，是标准（非 Joseph）形式。

```cpp
temp_matrix  = K * H;                // K · H
temp_matrix1 = temp_matrix * Pminus; // K · H · P⁻
P            = Pminus - temp_matrix1; // P⁺ = P⁻ - K·H·P⁻
```

---

## 完整更新流程 — `Update()`

```
Update()
├── Measure()                      // 读入观测 z 和控制 u
│   ├── (UseAutoAdjustment=true)
│   │   └── AdjustMeasurementMatrices()  // 自动筛选有效观测，压缩 H/R/z
│   └── (UseAutoAdjustment=false)
│       └── memcpy(z, MeasuredVector)
│
├── [Hook 0]  User_Func0_f        // 用户扩展点
├── UpdateXhatMinus()              // Eq1: 先验状态 x̂⁻
├── [Hook 1]  User_Func1_f        // 用户扩展点
├── UpdatePminus()                 // Eq2: 先验协方差 P⁻
├── [Hook 2]  User_Func2_f        // 用户扩展点
│
├── if (有有效观测)
│   ├── SetGain()                  // Eq3: 卡尔曼增益 K
│   ├── [Hook 3]  User_Func3_f
│   ├── UpdateXhat()               // Eq4: 后验状态 x̂
│   ├── [Hook 4]  User_Func4_f
│   └── UpdateCovariance()         // Eq5: 后验协方差 P
└── else (无有效观测, 自动模式下)
    └── x̂ = x̂⁻, P = P⁻            // 只保留预测

├── [Hook 5]  User_Func5_f
├── 协方差下限钳位: P[i][i] = max(P[i][i], StateMinVariance[i])
├── FilteredValue = xhat_data      // 输出滤波结果
├── [Hook 6]  User_Func6_f
└── return FilteredValue
```

---

## 自动观测调节模式

当 `UseAutoAdjustment = true` 时，`Measure()` 内部调用 `AdjustMeasurementMatrices()`：

1. **遍历** `MeasuredVector`，非零项视为有效观测
2. **压缩** 有效观测到 `z_data` 前段
3. **重建 H**：`H[validIdx][MeasurementMap[i] - 1] = MeasurementDegree[i]`
4. **重建 R**：`R[validIdx][validIdx] = MatR_DiagonalElements[i]`
5. **收缩** `H/R/K/z` 的行数到 `MeasurementValidNum`

> 适用场景：多传感器融合，部分传感器暂时不可用（值为 0）时自动降维。

---

## 设计要点

- **零堆分配**：所有矩阵数据存储在 `KalmanFilter` 对象内部的静态数组（`_storage_` 系列），通过 `BindStorage()` 绑定指针，适合裸机/RTOS 环境。
- **矩阵视图**：`Matrix` 结构不持有内存，只记录行列数和指针。
- **Gauss-Jordan 求逆**：`MatrixInverse()` 使用列主元消去，最小值判断阈值 `kPivotEpsilon = 1e-8`。
- **跳过标志**：`SkipEq1`~`SkipEq5` 配合 7 个 Hook 点 (`User_Func0_f`~`User_Func6_f`)，允许外部改为 EKF 等非线性实现。
- **协方差下限钳位**：`StateMinVariance` 防止滤波发散或锁定。
- **最大维度常量**：`kMaxStateSize=16`, `kMaxControlSize=8`, `kMaxMeasureSize=8`。
