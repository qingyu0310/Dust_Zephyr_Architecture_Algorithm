# PowerCtrl — 电机功率控制器

功率模型预测 → RLS 在线参数辨识 → 隶属度分配 → 受限力矩求解

## 概述

PowerCtrl 管理一组同类型电机（最多 4 个），基于功率模型预测各电机功耗，通过 RLS 在线辨识电机损耗参数，在总功率预算内做最优化分配并求解受限电流。主要用于 RoboMaster 等需要严格功率限制的场景。

## 依赖

- `alg::rls::RLS`（[identify/rls](../identify/rls/README.md)）
- `math.h`（sqrtf）

## 功率模型

```
P_in = K1 · τ² + K2 · ω² + τ · ω + K3
```

| 项 | 物理含义 | 说明 |
|----|----------|------|
| K1·τ² | 铜损 | 绕组电阻发热，与转矩平方成正比 |
| K2·ω² | 铁损 + 摩擦 | 磁滞/涡流损耗 + 轴承摩擦，与转速平方成正比 |
| τ·ω | 机械功率 | 实际对外做功，P_mech = τ·ω |
| K3 | 常数损耗 | 控制器待机功耗、传感器供电等固定损耗 |

**辨识目标**：RLS 以 `x = [Στ², Σω²]` 为输入、`y = P_meas - K3` 为期望输出，在线拟合 K1（铜损系数）和 K2（铁损系数）。`τ·ω` 和 K3 是已知或可直接计算的项，不作为辨识参数。

## 管线

```
每个控制周期 (1ms)

① SetMotorData(i, current, velocity, pidErr)
   ├── τ = current × torqueK
   ├── ω = velocity / 9.55      (RPM → rad/s)
   ├── τ², ω² 缓存
   └── pidErr 缓存（隶属度计算用）

② Predict()
   ├── P_pred[i] = K1·τ² + K2·ω² + τ·ω + K3
   ├── Στ², Σω²
   └── RLS Update(x=[Στ², Σω²], y=P_meas - K3)
       └── 更新 K1、K2

③ Allocate(budget)
   ├── 权重 K = 误差水平插值 (errLower~errUpper)
   ├── membership[i] = K·|err|/Σ|err| + (1-K)·|P_pred|/Σ|P_pred|
   ├── powerLimit[i] = membership[i] × budget
   ├── ΣP_pred ≤ budget → 电流直通
   └── ΣP_pred > budget → 解 K1·τ² + ω·τ + (K2·ω² + K3 - P_limit) = 0

④ GetLimitedCurrent(i) → 功率限制后的电流值
```

## 使用方式

### 实例化

两种构造方式：

```cpp
// 方式一：默认构造 + Init
alg::power_ctrl::PowerCtrl ctrl;
alg::power_ctrl::PowerCtrl::Config cfg{};
cfg.motorCount = 1;
cfg.rlsEnable = true;      // 开启 RLS 在线辨识
ctrl.Init(cfg);

// 方式二：构造时传入 Config
alg::power_ctrl::PowerCtrl ctrl({
    4,              // motorCount
    4.577e-5f,      // torqueK
    1.453e-7f,      // k1Init
    1.453e-7f,      // k2Init
    3.0f,           // k3
    50.0f,          // errUpper
    0.01f,          // errLower
    0.99999f,       // rlsLambda
    false           // rlsEnable
});
```

### Config 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| motorCount | 4 | 电机数量（≤ 4） |
| torqueK | 4.577e-5 | 电流 → 转矩系数（N·m/A），M3508 常数 |
| k1Init | 1.453e-7 | 铜损系数 K1 初始值 |
| k2Init | 1.453e-7 | 铁损系数 K2 初始值 |
| k3 | 3.0 | 常数损耗（W），从 RLS 目标中减掉 |
| errUpper | 50.0 | 隶属度上限阈值（|err| ≥ 此值时全按需求分配） |
| errLower | 0.01 | 隶属度下限阈值（|err| ≤ 此值时全按功率分配） |
| rlsLambda | 0.99999 | RLS 遗忘因子 |
| rlsEnable | false | RLS 在线辨识开关 |

### 控制循环

```cpp
void control_task(void*, void*, void*)
{
    for (;;)
    {
        /* ① 喂数据 */
        ctrl.SetMotorData(0, motor.GetNowCurrent(), motor.GetNowVelocity(), pidErr);

        /* ② 功率预测 + RLS */
        ctrl.SetMeasuredPower(powerMeter.GetPower());
        ctrl.Predict();

        /* ③ 功率分配 */
        ctrl.Allocate(45.0f);

        /* ④ 取限制后的电流 */
        float limited = ctrl.GetLimitedCurrent(0);

        k_msleep(1);
    }
}
```

### 辨识结果观察

开启 RLS 后，Predict() 内部自动更新 K1、K2：

```cpp
printf("K1=%e, K2=%e, Pred=%f, Meas=%f\n",
       ctrl.GetK1(), ctrl.GetK2(),
       ctrl.GetTotalPower(), measuredPower);
```

K1/K2 会从初始值逐渐收敛到真实值。收敛速度取决于 λ 和激励的丰富程度（转速/转矩变化范围）。

## RLS 在 PowerCtrl 中的工作方式

### 数据流

```
           ┌─────────────────┐
           │  SetMotorData   │
           │  SetMeasuredPower│
           └────────┬────────┘
                    ↓
           ┌─────────────────┐
           │   Predict()     │
           │                  │
           │  P_pred = K1·τ² │
           │        + K2·ω²  │
           │        + τ·ω    │
           │        + K3     │
           │                  │
           │  x = [Στ², Σω²] │
           │  y = P_meas-K3  │
           │       ↓         │
           │  RLS.Update(x,y)│
           │  K1=w[0], K2=w[1]│
           └─────────────────┘
```

### 为什么减去 K3

RLS 建模的是 `K1·τ² + K2·ω² = P_meas - K3 - τ·ω`。如果不减 K3，K3 作为常数偏置会吸收到 K1 或 K2 中，导致参数偏移。剩下 `τ·ω` 是真实的机械功率，不应被 RLS 吃掉，所以保留在残差中。

### 激励条件

RLS 需要 **持续激励** 才能收敛。具体来说，输入向量 `x = [τ², ω²]` 的两个分量需要独立变化。如果电机始终匀速（ω² 恒定）且转矩极小（τ² ≈ 0），K1 和 K2 无法同时辨识。建议的做法：

1. 在电机起停、变速、加减载过程中收集数据
2. 速度变化范围大 → K2 收敛快
3. 转矩变化范围大 → K1 收敛快
4. 固定 K3 作为已知量，减少辨识维度

## 关键边界条件

- **不超限**：总预测功率未超过预算时，电流直通，不做限制
- **超限求解**：对每个电机解 `K1·τ² + ω·τ + (K2·ω² + K3 - P_limit) = 0`，选择与原始力矩同号的根
- **电机数 0**：Allocate 直接返回
- **index 越界**：SetMotorData 和 GetLimitedCurrent 静默返回
- **delta < 0**：二次方程判别式负数时归零

## 典型配置示例

### 单电机辨识

```cpp
alg::power_ctrl::PowerCtrl ctrl;
alg::power_ctrl::PowerCtrl::Config cfg{};
cfg.motorCount = 1;
cfg.rlsEnable = true;
ctrl.Init(cfg);

// 循环中:
ctrl.SetMotorData(0, current, rpm, 0.0f);
ctrl.SetMeasuredPower(powerMeter.GetPower());
ctrl.Predict();
// K1, K2 逐步收敛
```

### 四电机底盘功率限制

```cpp
alg::power_ctrl::PowerCtrl chassisCtrl({
    4, 4.577e-5f, 1.453e-7f, 1.453e-7f, 3.5f,
    0.001f, 500.0f, 0.99999f, false   // 不开启 RLS，用固定 K1/K2
});

// 循环中:
for (int i = 0; i < 4; i++)
    chassisCtrl.SetMotorData(i, current[i], rpm[i], pidErr[i]);
chassisCtrl.Predict();
chassisCtrl.Allocate(45.0f);     // 裁判系统 45W
for (int i = 0; i < 4; i++)
    limitedCurrent[i] = chassisCtrl.GetLimitedCurrent(i);
```
