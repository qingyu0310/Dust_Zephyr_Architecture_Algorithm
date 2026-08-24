# algorithm/ 架构说明

> 更新日期：2026-08-21
> 依据文件：`framework/algorithm/Kconfig`、`framework/algorithm/CMakeLists.txt` 和当前源码目录。

## 职责

`framework/algorithm` 是纯计算层，提供缓冲区、控制器、滤波器、参数辨识和线性代数能力。它只接收数据、更新状态、输出结果，不直接操作设备树节点、Zephyr device 句柄、线程对象或 topic。

当前源码包含：

- `buffer/`：`BipBuffer` 双区连续缓冲、`RingBuf` FIFO。
- `controller/`：`Pid`、`PowerCtrl<N>`、协作式 `Timer`。
- `filter/`：一阶 LPF/HPF、标准 `Kalman<N,M>`、`ExtendedKalman<N,M,U>`、`QuaternionEkf`。
- `identify/`：`RLS<N>`、`MotorPlant`、`stability::WinStable` 和波形发生器。
- `math/eigen/`：内置 Eigen header-only 线性代数库。

## 边界

| 管 | 不管 |
| --- | --- |
| 数学计算、状态更新、配置参数、固定容量缓冲 | 不创建 `k_thread` |
| 纯数据输入输出、模板化静态内存 | 不持有 GPIO/PWM/UART/CAN/USB 等硬件句柄 |
| Kconfig 可裁剪的算法能力 | 不发布 topic、不解析设备协议、不编排业务流程 |

`Timer` 是协作式计数器，由外部周期调用 `Update()` 驱动；它不是硬件 timer。`DUST_CTL_EXECTIMER` 当前在 Kconfig/CMake 中保留，但当前源码未提供独立 `ExecTimer` 类，使用前需要补实现或清理该开关。

## 目录结构

```text
algorithm/
├── buffer/
│   ├── bipbuf/bipbuf.hpp
│   └── ringbuf/ringbuf.hpp
├── controller/
│   ├── pid/pid.hpp + pid.cpp
│   ├── power_ctrl/power_ctrl.hpp
│   └── timer/timer.hpp
├── filter/
│   ├── hpf/hpf.hpp
│   ├── lpf/lpf.hpp
│   ├── kalman/kalman.hpp + kalman_ekf.hpp + kalman_check.cpp
│   └── quaternion/quaternion.hpp + quaternion.cpp
├── identify/
│   ├── rls/rls.hpp
│   ├── motor/motorplant.hpp + motorplant.cpp
│   └── stability.hpp
├── math/eigen/
├── CMakeLists.txt
└── Kconfig
```

## 装配关系

所有算法开关默认 `n`，由上层业务开关通过 `select DUST_*` 拉起。`CMakeLists.txt` 按 `CONFIG_DUST_*` 把对应 include 目录和少量 `.cpp` 加入 `app`。

关键依赖以当前 `Kconfig` 为准：

| 符号 | 当前装配 |
| --- | --- |
| `DUST_FLT_KALMAN` | select `DUST_MATH_EIGEN`，编译 `kalman_check.cpp` |
| `DUST_FLT_KALMAN_EKF` | select `DUST_MATH_EIGEN`，共用 `kalman_check.cpp` |
| `DUST_FLT_QUATERNION` | select `DUST_FLT_KALMAN_EKF`，编译 `quaternion.cpp` |
| `DUST_ID_RLS` | select `DUST_MATH_EIGEN` |
| `DUST_ID_MOTOR_PLANT` | select `DUST_ID_RLS`，编译 `motorplant.cpp` |
| `DUST_MOD_CTL_POWER` | select `DUST_ID_RLS` 和 `DUST_FLT_LPF`；消耗 PID 目标/误差，但不实例化 PID，也不 select `DUST_CTL_PID` |

## 调用方

算法类通常由 `project/thread/` 的业务线程或 `framework/modules/` 的设备模块实例化。调用方负责提供采样周期、传感器数据、目标值、互斥保护和线程调度；算法层只保证计算语义，不保证调用时序。
