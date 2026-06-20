# Dust_Zephyr_Architecture_Algorithm

纯计算模块。提供控制器、滤波器、参数辨识等算法封装。

- **controller/** — PID 控制器、功率控制器、软件定时器
- **filter/** — 低通滤波、高通滤波、卡尔曼滤波、四元数姿态 EKF
- **identify/** — 递归最小二乘法（RLS）
- **tflm/** — TensorFlow Lite Micro 移植

所有算法零硬件依赖，仅数学计算，可独立于框架测试。

详见 [ARCHITECTURE.md](ARCHITECTURE.md)。
