/**
 * @file eso.hpp
 * @author qingyu
 * @brief 模板化 N 阶线性扩张状态观测器（LESO）—— 带宽参数化
 * @version 0.2
 * @date 2026-10-05
 *
 * @par 算法原理
 *      被控对象为 kN-1 阶:  y^(kN-1) = b0·u + f
 *      将总扰动 f 扩张为第 kN 个状态，观测器状态 z = [z1, z2, ..., z_kN]:
 *      - z1 ~ y, z2 ~ y', ..., z_{kN-1} ~ y^(kN-2)   被控量及其各阶导估计
 *      - z_kN ~ f                                     总扰动估计
 *
 *      连续域观测器方程:
 *      ż1     = z2     - β1·(z1 - y)
 *      ż2     = z3     - β2·(z1 - y)
 *      ...
 *      ż_kN-1 = z_kN   - β_kN-1·(z1 - y) + b0·u
 *      ż_kN   =          - β_kN·(z1 - y)
 *
 *      带宽参数化（Gao）: 全部极点配置到 -ωo，β_i = C(kN, i)·ωo^i，
 *      只需一个观测器带宽 ωo，避免逐项调参。
 *
 * @tparam kN  ESO 阶数（状态维数）= 被控对象阶数 + 1（一阶对象取 kN=2）
 *
 * @par 使用示例
 * @code
 *   // 一阶对象（速度环）: ω' = b0·u + f
 *   alg::observer::Eso<2> eso(100.0f, 0.001f, 50.0f);  // ωo=100 rad/s, dt=1 ms, b0=50
 *
 *   // 每控制周期
 *   eso.Update(omega_meas, u);              // 喂入测量值 + 控制量
 *   float omega_est = eso.GetEstimate();    // z1: 速度估计
 *   float f_est     = eso.GetDisturbance(); // z2: 总扰动估计
 * @endcode
 */

#pragma once

#pragma message "Compiling Algorithm/Observer/Eso"

#include <stdint.h>

namespace alg::observer {

/**
 * @brief N 阶线性扩张状态观测器（LESO）
 *
 *        离散实现（前向欧拉）:  z += dt·ż，Update 每个采样周期调用一次
 *        带宽参数化: 全部极点配置到 -ωo，β_i = C(kN, i)·ωo^i
 */
template <uint8_t kN>
class Eso final
{
public:
    static_assert(kN >= 2, "ESO 阶数至少为 2（被控对象至少一阶）");

    Eso() = default;

    /**
     * @brief 构造并初始化
     * @param wo  观测器带宽 ωo (rad/s)，越大收敛越快、抗噪越差
     * @param dt  采样周期 (s)
     * @param b0  控制增益（对象 y^(kN-1) = b0·u + f 中的 b0）
     */
    Eso(float wo, float dt, float b0 = 0.0f) { Init(wo, dt, b0); }

    /**
     * @brief（重新）配置观测器：按带宽算增益并清零状态
     * @param wo  观测器带宽 ωo (rad/s)
     * @param dt  采样周期 (s)
     * @param b0  控制增益
     */
    void Init(float wo, float dt, float b0 = 0.0f)
    {
        wo_ = wo;
        dt_ = dt;
        b0_ = b0;

        // β_i = C(kN, i)·ωo^i（i = 1..kN），把 (A − L·C) 的极点全配到 −ωo
        float woPow = 1.0f;
        for (uint8_t i = 0; i < kN; i++) {
            woPow *= wo;
            beta_[i] = static_cast<float>(Binomial(kN, i + 1)) * woPow;
        }

        for (uint8_t i = 0; i < kN; i++) z_[i] = 0.0f;
    }

    /**
     * @brief 执行一步观测更新
     * @param y  被控量测量值
     * @param u  本周期作用于对象的控制量
     * @return   总扰动估计 z_kN
     */
    float Update(float y, float u)
    {
        const float e = z_[0] - y;                  // 观测误差 z1 − y

        float dz[kN];
        for (uint8_t i = 0; i + 1 < kN; i++) {
            dz[i] = z_[i + 1] - beta_[i] * e;
        }
        dz[kN - 2] += b0_ * u;                      // 控制量注入 y^(kN-2) 通道
        dz[kN - 1]  = -beta_[kN - 1] * e;           // 总扰动通道

        for (uint8_t i = 0; i < kN; i++) {
            z_[i] += dt_ * dz[i];
        }

        return z_[kN - 1];
    }

    /**
     * @brief 第 i 个状态估计 z_{i+1}（0 基，i=0 为被控量估计）
	 *
     */
    float GetState(uint8_t i) const { return (i < kN) ? z_[i] : 0.0f; }
	
    /**
     * @brief 被控量估计 z1
	 *
     */
    float GetEstimate() 	  const { return z_[0]; }

    /**
     * @brief 总扰动估计 z_kN
	 *
     */
    float GetDisturbance() 	  const { return z_[kN - 1]; }

    /**
     * @brief 运行时修改控制增益 b0
     * @param b0  控制增益
     */
    void SetB0(float b0) { b0_ = b0; }

    /**
     * @brief 清零全部状态（保留带宽/周期/b0 配置）
     * 
     */
    void Reset()
    {
        for (uint8_t i = 0; i < kN; i++) z_[i] = 0.0f;
    }

private:
    float wo_ = 0.0f;                   // 观测器带宽 ωo (rad/s)
    float dt_ = 0.0f;                   // 采样周期 (s)
    float b0_ = 0.0f;                   // 控制增益

    float z_[kN]    {};                 // 状态 z1..z_kN（末位为总扰动估计）
    float beta_[kN] {};                 // 观测器增益 β1..β_kN（带宽参数化）

    /**
     * @brief 二项式系数 C(n, k)
     * @param n  上标
     * @param k  下标
     * @return   C(n, k)，k > n 时为 0
     */
    static uint32_t Binomial(uint32_t n, uint32_t k)
    {
        if (k > n) return 0;

        uint32_t r = 1;
        for (uint32_t i = 1; i <= k; i++) {
            r = r * (n - k + i) / i;
        }
        return r;
    }
};

} // namespace alg::observer
