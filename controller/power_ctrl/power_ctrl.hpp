/**
 * @file power_ctrl.hpp
 * @author qingyu
 * @brief 电机功率控制器 — 模型预测 / RLS 在线辨识 / 隶属度分配 / 力矩限幅
 * @version 0.2
 * @date 2026-06-27
 */

#pragma once

#pragma message "Compiling Algorithm/Controller/PowerCtrl"

#include <stdint.h>
#include <math.h>
#include "rls.hpp"

namespace alg::power_ctrl {

/**
 * @brief  单电机状态
 */
struct MotorState {
    float torque        = 0.0f;         // τ = current × torqueK
    float omega         = 0.0f;         // ω = rpm / 9.55 (rad/s)
    float torque2       = 0.0f;         // τ²
    float omega2        = 0.0f;         // ω²
    float powerPred     = 0.0f;         // P_in = K1·τ² + K2·ω² + τ·ω + K3
    float pidErr        = 0.0f;         // PID 误差（隶属度计算用）
    float powerLimit    = 0.0f;         // 该电机分配的功率预算
    float targetCurrent = 0.0f;         // PID 目标电流（串级输入）
    float currentOut    = 0.0f;         // 限幅后的输出电流
};

/**
 * @brief  单组电机功率控制器（如底盘 4 个电机）
 *
 * @tparam kMotorCount  电机数量（编译期分配）
 *
 * @par 管线（每控制周期调用一次）：
 *       @n 1. SetMotorData(i, current, velocity, pidErr)   — 喂入各电机数据
 *       @n 2. Predict()                                    — 功率预测 + RLS
 *       @n 3. Allocate(totalBudget)                        — 隶属度分配 + 限幅
 *       @n 4. GetLimitedCurrent(i)                         — 读取结果
 *
 * @par RLS：
 *       在线辨识 K1（铜损系数）和 K2（铁损系数）。
 *       K3 为固定常数损耗，不参与辨识。
 */
template <uint8_t kMotorCount>
class PowerCtrl final
{
public:
    struct Config {
        float   torqueK    = 4.577e-5f;         // 电流 → 转矩（M3508 常数）
        float   k1Init     = 1.453e-7f;         // RLS 初始 K1（铜损系数）
        float   k2Init     = 1.453e-7f;         // RLS 初始 K2（铁损系数）
        float   k3         = 3.0f;              // 固定常数损耗
        float   errUpper   = 50.0f;             // 隶属度：全需求阈值
        float   errLower   = 0.01f;             // 隶属度：全功率阈值
        float   rlsLambda  = 0.99999f;          // RLS 遗忘因子
        bool    rlsEnable  = false;             // false = 使用固定 K1/K2
        bool    tauOmegaEnable = true;          // 是否包含 τ·ω 项
    };

    PowerCtrl() : PowerCtrl(Config{}) {}

    /**
     * @brief 构造功率控制器，初始化 RLS
     * @param cfg 配置参数
     */
    explicit PowerCtrl(const Config& cfg)
        : cfg_(cfg)
        , k1_(cfg.k1Init)
        , k2_(cfg.k2Init)
        , rls_(cfg.rlsLambda, 1e-5f)
    {
        rlsInited_ = cfg.rlsEnable;

        if (rlsInited_) {
            float w[2] = { k1_, k2_ };
            rls_.SetWeights(w);
        }
    }

    /**
     * @brief 重配置功率控制器，清空状态
     * @param cfg 配置参数
     */
    void Init(const Config& cfg)
    {
        cfg_ = cfg;
        k1_  = cfg.k1Init;
        k2_  = cfg.k2Init;
        rlsInited_ = cfg.rlsEnable;

        if (rlsInited_) {
            float w[2] = { k1_, k2_ };
            rls_.SetWeights(w);
        }

        for (auto& m : motors_) m = MotorState{};
        for (auto& m : membership_) m = 0.0f;

        measuredPower_ = 0.0f;
    }

    /**
     * @brief 喂入单电机原始数据
     * @param idx     电机索引 [0, kMotorCount)
     * @param torque  实测转矩 (N·m)
     * @param omega   角速度 (rad/s)
     * @param pidErr  PID 误差（隶属度计算用）
     */
    void SetMotorData(uint8_t idx, float torque, float omega, float pidErr)
    {
        if (idx >= kMotorCount) return;

        auto& m = motors_[idx];
        m.torque   = torque;
        m.omega    = omega;
        m.torque2  = m.torque * m.torque;
        m.omega2   = m.omega  * m.omega;
        m.pidErr   = pidErr;
    }

    /**
     * @brief 设置 PID 目标电流（串级输入）
     * @param idx      电机索引
     * @param current  目标电流
     */
    void SetTarget(uint8_t idx, float current)
    {
        if (idx >= kMotorCount) return;
        motors_[idx].targetCurrent = current;
    }

    void SetMeasuredPower(float power) { measuredPower_ = power; }

    /**
     * @brief 功率预测 + RLS 在线更新
     *
     * 对每台电机计算 P = K1·τ² + K2·ω² + τ·ω(可选) + K3，
     * 再以 Στ²、Σω² 为输入，P_meas - K3 为期望输出，
     * 更新 K1（铜损系数）和 K2（铁损系数）。
     */
    void Predict()
    {
        float sumTorque2 = 0.0f;
        float sumOmega2  = 0.0f;

        for (uint8_t i = 0; i < kMotorCount; i++)
        {
            auto& m = motors_[i];

            // P = K1·τ² + K2·ω² + τ·ω(可选) + K3
            m.powerPred = k1_ * m.torque2 + k2_ * m.omega2 + cfg_.k3;
            if (cfg_.tauOmegaEnable) {
                m.powerPred += m.torque * m.omega;
            }

            sumTorque2 += m.torque2;
            sumOmega2  += m.omega2;
        }

        // RLS 在线更新 — 减掉已知的 K3，RLS 只管 K1·τ² + K2·ω²
        if (rlsInited_)
        {
            float x[2] = { sumTorque2, sumOmega2 };
            rls_.Update(x, measuredPower_ - cfg_.k3);
            k1_ = rls_.GetWeights()[0];
            k2_ = rls_.GetWeights()[1];
        }
    }

    /**
     * @brief 功率分配：隶属度计算 + 力矩限幅
     *
     * 总预测未超预算时直通 PID 目标电流；
     * 超限时逐电机解二次方程求受限转矩。
     * @param totalBudget 该组总功率预算 (W)
     */
    void Allocate(float totalBudget)
    {
        float sumAbsErr   = 0.0f;
        float sumPowerAbs = 0.0f;

        for (uint8_t i = 0; i < kMotorCount; i++) {
            sumAbsErr   += fabsf(motors_[i].pidErr);
            sumPowerAbs += fabsf(motors_[i].powerPred);
        }

        // 权重 K（误差水平插值）
        float k;
        if (sumAbsErr >= cfg_.errUpper) {
            k = 1.0f;
        } else if (sumAbsErr <= cfg_.errLower) {
            k = 0.0f;
        } else {
            float range = cfg_.errUpper - cfg_.errLower;
            k = (sumAbsErr - cfg_.errLower) / (range > 0.0f ? range : 1.0f);
        }

        // 隶属度 + 功率上限
        for (uint8_t i = 0; i < kMotorCount; i++)
        {
            float ratioErr   = (sumAbsErr   > 0.0f) ? fabsf(motors_[i].pidErr)     / sumAbsErr   : 0.0f;
            float ratioPower = (sumPowerAbs > 0.0f) ? fabsf(motors_[i].powerPred)  / sumPowerAbs : 0.0f;

            membership_[i] = k * ratioErr + (1.0f - k) * ratioPower;

            if (membership_[i] < 0.0f) membership_[i] = 0.0f;
            if (membership_[i] > 1.0f) membership_[i] = 1.0f;

            motors_[i].powerLimit = membership_[i] * totalBudget;
        }

        // 受限力矩求解
        float totalPred = GetTotalPower();

        if (totalPred <= totalBudget) {
            // 未超限：透传 PID 目标电流
            for (uint8_t i = 0; i < kMotorCount; i++) {
                motors_[i].currentOut = motors_[i].targetCurrent;
            }
            return;
        }

        // 超限：逐电机解方程，限制 PID 目标
        for (uint8_t i = 0; i < kMotorCount; i++)
        {
            const auto& m = motors_[i];

            float a = k1_;                                          // K1
            float b = m.omega;                                      // ω
            float c = k2_ * m.omega2 + cfg_.k3 - m.powerLimit;      // K2·ω² + K3 - P_limit

            bool positive = (m.targetCurrent >= 0.0f);
            float limitedTorque = SolveTorque(a, b, c, positive);

            float limitedCurrent = limitedTorque / cfg_.torqueK;
            motors_[i].currentOut = (fabsf(m.targetCurrent) <= fabsf(limitedCurrent)) ? m.targetCurrent : limitedCurrent;
        }
    }

    /**
     * @brief 取限幅后的输出电流
     * @param idx  电机索引
     * @return 限幅后的电流值，越界返回 0
     */
    float GetLimitedCurrent(uint8_t idx) const
    {
        return (idx < kMotorCount) ? motors_[idx].currentOut : 0.0f;
    }

    /**
     * @brief 所有电机总预测功率
     * @return 总功率 (W)
     */
    float GetTotalPower() const
    {
        float sum = 0.0f;
        for (uint8_t i = 0; i < kMotorCount; i++) {
            sum += motors_[i].powerPred;
        }
        return sum;
    }

    float GetK1() const { return k1_; }
    float GetK2() const { return k2_; }

private:
    Config cfg_;

    float k1_ = 0.0f;                           // 铜损系数
    float k2_ = 0.0f;                           // 铁损系数

    MotorState motors_[kMotorCount]{};          // 各电机状态（转矩/转速/功率等）
    float      membership_[kMotorCount]{};      // 各电机隶属度（功率分配系数 0~1）

    alg::rls::RLS<2> rls_;                      // RLS 在线辨识器
    bool rlsInited_ = false;

    float measuredPower_ = 0.0f;                // 功率计实测总功率

    /**
     * @brief 解二次方程 A·τ² + B·τ + C = 0，取与目标同号的根
     */
    static float SolveTorque(float a, float b, float c, bool positive)
    {
        if (a == 0.0f) {
            return (b != 0.0f) ? -c / b : 0.0f;
        }

        float delta = b * b - 4.0f * a * c;
        if (delta < 0.0f) delta = 0.0f;

        float sqrtDelta = sqrtf(delta);
        float t1 = (-b + sqrtDelta) / (2.0f * a);
        float t2 = (-b - sqrtDelta) / (2.0f * a);

        return positive ? t1 : t2;
    }
};

} // namespace alg::power_ctrl
