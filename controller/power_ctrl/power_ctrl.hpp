/**
 * @file power_ctrl.hpp
 * @author qingyu
 * @brief 电机功率控制器 — 模型预测 / RLS 在线辨识 / 隶属度分配 / 力矩限幅
 * @version 0.3
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
    float absOmega      = 0.0f;         // |ω| (rad/s) — 转速线性损耗回归量
    float powerPred     = 0.0f;         // P_in = k1·τ² + k2·|ω| + Kτ·ω·τ·ω + K3
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
 *       在线辨识 k1（τ² 铜损系数）、k2（转速线性损耗系数，|ω| 项）。
 *       K3（常数损耗）与 Kτ·ω（τ·ω 机械功率系数，默认 1.0）为已知固定项，不参与辨识。
 */
template <uint8_t kMotorCount>
class PowerCtrl final
{
public:
    struct Config {
        float   torqueK    = 4.577e-5f;         // 电流 → 转矩（M3508 常数）
        float   k1Init     = 1.453e-7f;         // RLS 初始 K1（铜损系数）
        float   k2Init     = 1.453e-7f;         // RLS 初始 k2（转速线性损耗系数，|ω| 项）
        float   kTauOmegaInit = 1.0f;           // τ·ω 耦合系数（固定，不参与 RLS 辨识）
        float   k3         = 3.0f;              // 固定常数损耗
        float   errUpper   = 50.0f;             // 隶属度：全需求阈值
        float   errLower   = 0.01f;             // 隶属度：全功率阈值
        float   rlsLambda  = 0.99999f;          // RLS 遗忘因子（改业务层传入 0.999，现在能生效）
        float   pInit      = 1e-5f;             // RLS 协方差初值（归一化后 1e-5 足够，可调）
        float   excMinAbsOmega = 0.0f;          // 激励门控：Σ|ω| 下限（0=不门控）
        float   excMinTau2   = 0.0f;            // 激励门控：Στ² 下限（0=不门控）
        float   powerMax     = 0.0f;            // 总功率钳制上限（0=不钳制；底盘传总预算 60）
        bool    fixK2        = false;           // true=固定 k2(|ω|) 只辨 k1(τ²)（单参数，防共线漂移）
        float   deadzonePower = 0.0f;           // RLS 更新死区：|P_meas|<此值跳过（0=不启用；港科大 5W）
        float   kFloor       = 1e-5f;           // 辨识系数下限钳位（防发散为负）
        bool    skipNegPower = false;           // true=预测功率为负（再生制动/停车）时跳过 RLS 更新
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
        , kTauOmega_(cfg.kTauOmegaInit)
        , rls_(cfg.rlsLambda, cfg.pInit)
        , rls1_(cfg.rlsLambda, cfg.pInit)
    {
        rlsInited_ = cfg.rlsEnable;

        if (rlsInited_) {
            float w[2] = { k1_, k2_ };                     // 去归一化：直接物理值
            rls_.SetWeights(w);
            float w1[1] = { k1_ };
            rls1_.SetWeights(w1);
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
        kTauOmega_ = cfg.kTauOmegaInit;
        rlsInited_ = cfg.rlsEnable;

        if (rlsInited_) {
            rls_.Reset(cfg.pInit);                          // 重置协方差 P 与权重（顺序：先 Reset）
            rls_.SetLambda(cfg.rlsLambda);                  // λ 生效（修复死配置）
            float w[2] = { k1_, k2_ };                      // 去归一化：直接物理值
            rls_.SetWeights(w);
            rls1_.Reset(cfg.pInit);                         // 单参数辨识器同参初始化
            rls1_.SetLambda(cfg.rlsLambda);
            float w1[1] = { k1_ };
            rls1_.SetWeights(w1);
        }

        for (auto& m : motors_) m = MotorState{};
        for (auto& m : membership_) m = 0.0f;

        measuredPower_ = 0.0f;
        powerValid_    = false;                             // 首帧保护：待功率计首帧
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
        m.absOmega = fabsf(m.omega);
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

    void SetMeasuredPower(float power, bool valid)
    {
        measuredPower_ = power;
        powerValid_    = valid;
    }

    /**
     * @brief 功率预测 + RLS 在线更新
     *
     * 对每台电机计算 P = k1·τ² + k2·|ω| + Kτ·ω·τ·ω(可选) + K3，
     * 再以 Στ²、Σ|ω| 为输入，(P_meas - K3 - Kτ·ω·Στ·ω) 为期望输出，
     * 更新 k1（τ² 铜损系数）、k2（转速线性损耗系数）。Kτ·ω 为固定系数，不参与辨识。
     */
    void Predict()
    {
        sumTorque2_  = 0.0f;
        sumAbsOmega_ = 0.0f;
        sumTauOmega_ = 0.0f;

        for (uint8_t i = 0; i < kMotorCount; i++)
        {
            auto& m = motors_[i];

            // P = k1·τ² + k2·|ω| + Kτ·ω·τ·ω(可选) + K3
            m.powerPred = k1_ * m.torque2 + k2_ * m.absOmega + cfg_.k3 / kMotorCount;   // K3 整车常数，按电机数分摊
            if (cfg_.tauOmegaEnable) {
                m.powerPred += kTauOmega_ * m.torque * m.omega;
                sumTauOmega_ += m.torque * m.omega;
            }

            sumTorque2_ += m.torque2;
            sumAbsOmega_ += m.absOmega;
        }

        // RLS 在线更新 — Kτ·ω 为已知固定系数（默认 1.0），只辨识 k1·τ² + k2·|ω|
        // 首帧保护：功率计有效帧到达前不更新（防 target=-K3-Στ·ω<0 清零窗口）
        // 三闸门（港科大防漂移）：
        //   ① 激励门控：Σ|ω|/Στ² 任一达标才更新（空载 P 不膨胀）
        //   ② 死区：|P_meas| < deadzonePower 跳过（小功率段 SNR 差，残差是噪声）
        //   ③ 负功率门控：skipNegPower 且 预测功率<0（再生制动/停车）跳过
        //      —— 制动时模型预测负功率、功率计读正（再生被吸收），残差被机械功率失配主导，
        //         会把 K1 顶飞（港科大用 estimatedPower<0 门控，不是 P_meas<0）
        if (rlsInited_ && powerValid_ &&
            (sumAbsOmega_ >= cfg_.excMinAbsOmega || sumTorque2_ >= cfg_.excMinTau2) &&
            fabsf(measuredPower_) >= cfg_.deadzonePower &&
            !(cfg_.skipNegPower && GetTotalPower() < 0.0f))
        {
            float target = measuredPower_ - cfg_.k3;
            if (cfg_.tauOmegaEnable) {
                target -= kTauOmega_ * sumTauOmega_;  // 减掉已知机械功率，残差 = 铜损 + 转速损耗
            }

            if (cfg_.fixK2) {
                // 单参数模式：固定 k2(|ω|)，只辨 k1(τ²)
                // target 减掉固定 k2 转速损耗，回归量只剩 Στ²，k1 无共线干扰
                target -= cfg_.k2Init * sumAbsOmega_;
                float x1[1] = { sumTorque2_ };
                rls1_.Update(x1, target);
                k1_ = fmaxf(rls1_.GetWeights()[0], cfg_.kFloor);   // 下限钳位防发散为负
                k2_ = cfg_.k2Init;                                  // 固定 k2(|ω|)
            } else {
                // 双参数模式：辨识 k1(τ²)、k2(|ω|)，无归一化
                float x2[2] = { sumTorque2_, sumAbsOmega_ };
                rls_.Update(x2, target);
                k1_ = fmaxf(rls_.GetWeights()[0], cfg_.kFloor);
                k2_ = fmaxf(rls_.GetWeights()[1], cfg_.kFloor);
            }
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

            float a = k1_;                                          				// k1（τ² 系数）
            float b = cfg_.tauOmegaEnable ? (kTauOmega_ * m.omega) : 0.0f;   		// Kτ·ω·ω（默认 1.0；禁用时无此项）
            float c = k2_ * m.absOmega + cfg_.k3 / kMotorCount - m.powerLimit;    	// k2·|ω| + K3/电机数 - P_limit

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

    /**
     * @brief 钳制后的总预测功率（±cfg_.powerMax，供显示/诊断）
     *
     * 参考 COD 哨兵：制动回馈使 τ·ω 为负、总功率预测可负，
     * 显示用钳制值防误导；限幅判断仍用 GetTotalPower() 原始值。
     * @return 总功率 (W)，限制在 ±powerMax 内（powerMax<=0 不钳制）
     */
    float GetTotalPowerClamped() const
    {
        float sum = GetTotalPower();
        if (cfg_.powerMax > 0.0f) {
            if (sum >  cfg_.powerMax) sum =  cfg_.powerMax;
            if (sum < -cfg_.powerMax) sum = -cfg_.powerMax;
        }
        return sum;
    }

    float GetK1() const { return k1_; }
    float GetK2() const { return k2_; }
    float GetKTauOmega() const { return kTauOmega_; }

    /** @brief 本次周期 RLS 回归量 Στ²（Predict 后有效，供诊断打印） */
    float GetSumTorque2()  const { return sumTorque2_; }
    /** @brief 本次周期 RLS 回归量 Σ|ω|（原始值，未归一化） */
    float GetSumAbsOmega() const { return sumAbsOmega_; }
    /** @brief 本次周期 Στ·ω（RLS target 扣减项） */
    float GetSumTauOmega() const { return sumTauOmega_; }

    /**
     * @brief 本次周期 RLS 目标值（target = P_meas − K3 − Στ·ω，Predict 后有效，供诊断打印）
     * @return RLS 输入目标功率 (W)
     */
    float GetTarget() const
    {
        float t = measuredPower_ - cfg_.k3;
        if (cfg_.tauOmegaEnable) {
            t -= kTauOmega_ * sumTauOmega_;
        }
        return t;
    }

private:
    Config cfg_;

    float k1_ = 0.0f;                           // τ² 铜损系数
    float k2_ = 0.0f;                           // 转速线性损耗系数（|ω| 项）
    float kTauOmega_ = 1.0f;                    // τ·ω 耦合系数（固定，不参与 RLS）

    MotorState motors_[kMotorCount]{};          // 各电机状态（转矩/转速/功率等）
    float      membership_[kMotorCount]{};      // 各电机隶属度（功率分配系数 0~1）

    alg::rls::RLS<2> rls_;                      // RLS 在线辨识器（双参数：辨 k1、k2）
    alg::rls::RLS<1> rls1_;                     // RLS 单参数辨识器（fixK2 模式：只辨 k1(τ²)）
    bool rlsInited_ = false;

    float measuredPower_ = 0.0f;                // 功率计实测总功率
    bool  powerValid_    = false;               // 功率计是否已收到有效帧（首帧保护）

    float sumTorque2_  = 0.0f;                  // 本次周期 Στ²（RLS 输入，供诊断打印）
    float sumAbsOmega_ = 0.0f;                  // 本次周期 Σ|ω|（RLS 输入，供诊断打印）
    float sumTauOmega_ = 0.0f;                  // 本次周期 Στ·ω（RLS target 扣减项）

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
