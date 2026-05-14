/**
 * @file power_ctrl.hpp
 * @author qingyu
 * @brief Motor power control — model prediction, RLS identification,
 *        membership allocation, and torque limiting
 * @version 0.1
 * @date 2026-04-30
 */

#pragma once

#include <stdint.h>
#include "rls.hpp"

namespace alg::power_ctrl {

inline constexpr uint8_t kMaxMotors = 4;

/**
 * @brief  Per-motor state
 */
struct MotorState {
    float torque    = 0.0f;         ///< τ = current × torqueK
    float omega     = 0.0f;         ///< ω = rpm / 9.55 (rad/s)
    float torque2   = 0.0f;         ///< τ²
    float omega2    = 0.0f;         ///< ω²
    float powerPred = 0.0f;         ///< P_in = K1·τ² + K2·ω² + τ·ω + K3
    float pidErr       = 0.0f;      ///< PID error (for membership)
    float powerLimit    = 0.0f;     ///< allocated power budget for this motor
    float targetCurrent = 0.0f;     ///< PID target current (cascade input)
    float currentOut    = 0.0f;     ///< limited output current
};

/**
 * @brief  Power controller for one motor group (e.g. 4 steer motors)
 *
 * @par Pipeline (call once per control cycle):
 *       @n 1. SetMotorData(i, current, velocity, pidErr)  — for each motor
 *       @n 2. Predict()                                    — run model + RLS
 *       @n 3. Allocate(totalBudget)                        — membership + limit
 *       @n 4. GetLimitedCurrent(i)                         — read result
 *
 * @par RLS:
 *       Internally identifies K1 (copper loss) and K2 (iron loss).
 *       K3 is a fixed constant.  The RLS runs when enabled.
 */
class PowerCtrl final
{
public:
    struct Config {
        uint8_t motorCount = 4;         ///< motors in this group (≤ kMaxMotors)
        float   torqueK    = 4.577e-5f; ///< current → torque (M3508 constant)
        float   k1Init     = 1.453e-7f; ///< RLS initial K1
        float   k2Init     = 1.453e-7f; ///< RLS initial K2
        float   k3         = 3.0f;      ///< fixed constant loss
        float   errUpper   = 50.0f;     ///< membership: full-need threshold
        float   errLower   = 0.01f;     ///< membership: no-need threshold
        float   rlsLambda  = 0.99999f;  ///< RLS forgetting factor
        bool    rlsEnable  = false;     ///< false = use fixed K1/K2
        bool    tauOmegaEnable = true;  ///< include τ·ω term in power model
    };

    PowerCtrl() : PowerCtrl(Config{}) {}
    explicit PowerCtrl(const Config& cfg);
    void Init(const Config& cfg);

    /* ---------- 输入 ---------- */

    /**
     * @brief  Feed raw data for one motor
     * @param idx       motor index [0, motorCount)
     * @param torque    measured torque (N·m)
     * @param omega     angular velocity (rad/s)
     * @param pidErr    PID position/speed error
     */
    void SetMotorData(uint8_t idx, float torque, float omega, float pidErr);

    /** @brief  设置 PID 目标电流（串级输入） */
    void SetTarget(uint8_t idx, float current);

    /** @brief  输入实测总功率（来自功率计），供 RLS 使用 */
    void SetMeasuredPower(float power) { measuredPower_ = power; }

    /* ---------- 处理 ---------- */

    /** @brief  Step ① + ②: power prediction + RLS update */
    void Predict();

    /**
     * @brief  Step ③ + ④: membership allocation + torque limit
     * @param totalBudget  total power budget for this group (W)
     */
    void Allocate(float totalBudget);

    /* ---------- 结果 ---------- */

    /** @brief  Limited output current for motor idx */
    float GetLimitedCurrent(uint8_t idx) const;

    /** @brief  Sum of predicted power across all motors */
    float GetTotalPower() const;

    /* ---------- RLS 参数 ---------- */

    float GetK1() const { return k1_; }
    float GetK2() const { return k2_; }

private:
    Config cfg_;

    float k1_;     ///< copper-loss coefficient
    float k2_;     ///< iron-loss coefficient

    MotorState motors_[kMaxMotors]{};
    float      membership_[kMaxMotors]{};

    alg::rls::RLS<2, 1> rls_;   ///< 2-input, 1-output online identifier
    bool rlsInited_ = false;

    float measuredPower_ = 0.0f;  ///< measured total power from power meter

    /* 二次方程求根: A·τ² + B·τ + C = 0 */
    static float SolveTorque(float a, float b, float c, bool positive);
};

} // namespace alg::power_ctrl
