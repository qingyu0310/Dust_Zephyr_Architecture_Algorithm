/**
 * @file pid.hpp
 * @author qingyu
 * @brief Positional PID controller with D-first, integral separation,
 *        variable integral speed, and D-term LPF
 * @version 0.1
 * @date 2026-04-29
 *
 * @copyright Copyright (c) 2026
 */

#pragma once

#include <stdint.h>

namespace alg::pid {

/// Derivative-on-measurement mode (avoids derivative kick on target change)
enum class DFirst : uint8_t {
    Disable = 0,  ///< D on error
    Enable,       ///< D on measurement (recommended for setpoint tracking)
};

/**
 * @brief Positional PID controller
 *
 * @par Features
 *      - Standard P/I/D + feedforward
 *      - Dead zone
 *      - Variable integral speed (变速积分)
 *      - Integral separation (积分分离)
 *      - Derivative on measurement (微分先行)
 *      - D-term low-pass filter
 *      - Output / integral clamping
 *      - Angle-mode with wraparound handling
 *
 * @par Usage
 * @code
 *   alg::pid::Pid pid({.kp = 3.0f, .ki = 0.1f, .kd = 0.05f, .outMax = 60.0f});
 *   float out = pid.Calc(target, now);
 * @endcode
 */
class Pid final
{
public:
    /**
     * @brief PID configuration (all fields have defaults)
     *
     * @note Use C++20 designated initializers for clarity:
     *       @c Pid::Config{.kp=3, .ki=0.1, .outMax=60}
     */
    struct Config {
        float kp               = 0.0f;
        float ki               = 0.0f;
        float kd               = 0.0f;
        float kf               = 0.0f;   ///< feedforward gain
        float iOutMax          = 0.0f;   ///< 0 = no integral limit
        float outMax           = 0.0f;   ///< 0 = no output limit
        float dt               = 0.001f; ///< control period (seconds)
        float deadZone         = 0.0f;
        float iVariableSpeedA  = 0.0f;   ///< variable-integral lower threshold
        float iVariableSpeedB  = 0.0f;   ///< variable-integral upper threshold
        float iSeparateThresh  = 0.0f;   ///< integral-separation threshold (0 = off)
        DFirst dFirst          = DFirst::Disable;
        float dLpfTau          = 0.0f;   ///< D LPF time constant (0 = no filter)
    };

    Pid() = default;
    explicit Pid(const Config& cfg) : cfg_(cfg) {}

    /** @brief (Re)configure all parameters */
    void Init(const Config& cfg) { cfg_ = cfg; }

    /**
     * @brief Shadow-set params (double-buffered, safe during runtime)
     * @details Writes to a shadow config. Calc() swaps it in atomically
     *          at the next cycle, so params never change mid-computation.
     */
    void SetShadow(const Config& cfg)
    {
        shadowCfg_     = cfg;
        shadowPending_ = true;
    }

    /**
     * @brief One-shot PID calculation
     * @param target  setpoint
     * @param now     feedback
     * @return        controller output
     */
    float Calc(float target, float now)
    {
        SetTarget(target);
        SetNow(now);
        return Calc();
    }

    /**
     * @brief Stateful PID calculation (use SetTarget / SetNow first)
     * @return controller output
     */
    float Calc();

    /**
     * @brief Angle-mode PID with wraparound handling (-π ~ π)
     * @return controller output
     */
    float CalcAngle();

    void SetKp(float v)             { cfg_.kp = v; }
    void SetKi(float v)             { cfg_.ki = v; }
    void SetKd(float v)             { cfg_.kd = v; }
    void SetKf(float v)             { cfg_.kf = v; }
    void SetIOutMax(float v)        { cfg_.iOutMax = v; }
    void SetOutMax(float v)         { cfg_.outMax = v; }
    void SetTarget(float v)         { target_ = v; }
    void SetNow(float v)            { now_ = v; }
    void SetIntegralError(float v)  { integralError_ = v; }

    const Config& GetConfig() const { return cfg_; }

    float GetKp()             const { return cfg_.kp; }
    float GetKi()             const { return cfg_.ki; }
    float GetKd()             const { return cfg_.kd; }
    float GetKf()             const { return cfg_.kf; }
    float GetOutMax()         const { return cfg_.outMax; }
    float GetIOutMax()        const { return cfg_.iOutMax; }
    float GetDt()             const { return cfg_.dt; }

    float GetOut()            const { return out_; }
    float GetTarget()         const { return target_; }
    float GetIntegralError()  const { return integralError_; }
    float GetError()          const { return preError_; }

private:

    Config cfg_ {};
    Config shadowCfg_ {};
    volatile bool shadowPending_ = false;

    float target_        = 0.0f;
    float now_           = 0.0f;
    float out_           = 0.0f;

    float preNow_        = 0.0f;
    float preTarget_     = 0.0f;
    float preOut_        = 0.0f;
    float preError_      = 0.0f;

    float integralError_ = 0.0f;
    float dLpfOutput_    = 0.0f;

    static void Clamp(float* v, float lo, float hi);

    /** @brief 公共 PID 计算（双缓冲 → 死区 → 变速积分 → PID → 前馈 → 限幅） */
    float CalcImpl(float error);
};

} // namespace alg::pid
