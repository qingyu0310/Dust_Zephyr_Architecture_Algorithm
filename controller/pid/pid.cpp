/**
 * @file pid.cpp
 * @author qingyu
 * @brief 位置式 PID 控制器实现：条件积分抗饱和、角度微分和实际周期
 * @version 0.2
 * @date 2026-10-31
 *
 * @copyright Copyright (c) 2026
 */
#pragma message "Compiling Algorithm/Controller/Pid"

#include "pid.hpp"
#include <cmath>
#include <algorithm>

namespace alg::pid {

namespace {

constexpr float kPi  = 3.141592653589793f;
constexpr float k2Pi = 6.283185307179586f;

/**
 * @brief 将任意圈数的角度差归一化到劣弧
 * @param value 角度差，单位 rad
 * @return [-pi, pi] 范围内的角度差
 */
float WrapAngle(float value)
{
    value = std::fmod(value, k2Pi);
    if (value > kPi) {
        value -= k2Pi;
    } else if (value < -kPi) {
        value += k2Pi;
    }
    return value;
}

}

/**
 * @brief 构造并初始化控制器
 * @param cfg PID 配置参数
 */
Pid::Pid(const Config& cfg)
{
    Init(cfg);
}

/**
 * @brief 重配参数并清除控制历史；在线调参可使用各增益 setter 保留历史
 * @param cfg PID 配置参数
 */
void Pid::Init(const Config& cfg)
{
    cfg_ = cfg;
    cfg_.iOutMax = std::fabs(cfg_.iOutMax);
    cfg_.outMax = std::fabs(cfg_.outMax);
    cfg_.deadZone = std::fabs(cfg_.deadZone);
    cfg_.iSpeedThreshLo = std::fabs(cfg_.iSpeedThreshLo);
    cfg_.iSpeedThreshHi = std::fabs(cfg_.iSpeedThreshHi);
    cfg_.iSeparateThresh = std::fabs(cfg_.iSeparateThresh);
    SetDt(cfg_.dt);
    Reset();
}

/**
 * @brief 清除控制历史，保留调用方设置的目标值和测量值
 */
void Pid::Reset()
{
    out_ = 0.0f;
    preOut_ = 0.0f;
    integralError_ = 0.0f;
    dLpfOutput_ = 0.0f;
    preNow_ = now_;
    preTarget_ = target_;
    preError_ = 0.0f;
    historyReady_ = false;
}

/**
 * @brief 更新周期并重算微分低通系数；无效周期由计算入口拒绝
 * @param dt 本次计算与上次计算的实际间隔，单位 s
 */
void Pid::SetDt(float dt)
{
    cfg_.dt = dt;
    dLpfAlpha_ = 0.0f;
    if (std::isfinite(dt) && dt > 0.0f && std::isfinite(cfg_.dLpfCutoffHz) && cfg_.dLpfCutoffHz > 0.0f) {
        const float tau = 1.0f / (2.0f * kPi * cfg_.dLpfCutoffHz);
        dLpfAlpha_ = tau / (tau + dt);
    }
}

/**
 * @brief 位置式 PID 计算
 * @return 控制器输出；输入或参数无效时清历史并返回零
 */
float Pid::Calc()
{
    return CalcImpl(target_ - now_, false);
}

/**
 * @brief 角度 PID，误差及相邻采样差均采用劣弧
 * @return 控制器输出
 */
float Pid::CalcAngle()
{
    return CalcImpl(WrapAngle(target_ - now_), true);
}

/**
 * @brief 公共 PID 计算，积分限幅在累加后生效，输出饱和不清空已有积分
 * @param error 本次误差
 * @param angle 是否按角度处理微分和目标增量
 * @return 本次控制输出
 */
float Pid::CalcImpl(float error, bool angle)
{
    if (!std::isfinite(error) || !std::isfinite(target_) || !std::isfinite(now_) || !std::isfinite(cfg_.dt) || cfg_.dt <= 0.0f || !std::isfinite(cfg_.kp) || !std::isfinite(cfg_.ki) || !std::isfinite(cfg_.kd) || !std::isfinite(cfg_.kf) || !std::isfinite(cfg_.iOutMax) || !std::isfinite(cfg_.outMax) || !std::isfinite(cfg_.deadZone) || !std::isfinite(cfg_.iSpeedThreshLo) || !std::isfinite(cfg_.iSpeedThreshHi) || !std::isfinite(cfg_.iSeparateThresh) || !std::isfinite(cfg_.dLpfCutoffHz) || !std::isfinite(integralError_)) {
        Reset();
        return 0.0f;
    }

    const bool variableIntegral = cfg_.iSpeedThreshLo != 0.0f || cfg_.iSpeedThreshHi != 0.0f;
    if (variableIntegral && cfg_.iSpeedThreshHi <= cfg_.iSpeedThreshLo) {
        Reset();
        return 0.0f;
    }

    // 死区只处理本次误差，不覆盖调用方的目标。
    float absError = std::fabs(error);
    if (absError <= cfg_.deadZone) {
        error = 0.0f;
    } else {
        error += (error > 0.0f) ? -cfg_.deadZone : cfg_.deadZone;
    }
    absError = std::fabs(error);

    // 保留原有变速积分上限清零和积分分离语义。
    float speedRatio = 1.0f;
    if (variableIntegral)
    {
        if (absError <= cfg_.iSpeedThreshLo) {
            speedRatio = 1.0f;
        } else if (absError >= cfg_.iSpeedThreshHi) {
            speedRatio = 0.0f;
            integralError_ = 0.0f;
        } else {
            speedRatio = (cfg_.iSpeedThreshHi - absError) / (cfg_.iSpeedThreshHi - cfg_.iSpeedThreshLo);
        }
    }
    const bool integralEnabled = cfg_.ki != 0.0f && (cfg_.iSeparateThresh == 0.0f || absError < cfg_.iSeparateThresh);
    if (!integralEnabled) {
        integralError_ = 0.0f;
    }

    // 首次计算及切换角度模式时不把默认历史值当成真实采样。
    const bool hasHistory = historyReady_ && angleMode_ == angle;
    float dRaw = 0.0f;
    if (hasHistory && cfg_.kd != 0.0f)
    {
        float delta = (cfg_.dFirst == DFirst::Enable) ? now_ - preNow_ : error - preError_;
        if (angle) {
            delta = WrapAngle(delta);
        }
        dRaw = cfg_.kd * delta / cfg_.dt;
        if (cfg_.dFirst == DFirst::Enable) {
            dRaw = -dRaw;
        }
    }
    if (!hasHistory) {
        dLpfOutput_ = 0.0f;
    }
    const float dOut = (dLpfAlpha_ > 0.0f) ? dRaw + dLpfAlpha_ * (dLpfOutput_ - dRaw) : dRaw;

    float targetDelta = target_ - preTarget_;
    if (angle) {
        targetDelta = WrapAngle(targetDelta);
    }
    const float pOut = cfg_.kp * error;
    const float fOut = cfg_.kf * targetDelta;

    float previousIntegral = integralError_;
    float candidateIntegral = previousIntegral;

    // 上下界始终为正负同一幅值，支持负 Ki；手工设置的积分也必须受限。
    float iClamp = 0.0f;
    if (cfg_.iOutMax > 0.0f && cfg_.ki != 0.0f) {
        iClamp = cfg_.iOutMax / std::fabs(cfg_.ki);
        previousIntegral = std::clamp(previousIntegral, -iClamp, iClamp);
        candidateIntegral = previousIntegral;
    }
    if (integralEnabled) {
        candidateIntegral += speedRatio * cfg_.dt * error;
    }
    if (cfg_.iOutMax > 0.0f && cfg_.ki != 0.0f) {
        candidateIntegral = std::clamp(candidateIntegral, -iClamp, iClamp);
    }

    // 条件积分：只拒绝继续推向饱和的增量，保留偏置并允许反向退饱和。
    const float candidateIOut = cfg_.ki * candidateIntegral;
    const float candidateOut = pOut + candidateIOut + dOut + fOut;
    const float integralDeltaOut = cfg_.ki * (candidateIntegral - previousIntegral);
    if (cfg_.outMax > 0.0f && ((candidateOut > cfg_.outMax && integralDeltaOut > 0.0f) || (candidateOut < -cfg_.outMax && integralDeltaOut < 0.0f))) {
        candidateIntegral = previousIntegral;
    }

    const float iOut = cfg_.ki * candidateIntegral;
    const float output = pOut + iOut + dOut + fOut;
    if (!std::isfinite(output) || !std::isfinite(candidateIntegral)) {
        Reset();
        return 0.0f;
    }

    integralError_ = candidateIntegral;
    out_ = (cfg_.outMax > 0.0f) ? std::clamp(output, -cfg_.outMax, cfg_.outMax) : output;
    dLpfOutput_ = dOut;
    preNow_ = now_;
    preTarget_ = target_;
    preOut_ = out_;
    preError_ = error;
    angleMode_ = angle;
    historyReady_ = true;
    return out_;
}

} // namespace alg::pid
