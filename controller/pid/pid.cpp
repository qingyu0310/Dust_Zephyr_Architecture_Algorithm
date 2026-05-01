/**
 * @file pid.cpp
 * @author qingyu
 * @brief Positional PID implementation
 * @version 0.1
 * @date 2026-04-29
 */

#include "pid.hpp"

namespace alg::pid {

namespace {
    constexpr float kPi  = 3.141592653589793f;
    constexpr float k2Pi = 6.283185307179586f;
}

/*  Clamp helper                                                       */

void Pid::Clamp(float* v, float lo, float hi)
{
    if (*v < lo) *v = lo;
    if (*v > hi) *v = hi;
}

/*  Stateful PID                                                       */

float Pid::Calc()
{
    /* ---- 0. 双缓冲切换 ---- */
    if (shadowPending_) {
        cfg_ = shadowCfg_;
        shadowPending_ = false;
    }

    /* ---------- 1. 误差 ---------- */
    float error    = target_ - now_;
    float absError = Abs(error);

    /* ---------- 2. 死区 ---------- */
    if (absError < cfg_.deadZone) {
        target_   = now_;
        error     = 0.0f;
        absError  = 0.0f;
    } else if (error > 0.0f) {
        error    -= cfg_.deadZone;
        absError  = Abs(error);
    } else if (error < 0.0f) {
        error    += cfg_.deadZone;
        absError  = Abs(error);
    }

    /* ---------- 3. 变速积分系数 ---------- */
    float speedRatio = 1.0f;
    if (cfg_.iVariableSpeedA != 0.0f || cfg_.iVariableSpeedB != 0.0f) {
        if (absError <= cfg_.iVariableSpeedA) {
            speedRatio = 1.0f;
        } else if (absError >= cfg_.iVariableSpeedB) {
            speedRatio = 0.0f;
            integralError_ = 0.0f;          // 大误差时清零防饱和
        } else {
            float denom = cfg_.iVariableSpeedB - cfg_.iVariableSpeedA;
            if (denom > 0.0f)
                speedRatio = (cfg_.iVariableSpeedB - absError) / denom;
        }
    }

    /* ---------- 4. P项 ---------- */
    const float pOut = cfg_.kp * error;

    /* ---------- 5. I项 ---------- */
    float iOut = 0.0f;

    // 积分限幅
    if (cfg_.iOutMax != 0.0f && cfg_.ki != 0.0f) {
        float iClamp = cfg_.iOutMax / cfg_.ki;
        Clamp(&integralError_, -iClamp, iClamp);
    }

    if (cfg_.iSeparateThresh == 0.0f || absError < cfg_.iSeparateThresh) {
        if (isnanf(error)) error = preError_;
        integralError_ += speedRatio * cfg_.dt * error;
        iOut = cfg_.ki * integralError_;
    } else {
        integralError_ = 0.0f;
    }

    /* ---------- 6. D项 ---------- */
    float dRaw = 0.0f;
    if (cfg_.dFirst == DFirst::Enable) {
        dRaw = -cfg_.kd * (now_ - preNow_) / cfg_.dt;
    } else {
        dRaw = cfg_.kd * (error - preError_) / cfg_.dt;
    }

    float dOut = dRaw;
    if (cfg_.dLpfTau > 0.0f) {
        float alpha = cfg_.dLpfTau / (cfg_.dLpfTau + cfg_.dt);
        dOut = alpha * dLpfOutput_ + (1.0f - alpha) * dRaw;
        dLpfOutput_ = dOut;
    }

    /* ---------- 7. 前馈 ---------- */
    const float fOut = cfg_.kf * (target_ - preTarget_);

    /* ---------- 8. 输出 ---------- */
    out_ = pOut + iOut + dOut + fOut;
    if (cfg_.outMax != 0.0f) {
        Clamp(&out_, -cfg_.outMax, cfg_.outMax);
    }

    /* ---------- 9. 状态保持 ---------- */
    preNow_    = now_;
    preTarget_ = target_;
    preOut_    = out_;
    preError_  = error;

    return out_;
}

/*  Angle PID (劣弧 + 防卡角)                                            */

float Pid::CalcAngle()
{
    /* ---- 0. 双缓冲切换 ---- */
    if (shadowPending_) {
        cfg_ = shadowCfg_;
        shadowPending_ = false;
    }

    /* ---------- 1. 劣弧角度差 ---------- */
    float error = target_ - now_;
    if (error > kPi)
        error -= k2Pi;
    else if (error < -kPi)
        error += k2Pi;

    float absError = Abs(error);

    /* ---------- 2. 死区 ---------- */
    if (absError < cfg_.deadZone) {
        target_   = now_;
        error     = 0.0f;
        absError  = 0.0f;
    } else if (error > 0.0f) {
        error    -= cfg_.deadZone;
        absError  = Abs(error);
    } else if (error < 0.0f) {
        error    += cfg_.deadZone;
        absError  = Abs(error);
    }

    /* ---------- 3. 变速积分 ---------- */
    float speedRatio = 1.0f;
    if (cfg_.iVariableSpeedA != 0.0f || cfg_.iVariableSpeedB != 0.0f) {
        if (absError <= cfg_.iVariableSpeedA) {
            speedRatio = 1.0f;
        } else if (absError >= cfg_.iVariableSpeedB) {
            speedRatio = 0.0f;
        } else {
            float denom = cfg_.iVariableSpeedB - cfg_.iVariableSpeedA;
            if (denom > 0.0f)
                speedRatio = (cfg_.iVariableSpeedB - absError) / denom;
        }
    }

    /* ---------- 4. P项 ---------- */
    const float pOut = cfg_.kp * error;

    /* ---------- 5. I项 ---------- */
    float iOut = 0.0f;

    if (cfg_.iOutMax != 0.0f) {
        float iClamp = cfg_.iOutMax / cfg_.ki;
        Clamp(&integralError_, -iClamp, iClamp);
    }

    // 防卡角: 输出饱和或误差反向时清零
    float absOut = Abs(out_);
    if ((cfg_.outMax != 0.0f && absOut >= cfg_.outMax) ||
        (preError_ > 0.0f && error < 0.0f) ||
        (preError_ < 0.0f && error > 0.0f))
    {
        integralError_ = 0.0f;
    }

    if (cfg_.iSeparateThresh == 0.0f || absError < cfg_.iSeparateThresh) {
        integralError_ += speedRatio * cfg_.dt * error;
        iOut = cfg_.ki * integralError_;
    } else {
        integralError_ = 0.0f;
    }

    /* ---------- 6. D项 ---------- */
    float dOut;
    if (cfg_.dFirst == DFirst::Enable) {
        dOut = -cfg_.kd * (now_ - preNow_) / cfg_.dt;
    } else {
        dOut = cfg_.kd * (error - preError_) / cfg_.dt;
    }

    /* ---------- 7. 前馈 ---------- */
    const float fOut = cfg_.kf * (target_ - preTarget_);

    /* ---------- 8. 输出 ---------- */
    out_ = pOut + iOut + dOut + fOut;
    if (cfg_.outMax != 0.0f) {
        Clamp(&out_, -cfg_.outMax, cfg_.outMax);
    }

    /* ---------- 9. 状态保持 ---------- */
    preNow_    = now_;
    preTarget_ = target_;
    preOut_    = out_;
    preError_  = error;

    return out_;
}

} // namespace alg::pid
