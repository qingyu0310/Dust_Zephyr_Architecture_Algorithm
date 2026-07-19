/**
 * @file stability.hpp
 * @author qingyu
 * @brief 辨识用稳定判据与扰动工具 — 内联实现，仅辨识阶段使用
 * @version 0.1
 * @date 2026-07-19
 */

#pragma once

#pragma message "Compiling Algorithm/Identify/Stability"

#include <cmath>
#include <cstdint>

namespace stability {

/**
 * @brief 滑动窗口稳定判据
 *
 * 积累 kWinSize 帧数据后，判断极差和斜率是否在限内。
 *
 * @tparam kWinSize  滑动窗口大小（帧数）
 */
template <uint32_t kWinSize>
class WinStable final
{
public:
    constexpr WinStable() = default;

    /**
     * @brief 推入一帧温度并判断是否稳定
     *
     * @param temp_c  当前温度 (°C)
     * @param dt_s    帧间隔 (s)，用于斜率计算
     * @param slope_limit  允许的最大斜率 (°C/s)
     * @param noise_limit  允许的最大极差 (°C)
     * @return true  窗口满且稳定
     */
    bool Check(float temp_c, float dt_s, float slope_limit = 0.01f, float noise_limit = 0.10f)
    {
        if (cnt_ < kWinSize) {
            Push(temp_c);
            return false;
        }

        float mn = buf_[0], mx = buf_[0];
        for (uint32_t i = 1; i < kWinSize; i++) {
            if (buf_[i] < mn) mn = buf_[i];
            if (buf_[i] > mx) mx = buf_[i];
        }

        const float slope = dt_s > 0.0f ? (temp_c - prev_) / dt_s : 0.0f;
        Push(temp_c);

        return std::abs(slope) <= slope_limit && (mx - mn) <= noise_limit;
    }

    void Reset()
    {
        head_ = 0;
        cnt_  = 0;
        prev_ = 0.0f;
    }

private:
    float buf_[kWinSize] {};
    float prev_     = 0.0f;
    uint32_t head_  = 0;
    uint32_t cnt_   = 0;

    void Push(float val)
    {
        buf_[head_] = val;
        head_ = (head_ + 1) % kWinSize;
        if (cnt_ < kWinSize) cnt_++;
        prev_ = val;
    }
};

/**
 * @brief 方波发生器：±amp
 */
class SquareGen final
{
public:
    float WaveGen(float amp, uint32_t period)
    {
        return (++cnt_ / period) % 2 ? amp : -amp;
    }
    void Reset() { cnt_ = 0; }
private:
    uint32_t cnt_ = 0;
};

/**
 * @brief 三角波发生器：-amp ~ +amp
 */
class TriangleGen final
{
public:
    float WaveGen(float amp, uint32_t period)
    {
        uint32_t t = cnt_++ % period;
        float v = 2.0f * static_cast<float>(t) / static_cast<float>(period);
        return (v < 1.0f) ? (v * amp) : ((2.0f - v) * amp);
    }
    void Reset() { cnt_ = 0; }
private:
    uint32_t cnt_ = 0;
};

/**
 * @brief 正弦波发生器：amp·sin(2π·freq·t)
 */
class SineGen final
{
public:
    float WaveGen(float amp, float freq_hz, float dt_s)
    {
        phase_ += freq_hz * dt_s;
        if (phase_ >= 1.0f) phase_ -= 1.0f;
        return amp * std::sin(2.0f * 3.14159265f * phase_);
    }
    void Reset() { phase_ = 0.0f; }
private:
    float phase_ = 0.0f;
};

} // namespace stability // namespace stability
