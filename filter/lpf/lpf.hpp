/**
 * @file lpf.hpp
 * @author qingyu
 * @brief First-order low-pass filter (cutoff-frequency based)
 * @version 0.1
 * @date 2026-04-30
 *
 * @par Usage
 * @code
 *   alg::filter::LowPassFilter lpf(100.0f, 0.001f);  // 100 Hz, 1 ms
 *   float out = lpf.Update(raw);
 * @endcode
 */

#pragma once

namespace alg::filter {

/**
 * @brief First-order RC low-pass filter
 *
 *        out[n] = alpha * in[n] + (1-alpha) * out[n-1]
 *
 *        alpha = dt / (dt + RC),  RC = 1 / (2π · fc)
 */
class LowPassFilter final
{
public:
    LowPassFilter() = default;

    /**
     * @brief Construct and initialize
     * @param cutoffHz  cutoff frequency in Hz (≤ 0 = passthrough)
     * @param dt        sampling period in seconds
     */
    LowPassFilter(float cutoffHz, float dt) { Init(cutoffHz, dt); }

    /**
     * @brief (Re)configure filter
     * @param cutoffHz  cutoff frequency in Hz (≤ 0 = passthrough)
     * @param dt        sampling period in seconds
     */
    void Init(float cutoffHz, float dt)
    {
        if (cutoffHz <= 0.0f) {
            alpha_ = 1.0f;      // passthrough
        } else {
            float rc = 1.0f / (k2Pi * cutoffHz);
            alpha_ = dt / (dt + rc);
        }
        output_ = 0.0f;
        initialized_ = true;
    }

    /**
     * @brief Update filter with new sample
     * @param input  raw input value
     * @return       filtered output
     */
    float Update(float input)
    {
        if (!initialized_) return input;
        output_ = alpha_ * input + (1.0f - alpha_) * output_;
        return output_;
    }

    float GetOutput() const { return output_; }
    void  Reset(float val = 0.0f) { output_ = val; }

private:
    static constexpr float k2Pi = 6.283185307179586f;

    float alpha_       = 0.0f;
    float output_      = 0.0f;
    bool  initialized_ = false;
};

} // namespace alg::filter
