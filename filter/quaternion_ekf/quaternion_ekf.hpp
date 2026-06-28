/**
 * @file quaternion_ekf.hpp
 * @author qingyu
 * @brief 基于四元数的 IMU 姿态扩展 Kalman 滤波器
 * @version 0.2
 * @date 2026-06-28
 */

#pragma once
#pragma message "Compiling Algorithm/Filter/QuaternionEkf"

#include "kalman_ekf.hpp"
#include <stdint.h>

struct Sample;

namespace alg::attitude {

/**
 * @brief 基于四元数的 IMU 姿态扩展 Kalman 滤波器。
 *
 * 状态向量维度为 6：
 * - q0, q1, q2, q3: 姿态四元数。
 * - gyro_bias_x, gyro_bias_y: 陀螺 x/y 零偏。
 *
 * 观测量维度为 3，使用归一化后的加速度方向作为重力方向观测。
 */
class QuaternionEkf final
{
public:
    struct Config {
        float process_noise_quaternion = 10.0f;
        float process_noise_bias       = 0.001f;
        float measurement_noise        = 1000000.0f;
        float fading_factor            = 1.0f;
        float accel_lpf_coefficient    = 0.0f;
        float chi_square_threshold     = 1e-8f;
        float stable_gyro_threshold    = 0.3f;
        float stable_accel_reference   = 9.8f;
        float stable_accel_tolerance   = 0.5f;
        float bias_correction_limit    = 1e-2f;
        float bias_covariance_limit    = 10000.0f;
        uint32_t divergence_limit      = 50;
    };

    struct State {
        bool Initialized  = false;
        bool ConvergeFlag = false;
        bool StableFlag   = false;

        uint64_t ErrorCount  = 0;
        uint64_t UpdateCount = 0;

        float q[4]                 {1.0f, 0.0f, 0.0f, 0.0f};
        float GyroBias[3]          {0.0f, 0.0f, 0.0f};
        float Gyro[3]              {0.0f, 0.0f, 0.0f};
        float Accel[3]             {0.0f, 0.0f, 0.0f};
        float OrientationCosine[3] {0.0f, 0.0f, 0.0f};

        float accLPFcoef        = 0.0f;
        float gyro_norm         = 0.0f;
        float accl_norm         = 0.0f;
        float AdaptiveGainScale = 1.0f;

        float Roll = 0.0f;
        float Pitch = 0.0f;
        float Yaw = 0.0f;
        float YawTotalAngle = 0.0f;

        float Q1 = 10.0f;
        float Q2 = 0.001f;
        float R  = 1000000.0f;
        float dt = 0.0f;
        float ChiSquare = 0.0f;
        float ChiSquareTestThreshold = 1e-8f;
        float lambda = 1.0f;

        int16_t YawRoundCount = 0;
        float   YawAngleLast  = 0.0f;
    };

    QuaternionEkf() = default;

    void Init(const Config& config);
    void Update(const Sample& sample);

    const State& GetState() const { return state_; }
    State&       GetState()       { return state_; }
    const filter::KalmanFilter& Filter() const { return filter_; }
    filter::KalmanFilter&       Filter()       { return filter_; }

private:
    static constexpr uint8_t kStateSize       = 6;
    static constexpr uint8_t kMeasurementSize = 3;

    static void ObserveHook(filter::KalmanFilter *kf);
    static void LinearizeHook(filter::KalmanFilter *kf);
    static void SetHHook(filter::KalmanFilter *kf);
    static void UpdatePosteriorHook(filter::KalmanFilter *kf);

    static float FastInverseSqrt(float x);

    void Observe();
    void LinearizeAndFade();
    void SetObservationJacobian();
    void UpdatePosterior();

    Config config_;
    State  state_ {};
    filter::KalmanFilter filter_ {};
    float  chi_square_data_[1] = {};  // 卡方检验结果 1×1
};

QuaternionEkf& GlobalQuaternionEkf();

} // namespace alg::attitude

using QEKF_INS_t = alg::attitude::QuaternionEkf::State;

extern QEKF_INS_t QEKF_INS;
extern float      chiSquare;
extern float      ChiSquareTestThreshold;

void IMU_QuaternionEKF_Init(float process_noise1, float process_noise2,
                            float measure_noise, float lambda, float lpf);
void IMU_QuaternionEKF_Update(float gx, float gy, float gz,
                              float ax, float ay, float az, float dt);
