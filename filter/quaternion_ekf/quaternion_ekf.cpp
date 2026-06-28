/**
 * @file quaternion_ekf.cpp
 * @author qingyu
 * @brief 基于四元数的 IMU 姿态扩展 Kalman 滤波器 — Eigen 版本
 * @version 0.2
 * @date 2026-06-28
 */
#pragma message "Compiling Algorithm/Filter/QuaternionEkf"

#include "quaternion_ekf.hpp"
#include "imu.hpp"

#include <zephyr/kernel.h>
#include <math.h>
#include <string.h>
#include <Eigen/Dense>

namespace {

constexpr float kBaseTransition[36] {
     1.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,
     0.0f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f,
     0.0f,  0.0f,  1.0f,  0.0f,  0.0f,  0.0f,
     0.0f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,
     0.0f,  0.0f,  0.0f,  0.0f,  1.0f,  0.0f,
     0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  1.0f,
};

constexpr float kInitialCovariance[36] {
     100000.0f, 0.1f,      0.1f,       0.1f,       0.1f,    0.1f,
     0.1f,      100000.0f, 0.1f,       0.1f,      0.1f,    0.1f,
     0.1f,      0.1f,      100000.0f,  0.1f,       0.1f,    0.1f,
     0.1f,      0.1f,      0.1f,       100000.0f,  0.1f,    0.1f,
     0.1f,      0.1f,      0.1f,       0.1f,       100.0f,  0.1f,
     0.1f,      0.1f,      0.1f,       0.1f,       0.1f,    100.0f,
};

constexpr float kRadToDeg = 57.295779513f;
constexpr float kYawBiasStaticThreshold = 0.02f;

float Clamp(float value, float min_value, float max_value)
{
    if (value < min_value) return min_value;
    if (value > max_value) return max_value;
    return value;
}

} // namespace

namespace alg::attitude {

static void SyncLegacyExports(const QuaternionEkf::State& state)
{
    ::QEKF_INS = state;
    ::chiSquare = state.ChiSquare;
    ::ChiSquareTestThreshold = state.ChiSquareTestThreshold;
}

void QuaternionEkf::Init(const Config& config)
{
    config_ = config;
    if (config_.fading_factor > 1.0f) config_.fading_factor = 1.0f;

    state_ = {};
    state_.Initialized  = true;
    state_.Q1           = config_.process_noise_quaternion;
    state_.Q2           = config_.process_noise_bias;
    state_.R            = config_.measurement_noise;
    state_.lambda       = config_.fading_factor;
    state_.accLPFcoef   = config_.accel_lpf_coefficient;
    state_.ChiSquareTestThreshold = config_.chi_square_threshold;

    (void)filter_.Init(kStateSize, 0, kMeasurementSize);
    filter_.SetUserData(this);
    filter_.User_Func0_f = ObserveHook;
    filter_.User_Func1_f = LinearizeHook;
    filter_.User_Func2_f = SetHHook;
    filter_.User_Func3_f = UpdatePosteriorHook;
    filter_.SkipEq3 = true;
    filter_.SkipEq4 = true;

    filter_.xhat_data[0] = 1.0f;
    filter_.xhat_data[1] = 0.0f;
    filter_.xhat_data[2] = 0.0f;
    filter_.xhat_data[3] = 0.0f;

    memcpy(filter_.F_data, kBaseTransition, sizeof(kBaseTransition));
    memcpy(filter_.P_data, kInitialCovariance, sizeof(kInitialCovariance));

    SyncLegacyExports(state_);
}

void QuaternionEkf::Update(const Sample& sample)
{
    if (!state_.Initialized) Init(Config{});

    state_.dt = sample.dt;

    state_.Gyro[0] = sample.gyro[0] - state_.GyroBias[0];
    state_.Gyro[1] = sample.gyro[1] - state_.GyroBias[1];
    state_.Gyro[2] = sample.gyro[2] - state_.GyroBias[2];

    const float half_gx_dt = 0.5f * state_.Gyro[0] * state_.dt;
    const float half_gy_dt = 0.5f * state_.Gyro[1] * state_.dt;
    const float half_gz_dt = 0.5f * state_.Gyro[2] * state_.dt;

    memcpy(filter_.F_data, kBaseTransition, sizeof(kBaseTransition));
    filter_.F_data[1]  = -half_gx_dt;
    filter_.F_data[2]  = -half_gy_dt;
    filter_.F_data[3]  = -half_gz_dt;
    filter_.F_data[6]  =  half_gx_dt;
    filter_.F_data[8]  =  half_gz_dt;
    filter_.F_data[9]  = -half_gy_dt;
    filter_.F_data[12] =  half_gy_dt;
    filter_.F_data[13] = -half_gz_dt;
    filter_.F_data[15] =  half_gx_dt;
    filter_.F_data[18] =  half_gz_dt;
    filter_.F_data[19] =  half_gy_dt;
    filter_.F_data[20] = -half_gx_dt;

    if (state_.UpdateCount == 0) {
        state_.Accel[0] = sample.accel[0];
        state_.Accel[1] = sample.accel[1];
        state_.Accel[2] = sample.accel[2];
    }

    const float accel_filter_den = state_.dt + state_.accLPFcoef;
    if (accel_filter_den > 0.0f) {
        state_.Accel[0] = state_.Accel[0] * state_.accLPFcoef / accel_filter_den + sample.accel[0] * state_.dt / accel_filter_den;
        state_.Accel[1] = state_.Accel[1] * state_.accLPFcoef / accel_filter_den + sample.accel[1] * state_.dt / accel_filter_den;
        state_.Accel[2] = state_.Accel[2] * state_.accLPFcoef / accel_filter_den + sample.accel[2] * state_.dt / accel_filter_den;
    } else {
        state_.Accel[0] = sample.accel[0];
        state_.Accel[1] = sample.accel[1];
        state_.Accel[2] = sample.accel[2];
    }

    const float accel_inv_norm = FastInverseSqrt(state_.Accel[0] * state_.Accel[0] +
                                                  state_.Accel[1] * state_.Accel[1] +
                                                  state_.Accel[2] * state_.Accel[2]);
    for (uint8_t i = 0; i < kMeasurementSize; i++) {
        filter_.MeasuredVector[i] = state_.Accel[i] * accel_inv_norm;
    }

    const float gyro_inv_norm = FastInverseSqrt(state_.Gyro[0] * state_.Gyro[0] +
                                                 state_.Gyro[1] * state_.Gyro[1] +
                                                 state_.Gyro[2] * state_.Gyro[2]);
    state_.gyro_norm = gyro_inv_norm > 0.0f ? 1.0f / gyro_inv_norm : 0.0f;
    state_.accl_norm = accel_inv_norm > 0.0f ? 1.0f / accel_inv_norm : 0.0f;

    if (state_.gyro_norm < config_.stable_gyro_threshold &&
        state_.accl_norm > config_.stable_accel_reference - config_.stable_accel_tolerance &&
        state_.accl_norm < config_.stable_accel_reference + config_.stable_accel_tolerance) {
        state_.StableFlag = true;
    } else {
        state_.StableFlag = false;
    }

    filter_.Q_data[0]  = state_.Q1 * state_.dt;
    filter_.Q_data[7]  = state_.Q1 * state_.dt;
    filter_.Q_data[14] = state_.Q1 * state_.dt;
    filter_.Q_data[21] = state_.Q1 * state_.dt;
    filter_.Q_data[28] = state_.Q2 * state_.dt;
    filter_.Q_data[35] = state_.Q2 * state_.dt;

    filter_.R_data[0] = state_.R;
    filter_.R_data[4] = state_.R;
    filter_.R_data[8] = state_.R;

    {
        const uint32_t _kf_start = k_cycle_get_32();
        filter_.Update();
        const uint32_t _kf_end = k_cycle_get_32();

        static uint32_t _kf_cnt = 0;
        if (++_kf_cnt % 1000 == 0) {
            const uint64_t _kf_ns = k_cyc_to_ns_floor64(_kf_end - _kf_start);
            printk("[KF] kalman filter: %.3f ms\n", (double)_kf_ns / 1000000.0);
        }
    }

    state_.q[0] = filter_.FilteredValue[0];
    state_.q[1] = filter_.FilteredValue[1];
    state_.q[2] = filter_.FilteredValue[2];
    state_.q[3] = filter_.FilteredValue[3];
    state_.GyroBias[0] = filter_.FilteredValue[4];
    state_.GyroBias[1] = filter_.FilteredValue[5];

    if (state_.StableFlag &&
        fabsf(state_.Gyro[0]) < config_.stable_gyro_threshold &&
        fabsf(state_.Gyro[1]) < config_.stable_gyro_threshold &&
        fabsf(state_.Gyro[2]) < kYawBiasStaticThreshold) {
        const float bias_limit = config_.bias_correction_limit * state_.dt;
        const float bias_error = sample.gyro[2] - state_.GyroBias[2];
        state_.GyroBias[2] += Clamp(bias_error, -bias_limit, bias_limit);
    }

    state_.Gyro[0] = sample.gyro[0] - state_.GyroBias[0];
    state_.Gyro[1] = sample.gyro[1] - state_.GyroBias[1];
    state_.Gyro[2] = sample.gyro[2] - state_.GyroBias[2];

    state_.Yaw  = atan2f(2.0f * (state_.q[0] * state_.q[3] + state_.q[1] * state_.q[2]),
                         2.0f * (state_.q[0] * state_.q[0] + state_.q[1] * state_.q[1]) - 1.0f) * kRadToDeg;
    state_.Roll = atan2f(2.0f * (state_.q[0] * state_.q[1] + state_.q[2] * state_.q[3]),
                         2.0f * (state_.q[0] * state_.q[0] + state_.q[3] * state_.q[3]) - 1.0f) * kRadToDeg;
    const float pitch_input = Clamp(-2.0f * (state_.q[1] * state_.q[3] - state_.q[0] * state_.q[2]), -1.0f, 1.0f);
    state_.Pitch = asinf(pitch_input) * kRadToDeg;

    if (state_.Yaw - state_.YawAngleLast > 180.0f) {
        state_.YawRoundCount--;
    } else if (state_.Yaw - state_.YawAngleLast < -180.0f) {
        state_.YawRoundCount++;
    }
    state_.YawTotalAngle = 360.0f * state_.YawRoundCount + state_.Yaw;
    state_.YawAngleLast  = state_.Yaw;
    state_.UpdateCount++;

    SyncLegacyExports(state_);
}

void QuaternionEkf::Observe()
{
    state_.ChiSquare = chi_square_data_[0];
}

void QuaternionEkf::LinearizeAndFade()
{
    float q0 = filter_.xhatminus_data[0];
    float q1 = filter_.xhatminus_data[1];
    float q2 = filter_.xhatminus_data[2];
    float q3 = filter_.xhatminus_data[3];

    const float q_inv_norm = FastInverseSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    if (q_inv_norm > 0.0f) {
        for (uint8_t i = 0; i < 4; i++) {
            filter_.xhatminus_data[i] *= q_inv_norm;
        }
    } else {
        filter_.xhatminus_data[0] = 1.0f;
        filter_.xhatminus_data[1] = 0.0f;
        filter_.xhatminus_data[2] = 0.0f;
        filter_.xhatminus_data[3] = 0.0f;
    }

    q0 = filter_.xhatminus_data[0];
    q1 = filter_.xhatminus_data[1];
    q2 = filter_.xhatminus_data[2];
    q3 = filter_.xhatminus_data[3];

    filter_.F_data[4]  =  q1 * state_.dt * 0.5f;
    filter_.F_data[5]  =  q2 * state_.dt * 0.5f;
    filter_.F_data[10] = -q0 * state_.dt * 0.5f;
    filter_.F_data[11] =  q3 * state_.dt * 0.5f;
    filter_.F_data[16] = -q3 * state_.dt * 0.5f;
    filter_.F_data[17] = -q0 * state_.dt * 0.5f;
    filter_.F_data[22] =  q2 * state_.dt * 0.5f;
    filter_.F_data[23] = -q1 * state_.dt * 0.5f;

    filter_.P_data[28] /= state_.lambda;
    filter_.P_data[35] /= state_.lambda;

    if (filter_.P_data[28] > config_.bias_covariance_limit)
        filter_.P_data[28] = config_.bias_covariance_limit;
    if (filter_.P_data[35] > config_.bias_covariance_limit)
        filter_.P_data[35] = config_.bias_covariance_limit;
}

void QuaternionEkf::SetObservationJacobian()
{
    const float dq0 = 2.0f * filter_.xhatminus_data[0];
    const float dq1 = 2.0f * filter_.xhatminus_data[1];
    const float dq2 = 2.0f * filter_.xhatminus_data[2];
    const float dq3 = 2.0f * filter_.xhatminus_data[3];

    memset(filter_.H_data, 0, sizeof(float) * filter_.zSize * filter_.xhatSize);

    // H = dh/dq for gravity direction observation model
    filter_.H_data[0]  = -dq2;
    filter_.H_data[1]  =  dq3;
    filter_.H_data[2]  = -dq0;
    filter_.H_data[3]  =  dq1;

    filter_.H_data[6]  =  dq1;
    filter_.H_data[7]  =  dq0;
    filter_.H_data[8]  =  dq3;
    filter_.H_data[9]  =  dq2;

    filter_.H_data[12] =  dq0;
    filter_.H_data[13] = -dq1;
    filter_.H_data[14] = -dq2;
    filter_.H_data[15] =  dq3;
}

void QuaternionEkf::UpdatePosterior()
{
    constexpr uint8_t nx = 6;
    constexpr uint8_t nz = 3;

    auto fallback_to_prediction = [this]() {
        memcpy(filter_.xhat_data, filter_.xhatminus_data, sizeof(float) * nx);
        memcpy(filter_.P_data, filter_.Pminus_data, sizeof(float) * nx * nx);
        filter_.SkipEq5 = false;
    };

    const float q0 = filter_.xhatminus_data[0];
    const float q1 = filter_.xhatminus_data[1];
    const float q2 = filter_.xhatminus_data[2];
    const float q3 = filter_.xhatminus_data[3];

    auto H_m      = Eigen::Map<Eigen::MatrixXf>(filter_.H_data, nz, nx);
    auto Pminus_m = Eigen::Map<Eigen::MatrixXf>(filter_.Pminus_data, nx, nx);
    auto R_m      = Eigen::Map<Eigen::MatrixXf>(filter_.R_data, nz, nz);
    auto z_m      = Eigen::Map<Eigen::VectorXf>(filter_.z_data, nz);
    auto HT_m     = Eigen::Map<Eigen::MatrixXf>(filter_.HT_data, nx, nz);

    HT_m = H_m.transpose();

    // PHT = Pminus * H^T (固定尺寸 6×3，Eigen 编译期展开)
    Eigen::Matrix<float, 6, 3> PHT = Pminus_m * HT_m;

    // S = H * PHT + R (固定 3×3)
    Eigen::Matrix<float, 3, 3> S = H_m * PHT + R_m;

    // 先验四元数预测的加速度方向 h(xhatminus)
    Eigen::Vector3f h_pred;
    h_pred(0) = 2.0f * (q1 * q3 - q0 * q2);
    h_pred(1) = 2.0f * (q0 * q1 + q2 * q3);
    h_pred(2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

    for (uint8_t i = 0; i < 3; i++) {
        state_.OrientationCosine[i] = acosf(Clamp(fabsf(h_pred(i)), 0.0f, 1.0f));
    }

    // 观测残差
    Eigen::Vector3f innov = z_m - h_pred;

    // 固定 3×3，Eigen 用 cofactor 公式求逆（比 LLT 快）
    Eigen::Matrix3f S_inv = S.inverse();

    // 卡方值 = innov^T * S^-1 * innov
    float chi2 = innov.dot(S_inv * innov);
    chi_square_data_[0] = chi2;
    state_.ChiSquare = chi2;

    if (chi2 < 0.5f * state_.ChiSquareTestThreshold) {
        state_.ConvergeFlag = true;
    }

    if (chi2 > state_.ChiSquareTestThreshold && state_.ConvergeFlag) {
        if (state_.StableFlag) {
            state_.ErrorCount++;
        } else {
            state_.ErrorCount = 0;
        }

        if (state_.ErrorCount > config_.divergence_limit) {
            state_.ConvergeFlag = false;
            filter_.SkipEq5 = false;
        } else {
            memcpy(filter_.xhat_data, filter_.xhatminus_data, sizeof(float) * nx);
            memcpy(filter_.P_data, filter_.Pminus_data, sizeof(float) * nx * nx);
            filter_.SkipEq5 = true;
            return;
        }
    } else {
        if (chi2 > 0.1f * state_.ChiSquareTestThreshold && state_.ConvergeFlag) {
            state_.AdaptiveGainScale =
                (state_.ChiSquareTestThreshold - chi2) / (0.9f * state_.ChiSquareTestThreshold);
        } else {
            state_.AdaptiveGainScale = 1.0f;
        }
        state_.ErrorCount = 0;
        filter_.SkipEq5 = false;
    }

    // K = PHT * S^-1
    auto K_m = Eigen::Map<Eigen::MatrixXf>(filter_.K_data, nx, nz);
    K_m = PHT * S_inv;

    K_m *= state_.AdaptiveGainScale;

    // 根据姿态方向削弱零偏通道增益
    for (uint8_t i = 4; i < 6; i++) {
        for (uint8_t j = 0; j < 3; j++) {
            filter_.K_data[i * 3 + j] *= state_.OrientationCosine[i - 4] / 1.5707963f;
        }
    }

    // 卡尔曼修正: xhat = xhatminus + K * innov (固定 6×1)
    Eigen::Matrix<float, 6, 1> correction = K_m * innov;

    if (state_.ConvergeFlag) {
        for (uint8_t i = 4; i < 6; i++) {
            const float limit = config_.bias_correction_limit * state_.dt;
            if (correction(i) > limit)  correction(i) = limit;
            if (correction(i) < -limit) correction(i) = -limit;
        }
    }

    Eigen::Map<Eigen::VectorXf>(filter_.xhat_data, nx) =
        Eigen::Map<Eigen::VectorXf>(filter_.xhatminus_data, nx) + correction;
}

void QuaternionEkf::ObserveHook(filter::KalmanFilter *kf)
{
    static_cast<QuaternionEkf *>(kf->UserData())->Observe();
}

void QuaternionEkf::LinearizeHook(filter::KalmanFilter *kf)
{
    static_cast<QuaternionEkf *>(kf->UserData())->LinearizeAndFade();
}

void QuaternionEkf::SetHHook(filter::KalmanFilter *kf)
{
    static_cast<QuaternionEkf *>(kf->UserData())->SetObservationJacobian();
}

void QuaternionEkf::UpdatePosteriorHook(filter::KalmanFilter *kf)
{
    static_cast<QuaternionEkf *>(kf->UserData())->UpdatePosterior();
}

float QuaternionEkf::FastInverseSqrt(float x)
{
    if (x <= 0.0f) return 0.0f;

    const float half_x = 0.5f * x;
    union { float f; uint32_t i; } value = { x };
    value.i = 0x5f375a86u - (value.i >> 1);
    value.f = value.f * (1.5f - (half_x * value.f * value.f));
    return value.f;
}

QuaternionEkf& GlobalQuaternionEkf()
{
    static QuaternionEkf instance {};
    return instance;
}

} // namespace alg::attitude

QEKF_INS_t QEKF_INS {};
float      chiSquare = 0.0f;
float      ChiSquareTestThreshold = 0.0f;

void IMU_QuaternionEKF_Init(float process_noise1, float process_noise2, float measure_noise, float lambda, float lpf)
{
    alg::attitude::QuaternionEkf::Config config {};
    config.process_noise_quaternion = process_noise1;
    config.process_noise_bias       = process_noise2;
    config.measurement_noise        = measure_noise;
    config.fading_factor            = lambda;
    config.accel_lpf_coefficient    = lpf;

    auto& ekf = alg::attitude::GlobalQuaternionEkf();
    ekf.Init(config);
}

void IMU_QuaternionEKF_Update(float gx, float gy, float gz,
                              float ax, float ay, float az, float dt)
{
    Sample sample {};
    sample.gyro[0] = gx;
    sample.gyro[1] = gy;
    sample.gyro[2] = gz;
    sample.accel[0] = ax;
    sample.accel[1] = ay;
    sample.accel[2] = az;
    sample.dt = dt;

    auto& ekf = alg::attitude::GlobalQuaternionEkf();
    ekf.Update(sample);
}
