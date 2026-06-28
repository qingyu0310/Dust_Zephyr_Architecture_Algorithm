/**
 * @file quaternion.cpp
 * @author qingyu
 * @brief 基于四元数的 IMU 姿态扩展 Kalman 滤波器实现
 * @version 0.1
 * @date 2026-06-28
 */

#pragma message "Compiling Algorithm/Filter/QuaternionEkf"

#include "quaternion.hpp"
#include "imu.hpp"
#include <zephyr/kernel.h>
#include <math.h>
#include <algorithm>
#include <Eigen/Dense>

namespace {

constexpr float kRadToDeg = 57.295779513f;
constexpr float kYawBiasStaticThreshold = 0.02f;

} // namespace

namespace alg::attitude {

/**
 * @brief 快速倒数平方根
 * @param x  输入
 * @return 1/√x
 */
float QuaternionEkf::FastInverseSqrt(float x)
{
    if (x <= 0.0f) return 0.0f;

    const float half_x = 0.5f * x;
    union { float f; uint32_t i; } value = { x };
    value.i = 0x5f375a86u - (value.i >> 1);
    value.f = value.f * (1.5f - (half_x * value.f * value.f));
    return value.f;
}

/**
 * @brief 单位化 EKF 状态中的四元数
 */
void QuaternionEkf::NormalizeQuaternion(EkfSt& x)
{
    const float q_inv = FastInverseSqrt(x(0) * x(0) + x(1) * x(1) + x(2) * x(2) + x(3) * x(3));
    if (q_inv <= 0.0f) {
        x(0) = 1.0f;
        x(1) = 0.0f;
        x(2) = 0.0f;
        x(3) = 0.0f;
        return;
    }

    x(0) *= q_inv;
    x(1) *= q_inv;
    x(2) *= q_inv;
    x(3) *= q_inv;
}

/**
 * @brief EKF 系统模型回调
 *        x⁻ = f(x, u),  F = ∂f/∂x
 *
 * 四元数运动学积分（一阶龙格-库塔）：
 *   q' = q + (dt/2) · Ω(ω) · q
 *   零偏随机游走：bg' = bg
 *
 * 控制输入 u = [ω·dt, dt]，包含 dt 信息。
 *
 * @param x_in  当前状态 [q0..q3, bgx, bgy]
 * @param u_in  控制输入 [ωx·dt, ωy·dt, ωz·dt, dt]
 * @param x_out 先验状态
 * @param F_out Jacobian 6×6
 */
void QuaternionEkf::SystemFunc(const EkfSt& x_in, const EKF::Ctrl& u_in, EkfSt& x_out, EkfCov& F_out)
{
    const float q0 = x_in(0), q1 = x_in(1), q2 = x_in(2), q3 = x_in(3);

    const float hx = 0.5f * u_in(0);  // 1/2 ωx·dt
    const float hy = 0.5f * u_in(1);  // 1/2 ωy·dt
    const float hz = 0.5f * u_in(2);  // 1/2 ωz·dt
    const float hdt = 0.5f * u_in(3); // 1/2 dt

    // RK1 四元数积分: q' = q + 0.5 * Omega(w) * q * dt。
    x_out(0) = q0 + (-q1 * hx - q2 * hy - q3 * hz);
    x_out(1) = q1 + ( q0 * hx + q2 * hz - q3 * hy);
    x_out(2) = q2 + ( q0 * hy - q1 * hz + q3 * hx);
    x_out(3) = q3 + ( q0 * hz + q1 * hy - q2 * hx);
    NormalizeQuaternion(x_out);

    // 零偏随机游走
    x_out(4) = x_in(4);
    x_out(5) = x_in(5);

    // Jacobian F：左上 4x4 是姿态积分，右上 4x2 是 bgx/bgy 耦合项。
    F_out.setIdentity();
    F_out(0, 1) = -hx;   F_out(0, 2) = -hy;   F_out(0, 3) = -hz;
    F_out(1, 0) =  hx;   F_out(1, 2) =  hz;   F_out(1, 3) = -hy;
    F_out(2, 0) =  hy;   F_out(2, 1) = -hz;   F_out(2, 3) =  hx;
    F_out(3, 0) =  hz;   F_out(3, 1) =  hy;   F_out(3, 2) = -hx;

    F_out(0, 4) =  q1 * hdt;   F_out(0, 5) =  q2 * hdt;
    F_out(1, 4) = -q0 * hdt;   F_out(1, 5) =  q3 * hdt;
    F_out(2, 4) = -q3 * hdt;   F_out(2, 5) = -q0 * hdt;
    F_out(3, 4) =  q2 * hdt;   F_out(3, 5) = -q1 * hdt;
}

/**
 * @brief EKF 观测模型回调
 *        z_pred = h(x⁻),  H = ∂h/∂x
 *
 * 由四元数预测重力方向（归一化加速度应指向的向量）：
 *   g(q) = [2·(q1·q3 - q0·q2),
 *           2·(q0·q1 + q2·q3),
 *           q0² - q1² - q2² + q3²]
 *
 * @param x_in  先验状态 x⁻
 * @param z_out 预测观测值（3×1 归一化重力方向）
 * @param H_out Jacobian H（3×6，后两列对零偏为 0）
 */
void QuaternionEkf::ObserveFunc(const EkfSt& x_in, EKF::Obs& z_out, EKF::ObsMat& H_out)
{
    const float q0 = x_in(0), q1 = x_in(1), q2 = x_in(2), q3 = x_in(3);

    // 重力方向预测
    z_out(0) = 2.0f * (q1 * q3 - q0 * q2);
    z_out(1) = 2.0f * (q0 * q1 + q2 * q3);
    z_out(2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

    // 观测 Jacobian H = ∂g/∂q   (3×4)，零偏列补零
    const float dq0 = 2.0f * q0, dq1 = 2.0f * q1;
    const float dq2 = 2.0f * q2, dq3 = 2.0f * q3;

    H_out.setZero();
    H_out(0, 0) = -dq2;   H_out(0, 1) =  dq3;   H_out(0, 2) = -dq0;   H_out(0, 3) =  dq1;
    H_out(1, 0) =  dq1;   H_out(1, 1) =  dq0;   H_out(1, 2) =  dq3;   H_out(1, 3) =  dq2;
    H_out(2, 0) =  dq0;   H_out(2, 1) = -dq1;   H_out(2, 2) = -dq2;   H_out(2, 3) =  dq3;
    // H(0:2, 4:5) = 0（零偏不影响重力方向预测）
}

/**
 * @brief 初始化滤波器
 * @param config 配置参数
 */
void QuaternionEkf::Init(const Config& config)
{
    config_ = config;
    if (config_.lambda > 1.0f) config_.lambda = 1.0f;

    state_ = {};
    state_.init    = true;
    state_.Qq      = config_.Qq;
    state_.Qb      = config_.Qb;
    state_.R       = config_.R;
    state_.lambda  = config_.lambda;
    state_.alpha   = config_.alpha;
    state_.chi2_th = config_.chi2_th;

    // EKF 回调
    ekf_.SetSystemFunc(SystemFunc);
    ekf_.SetObserveFunc(ObserveFunc);

    // 初始过程噪声（dt 缩放由 Update 每周期处理）
    EKF::Cov Q = EKF::Cov::Zero();
    Q(0,0) = Q(1,1) = Q(2,2) = Q(3,3) = config_.Qq;
    Q(4,4) = Q(5,5) = config_.Qb;

    // 观测噪声
    EKF::ObsCov R = EKF::ObsCov::Identity() * config_.R;

    // 初始协方差
    EKF::Cov P0;
    P0 << 100000.0f, 0.1f,      0.1f,      0.1f,      0.1f,     0.1f,
          0.1f,      100000.0f, 0.1f,      0.1f,      0.1f,     0.1f,
          0.1f,      0.1f,      100000.0f, 0.1f,      0.1f,     0.1f,
          0.1f,      0.1f,      0.1f,      100000.0f, 0.1f,     0.1f,
          0.1f,      0.1f,      0.1f,      0.1f,      100.0f,   0.1f,
          0.1f,      0.1f,      0.1f,      0.1f,      0.1f,     100.0f;

    ekf_.Init(Q, R, P0);

    // 初始状态 [1, 0, 0, 0, 0, 0]（单位四元数 + 零偏归零）
    EKF::State x0;
    x0 << 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
    ekf_.SetState(x0);

    state_.q[0] = 1.0f;
}

/**
 * @brief 姿态更新
 * @param sample IMU 采样
 */
void QuaternionEkf::Update(const Sample& sample)
{
    if (!state_.init) Init(Config{});

    const float dt = sample.dt;
    if (dt <= 0.0f) return;
    state_.dt = dt;

    // 去偏角速度
    state_.w[0] = sample.gyro[0] - state_.bg[0];
    state_.w[1] = sample.gyro[1] - state_.bg[1];
    state_.w[2] = sample.gyro[2] - state_.bg[2];

    // 控制输入 = [ω·dt, dt]
    EKF::Ctrl u;
    u << state_.w[0] * dt, state_.w[1] * dt, state_.w[2] * dt, dt;
    ekf_.SetU(u);

    // 过程噪声 Q（dt 缩放）
    EKF::Cov Q = EKF::Cov::Zero();
    Q(0,0) = Q(1,1) = Q(2,2) = Q(3,3) = state_.Qq * dt;
    Q(4,4) = Q(5,5) = state_.Qb * dt;
    ekf_.SetQ(Q);

    // 加速度 LPF
    if (state_.upd_cnt == 0) {
        state_.a[0] = sample.accel[0];
        state_.a[1] = sample.accel[1];
        state_.a[2] = sample.accel[2];
    }

    const float den = dt + state_.alpha;
    if (den > 0.0f) 
    {
        const float lpf_a = state_.alpha / den;
        const float lpf_b = dt / den;
        state_.a[0] = state_.a[0] * lpf_a + sample.accel[0] * lpf_b;
        state_.a[1] = state_.a[1] * lpf_a + sample.accel[1] * lpf_b;
        state_.a[2] = state_.a[2] * lpf_a + sample.accel[2] * lpf_b;
    } 
    else 
    {
        state_.a[0] = sample.accel[0];
        state_.a[1] = sample.accel[1];
        state_.a[2] = sample.accel[2];
    }

    // 归一化加速度 → 观测 z
    const float a_inv = FastInverseSqrt(state_.a[0] * state_.a[0] +
        state_.a[1] * state_.a[1] + state_.a[2] * state_.a[2]);

    EKF::Obs z;
    z << state_.a[0] * a_inv, state_.a[1] * a_inv, state_.a[2] * a_inv;
    ekf_.SetZ(z);

    state_.a_norm = a_inv > 0.0f ? 1.0f / a_inv : 0.0f;

    // 稳定性检测
    {
        const float w_inv = FastInverseSqrt(state_.w[0] * state_.w[0] + 
            state_.w[1] * state_.w[1] + state_.w[2] * state_.w[2]);
        state_.w_norm = w_inv > 0.0f ? 1.0f / w_inv : 0.0f;

        state_.stable =
            state_.w_norm < config_.w_stable_th &&
            state_.a_norm > config_.a_ref - config_.a_tol &&
            state_.a_norm < config_.a_ref + config_.a_tol;
    }

    // EKF 预测 + 渐消 + 校正
    ekf_.Predict();

    // 渐消因子：膨胀零偏协方差
    {
        auto& P = ekf_.PMinus();
        P(4,4) = std::min(P(4,4) / state_.lambda, config_.Pb_limit);
        P(5,5) = std::min(P(5,5) / state_.lambda, config_.Pb_limit);
    }

    ekf_.SetGainScale(state_.k_scale);
    ekf_.Update();
    state_.chi2 = ekf_.GetChi2();

    // 卡方门控
    {
        const float chi2 = state_.chi2;
        const float th   = state_.chi2_th;

        if (chi2 < 0.5f * th) {
            state_.converg = true;
        }

        if (chi2 > th && state_.converg) 
        {
            if (state_.stable) {
                state_.err_cnt++;
            } else {
                state_.err_cnt = 0;
            }

            if (state_.err_cnt > config_.div_limit) {
                state_.converg = false;
                state_.err_cnt = 0;
            } else {
                ekf_.FallbackToPrediction();
            }
        } else {
            if (state_.converg && chi2 > 0.1f * th) {
                state_.k_scale = (th - chi2) / (0.9f * th);
                state_.k_scale = std::clamp(state_.k_scale, 0.0f, 1.0f);
            } else {
                state_.k_scale = 1.0f;
            }
            state_.err_cnt = 0;
        }
    }

    // 提取结果
    {
        auto x = ekf_.GetX();
        NormalizeQuaternion(x);
        ekf_.SetState(x);

        state_.q[0]  = x(0);
        state_.q[1]  = x(1);
        state_.q[2]  = x(2);
        state_.q[3]  = x(3);
        state_.bg[0] = x(4);
        state_.bg[1] = x(5);
    }

    // Z 轴零偏修正（静态时）
    if (state_.stable &&
        std::abs(state_.w[0]) < config_.w_stable_th &&
        std::abs(state_.w[1]) < config_.w_stable_th &&
        std::abs(state_.w[2]) < kYawBiasStaticThreshold)
    {
        const float bias_limit = config_.bias_limit * dt;
        const float bias_error = sample.gyro[2] - state_.bg[2];
        state_.bg[2] += std::clamp(bias_error, -bias_limit, bias_limit);
    }

    // 重新去偏（影响欧拉角计算用的角速度）
    state_.w[0] = sample.gyro[0] - state_.bg[0];
    state_.w[1] = sample.gyro[1] - state_.bg[1];
    state_.w[2] = sample.gyro[2] - state_.bg[2];

    // 欧拉角
    {
        const float q0 = state_.q[0], q1 = state_.q[1];
        const float q2 = state_.q[2], q3 = state_.q[3];

        state_.yaw   = atan2f(2.0f * (q0 * q3 + q1 * q2),
                                    2.0f * (q0 * q0 + q1 * q1) - 1.0f) * kRadToDeg;
        state_.roll  = atan2f(2.0f * (q0 * q1 + q2 * q3),
                              2.0f * (q0 * q0 + q3 * q3) - 1.0f) * kRadToDeg;

        const float pitch_input = std::clamp(2.0f * (q0 * q2 - q1 * q3), -1.0f, 1.0f);
        state_.pitch = asinf(pitch_input) * kRadToDeg;

        // 偏航绕圈计数
        const float yaw_diff = state_.yaw - state_.yaw_prev;
        if (yaw_diff > 180.0f)       state_.yaw_rnd--;
        else if (yaw_diff < -180.0f) state_.yaw_rnd++;
        state_.yaw_sum = 360.0f * state_.yaw_rnd + state_.yaw;
        state_.yaw_prev = state_.yaw;
    }

    state_.upd_cnt++;
}

} // namespace alg::attitude
