/**
 * @file kalman.hpp
 * @author qingyu
 * @brief 模板化通用线性 Kalman 滤波器
 * @version 0.1
 * @date 2026-06-28
 *
 * @par 使用示例（弹簧阻尼系统）
 * @code
 *   //   x[0] = pos          位置
 *   //   x[1] = vel          速度
 *   //   z[0] = pos_meas     位置测量
 *   //
 *   //   A = [1       dt    ]
 *   //       [-k·dt  1-c·dt ]
 *   //
 *   //   H = [1  0]
 *
 *   constexpr float k  = 10.0f;
 *   constexpr float c  = 0.5f;
 *   constexpr float dt = 0.01f;
 *
 *   Eigen::Matrix2f A;
 *   A << 1.0f, dt, -k * dt, 1.0f - c * dt;
 *
 *   Eigen::RowVector2f H;
 *   H << 1.0f, 0.0f;
 *
 *   Kalman<2, 1> kf;
 *   kf.Init(A, H,
 *           Eigen::Matrix2f::Identity() * 0.01f,
 *           Eigen::Matrix<float, 1, 1>::Identity() * 0.1f);
 *
 *   // 每控制周期
 *   Eigen::Matrix<float, 1, 1> z;
 *   z(0) = sensor.Read();
 *   kf.SetZ(z);                 // 设观测值
 *   // kf.SetU(u);              // 设控制输入（nu=0 时不需要）
 *   kf.Predict();
 *   kf.Update();
 *
 *   float pos = kf.GetX()(0);
 *   float vel = kf.GetX()(1);
 * @endcode
 */

#pragma once

#pragma message "Compiling Algorithm/Filter/Kalman"

#include <Eigen/Dense>

namespace alg::filter {

template <uint8_t nx, uint8_t nz, uint8_t nu = 0>
class Kalman
{
public:
    using State   = Eigen::Matrix<float, nx, 1>;    // n×1  后验状态 x
    using Cov     = Eigen::Matrix<float, nx, nx>;   // n×n  协方差矩阵 P
    using Obs     = Eigen::Matrix<float, nz, 1>;    // m×1  观测值 z
    using ObsCov  = Eigen::Matrix<float, nz, nz>;   // m×m  观测噪声 R
    using Ctrl    = Eigen::Matrix<float, nu, 1>;    // p×1  控制输入 u
    using CtrlMat = Eigen::Matrix<float, nx, nu>;   // n×p  控制矩阵 B
    using ObsMat  = Eigen::Matrix<float, nz, nx>;   // m×n  观测矩阵 H
    using Gain    = Eigen::Matrix<float, nx, nz>;   // n×m  卡尔曼增益 K

    inline void SetZ(const Obs& z)        {  z_ = z;  }
    inline void SetU(const Ctrl& u)       {  u_ = u;  }
    inline void SetState(const State& x0) {  x_ = x0;  x_minus_ = x0; }
    inline void Reset() {  x_.setZero();  x_minus_.setZero();  P_ = Cov::Identity();  P_minus_ = P_;  }

    inline const State& GetX() const { return x_; }
    inline const Cov&   GetP() const { return P_; }
    inline const Gain&  GetK() const { return K_; }

    void Init(const Cov&     A,
              const ObsMat&  H  = ObsMat::Zero(),
              const Cov&     Q  = Cov   ::Zero(),
              const ObsCov&  R  = ObsCov::Zero(),
              const Cov&     P0 = Cov   ::Identity(),
              const State&   P_min_diag = State::Zero(),
              const CtrlMat& B  = CtrlMat::Zero())
    {
        A_  = A;
        AT_ = A_.transpose();
        B_  = B;
        H_  = H;
        HT_ = H_.transpose();
        Q_  = Q;
        R_  = R;
        P_  = P0;
        P_min_diag_ = P_min_diag;
        z_.setZero();
        u_.setZero();
        h_.setZero();
        x_.setZero();
        x_minus_.setZero();
        P_minus_ = P_;
        K_.setZero();
        S_.setZero();
    }

    void Predict()
    {
        x_minus_.noalias() = A_ * x_;
        if constexpr (nu > 0) {
            x_minus_.noalias() += B_ * u_;
        }
        P_minus_.noalias() = A_ * P_ * AT_;
        P_minus_ += Q_;
        MakeSymmetric(P_minus_);
    }

    void Predict(const Ctrl& u)
    {
        u_ = u;
        Predict();
    }

    void Update()
    {
        // 创新协方差 S = H * P_minus * H^T + R
        S_.noalias() = H_ * P_minus_ * HT_;
        S_ += R_;
        MakeSymmetric(S_);

        // 卡尔曼增益 K = P_minus * H^T * S^-1
        auto solver = S_.ldlt();
        if (solver.info() != Eigen::Success) return;
        K_.noalias() = P_minus_ * HT_;
        K_ = solver.solve(K_.transpose()).transpose();

        // 状态修正 x = x_minus + K * (z - H * x_minus)
        h_.noalias() = H_ * x_minus_;
        x_ = x_minus_ + K_ * (z_ - h_);

        // 协方差修正 P = (I - K*H) * Pminus
        P_.noalias() = P_minus_ - K_ * H_ * P_minus_;
        MakeSymmetric(P_);

        ClampCovariance();
    }

    void Update(const Obs& z)
    {
        z_ = z;
        Update();
    }

private:
    Cov     A_;         // 状态转移  x_minus = A_ * x_
    Cov     AT_;        // A_^T
    CtrlMat B_;         // 控制输入  n×p
    Cov     Q_;         // 过程噪声
    ObsCov  R_;         // 观测噪声
    ObsMat  H_;         // 观测矩阵
    Eigen::Matrix<float, nx, nz> HT_;   // H_^T

    Obs     z_;         // 实际观测
    Ctrl    u_;         // 控制输入

    Obs     h_;         // 预测观测  h_ = H_ * x_minus_

    State   x_;         // 后验状态
    State   x_minus_;   // 先验状态
    Cov     P_;         // 后验协方差
    Cov     P_minus_;   // 先验协方差
    Gain    K_;         // 卡尔曼增益

    ObsCov  S_;         // 创新协方差  H_ * P_minus_ * H_^T + R_

    State P_min_diag_;

    template<typename Mat>
    static inline void MakeSymmetric(Mat& M)
    {
        M.template triangularView<Eigen::Upper>() = M.transpose();
    }

    void ClampCovariance()
    {
        P_.diagonal() = P_.diagonal().cwiseMax(P_min_diag_);
    }
};

} // namespace alg::filter
