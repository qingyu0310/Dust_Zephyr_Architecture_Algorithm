/**
 * @file lqr.hpp
 * @author qingyu
 * @brief 模板化离散 LQR 控制器
 * @version 0.1
 * @date 2026-09-25
 *
 * @par 使用示例
 * @code
 *   alg::lqr::Lqr<2, 1> lqr;
 *
 *   Eigen::Matrix2f A;
 *   A << 1.0f, dt,
 *        0.0f, 1.0f;
 *
 *   Eigen::Matrix<float, 2, 1> B;
 *   B << 0.5f * dt * dt,
 *        dt;
 *
 *   Eigen::Matrix2f Q = Eigen::Matrix2f::Identity();
 *   Eigen::Matrix<float, 1, 1> R;
 *   R << 0.1f;
 *
 *   lqr.Init(A, B, Q, R);
 *   lqr.Solve();
 *
 *   Eigen::Vector2f x;
 *   x << pos, vel;
 *   float u = lqr.Calc(x)(0);
 * @endcode
 */

#pragma once

#pragma message "Compiling Algorithm/Controller/Lqr"

#include <stdint.h>
#include <Eigen/Dense>

namespace alg::lqr {

/**
 * @brief 固定维度离散 LQR 控制器
 *
 * 求解离散代数 Riccati 方程：
 *   P = A^T P A - A^T P B (R + B^T P B)^-1 B^T P A + Q
 *
 * 反馈增益：
 *   K = (R + B^T P B)^-1 B^T P A
 *
 * 控制律：
 *   u = -K * (x - ref) + uff
 *
 * @tparam nx  状态维数
 * @tparam nu  控制输入维数
 */
template <uint8_t nx, uint8_t nu>
class Lqr final
{
public:
    using State    = Eigen::Matrix<float, nx, 1>;
    using Control  = Eigen::Matrix<float, nu, 1>;
    using StateMat = Eigen::Matrix<float, nx, nx>;
    using CtrlMat  = Eigen::Matrix<float, nx, nu>;
    using CostQ    = Eigen::Matrix<float, nx, nx>;
    using CostR    = Eigen::Matrix<float, nu, nu>;
    using Gain     = Eigen::Matrix<float, nu, nx>;

    struct Config {
        uint16_t maxIter = 80;      // Riccati 最大迭代次数
        float eps        = 1e-5f;   // 收敛阈值，比较 P 的最大元素变化
    };

    Lqr() = default;

    Lqr(const StateMat& A, const CtrlMat& B, const CostQ& Q, const CostR& R, const Config& cfg = Config{})
    {
        Init(A, B, Q, R, cfg);
    }

    /**
     * @brief 初始化系统模型和权重
     * @param A    离散状态转移矩阵
     * @param B    离散控制输入矩阵
     * @param Q    状态权重，需半正定
     * @param R    控制权重，需正定
     * @param cfg  迭代配置
     */
    void Init(const StateMat& A, const CtrlMat& B, const CostQ& Q, const CostR& R, const Config& cfg = Config{})
    {
        A_  = A;
        AT_ = A_.transpose();
        B_  = B;
        BT_ = B_.transpose();
        Q_  = Q;
        R_  = R;
        P_  = Q_;
        K_  .setZero();
        ref_.setZero();
        uff_.setZero();
        cfg_ 	 = cfg;
        solved_  = false;
        iter_ 	 = 0;
        lastErr_ = 0.0f;
    }

    /**
     * @brief 迭代求解离散 LQR 增益
     * @return true 表示收敛并更新 K；false 表示矩阵分解失败或未在 maxIter 内收敛
     */
    bool Solve()
    {
        StateMat P = P_;

        for (uint16_t i = 0; i < cfg_.maxIter; i++) 
		{
            CostR S;
            S.noalias() = R_ + BT_ * P * B_;
            auto solver = S.ldlt();
            if (solver.info() != Eigen::Success) {
                solved_ = false;
                return false;
            }

            const Gain K = solver.solve(BT_ * P * A_);

            StateMat nextP;
            nextP.noalias() = AT_ * P * A_;
            nextP.noalias() -= AT_ * P * B_ * K;
            nextP += Q_;
            MakeSymmetric(nextP);

            lastErr_ = (nextP - P).cwiseAbs().maxCoeff();
            P = nextP;
            iter_ = static_cast<uint16_t>(i + 1U);

            if (lastErr_ <= cfg_.eps) {
                P_ = P;
                K_ = K;
                solved_ = true;
                return true;
            }
        }

        P_ = P;
        UpdateGain();
        solved_ = false;
        return false;
    }

    /**
     * @brief 直接设置反馈增益，适合离线整定后烧录固定 K
     * @param K  控制反馈增益
     */
    void SetGain(const Gain& K)
    {
        K_ = K;
        solved_ = true;
    }

    void SetRef(const State& ref) { ref_ = ref; }
    void SetFeedForward(const Control& uff) { uff_ = uff; }

    /**
     * @brief 计算控制量
     * @param x  当前状态
     * @return 控制输入 u
     */
    Control Calc(const State& x) const
    {
        Control u;
        u.noalias() = -K_ * (x - ref_);
        u += uff_;
        return u;
    }

    /**
     * @brief 计算控制量，同时传入参考状态
     */
    Control Calc(const State& x, const State& ref)
    {
        ref_ = ref;
        return Calc(x);
    }

    const Gain& 	GetK() 		 const { return K_; }
    const StateMat& GetP() 		 const { return P_; }
    const State& 	GetRef() 	 const { return ref_; }
    bool     		IsSolved() 	 const { return solved_; }
    uint16_t 		GetIter() 	 const { return iter_; }
    float 			GetLastErr() const { return lastErr_; }

private:
    Config 		cfg_{};

    StateMat 	A_{};
    StateMat 	AT_{};
    CtrlMat 	B_{};
    Eigen::Matrix<float, nu, nx> BT_{};
    CostQ 		Q_{};
    CostR 		R_{};
    StateMat 	P_{};
    Gain 		K_{};

    State 		ref_{};
    Control 	uff_{};

    bool 		solved_  = false;
    uint16_t 	iter_ 	 = 0;
    float 		lastErr_ = 0.0f;

    void UpdateGain()
    {
        CostR S;
        S.noalias() = R_ + BT_ * P_ * B_;
        auto solver = S.ldlt();
        if (solver.info() != Eigen::Success) return;

        K_ = solver.solve(BT_ * P_ * A_);
    }

    template<typename Mat>
    static inline void MakeSymmetric(Mat& M)
    {
        M = 0.5f * (M + M.transpose()).eval();
    }
};

} // namespace alg::lqr
