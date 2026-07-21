/**
 * @file rls.hpp
 * @author qingyu
 * @brief 递归最小二乘(RLS)自适应滤波器 — 模板版本，编译期分配内存
 * @version 0.3
 * @date 2026-06-27
 *
 * @par 算法原理
 *      标准 SISO-RLS，符号约定:
 *      - x_k: n×1 待估参数（权重）
 *      - α_k: 1×n 输入行向量
 *      - β_k: 1×1 实际输出（标量）
 *      - P_k: n×n 协方差矩阵
 *      - λ:   遗忘因子
 *
 *      递归步骤:
 *      ① β̂_k = α_k · x_{k-1}                           预测输出
 *      ② e_k  = β_k - β̂_k                              先验误差
 *      ③ K_k  = P_{k-1}·α_kᵀ / (λ + α_k·P_{k-1}·α_kᵀ)  增益
 *      ④ x_k  = x_{k-1} + K_k · e_k                    参数更新
 *      ⑤ P_k  = (P_{k-1} - K_k·α_k·P_{k-1}) / λ        协方差更新
 *
 * @tparam n  输入/参数维度（编译期分配）
 *
 * @par 使用示例
 * @code
 *   // 2 输入，遗忘因子 0.99，初始协方差 10.0
 *   alg::rls::RLS<2> filter(0.99f, 10.0f);
 *
 *   float α[2] = {1.0f, 2.0f};
 *   float β    = 5.0f;
 *   filter.Update(α, β);
 *   float k1 = filter.GetWeights()[0];
 *   float k2 = filter.GetWeights()[1];
 * @endcode
 */

#pragma once

#pragma message "Compiling Algorithm/Identify/Rls"

#include <stdint.h>
#include <string.h>
#include <Eigen/Dense>

namespace alg::rls {

template <uint8_t n>
class RLS final
{
public:
    /**
     * @brief 构造 RLS 辨识器
     * @param lambda  遗忘因子 (0 < λ ≤ 1)，越小对过去遗忘越快
     * @param p_init  协方差矩阵 P 的初始对角值
     *
     * @note P 初始值建议:
     *       - SNR 较高 (信号可信): p_init = 1~10
     *       - SNR 较低 (噪声大):   p_init = 100~1000
     */
    RLS(float lambda, float p_init) : lambda_(lambda)
    {
        for (uint16_t i = 0; i < n * n; i += (n + 1)) {
            P_k[i] = p_init;
        }
    }

    /**
     * @brief 执行一次 RLS 迭代
     * @param alpha  输入行向量 α_k，长度 n
     * @param beta   期望输出标量 β_k
     *
     * @par 算法步骤
     *      ① β̂_k  = α_k · x_k              预测输出
     *      ② e_k  = β_k - β̂_k              先验误差
     *      ③ K_k  = P·α_kᵀ / (λ + α·P·αᵀ)  卡尔曼增益
     *      ④ x_k += K_k · e_k              参数更新
     *      ⑤ P_k  = (P - K·α·P) / λ        协方差更新
     */
    /**
     * @brief 执行一次 RLS 迭代
     * @param alpha           输入行向量 α_k，长度 n
     * @param beta            期望输出标量 β_k
     * @param lambda_override  可选：覆盖内部 lambda_，-1 表示使用默认值
     *
     * @note lambda_override 用于异步帧场景：每帧传入 exp(-Δt / T_forget)
     *       约束: 0 < lambda ≤ 1，否则协方差更新可能发散
     */
    void Update(const float (&alpha)[n], float beta, float lambda_override = -1.0f)
    {
        // 拷贝输入
        memcpy(alpha_k, alpha, n * sizeof(float));
        beta_k = beta;

        // 确定本次使用的 lambda
        float lambda = (lambda_override > 0.0f) ? lambda_override : lambda_;

        // Eigen::Map 包装内部 float 数组，零拷贝
        Eigen::Map<Eigen::VectorXf>            x(x_k, n);       // n×1 待估参数
        Eigen::Map<const Eigen::RowVectorXf>   a(alpha_k, n);   // 1×n 输入行向量
        Eigen::Map<Eigen::MatrixXf>            P(P_k, n, n);    // n×n 协方差矩阵

        // 预测 beta_hat_k = alpha_k . x_k
        beta_hat_k = a.dot(x);

        // 误差 e_k = beta_k - beta_hat_k
        e_k = beta_k - beta_hat_k;

        // 增益: K = P . alpha_kT / (lambda + alpha_k . P . alpha_kT)
        Eigen::VectorXf K = P * a.transpose();
        float denom = lambda + a.dot(K);
        if (denom < 1e-12f) denom = 1e-12f;
        K /= denom;

        // 参数更新 x_k += K . e_k
        x += K * e_k;

        // 协方差更新 P = (P - K . alpha . P) / lambda
        P = (P - K * a * P) / lambda;
        P = P.selfadjointView<Eigen::Lower>();
    }

    const float* GetWeights() const { return x_k; }
    void  SetWeights( const float (&w)[n] ) { memcpy(x_k, w, n * sizeof(float)); }
    float GetOutput() const { return beta_hat_k; }
    float GetError()  const { return e_k; }

    /**
     * @brief 重置全部状态（权重清零、协方差重新初始化），用于多轮辨识
     */
    void Reset(float p_init)
    {
        memset(x_k,     0, sizeof(x_k));
        memset(alpha_k, 0, sizeof(alpha_k));
        memset(P_k,     0, sizeof(P_k));

        beta_k     = 0.0f;
        beta_hat_k = 0.0f;
        e_k        = 0.0f;

        for (uint16_t i = 0; i < n * n; i += n + 1) {
            P_k[i] = p_init;
        }
    }

private:
    const float lambda_ = 0.99999f;

    float x_k[n]     {};
    float alpha_k[n] {};
    float beta_k    = 0;
    float P_k[n * n] {};

    float beta_hat_k = 0;
    float e_k        = 0;
};

}; // namespace alg::rls
