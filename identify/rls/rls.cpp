/**
 * @file rls.cpp
 * @author qingyu
 * @brief RLS 递归最小二乘实现
 * @version 0.3
 * @date 2026-06-27
 */
#pragma message "Compiling Identify/Rls"

#include "rls.hpp"
#include <string.h>

namespace alg::rls {

/**
 * @brief 构造 RLS 辨识器
 * @param lambda  遗忘因子 (0 < λ ≤ 1)，越小对过去遗忘越快
 * @param p_init  协方差矩阵 P 的初始对角值
 *
 * @note P 初始值建议:
 *       - SNR 较高 (信号可信): p_init = 1~10
 *       - SNR 较低 (噪声大):   p_init = 100~1000
 */
template <uint8_t n>
RLS<n>::RLS(float lambda, float p_init) : lambda_(lambda)
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
 *      ① β̂_k = α_k · x_k         预测输出
 *      ② e_k  = β_k - β̂_k        先验误差
 *      ③ K_k  = P·α_kᵀ / (λ + α·P·αᵀ)  卡尔曼增益
 *      ④ x_k += K_k · e_k        参数更新
 *      ⑤ P_k  = (P - K·α·P) / λ  协方差更新
 */
template <uint8_t n>
void RLS<n>::Update(const float (&alpha)[n], float beta)
{
    // 拷贝输入
    memcpy(alpha_k, alpha, n * sizeof(float));
    beta_k = beta;

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
    float denom = lambda_ + (a * P * a.transpose())(0);
    K /= denom;

    // 参数更新 x_k += K . e_k
    x += K * e_k;

    // 协方差更新 P = (P - K . alpha . P) / lambda
    P = (P - K * a * P) / lambda_;
}

}; // namespace alg::rls
