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
 *      ① β̂_k = α_k · x_{k-1}                      预测输出
 *      ② e_k  = β_k - β̂_k                         先验误差
 *      ③ K_k  = P_{k-1}·α_kᵀ / (λ + α_k·P_{k-1}·α_kᵀ)  增益
 *      ④ x_k  = x_{k-1} + K_k · e_k                参数更新
 *      ⑤ P_k  = (P_{k-1} - K_k·α_k·P_{k-1}) / λ    协方差更新
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

#include <stdint.h>
#include <Eigen/Dense>

namespace alg::rls {

template <uint8_t n>
class RLS final
{
public:
    RLS(float lambda, float p_init);
    void Update(const float (&alpha)[n], float beta);

    const float* GetWeights() const { return x_k; }
    void SetWeights(const float (&w)[n]) { memcpy(x_k, w, n * sizeof(float)); }
    float GetOutput() const { return beta_hat_k; }
    float GetError()   const { return e_k; }

private:
    const float lambda_;

    float x_k[n] {};
    float alpha_k[n] {};
    float beta_k = 0;
    float P_k[n * n] {};

    float beta_hat_k = 0;
    float e_k        = 0;
};

}; // namespace alg::rls
