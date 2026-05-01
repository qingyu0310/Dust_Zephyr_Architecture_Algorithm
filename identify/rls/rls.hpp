/**
 * @file rls.hpp
 * @author qingyu
 * @brief 递归最小二乘(RLS)自适应滤波器
 * @version 0.1
 * @date 2026-04-28
 *
 * @copyright Copyright (c) 2026
 *
 * @par 算法原理
 *      RLS 通过最小化指数加权误差平方和，递推更新权重向量:
 *      - E(n) = Y(n) - U(n)                        误差
 *      - K(n) = P(n-1)·X(n) / [λ + Xᵀ(n)·P(n-1)·X(n)]  卡尔曼增益
 *      - W(n) = W(n-1) + K(n)·E(n)                 权重更新
 *      - P(n) = [P(n-1) - K(n)·Xᵀ(n)·P(n-1)] / λ  协方差更新
 *
 * @par 使用示例
 * @code
 *   // 3 输入 1 输出，遗忘因子 0.99，初始协方差 10.0
 *   alg::rls::RLS filter(3, 1, 0.99f, 10.0f);
 *
 *   float input[3] = {1.0f, 2.0f, 3.0f};
 *   float desired = 5.0f;
 *   filter.Update(input, &desired);
 *
 *   printk("误差: %f\n", filter.GetError()[0]);
 * @endcode
 */

#pragma once

#include <stdint.h>

namespace alg::rls {

class RLS
{
public:
    /**
     * @brief 构造函数，分配内存并初始化
     * @param x_size  输入向量维度
     * @param y_size  输出向量维度
     * @param lambda  遗忘因子 (0 < λ ≤ 1)，越小对过去遗忘越快
     * @param p_init  协方差矩阵 P 的初始对角值
     *
     * @note P 初始值建议:
     *       - SNR 较高 (信号可信): P_init = 1~10
     *       - SNR 较低 (噪声大):   P_init = 100~1000
     */
    RLS(uint8_t x_size, uint8_t y_size, float lambda, float p_init);

    ~RLS();

    /// 禁止拷贝
    RLS(const RLS&) = delete;
    RLS& operator=(const RLS&) = delete;

    /**
     * @brief 执行一次 RLS 迭代
     * @param input    输入向量 X，长度 x_size
     * @param desired  期望输出 Y，长度 y_size
     *
     * @note 调用前需保证 input / desired 有效；内部自动更新权重与协方差
     */
    void Update(const float* input, const float* desired);

    /** @brief 获取当前权重向量 W (长度 x_size) */
    const float* GetWeights() const { return w_; }

    /** @brief 设置初始权重 W (长度 x_size) */
    void SetWeights(const float* w);

    /** @brief 获取模型估计输出 U (长度 y_size) */
    const float* GetOutput()   const { return u_; }

    /** @brief 获取当前误差 E = Y - U (长度 y_size) */
    const float* GetError()    const { return e_; }

    /** @brief 获取输入维度 */
    uint8_t GetXSize() const { return x_size_; }

    /** @brief 获取输出维度 */
    uint8_t GetYSize() const { return y_size_; }

private:
    /* ---------- 矩阵运算辅助 (私有，避免外部依赖) ---------- */

    /** C[m×n] = A[m×k] × B[k×n] */
    void Multiply(float* c, const float* a, const float* b,
                  uint8_t m, uint8_t k, uint8_t n) const;

    /** C[m×n] = A[m×n] + B[m×n] */
    void Add(float* c, const float* a, const float* b,
             uint8_t rows, uint8_t cols) const;

    /** C[m×n] = A[m×n] - B[m×n] */
    void Subtract(float* c, const float* a, const float* b,
                  uint8_t rows, uint8_t cols) const;

    /** B[n×m] = Aᵀ[m×n] */
    void Transpose(float* b, const float* a,
                   uint8_t rows, uint8_t cols) const;

    /* ---------- 状态 ---------- */

    uint8_t  x_size_;       ///< 输入向量维度
    uint8_t  y_size_;       ///< 输出向量维度
    uint16_t p_size_;       ///< P 矩阵元素个数 (x_size_²)
    float    lambda_;       ///< 遗忘因子

    /* 所有矩阵以 flat float 数组存储，按列优先 (column-major) */

    float* w_;         ///< 权重向量         [x_size_ × 1]
    float* x_;         ///< 当前输入向量     [x_size_ × 1]
    float* xt_;        ///< 输入转置         [1 × x_size_]
    float* p_;         ///< 协方差矩阵       [x_size_ × x_size_]
    float* y_;         ///< 期望输出         [y_size_ × 1]
    float* u_;         ///< 模型估计输出     [y_size_ × 1]
    float* e_;         ///< 误差             [y_size_ × 1]
    float* k_;         ///< 卡尔曼增益       [x_size_ × 1]
    float* kn_;        ///< 增益分子 P·X     [x_size_ × 1]
    float  kd_;        ///< 增益分母标量 λ + Xᵀ·P·X
    float* cache_[3];  ///< 3 块临时空间     [x_size_ × x_size_]
    float* output_;    ///< 权重更新暂存     [x_size_ × 1]
};

}; // namespace alg::rls
