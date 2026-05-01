/**
 * @file rls.cpp
 * @author qingyu
 * @brief 递归最小二乘(RLS)自适应滤波器实现
 * @version 0.1
 * @date 2026-04-28
 *
 * @copyright Copyright (c) 2026
 */

#include "rls.hpp"
#include <string.h>

namespace alg::rls {

RLS::RLS(uint8_t x_size, uint8_t y_size, float lambda, float p_init)
    : x_size_(x_size)
    , y_size_(y_size)
    , p_size_(static_cast<uint16_t>(x_size) * x_size)
    , lambda_(lambda)
{
    /* 分配所有 flat 数组 */
    w_       = new float[x_size_ * y_size_]();
    x_       = new float[x_size_]();
    xt_      = new float[x_size_]();
    p_       = new float[p_size_]();
    y_       = new float[y_size_]();
    u_       = new float[y_size_]();
    e_       = new float[y_size_]();
    k_       = new float[x_size_]();
    kn_      = new float[x_size_]();
    kd_      = 0.0f;
    cache_[0] = new float[p_size_]();
    cache_[1] = new float[p_size_]();
    cache_[2] = new float[p_size_]();
    output_  = new float[x_size_]();

    /* P 初始化为对角阵 P_init × I */
    for (uint16_t i = 0; i < p_size_; i += (x_size_ + 1)) {
        p_[i] = p_init;
    }
}

RLS::~RLS()
{
    delete[] w_;
    delete[] x_;
    delete[] xt_;
    delete[] p_;
    delete[] y_;
    delete[] u_;
    delete[] e_;
    delete[] k_;
    delete[] kn_;
    delete[] cache_[0];
    delete[] cache_[1];
    delete[] cache_[2];
    delete[] output_;
}

void RLS::Update(const float* input, const float* desired)
{
    /* 拷贝输入 / 期望输出到内部缓冲区 */
    memcpy(x_, input,  x_size_ * sizeof(float));
    memcpy(y_, desired, y_size_ * sizeof(float));

    const uint8_t  m = x_size_;   // 简写
    const uint16_t n = p_size_;

    /* ---- 1. 模型输出 U = Wᵀ · X  [y_size_ × 1] ---- */
    for (uint8_t i = 0; i < y_size_; i++) {
        u_[i] = 0.0f;
        for (uint8_t j = 0; j < x_size_; j++) {
            u_[i] += w_[j * y_size_ + i] * x_[j];
        }
    }

    /* ---- 2. 输入转置 XT = Xᵀ  [1 × m] ---- */
    Transpose(xt_, x_, m, 1);

    /* ---- 3. 误差 E = Y - U  [y_size_ × 1] ---- */
    Subtract(e_, y_, u_, y_size_, 1);

    /* ---- 3. 增益分子 KN = P · X  [m × 1] ---- */
    Multiply(kn_, p_, x_, m, m, 1);

    /* ---- 4. C1 = KN · XT = P·X·Xᵀ  [m × m] ---- */
    Multiply(cache_[1], kn_, xt_, m, 1, m);

    /* ---- 5. C2 = C1 · P = P·X·Xᵀ·P  [m × m] ---- */
    Multiply(cache_[2], cache_[1], p_, m, m, m);

    /* ---- 6. C0 = XT · P  [1 × m] ---- */
    Multiply(cache_[0], xt_, p_, 1, m, m);

    /* ---- 7. C1 = C0 · X = Xᵀ·P·X  [1 × 1] (标量) ---- */
    Multiply(cache_[1], cache_[0], x_, 1, m, 1);
    // cache_[1][0] 此时 = Xᵀ·P·X

    /* ---- 8. C0 = λ + Xᵀ·P·X  [1 × 1] ---- */
    cache_[0][0] = lambda_ + cache_[1][0];

    /* ---- 9. 增益分母 KD = 1 / (λ + Xᵀ·P·X)  [1 × 1] ---- */
    kd_ = 1.0f / cache_[0][0];

    /* ---- 10. K = KN · KD = P·X / (λ + Xᵀ·P·X)  [m × 1] ---- */
    for (uint8_t i = 0; i < m; i++) {
        k_[i] = kn_[i] * kd_;
    }

    /* ---- 11. C0 = K · E  [m × y_size_] ---- */
    Multiply(cache_[0], k_, e_, m, 1, y_size_);
    // cache_[0][i] = K[i] · E[0]  (y_size_ = 1 时等价于标量乘)

    /* ---- 12. W = W + K·E  [m × 1] ---- */
    Add(output_, w_, cache_[0], m, 1);

    /* 写回权重 */
    memcpy(w_, output_, m * sizeof(float));

    /* ---- 13. 协方差更新 ---- */
    // C0 = C2 · KD = P·X·Xᵀ·P / (λ + Xᵀ·P·X)  [m × m]
    for (uint16_t i = 0; i < n; i++) {
        cache_[0][i] = cache_[2][i] * kd_;
    }

    // C1 = P - C0  [m × m]
    Subtract(cache_[1], p_, cache_[0], m, m);

    // P = (1/λ) · C1  [m × m]
    const float inv_lambda = 1.0f / lambda_;
    for (uint16_t i = 0; i < n; i++) {
        p_[i] = inv_lambda * cache_[1][i];
    }
}

void RLS::SetWeights(const float* w)
{
    memcpy(w_, w, x_size_ * y_size_ * sizeof(float));
}

void RLS::Multiply(float* c, const float* a, const float* b,
                   uint8_t m, uint8_t k, uint8_t n) const
{
    /* C[m×n] = A[m×k] × B[k×n]  (列优先) */
    for (uint8_t col = 0; col < n; col++) {
        for (uint8_t row = 0; row < m; row++) {
            float sum = 0.0f;
            for (uint8_t i = 0; i < k; i++) {
                sum += a[i * m + row] * b[col * k + i];
            }
            c[col * m + row] = sum;
        }
    }
}

void RLS::Add(float* c, const float* a, const float* b,
              uint8_t rows, uint8_t cols) const
{
    const uint16_t len = static_cast<uint16_t>(rows) * cols;
    for (uint16_t i = 0; i < len; i++) {
        c[i] = a[i] + b[i];
    }
}

void RLS::Subtract(float* c, const float* a, const float* b,
                   uint8_t rows, uint8_t cols) const
{
    const uint16_t len = static_cast<uint16_t>(rows) * cols;
    for (uint16_t i = 0; i < len; i++) {
        c[i] = a[i] - b[i];
    }
}

void RLS::Transpose(float* b, const float* a,
                    uint8_t rows, uint8_t cols) const
{
    for (uint8_t r = 0; r < rows; r++) {
        for (uint8_t c = 0; c < cols; c++) {
            b[r * cols + c] = a[c * rows + r];
        }
    }
}

}; // namespace alg::rls
