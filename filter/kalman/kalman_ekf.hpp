/**
 * @file kalman_ekf.hpp
 * @author qingyu
 * @brief 通用线性 Kalman 滤波器 — 使用 Eigen 进行矩阵运算
 * @version 0.2
 * @date 2026-06-28
 */

#pragma once

#include <stdint.h>
#include <Eigen/Dense>

namespace alg::filter {

/**
 * @brief 通用线性 Kalman 滤波器。
 *
 * 状态向量:
 * - xhat: 当前后验估计。
 * - xhatminus: 当前先验预测。
 *
 * 常用矩阵:
 * - F: 状态转移矩阵。
 * - B: 控制输入矩阵，uSize 为 0 时不参与计算。
 * - H: 观测矩阵。
 * - Q/R: 过程噪声和观测噪声。
 * - P/Pminus: 后验/先验协方差。
 * - K: Kalman 增益。
 *
 * 所有矩阵以固定大小 float 数组存储，通过 Eigen::Map 运算。
 * 上层可直接读写 xhat_data / F_data / P_data 等原始指针。
 */
class KalmanFilter final
{
public:
    using Hook = void (*)(KalmanFilter *kf);

    static constexpr uint8_t kMaxStateSize   = 16;
    static constexpr uint8_t kMaxControlSize = 8;
    static constexpr uint8_t kMaxMeasureSize = 8;

    KalmanFilter();

    bool   Init(uint8_t xhat_size, uint8_t u_size, uint8_t z_size);
    void   Measure();
    void   UpdateXhatMinus();
    void   UpdatePminus();
    void   SetGain();
    void   UpdateXhat();
    void   UpdateCovariance();
    float *Update();

    void  SetUserData(void *user_data) { user_data_ = user_data; }
    void *UserData() const { return user_data_; }

    /* 对外暴露指针，兼容旧工程调参/Hook 写法。 */
    float   *FilteredValue         = nullptr;
    float   *MeasuredVector        = nullptr;
    float   *ControlVector         = nullptr;
    uint8_t  xhatSize              = 0;
    uint8_t  uSize                 = 0;
    uint8_t  zSize                 = 0;
    bool     UseAutoAdjustment     = false;
    uint8_t  MeasurementValidNum   = 0;
    uint8_t *MeasurementMap        = nullptr;
    float   *MeasurementDegree     = nullptr;
    float   *MatR_DiagonalElements = nullptr;
    float   *StateMinVariance      = nullptr;

    bool SkipEq1 = false;
    bool SkipEq2 = false;
    bool SkipEq3 = false;
    bool SkipEq4 = false;
    bool SkipEq5 = false;

    /* 原始数据指针，用于直接读写矩阵元素（兼容 quaternion_ekf）。 */
    float *xhat_data         = nullptr;
    float *xhatminus_data    = nullptr;
    float *u_data            = nullptr;
    float *z_data            = nullptr;
    float *P_data            = nullptr;
    float *Pminus_data       = nullptr;
    float *F_data            = nullptr;
    float *FT_data           = nullptr;
    float *B_data            = nullptr;
    float *H_data            = nullptr;
    float *HT_data           = nullptr;
    float *Q_data            = nullptr;
    float *R_data            = nullptr;
    float *K_data            = nullptr;
    float *S_data            = nullptr;
    float *temp_matrix_data  = nullptr;
    float *temp_matrix_data1 = nullptr;
    float *temp_vector_data  = nullptr;
    float *temp_vector_data1 = nullptr;

    /* Update() 各阶段的扩展点。 */
    Hook User_Func0_f = nullptr;
    Hook User_Func1_f = nullptr;
    Hook User_Func2_f = nullptr;
    Hook User_Func3_f = nullptr;
    Hook User_Func4_f = nullptr;
    Hook User_Func5_f = nullptr;
    Hook User_Func6_f = nullptr;

private:
    static constexpr uint16_t kMaxStateSquare   = kMaxStateSize * kMaxStateSize;
    static constexpr uint16_t kMaxStateControl  = kMaxStateSize * kMaxControlSize;
    static constexpr uint16_t kMaxStateMeasure  = kMaxStateSize * kMaxMeasureSize;
    static constexpr uint16_t kMaxMeasureState  = kMaxMeasureSize * kMaxStateSize;
    static constexpr uint16_t kMaxMeasureSquare = kMaxMeasureSize * kMaxMeasureSize;
    static constexpr uint8_t  kMaxVectorSize    = kMaxStateSize;

    void BindStorage();
    void ResetStorage();
    void AdjustMeasurementMatrices();

    void *user_data_ = nullptr;

    uint8_t measurement_map_storage_[kMaxMeasureSize] {};
    float   measurement_degree_storage_[kMaxMeasureSize] {};
    float   mat_r_diagonal_storage_[kMaxMeasureSize] {};
    float   state_min_variance_storage_[kMaxStateSize] {};
    uint8_t temp_index_storage_[kMaxMeasureSize] {};

    float filtered_value_storage_[kMaxStateSize] {};
    float measured_vector_storage_[kMaxMeasureSize] {};
    float control_vector_storage_[kMaxControlSize] {};

    float xhat_data_storage_[kMaxStateSize] {};
    float xhatminus_data_storage_[kMaxStateSize] {};
    float u_data_storage_[kMaxControlSize] {};
    float z_data_storage_[kMaxMeasureSize] {};
    float p_data_storage_[kMaxStateSquare] {};
    float pminus_data_storage_[kMaxStateSquare] {};
    float f_data_storage_[kMaxStateSquare] {};
    float ft_data_storage_[kMaxStateSquare] {};
    float b_data_storage_[kMaxStateControl] {};
    float h_data_storage_[kMaxMeasureState] {};
    float ht_data_storage_[kMaxStateMeasure] {};
    float q_data_storage_[kMaxStateSquare] {};
    float r_data_storage_[kMaxMeasureSquare] {};
    float k_data_storage_[kMaxStateMeasure] {};
    float s_data_storage_[kMaxStateSquare] {};
    float temp_matrix_data_storage_[kMaxStateSquare] {};
    float temp_matrix1_data_storage_[kMaxStateSquare] {};
    float temp_vector_data_storage_[kMaxVectorSize] {};
    float temp_vector1_data_storage_[kMaxVectorSize] {};
};

} // namespace alg::filter
