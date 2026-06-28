/**
 * @file kalman_ekf.cpp
 * @author qingyu
 * @brief 通用线性 Kalman 滤波器实现 — Eigen 版本
 * @version 0.2
 * @date 2026-06-28
 */
 
#pragma message "Compiling Algorithm/Filter/KalmanEkf"

#include "kalman_ekf.hpp"
#include <string.h>
#include <zephyr/kernel.h>

namespace alg::filter {

KalmanFilter::KalmanFilter()
{
    BindStorage();
}

bool KalmanFilter::Init(uint8_t xhat_size, uint8_t u_size, uint8_t z_size)
{
    if (xhat_size == 0 || xhat_size > kMaxStateSize ||
        u_size > kMaxControlSize || z_size == 0 || z_size > kMaxMeasureSize) {
        return false;
    }

    xhatSize = xhat_size;
    uSize    = u_size;
    zSize    = z_size;

    UseAutoAdjustment   = false;
    MeasurementValidNum = 0;
    SkipEq1 = false;
    SkipEq2 = false;
    SkipEq3 = false;
    SkipEq4 = false;
    SkipEq5 = false;

    User_Func0_f = nullptr;
    User_Func1_f = nullptr;
    User_Func2_f = nullptr;
    User_Func3_f = nullptr;
    User_Func4_f = nullptr;
    User_Func5_f = nullptr;
    User_Func6_f = nullptr;

    ResetStorage();
    return true;
}

void KalmanFilter::Measure()
{
    if (UseAutoAdjustment) {
        AdjustMeasurementMatrices();
    } else {
        memcpy(z_data, MeasuredVector, sizeof(float) * zSize);
        memset(MeasuredVector, 0, sizeof(float) * zSize);
    }

    if (uSize > 0) {
        memcpy(u_data, ControlVector, sizeof(float) * uSize);
    }
}

void KalmanFilter::UpdateXhatMinus()
{
    if (SkipEq1) return;

    const uint8_t nx = xhatSize;
    auto F_m = Eigen::Map<Eigen::MatrixXf>(F_data, nx, nx);
    auto x_m = Eigen::Map<Eigen::VectorXf>(xhat_data, nx);

    if (uSize > 0) {
        auto B_m = Eigen::Map<Eigen::MatrixXf>(B_data, nx, uSize);
        auto u_m = Eigen::Map<Eigen::VectorXf>(u_data, uSize);
        Eigen::Map<Eigen::VectorXf>(xhatminus_data, nx) = F_m * x_m + B_m * u_m;
    } else {
        Eigen::Map<Eigen::VectorXf>(xhatminus_data, nx) = F_m * x_m;
    }
}

void KalmanFilter::UpdatePminus()
{
    if (SkipEq2) return;

    const uint8_t nx = xhatSize;
    auto F_m = Eigen::Map<Eigen::MatrixXf>(F_data, nx, nx);
    auto P_m = Eigen::Map<Eigen::MatrixXf>(P_data, nx, nx);
    auto Q_m = Eigen::Map<Eigen::MatrixXf>(Q_data, nx, nx);

    Eigen::Map<Eigen::MatrixXf>(FT_data, nx, nx) = F_m.transpose();
    Eigen::Map<Eigen::MatrixXf>(Pminus_data, nx, nx) = F_m * P_m * Eigen::Map<Eigen::MatrixXf>(FT_data, nx, nx) + Q_m;
    Eigen::Map<Eigen::MatrixXf>(Pminus_data, nx, nx) = Eigen::Map<Eigen::MatrixXf>(Pminus_data, nx, nx).selfadjointView<Eigen::Lower>();
}

void KalmanFilter::SetGain()
{
    if (SkipEq3) return;

    const uint8_t nx = xhatSize;
    const uint8_t nz = zSize;
    auto Pminus_m = Eigen::Map<Eigen::MatrixXf>(Pminus_data, nx, nx);
    auto H_m      = Eigen::Map<Eigen::MatrixXf>(H_data, nz, nx);
    auto R_m      = Eigen::Map<Eigen::MatrixXf>(R_data, nz, nz);

    Eigen::Map<Eigen::MatrixXf> HT_m(HT_data, nx, nz);
    HT_m = H_m.transpose();

    // S = H * Pminus * H^T + R
    Eigen::Map<Eigen::MatrixXf>(S_data, nz, nz) = (H_m * Pminus_m * HT_m + R_m).eval();
    Eigen::Map<Eigen::MatrixXf>(S_data, nz, nz) = Eigen::Map<Eigen::MatrixXf>(S_data, nz, nz).selfadjointView<Eigen::Lower>();

    // K = Pminus * H^T * S^-1 (LLT)
    auto K_m = Eigen::Map<Eigen::MatrixXf>(K_data, nx, nz);
    K_m = Pminus_m * HT_m * Eigen::Map<Eigen::MatrixXf>(S_data, nz, nz).llt().solve(Eigen::MatrixXf::Identity(nz, nz));
}

void KalmanFilter::UpdateXhat()
{
    if (SkipEq4) return;

    const uint8_t nx = xhatSize;
    const uint8_t nz = zSize;
    auto H_m      = Eigen::Map<Eigen::MatrixXf>(H_data, nz, nx);
    auto xminus_m = Eigen::Map<Eigen::VectorXf>(xhatminus_data, nx);
    auto z_m      = Eigen::Map<Eigen::VectorXf>(z_data, nz);
    auto K_m      = Eigen::Map<Eigen::MatrixXf>(K_data, nx, nz);

    Eigen::Map<Eigen::VectorXf>(temp_vector_data, nz)  = H_m * xminus_m;
    Eigen::Map<Eigen::VectorXf>(temp_vector_data1, nz) = z_m - Eigen::Map<Eigen::VectorXf>(temp_vector_data, nz);
    Eigen::Map<Eigen::VectorXf>(xhat_data, nx) = xminus_m + K_m * Eigen::Map<Eigen::VectorXf>(temp_vector_data1, nz);
}

void KalmanFilter::UpdateCovariance()
{
    if (SkipEq5) return;

    const uint8_t nx = xhatSize;
    Eigen::Map<Eigen::MatrixXf>(temp_matrix_data, nx, nx) = (Eigen::Map<Eigen::MatrixXf>(K_data, nx, zSize) * Eigen::Map<Eigen::MatrixXf>(H_data, zSize, nx) * Eigen::Map<Eigen::MatrixXf>(Pminus_data, nx, nx)).eval();
    Eigen::Map<Eigen::MatrixXf>(P_data, nx, nx) = Eigen::Map<Eigen::MatrixXf>(Pminus_data, nx, nx) - Eigen::Map<Eigen::MatrixXf>(temp_matrix_data, nx, nx);
    Eigen::Map<Eigen::MatrixXf>(P_data, nx, nx) = Eigen::Map<Eigen::MatrixXf>(P_data, nx, nx).selfadjointView<Eigen::Lower>();
}

float *KalmanFilter::Update()
{
    static uint32_t cnt = 0;
    bool log = (++cnt % 1000 == 0);
    uint32_t t0, t1, t2, t3, t4, t5, t6, t7, t8;

    t0 = k_cycle_get_32();
    Measure();
    if (User_Func0_f != nullptr) User_Func0_f(this);
    t1 = k_cycle_get_32();

    UpdateXhatMinus();
    if (User_Func1_f != nullptr) User_Func1_f(this);
    t2 = k_cycle_get_32();

    UpdatePminus();
    if (User_Func2_f != nullptr) User_Func2_f(this);
    t3 = k_cycle_get_32();

    if (MeasurementValidNum != 0 || !UseAutoAdjustment) {
        SetGain();
        if (User_Func3_f != nullptr) User_Func3_f(this);
        t4 = k_cycle_get_32();

        UpdateXhat();
        if (User_Func4_f != nullptr) User_Func4_f(this);
        t5 = k_cycle_get_32();

        UpdateCovariance();
    } else {
        memcpy(xhat_data, xhatminus_data, sizeof(float) * xhatSize);
        memcpy(P_data, Pminus_data, sizeof(float) * xhatSize * xhatSize);
        t4 = t5 = k_cycle_get_32();
    }
    t6 = k_cycle_get_32();

    if (User_Func5_f != nullptr) User_Func5_f(this);
    t7 = k_cycle_get_32();

    for (uint8_t i = 0; i < xhatSize; i++) {
        const uint16_t idx = i * xhatSize + i;
        if (P_data[idx] < StateMinVariance[i]) {
            P_data[idx] = StateMinVariance[i];
        }
    }
    memcpy(FilteredValue, xhat_data, sizeof(float) * xhatSize);
    if (User_Func6_f != nullptr) User_Func6_f(this);
    t8 = k_cycle_get_32();

    if (log) {
        printk("[KF] Measure+F0=%uus  xhat-= %uus  P-= %uus  SetGain+F3=%uus  xhat+= %uus  P= %uus  F5=%uus  final=%uus  total=%uus\n",
            (unsigned)(k_cyc_to_ns_floor64(t1 - t0) / 1000),
            (unsigned)(k_cyc_to_ns_floor64(t2 - t1) / 1000),
            (unsigned)(k_cyc_to_ns_floor64(t3 - t2) / 1000),
            (unsigned)(k_cyc_to_ns_floor64(t4 - t3) / 1000),
            (unsigned)(k_cyc_to_ns_floor64(t5 - t4) / 1000),
            (unsigned)(k_cyc_to_ns_floor64(t6 - t5) / 1000),
            (unsigned)(k_cyc_to_ns_floor64(t7 - t6) / 1000),
            (unsigned)(k_cyc_to_ns_floor64(t8 - t7) / 1000),
            (unsigned)(k_cyc_to_ns_floor64(t8 - t0) / 1000));
    }

    return FilteredValue;
}

void KalmanFilter::BindStorage()
{
    MeasurementMap        = measurement_map_storage_;
    MeasurementDegree     = measurement_degree_storage_;
    MatR_DiagonalElements = mat_r_diagonal_storage_;
    StateMinVariance      = state_min_variance_storage_;
    FilteredValue         = filtered_value_storage_;
    MeasuredVector        = measured_vector_storage_;
    ControlVector         = control_vector_storage_;

    xhat_data         = xhat_data_storage_;
    xhatminus_data    = xhatminus_data_storage_;
    u_data            = u_data_storage_;
    z_data            = z_data_storage_;
    P_data            = p_data_storage_;
    Pminus_data       = pminus_data_storage_;
    F_data            = f_data_storage_;
    FT_data           = ft_data_storage_;
    B_data            = b_data_storage_;
    H_data            = h_data_storage_;
    HT_data           = ht_data_storage_;
    Q_data            = q_data_storage_;
    R_data            = r_data_storage_;
    K_data            = k_data_storage_;
    S_data            = s_data_storage_;
    temp_matrix_data  = temp_matrix_data_storage_;
    temp_matrix_data1 = temp_matrix1_data_storage_;
    temp_vector_data  = temp_vector_data_storage_;
    temp_vector_data1 = temp_vector1_data_storage_;
}

void KalmanFilter::ResetStorage()
{
    memset(measurement_map_storage_, 0, sizeof(measurement_map_storage_));
    memset(measurement_degree_storage_, 0, sizeof(measurement_degree_storage_));
    memset(mat_r_diagonal_storage_, 0, sizeof(mat_r_diagonal_storage_));
    memset(state_min_variance_storage_, 0, sizeof(state_min_variance_storage_));
    memset(temp_index_storage_, 0, sizeof(temp_index_storage_));

    memset(filtered_value_storage_, 0, sizeof(filtered_value_storage_));
    memset(measured_vector_storage_, 0, sizeof(measured_vector_storage_));
    memset(control_vector_storage_, 0, sizeof(control_vector_storage_));

    memset(xhat_data_storage_, 0, sizeof(xhat_data_storage_));
    memset(xhatminus_data_storage_, 0, sizeof(xhatminus_data_storage_));
    memset(u_data_storage_, 0, sizeof(u_data_storage_));
    memset(z_data_storage_, 0, sizeof(z_data_storage_));
    memset(p_data_storage_, 0, sizeof(p_data_storage_));
    memset(pminus_data_storage_, 0, sizeof(pminus_data_storage_));
    memset(f_data_storage_, 0, sizeof(f_data_storage_));
    memset(ft_data_storage_, 0, sizeof(ft_data_storage_));
    memset(b_data_storage_, 0, sizeof(b_data_storage_));
    memset(h_data_storage_, 0, sizeof(h_data_storage_));
    memset(ht_data_storage_, 0, sizeof(ht_data_storage_));
    memset(q_data_storage_, 0, sizeof(q_data_storage_));
    memset(r_data_storage_, 0, sizeof(r_data_storage_));
    memset(k_data_storage_, 0, sizeof(k_data_storage_));
    memset(s_data_storage_, 0, sizeof(s_data_storage_));
    memset(temp_matrix_data_storage_, 0, sizeof(temp_matrix_data_storage_));
    memset(temp_matrix1_data_storage_, 0, sizeof(temp_matrix1_data_storage_));
    memset(temp_vector_data_storage_, 0, sizeof(temp_vector_data_storage_));
    memset(temp_vector1_data_storage_, 0, sizeof(temp_vector1_data_storage_));
}

void KalmanFilter::AdjustMeasurementMatrices()
{
    MeasurementValidNum = 0;

    memcpy(z_data, MeasuredVector, sizeof(float) * zSize);
    memset(MeasuredVector, 0, sizeof(float) * zSize);

    memset(R_data, 0, sizeof(float) * zSize * zSize);
    memset(H_data, 0, sizeof(float) * xhatSize * zSize);

    for (uint8_t i = 0; i < zSize; i++) 
    {
        if (z_data[i] != 0.0f) 
        {
            z_data[MeasurementValidNum] = z_data[i];
            temp_index_storage_[MeasurementValidNum] = i;
            H_data[xhatSize * MeasurementValidNum + MeasurementMap[i] - 1] = MeasurementDegree[i];
            MeasurementValidNum++;
        }
    }

    for (uint8_t i = 0; i < MeasurementValidNum; i++) {
        R_data[i * MeasurementValidNum + i] = MatR_DiagonalElements[temp_index_storage_[i]];
    }
}

} // namespace alg::filter
