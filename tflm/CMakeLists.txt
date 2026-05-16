if(CONFIG_TFLM)

# ============================================================
# TFLM 编译宏 — 禁用动态内存分配、C++ exception 等
# ============================================================
zephyr_library_compile_definitions(
    TF_LITE_STATIC_MEMORY
    TF_LITE_USE_GLOBAL_CMATH_FUNCTIONS
    TF_LITE_USE_GLOBAL_MIN
    TF_LITE_USE_GLOBAL_MAX
    TF_LITE_DISABLE_X86_X87_CUSTOM_OP
    DUNNESS_DISABLE_CORE_NAN_CHECK
    TF_LITE_DISABLE_SCODE_CHECK
)

# ============================================================
# FlatBuffers | gemmlowp | ruy 头文件路径
# ============================================================
zephyr_library_include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/flatbuffers/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/gemmlowp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/ruy
)

# ============================================================
# TFLM 内部 include 路径（不泄露到 app）
# ============================================================
zephyr_library_include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/micro
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/micro/arena_allocator
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/micro/memory_planner
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/micro/tflite_bridge
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/micro/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/c
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/core/c
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/core
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/core/api
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/kernels/internal
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/kernels/internal/reference
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/kernels/internal/optimized
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/schema
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/compiler/mlir/lite/schema
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/compiler/mlir/lite/kernels/internal
    ${CMAKE_CURRENT_SOURCE_DIR}/tensorflow/lite/micro/compression
    ${CMAKE_CURRENT_SOURCE_DIR}/signal/micro/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/port
)

# ============================================================
# TFLM 框架核心
# ============================================================
zephyr_library_sources(
    tensorflow/lite/micro/micro_interpreter.cc
    tensorflow/lite/micro/micro_interpreter_context.cc
    tensorflow/lite/micro/micro_interpreter_graph.cc
    tensorflow/lite/micro/micro_allocator.cc
    tensorflow/lite/micro/micro_allocation_info.cc
    tensorflow/lite/micro/micro_op_resolver.cc
    tensorflow/lite/micro/micro_context.cc
    tensorflow/lite/micro/micro_log.cc
    tensorflow/lite/micro/micro_utils.cc
    tensorflow/lite/micro/micro_time.cc
    tensorflow/lite/micro/memory_helpers.cc
    tensorflow/lite/micro/flatbuffer_utils.cc
    tensorflow/lite/micro/micro_resource_variable.cc
    tensorflow/lite/micro/micro_profiler.cc
    tensorflow/lite/micro/system_setup.cc
    tensorflow/lite/micro/mock_micro_graph.cc
    tensorflow/lite/micro/test_helpers.cc
    tensorflow/lite/micro/test_helper_custom_ops.cc
    tensorflow/lite/micro/recording_micro_allocator.cc
)

# ============================================================
# Arena 分配器
# ============================================================
zephyr_library_sources(
    tensorflow/lite/micro/arena_allocator/single_arena_buffer_allocator.cc
    tensorflow/lite/micro/arena_allocator/persistent_arena_buffer_allocator.cc
    tensorflow/lite/micro/arena_allocator/non_persistent_arena_buffer_allocator.cc
    tensorflow/lite/micro/arena_allocator/recording_single_arena_buffer_allocator.cc
)

# ============================================================
# 内存规划器
# ============================================================
zephyr_library_sources(
    tensorflow/lite/micro/memory_planner/greedy_memory_planner.cc
    tensorflow/lite/micro/memory_planner/linear_memory_planner.cc
    tensorflow/lite/micro/memory_planner/non_persistent_buffer_planner_shim.cc
)

# ============================================================
# TFLite Bridge
# ============================================================
zephyr_library_sources(
    tensorflow/lite/micro/tflite_bridge/flatbuffer_conversions_bridge.cc
    tensorflow/lite/micro/tflite_bridge/micro_error_reporter.cc
)

# ============================================================
# Core C API / Core API / Schema
# ============================================================
zephyr_library_sources(
    tensorflow/lite/core/c/common.cc
    tensorflow/lite/core/api/flatbuffer_conversions.cc
    tensorflow/lite/core/api/tensor_utils.cc
    tensorflow/compiler/mlir/lite/schema/schema_utils.cc
    tensorflow/compiler/mlir/lite/core/api/error_reporter.cc
)

# ============================================================
# Kernels 通用支持
# ============================================================
zephyr_library_sources(
    tensorflow/lite/kernels/kernel_util.cc
    tensorflow/lite/micro/kernels/kernel_runner.cc
    tensorflow/lite/micro/kernels/kernel_util.cc
    tensorflow/lite/micro/kernels/micro_tensor_utils.cc
)

# ============================================================
# Kernels Internal
# ============================================================
zephyr_library_sources(
    tensorflow/lite/kernels/internal/common.cc
    tensorflow/lite/kernels/internal/portable_tensor_utils.cc
    tensorflow/lite/kernels/internal/quantization_util.cc
    tensorflow/lite/kernels/internal/runtime_shape.cc
    tensorflow/lite/kernels/internal/tensor_ctypes.cc
    tensorflow/lite/kernels/internal/tensor_utils.cc
    tensorflow/lite/kernels/internal/reference/comparisons.cc
    tensorflow/lite/kernels/internal/reference/portable_tensor_utils.cc
)

# ============================================================
# TFLM 算子
# ============================================================
zephyr_library_sources(
    tensorflow/lite/micro/kernels/fully_connected.cc
    tensorflow/lite/micro/kernels/fully_connected_common.cc
    tensorflow/lite/micro/kernels/activations.cc
    tensorflow/lite/micro/kernels/activations_common.cc
    tensorflow/lite/micro/kernels/tanh.cc
)

# ============================================================
# 平台适配（Zephyr 通用）
# ============================================================
zephyr_library_sources(
    port/debug_log.cc
)

# ============================================================
# TFLM 封装 API（始终编译）
# ============================================================
zephyr_library_sources(
    tflm.cpp
)

# ============================================================
# 模型选择（TRD_TFLM select 对应模型）
# ============================================================
if(CONFIG_TFLM_MODEL_CHASSIS)
    zephyr_library_sources(
        models/chassis/chassis_model_data.cc
    )
endif()

# ============================================================
# 暴露公共 API 头文件路径给整个项目（只暴露这一个路径）
# ============================================================
zephyr_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

endif() # CONFIG_TFLM
