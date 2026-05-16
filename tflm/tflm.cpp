/**
 * @file tflm.cpp
 * @author qingyu
 * @brief
 * @version 0.1
 * @date 2026-05-16
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "tflm.hpp"

#include <zephyr/kernel.h>

#include "tensorflow/lite/micro/micro_interpreter.h"
#include "tensorflow/lite/micro/micro_mutable_op_resolver.h"
#include "tensorflow/lite/micro/system_setup.h"
#include "tensorflow/lite/schema/schema_generated.h"

namespace tflm {

void init()
{
    tflite::InitializeTarget();
}

void print_info() {
    printk("TFLM initialized\n");
}

}  // namespace tflm
