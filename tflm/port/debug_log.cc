/* Copyright 2025 The TensorFlow Authors. All Rights Reserved.
   SPDX-License-Identifier: Apache-2.0 */

#include "tensorflow/lite/micro/debug_log.h"

#include <zephyr/kernel.h>
#include <cstdarg>
#include <cstdio>

extern "C" void DebugLog(const char* format, va_list args) {
#ifndef TF_LITE_STRIP_ERROR_STRINGS
  char buffer[256];
  vsnprintf(buffer, sizeof(buffer), format, args);
  printk("%s", buffer);
#endif
}

extern "C" int DebugVsnprintf(char* buffer, size_t buf_size,
                               const char* format, va_list vlist) {
  return vsnprintf(buffer, buf_size, format, vlist);
}
