#pragma once

#ifdef FMT_FORMAT_WORKAROUND
#include <fmt/core.h>
namespace std {
  using namespace fmt;
}
#else
#include <format>
#endif

