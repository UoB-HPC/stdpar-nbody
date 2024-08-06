#pragma once

#if defined(FMT_FORMAT_WORKAROUND) || (__cplusplus < 202302L)
  #define FMT_HEADER_ONLY
  #include <fmt/core.h>
namespace std {
using namespace fmt;
}
#else
  #include <format>
#endif
