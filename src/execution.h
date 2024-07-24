#pragma once
#include <execution>

#if defined(SEQUENTIAL)
inline constexpr auto par       = std::execution::seq;
inline constexpr auto par_unseq = std::execution::seq;
#else
inline constexpr auto par       = std::execution::par;
inline constexpr auto par_unseq = std::execution::par_unseq;
#endif
