#pragma once
#include <execution>

#if defined(SEQUENTIAL)
inline constexpr auto par       = std::execution::seq;
inline constexpr auto par_unseq = std::execution::seq;
#elif defined(UNSAFE_PAR_UNSEQ)
  #pragma message "This produces hangs and or invalid results on certain devices"
inline constexpr auto par       = std::execution::par_unseq;
inline constexpr auto par_unseq = std::execution::par_unseq;
#else
inline constexpr auto par       = std::execution::par;
inline constexpr auto par_unseq = std::execution::par_unseq;
#endif
