#pragma once
#include <ranges>

// Workaround for two nvc++ bugs:
// - Not supporting 128-bit integers yet on GPUs,
//   used by std::views::iota(64-bit int) overflow checking.
// - Link failure with OpenMP on CPU due to unknown issue yet.
#if defined(_NVHPC_STDPAR_GPU) || defined(_NVHPC_STDPAR_MULTICORE)
  #include <thrust/iterator/counting_iterator.h>
#endif

template <typename T>
constexpr auto counting_iterator(T v = T(0)) {
#if defined(_NVHPC_STDPAR_GPU) || defined(_NVHPC_STDPAR_MULTICORE)
  return thrust::counting_iterator<T>(T(v));
#else
  return std::views::iota(T(v)).begin();
#endif
}
