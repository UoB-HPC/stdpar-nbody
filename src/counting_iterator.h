#ifndef COUNTING_ITERATOR_H
#define COUNTING_ITERATOR_H

#include <ranges>

// Workaround for nvc++ not supporting 128-bit integers yet
// used by std::views::iota(64-bit int) overflow checking:
#if defined(_NVHPC_STDPAR_GPU)
#include <thrust/iterator/counting_iterator.h>
#endif

template <typename T>
auto counting_iterator(T v = T(0)) {
  #if defined(_NVHPC_STDPAR_GPU)
      return thrust::counting_iterator<T>(T(v));
  #else
      return std::views::iota(T(v)).begin();
  #endif
}

#endif // COUNTING_ITERATOR_H
