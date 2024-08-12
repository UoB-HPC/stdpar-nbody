#pragma once
#include "vec.h"

template <typename T, dim_t N>
struct monopole {
  vec<T, N + 1> data;

  monopole(T mass, vec<T, N> x) {
    for (dim_t i = 0; i < N; ++i) data[i] = x[i];
    data[N] = mass;
  }

  T mass() { return data[N]; }
  vec<T, N> x() {
    vec<T, N> x;
    for (dim_t i = 0; i < N; ++i) x[i] = m[i];
    return x;
  }
};
