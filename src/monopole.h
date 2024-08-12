#pragma once
#include <type_traits>
#include <utility>

#include "vec.h"

template <typename T, dim_t N>
struct monopole {
  alignas(N % 2 == 0? alignof(vec<T, N>) : alignof(T) * (N+1)) vec<T, N> x_;
  T m_;

  monopole()                            = default;
  monopole(monopole const &)            = default;
  monopole(monopole &&)                 = default;
  monopole &operator=(monopole const &) = default;
  monopole &operator=(monopole &&)      = default;

  monopole(T m, vec<T, N> x) {
    for (dim_t i = 0; i < N; ++i) x_[i] = x[i];
    m_ = m;
  }
  T& mass() { return m_; }
  T const& mass() const { return m_; }
  vec<T, N>& x() { return x_; }
  vec<T, N> const& x() const { return x_; }
};

namespace std {

template <typename T, dim_t N>
struct tuple_size<::monopole<T, N>> : integral_constant<size_t, 2> {};

template <size_t Index, typename T, dim_t N>
struct tuple_element<Index, ::monopole<T, N>> : tuple_element<Index, tuple<T, vec<T, N>>> {};

}  // namespace std

template <std::size_t Index, typename T, dim_t N>
std::tuple_element_t<Index, monopole<T, N>> get(monopole<T, N> m) {
  if constexpr (Index == 0) return m.mass();
  if constexpr (Index == 1) return m.x();
}
