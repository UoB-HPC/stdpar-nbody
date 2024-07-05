#pragma once
// Eventually we'd replace this functionality with `std::simd`.
#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>

using dim_t = std::uint32_t;

template <dim_t N>
inline constexpr dim_t child_count = 2 * child_count<N - 1>;

template <>
inline constexpr dim_t child_count<1> = 2;

template <typename T, int N>
struct vec {
  // ignore alignment for odd dimensions
  alignas(N % 2 == 0 ? alignof(T) * N : alignof(T)) T data[N];

  constexpr T& operator[](int i) { return data[i]; }
  constexpr T operator[](int i) const { return data[i]; }

  static constexpr vec<T, N> splat(T v) {
    vec<T, N> o;
    for (int i = 0; i < N; ++i) o[i] = v;
    return o;
  }
  static constexpr vec<T, N> from(T x, T y) {
    vec<T, N> v;
    v[0] = x;
    v[1] = y;
    return v;
  }
};

template <typename T>
void write(vec<T, 2> v, T& x, T& y) {
  x = v[0];
  y = v[1];
}

////////////////////////////////////////////////////////////////////////////////
// Arithmetic ops:

template <typename T, int N>
constexpr bool operator==(vec<T, N> const & a, vec<T, N> const & b) {
  bool is_equal = true;
  for (int i = 0; i < N; ++i) is_equal = is_equal && (a[i] == b[i]);
  return is_equal;
}

template <typename T, int N>
constexpr vec<T, N>& operator+=(vec<T, N>& a, vec<T, N> b) {
  for (int i = 0; i < N; ++i) a[i] += b[i];
  return a;
}

template <typename T, int N>
constexpr vec<T, N> operator+(vec<T, N> a, vec<T, N> b) {
  a += b;
  return a;
}

template <typename T, int N>
constexpr vec<T, N>& operator-=(vec<T, N>& a, vec<T, N> b) {
  for (int i = 0; i < N; ++i) a[i] = a[i] - b[i];
  return a;
}

template <typename T, int N>
constexpr vec<T, N> operator-(vec<T, N> a, vec<T, N> b) {
  a -= b;
  return a;
}

template <typename T, int N>
constexpr vec<T, N>& operator*=(vec<T, N>& a, T s) {
  for (int i = 0; i < N; ++i) a[i] *= s;
  return a;
}

template <typename T, int N>
constexpr vec<T, N> operator*(T s, vec<T, N> b) {
  b *= s;
  return b;
}
template <typename T, int N>
constexpr vec<T, N> operator*(vec<T, N> a, T s) {
  a *= s;
  return a;
}

template <typename T, int N>
constexpr vec<T, N>& operator/=(vec<T, N>& a, T s) {
  for (int i = 0; i < N; ++i) a[i] /= s;
  return a;
}

template <typename T, int N>
constexpr vec<T, N> operator/(T s, vec<T, N> b) {
  b /= s;
  return b;
}
template <typename T, int N>
constexpr vec<T, N> operator/(vec<T, N> a, T s) {
  a /= s;
  return a;
}

////////////////////////////////////////////////////////////////////////////////
// Horizontal reductions

template <typename T>
T gmin(T a, T b) {
  if constexpr (std::is_floating_point_v<std::decay_t<T>>) return std::fmin(a, b);
  else return std::min(a, b);
}

template <typename T, int N>
constexpr T min(vec<T, N> v) {
  T r = v[0];
  for (int i = 1; i < N; ++i) r = gmin(r, v[i]);
  return r;
}

template <typename T>
T gmax(T a, T b) {
  if constexpr (std::is_floating_point_v<std::decay_t<T>>) return std::fmax(a, b);
  else return std::max(a, b);
}

template <typename T, int N>
constexpr T max(vec<T, N> v) {
  T r = v[0];
  for (int i = 1; i < N; ++i) r = gmax(r, v[i]);
  return r;
}

////////////////////////////////////////////////////////////////////////////////
// Vector math ops:

// Fused Multiply Add: s * a + b
template <typename T, int N>
constexpr vec<T, N> fma(T s, vec<T, N> a, vec<T, N> b) {
  for (int i = 0; i < N; ++i) a[i] = std::fma(s, a[i], b[i]);
  return a;
}

// Fused Multiply Add: a * b + c
template <typename T, int N>
constexpr vec<T, N> fma(vec<T, N> a, vec<T, N> b, vec<T, N> c) {
  for (int i = 0; i < N; ++i) a[i] = std::fma(a[i], b[i], c[i]);
  return a;
}

// L2-Norm^2
template <typename T, int N>
constexpr T l2norm2(vec<T, N> v) {
  T tmp = T(0.);
  for (int i = 0; i < N; ++i) tmp += v[i] * v[i];
  return tmp;
}

// L2-Norm
template <typename T, int N>
constexpr T l2norm(vec<T, N> v) {
  return std::sqrt(l2norm2(v));
}

// Distance^2
template <typename T, int N>
constexpr T dist2(vec<T, N> a, vec<T, N> b) {
  T tmp = T(0.);
  for (int i = 0; i < N; ++i) {
    T di = std::abs(a[i] - b[i]);
    tmp += di * di;
  }
  return tmp;
}

// Distance
template <typename T, int N>
constexpr T dist(vec<T, N> a, vec<T, N> b) {
  return std::sqrt(dist2(a, b)) + std::numeric_limits<T>::epsilon();
}

// Distance^3
template <typename T, int N>
constexpr T dist3(vec<T, N> a, vec<T, N> b) {
  return std::pow(dist2(a, b), T(3.) / T(2.)) + std::numeric_limits<T>::epsilon();
}
