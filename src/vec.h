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

template <typename T, dim_t N>
struct vec {
  // ignore alignment for odd dimensions
  alignas(N % 2 == 0 ? alignof(T) * N : alignof(T)) T data[N];

  constexpr T& operator[](int i) { return data[i]; }
  constexpr T operator[](int i) const { return data[i]; }

  static constexpr vec<T, N> splat(T v) {
    vec<T, N> o;
    for (dim_t i = 0; i < N; ++i) o[i] = v;
    return o;
  }
  static constexpr vec<T, N> from(T x, T y) {
    vec<T, N> v;
    v[0] = x;
    v[1] = y;
    return v;
  }

  bool all()
    requires(std::is_same_v<T, bool>)
  {
    bool r = true;
    for (dim_t i = 0; i < N; ++i) r &= (*this)[i];
    return r;
  }
};

template <typename T>
void write(vec<T, 2> v, T& x, T& y) {
  x = v[0];
  y = v[1];
}

////////////////////////////////////////////////////////////////////////////////
// Arithmetic ops:

template <typename T, dim_t N>
constexpr vec<T, N>& operator+=(vec<T, N>& a, vec<T, N> b) {
  for (dim_t i = 0; i < N; ++i) a[i] += b[i];
  return a;
}

template <typename T, dim_t N>
constexpr vec<T, N> operator+(vec<T, N> a, vec<T, N> b) {
  a += b;
  return a;
}

template <typename T, dim_t N>
constexpr vec<T, N>& operator-=(vec<T, N>& a, vec<T, N> b) {
  for (dim_t i = 0; i < N; ++i) a[i] = a[i] - b[i];
  return a;
}

template <typename T, dim_t N>
constexpr vec<T, N> operator-(vec<T, N> a, vec<T, N> b) {
  a -= b;
  return a;
}

template <typename T, dim_t N>
constexpr vec<T, N>& operator*=(vec<T, N>& a, T s) {
  for (dim_t i = 0; i < N; ++i) a[i] *= s;
  return a;
}

template <typename T, dim_t N>
constexpr vec<T, N> operator*(T s, vec<T, N> b) {
  b *= s;
  return b;
}
template <typename T, dim_t N>
constexpr vec<T, N> operator*(vec<T, N> a, T s) {
  a *= s;
  return a;
}

template <typename T, dim_t N>
constexpr vec<T, N>& operator/=(vec<T, N>& a, T s) {
  for (dim_t i = 0; i < N; ++i) a[i] /= s;
  return a;
}

template <typename T, dim_t N>
constexpr vec<T, N> operator/(T s, vec<T, N> b) {
  b /= s;
  return b;
}
template <typename T, dim_t N>
constexpr vec<T, N> operator/(vec<T, N> a, T s) {
  a /= s;
  return a;
}
template <typename T, dim_t N>
constexpr vec<T, N> operator/(vec<T, N> a, vec<T, N> b) {
  for (dim_t i = 0; i < N; ++i) a[i] /= b[i];
  return a;
}

////////////////////////////////////////////////////////////////////////////////
// Relational ops

template <typename T, dim_t N>
constexpr bool operator==(vec<T, N> const & a, vec<T, N> const & b) {
  bool is_equal = true;
  for (dim_t i = 0; i < N; ++i) is_equal = is_equal && (a[i] == b[i]);
  return is_equal;
}

template <typename T, dim_t N>
constexpr vec<bool, N> operator<(vec<T, N> a, vec<T, N> b) {
  vec<bool, N> v;
  for (dim_t i = 0; i < N; ++i) v[i] = a[i] < b[i];
  return v;
}

template <typename T, dim_t N>
constexpr vec<bool, N> operator<(vec<T, N> a, T s) {
  vec<bool, N> v;
  for (dim_t i = 0; i < N; ++i) v[i] = a[i] < s;
  return v;
}

////////////////////////////////////////////////////////////////////////////////
// Horizontal reductions

template <typename T>
T gmin(T a, T b) {
  if constexpr (std::is_floating_point_v<std::decay_t<T>>) return std::fmin(a, b);
  else return std::min(a, b);
}

template <typename T, dim_t N>
constexpr T min(vec<T, N> v) {
  T r = v[0];
  for (dim_t i = 1; i < N; ++i) r = gmin(r, v[i]);
  return r;
}

template <typename T>
T gmax(T a, T b) {
  if constexpr (std::is_floating_point_v<std::decay_t<T>>) return std::fmax(a, b);
  else return std::max(a, b);
}

template <typename T, dim_t N>
constexpr T max(vec<T, N> v) {
  T r = v[0];
  for (dim_t i = 1; i < N; ++i) r = gmax(r, v[i]);
  return r;
}

////////////////////////////////////////////////////////////////////////////////
// Vector math ops:

// Fused Multiply Add: s * a + b
template <typename T, dim_t N>
constexpr vec<T, N> fma(T s, vec<T, N> a, vec<T, N> b) {
  for (dim_t i = 0; i < N; ++i) a[i] = std::fma(s, a[i], b[i]);
  return a;
}

// Fused Multiply Add: a * b + c
template <typename T, dim_t N>
constexpr vec<T, N> fma(vec<T, N> a, vec<T, N> b, vec<T, N> c) {
  for (dim_t i = 0; i < N; ++i) a[i] = std::fma(a[i], b[i], c[i]);
  return a;
}

template <typename T, dim_t N>
constexpr vec<T, N> max(vec<T, N> a, vec<T, N> b) {
  for (dim_t i = 0; i < N; ++i) a[i] = gmax(a[i], b[i]);
  return a;
}

template <typename T, dim_t N>
constexpr vec<T, N> min(vec<T, N> a, vec<T, N> b) {
  for (dim_t i = 0; i < N; ++i) a[i] = gmin(a[i], b[i]);
  return a;
}

template <typename T, dim_t N>
constexpr vec<T, N> abs(vec<T, N> a) {
  for (dim_t i = 0; i < N; ++i) a[i] = std::abs(a[i]);
  return a;
}

template <typename T, dim_t N>
constexpr T sum(vec<T, N> v) {
  T tmp = T(0.);
  for (dim_t i = 0; i < N; ++i) tmp += v[i];
  return tmp;
}

template <typename T, dim_t N>
constexpr T avg(vec<T, N> v) {
  return sum(v) * (T(1) / T(N));
}

// L2-Norm^2
template <typename T, dim_t N>
constexpr T l2norm2(vec<T, N> v) {
  T tmp = T(0.);
  for (dim_t i = 0; i < N; ++i) tmp += v[i] * v[i];
  return tmp;
}

// L2-Norm
template <typename T, dim_t N>
constexpr T l2norm(vec<T, N> v) {
  return std::sqrt(l2norm2(v));
}

// Distance^2
template <typename T, dim_t N>
constexpr T dist2(vec<T, N> a, vec<T, N> b) {
  T tmp = T(0.);
  for (dim_t i = 0; i < N; ++i) {
    T di = std::abs(a[i] - b[i]);
    tmp += di * di;
  }
  return tmp;
}

// Distance
template <typename T, dim_t N>
constexpr T dist(vec<T, N> a, vec<T, N> b) {
  return std::sqrt(dist2(a, b)) + std::numeric_limits<T>::epsilon();
}

// Distance^3
template <typename T, dim_t N>
constexpr T dist3(vec<T, N> a, vec<T, N> b) {
  return std::pow(dist2(a, b), T(3.) / T(2.)) + std::numeric_limits<T>::epsilon();
}

////////////////////////////////////////////////////////////////////////////////
// Conversion ops
template <typename U, typename T, dim_t N>
vec<U, N> cast(vec<T, N> v) {
  vec<U, N> r;
  for (dim_t i = 0; i < N; ++i) r[i] = static_cast<U>(v[i]);
  return r;
}

////////////////////////////////////////////////////////////////////////////////
// Hilbert Curve ops

template <dim_t N>
constexpr uint64_t interleave_bits(vec<uint32_t, N> x) {
  if constexpr (N == 2) {
    auto bit_split = [](uint64_t x) {
      x = (x | x << 16) & uint64_t(0xffff0000ffff);
      x = (x | x << 8) & uint64_t(0xff00ff00ff00ff);
      x = (x | x << 4) & uint64_t(0xf0f0f0f0f0f0f0f);
      x = (x | x << 2) & uint64_t(0x3333333333333333);
      x = (x | x << 1) & uint64_t(0x5555555555555555);
      return x;
    };
    return bit_split(x[1]) | (bit_split(x[0]) << 1);
  } else if constexpr (N == 3) {
    auto bit_split = [](uint64_t x) {
      x &= uint64_t(0x1fffff);
      x = (x | x << 32) & uint64_t(0x1f00000000ffff);
      x = (x | x << 16) & uint64_t(0x1f0000ff0000ff);
      x = (x | x << 8) & uint64_t(0x100f00f00f00f00f);
      x = (x | x << 4) & uint64_t(0x10c30c30c30c30c3);
      x = (x | x << 2) & uint64_t(0x1249249249249249);
      return x;
    };
    return bit_split(x[2]) | (bit_split(x[1]) << 1) | (bit_split(x[0]) << 2);
  } else {
    printf("unimplemented interleave_bits!\n");
    abort();
  }
}

/// Computes the Hilbert index of coordinate 'x'.
/// See:
///   John Skilling, Programming the Hilbert curve,
///   AIP Conference Proceedings 707, 381 (2004), doi: 10.1063/1.1751381
template <dim_t N>
uint64_t hilbert(vec<uint32_t, N> x) {
  if constexpr (N == 2) {
    constexpr int32_t n = 2, bits = 32;
    constexpr uint32_t M = 1U << (bits - 1);

    for (uint32_t Q = M; Q > 1; Q >>= 1) {  // Inverse undo
      uint32_t P = Q - 1;
      for (uint32_t i = 0; i < n; ++i) {
        if ((x[i] & Q) != 0) {
          x[0] ^= P;  // invert
        } else {
          uint32_t t = (x[0] ^ x[i]) & P;
          x[0] ^= t;
          x[i] ^= t;
        }
      }  // exchange
    }

    // Gray encode
    for (uint32_t i = 1; i < n; ++i) x[i] ^= x[i - 1];
    uint32_t t = 0;
    for (uint32_t Q = M; Q > 1; Q >>= 1)
      if ((x[n - 1] & Q) != 0) t ^= Q - 1;
    for (uint32_t i = 0; i < n; ++i) x[i] ^= t;
    
    return interleave_bits(x);
  
  } else if constexpr (N == 3) {
    constexpr int32_t n = 2, bits = 21;
    constexpr uint32_t M = 1U << (bits - 1);
    // Inverse undo
    for (uint32_t Q = M; Q > 1; Q >>= 1) {
      uint32_t P = Q - 1;
      for (uint32_t i = 0; i < n; ++i) {
        if ((x[i] & Q) != 0) {
          x[0] ^= P;  // invert
        } else {
          uint32_t t = (x[0] ^ x[i]) & P;
          x[0] ^= t;
          x[i] ^= t;
        }
      }
    }  // exchange

    // Gray encode
    for (uint32_t i = 1; i < n; ++i) x[i] ^= x[i - 1];
    uint32_t t = 0;
    for (uint32_t Q = M; Q > 1; Q >>= 1)
      if ((x[n - 1] & Q) != 0) t ^= Q - 1;
    for (uint32_t i = 0; i < n; ++i) x[i] ^= t;

    return interleave_bits(x);
  } else {
    printf("unimplemented interleave_bits!\n");
    abort();
  }
}

////////////////////////////////////////////////////////////////////////////////
// Output

template <typename Os, typename T, dim_t N>
Os& operator<<(Os& os, vec<T, N> const & v) {
  os << "(";
  for (dim_t i = 0; i < N; ++i) {
    os << v[i];
    if (i != N - 1) os << ",";
  }
  os << ")";
  return os;
}

////////////////////////////////////////////////////////////////////////////////
// Axis-Aligned Bounding Box

struct from_bounds_t {};
inline constexpr from_bounds_t from_bounds;

struct from_points_t {};
inline constexpr from_points_t from_points;

template <typename T, dim_t N>
struct aabb {
  vec<T, N> xmin = vec<T, N>::splat(T(0));
  vec<T, N> xmax = vec<T, N>::splat(T(0));

  aabb() = default;
  aabb(from_bounds_t, vec<T, N> xmin, vec<T, N> xmax) : xmin(xmin), xmax(xmax) {}
  aabb(from_points_t, vec<T, N> p0) {
    constexpr auto tol = vec<T, N>::splat(std::numeric_limits<T>::epsilon() * 10.);
    xmin               = p0 - tol;
    xmax               = p0 + tol;
  }
  aabb(from_points_t, vec<T, N> p0, vec<T, N> p1) {
    constexpr auto tol = vec<T, N>::splat(std::numeric_limits<T>::epsilon() * 10.);
    xmin               = min(p0, p1) - tol;
    xmax               = max(p0, p1) + tol;
  }

  vec<T, N> lengths() { return xmax - xmin; }
};

template <typename T, dim_t N>
aabb<T, N> merge(aabb<T, N> a, aabb<T, N> b) {
  return aabb(from_bounds, min(a.xmin, b.xmin), max(a.xmax, b.xmax));
}

template <typename Os, typename T, dim_t N>
Os& operator<<(Os& os, aabb<T, N> const & b) {
  os << "{min:" << b.xmin << ",max:" << b.xmax << "}";
  return os;
}

////////////////////////////////////////////////////////////////////////////////
// Other math

/// Integer power of two:
template <typename T>
constexpr T ipow2(T e)
  requires(std::is_integral_v<T>)
{
  return T(1) << e;
}
