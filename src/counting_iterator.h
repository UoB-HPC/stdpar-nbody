#pragma once
#include <limits>
#include <ranges>

// Workaround for two nvc++ bugs:
// - Not supporting 128-bit integers yet on GPUs,
//   used by std::views::iota(64-bit int) overflow checking.
// - Link failure with OpenMP on CPU due to unknown issue yet.
#if defined(_NVHPC_STDPAR_GPU) || defined(_NVHPC_STDPAR_MULTICORE)
  #include <thrust/iterator/counting_iterator.h>
#endif

// Workaround for compilers not supporting views::iota properly. Originally
// from BabelStream.
//
// A lightweight counting iterator which will be used by the STL algorithms
// NB: C++ <= 17 doesn't have this built-in, and it's only added later in ranges-v3 (C++2a) which this
// implementation doesn't target
template <typename N>
struct ranged {
  struct iterator {
    using difference_type   = ptrdiff_t;
    using value_type        = N;
    using pointer           = N const *;
    using reference         = N;
    using iterator_category = std::random_access_iterator_tag;

    N i_{};

    // XXX This is not part of the iterator spec, it gets picked up by oneDPL if enabled.
    // Without this, the DPL SYCL backend collects the iterator data on the host and copies to the device.
    // This type is unused for any nother STL impl.
    using is_passed_directly = std::true_type;

    reference operator*() const { return i_; }
    iterator &operator++() {
      ++i_;
      return *this;
    }
    iterator operator++(int) {
      iterator copy(*this);
      ++i_;
      return copy;
    }

    iterator &operator--() {
      --i_;
      return *this;
    }
    iterator operator--(int) {
      iterator copy(*this);
      --i_;
      return copy;
    }

    iterator &operator+=(difference_type by) {
      i_ += by;
      return *this;
    }
    iterator &operator-=(difference_type by) {
      i_ -= by;
      return *this;
    }

    value_type operator[](difference_type const &i) const { return i_ + i; }

    difference_type operator-(iterator const &it) const { return i_ - it.i_; }
    iterator operator+(difference_type const v) const { return iterator(i_ + v); }
    iterator operator-(difference_type const v) const { return iterator(i_ - v); }
    friend iterator operator+(difference_type const v, iterator const it) { return it + v; }
    friend iterator operator-(difference_type const v, iterator const it) { return it - v; }

    bool operator==(iterator const &other) const { return i_ == other.i_; }
    bool operator!=(iterator const &other) const { return !(*this == other); }
    bool operator<(iterator const &other) const { return i_ < other.i_; }
    bool operator<=(iterator const &other) const { return i_ <= other.i_; }
    bool operator>(iterator const &other) const { return i_ > other.i_; }
    bool operator>=(iterator const &other) const { return i_ >= other.i_; }

    iterator()                            = default;
    iterator(iterator const &)            = default;
    iterator(iterator &&)                 = default;
    iterator &operator=(iterator const &) = default;
    iterator &operator=(iterator &&)      = default;
    explicit iterator(N start) : i_(start) {}
  };

  [[nodiscard]] iterator begin() const { return begin_; }
  [[nodiscard]] iterator end() const { return end_; }
  ranged(N begin, N end) : begin_(begin), end_(end) {}
  iterator begin_;
  iterator end_;
};

static_assert(std::random_access_iterator<ranged<unsigned int>::iterator>);

template <typename T>
constexpr auto counting_iterator(T v = T(0)) {
#if defined(_NVHPC_STDPAR_GPU) || defined(_NVHPC_STDPAR_MULTICORE)
  return thrust::counting_iterator<T>(T(v));
#elif defined(__clang__)
  return ranged<T>{v, std::numeric_limits<T>::max()}.begin();
#else
  return std::views::iota(T(v)).begin();
#endif
}
