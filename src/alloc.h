#pragma once
#include <cstdint>
#include <cstdlib>
#include <memory>

template <typename T>
void alloc(T*& ptr, std::size_t n) {
#ifdef MEM_ALIGNED
  ptr = (T*)std::aligned_alloc(/*2MiB*/ 2*1024*1024, sizeof(T) * n);
#else
  ptr = new T[n];
#endif  
}

template <typename T>
void dealloc(T* ptr, std::size_t n) {
  (void)n;
#ifdef MEM_ALIGNED
  std::free(ptr);
#else
  delete[] ptr;
#endif
}

template <typename T>
struct allocator {
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t ;
  using propagate_on_container_move_assignment  = std::true_type;
  
  T* allocate(std::size_t n, const void* hint = 0) {
    (void)hint;
    T* ptr;
    alloc(ptr, n);
    return ptr;
  }
  void deallocate(T* ptr, std::size_t n) {
    dealloc(ptr, n);
  }
  constexpr size_type max_size() const { return size_type(std::numeric_limits<unsigned int>::max() / sizeof(T)); }
  void construct(T* p, const T& value) { std::construct_at(p, value); }
  void destroy(T* p) { std::destroy_at(p); }
};

template <typename T, typename U>
constexpr bool operator==(allocator<T> const &, allocator<U> const&) noexcept { return true; }

template <typename T, typename U>
constexpr bool operator!=(allocator<T> const&, allocator<U> const&) noexcept { return false; }
