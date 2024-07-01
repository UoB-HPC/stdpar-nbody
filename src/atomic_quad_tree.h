#ifndef ATOMIC_QUAD_TREE_H
#define ATOMIC_QUAD_TREE_H

#include "atomic.h"
#include <cassert>
#include <vector>

#include "system.h"


template<typename T, typename Index_t, dim_t N>
class AtomicQuadTree {
public:
    Index_t capacity;
    T root_side_length;
    vec<T, N> root_x;
    Index_t* first_child;
    Index_t* parent;
    // bump ptr used to keep track of allocated nodes
    atomic<Index_t>* bump_allocator;

    T* total_masses;
    vec<T, N>* centre_masses;

    // used for mass calc
    atomic<Index_t>* child_mass_complete;  // stores number of children that have correct mass

    static constexpr Index_t empty = std::numeric_limits<Index_t>::max();
    static constexpr Index_t body = empty - 1;
    static constexpr Index_t locked = body - 1;

    Index_t next_node(Index_t i) const noexcept {
      // Root node:
      if (i == 0) return empty;
      // Sibling group
      Index_t sg = (i - 1) / child_count<N>;
      // Child within sibling group
      Index_t cp = (i - 1) % child_count<N>;
      return cp == (child_count<N> - 1) ? parent[sg] : i + 1;
    }

    // Index of the "sibling group" of i:
    Index_t sg(Index_t i) {
      return i == 0 ? 0 : (i - 1) / child_count<N>;
    }

    void clear(Index_t i) {
      if (i == 0) bump_allocator->store(1, memory_order_relaxed);
      first_child[i] = empty;
      parent[sg(i)] = empty;
      total_masses[i] = T(0);
      centre_masses[i] = vec<T, N>::splat(0);
      child_mass_complete[i].store(0, memory_order_relaxed);
    }

    static AtomicQuadTree alloc(size_t size) {
      AtomicQuadTree qt;
      qt.capacity = size;
      qt.first_child = new Index_t[size];
      qt.parent = new Index_t[1 + size/4];
      qt.bump_allocator = new atomic<Index_t>(1);

      qt.total_masses = new T[size];
      qt.centre_masses = new vec<T, N>[size];

      qt.child_mass_complete = new atomic<Index_t>[size];
      return qt;
    }

    static void dealloc(AtomicQuadTree* qt) {
      delete[] qt->first_child;
      delete[] qt->parent;
      delete[] qt->total_masses;
      delete[] qt->centre_masses;
      delete[] qt->child_mass_complete;
      delete qt->bump_allocator;
    }
};

#endif //ATOMIC_QUAD_TREE_H
