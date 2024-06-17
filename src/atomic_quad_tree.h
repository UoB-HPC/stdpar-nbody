#ifndef ATOMIC_QUAD_TREE_H
#define ATOMIC_QUAD_TREE_H

#include "atomic.h"
#include <cassert>
#include <vector>

#include "system.h"

// define sentinel value
template<typename Index_t>
constexpr Index_t is_leaf = static_cast<Index_t>(-1);

// phases that nodes in the tree can go through
enum class NodeStatus : uint32_t {
    Locked,
    EmptyLeaf,
    FullLeaf,
    NotLeaf,
};
static_assert(atomic<NodeStatus>::is_always_lock_free);

// Quad tree using Class of Vectors (CoV) (i.e. Structure of Arrays)
template<typename T, typename Index_t>
class AtomicQuadTree {
public:
    Index_t capacity;
    T root_side_length;
    vec<T, 2> root_x;
    Index_t* first_child;
    Index_t* next_nodes;
    Index_t* parent;
    atomic<NodeStatus>* node_status;
    // bump ptr used to keep track of allocated nodes
    atomic<Index_t>* bump_allocator;

    T* total_masses;
    vec<T, 2>* centre_masses;

    // used for mass calc
    atomic<Index_t>* child_mass_complete;  // stores number of children that have correct mass

    void clear(Index_t i) {
      if (i == 0) bump_allocator->store(1, memory_order_relaxed);
      first_child[i] = is_leaf<Index_t>;
      next_nodes[i] = is_leaf<Index_t>;
      parent[i] = is_leaf<Index_t>;
      node_status[i].store(NodeStatus::EmptyLeaf, memory_order_relaxed);
      total_masses[i] = T(0);
      centre_masses[i] = vec<T, 2>::splat(0);
      child_mass_complete[i].store(0, memory_order_relaxed);
    }

    static AtomicQuadTree alloc(size_t size) {
      AtomicQuadTree qt;
      qt.capacity = size;
      qt.first_child = new Index_t[size];
      qt.next_nodes = new Index_t[size];
      qt.parent = new Index_t[size];
      qt.node_status = new atomic<NodeStatus>[size];
      qt.bump_allocator = new atomic<Index_t>(1);

      qt.total_masses = new T[size];
      qt.centre_masses = new vec<T, 2>[size];

      qt.child_mass_complete = new atomic<Index_t>[size];
      return qt;
    }

    static void dealloc(AtomicQuadTree* qt) {
      delete qt->first_child[];
      delete qt->next_nodes[];
      delete qt->parent[];
      delete qt->node_status[];
      delete qt->total_masses[];
      delete qt->centre_masses[];
      delete qt->child_mass_complete[];
      delete qt->bump_allocator;
    }
};

#endif //ATOMIC_QUAD_TREE_H
