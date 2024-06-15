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

// container class for lambdas
template<typename T, typename Index_t>
class ConstAtomicQuadTreeContainer {
public:
    T const root_side_length;
    // T const root_x;
    // T const root_y;
    T const * const total_masses;
    T const * const centre_masses_x;
    T const * const centre_masses_y;
    Index_t const * const first_child;
    Index_t const * const next_nodes;
    // atomic<NodeStatus> const * const node_status;

    // atomic<int32_t> const * const leaf_count;
    // atomic<Index_t> const * const child_mass_complete;
    Index_t const * const parent;

    // only pointer field not a vector
    // atomic<Index_t> const * const bump_allocator;
};

// container class for lambdas
template<typename T, typename Index_t>
class AtomicQuadTreeContainer {
public:
    T root_side_length;
    T root_x;
    T root_y;
    T * const total_masses;
    T * const centre_masses_x;
    T * const centre_masses_y;
    Index_t * const first_child;
    Index_t * const next_nodes;
    atomic<NodeStatus> * const node_status;

    atomic<int32_t> * const leaf_count;
    atomic<Index_t> * const child_mass_complete;
    Index_t * const parent;

    // only a pointer field not a vector
    atomic<Index_t> * const bump_allocator;

    auto to_const() {
        return ConstAtomicQuadTreeContainer<T, Index_t>{
                root_side_length, // root_x, root_y,
                total_masses, centre_masses_x, centre_masses_y,
                first_child, next_nodes, // node_status,
                //leaf_count, child_mass_complete,
                parent,
                //bump_allocator
        };
    }
};

// Quad tree using Class of Vectors (CoV) (i.e. Structure of Arrays)
template<typename T, typename Index_t>
class AtomicQuadTree {
public:
    T root_side_length;
    T root_x;
    T root_y;
    std::vector<T> total_masses;
    std::vector<T> centre_masses_x;
    std::vector<T> centre_masses_y;
    std::vector<Index_t> first_child;
    std::vector<Index_t> next_nodes;
    std::vector<atomic<NodeStatus>> node_status;

    // used for mass calc
    std::vector<atomic<int32_t>> leaf_count;  // stores total number of sub leaves a node has
    std::vector<atomic<Index_t>> child_mass_complete;  // stores number of children that have correct mass
    std::vector<Index_t> parent;

    // bump ptr used to keep track of allocated nodes
    std::unique_ptr<atomic<Index_t>> bump_allocator = std::make_unique<atomic<Index_t>>(1);

    explicit AtomicQuadTree(size_t size)
        : total_masses(size, 0), centre_masses_x(size, 0), centre_masses_y(size, 0),
          first_child(size, is_leaf<Index_t>), next_nodes(size, is_leaf<Index_t>), node_status(size),
          leaf_count(size), child_mass_complete(size), parent(size, is_leaf<Index_t>){

    }

    auto get_container() {
        return AtomicQuadTreeContainer<T, Index_t>{
            root_side_length, root_x, root_y,
            total_masses.data(), centre_masses_x.data(), centre_masses_y.data(),
            first_child.data(), next_nodes.data(), node_status.data(),
            leaf_count.data(), child_mass_complete.data(), parent.data(),
            bump_allocator.get()
        };
    }
};

#endif //ATOMIC_QUAD_TREE_H
