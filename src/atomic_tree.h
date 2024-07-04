#pragma once

#include "atomic.h"
#include <cassert>
#include <vector>

#include "system.h"

template<typename T, dim_t N, typename Index = std::uint32_t>
class atomic_tree {
public:
    Index capacity;  //< Tree node capacity
    // Bounding box:
    T root_side_length;
    vec<T, N> root_x;
    // Tree connectivity:
    Index* first_child; //< First child of a node (all siblings are stored contiguously)
    Index* parent;      //< Parent of a sibling group (nodes that share the same parent)
    // Atomic offset to next free child_group (bump allocator):
    atomic<Index>* next_free_child_group; 

    // Node data:
    T* total_masses;
    vec<T, N>* centre_masses;

    // Helper latch per node for the last thread to proceed to the parent during the
    // tree computation of tree node centroids and masses:
    atomic<Index>* child_mass_complete;

    // Sentinel values used in "first_child" array to indicate whether:
    static constexpr Index empty = std::numeric_limits<Index>::max(); //< Empty tree leaf node
    static constexpr Index body = empty - 1;  //< Tree leaf node contains body
    static constexpr Index locked = body - 1; //< Tree leaf node is concurrently being split by another thread ("locked")

    // Allocates tree with capacity to hold `capacity` nodes:
    static atomic_tree alloc(size_t capacity) {
      atomic_tree qt;
      qt.capacity = capacity;
      qt.first_child = new Index[capacity];
      qt.parent = new Index[1 + capacity/4];
      qt.next_free_child_group = new atomic<Index>(1);

      qt.total_masses = new T[capacity];
      qt.centre_masses = new vec<T, N>[capacity];

      qt.child_mass_complete = new atomic<Index>[capacity];
      return qt;
    }

    // Deallocates the tree
    static void dealloc(atomic_tree* qt) {
      delete[] qt->first_child;
      delete[] qt->parent;
      delete[] qt->total_masses;
      delete[] qt->centre_masses;
      delete[] qt->child_mass_complete;
      delete qt->next_free_child_group;
    }

    // Computes which node to traverse next after `i` in a bottom to top (child to root) tree traversal:
    Index next_node(Index i) const noexcept {
      // Root node:
      if (i == 0) return empty;
      // Sibling group
      Index sg = (i - 1) / child_count<N>;
      // Child within sibling group
      Index cp = (i - 1) % child_count<N>;
      return cp == (child_count<N> - 1) ? parent[sg] : i + 1;
    }

    // Index of the "sibling group" (all nodes that share the same parent node) of `i`:
    Index sg(Index i) {
      return i == 0 ? 0 : (i - 1) / child_count<N>;
    }

    // Resets tree node `i`
    void clear(Index i) {
      if (i == 0) next_free_child_group->store(1, memory_order_relaxed);
      first_child[i] = empty;
      parent[sg(i)] = empty;
      total_masses[i] = T(0);
      centre_masses[i] = vec<T, N>::splat(0);
      child_mass_complete[i].store(0, memory_order_relaxed);
    }

    // Resets all nodes in the tree to ready it for the next iteration:
    void clear(System<T, N>& system, Index last_node) {
      auto r = system.body_indices();
      std::for_each_n(
        std::execution::par_unseq,
        r.begin(), last_node,
        [tree=*this] (auto tree_index) mutable { tree.clear(tree_index); });
    }

    // Compute bounds of the root node of the tree
    // (finds min/max xy-coord of all bodies in `system`):
    void compute_bounds(System<T, N>& system) {
      auto r = system.body_indices();
      auto [min_size, max_size] = std::transform_reduce(
        std::execution::par_unseq,
        r.begin(), r.end(),
        std::make_tuple<T, T>(0, 0),
        [] (auto lhs, auto rhs) -> std::tuple<T, T> {
          return {gmin(std::get<0>(lhs), std::get<0>(rhs)), gmax(std::get<1>(lhs), std::get<1>(rhs))};
        },
        [s=system.state()] (auto i) -> std::tuple<T, T> {
          return {min(s.x[i]), max(s.x[i])};
        }
      );

      // adjust boundary
      max_size += 1;
      min_size -= 1;

      T divide = (max_size + min_size) / static_cast<T>(2);

      // add root node to tree
      root_side_length = max_size - min_size;
      root_x = vec<T, N>::splat(divide);
      first_child[0] = empty;
    }

    // Inserts body with `mass` at position `pos`:
    void insert(T mass, vec<T, N> pos) {
      Index tree_index = 0;  // insert into root
      vec<T, N> divide = root_x;
      T side_length = root_side_length;

      while (true) {
        // find current node status
        atomic_ref<Index> fc{first_child[tree_index]};
        auto status = fc.load(memory_order_acquire);
        if (status != empty && status != body && status != locked) {
            // If the node has children
            T half_length = side_length / static_cast<T>(4);  // / 2 is for new quad length, then / 2 is for half length

            // calculate correct "hyperant" that pos falls into
            Index child_pos = 0;
            Index level = 1;
            for (dim_t i = 0; i < N; i++) {
                child_pos += level * (pos[i] > divide[i]);
                level *= 2;
                divide[i] += (2 * (pos[i] > divide[i]) - 1) * half_length;
            }
            tree_index = status + child_pos; // status is the first child of tree_index
            side_length /= static_cast<T>(2);
        } else if (status == empty && fc.compare_exchange_weak(status, locked, memory_order_acquire)) {
            // compare_exchange_weak is fine because if it fails spuriously, we'll retry again
            total_masses[tree_index] = mass;
            centre_masses[tree_index] = pos;
            fc.store(body, memory_order_release);
            break;
        } else if (status == body && fc.compare_exchange_weak(status, locked, memory_order_acquire)) {
            // compare_exchange_weak is fine because if it fails spuriously, we'll retry again

            // create children
            Index first_child_index = next_free_child_group->fetch_add(child_count<N>, memory_order_relaxed);
            parent[sg(first_child_index)] = tree_index;

            // evict body at current index and insert into children keeping node locked
            auto p_x = centre_masses[tree_index];
            Index child_pos = 0;
            Index level = 1;
            for (dim_t i = 0; i < N; i++) {
                child_pos += level * (p_x[i] > divide[i]);
                level *= 2;
            }
            Index evicted_index = first_child_index + child_pos;
            centre_masses[evicted_index] = p_x;
            total_masses[evicted_index] = total_masses[tree_index];
            first_child[evicted_index] = body;

            // release node and continue to try to insert body
            fc.store(first_child_index, memory_order_release);
        }
      }
    }

    // Inserts all bodies in system into the tree:
    void insert(System<T, N>& system) {
      auto r = system.body_indices();
      std::for_each(
        std::execution::par,
        r.begin(), r.end(),
        [s=system.state(),tree=*this] (Index i) mutable { tree.insert(s.m[i], s.x[i]); }
      );
    }

    // Computes tree node centroids and masses that depend on the body at the leaf node `i`:
    void compute_tree(Index tree_index) {
      // If this node is not a leaf node with a body, we are done:
      if (first_child[tree_index] != body)
	return;

      // Accumulate masses up to the root
      do {
        // move up to parent
        tree_index = parent[sg(tree_index)];
        if (tree_index == empty) break;  // if reached root

        // No thread will be arriving from siblings that are empty leaves,
        // so count those:
        Index local_leaf_count = 0;
        auto leaf_child_index = first_child[tree_index];
        for (dim_t i = 0; i < child_count<N>; ++i)
          local_leaf_count += static_cast<Index>(first_child[leaf_child_index + i] == empty);
        Index expected_count = child_count<N> - 1 - local_leaf_count;

        // Arrive at parent releasing previous masses accumulated, and acquiring masses accumulated by other threads:
        if (child_mass_complete[tree_index].fetch_add(1, memory_order_acq_rel) != expected_count)
            break;

        // pick up all child masses
        T m = 0;
        auto x = vec<T, N>::splat(0);
        for (dim_t i = 0; i < child_count<N>; i++) {
            auto child_index = first_child[tree_index] + i;
            auto child_total_mass = total_masses[child_index];
            m += child_total_mass;
            x += child_total_mass * centre_masses[child_index];
        }
        x /= m;

        total_masses[tree_index] = m;
        centre_masses[tree_index] = x;
      } while(true);
    }

    // Compute tree nodes centroids and masses for all bodies in system:
    void compute_tree(System<T, N>& system) {
      auto r = system.body_indices();
      std::for_each_n(
        std::execution::par,
        r.begin(), capacity,
        [tree=*this] (auto i) mutable { tree.compute_tree(i); }
      );
    }

    // Compute force at position `x` using cut-off parameter `theta`:
    vec<T, N> compute_force(vec<T, N> x, T const theta) const {
      Index tree_index = 0;
      T side_length = root_side_length;
      auto a = vec<T, N>::splat(0);

      // stackless tree traversal
      bool came_forwards = true;
      while (tree_index != empty) {
        Index next_node_index = next_node(tree_index);
        if (came_forwards) {  // child or sibling node
            vec<T, N> xj = centre_masses[tree_index];
            // check if below threshold
            auto fc = first_child[tree_index];
            if (fc == empty || fc == body || side_length / dist(x, xj) < theta) {
                T mj = total_masses[tree_index];
                a += mj * (xj - x)/ dist3(x, xj);
            } else {  // visit children
                next_node_index = first_child[tree_index];
                side_length /= static_cast<T>(2);
            }
        }
        // this works as children are always allocated after their parent
        came_forwards = next_node_index > tree_index;
        side_length *= came_forwards ? static_cast<T>(1) : static_cast<T>(2);
        tree_index = next_node_index;
      }

      return a;
    }

    // Compute forces of all bodies in `system` using cut-off parameter `theta`:
    void compute_force(System<T, N>& system, T const theta) {
      auto r = system.body_indices();
      std::for_each(
        std::execution::par_unseq,
        r.begin(), r.end(),
        [s=system.state(),theta, tree=*this] (Index i) {
            s.a[i] = s.c * tree.compute_force(s.x[i], theta);
        }
      );
    }
};

