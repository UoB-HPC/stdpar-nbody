#pragma once

#include <algorithm>
#include <cassert>
#include <ranges>
#include <vector>

#include "alloc.h"
#include "arguments.h"
#include "atomic.h"
#include "execution.h"
#include "format.h"
#include "saving.h"
#include "system.h"
#include "timer.h"

template <typename T, dim_t N, typename Index = std::uint32_t>
struct octree {
  Index capacity;  //< Tree node capacity
  // Bounding box:
  T root_side_length;
  vec<T, N> root_x;
  // Tree connectivity:
  mutable Index* first_child;  //< First child of a node (all siblings are stored contiguously)
  mutable Index* parent;       //< Parent of a sibling group (nodes that share the same parent)
  // Atomic offset to next free child_group (bump allocator):
  mutable atomic<Index>* next_free_child_group;

  // Node data:
  mutable monopole<T, N>* m;

  // Helper latch per node for the last thread to proceed to the parent during
  // the tree computation of tree node centroids and masses:
  mutable atomic<Index>* child_mass_complete;

  // Sentinel values used in "first_child" array to indicate whether:
  static constexpr Index empty  = std::numeric_limits<Index>::max();  //< Empty leaf node
  static constexpr Index body   = empty - 1;                          //< Tree leaf node contains body
  static constexpr Index locked = body - 1;                           //< Being modified by another thread

  // Allocates tree with capacity to hold `capacity` nodes:
  static octree alloc(size_t capacity) {
    octree qt;
    qt.capacity = capacity;
    ::alloc(qt.first_child, capacity);
    ::alloc(qt.parent, 1 + capacity / child_count<N>);
    ::alloc(qt.next_free_child_group, 1);
    ::alloc(qt.m, capacity);
    ::alloc(qt.child_mass_complete, capacity);
    return qt;
  }

  // Deallocates the tree
  static void dealloc(octree* qt) {
    dealloc(qt->first_child, qt->capacity);
    dealloc(qt->parent, 1 + qt->capacity / child_count<N>);
    dealloc(qt->next_free_child_group, 1);
    dealloc(qt->m, qt->capacity);
    dealloc(qt->child_mass_complete, qt->capacity);
  }

  // Computes which node to traverse next after `i` in a child to root traversal:
  Index next_node(Index i) const noexcept {
    // Root node:
    if (i == 0) return empty;
    // Sibling group
    Index sg = (i - 1) / child_count<N>;
    // Child within sibling group
    Index cp = (i - 1) % child_count<N>;
    return cp == (child_count<N> - 1) ? parent[sg] : i + 1;
  }

  // Index of "sibling group" (nodes sharing same parent) of `i`:
  static constexpr Index sg(Index i) { return i == 0 ? 0 : (i - 1) / child_count<N>; }

  // Resets tree node `i`
  void clear(Index i) const {
    if (i == 0) next_free_child_group->store(1, memory_order_relaxed);
    first_child[i] = empty;
    parent[sg(i)]  = empty;
    m[i]           = monopole(T(0), vec<T, N>::splat(0));
    child_mass_complete[i].store(0, memory_order_relaxed);
  }

  // Resets all nodes in the tree to ready it for the next iteration:
  void clear(System<T, N>& system, Index last_node) {
    auto r = system.body_indices();
    std::for_each_n(par_unseq, r.begin(), last_node, [tree = *this](auto tree_index) { tree.clear(tree_index); });
  }

  // Compute bounds of the root node of the tree
  // (finds min/max xy-coord of all bodies in `system`):
  void compute_bounds(System<T, N>& system) {
    auto r                    = system.body_indices();
    auto [min_size, max_size] = std::transform_reduce(
     par_unseq, r.begin(), r.end(), std::make_tuple<T, T>(0, 0),
     [](auto lhs, auto rhs) -> std::tuple<T, T> {
       return {gmin(std::get<0>(lhs), std::get<0>(rhs)), gmax(std::get<1>(lhs), std::get<1>(rhs))};
     },
     [s = system.state()](auto i) -> std::tuple<T, T> { return {min(s.x[i]), max(s.x[i])}; });

    // adjust boundary
    max_size += 1;
    min_size -= 1;

    T divide = (max_size + min_size) / static_cast<T>(2);

    // add root node to tree
    root_side_length = max_size - min_size;
    root_x           = vec<T, N>::splat(divide);
    first_child[0]   = empty;
  }

  // Inserts body with `mass` at position `pos`:
  void insert(T mass, vec<T, N> pos) const {
    Index tree_index = 0;  // insert into root
    vec<T, N> divide = root_x;
    T side_length    = root_side_length;

    while (true) {
      atomic_ref<Index> fc{first_child[tree_index]};
      auto status = fc.load(memory_order_acquire);
      if (status != empty && status != body && status != locked) {
        // If the node has children, traverse to the children:

        // Here /4: /2 is for new quad length, then / 2 is for half length
        T half_length = side_length / static_cast<T>(4);

        // calculate correct "hyperant" that pos falls into
        Index child_pos = 0;
        Index level     = 1;
        for (dim_t i = 0; i < N; i++) {
          child_pos += level * (pos[i] > divide[i]);
          level *= 2;
          divide[i] += (2 * (pos[i] > divide[i]) - 1) * half_length;
        }
        tree_index = status + child_pos;  // status is the first child of tree_index
        side_length /= static_cast<T>(2);
        continue;
      } else if (status == empty && fc.compare_exchange_weak(status, locked, memory_order_acquire)) {
        // If the node is empty and we locked it: insert body, unlock it, and done.
        // compare_exchange_weak suffices: if it fails spuriously, we'll retry again
        m[tree_index] = monopole(mass, pos);
        fc.store(body, memory_order_release);
        break;
      } else if (status == body && fc.compare_exchange_weak(status, locked, memory_order_acquire)) {
        // If node has body and we locked it: split node and unlock it.
        // We'll then continue traversing to try insert our current body.
        // compare_exchange_weak suffices: if it fails spuriously, we'll retry again

        // create children
        Index first_child_index       = next_free_child_group->fetch_add(child_count<N>, memory_order_relaxed);
        parent[sg(first_child_index)] = tree_index;

        // evict body at current index and insert into children keeping node locked
        auto [p_m, p_x] = m[tree_index];
        Index child_pos = 0;
        Index level     = 1;
        for (dim_t i = 0; i < N; i++) {
          child_pos += level * (p_x[i] > divide[i]);
          level *= 2;
        }
        Index evicted_index        = first_child_index + child_pos;
        m[evicted_index]           = monopole(p_m, p_x);
        first_child[evicted_index] = body;

        // release node and continue to try to insert body
        fc.store(first_child_index, memory_order_release);
        continue;
      }
      // Otherwise, we failed to lock the node because some other thread was
      // modifying it, so we try again.
    }
  }

  // Inserts all bodies in system into the tree:
  void insert(System<T, N>& system) {
    auto r = system.body_indices();
    std::for_each(par, r.begin(), r.end(),
                  [s = system.state(), tree = *this](Index i) { tree.insert(s.m[i], s.x[i]); });
  }

  // Computes tree node centroids and masses that depend on the body at the leaf node `i`:
  void compute_tree(Index tree_index) const {
    // If this node is not a leaf node with a body, we are done:
    if (first_child[tree_index] != body) return;

    // Accumulate masses up to the root
    do {
      // move up to parent
      tree_index = parent[sg(tree_index)];
      if (tree_index == empty) break;  // if reached root

      // No thread will be arriving from siblings that are empty leaves, so count those:
      Index local_leaf_count = 0;
      auto leaf_child_index  = first_child[tree_index];
      for (dim_t i = 0; i < child_count<N>; ++i)
        local_leaf_count += static_cast<Index>(first_child[leaf_child_index + i] == empty);
      Index expected_count = child_count<N> - 1 - local_leaf_count;

      // Arrive at parent releasing previous masses accumulated,
      // and acquiring masses accumulated by other threads:
      if (child_mass_complete[tree_index].fetch_add(1, memory_order_acq_rel) != expected_count) break;

      // pick up all child masses
      T m    = 0;
      auto x = vec<T, N>::splat(0);
      for (dim_t i = 0; i < child_count<N>; i++) {
        auto child_index        = first_child[tree_index] + i;
        auto [child_m, child_x] = this->m[child_index];
        m += child_m;
        x += child_m * child_x;
      }
      x /= m;

      this->m[tree_index] = monopole(m, x);
    } while (true);
  }

  // Compute tree nodes centroids and masses for all bodies in system:
  void compute_tree(System<T, N>& system) {
    auto r = system.body_indices();
    std::for_each_n(par, r.begin(), capacity, [tree = *this](auto i) { tree.compute_tree(i); });
  }

  // Compute force at position `x` using cut-off parameter `theta`:
  vec<T, N> compute_force(vec<T, N> x, T const theta) const {
    Index tree_index = 0;
    T side_length    = root_side_length;
    auto a           = vec<T, N>::splat(0);

    // stackless tree traversal
    bool came_forwards = true;
    while (tree_index != empty) {
      Index next_node_index = next_node(tree_index);
      if (came_forwards) {  // child or sibling node
        auto [mj, xj] = m[tree_index];
        // check if below threshold
        auto fc = first_child[tree_index];
        auto dx = dist(x, xj);
        if (fc == empty || fc == body || side_length / dx < theta) {
          a += mj * (xj - x) / (dx * dx * dx);
        } else {  // visit children
          next_node_index = first_child[tree_index];
          side_length /= static_cast<T>(2);
        }
      }
      // Relies on children allocated after their parent
      came_forwards = next_node_index > tree_index;
      side_length *= came_forwards ? static_cast<T>(1) : static_cast<T>(2);
      tree_index = next_node_index;
    }

    return a;
  }

  // Compute forces of all bodies in `system` using cut-off parameter `theta`:
  void compute_force(System<T, N>& system, T const theta) {
    auto r = system.body_indices();
    std::for_each(par_unseq, r.begin(), r.end(), [s = system.state(), theta, tree = *this](Index i) {
      s.a[i] = s.c * tree.compute_force(s.x[i], theta);
    });
  }
};

template <typename T, dim_t N>
void run_octree(System<T, N>& system, Arguments arguments) {
  Saver<T, N> saver(arguments);
  saver.save_all(system);

  // Benchmarking output
  if (arguments.csv_total) {
    if (arguments.print_state) abort();
    if (arguments.print_info) abort();
    if (arguments.save_pos) abort();
    if (arguments.save_energy) abort();
  }
  if (arguments.csv_total || arguments.csv_detailed) {
    std::cout << "algorithm,dim,precision,nsteps,nbodies,total [s]";
    if (arguments.csv_detailed)
      std::cout << ",force [s],accel [s],clear [s],bbox [s],insert [s],multipoles [s],force approx [s]";
    std::cout << "\n";
  }

  // init tree structure
  auto tree = octree<T, N>::alloc(system.max_tree_node_size);
  if (arguments.print_info) std::cout << "Tree init complete\n";

  auto dt_force     = dur_t(0);
  auto dt_accel     = dur_t(0);
  auto dt_clear     = dur_t(0);
  auto dt_bbox      = dur_t(0);
  auto dt_insert    = dur_t(0);
  auto dt_monopoles = dur_t(0);
  auto dt_fapprox   = dur_t(0);
  auto dt_total     = dur_t(0);

  if (arguments.csv_detailed) {
    dt_total = time([&] {
      for (size_t step = 0; step < arguments.steps; step++) {
        dt_force += time([&] {
          dt_clear += time([&] {
            tree.clear(system, (step == 0) ? tree.capacity : tree.next_free_child_group->load(memory_order_relaxed));
          });
          dt_bbox += time([&] { tree.compute_bounds(system); });
          dt_insert += time([&] { tree.insert(system); });
          dt_monopoles += time([&] { tree.compute_tree(system); });
          dt_fapprox += time([&] { tree.compute_force(system, static_cast<T>(arguments.theta)); });
        });

        dt_accel += time([&] { system.accelerate_step(); });

        if (arguments.print_info) {
          std::cout << std::format("Tree size: {}\n", tree.next_free_child_group->load());
          std::cout << std::format("Total mass: {: .5f}\n", tree.m[0].mass());
        }
        saver.save_all(system);
      }
    });
  } else {
    auto kernels = [&](size_t step) {
      tree.clear(system, (step == 0) ? tree.capacity : tree.next_free_child_group->load(memory_order_relaxed));
      tree.compute_bounds(system);
      tree.insert(system);
      tree.compute_tree(system);
      tree.compute_force(system, static_cast<T>(arguments.theta));
      system.accelerate_step();
    };
    for (size_t step = 0; step < arguments.warmup_steps; step++) kernels(step);
    dt_total = time([&] {
      for (size_t step = arguments.warmup_steps; step < arguments.steps; step++) kernels(step);
    });
    arguments.steps -= arguments.warmup_steps;
  }

  if (arguments.csv_detailed || arguments.csv_total) {
    std::cout << std::format("{},{},{},{},{},{:.2f}", "octree", N, sizeof(T) * 8, arguments.steps, system.size,
                             dt_total.count());

    if (arguments.csv_detailed) {
      std::cout << std::format(",{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}", dt_force.count(), dt_accel.count(),
                               dt_clear.count(), dt_bbox.count(), dt_insert.count(), dt_monopoles.count(),
                               dt_fapprox.count());
    }
    std::cout << "\n";
  }
}
