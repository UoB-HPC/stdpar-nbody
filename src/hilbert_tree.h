#pragma once
#include <algorithm>
#include <bit>
#include <execution>
#include <ranges>

#include "arguments.h"
#include "atomic_tree.h"
#include "saving.h"
#include "system.h"

using clock_timer = std::chrono::steady_clock;

/// Computes the bounding box of the grid.
template <typename T, dim_t N>
aabb<T, N> bounding_box(std::span<vec<T, N>> xs) {
  return std::transform_reduce(
   std::execution::par_unseq, xs.begin(), xs.end(), aabb<T, N>(from_points, xs[0]),
   [](auto a, auto b) { return merge(a, b); }, [](auto a) { return aabb<T, N>(from_points, a); });
}

/// Sorts bodies along the Hilbert curve.
template <typename T, dim_t N>
void hilbert_sort(System<T, N>& system, aabb<T, N> bbox) {
  if (N != 2 && N != 3) {
    std::cerr << "Hilbert sort is currently only implemented for 2D and 3D" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Create a Cartesian grid to bucket bodies into
  uint32_t hilbert_cells_per_dim = N == 2 ? 0xffffffff : (N == 3) ? 0x1fffff : 0;  // 21 bits set256;
  vec<T, N> grid_cell_size       = bbox.lengths() / (T)hilbert_cells_per_dim;

  // Compute the Hilbert index for each body in the Cartesian grid
  auto bids = system.body_indices();
  std::vector<uint32_t> hilbert_ids(system.size);
  std::for_each(std::execution::par_unseq, bids.begin(), bids.end(),
                [hids = hilbert_ids.data(), x = system.x.data(), mins = bbox.xmin, grid_cell_size](auto idx) {
                  // Bucket the body into a Cartesian grid cell:
                  vec<uint32_t, N> cell_idx = cast<uint32_t>((x[idx] - mins) / grid_cell_size);
                  // Compute the Hilber index of the cell and assign it to the body:
                  hids[idx] = hilbert(cell_idx);
                });

  // Sort the body ids according to the hilbert id
#if defined(__clang__) || defined(__NVCOMPILER)
  // clang & nvc++ struggle with zip_view
  // Workaround: copy everything to a vector of tuples (allocated only once), sort that, then copy things back
  // TODO: sort an array of keys and then apply a permutation in O(N) time and O(1) storage
  // (instead of O(N) time and O(N) storage).
  static std::vector<std::tuple<uint64_t, vec<T, 3>, T, vec<T, 3>, vec<T, 3>, vec<T, 3>>> tmp(system.size);
  std::for_each(std::execution::par_unseq, tmp.begin(), tmp.end(),
                [tmp = tmp.data(), hid = hilbert_ids.data(), x = system.x.data(), m = system.m.data(),
                 v = system.v.data(), a = system.a.data(), ao = system.ao.data()](auto& e) {
                  auto idx = &e - tmp;
                  e        = std::make_tuple(hid[idx], x[idx], m[idx], v[idx], a[idx], ao[idx]);
                });
  // sort by hilbert id
  std::sort(std::execution::par_unseq, tmp.begin(), tmp.end(),
            [](auto const & a, auto const & b) { return std::get<0>(a) < std::get<0>(b); });
  // copy back
  std::for_each(std::execution::par_unseq, tmp.begin(), tmp.end(),
                [tmp = tmp.data(), hid = hilbert_ids.data(), x = system.x.data(), m = system.m.data(),
                 v = system.v.data(), a = system.a.data(), ao = system.ao.data()](auto& e) {
                  auto idx = &e - tmp;
                  x[idx]   = std::get<1>(e);
                  m[idx]   = std::get<2>(e);
                  v[idx]   = std::get<3>(e);
                  a[idx]   = std::get<4>(e);
                  ao[idx]  = std::get<5>(e);
                });

#else
  auto r = std::views::zip(hilbert_ids, system.x, system.m, system.v, system.a, system.ao);
  std::sort(std::execution::par_unseq, r.begin(), r.end(),
            [hids = hilbert_ids.data()](auto a, auto b) { return std::get<0>(a) < std::get<0>(b); });
#endif
}

template <typename T, dim_t N>
struct bvh {
  uint32_t deep_level;
  aabb<T, N>* b;  //< Axis-Aligned Bounding-Boxes
  T* m;           //< Monopole Masses
  vec<T, N>* x;   //< Monopole Centers of Mass

  // The number of nodes at each level is "2^l" because each tree level is fully refined:
  static constexpr uint32_t nnodes_at_level(uint32_t l) { return ipow2(l); }

  // The total number of nodes is "sum over l in [0, nlevels) of 2^l" which is just "2^(nlevels+1)-1"
  static constexpr uint32_t nnodes_until_level(uint32_t l) { return ipow2(l) - 1; }

  // The parent of node `n` which is at level `l`:
  static constexpr uint64_t parent(uint64_t n, uint32_t l) {
    if (l == 0) return uint64_t(0);
    auto b = nnodes_until_level(l);
    auto o = n - b;
    return uint64_t(nnodes_until_level(l - 1) + (o / 2));
  }

  // The left child of node `n` which is at level `l`:
  static constexpr uint64_t left_child(uint64_t n, uint32_t l) {
    auto first = nnodes_until_level(l);
    auto count = nnodes_at_level(l);
    return (n - first) * 2 + first + count;
  }

  // Range of nodes at level `l`:
  static constexpr auto nodes(uint32_t l) {
    uint64_t first = nnodes_until_level(l);
    uint64_t count = nnodes_at_level(l);
    uint64_t last  = first + count;
    return std::views::iota((uint32_t)first, (uint32_t)last);
  }

  // Allocate the bvh:
  static bvh alloc(System<T, N> const & system) {
    // #leafs is the smallest power of two larger than #bodies:
    uint64_t nleafs = std::bit_ceil(system.size);

    // #levels is the number of trailing zeros of #leafs plus one.
    // We already have the leaf level built (its just the bodies),
    // so we do not include it here:
    uint32_t nlevels = std::countr_zero(nleafs);

    uint64_t nnodes = nnodes_until_level(nlevels);
    return bvh<T, N>{
     .deep_level = nlevels - 1, .b = new aabb<T, N>[nnodes], .m = new T[nnodes], .x = new vec<T, N>[nnodes]};
  }

  // Deallocate the bvh:
  static void dealloc(bvh t) {
    delete[] t.b;
    delete[] t.m;
    delete[] t.x;
  }

  // Build the bvh:
  void build_tree(System<T, N>& system) {
    uint64_t nbodies = system.size;
    // We build the deepest BVH level from the bodies:
    {
      auto ids   = nodes(deep_level);
      auto first = ids[0];
      std::for_each(std::execution::par_unseq, ids.begin(), ids.end(), [=, s = system.state(), *this](uint64_t i) {
        auto li = i - first;
        auto bl = (li * 2);
        auto br = bl + 1;
        if (bl >= nbodies) {
          m[i] = 0.;  // If the mass of a node is zero, that node is "dead"
          return;
        }

        if (br >= nbodies) {
          m[i] = s.m[bl];
          x[i] = s.x[bl];
          b[i] = aabb<T, N>(from_points, s.x[bl]);
        } else {
          m[i] = s.m[bl] + s.m[br];
          x[i] = s.m[bl] * s.x[bl] + s.m[br] * s.x[br];
          x[i] /= m[i];
          b[i] = aabb<T, N>(from_points, s.x[bl], s.x[br]);
        }
      });
    }

    // We then recursively build the remainin levels:
    for (int32_t l = deep_level - 1; l >= 0; --l) {
      auto ids   = nodes(l);
      auto first = ids[0];
      auto count = std::ranges::size(ids);
      std::for_each(std::execution::par_unseq, ids.begin(), ids.end(), [=, *this](uint64_t i) {
        auto li = i - first;
        auto bl = (li * 2) + first + count;
        auto br = bl + 1;

        auto ibl = m[bl] != 0.;
        auto ibr = m[br] != 0.;

        if (!ibl) {
          m[i] = 0.;
          return;
        }

        if (!ibr) {
          m[i] = m[bl];
          x[i] = x[bl];
          b[i] = b[bl];
        } else {
          m[i] = m[bl] + m[br];
          x[i] = m[bl] * x[bl] + m[br] * x[br];
          x[i] /= m[i];
          b[i] = merge(b[bl], b[br]);
        }
      });
    }
  }

  // Compute the force for each body
  void compute_force(System<T, N>& system, T theta) {
    uint64_t nbodies         = system.size;
    constexpr uint64_t empty = std::numeric_limits<uint64_t>::max();
    auto ids                 = system.body_indices();
    std::for_each(std::execution::par_unseq, ids.begin(), ids.end(), [=, s = system.state(), *this](uint64_t i) {
      auto xs = s.x[i];

      uint64_t tree_index = 0;
      auto a              = vec<T, N>::splat(0);
      uint32_t level      = 0;

      // stackless tree traversal
      bool came_forwards = true;
      while (tree_index != empty) {
        auto next_node = [&](uint64_t i) {
          if (i == 0) return empty;  // Back at the root? Done!
          // If left child, go to right child; otherwise go to parent.
          return uint64_t((i - 1) % 2 ? parent(i, level) : i + 1);
        };
        uint64_t next_node_index = next_node(tree_index);
        if (came_forwards) {  // child or sibling node
          vec<T, N> xj = x[tree_index];
          T mj         = m[tree_index];
          if (mj == 0.) {
            // dead node, nothing to do
          } else if ((dist(xs, xj) / b[tree_index].lengths() < theta).all()) {
            // below threshold
            a += mj * (xj - xs) / dist3(xs, xj);
          } else if (level == deep_level) {
            // force with other bodies (not with itself)
            auto f        = [&](auto idx) { a += s.m[idx] * (s.x[idx] - xs) / dist3(xs, s.x[idx]); };
            uint64_t bidx = 2 * (tree_index - nnodes_until_level(deep_level));
            if (bidx < nbodies && bidx != i) f(bidx);
            bidx = bidx + 1;
            if (bidx < nbodies && bidx != i) f(bidx);
          } else {
            // go to children
            next_node_index = left_child(tree_index, level);
            level++;
          }
        }

        // Relies on children allocated after their parent
        came_forwards = next_node_index > tree_index;
        level         = came_forwards ? level : level - 1;
        tree_index    = next_node_index;
      }

      s.a[i] = s.c * a;
    });
  }
};

template <typename T, dim_t N>
void run_hilbert_binary_tree(System<T, N>& system, Arguments arguments) {
  Saver<T, N> saver(arguments);
  saver.save_all(system);
  T theta = arguments.theta;

  // Allocate the bvh:
  auto tree = bvh<T, N>::alloc(system);

  for (size_t step = 0; step < arguments.steps; step++) {
    auto s0 = clock_timer::now();
    // Sort bodies along Hilbert curve:
    hilbert_sort(system, bounding_box(std::span{system.x}));

    auto s1 = clock_timer::now();
    // Build the bvh and monopoles:
    tree.build_tree(system);

    auto s2 = clock_timer::now();
    // Compute the force on each body using Barnes-Hut:
    tree.compute_force(system);

    auto s3 = clock_timer::now();
    // Apply acceleration
    system.accelerate_step();
    auto s4 = clock_timer::now();

    if (arguments.print_info) {
      using dur_t = std::chrono::duration<double, std::milli>;
      std::cout << std::format("Timings:\n- Sort Bodies {:.2f} ms\n- Build Tree {:.2f} ms\n- Calc force "
                               "{:.2f} ms\n- Calc acceleration {:.2f} ms",
                               dur_t(s1 - s0).count(), dur_t(s2 - s1).count(), dur_t(s3 - s2).count(),
                               dur_t(s4 - s3).count())
                << std::endl;
      std::cout << std::format("Total mass: {: .5f}\n", tree.m[0]);
    }
    saver.save_all(system);
  }

  // Deallocate the bvh:
  bvh<T, N>::dealloc(tree);
}
