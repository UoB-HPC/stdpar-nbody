#pragma once
#include <algorithm>
#include <bit>
#include <ranges>

#include "arguments.h"
#include "atomic_tree.h"
#include "execution.h"
#include "saving.h"
#include "system.h"

using clock_timer = std::chrono::steady_clock;

/// Computes the bounding box of the grid.
template <typename T, dim_t N>
aabb<T, N> bounding_box(std::span<vec<T, N>> xs) {
  return std::transform_reduce(
   par_unseq, xs.begin(), xs.end(), aabb<T, N>(from_points, vec<T, N>::splat(0.)),
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
  static std::vector<uint64_t> hilbert_ids(system.size);
  std::for_each(par_unseq, bids.begin(), bids.end(),
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
  static std::vector<std::tuple<uint64_t, vec<T, N>, T, vec<T, N>, vec<T, N>, vec<T, N>>> tmp(system.size);
  std::for_each(par_unseq, tmp.begin(), tmp.end(),
                [tmp = tmp.data(), hid = hilbert_ids.data(), x = system.x.data(), m = system.m.data(),
                 v = system.v.data(), a = system.a.data(), ao = system.ao.data()](auto& e) {
                  auto idx = &e - tmp;
                  e        = std::make_tuple(hid[idx], x[idx], m[idx], v[idx], a[idx], ao[idx]);
                });
  // sort by hilbert id
  std::sort(par_unseq, tmp.begin(), tmp.end(),
            [](auto const & a, auto const & b) { return std::get<0>(a) < std::get<0>(b); });
  // copy back
  std::for_each(par_unseq, tmp.begin(), tmp.end(),
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
  std::sort(par_unseq, r.begin(), r.end(), [](auto a, auto b) { return std::get<0>(a) < std::get<0>(b); });
#endif
}

template <typename T, dim_t N>
struct bvh {
  using node_t  = uint32_t;
  using level_t = uint32_t;
  level_t last_level;  //< Last level in the tree, i.e., leaf-level - 1 (leaf level is just bodies).
  aabb<T, N>* b;       //< Axis-Aligned Bounding-Boxes
  T* m;                //< Monopole Masses
  vec<T, N>* x;        //< Monopole Centers of Mass

  // The number of nodes at each level is "2^l" because each tree level is fully refined:
  static constexpr level_t nnodes_at_level(level_t l) { return ipow2(l); }

  // The total number of nodes is "sum over l in [0, nlevels) of 2^l" which is just "2^(nlevels+1)-1"
  static constexpr level_t nnodes_until_level(level_t l) { return ipow2(l) - 1; }

  // The parent of node `n` which is at level `l`:
  static constexpr node_t parent(node_t n, level_t l) {
    if (l == 0) return node_t(0);
    auto b = nnodes_until_level(l);
    auto o = n - b;
    return node_t(nnodes_until_level(l - 1) + (o / 2));
  }

  // The left child of node `n` which is at level `l`:
  static constexpr node_t left_child(node_t n, level_t l) {
    auto first = nnodes_until_level(l);
    auto count = nnodes_at_level(l);
    return (n - first) * 2 + first + count;
  }

  // Range of nodes at level `l`:
  static constexpr auto nodes(level_t l) {
    node_t first = nnodes_until_level(l);
    node_t count = nnodes_at_level(l);
    node_t last  = first + count;
    return std::views::iota((uint32_t)first, (uint32_t)last);
  }

  // Allocate the bvh:
  static bvh alloc(System<T, N> const & system) {
    // #leafs is the smallest power of two larger than #bodies:
    node_t nleafs = std::bit_ceil(system.size);

    // #levels is the number of trailing zeros of #leafs plus one.
    // We already have the leaf level built (its just the bodies),
    // so we do not include it here:
    level_t nlevels = std::countr_zero(nleafs);

    node_t nnodes = nnodes_until_level(nlevels);
    return bvh<T, N>{
     .last_level = nlevels - 1, .b = new aabb<T, N>[nnodes], .m = new T[nnodes], .x = new vec<T, N>[nnodes]};
  }

  // Deallocate the bvh:
  static void dealloc(bvh t) {
    delete[] t.b;
    delete[] t.m;
    delete[] t.x;
  }

  // Build the bvh:
  void build_tree(System<T, N>& system) {
    node_t nbodies = system.size;
    // We build the deepest BVH level from the bodies:
    {
      auto ids   = nodes(last_level);
      auto first = ids[0];
      std::for_each(par_unseq, ids.begin(), ids.end(), [=, s = system.state(), *this](node_t i) {
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
    for (int32_t l = last_level - 1; l >= 0; --l) {
      auto ids   = nodes(l);
      auto first = ids[0];
      auto count = std::ranges::size(ids);
      std::for_each(par_unseq, ids.begin(), ids.end(), [=, *this](node_t i) {
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
    node_t nbodies         = system.size;
    constexpr node_t empty = std::numeric_limits<node_t>::max();
    auto ids               = system.body_indices();
    std::for_each(par_unseq, ids.begin(), ids.end(), [=, s = system.state(), *this](node_t i) {
      auto xs = s.x[i];

      node_t tree_index = 0;
      auto a            = vec<T, N>::splat(0);
      level_t level     = 0;

      // stackless tree traversal
      bool came_forwards = true;
      while (tree_index != empty) {
        auto next_node = [&](node_t i) {
          if (i == 0) return empty;  // Back at the root? Done!
          // If left child, go to right child; otherwise go to parent.
          return node_t((i - 1) % 2 ? parent(i, level) : i + 1);
        };
        node_t next_node_index = next_node(tree_index);
        if (came_forwards) {  // child or sibling node
          vec<T, N> xj = x[tree_index];
          T mj         = m[tree_index];
          if (mj == 0.) {
            // dead node, nothing to do
          } else if ((dist(xs, xj) / b[tree_index].lengths() < theta).all()) {
            // below threshold
            a += mj * (xj - xs) / dist3(xs, xj);
          } else if (level == last_level) {
            // force with other bodies (not with itself)
            auto f      = [&](auto idx) { a += s.m[idx] * (s.x[idx] - xs) / dist3(xs, s.x[idx]); };
            node_t bidx = 2 * (tree_index - nnodes_until_level(last_level));
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

  // Benchmarking output
  if (arguments.csv_detailed || arguments.csv_total) {
    if (arguments.print_state) abort();
    if (arguments.print_info) abort();
    if (arguments.save_pos) abort();
    if (arguments.save_energy) abort();
    std::cout << "algorithm,dim,precision,nsteps,nbodies,total [s]";
    if (arguments.csv_detailed) std::cout << ",force [s],accel [s],bbox [s],sort [s],multipoles [s],force approx [s]";
    std::cout << "\n";
  }

  // Allocate the bvh:
  auto tree = bvh<T, N>::alloc(system);

  auto dt_force      = dur_t(0);
  auto dt_accel      = dur_t(0);
  auto dt_bbox       = dur_t(0);
  auto dt_sort       = dur_t(0);
  auto dt_multipoles = dur_t(0);
  auto dt_fapprox    = dur_t(0);
  auto dt_total      = dur_t(0);
  if (arguments.csv_detailed) {
    dt_total = time([&] {
      for (size_t step = 0; step < arguments.steps; step++) {
        dt_force += time([&] {
          // Bounding box
          aabb<T, N> bbox;
          dt_bbox += time([&] { bbox = bounding_box(std::span{system.x}); });

          // Sort bodies along Hilbert curve:
          dt_sort += time([&] { hilbert_sort(system, bbox); });

          // Build the bvh and monopoles:
          dt_multipoles += time([&] { tree.build_tree(system); });

          // Compute the force on each body using Barnes-Hut:
          dt_fapprox += time([&] { tree.compute_force(system, theta); });
        });

        // Apply acceleration
        dt_accel += time([&] { system.accelerate_step(); });

        if (arguments.print_info) std::cout << std::format("Total mass: {: .5f}\n", tree.m[0]);
        saver.save_all(system);
      }
    });
  } else {
    dt_total = time([&] {
      for (size_t step = 0; step < arguments.steps; step++) {
        // Bounding box
        aabb<T, N> bbox = bounding_box(std::span{system.x});

        // Sort bodies along Hilbert curve:
        hilbert_sort(system, bbox);

        // Build the bvh and monopoles:
        tree.build_tree(system);

        // Compute the force on each body using Barnes-Hut:
        tree.compute_force(system, theta);

        // Apply acceleration
        system.accelerate_step();
      }
    });
  }

  if (arguments.csv_detailed || arguments.csv_total) {
    std::cout << std::format("{},{},{},{},{},{:.2f}", "hilbert-tree", N, sizeof(T) * 8, arguments.steps, system.size,
                             dt_total.count());

    if (arguments.csv_detailed) {
      std::cout << std::format(",{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}", dt_force.count(), dt_accel.count(),
                               dt_bbox.count(), dt_sort.count(), dt_multipoles.count(), dt_fapprox.count());
    }
    std::cout << "\n";
  }

  // Deallocate the bvh:
  bvh<T, N>::dealloc(tree);
}
