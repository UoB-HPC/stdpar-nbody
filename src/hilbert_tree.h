#pragma once
#include <algorithm>
#include <bit>
#include <ranges>

#include "arguments.h"
#include "atomic_tree.h"
#include "counting_iterator.h"
#include "execution.h"
#include "format.h"
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
#if defined(__NVCOMPILER)
  // Workaround for nvc++: we can use Thrust zip_iterator, which predates zip_view, but provides the same functionality,
  // and works just fine:
  auto b = thrust::make_zip_iterator(hilbert_ids.begin(), system.x.begin(), system.m.begin(), system.v.begin(),
                                     system.a.begin(), system.ao.begin());
  auto e = thrust::make_zip_iterator(hilbert_ids.end(), system.x.end(), system.m.end(), system.v.end(), system.a.end(),
                                     system.ao.end());
  std::sort(par_unseq, b, e, [](auto a, auto b) { return thrust::get<0>(a) < thrust::get<0>(b); });
#elif defined(__clang__) || (__cplusplus < 202302L)
  // clang & nvc++ struggle with zip_view
  // Workaround for clang: sort pair of (hilbert index, original index), then copy things back
  // TODO: sort an array of keys and then apply a permutation in O(N) time and O(1) storage
  // (instead of O(N) time and O(N) storage).

  static std::vector<std::pair<uint64_t, std::size_t>> hilbert_index_map(system.size);
  std::for_each_n(par_unseq, counting_iterator<std::size_t>(0), system.size,
                  [hmap = hilbert_index_map.data(), hids = hilbert_ids.data()](std::size_t idx) {
                    hmap[idx] = std::make_pair(hids[idx], idx);
                  });
  // sort by hilbert id
  std::sort(par_unseq, hilbert_index_map.begin(), hilbert_index_map.end(),
            [](auto const & a, auto const & b) { return a.first < b.first; });

  // create temp copy of system so that we don't get race conditions when
  // rearranging values in the next step
  static std::vector<std::tuple<vec<T, N>, T, vec<T, N>, vec<T, N>, vec<T, N>>> tmp_system(system.size);
  std::for_each_n(par_unseq, counting_iterator<std::size_t>(0), system.size,
                  [tmp_sys = tmp_system.data(), x = system.x.data(), m = system.m.data(), v = system.v.data(),
                   a = system.a.data(), ao = system.ao.data()](std::size_t idx) {
                    tmp_sys[idx] = std::make_tuple(x[idx], m[idx], v[idx], a[idx], ao[idx]);
                  });

  // copy back
  std::for_each_n(par_unseq, counting_iterator<std::size_t>(0), system.size,
                  [tmp_sys = tmp_system.data(), hmap = hilbert_index_map.data(), x = system.x.data(),
                   m = system.m.data(), v = system.v.data(), a = system.a.data(), ao = system.ao.data()](auto idx) {
                    std::size_t original_index = hmap[idx].second;
                    auto e = tmp_sys[original_index];
                    x[idx] = std::get<0>(e);
                    m[idx] = std::get<1>(e);
                    v[idx] = std::get<2>(e);
                    a[idx] = std::get<3>(e);
                    ao[idx] = std::get<4>(e);
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

  // returns number of leaves contained in a node of the given level.
  static constexpr auto ncontained_leaves_at_level(level_t l, level_t nlevels) {
    return 1 << (nlevels - l);
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

  static constexpr bool can_approximate(vec<T, N> xs, vec<T, N> xj, aabb<T, N> b, T theta) {
    auto lengths = b.lengths();
    T max_length = lengths[0];
    for (dim_t i = 1; i < N; ++i)
      if (lengths[i] > max_length) max_length = lengths[i];

    return (max_length * max_length) / dist2(xs, xj) < theta * theta;
  }

  // Compute the force for each body
  void compute_force(System<T, N>& system, T theta) {
    node_t nbodies         = system.size;
    auto ids               = system.body_indices();
    std::for_each(par_unseq, ids.begin(), ids.end(), [=, s = system.state(), *this](node_t i) {
      auto xs = s.x[i];

      node_t tree_index = 0;
      auto a            = vec<T, N>::splat(0);
      level_t level     = 0;
      level_t leaf_level = last_level + 1;
      // this only refers to levels in the tree, so leaf level is not
      // included.
      level_t nlevels   = last_level + 1;

      // stackless tree traversal
      std::size_t num_covered_particles = 0;
      while (num_covered_particles < s.sz) {
        node_t next_node_index = 0;
        level_t next_level = level;

        auto force_ascend_right = [&] {
          next_node_index = parent(tree_index, level) + 1;
          next_level = level - 1;
        };

        auto ascend_right = [&] {
          // If left child, go to right child; otherwise go to right uncle.
          if(node_t((tree_index - 1) % 2)) {
            force_ascend_right();
          } else {
            next_node_index = tree_index + 1;
          }
        };
        // descend right away
        auto descend_directly = [&] {
          next_node_index = left_child(tree_index, level);
          next_level = level + 1;
        };

        if (level == leaf_level) {
          // force with other bodies (not with itself)
          node_t bidx = tree_index - nnodes_until_level(leaf_level);

          for(int k = 0; k < 2; ++k) {
            if (bidx < nbodies && bidx != i){
              vec<T, N> xj = s.x[bidx];
              T mj         = s.m[bidx];

              a += mj * (xj - xs) / dist3(xs, xj);
            }
            ++bidx;
          }
          num_covered_particles += 2;

          force_ascend_right();
        } else {

          vec<T, N> xj = x[tree_index];
          T mj         = m[tree_index];

          if (can_approximate(xs, xj, b[tree_index], theta)) {
            // below threshold
            a += mj * (xj - xs) / dist3(xs, xj);
            num_covered_particles += ncontained_leaves_at_level(level, nlevels);

            ascend_right();
          } else {
            // go to children
            descend_directly();
          }
        }

        level         = next_level;
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
  if (arguments.csv_total) {
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

  auto dt_force     = dur_t(0);
  auto dt_accel     = dur_t(0);
  auto dt_bbox      = dur_t(0);
  auto dt_sort      = dur_t(0);
  auto dt_monopoles = dur_t(0);
  auto dt_fapprox   = dur_t(0);
  auto dt_total     = dur_t(0);
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
          dt_monopoles += time([&] { tree.build_tree(system); });

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
    auto kernels = [&] {
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
    };
    for (size_t step = 0; step < arguments.warmup_steps; step++) kernels();
    dt_total = time([&] {
      for (size_t step = arguments.warmup_steps; step < arguments.steps; step++) kernels();
    });
    arguments.steps -= arguments.warmup_steps;
  }

  if (arguments.csv_detailed || arguments.csv_total) {
    std::cout << std::format("{},{},{},{},{},{:.2f}", "hilbert-tree", N, sizeof(T) * 8, arguments.steps, system.size,
                             dt_total.count());

    if (arguments.csv_detailed) {
      std::cout << std::format(",{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}", dt_force.count(), dt_accel.count(),
                               dt_bbox.count(), dt_sort.count(), dt_monopoles.count(), dt_fapprox.count());
    }
    std::cout << "\n";
  }

  // Deallocate the bvh:
  bvh<T, N>::dealloc(tree);
}
