#pragma once
#include <algorithm>
#include <ranges>

#include "arguments.h"
#include "atomic_tree.h"
#include "execution.h"
#include "saving.h"
#include "system.h"
#include "timer.h"
#include "format.h"

template <typename T, dim_t N>
void run_barnes_hut(System<T, N>& system, Arguments arguments) {
  Saver<T, N> saver(arguments);
  saver.save_all(system);

  // Benchmarking output
  if (arguments.csv_detailed || arguments.csv_total) {
    if (arguments.print_state) abort();
    if (arguments.print_info) abort();
    if (arguments.save_pos) abort();
    if (arguments.save_energy) abort();
    std::cout << "algorithm,dim,precision,nsteps,nbodies,total [s]";
    if (arguments.csv_detailed)
      std::cout << ",force [s],accel [s],clear [s],bbox [s],insert [s],multipoles [s],force approx [s]";
    std::cout << "\n";
  }

  // init tree structure
  auto tree = atomic_tree<T, N>::alloc(system.max_tree_node_size);
  if (arguments.print_info) std::cout << "Tree init complete\n";

  auto dt_force      = dur_t(0);
  auto dt_accel      = dur_t(0);
  auto dt_clear      = dur_t(0);
  auto dt_bbox       = dur_t(0);
  auto dt_insert     = dur_t(0);
  auto dt_monopoles = dur_t(0);
  auto dt_fapprox    = dur_t(0);
  auto dt_total      = dur_t(0);

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
          std::cout << std::format("Total mass: {: .5f}\n", tree.total_masses[0]);
        }
        saver.save_all(system);
      }
    });
  } else {
    dt_total = time([&] {
      for (size_t step = 0; step < arguments.steps; step++) {
        tree.clear(system, (step == 0) ? tree.capacity : tree.next_free_child_group->load(memory_order_relaxed));
        tree.compute_bounds(system);
        tree.insert(system);
        tree.compute_tree(system);
        tree.compute_force(system, static_cast<T>(arguments.theta));
        system.accelerate_step();
      }
    });
  }

  if (arguments.csv_detailed || arguments.csv_total) {
    std::cout << std::format("{},{},{},{},{},{:.2f}", "barnes-hut", N, sizeof(T) * 8, arguments.steps, system.size,
                             dt_total.count());

    if (arguments.csv_detailed) {
      std::cout << std::format(",{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}", dt_force.count(), dt_accel.count(),
                               dt_clear.count(), dt_bbox.count(), dt_insert.count(), dt_monopoles.count(),
                               dt_fapprox.count());
    }
    std::cout << "\n";
  }
}
