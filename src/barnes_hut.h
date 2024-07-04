#pragma once
#include <algorithm>
#include <execution>
#include <ranges>

#include "arguments.h"
#include "atomic_tree.h"
#include "saving.h"
#include "system.h"

using clock_timer = std::chrono::steady_clock;

template <typename T, dim_t N, typename Index>
void barnes_hut_step(System<T, N>& system, Arguments arguments, atomic_tree<T, N, Index> tree, bool first) {
  auto start_timer = clock_timer::now();

  tree.clear(system, first ? tree.capacity : tree.next_free_child_group->load(memory_order_relaxed));
  tree.compute_bounds(system);
  tree.insert(system);
  auto built_tree_timer = clock_timer::now();

  tree.compute_tree(system);
  auto calc_mass_timer = clock_timer::now();

  tree.compute_force(system, static_cast<T>(arguments.theta));
  auto force_timer = clock_timer::now();

  system.accelerate_step();
  auto accel_timer = clock_timer::now();

  if (arguments.print_info) {
    using dur_t = std::chrono::duration<double, std::milli>;
    std::cout << std::format("Timings:\n- Build Tree {:.2f} ms\n- Calc mass {:.2f} ms\n- Calc force "
                             "{:.2f} ms\n- Calc acceleration {:.2f} ms",
                             dur_t(built_tree_timer - start_timer).count(),
                             dur_t(calc_mass_timer - built_tree_timer).count(),
                             dur_t(force_timer - calc_mass_timer).count(), dur_t(accel_timer - force_timer).count())
              << std::endl;
    std::cout << std::format("Tree size: {}\n", tree.next_free_child_group->load());
    std::cout << std::format("Total mass: {: .5f}\n", tree.total_masses[0]);
  }
}

template <typename T, dim_t N>
void run_barnes_hut(System<T, N>& system, Arguments arguments) {
  Saver<T, N> saver(arguments);
  saver.save_all(system);

  // init tree structure
  auto tree = atomic_tree<T, N>::alloc(system.max_tree_node_size);
  if (arguments.print_info) std::cout << "Tree init complete\n";
  for (size_t step = 0; step < arguments.steps; step++) {
    barnes_hut_step(system, arguments, tree, step == 0);
    saver.save_all(system);
  }
}
