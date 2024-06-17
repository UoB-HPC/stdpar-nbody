#ifndef BARNES_HUT_H
#define BARNES_HUT_H

#include <algorithm>
#include <execution>
#include <ranges>

#include "arguments.h"
#include "atomic_quad_tree.h"
#include "kernels.h"
#include "system.h"
#include "saving.h"

using clock_timer = std::chrono::steady_clock;

template<typename T, typename Index_t>
void barnes_hut_step(System<T>& system, Arguments arguments, AtomicQuadTree<T, Index_t> tree) {
    auto start_timer = clock_timer::now();

    clear_tree(system, tree);
    compute_bounded_atomic_quad_tree(system, tree);
    build_atomic_tree(system, tree);
    auto built_tree_timer = clock_timer::now();

    calc_mass_atomic_tree(system, tree);
    auto calc_mass_timer = clock_timer::now();

    calc_force_atomic_tree(system, tree, static_cast<T>(arguments.theta));
    auto force_timer = clock_timer::now();

    system.accelerate_step();
    auto accel_timer = clock_timer::now();

    if (arguments.print_info) {
        using dur_t = std::chrono::duration<double, std::milli>;
        std::cout << std::format("Timings:\n- Build Tree {:.2f} ms\n- Calc mass {:.2f} ms\n- Calc force {:.2f} ms\n- Calc acceleration {:.2f} ms",
                                 dur_t(built_tree_timer - start_timer).count(),
                                 dur_t(calc_mass_timer - built_tree_timer).count(),
                                 dur_t(force_timer - calc_mass_timer).count(),
                                 dur_t(accel_timer - force_timer).count())
                  << std::endl;
        std::cout << std::format("Tree size: {}\n", tree.bump_allocator->load());
        std::cout << std::format("Total mass: {: .5f}\n", tree.total_masses[0]);
    }
}

template<typename T>
void run_barnes_hut(System<T>& system, Arguments arguments) {
    using Index_t = std::uint32_t;

    Saver<T> saver(arguments);
    saver.save_points(system);

    // init tree structure
    auto tree = AtomicQuadTree<T, Index_t>::alloc(system.max_tree_node_size);
    if (arguments.print_info) {
        std::cout << "Tree init complete\n";
    }
    for (size_t step = 0; step < arguments.steps; step++) {
        barnes_hut_step<T, Index_t>(system, arguments, tree);
	saver.save_points(system);
    }
}

#endif //BARNES_HUT_H
