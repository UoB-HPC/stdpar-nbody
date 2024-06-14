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
void barnes_hut_step(System<T>& system, Arguments arguments, AtomicQuadTreeContainer<T, Index_t> tree) {
    auto start_timer = clock_timer::now();

    clear_tree(system, tree);
    compute_bounded_atomic_quad_tree(system, tree);
    build_atomic_tree(system, tree);
    auto built_tree_timer = clock_timer::now();

    calc_mass_atomic_tree(system, tree);
    auto calc_mass_timer = clock_timer::now();

    calc_force_atomic_tree(system, tree.to_const(), static_cast<T>(arguments.theta));
    auto force_timer = clock_timer::now();

    system.accelerate_step();
    auto accel_timer = clock_timer::now();

    if (arguments.print_info) {
        std::cout << std::format("Timings:\n- Build Tree {:.2f} ms\n- Calc mass {:.2f} ms\n- Calc force {:.2f} ms\n- Calc acceleration {:.2f} ms",
                                 std::chrono::duration<double, std::milli>(built_tree_timer - start_timer).count(),
                                 std::chrono::duration<double, std::milli>(calc_mass_timer - built_tree_timer).count(),
                                 std::chrono::duration<double, std::milli>(force_timer - calc_mass_timer).count(),
                                 std::chrono::duration<double, std::milli>(accel_timer - force_timer).count())
                  << std::endl;
        std::cout << std::format("Tree size: {}\n", tree.bump_allocator->load());
        std::cout << std::format("Total mass: {: .5f}\n", tree.total_masses[0]);
    }
}

template<typename T>
void barnes_hut_run(System<T>& system, Arguments arguments) {
    using Index_t = std::uint32_t;

    // enable saving
    std::optional<Saver<T>> saver = {};
    if (arguments.save_output) {
        saver = Saver<T>(arguments);
        saver.value().save_points(system);
    }

    // init tree structure
    auto vector_tree = AtomicQuadTree<T, Index_t>(system.max_tree_node_size);
    auto tree = vector_tree.get_container();
    if (arguments.print_info) {
        std::cout << "Tree init complete\n";
    }
    for (auto step = 0; step < arguments.steps; step++) {
        barnes_hut_step<T, Index_t>(system, arguments, tree);

        if (arguments.save_output) {
            saver.value().save_points(system);
        }
    }

    if (arguments.save_output) {
        saver.value().finish();
    }
}

#endif //BARNES_HUT_H
