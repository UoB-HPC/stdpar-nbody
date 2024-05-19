#ifndef ALL_PAIRS_H
#define ALL_PAIRS_H

#include <execution>
#include <ranges>

#include "arguments.h"
#include "system.h"
#include "saving.h"

template<typename T>
void run_all_pairs_step(System<T>& system, Arguments arguments) {
    // enable saving
    std::optional<Saver<T>> saver = {};
    if (arguments.save_output) {
        saver = Saver<T>(arguments);
        saver.value().save_points(system);
    }

    // all pairs algorithm time step
    for (auto step = 0; step < arguments.steps; step++) {
        // force step
        std::for_each_n(
            std::execution::par_unseq,
            std::begin(system.index), system.size,
            [
                p_xs=system.positions_x.data(), p_ys=system.positions_y.data(),
                constant=system.constant, masses=system.masses.data(), size=system.size,
                a_xs=system.accel_x.data(), a_ys=system.accel_y.data()
            ] (uint32_t i) {
                T accel_x_i = 0;
                T accel_y_i = 0;
                T position_x_i = p_xs[i];
                T position_y_i = p_ys[i];
                for (auto j = 0; j < size; j++) {
                    T mass_j = masses[j];
                    T position_x_j = p_xs[j];
                    T position_y_j = p_ys[j];
                    T cube_l2_norm = calc_cube_l2<T>(position_x_i, position_y_i, position_x_j, position_y_j);
                    accel_x_i += mass_j * (position_x_j - position_x_i) / cube_l2_norm;
                    accel_y_i += mass_j * (position_y_j - position_y_i) / cube_l2_norm;
                }
                a_xs[i] = constant * accel_x_i;
                a_ys[i] = constant * accel_y_i;
        });

        // position update step
        system.accelerate_step();

        // save positions
        if (arguments.save_output) {
            saver.value().save_points(system);
        }
    }

    // close saving file
    if (arguments.save_output) {
        saver.value().finish();
    }
}

#endif //ALL_PAIRS_H
