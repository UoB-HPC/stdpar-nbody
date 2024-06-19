#ifndef ALL_PAIRS_H
#define ALL_PAIRS_H

#include <execution>
#include <ranges>

#include "arguments.h"
#include "system.h"
#include "saving.h"

template<typename T>
void run_all_pairs_step(System<T>& system, Arguments arguments) {
    Saver<T> saver(arguments);
    saver.save_points(system);

    // all pairs algorithm time step
    for (size_t step = 0; step < arguments.steps; step++) {
        // force step
        auto r = system.body_indices();
        std::for_each(
            std::execution::par_unseq,
            r.begin(), r.end(),
            [s = system.state()] (auto i) {
                auto ai = vec<T, 2>::splat(0);
                auto pi = s.x[i];
                for (typename System<T>::index_t j = 0; j < s.sz; j++) {
                    auto pj = s.x[j];
                    ai += s.m[j] * (pj - pi) / dist3(pi, pj);
                }
                s.a[i] = s.c * ai;
        });

        // position update step
        system.accelerate_step();

        // save positions
	saver.save_points(system);
    }
}

#endif //ALL_PAIRS_H
