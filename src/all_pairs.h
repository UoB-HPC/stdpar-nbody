#ifndef ALL_PAIRS_H
#define ALL_PAIRS_H

#include <execution>
#include <ranges>
#include "arguments.h"
#include "atomic.h"
#include "system.h"
#include "saving.h"
#include "counting_iterator.h"

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


template<typename T>
void run_all_pairs_collapsed_step(System<T>& system, Arguments arguments) {
    Saver<T> saver(arguments);
    saver.save_points(system);

    // all pairs algorithm time step
    for (size_t step = 0; step < arguments.steps; step++) {
        // force step
        auto it = counting_iterator<uint64_t>(0);
        std::for_each_n(
            std::execution::par_unseq,
            it, system.size * system.size,
            [s = system.state()] (auto p) {
                auto j = p / s.sz;
	        auto i = p % s.sz;
		if (i == j) {
		  // Reset to zero by subtracting old acceleration:
		  atomic_ref<T>{s.a[i][0]}.fetch_sub(s.ao[i][0], memory_order_relaxed);
		  atomic_ref<T>{s.a[i][1]}.fetch_sub(s.ao[i][1], memory_order_relaxed);
		  return;
		}

		auto pi = s.x[i];
		auto pj = s.x[j];
		auto a = s.c * s.m[j] * (pj - pi) / dist3(pi, pj);
		atomic_ref<T>{s.a[i][0]}.fetch_add(a[0], memory_order_relaxed);
		atomic_ref<T>{s.a[i][1]}.fetch_add(a[1], memory_order_relaxed);
        });
        // position update step
        system.accelerate_step();

        // save positions
	saver.save_points(system);
    }
}

#endif //ALL_PAIRS_H
