#pragma once

#include <ranges>

#include "arguments.h"
#include "atomic.h"
#include "counting_iterator.h"
#include "execution.h"
#include "saving.h"
#include "system.h"
#include "timer.h"
#include "format.h"

template <typename T, dim_t N>
void all_pairs_force(System<T, N>& system) {
  auto r = system.body_indices();
  std::for_each(par_unseq, r.begin(), r.end(), [s = system.state()](auto i) {
    auto ai = vec<T, N>::splat(0);
    auto pi = s.x[i];
    for (typename System<T, N>::index_t j = 0; j < s.sz; j++) {
      auto pj = s.x[j];
      ai += s.m[j] * (pj - pi) / dist3(pi, pj);
    }
    s.a[i] = s.c * ai;
  });
}

template <typename T, dim_t N>
void all_pairs_collapsed_force(System<T, N>& system) {

  auto it = counting_iterator<uint64_t>(0);
  std::for_each_n(par_unseq, it, system.size * system.size, [s = system.state()](auto p) {
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
    auto a  = s.c * s.m[j] * (pj - pi) / dist3(pi, pj);
    atomic_ref<T>{s.a[i][0]}.fetch_add(a[0], memory_order_relaxed);
    atomic_ref<T>{s.a[i][1]}.fetch_add(a[1], memory_order_relaxed);
  });
}

template <typename T, dim_t N, typename Force>
void run_all_pairs(System<T, N>& system, Arguments arguments, char const * name, Force&& f) {
  Saver<T, N> saver(arguments);
  saver.save_all(system);

  // Benchmarking output
  if (arguments.csv_detailed || arguments.csv_total) {
    if (arguments.print_state) abort();
    if (arguments.print_info) abort();
    if (arguments.save_pos) abort();
    if (arguments.save_energy) abort();
    std::cout << "algorithm,dim,precision,nsteps,nbodies,total [s]";
    if (arguments.csv_detailed) std::cout << ",force [s],accel [s]";
    std::cout << "\n";
  }

  auto dt_force = dur_t(0);
  auto dt_accel = dur_t(0);
  auto dt_total = dur_t(0);

  if (arguments.csv_detailed) {
    dt_total = time([&] {
      // all pairs algorithm time step
      for (size_t step = 0; step < arguments.steps; step++) {
        // force step
        dt_force += time([&] { f(system); });
        // position update step
        dt_accel += time([&] { system.accelerate_step(); });
        // save positions
        saver.save_all(system);
      }
    });
  } else {
    dt_total = time([&] {
      // all pairs algorithm time step
      for (size_t step = 0; step < arguments.steps; step++) {
        // force step
        f(system);
        // position update step
        system.accelerate_step();
      }
    });
  }

  if (arguments.csv_detailed || arguments.csv_total) {
    std::cout << std::format("{},{},{},{},{},{:.2f}", name, N, sizeof(T) * 8, arguments.steps, system.size,
                             dt_total.count());
    if (arguments.csv_detailed) std::cout << std::format(",{:.2f},{:.2f}", dt_force.count(), dt_accel.count());
    std::cout << "\n";
  }
}

template <typename T, dim_t N>
void run_all_pairs_step(System<T, N>& system, Arguments arguments) {
  run_all_pairs(system, arguments, "all-pairs", [](auto& system) { all_pairs_force(system); });
}

template <typename T, dim_t N>
void run_all_pairs_collapsed_step(System<T, N>& system, Arguments arguments) {
  run_all_pairs(system, arguments, "all-pairs-collapsed", [](auto& system) { all_pairs_collapsed_force(system); });
}
