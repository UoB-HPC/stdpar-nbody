#include <algorithm>
#include <chrono>
#include <iostream>

#ifndef DIM_SIZE
  #error Must specify spatial dimensions by compiling with -DDIM_SIZE=2 or -DDIM_SIZE=3 .
#endif

using clock_timer = std::chrono::steady_clock;

#include "all_pairs.h"
#include "arguments.h"
#include "barnes_hut.h"
#include "hilbert_tree.h"
#include "models.h"

template <typename T, dim_t N>
using sim_func_t = std::function<void(System<T, N>&, Arguments)>;

template <typename T, dim_t N>
void run_simulation(Arguments arguments, System<T, N>& system, sim_func_t<T, N> sim_algo) {
  if (arguments.print_state) {
    std::cout << "Starting state:" << std::endl;
    system.print();
  }
  std::cout << "Starting simulation" << std::endl;
  auto start = clock_timer::now();

  sim_algo(system, arguments);

  auto end = clock_timer::now();

  if (arguments.print_state) {
    std::cout << "Final state:" << std::endl;
    system.print();
  }
  std::cout << std::format("Done simulation\nTotal time: {:.2f} ms\n",
                           std::chrono::duration<double, std::milli>(end - start).count());
}

template <typename T, dim_t N>
auto run_precision(Arguments arguments) -> void {
  // arguments maybe modified by these functions
  auto system = [&arguments] {
    switch (arguments.simulation_type) {
      case SimulationType::Plummer: return build_plummer_model<T, N>(arguments);
      case SimulationType::Uniform: return build_uniform_model<T, N>(arguments);
      case SimulationType::Galaxy: return build_galaxy_model<T, N>(arguments);
      case SimulationType::Load: {
        auto system    = Saver<T, N>::load_system(arguments.load_input.value());
        arguments.size = system.size;
        return system;
      }
      default: throw std::runtime_error("Unknown simulation type");
    }
  }();
  switch (arguments.simulation_algo) {
    case SimulationAlgo::AllPairs: return run_simulation<T, N>(arguments, system, run_all_pairs_step<T, N>);
    case SimulationAlgo::AllPairsCollapsed:
      return run_simulation<T, N>(arguments, system, run_all_pairs_collapsed_step<T, N>);
    case SimulationAlgo::BarnesHut: return run_simulation<T, N>(arguments, system, run_barnes_hut<T, N>);
    case SimulationAlgo::HilbertTree: return run_simulation<T, N>(arguments, system, run_hilbert_binary_tree<T, N>);
  }
}

auto main(int argc, char* argv[]) -> int {
  auto arguments = parse_args(std::vector<std::string>(argv + 1, argv + argc));

  if (arguments.single_precision) run_precision<float, DIM_SIZE>(arguments);
  else run_precision<double, DIM_SIZE>(arguments);

  return EXIT_SUCCESS;
}
