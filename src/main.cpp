#include <algorithm>
#include <chrono>
#include <iostream>

using clock_timer = std::chrono::steady_clock;

#include "arguments.h"
#include "all_pairs.h"
#include "barnes_hut.h"
#include "models.h"

template<typename T>
using sim_func_t = std::function<void(System<T>&, Arguments)>;

template<typename T>
void run_simulation(Arguments arguments, System<T>& system, sim_func_t<T> sim_algo) {
    if (arguments.print_state) {
        std::cout << "Starting state:" << std::endl;
        system.print();
    }
    std::cout << "Starting simulation\n" << std::endl;
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

template<typename T>
auto run_precision(Arguments arguments) -> void {
    auto system = [arguments] {
        switch (arguments.simulation_type) {
            case SimulationType::Plummer:
                return build_plummer_model<T>(arguments);
            case SimulationType::Solar:
                return build_solar_system<T>();
            case SimulationType::Galaxy:
                return build_galaxy_model<T>(arguments);
            default:
                throw std::runtime_error("Unknown simulation type");
        }
    }();
    auto sim_algo = arguments.barnes_hut ? barnes_hut_run<T> : run_all_pairs_step<T>;
    run_simulation<T>(arguments, system, sim_algo);
}


auto main(int argc, char* argv[]) -> int {
    auto arguments = parse_args(std::vector<std::string>(argv + 1, argv + argc));

    if (arguments.single_precision) {
        run_precision<float>(arguments);
    } else {
        run_precision<double>(arguments);
    }

    return EXIT_SUCCESS;
}
