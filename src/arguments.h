#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <iostream>
#include <string>
#include <vector>

enum class SimulationType {
    Uniform,
    Plummer,
    Galaxy,
};

enum class SimulationAlgo {
    BarnesHut,
    AllPairs,
    AllPairsCollapsed,
};

struct Arguments {
    std::size_t size = 1'000;
    std::size_t steps = 1;
    bool single_precision = true;
    SimulationType simulation_type = SimulationType::Uniform;
    SimulationAlgo simulation_algo = SimulationAlgo::BarnesHut;
    bool print_state = false;
    bool print_info = false;
    double theta = 0.5;
    bool save_pos = false;
    bool save_energy = false;
};

auto parse_args(std::vector<std::string>&& args) {
    auto arguments = Arguments{};
    size_t arg_index = 0;

    for (; arg_index < args.size(); arg_index++) {
        if (args[arg_index] == "-n") {
            arg_index += 1;
            auto size = std::stoi(args[arg_index]);
            arguments.size = size;
        } else if (args[arg_index] == "-s") {
            arg_index += 1;
            auto steps = std::stoi(args[arg_index]);
            arguments.steps = steps;
        } else if (args[arg_index] == "--theta") {
            arg_index += 1;
            arguments.theta = std::stod(args[arg_index]);
        } else if (args[arg_index] == "--double") {
            arguments.single_precision = false;
        } else if (args[arg_index] == "--all-pairs") {
            arguments.simulation_algo = SimulationAlgo::AllPairs;
        } else if (args[arg_index] == "--all-pairs-collapsed") {
            arguments.simulation_algo = SimulationAlgo::AllPairsCollapsed;
        } else if (args[arg_index] == "--plummer") {
            arguments.simulation_type = SimulationType::Plummer;
        } else if (args[arg_index] == "--galaxy") {
            arguments.simulation_type = SimulationType::Galaxy;
        } else if (args[arg_index] == "--print-state") {
            arguments.print_state = true;
        } else if (args[arg_index] == "--print-info") {
            arguments.print_info = true;
        } else if (args[arg_index] == "--save-pos") {
            arguments.save_pos = true;
        } else if (args[arg_index] == "--save-energy") {
            arguments.save_energy = true;
        } else if (args[arg_index] == "--help" || args[arg_index] == "-h") {
            std::cout << ("Help:\n"
                          "-n size\t\tNumber of particles to simulate\n"
                          "-s steps\t\tNumber of steps to run simulation for\n"
                          "--theta t\t\tTheta threshold parameter to use in Barnes-Hut\n"
                          "--double\t\tUse double precision floating point (default is single precision)\n"
                          "--all-pairs\t\tUse all pairs simulation algorithm (default is barnes-hut)\n"
                          "--all-pairs-collapsed\t\tUse collapsed all pairs simulation algorithm (default is barnes-hut)\n"
                          "--plummer\t\tUse plummer model (D=3 only, default is uniform)\n"
                          "--galaxy\t\tUse galaxy colliding model (D=2 only, default is uniform)\n"
                          "--print-state\t\tPrint the initial and final state of the simulation\n"
                          "--print-info\t\tPrint info every timestep\n"
                          "--save-pos\t\tSave positions every timestep to positions.bin\n"
                          "--save-energy\t\tSave kinetic and gravitational energy every timestep to energy.bin\n"
                          "--help\t\tDisplay this help message and quit\n"
            );
            std::exit(0);
        } else {
            std::cout << std::format("Unknown argument: '{}'\n", args[arg_index]);
            std::exit(1);
        }
    }

    return arguments;
}

#endif //ARGUMENTS_H
