#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <iostream>
#include <string>
#include <vector>

enum class SimulationType {
    Solar,
    Plummer,
    Galaxy,
};

struct Arguments {
    std::size_t size = 1'000;
    std::size_t steps = 1;
    bool single_precision = true;
    bool barnes_hut = true;
    SimulationType simulation_type = SimulationType::Plummer;
    bool print_state = false;
    bool print_info = false;
    double theta = 0.5;
    bool save_output = false;
};

auto parse_args(std::vector<std::string>&& args) {
    auto arguments = Arguments{};
    auto arg_index = 0;

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
            arguments.barnes_hut = false;
        } else if (args[arg_index] == "--solar") {
            arguments.simulation_type = SimulationType::Solar;
        } else if (args[arg_index] == "--galaxy") {
            arguments.simulation_type = SimulationType::Galaxy;
        } else if (args[arg_index] == "--print-state") {
            arguments.print_state = true;
        } else if (args[arg_index] == "--print-info") {
            arguments.print_info = true;
        } else if (args[arg_index] == "--save") {
            arguments.save_output = true;
        } else if (args[arg_index] == "--help" || args[arg_index] == "-h") {
            std::cout << ("Help:\n"
                          "-n size\t\tNumber of particles to simulate\n"
                          "-s steps\t\tNumber of steps to run simulation for\n"
                          "--theta t\t\tTheta threshold parameter to use in Barnes-Hut\n"
                          "--double\t\tUse double precision floating point (default is single precision)\n"
                          "--all-pairs\t\tUse all pairs simulation (default is barnes-hut)\n"
                          "--solar\t\tUse solar system planet distribution (ignores size, steps are in days, uses double precision, default is plummer distribution)\n"
                          "--galaxy\t\tUse galaxy colliding model (default is plummer distribution)\n"
                          "--print-state\t\tPrint the initial and final state of the simulation\n"
                          "--print-info\t\tPrint info every timestep\n"
                          "--help\t\tDisplay this help message and quit\n"
            );
            std::exit(0);
        } else {
            std::cout << std::format("Unknown argument: '{}'\n", args[arg_index]);
            std::exit(1);
        }
    }

    // adjust arguments for solar simulation
    if (arguments.simulation_type == SimulationType::Solar) {
        arguments.size = 9;
        arguments.steps *= 24;
        // double precision is required to work
        arguments.single_precision = false;
    }

    return arguments;
}

#endif //ARGUMENTS_H
