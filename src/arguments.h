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
    std::size_t size = 1'00'0;
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

inline auto parse_args(std::vector<std::string>&& args) {
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
	} else if (args[arg_index] == "--precision") {
            arg_index += 1;
	    if (args[arg_index] == "float") {
	      arguments.single_precision = true;
	    } else if (args[arg_index] == "double") {
	      arguments.single_precision = false;
	    } else {
	      std::cerr << "Unknown precision: \"" << args[arg_index] << "\"." << std::endl;
	      std::cerr << "Options are: double, float (default)." << std::endl;
	      std::exit(EXIT_FAILURE);
	    }
        } else if (args[arg_index] == "--double") {
            arguments.single_precision = false;
	} else if (args[arg_index] == "--algorithm") {
            arg_index += 1;
	    if (args[arg_index] == "all-pairs") {
	      arguments.simulation_algo = SimulationAlgo::AllPairs;
	    } else if (args[arg_index] == "all-pairs-collapsed") {
	      arguments.simulation_algo = SimulationAlgo::AllPairsCollapsed;
	    } else if (args[arg_index] == "barnes-hut") {
	      arguments.simulation_algo = SimulationAlgo::BarnesHut;
	    } else {
	      std::cerr << "Unknown algorithm: \"" << args[arg_index] << "\"." << std::endl;
	      std::cerr << "Options are: all-pairs, all-pairs-collapsed, barnes-hut (default)." << std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	} else if (args[arg_index] == "--workload") {
            arg_index += 1;
	    if (args[arg_index] == "plummer") {
	      arguments.simulation_type = SimulationType::Plummer;
	    } else if (args[arg_index] == "galaxy") {
	      arguments.simulation_type = SimulationType::Galaxy;
	    } else if (args[arg_index] == "uniform") {
	      arguments.simulation_type = SimulationType::Uniform;
	    } else {
	      std::cerr << "Unknown workload: \"" << args[arg_index] << "\"." << std::endl;
	      std::cerr << "Options are: plummer, galaxy, uniform (default)." << std::endl;
	      std::exit(EXIT_FAILURE);
	    }
        } else if (args[arg_index] == "--print-state") {
            arguments.print_state = true;
        } else if (args[arg_index] == "--print-info") {
            arguments.print_info = true;
        } else if (args[arg_index] == "--save") {
	  arg_index += 1;
	  if (args[arg_index] == "pos") {
            arguments.save_pos = true;
	  } else if (args[arg_index] == "energy") {
            arguments.save_energy = true;
	  } else if (args[arg_index] == "all") {
	    arguments.save_pos = true;
	    arguments.save_energy = true;
	  } else if (args[arg_index] == "none") {
	    arguments.save_pos = false;
	    arguments.save_energy = false;
	  } else {
	    std::cerr << "Unknown save options: \"" << args[arg_index] << "\"." << std::endl;
	    std::cerr << "Options are: pos, energy, all, none (default)." << std::endl;
	    std::exit(EXIT_FAILURE);
	  }
        } else if (args[arg_index] == "--help" || args[arg_index] == "-h") {
            std::cout << ("Help:\n"
                          "-n size\t\tNumber of particles to simulate\n"
                          "-s steps\t\tNumber of steps to run simulation for\n"
                          "--theta t\t\tTheta threshold parameter to use in Barnes-Hut\n"
                          "--precision double|float(default)\t\tSelects floating-point precision\n"
                          "--algorithm all-pairs|all-pairs-collapsed|barnes-hut(default)<algo>\t\tSelects simulation algorithm\n"
                          "--workload plummer|galaxy|uniform(default)\t\tSelects workload\n"
                          "--print-state\t\tPrint the initial and final state of the simulation\n"
                          "--print-info\t\tPrint info every timestep\n"
                          "--save pos|energy|all|none(default) \t\tSelects what data to save every timestep\n"
                          "--help\t\tDisplay this help message and quit\n"
            );
            std::exit(EXIT_SUCCESS);
        } else {
            std::cout << std::format("Unknown argument: '{}'\n", args[arg_index]);
            std::exit(1);
        }
    }

    return arguments;
}

#endif //ARGUMENTS_H
