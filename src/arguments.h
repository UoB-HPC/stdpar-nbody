#pragma once
#include <iostream>
#include <optional>
#include <string>
#include <vector>

enum class SimulationType {
  Uniform,
  Plummer,
  Galaxy,
  Load,
};

enum class SimulationAlgo {
  AllPairs,
  AllPairsCollapsed,
  BarnesHut,
  HilbertTree,
};

struct Arguments {
  std::size_t size                      = 1'000;
  std::size_t steps                     = 1;
  bool single_precision                 = true;
  SimulationType simulation_type        = SimulationType::Uniform;
  SimulationAlgo simulation_algo        = SimulationAlgo::BarnesHut;
  bool print_state                      = false;
  bool print_info                       = false;
  double theta                          = 0.5;
  bool save_pos                         = false;
  bool save_energy                      = false;
  std::optional<std::string> load_input = {};
};

inline auto parse_args(std::vector<std::string>&& args) {
  using namespace std;
  auto arguments   = Arguments{};
  size_t arg_index = 0;

  for (; arg_index < args.size(); arg_index++) {
    if (args[arg_index] == "-n") {
      arg_index += 1;
      auto size      = stoi(args[arg_index]);
      arguments.size = size;
    } else if (args[arg_index] == "-s") {
      arg_index += 1;
      auto steps      = stoi(args[arg_index]);
      arguments.steps = steps;
    } else if (args[arg_index] == "--theta") {
      arg_index += 1;
      arguments.theta = stod(args[arg_index]);
    } else if (args[arg_index] == "--precision") {
      arg_index += 1;
      if (args[arg_index] == "float") {
        arguments.single_precision = true;
      } else if (args[arg_index] == "double") {
        arguments.single_precision = false;
      } else {
        cerr << "Unknown precision: \"" << args[arg_index] << "\"." << endl;
        cerr << "Options are: double, float (default)." << endl;
        exit(EXIT_FAILURE);
      }
    } else if (args[arg_index] == "--algorithm") {
      arg_index += 1;
      if (args[arg_index] == "all-pairs") {
        arguments.simulation_algo = SimulationAlgo::AllPairs;
      } else if (args[arg_index] == "all-pairs-collapsed") {
        arguments.simulation_algo = SimulationAlgo::AllPairsCollapsed;
      } else if (args[arg_index] == "barnes-hut") {
        arguments.simulation_algo = SimulationAlgo::BarnesHut;
      } else if (args[arg_index] == "hilbert-tree") {
        arguments.simulation_algo = SimulationAlgo::HilbertTree;
      } else {
        cerr << "Unknown algorithm: \"" << args[arg_index] << "\"." << endl;
        cerr << "Options are: all-pairs, all-pairs-collapsed, barnes-hut (default)." << endl;
        exit(EXIT_FAILURE);
      }
    } else if (args[arg_index] == "--workload") {
      arg_index += 1;
      if (args[arg_index] == "plummer") {
        arguments.simulation_type = SimulationType::Plummer;
      } else if (args[arg_index] == "galaxy") {
        arguments.simulation_type = SimulationType::Galaxy;
      } else if (args[arg_index] == "uniform") {
        arguments.simulation_type = SimulationType::Uniform;
      } else if (args[arg_index] == "load") {
        arg_index += 1;
        arguments.load_input      = {args[arg_index]};
        arguments.simulation_type = SimulationType::Load;
      } else {
        cerr << "Unknown workload: \"" << args[arg_index] << "\"." << endl;
        cerr << "Options are: plummer, galaxy, uniform (default)." << endl;
        exit(EXIT_FAILURE);
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
        arguments.save_pos    = true;
        arguments.save_energy = true;
      } else if (args[arg_index] == "none") {
        arguments.save_pos    = false;
        arguments.save_energy = false;
      } else {
        cerr << "Unknown save options: \"" << args[arg_index] << "\"." << endl;
        cerr << "Options are: pos, energy, all, none (default)." << endl;
        exit(EXIT_FAILURE);
      }
    } else if (args[arg_index] == "--help" || args[arg_index] == "-h") {
      cout << ("Help:\n"
               "-n size\t\tNumber of particles to simulate\n"
               "-s steps\t\tNumber of steps to run simulation for\n"
               "--theta t\t\tTheta threshold parameter to use in Barnes-Hut\n"
               "--precision double|float(default)\t\tSelects floating-point "
               "precision\n"
               "--algorithm "
               "all-pairs|all-pairs-collapsed|barnes-hut(default)<algo>"
               "\t\tSelects simulation algorithm\n"
               "--workload plummer|galaxy|uniform(default)|load "
               "<file.bin>\t\tSelects workload\n"
               "--print-state\t\tPrint the initial and final state of the "
               "simulation\n"
               "--print-info\t\tPrint info every timestep\n"
               "--save pos|energy|all|none(default) \t\tSelects what data to "
               "save every timestep\n"
               "--help\t\tDisplay this help message and quit\n");
      exit(EXIT_SUCCESS);
    } else {
      cout << format("Unknown argument: '{}'\n", args[arg_index]);
      exit(1);
    }
  }

  return arguments;
}
