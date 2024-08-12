#pragma once
#include <fstream>

#include "arguments.h"
#include "format.h"
#include "system.h"

template <typename T, dim_t N>
class Saver {
 public:
  explicit Saver(Arguments arguments)
   : is_save_pos(arguments.save_pos),
     is_save_energy(arguments.save_energy),
     simulation_size(arguments.size),
     timestep_count(arguments.steps) {
    if (is_save_pos) setup_position_file();
    if (is_save_energy) setup_energy_file();
  }

  auto save_all(System<T, N>& system) -> void {
    save_points(system);
    save_energy(system);
  }

  static auto load_system(std::string filename) -> System<T, N> {
    using int_t   = std::uint32_t;
    using float_t = float;
    auto in_file  = std::ifstream(filename, std::ios::binary);

    // read in state parameters
    int_t system_size;
    in_file.read(reinterpret_cast<char *>(&system_size), sizeof(system_size));
    int_t system_dim;
    in_file.read(reinterpret_cast<char *>(&system_dim), sizeof(system_dim));

    float_t time_step;
    in_file.read(reinterpret_cast<char *>(&time_step), sizeof(time_step));
    float_t constant;
    in_file.read(reinterpret_cast<char *>(&constant), sizeof(constant));

    if (system_dim != N) {
      throw std::runtime_error(
       std::format("This version is compiled with D={}, but the file provided is D={}", N, system_dim));
    }

    // read in state (mass, position, velocity)
    auto planet_size = 1 + 2 * system_dim;
    auto data_size   = system_size * planet_size;
    auto data        = std::vector<float_t>(data_size);
    in_file.read(reinterpret_cast<char *>(data.data()), data_size * sizeof(float_t));

    // create system
    auto system = System<T, N>(system_size, time_step, constant);

    // add points to system
    for (std::size_t i = 0; i < system_size; i++) {
      T mass   = data[i * planet_size];
      auto pos = vec<T, N>::splat(0);
      auto vel = vec<T, N>::splat(0);
      for (int_t j = 0; j < system_dim; j++) {
        pos[j] = data[i * planet_size + 1 + j];
        vel[j] = data[i * planet_size + 1 + N + j];
      }
      system.add_point(mass, pos, vel);
    }

    return system;
  }

  ~Saver() {
    if (is_save_pos) pos_out_file.close();
    if (is_save_energy) energy_out_file.close();
  }

 private:
  bool is_save_pos    = false;
  bool is_save_energy = false;
  std::uint32_t simulation_size;
  std::uint32_t timestep_count;
  std::uint32_t data_size = sizeof(T);

  std::ofstream pos_out_file;
  std::ofstream energy_out_file;

  auto setup_position_file() -> void {
    pos_out_file = std::ofstream("positions.bin", std::ios::out | std::ios::binary);

    // store position file parameters
    // number of bodies
    pos_out_file.write(reinterpret_cast<char const *>(&simulation_size), sizeof(simulation_size));
    // number of timesteps
    pos_out_file.write(reinterpret_cast<char const *>(&timestep_count), sizeof(timestep_count));
    // size (in bytes) of data type used to store positions
    pos_out_file.write(reinterpret_cast<char const *>(&data_size), sizeof(data_size));
    // dimension
    std::uint32_t dim = N;
    pos_out_file.write(reinterpret_cast<char const *>(&dim), sizeof(dim));
  }

  auto setup_energy_file() -> void {
    energy_out_file = std::ofstream("energy.bin", std::ios::out | std::ios::binary);

    // store position file parameters
    // number of timesteps
    energy_out_file.write(reinterpret_cast<char const *>(&timestep_count), sizeof(timestep_count));
    // size (in bytes) of data type used to store energy values
    energy_out_file.write(reinterpret_cast<char const *>(&data_size), sizeof(data_size));
  }

  auto save_points(System<T, N>& system) -> void {
    if (!is_save_pos) return;

    static vector<vec<T, N>> xs(simulation_size, vec<T, N>::splat(T(0)));
    auto it = counting_iterator<uint32_t>(0);
    std::for_each_n(par_unseq, it, simulation_size, [s = system.state(), xs = xs.data()](uint32_t i) {
      xs[i] = s.p[i].x();
    });

    pos_out_file.write(reinterpret_cast<char const *>(xs.data()), simulation_size * data_size * std::uint32_t(N));
  }

  auto save_energy(System<T, N> const &system) -> void {
    if (!is_save_energy) return;

    auto [kinetic, grav] = system.calc_energies();
    energy_out_file.write(reinterpret_cast<char const *>(&kinetic), sizeof(kinetic));
    energy_out_file.write(reinterpret_cast<char const *>(&grav), sizeof(grav));
  }
};
