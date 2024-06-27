#ifndef SAVING_H
#define SAVING_H

#include <fstream>

#include "arguments.h"
#include "system.h"

template<typename T>
class Saver {
public:
    explicit Saver(Arguments arguments)
    : is_save_pos(arguments.save_pos), is_save_energy(arguments.save_energy), simulation_size(arguments.size), timestep_count(arguments.steps) {

        if (is_save_pos) setup_position_file();
        if (is_save_energy) setup_energy_file();
    }

    auto save_all(System<T> const & system) -> void {
        save_points(system);
        save_energy(system);
    }

    ~Saver() {
        if (is_save_pos) pos_out_file.close();
        if (is_save_energy) energy_out_file.close();
    }
private:
    bool is_save_pos = false;
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
        pos_out_file.write(reinterpret_cast<const char*>(&simulation_size), sizeof(simulation_size));
        // number of timesteps
        pos_out_file.write(reinterpret_cast<const char*>(&timestep_count), sizeof(timestep_count));
        // size (in bytes) of data type used to store positions
        pos_out_file.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size));
    }

    auto setup_energy_file() -> void {
        energy_out_file = std::ofstream("energy.bin", std::ios::out | std::ios::binary);

        // store position file parameters
        // number of timesteps
        energy_out_file.write(reinterpret_cast<const char*>(&timestep_count), sizeof(timestep_count));
        // size (in bytes) of data type used to store energy values
        energy_out_file.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size));
    }

    auto save_points(System<T> const & system) -> void {
        if (!is_save_pos) return;

        pos_out_file.write(reinterpret_cast<const char*>(system.x.data()),
                           simulation_size * data_size * std::uint32_t(2));

    }

    auto save_energy(System<T> const & system) -> void {
        if (!is_save_energy) return;

        auto [kinetic, grav] = system.calc_energies();
        energy_out_file.write(reinterpret_cast<const char*>(&kinetic), sizeof(kinetic));
        energy_out_file.write(reinterpret_cast<const char*>(&grav), sizeof(grav));
    }
};

#endif //SAVING_H
