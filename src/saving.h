#ifndef SAVING_H
#define SAVING_H

#include <fstream>

#include "arguments.h"
#include "system.h"

template<typename T>
class Saver {
public:
    explicit Saver(Arguments arguments) : simulation_size(arguments.size), timestep_count(arguments.steps) {
        out_file = std::ofstream("positions.bin", std::ios::out | std::ios::binary);

        // store position file parameters
        // number of bodies
        out_file.write(reinterpret_cast<const char*>(&simulation_size), sizeof(simulation_size));
        // number of timesteps
        out_file.write(reinterpret_cast<const char*>(&timestep_count), sizeof(timestep_count));
        // size (in bytes) of data type used to store positions
        out_file.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size));
    }

    auto save_points(System<T> const & system) -> void {
        out_file.write(reinterpret_cast<const char*>(system.positions_x.data()), simulation_size * data_size);
        out_file.write(reinterpret_cast<const char*>(system.positions_y.data()), simulation_size * data_size);
    }

    auto finish() -> void {
        out_file.close();
    }
private:
    std::uint32_t simulation_size;
    std::uint32_t timestep_count;
    std::uint32_t data_size = sizeof(T);

    std::ofstream out_file;
};

#endif //SAVING_H
