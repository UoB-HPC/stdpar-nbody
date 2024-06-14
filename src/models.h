#ifndef MODELS_H
#define MODELS_H

#include "system.h"

template<typename T>
auto build_solar_system() {
    auto timestep = 60 * 60;  // an hour in seconds

    auto system = System<T>(9, timestep, static_cast<T>(6.674e-11), 100);

    auto total_mass = 1.989e30;

    // The Sun
    system.add_point(1.989e30, 0, total_mass);

    // https://nssdc.gsfc.nasa.gov/planetary/factsheet/

    // Mercury
    system.add_point(0.33e24, 57.9e9, total_mass);
    // Venus
    system.add_point(4.87e24, 108.2e9, total_mass);
    // Earth
    system.add_point(5.97e24, 149.6e9, total_mass);
    // Mars
    system.add_point(0.642e24, 228e9, total_mass);
    // Jupiter
    system.add_point(1898e24, 778.5e9, total_mass);
    // Saturn
    system.add_point(568e24, 1432e9, total_mass);
    // Uranus
    system.add_point(86.8e24, 2867e9, total_mass);
    // Neptune
    system.add_point(102e24, 4515e9, total_mass);

    return system;
}

template<typename T>
auto build_plummer_model(Arguments arguments) {
    auto number_of_particles = arguments.size;

    auto system = System<T>(number_of_particles, 1000, static_cast<T>(6.674e-11), std::max<std::size_t>(4 * number_of_particles, 1000));

    for (size_t _ = 0; _ < number_of_particles; _++) {
        auto mass = 1.0 / static_cast<T>(number_of_particles);
        auto radius = 1.0 / std::sqrt(std::pow(system.unit_dis(system.gen), -2.0 / 3.0) - 1);
        auto pos_angle = system.angle_dis(system.gen);

        auto const x = radius * std::sin(pos_angle);
        auto const y = radius * std::cos(pos_angle);

        // rejection sampling to determine velocity
        auto q = 0.0;
        auto g = 0.1;
        while (g > q * q * std::pow(1.0 - q * q, 3.5)) {
            q = system.unit_dis(system.gen);
            g = 0.1 * system.unit_dis(system.gen);
        }
        auto velocity = q * std::numbers::sqrt2 * std::pow(radius * radius + 1, -0.25);
        auto velocity_angle = system.angle_dis(system.gen);
        auto v_x = velocity * std::sin(velocity_angle);
        auto v_y = velocity * std::cos(velocity_angle);

        system.add_point(mass, x, y, v_x, v_y);
    }

    return system;
}

template<typename T>
auto create_circular_orbit(System<T>& system, std::size_t number_of_particles, T total_mass, T orbit_mass, T centre_x, T centre_y) {
    for (size_t _ = 0; _ < number_of_particles; _++) {
        auto mass = orbit_mass / static_cast<T>(number_of_particles);
        auto radius = 30 + 20 * system.unit_dis(system.gen);

        // calculate orbit start position
        auto const angle = system.angle_dis(system.gen);
        auto const x = radius * std::sin(angle);
        auto const y = radius * std::cos(angle);

        // assume orbit around centre of system
        auto const velocity = std::sqrt(system.constant * total_mass / (radius + std::numeric_limits<T>::epsilon()));
        auto const norm = std::sqrt(x * x + y * y) + std::numeric_limits<T>::epsilon();
        auto v_x = -y / norm * velocity;
        auto v_y =  x / norm * velocity;

        // adjust position based on centre
        system.add_point(mass, x + centre_x, y + centre_y, v_x, v_y);
    }
}

template<typename T>
auto build_galaxy_model(Arguments arguments) {
    // create model of two spinning 'galaxies' that collide

    auto galaxy_number_of_particles = arguments.size / 2.0;
    auto system = System<T>(2 * galaxy_number_of_particles, 1e1, 1e-4, 10 * galaxy_number_of_particles);

    auto centre_mass = 1e4;

    auto offset = 100.0;

    system.add_point(centre_mass, -offset, offset / 2, 0, 0);
    create_circular_orbit<T>(
        system, galaxy_number_of_particles - 1,
        centre_mass + 1, 1,
        -offset, offset / 2
    );

    centre_mass /= 10;
    system.add_point(centre_mass, offset, -offset / 2, 0, 0);
    create_circular_orbit<T>(
        system, galaxy_number_of_particles - 1,
        centre_mass + 1, 1,
        offset, -offset / 2
    );

    return system;
}

#endif //MODELS_H
