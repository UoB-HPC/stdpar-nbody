#pragma once
#include <format>
#include <stdexcept>
#include "arguments.h"
#include "system.h"

template<auto Value, auto Expected>
concept is_equal = (Value == Expected);

template<typename T, dim_t N>
auto build_uniform_model(Arguments arguments) -> System<T, N> {
    auto number_of_particles = arguments.size;
    auto system = System<T, N>(number_of_particles, 1e-1, 1);
    for (std::size_t _ = 0; _ < number_of_particles; _++) {
        T const mass = 1.0 / static_cast<T>(number_of_particles);
        auto pos = vec<T, N>::splat(0);
        auto vel = vec<T, N>::splat(0);
        for (std::size_t i = 0; i < N; i++) {
            pos[i] = system.sym_dis(system.gen);
            vel[i] = system.sym_dis(system.gen);
        }
        system.add_point(mass, pos, vel);
    }

    return system;
}

template<typename T, dim_t N>
requires is_equal<N, 3>
auto build_plummer_model(Arguments arguments) -> System<T, N> {
    auto number_of_particles = arguments.size;

    auto system = System<T, N>(number_of_particles, 1, static_cast<T>(6.674e-11));

    for (size_t _ = 0; _ < number_of_particles; _++) {
        T const mass = 1.0 / static_cast<T>(number_of_particles);
        T const radius = 1.0 / std::sqrt(std::pow(system.unit_dis(system.gen), -2.0 / 3.0) - 1);
        T const p_theta = std::acos(system.sym_dis(system.gen));
        T const p_phi = system.angle_dis(system.gen);

        auto pos = radius * vec<T, N>{{
            std::sin(p_theta) * std::cos(p_phi),
            std::sin(p_theta) * std::sin(p_phi),
            std::cos(p_theta)
        }};

        // rejection sampling to determine velocity
        T q = 0.0;
        T g = 0.1;
        while (g > q * q * std::pow(1.0 - q * q, 3.5)) {
            q = system.unit_dis(system.gen);
            g = 0.1 * system.unit_dis(system.gen);
        }
        T const velocity_norm = q * std::numbers::sqrt2 * std::pow(radius * radius + 1, -0.25);
        T const v_theta = std::acos(system.sym_dis(system.gen));
        T const v_phi = system.angle_dis(system.gen);

        auto velocity = velocity_norm * vec<T, N>{{
            std::sin(v_theta) * std::cos(v_phi),
            std::sin(v_theta) * std::sin(v_phi),
            std::cos(v_theta)
        }};

        system.add_point(mass, pos, velocity);
    }

    return system;
}

template<typename T, dim_t N>
auto build_plummer_model(Arguments) -> System<T, N> {
    throw std::runtime_error(std::format("Cannot build Plummer model for D={}", N));
}

template<typename T>
auto rotate_vec(std::array<std::array<T, 3>, 3> const & matrix, vec<T, 3> v) -> vec<T, 3> {
    auto result = vec<T, 3>::splat(0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i] += matrix[i][j] * v[j];
        }
    }
    return result;
}

template<typename T, dim_t N>
auto create_circular_orbit(System<T, N>& system, std::size_t number_of_particles, T total_mass, T orbit_mass, vec<T, N> centre_pos) {
    for (size_t _ = 0; _ < number_of_particles; _++) {
        T const mass = orbit_mass / static_cast<T>(number_of_particles);
        T const radius = 30 + 20 * system.unit_dis(system.gen);

        // calculate orbit start position
        T const angle = system.angle_dis(system.gen);
        auto pos = radius * vec<T, N>{{
            std::sin(angle),
            std::cos(angle)
        }};

        // assume orbit around centre of system
        T const velocity_norm = std::sqrt(system.constant * total_mass / (radius + std::numeric_limits<T>::epsilon()));
        auto velocity = velocity_norm / (l2norm(pos) + std::numeric_limits<T>::epsilon()) * vec<T, N>{{
            -pos[1],
             pos[0]
        }};

        if constexpr (N == 3) {
            pos[2] = 10 * system.sym_dis(system.gen);
            velocity[2] = 0.00001 * system.sym_dis(system.gen);

            // roate the positions so they look better in 3d
            std::array<std::array<T, 3>, 3> const matrix = {{
                {0.0, -1.0, 0.0},
                {0.9, 0.0, 0.5},
                {0.5, 0.0, 0.9}
            }};

            pos = rotate_vec(matrix, pos);
            velocity = rotate_vec(matrix, velocity);
        }

        // adjust position based on centre
        system.add_point(mass, pos + centre_pos, velocity);
    }
}

template<typename T, dim_t N>
requires is_equal<N, 2> || is_equal<N, 3>
auto build_galaxy_model(Arguments arguments) -> System<T, N> {
    // create model of two spinning 'galaxies' that collide

    auto galaxy_number_of_particles = arguments.size / 2.0;
    auto system = System<T, N>(2 * galaxy_number_of_particles, 1e1, 1e-4);

    T centre_mass = 1e4;

    T offset = 100.0;

    // galaxy 1
    auto offset_pos = offset * vec<T, N>{{-1, 1 / 2.0}};
    system.add_point(centre_mass, offset_pos, vec<T, N>::splat(0));
    create_circular_orbit<T, N>(
        system, galaxy_number_of_particles - 1,
        centre_mass + 1, 1,
        offset_pos
    );

    // galaxy 2
    centre_mass /= 10;
    offset_pos = offset * vec<T, N>{{1, -1 / 2.0}};
    system.add_point(centre_mass, offset_pos, vec<T, N>::splat(0));
    create_circular_orbit<T, N>(
        system, galaxy_number_of_particles - 1,
        centre_mass + 1, 1,
        offset_pos
    );

    return system;
}

template<typename T, dim_t N>
auto build_galaxy_model(Arguments) -> System<T, N> {
    throw std::runtime_error(std::format("Cannot build Galaxy model for D={}", N));
}
