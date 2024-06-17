#ifndef SYSTEM_H
#define SYSTEM_H

#include <cmath>
#include <functional>
#include <random>
#include <vector>

template<typename T>
auto calc_cube_l2(T x_i, T y_i, T x_j, T y_j) -> T {
    T total = std::pow<T>(std::abs(x_i - x_j), static_cast<T>(2)) + std::pow<T>(std::abs(y_i - y_j), static_cast<T>(2));
    return std::pow<T>(total, static_cast<T>(3.0) / static_cast<T>(2.0)) + std::numeric_limits<T>::epsilon();  // to allow division by 0
}

template<typename T>
auto calc_l2_norm(T x_i, T y_i, T x_j, T y_j) -> T {
    T total = std::pow<T>(std::abs(x_i - x_j), static_cast<T>(2)) + std::pow<T>(std::abs(y_i - y_j), static_cast<T>(2));
    return std::pow<T>(total, static_cast<T>(1.0) / static_cast<T>(2.0)) + std::numeric_limits<T>::epsilon();  // to allow division by 0
}

template<typename T>
class System {
public:
    using index_t = uint32_t;
    index_t const size;
    index_t const max_tree_node_size;
    T const time_step;
    T const constant;
    std::vector<T> masses;
    std::vector<T> positions_x;
    std::vector<T> positions_y;
    std::vector<T> velocities_x;
    std::vector<T> velocities_y;
    std::vector<T> accel_x;
    std::vector<T> accel_y;
    std::vector<T> accel_old_x;
    std::vector<T> accel_old_y;

    // random generation
    std::mt19937 gen{42};  // fix random generation
    std::uniform_real_distribution<> angle_dis{0, 2 * std::numbers::pi};
    std::uniform_real_distribution<> unit_dis{0, 1};

    System(index_t size, T time_step, T constant, index_t max_tree_node_size):
        size(size), max_tree_node_size(max_tree_node_size),
	time_step(time_step), constant(constant), masses(size),
        positions_x(size), positions_y(size),
        velocities_x(size), velocities_y(size),
        accel_x(size), accel_y(size),
        accel_old_x(size), accel_old_y(size)
    {}

    auto body_indices() { return std::views::iota(index_t(0), size); }

    void accelerate_step() {
        // performs leap frog integration
        auto r = body_indices();
        std::for_each(
            std::execution::par_unseq,
	    r.begin(), r.end(),
            [
                masses=masses.data(), accel_x=accel_x.data(), accel_y=accel_y.data(),
                accel_old_x=accel_old_x.data(), accel_old_y=accel_old_y.data(),
                positions_x=positions_x.data(), positions_y=positions_y.data(),
                velocities_x=velocities_x.data(), velocities_y=velocities_y.data(), time_step=time_step
            ] (auto i) {
                positions_x[i] += velocities_x[i] * time_step + static_cast<T>(0.5) * accel_old_x[i] * time_step * time_step;
                positions_y[i] += velocities_y[i] * time_step + static_cast<T>(0.5) * accel_old_y[i] * time_step * time_step;

                velocities_x[i] += static_cast<T>(0.5) * (accel_x[i] + accel_old_x[i]) * time_step;
                velocities_y[i] += static_cast<T>(0.5) * (accel_y[i] + accel_old_y[i]) * time_step;

                accel_old_x[i] = accel_x[i];
                accel_old_y[i] = accel_y[i];
        });
    }

    auto calc_energies() const {
        auto r = body_indices();
        T kinetic_engery = static_cast<T>(0.5) * std::transform_reduce(
            std::execution::par_unseq,
            r.begin(), r.end(),
            static_cast<T>(0), std::plus<T>{},
            [masses=masses.data(), velocities_x=velocities_x.data(), velocities_y=velocities_y.data()] (auto i) {
                return masses[i] * (velocities_x[i] * velocities_x[i] + velocities_y[i] * velocities_y[i]);
            }
        );
        T gravitational_energy = -static_cast<T>(0.5) * constant * std::transform_reduce(
            std::execution::par_unseq,
            r.begin(), r.end(),
            static_cast<T>(0), std::plus<T>{},
            [size=size, masses=masses.data(), pos_x=positions_x.data(), pos_y=positions_y.data()] (auto i) {
                T total = 0;
                T mass_i = masses[i];
                T pos_x_i = pos_x[i];
                T pos_y_i = pos_y[i];
                for (index_t j = 0; j < size; j++) {
                    if (j != i) {
                        total += mass_i * masses[j] / calc_l2_norm<T>(pos_x_i, pos_y_i, pos_x[j], pos_y[j]);
                    }
                }
                return total;
            }
        );

        std::cout << std::format("Kinetic Enegy: {}\n", kinetic_engery)
                  << std::format("Gravitational Potential Enegy: {}\n", gravitational_energy)
                  << std::format("Total Energy: {}\n", kinetic_engery + gravitational_energy);
    }

    void add_point(T mass, T distance, T total_mass) {
        auto const angle = angle_dis(gen);
        auto const x = distance * std::sin(angle);
        auto const y = distance * std::cos(angle);

        // assume orbit around centre of system
        auto const velocity = std::sqrt(constant * total_mass / (distance + std::numeric_limits<T>::epsilon()));
        auto const norm = std::sqrt(x * x + y * y) + std::numeric_limits<T>::epsilon();
        auto v_x = -y / norm * velocity;
        auto v_y =  x / norm * velocity;

        add_point(mass, x, y, v_x, v_y);
    }

    void add_point(T mass, T p_x, T p_y, T v_x, T v_y) {
        auto index = next_point;
        next_point += 1;

        masses[index] = mass;
        positions_x[index] = p_x;
        positions_y[index] = p_y;
        velocities_x[index] = v_x;
        velocities_y[index] = v_y;
    }

    void print() const {
        for (size_t i = 0; i < size; i++) {
            std::cout << std::format(
                "{:02}: m={: .3e}, p=({: .3e}, {: .3e}), v=({: .3e}, {: .3e}), f=({: .3e}, {: .3e})",
                i, masses[i], positions_x[i], positions_y[i],
                velocities_x[i], velocities_y[i], accel_x[i], accel_y[i]
            ) << std::endl;
        }
    }
private:
    std::size_t next_point = 0;
};

#endif //SYSTEM_H
