#ifndef SYSTEM_H
#define SYSTEM_H

#include <cmath>
#include <functional>
#include <random>
#include <vector>
#include "vec.h"

template<typename T>
class System {
public:
    using index_t = uint32_t;
    index_t const size;
    index_t const max_tree_node_size;
    T const dt;
    T const constant;
    std::vector<T> m;
    std::vector<vec<T, 2>> x, v, a, ao;

    // random generation
    std::mt19937 gen{42};  // fix random generation
    std::uniform_real_distribution<> angle_dis{0, 2 * std::numbers::pi};
    std::uniform_real_distribution<> unit_dis{0, 1};

    System(index_t size, T time_step, T constant, index_t max_tree_node_size):
        size(size), max_tree_node_size(max_tree_node_size),
	dt(time_step), constant(constant), m(size),
        x(size), v(size), a(size), ao(size)
    {}

    auto body_indices() const { return std::views::iota(index_t(0), size); }

    // Helper to make it easier to access all state from parallel algorithms
    struct state_t {
      T* m;
      vec<T, 2>* x, *v, *a, *ao;
      T dt, c;
      index_t sz;
    };
    state_t state() {
      return { .m = m.data(), .x = x.data(), .v = v.data(), .a = a.data(), .ao = ao.data(),
	       .dt = dt, .c = constant, .sz = size };
    }
    
    void accelerate_step() {
        // performs leap frog integration
        auto r = body_indices();
        std::for_each(
            std::execution::par_unseq,
	    r.begin(), r.end(),
            [s = state()] (auto i) {
		s.x[i] += s.dt * s.v[i] + T(0.5) * s.dt * s.dt * s.ao[i];
                s.v[i] += T(0.5) * s.dt * (s.a[i] + s.ao[i]);
                s.ao[i] = s.a[i];
        });
    }

    auto calc_energies() const -> std::tuple<T, T> {
        auto r = body_indices();
        T kinetic_energy = static_cast<T>(0.5) * std::transform_reduce(
            std::execution::par_unseq,
            r.begin(), r.end(),
            static_cast<T>(0), std::plus<T>{},
            [m=m.data(), v=v.data()] (auto i) { return m[i] * l2norm2(v[i]); }
        );
        T gravitational_energy = -static_cast<T>(0.5) * constant * std::transform_reduce(
            std::execution::par_unseq,
            r.begin(), r.end(),
            static_cast<T>(0), std::plus<T>{},
            [m=m.data(), x=x.data(), sz=size] (auto i) {
                T total = 0;
                T mi = m[i];
                auto xi = x[i];
                for (index_t j = 0; j < sz; j++) {
                    if (j != i) total += mi * m[j] / dist(xi, x[j]);
                }
                return total;
            }
        );

        return { kinetic_energy, gravitational_energy };
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
        auto i = next_point;
        next_point += 1;

        m[i] = mass;
        x[i] = vec<T, 2>{{p_x, p_y}};
        v[i] = vec<T, 2>{{v_x, v_y}};
    }

    void print() const {
        for (size_t i = 0; i < size; i++) {
            std::cout << std::format(
                "{:02}: m={: .3e}, p=({: .3e}, {: .3e}), v=({: .3e}, {: .3e}), f=({: .3e}, {: .3e})",
                i, m[i], x[i][0], x[i][1],
                v[i][0], v[i][1], a[i][0], a[i][1]
            ) << std::endl;
        }
    }
private:
    std::size_t next_point = 0;
};

#endif //SYSTEM_H
