#pragma once
#include <cmath>
#include <functional>
#include <random>
#include <vector>
#include "vec.h"

template<typename T, dim_t N>
class System {
public:
    using index_t = uint32_t;
    index_t const size;
    index_t const max_tree_node_size;
    T const dt;
    T const constant;
    std::vector<T> m;
    std::vector<vec<T, N>> x, v, a, ao;

    // random generation
    std::mt19937 gen{42};  // fix random generation
    std::uniform_real_distribution<> angle_dis{0, 2 * std::numbers::pi};
    std::uniform_real_distribution<> unit_dis{0, 1};
    std::uniform_real_distribution<> sym_dis{-1, 1};

    System(index_t size, T time_step, T constant):
        size(size), max_tree_node_size(std::max<std::size_t>(child_count<N> * size, 1000)),
	dt(time_step), constant(constant), m(size),
        x(size), v(size), a(size), ao(size)
    {}

    auto body_indices() const { return std::views::iota(index_t(0), size); }

    // Helper to make it easier to access all state from parallel algorithms
    struct state_t {
      T* m;
      vec<T, N>* x, *v, *a, *ao;
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

    void add_point(T mass, vec<T, N> pos, vec<T, N> vel) {
        auto i = next_point;
        next_point += 1;

        m[i] = mass;
        x[i] = pos;
        v[i] = vel;
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
