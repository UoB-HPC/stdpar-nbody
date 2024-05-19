#ifndef QUAD_TREE_H
#define QUAD_TREE_H
#include <vector>

#include "system.h"

template<typename Index_t>
constexpr Index_t is_leaf = static_cast<Index_t>(-1);

// Quad tree using Class of Vectors (CoV) (i.e. Structure of Arrays)
template<typename T, typename Index_t>
class VectorQuadTree {
public:
    std::vector<T> side_lengths;
    std::vector<T> divides_x;
    std::vector<T> divides_y;
    std::vector<T> total_masses;
    std::vector<T> centre_masses_x;
    std::vector<T> centre_masses_y;
    std::vector<Index_t> first_child;
    std::vector<Index_t> next_nodes;
    std::vector<Index_t> particle_indicies;
    std::vector<T> particle_pos_x;
    std::vector<T> particle_pos_y;
    std::vector<int> is_full;  // if has a particle or children
    System<T> const& system;

    explicit VectorQuadTree(System<T>& system) : system(system) {}

    auto add_node(T side_length, T divide_x, T divide_y, Index_t next_node) {
        side_lengths.push_back(side_length);
        divides_x.push_back(divide_x);
        divides_y.push_back(divide_y);
        total_masses.push_back(0);
        centre_masses_x.push_back(0);
        centre_masses_y.push_back(0);
        first_child.push_back(is_leaf<Index_t>);
        next_nodes.push_back(next_node);
        particle_indicies.push_back(0);
        particle_pos_x.push_back(0);
        particle_pos_y.push_back(0);
        is_full.push_back(false);

        size += 1;
    }

    void insert(Index_t tree_index, Index_t particle_index, T pos_x, T pos_y) {
        if (!is_full[tree_index]) {
            is_full[tree_index] = true;
            particle_pos_x[tree_index] = pos_x;
            particle_pos_y[tree_index] = pos_y;
            particle_indicies[tree_index] = particle_index;
        } else if (first_child[tree_index] == is_leaf<Index_t>) {  // create children and insert
            T quad_length = side_lengths[tree_index] / static_cast<T>(2);
            T half_length = quad_length / static_cast<T>(2);
            first_child[tree_index] = size;
            add_node(quad_length, divides_x[tree_index] - half_length, divides_y[tree_index] + half_length, size + 1);
            add_node(quad_length, divides_x[tree_index] - half_length, divides_y[tree_index] - half_length, size + 1);
            add_node(quad_length, divides_x[tree_index] + half_length, divides_y[tree_index] + half_length, size + 1);
            add_node(quad_length, divides_x[tree_index] + half_length, divides_y[tree_index] - half_length, tree_index);  // link back to parent

            insert(tree_index, particle_indicies[tree_index], particle_pos_x[tree_index], particle_pos_y[tree_index]);
            insert(tree_index, particle_index, pos_x, pos_y);
        } else {  // recurse on quadrant
            // work out which quadrant
            if (pos_x < divides_x[tree_index]) {
                // left
                if (pos_y > divides_y[tree_index]) {
                    // top
                    insert(first_child[tree_index] + 0, particle_index, pos_x, pos_y);
                } else {
                    // bottom
                    insert(first_child[tree_index] + 1, particle_index, pos_x, pos_y);
                }
            } else {
                // right
                if (pos_y > divides_y[tree_index]) {
                    // top
                    insert(first_child[tree_index] + 2, particle_index, pos_x, pos_y);
                } else {
                    // bottom
                    insert(first_child[tree_index] + 3, particle_index, pos_x, pos_y);
                }
            }
        }
    }

    // traversed the tree and calculates the centre of mass and total mass at each node
    auto calc_mass(Index_t tree_index) -> void{
        if (!is_full[tree_index]) {
            // empty node - do nothing!
        } else if (first_child[tree_index] == is_leaf<Index_t>) {
            total_masses[tree_index] = system.masses[particle_indicies[tree_index]];
            centre_masses_x[tree_index] = particle_pos_x[tree_index];
            centre_masses_y[tree_index] = particle_pos_y[tree_index];
        } else {
            // recursively calculate
            for (auto i = 0; i < 4; i++) {
                auto child_index = first_child[tree_index] + i;
                calc_mass(child_index);
                auto child_total_mass = total_masses[child_index];
                total_masses[tree_index] += child_total_mass;
                centre_masses_x[tree_index] += child_total_mass * centre_masses_x[child_index];
                centre_masses_y[tree_index] += child_total_mass * centre_masses_y[child_index];
            }
            centre_masses_x[tree_index] /= total_masses[tree_index];
            centre_masses_y[tree_index] /= total_masses[tree_index];
        }
    }

private:
    Index_t size = 0;
};

template<typename T, typename Index_t>
auto compute_bounded_quad_tree(System<T>& system) {
    T max_size = static_cast<T>(0);
    T min_size = static_cast<T>(0);
    for (auto index = 0; index < system.size; index++) {
        max_size = std::max<T>(max_size, std::max<T>(system.positions_x[index], system.positions_y[index]));
        min_size = std::min<T>(min_size, std::min<T>(system.positions_x[index], system.positions_y[index]));
    }

    max_size += 1;
    min_size -= 1;

    T divide = (max_size + min_size) / static_cast<T>(2);

    VectorQuadTree<T, Index_t> tree{system};
    // insert root node
    tree.add_node(max_size - min_size, divide, divide, is_leaf<Index_t>);
    return tree;
}

#endif //QUAD_TREE_H
