#ifndef HPC_SENDERS_KERNELS_H
#define HPC_SENDERS_KERNELS_H

#include "atomic_quad_tree.h"

// raw kernels

template<typename T, typename Index_t>
void atomic_calc_mass(AtomicQuadTreeContainer<T, Index_t> tree) {
    Index_t tree_index = 0;

    // navigate to leaf
    while (tree.node_status[tree_index].load(memory_order_acquire) != NodeStatus::FullLeaf) {
        // work out which child to go to ...
        auto child_index = tree.first_child[tree_index];
        if      (tree.leaf_count[child_index + 0].fetch_sub(1, memory_order_relaxed) > 0) tree_index = child_index + 0;
        else if (tree.leaf_count[child_index + 1].fetch_sub(1, memory_order_relaxed) > 0) tree_index = child_index + 1;
        else if (tree.leaf_count[child_index + 2].fetch_sub(1, memory_order_relaxed) > 0) tree_index = child_index + 2;
        else if (tree.leaf_count[child_index + 3].fetch_sub(1, memory_order_relaxed) > 0) tree_index = child_index + 3;
    }

    // sum child masses
    do {
        // move up to parent
        tree_index = tree.parent[tree_index];
        if (tree_index == is_leaf<Index_t>) break;  // if reached root

        // count how many leafs this node has
        auto leaf_child_index = tree.first_child[tree_index];
        uint32_t local_leaf_count =   static_cast<uint32_t>(tree.node_status[leaf_child_index + 0].load(memory_order_acquire) == NodeStatus::EmptyLeaf)
                                    + static_cast<uint32_t>(tree.node_status[leaf_child_index + 1].load(memory_order_acquire) == NodeStatus::EmptyLeaf)
                                    + static_cast<uint32_t>(tree.node_status[leaf_child_index + 2].load(memory_order_acquire) == NodeStatus::EmptyLeaf)
                                    + static_cast<uint32_t>(tree.node_status[leaf_child_index + 3].load(memory_order_acquire) == NodeStatus::EmptyLeaf);
        if (tree.child_mass_complete[tree_index].fetch_add(1, memory_order_relaxed) != 4 - 1 - local_leaf_count) break;  // not this threads job

        // reset masses as node contains old information
        tree.total_masses[tree_index] = 0;
        tree.centre_masses_x[tree_index] = 0;
        tree.centre_masses_y[tree_index] = 0;

        // pick up all child masses
        for (auto i = 0; i < 4; i++) {
            auto child_index = tree.first_child[tree_index] + i;
            auto child_total_mass = tree.total_masses[child_index];
            tree.total_masses[tree_index] += child_total_mass;
            tree.centre_masses_x[tree_index] += child_total_mass * tree.centre_masses_x[child_index];
            tree.centre_masses_y[tree_index] += child_total_mass * tree.centre_masses_y[child_index];
        }
        tree.centre_masses_x[tree_index] /= tree.total_masses[tree_index];
        tree.centre_masses_y[tree_index] /= tree.total_masses[tree_index];
    } while(true);
}

template<typename T, typename Index_t>
void atomic_insert(
        // particle data
        T mass, T pos_x, T pos_y,
        // tree data
        AtomicQuadTreeContainer<T, Index_t> tree
) {
    Index_t tree_index = 0;  // insert into root
    T divide_x = tree.root_x;
    T divide_y = tree.root_y;
    T side_length = tree.root_side_length;

    while (true) {
        // find current node status
        auto local_node_status = tree.node_status[tree_index].load(memory_order_acquire);
        if (local_node_status == NodeStatus::NotLeaf) {  // if the node has children
            T half_length = side_length / static_cast<T>(4);  // / 2 is for new quad length, then / 2 is for half length

            if (pos_x < divide_x)  {  // left
                divide_x -= half_length;
                if (pos_y > divide_y) {  // top left
                    divide_y += half_length;
                    tree_index = tree.first_child[tree_index] + 0;
                }
                else {  // bottom left
                    divide_y -= half_length;
                    tree_index = tree.first_child[tree_index] + 1;
                }
            }
            else { // right
                divide_x += half_length;
                if (pos_y > divide_y) {  // top right
                    divide_y += half_length;
                    tree_index = tree.first_child[tree_index] + 2;
                }
                else {  // bottom right
                    divide_y -= half_length;
                    tree_index = tree.first_child[tree_index] + 3;
                }
            }

            tree.leaf_count[tree_index].fetch_add(1, memory_order_relaxed);  // count needed for mass traversal
            side_length /= static_cast<T>(2);
        } else if (local_node_status == NodeStatus::EmptyLeaf && tree.node_status[tree_index].compare_exchange_weak(local_node_status, NodeStatus::Locked, memory_order_acquire, memory_order_relaxed)) {
            tree.total_masses[tree_index] = mass;
            tree.centre_masses_x[tree_index] = pos_x;
            tree.centre_masses_y[tree_index] = pos_y;
            tree.node_status[tree_index].store(NodeStatus::FullLeaf, memory_order_release);
            break;
        } else if (local_node_status == NodeStatus::FullLeaf && tree.node_status[tree_index].compare_exchange_weak(local_node_status, NodeStatus::Locked, memory_order_acquire, memory_order_relaxed)) {
            // create children
            Index_t first_child_index = tree.bump_allocator->fetch_add(4, memory_order_relaxed);
            tree.first_child[tree_index] = first_child_index;
            T half_length = side_length / static_cast<T>(4);  // / 2 is for new quad length, then / 2 is for half length

            tree.next_nodes[first_child_index + 0] = first_child_index + 1;
            tree.next_nodes[first_child_index + 1] = first_child_index + 2;
            tree.next_nodes[first_child_index + 2] = first_child_index + 3;
            tree.next_nodes[first_child_index + 3] = tree_index;  // link back to parent

            tree.parent[first_child_index + 0] = tree_index;
            tree.parent[first_child_index + 1] = tree_index;
            tree.parent[first_child_index + 2] = tree_index;
            tree.parent[first_child_index + 3] = tree_index;

            // relaxed as values released with a later atomic
            tree.node_status[first_child_index + 0].store(NodeStatus::EmptyLeaf, memory_order_relaxed);
            tree.node_status[first_child_index + 1].store(NodeStatus::EmptyLeaf, memory_order_relaxed);
            tree.node_status[first_child_index + 2].store(NodeStatus::EmptyLeaf, memory_order_relaxed);
            tree.node_status[first_child_index + 3].store(NodeStatus::EmptyLeaf, memory_order_relaxed);
            // end of children creation

            // evict body at current index and insert into children keeping node locked
            T p_x = tree.centre_masses_x[tree_index];
            T p_y = tree.centre_masses_y[tree_index];
            Index_t evicted_index = first_child_index + 2 * static_cast<Index_t>(p_x >= divide_x) + 1 * static_cast<Index_t>(p_y <= divide_y);
            tree.centre_masses_x[evicted_index] = p_x;
            tree.centre_masses_y[evicted_index] = p_y;
            tree.total_masses[evicted_index] = tree.total_masses[tree_index];
            tree.node_status[evicted_index].store(NodeStatus::FullLeaf, memory_order_relaxed);
            tree.leaf_count[evicted_index].fetch_add(1, memory_order_relaxed);

            // release node and continue to try to insert body
            tree.node_status[tree_index].store(NodeStatus::NotLeaf, memory_order_release);
        }
    }
}

template<typename T, typename Index_t>
auto bh_calc_force(T const particle_pos_x, T const particle_pos_y, T const theta,
                   ConstAtomicQuadTreeContainer<T, Index_t> const tree
) -> std::tuple<T, T> {

    Index_t tree_index = 0;
    T side_length = tree.root_side_length;
    T particle_accel_x = 0;
    T particle_accel_y = 0;

    // stackless tree traversal
    bool came_forwards = true;
    while (tree_index != is_leaf<Index_t>) {
        Index_t next_node_index = tree.next_nodes[tree_index];
        if (came_forwards) {  // child or sibling node
            T position_x_j = tree.centre_masses_x[tree_index];
            T position_y_j = tree.centre_masses_y[tree_index];
            T distance = calc_l2_norm<T>(particle_pos_x, particle_pos_y, position_x_j, position_y_j);

            // check if below threshold
            if (tree.first_child[tree_index] == is_leaf<Index_t> || side_length / distance < theta) {
                T mass_j = tree.total_masses[tree_index];
                T cube_l2_norm = calc_cube_l2<T>(particle_pos_x, particle_pos_y, position_x_j, position_y_j);
                particle_accel_x += mass_j * (position_x_j - particle_pos_x) / cube_l2_norm;
                particle_accel_y += mass_j * (position_y_j - particle_pos_y) / cube_l2_norm;
            } else {  // visit children
                next_node_index = tree.first_child[tree_index];
                side_length /= static_cast<T>(2);
            }
        }
        // this works as children are always allocated after their parent
        came_forwards = next_node_index > tree_index;
        side_length *= came_forwards ? static_cast<T>(1) : static_cast<T>(2);
        tree_index = next_node_index;
    }

    return {particle_accel_x, particle_accel_y};
}

// launch kernels

template<typename T, typename Index_t>
auto build_atomic_tree(System<T>& system, AtomicQuadTreeContainer<T, Index_t> tree) {
    std::for_each_n(
        std::execution::par,
        std::begin(system.index), system.size,
        [
            masses=system.masses.data(),
            positions_x=system.positions_x.data(), positions_y=system.positions_y.data(),
            tree
        ] (Index_t particle_index) {
            atomic_insert<T, Index_t>(
                // particle data
                masses[particle_index], positions_x[particle_index], positions_y[particle_index],
                // tree data
                tree
            );
        }
    );
}


template<typename T, typename Index_t>
auto calc_mass_atomic_tree(System<T>& system, AtomicQuadTreeContainer<T, Index_t> tree) {
    std::for_each_n(
        std::execution::par,
        std::begin(system.index), system.size,
        [tree] (auto _) {
            atomic_calc_mass(tree);
        }
    );
}

template<typename T, typename Index_t>
auto calc_force_atomic_tree(System<T>& system, ConstAtomicQuadTreeContainer<T, Index_t> const tree, T const theta) {
    std::for_each_n(
        std::execution::par_unseq,
        std::begin(system.index), system.size,
        [
            p_xs=system.positions_x.data(), p_ys=system.positions_y.data(), constant=system.constant,
            masses=system.masses.data(), s_f_x=system.accel_x.data(), s_f_y=system.accel_y.data(),
            theta, tree
        ] (Index_t particle_index) {
            T const p_x = p_xs[particle_index];
            T const p_y = p_ys[particle_index];

            auto [a_x, a_y] = bh_calc_force<T, Index_t>(p_x, p_y, theta, tree);

            s_f_x[particle_index] = constant * a_x;
            s_f_y[particle_index] = constant * a_y;
        }
    );
}

template<typename T, typename Index_t>
auto compute_bounded_atomic_quad_tree(System<T>& system, AtomicQuadTreeContainer<T, Index_t>& tree){
    // find the minimum and maximum xy co-ordinate
    auto [min_size, max_size] = std::transform_reduce(
        std::execution::par_unseq,
        std::begin(system.index), std::next(std::begin(system.index), system.size),
        std::make_tuple<T, T>(0, 0),
        [] (auto lhs, auto rhs) -> std::tuple<T, T> {
            return {std::min<T>(std::get<0>(lhs), std::get<0>(rhs)), std::max<T>(std::get<1>(lhs), std::get<1>(rhs))};
        },
        [pos_x=system.positions_x.data(), pos_y=system.positions_y.data()] (auto index) -> std::tuple<T, T> {
            return {std::min<T>(pos_x[index], pos_y[index]), std::max<T>(pos_x[index], pos_y[index])};
        }
    );

    // adjust boundary
    max_size += 1;
    min_size -= 1;

    T divide = (max_size + min_size) / static_cast<T>(2);

    // add root node to tree
    tree.root_side_length = max_size - min_size;
    tree.root_x = divide;
    tree.root_y = divide;
    tree.next_nodes[0] = is_leaf<Index_t>;
    tree.node_status[0].store(NodeStatus::EmptyLeaf, memory_order_relaxed);
}

template<typename T, typename Index_t>
auto clear_tree(System<T>& system, AtomicQuadTreeContainer<T, Index_t> tree) {
    // clear the tree, ready for next iteration
    std::for_each_n(
        std::execution::par_unseq,
        std::begin(system.index), tree.bump_allocator->load(memory_order_acquire),
        [tree] (auto tree_index) {
            if (tree_index == 0) {
                tree.bump_allocator->store(1, memory_order_relaxed);
            }

            tree.total_masses[tree_index] = 0;
            tree.centre_masses_x[tree_index] = 0;
            tree.centre_masses_y[tree_index] = 0;
            tree.first_child[tree_index] = is_leaf<Index_t>;
            tree.next_nodes[tree_index] = is_leaf<Index_t>;

            tree.leaf_count[tree_index].store(0, memory_order_relaxed);
            tree.child_mass_complete[tree_index].store(0, memory_order_relaxed);
            tree.parent[tree_index] = is_leaf<Index_t>;
            tree.node_status[tree_index].store(NodeStatus::EmptyLeaf, memory_order_release);
    });
}

#endif //HPC_SENDERS_KERNELS_H
