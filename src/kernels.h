#ifndef HPC_SENDERS_KERNELS_H
#define HPC_SENDERS_KERNELS_H

#include "atomic_quad_tree.h"

// raw kernels
template<typename T, typename Index_t>
void atomic_calc_mass(AtomicQuadTree<T, Index_t> tree, Index_t tree_index) {
    // If this node is not a leaf node with a body, we are done:
    if (tree.first_child[tree_index] != tree.body)
      return;

    // Accumulate masses up to the root
    do {
        // move up to parent
        tree_index = tree.parent[tree_index];
        if (tree_index == tree.empty) break;  // if reached root

	// No thread will be arriving from siblings that are empty leaves,
	// so count those:
	uint32_t local_leaf_count = 0;
	auto leaf_child_index = tree.first_child[tree_index];
	for (int i = 0; i < 4; ++i)
	  local_leaf_count += static_cast<uint32_t>(tree.first_child[leaf_child_index + i] == tree.empty);
	uint32_t expected_count = 4 - 1 - local_leaf_count;

	// Arrive at parent releasing previous masses accumulated, and acquiring masses accumulated by other threads:
        if (tree.child_mass_complete[tree_index].fetch_add(1, memory_order_acq_rel) != expected_count)
	  break;

        // pick up all child masses
	T mass = 0;
	auto centre_masses = vec<T, 2>::splat(0);
        for (auto i = 0; i < 4; i++) {
            auto child_index = tree.first_child[tree_index] + i;
            auto child_total_mass = tree.total_masses[child_index];
            mass += child_total_mass;
            centre_masses += child_total_mass * tree.centre_masses[child_index];
        }
        centre_masses /= mass;

        tree.total_masses[tree_index] = mass;
        tree.centre_masses[tree_index] = centre_masses;
    } while(true);
}

template<typename T, typename Index_t>
void atomic_insert(T mass, vec<T, 2> pos, AtomicQuadTree<T, Index_t> tree) {
    Index_t tree_index = 0;  // insert into root
    vec<T, 2> divide = tree.root_x;
    T side_length = tree.root_side_length;

    while (true) {
        // find current node status
        atomic_ref<Index_t> fc{tree.first_child[tree_index]};
        auto status = fc.load(memory_order_acquire);
        if (status != tree.empty && status != tree.body && status != tree.locked) {
	    // If the node has children
            T half_length = side_length / static_cast<T>(4);  // / 2 is for new quad length, then / 2 is for half length

	    Index_t child_pos;
            if (pos[0] < divide[0])  {  // left
                divide[0] -= half_length;
                if (pos[1] > divide[1]) {  // top left
                    divide[1] += half_length;
                    child_pos = 0;
                } else {  // bottom left
                    divide[1] -= half_length;
                    child_pos = 1;
                }
            } else { // right
                divide[0] += half_length;
                if (pos[1] > divide[1]) {  // top right
                    divide[1] += half_length;
                    child_pos = 2;
                } else {  // bottom right
                    divide[1] -= half_length;
                    child_pos = 3;
                }
            }
	    tree_index = tree.first_child[tree_index] + child_pos;
            side_length /= static_cast<T>(2);
        } else if (status == tree.empty && fc.compare_exchange_weak(status, tree.locked, memory_order_acquire)) {
            // compare_exchange_weak is fine because if it fails spuriously, we'll retry again
            tree.total_masses[tree_index] = mass;
            tree.centre_masses[tree_index] = pos;
            fc.store(tree.body, memory_order_release);
            break;
        } else if (status == tree.body && fc.compare_exchange_weak(status, tree.locked, memory_order_acquire)) {
            // compare_exchange_weak is fine because if it fails spuriously, we'll retry again

	    // create children
            Index_t first_child_index = tree.bump_allocator->fetch_add(4, memory_order_relaxed);
	    for (int i = 0; i < 4; ++i) {
	      tree.parent[first_child_index + i] = tree_index;
	      tree.first_child[first_child_index + i] = tree.empty;
	    }

            // evict body at current index and insert into children keeping node locked
            auto p_x = tree.centre_masses[tree_index];
            Index_t evicted_index = first_child_index + 2 * static_cast<Index_t>(p_x[0] >= divide[0]) + 1 * static_cast<Index_t>(p_x[1] <= divide[1]);
            tree.centre_masses[evicted_index] = p_x;
            tree.total_masses[evicted_index] = tree.total_masses[tree_index];
            tree.first_child[evicted_index] = tree.body;

            // release node and continue to try to insert body
            fc.store(first_child_index, memory_order_release);
        }
    }
}

template<typename T, typename Index_t>
vec<T, 2> bh_calc_force(vec<T, 2> x, T const theta, AtomicQuadTree<T, Index_t> const tree) {
    Index_t tree_index = 0;
    T side_length = tree.root_side_length;
    auto a = vec<T, 2>::splat(0);

    // stackless tree traversal
    bool came_forwards = true;
    while (tree_index != tree.empty) {
        Index_t next_node_index = tree.next_node(tree_index);
        if (came_forwards) {  // child or sibling node
            vec<T, 2> xj = tree.centre_masses[tree_index];
            // check if below threshold
	    auto fc = tree.first_child[tree_index];
            if (fc == tree.empty || fc == tree.body || side_length / dist(x, xj) < theta) {
                T mj = tree.total_masses[tree_index];
                a += mj * (xj - x)/ dist3(x, xj);
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

    return a;
}

// launch kernels

template<typename T, typename Index_t>
auto build_atomic_tree(System<T>& system, AtomicQuadTree<T, Index_t> tree) {
    auto r = system.body_indices();
    std::for_each(
        std::execution::par,
        r.begin(), r.end(),
        [s=system.state(),tree] (Index_t i) {
            atomic_insert<T, Index_t>(s.m[i], s.x[i], tree);
        }
    );
}


template<typename T, typename Index_t>
auto calc_mass_atomic_tree(System<T>& system, AtomicQuadTree<T, Index_t> tree) {
    auto r = system.body_indices();
    std::for_each_n(
        std::execution::par,
        r.begin(), tree.capacity,
        [tree] (auto i) {
	  atomic_calc_mass(tree, i);
        }
    );
}

template<typename T, typename Index_t>
auto calc_force_atomic_tree(System<T>& system, AtomicQuadTree<T, Index_t> const tree, T const theta) {
    auto r = system.body_indices();
    std::for_each(
        std::execution::par_unseq,
        r.begin(), r.end(),
        [s=system.state(),theta, tree] (Index_t i) {
            s.a[i] = s.c * bh_calc_force<T, Index_t>(s.x[i], theta, tree);
        }
    );
}

template<typename T, typename Index_t>
auto compute_bounded_atomic_quad_tree(System<T>& system, AtomicQuadTree<T, Index_t>& tree){
    // find the minimum and maximum xy co-ordinate
    auto r = system.body_indices();
    auto [min_size, max_size] = std::transform_reduce(
        std::execution::par_unseq,
        r.begin(), r.end(),
        std::make_tuple<T, T>(0, 0),
        [] (auto lhs, auto rhs) -> std::tuple<T, T> {
	  return {gmin(std::get<0>(lhs), std::get<0>(rhs)), gmax(std::get<1>(lhs), std::get<1>(rhs))};
        },
        [s=system.state()] (auto i) -> std::tuple<T, T> {
	  return {min(s.x[i]), max(s.x[i])};
        }
    );

    // adjust boundary
    max_size += 1;
    min_size -= 1;

    T divide = (max_size + min_size) / static_cast<T>(2);

    // add root node to tree
    tree.root_side_length = max_size - min_size;
    tree.root_x = vec<T, 2>::splat(divide);
    tree.first_child[0] = tree.empty;
}

template<typename T, typename Index_t>
auto clear_tree(System<T>& system, AtomicQuadTree<T, Index_t> tree, Index_t last_node) {
    // clear the tree, ready for next iteration
    auto r = system.body_indices();
    std::for_each_n(
        std::execution::par_unseq,
        r.begin(), last_node,
        [tree] (auto tree_index) mutable { tree.clear(tree_index); });
}

#endif //HPC_SENDERS_KERNELS_H
