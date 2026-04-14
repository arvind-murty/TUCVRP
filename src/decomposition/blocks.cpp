#include "internal.hpp"

#include <stdexcept>
#include <vector>

namespace tucvrp {

namespace {

std::vector<bool> component_membership(const RootedTreeData& rooted_tree, const Component& component) {
    std::vector<bool> in_component(rooted_tree.parent.size(), false);
    for (const int v : component.vertices) {
        // Components already store their explicit vertex sets, so convert that list into a fast
        // membership bitmap for the block-decomposition helpers below.
        in_component[v] = true;
    }
    return in_component;
}

// Build the subtree T_U spanning the vertices in U from Section 4.1:
// the big terminals in the component, the component root, and possibly the component exit.
std::vector<bool> spanning_subtree_for_block_decomposition(const RootedTreeData& rooted_tree,
                                                           const Component& component,
                                                           const std::vector<bool>& in_component,
                                                           double gamma0) {
    std::vector<bool> in_u(rooted_tree.parent.size(), false);
    // By definition, the component root always belongs to U.
    in_u[component.root] = true;
    if (component.exit != -1) {
        // Internal components also include their exit vertex in U.
        in_u[component.exit] = true;
    }

    for (const int v : component.vertices) {
        // Big terminals are exactly the terminals of this component whose demand exceeds Gamma_0.
        if (v != component.root && v != component.exit && rooted_tree.is_terminal(v) &&
            rooted_tree.demands[v] > gamma0) {
            in_u[v] = true;
        }
    }

    std::vector<bool> in_tu(rooted_tree.parent.size(), false);
    for (const int start : component.vertices) {
        if (!in_u[start]) {
            continue;
        }
        // The subtree spanning U is the union of all root-to-u paths for u in U.
        // Since the input is already rooted, each such path is recovered by repeatedly following parent pointers.
        for (int v = start; v != component.root; v = rooted_tree.parent[v]) {
            if (v == -1 || !in_component[v]) {
                throw std::logic_error("component root is not an ancestor of a block key vertex");
            }
            in_tu[v] = true;
        }
        // Ensure the component root itself is marked once at the end of every path.
        in_tu[component.root] = true;
    }

    return in_tu;
}

std::vector<bool> block_key_vertices(const RootedTreeData& rooted_tree,
                                     const Component& component,
                                     const std::vector<bool>& in_component,
                                     const std::vector<bool>& in_tu,
                                     double gamma0) {
    std::vector<bool> is_block_key(rooted_tree.parent.size(), false);

    for (const int v : component.vertices) {
        if (!in_tu[v]) {
            // Vertices outside T_U do not participate in the block split.
            continue;
        }

        // A vertex belongs to U if it is the component root, the component exit,
        // or a big terminal in this component.
        bool in_u = (v == component.root || v == component.exit);
        if (rooted_tree.is_terminal(v) && v != component.root && v != component.exit &&
            rooted_tree.demands[v] > gamma0) {
            in_u = true;
        }

        int tu_child_count = 0;
        for (const int child : rooted_tree.children[v]) {
            if (in_component[child] && in_tu[child]) {
                // Only children that stay inside both the component and T_U count toward
                // the "two children in T_U" condition from Section 4.1.
                ++tu_child_count;
            }
        }

        // Section 4.1: a key vertex is either in U or has two children in T_U.
        if (in_u || tu_child_count >= 2) {
            is_block_key[v] = true;
        }
    }

    return is_block_key;
}

double block_demand(const RootedTreeData& rooted_tree,
                    const std::vector<int>& vertices,
                    int root,
                    int exit) {
    double demand = 0.0;
    for (const int v : vertices) {
        // By definition, only terminals strictly inside the block contribute to its demand.
        // So terminals at the root or the exit are excluded.
        if (v != root && v != exit && rooted_tree.is_terminal(v)) {
            demand += rooted_tree.demands[v];
        }
    }
    return demand;
}

void decompose_component_into_blocks(TreeDecomposition& decomposition,
                                     const RootedTreeData& rooted_tree,
                                     int component_id,
                                     double epsilon,
                                     double gamma) {
    const Component& component = decomposition.components[component_id];
    const double gamma0 = decomposition_detail::big_terminal_threshold(epsilon, gamma);
    // Precompute the three ingredients from Section 4.1:
    // 1. which vertices belong to the component,
    // 2. which vertices belong to the subtree T_U spanning U,
    // 3. which vertices of T_U are block key vertices.
    const std::vector<bool> in_component = component_membership(rooted_tree, component);
    const std::vector<bool> in_tu =
        spanning_subtree_for_block_decomposition(rooted_tree, component, in_component, gamma0);
    const std::vector<bool> is_block_key =
        block_key_vertices(rooted_tree, component, in_component, in_tu, gamma0);

    for (const int v : component.vertices) {
        if (!is_block_key[v]) {
            // Blocks start only from key vertices; non-key vertices are absorbed into a neighboring block.
            continue;
        }

        for (const int child : rooted_tree.children[v]) {
            if (!in_component[child]) {
                // Ignore edges leaving the current component.
                continue;
            }

            // Each child edge out of a key vertex starts one candidate block.
            // We walk downward until we hit another key vertex; everything in between belongs to this block.
            std::vector<int> block_vertices{v};
            int exit = -1;
            std::vector<int> stack{child};
            while (!stack.empty()) {
                const int u = stack.back();
                stack.pop_back();
                block_vertices.push_back(u);

                if (is_block_key[u]) {
                    // The first key vertex reached below v is the unique exit of this block.
                    // Because T_U is a tree, encountering a second distinct exit would violate the definition.
                    if (exit != -1 && exit != u) {
                        throw std::logic_error("block decomposition found multiple exit vertices");
                    }
                    exit = u;
                    continue;
                }

                for (auto it = rooted_tree.children[u].rbegin(); it != rooted_tree.children[u].rend(); ++it) {
                    if (in_component[*it]) {
                        // Continue the block through all non-key descendants that still stay in the component.
                        stack.push_back(*it);
                    }
                }
            }

            // The collected vertices form a maximally connected subgraph in which each key vertex
            // appears only as a boundary vertex of degree 1, exactly as required by Section 4.1.
            decomposition_detail::append_block(
                decomposition,
                component_id,
                v,
                exit,
                block_demand(rooted_tree, block_vertices, v, exit),
                std::move(block_vertices));
        }
    }
}

}  // namespace

void DecompositionBuilder::decompose_components_into_blocks(TreeDecomposition& decomposition,
                                                            const RootedTreeData& rooted_tree,
                                                            double epsilon) {
    if (epsilon <= 0.0 || epsilon >= 1.0) {
        throw std::invalid_argument("decompose_components_into_blocks requires epsilon in (0, 1)");
    }

    decomposition.blocks.clear();
    for (auto& component : decomposition.components) {
        component.block_ids.clear();
    }

    const double gamma = 12.0 / epsilon;
    for (int component_id = 0; component_id < static_cast<int>(decomposition.components.size()); ++component_id) {
        decompose_component_into_blocks(decomposition, rooted_tree, component_id, epsilon, gamma);
    }
}

}  // namespace tucvrp
