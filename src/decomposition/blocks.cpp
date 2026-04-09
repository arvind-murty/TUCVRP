#include "internal.hpp"

#include <stdexcept>
#include <vector>

namespace tucvrp {

namespace {

void decompose_component_into_blocks(TreeDecomposition& decomposition,
                                     const RootedTreeData& rooted_tree,
                                     int component_id,
                                     double epsilon,
                                     double gamma) {
    using namespace decomposition_detail;

    const Component& component = decomposition.components[component_id];
    const double gamma0 = big_terminal_threshold(epsilon, gamma);
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
            append_block(
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
