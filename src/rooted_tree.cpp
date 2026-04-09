#include "tucvrp/rooted_tree.hpp"

#include <stdexcept>

namespace tucvrp {

// Check terminal membership through the cached vertex-indexed terminal flags.
bool RootedTreeData::is_terminal(int vertex) const {
    return vertex >= 0 && vertex < static_cast<int>(terminal_flags.size()) && terminal_flags[vertex];
}

// Build parent, child, and distance annotations for the rooted version of an instance.
RootedTreeData RootedTreeBuilder::build(const Instance& instance) {
    instance.validate();

    RootedTreeData rooted_tree;
    rooted_tree.depot = instance.depot();
    rooted_tree.vertices = instance.vertices();
    rooted_tree.parent = instance.parent_array();
    rooted_tree.distances_from_depot = instance.distances_from_depot();
    rooted_tree.subtree_terminal_counts = instance.subtree_terminal_counts();
    rooted_tree.demands.assign(instance.vertex_count(), 0.0);
    rooted_tree.children.resize(instance.vertex_count());
    rooted_tree.terminal_flags.assign(instance.vertex_count(), false);

    for (int v : rooted_tree.vertices) {
        for (const auto& edge : instance.neighbors(v)) {
            if (edge.to != rooted_tree.parent[v]) {
                rooted_tree.children[v].push_back(edge.to);
            }
        }
    }

    rooted_tree.terminal_vertices.reserve(instance.terminals().size());
    for (const auto& terminal : instance.terminals()) {
        rooted_tree.terminal_vertices.push_back(terminal.vertex);
        rooted_tree.demands[terminal.vertex] = terminal.demand;
        rooted_tree.terminal_flags[terminal.vertex] = true;
    }

    return rooted_tree;
}

}  // namespace tucvrp
