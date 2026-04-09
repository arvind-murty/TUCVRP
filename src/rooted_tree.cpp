#include "tucvrp/rooted_tree.hpp"

#include <stdexcept>

namespace tucvrp {

// Check terminal membership through the cached terminal-distance map.
bool RootedTreeData::is_terminal(int vertex) const { return terminal_distances.contains(vertex); }

// Build parent, child, and distance annotations for the rooted version of an instance.
RootedTreeData RootedTreeBuilder::build(const Instance& instance) {
    instance.validate();

    RootedTreeData rooted_tree;
    rooted_tree.depot = instance.depot();
    rooted_tree.vertices = instance.vertices();
    rooted_tree.parent = instance.parent_array();
    rooted_tree.distances_from_depot = instance.distances_from_depot();
    rooted_tree.children.resize(instance.vertex_count());
    rooted_tree.terminal_distances = instance.terminal_distances();

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
    }

    return rooted_tree;
}

}  // namespace tucvrp
