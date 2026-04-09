#pragma once

#include "tucvrp/instance.hpp"

#include <vector>

namespace tucvrp {

// Cached rooted-tree annotations derived from an instance.
struct RootedTreeData {
    // Depot/root of the rooted tree.
    int depot = -1;
    // Active vertices in the rooted tree.
    std::vector<int> vertices;
    // Terminal vertices in the same order as the instance terminals list.
    std::vector<int> terminal_vertices;
    // Parent of each vertex, with -1 at the depot.
    std::vector<int> parent;
    // Rooted children of each vertex.
    std::vector<std::vector<int>> children;
    // Depot distance for every vertex id.
    std::vector<double> distances_from_depot;
    // Number of terminals in the subtree rooted at each vertex id.
    std::vector<int> subtree_terminal_counts;
    // Terminal membership for each vertex id.
    std::vector<bool> terminal_flags;

    // Return whether a vertex is a terminal in this rooted view.
    [[nodiscard]] bool is_terminal(int vertex) const;
};

class RootedTreeBuilder {
  public:
    // Build rooted-tree annotations from an instance.
    static RootedTreeData build(const Instance& instance);
};

}  // namespace tucvrp
