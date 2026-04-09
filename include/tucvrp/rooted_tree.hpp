#pragma once

#include "tucvrp/instance.hpp"

#include <unordered_map>
#include <vector>

namespace tucvrp {

struct RootedTreeData {
    int depot = -1;
    std::vector<int> vertices;
    std::vector<int> terminal_vertices;
    std::vector<int> parent;
    std::vector<std::vector<int>> children;
    std::vector<double> distance_from_depot;
    std::unordered_map<int, double> terminal_distances;

    // Return whether a vertex is a terminal in this rooted view.
    [[nodiscard]] bool is_terminal(int vertex) const;
};

class RootedTreeBuilder {
  public:
    // Build rooted-tree annotations from an instance.
    static RootedTreeData build(const Instance& instance);
};

}  // namespace tucvrp
