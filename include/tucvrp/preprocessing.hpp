#pragma once

#include "tucvrp/instance.hpp"

#include <vector>

namespace tucvrp {

struct BoundedDistanceStats {
    double min_distance = 0.0;
    double max_distance = 0.0;
    bool bounded = false;
};

struct Component {
    int root = -1;
    int exit = -1;
    std::vector<int> vertices;
};

class Preprocessor {
  public:
    // Compute terminal distance statistics and check the paper's bounded-distance condition.
    static BoundedDistanceStats bounded_distance_stats(const Instance& instance, double epsilon);
    // Compute the standard edge-load lower bound from subtree demands.
    static double edge_load_lower_bound(const Instance& instance);
    // Replace high-degree branching by zero-cost auxiliary nodes so every vertex has at most two children.
    static Instance make_tree_binary(const Instance& instance);
};

}  // namespace tucvrp
