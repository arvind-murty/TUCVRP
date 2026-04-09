#pragma once

#include "tucvrp/instance.hpp"

#include <vector>

namespace tucvrp {

struct BoundedDistanceStats {
    double min_distance = 0.0;
    double max_distance = 0.0;
    bool bounded = false;
};

class Preprocessor {
  public:
    // Compute terminal distance statistics and check the paper's bounded-distance condition.
    static BoundedDistanceStats bounded_distance_stats(const Instance& instance, double epsilon);
    // Compute the standard edge-load lower bound from subtree demands.
    static double edge_load_lower_bound(const Instance& instance);
    // Apply preprocessing so the output tree is binary and every terminal is a leaf.
    static Instance make_binary_leaf_tree(const Instance& instance);
};

}  // namespace tucvrp
