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
    static BoundedDistanceStats bounded_distance_stats(const Instance& instance, double epsilon);
    static double edge_load_lower_bound(const Instance& instance);
    static Instance binarize_tree(const Instance& instance);
};

}  // namespace tucvrp
