#pragma once

#include "tucvrp/instance.hpp"

#include <vector>

namespace tucvrp {

struct Tour {
    std::vector<int> terminals;
    double demand = 0.0;
    double cost = 0.0;
};

struct SolveResult {
    double cost = 0.0;
    std::vector<Tour> tours;
};

class ExactSolver {
  public:
    static SolveResult solve(const Instance& instance);
};

}  // namespace tucvrp
