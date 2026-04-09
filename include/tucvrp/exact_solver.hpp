#pragma once

#include "tucvrp/instance.hpp"
#include "tucvrp/solution.hpp"

namespace tucvrp {

class ExactSolver {
  public:
    // Solve small instances exactly by partitioning terminals into feasible tours.
    static SolveResult solve(const Instance& instance);
};

}  // namespace tucvrp
