#pragma once

#include "tucvrp/instance.hpp"
#include "tucvrp/solution.hpp"

namespace tucvrp {

// Deterministic adaptation of the Labbé-Laporte-Mercure tree heuristic to the
// codebase's current model: one vehicle type, capacity 1, and route-length-only objective.
class LabbeApproxSolver {
  public:
    // Solve the current instance by repeatedly eliminating leaves in the original rooted tree.
    static SolveResult solve(const Instance& instance);
};

}  // namespace tucvrp
