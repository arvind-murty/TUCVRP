#pragma once

#include "tucvrp/decomposition.hpp"
#include "tucvrp/instance.hpp"
#include "tucvrp/rooted_tree.hpp"
#include "tucvrp/solution.hpp"

#include <string>

namespace tucvrp {

// Tunable parameters for the paper's (1.5 + epsilon)-approximation.
struct OnePointFiveApproxParams {
    double epsilon = 0.1;
};

// Main entry point for the paper algorithm under development.
class OnePointFiveApproxSolver {
  public:
    // Solve using the paper algorithm with default approximation parameters.
    static SolveResult solve(const Instance& instance);
    // Solve using the paper algorithm with explicit approximation parameters.
    static SolveResult solve(const Instance& instance, const OnePointFiveApproxParams& params);
    // Solve the bounded-distance subproblem after the outer reduction step.
    static SolveResult solve_bounded_distance(const Instance& instance, const OnePointFiveApproxParams& params);
};

}  // namespace tucvrp
