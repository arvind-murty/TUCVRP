#include "tucvrp/solver.hpp"

#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace tucvrp {

SolveResult ExactSolver::solve(const Instance& instance) {
    instance.validate();
    const auto& terminals = instance.terminals();
    const int m = static_cast<int>(terminals.size());
    if (m == 0) {
        return {};
    }
    if (m > 20) {
        throw std::invalid_argument("exact solver is intended only for small instances (<= 20 terminals)");
    }

    const int mask_count = 1 << m;
    std::vector<double> subset_demand(mask_count, 0.0);
    std::vector<double> subset_cost(mask_count, 0.0);
    std::vector<std::vector<int>> subset_terminals(mask_count);

    for (int mask = 1; mask < mask_count; ++mask) {
        int bit = std::countr_zero(static_cast<unsigned int>(mask));
        int prev = mask & (mask - 1);
        subset_demand[mask] = subset_demand[prev] + terminals[bit].demand;
        subset_terminals[mask] = subset_terminals[prev];
        subset_terminals[mask].push_back(terminals[bit].vertex);
        if (subset_demand[mask] <= 1.0 + 1e-9) {
            subset_cost[mask] = instance.steiner_cost_for_terminal_subset(subset_terminals[mask]);
        } else {
            subset_cost[mask] = std::numeric_limits<double>::infinity();
        }
    }

    std::vector<double> dp(mask_count, std::numeric_limits<double>::infinity());
    std::vector<int> choice(mask_count, -1);
    dp[0] = 0.0;

    for (int mask = 1; mask < mask_count; ++mask) {
        for (int sub = mask; sub > 0; sub = (sub - 1) & mask) {
            if (!std::isfinite(subset_cost[sub])) {
                continue;
            }
            double candidate = dp[mask ^ sub] + subset_cost[sub];
            if (candidate < dp[mask]) {
                dp[mask] = candidate;
                choice[mask] = sub;
            }
        }
    }

    SolveResult result;
    result.cost = dp[mask_count - 1];

    int mask = mask_count - 1;
    while (mask != 0) {
        int sub = choice[mask];
        if (sub <= 0) {
            throw std::logic_error("failed to reconstruct exact solution");
        }
        Tour tour;
        tour.cost = subset_cost[sub];
        for (int i = 0; i < m; ++i) {
            if (sub & (1 << i)) {
                tour.terminals.push_back(terminals[i].vertex);
                tour.demand += terminals[i].demand;
            }
        }
        result.tours.push_back(tour);
        mask ^= sub;
    }

    return result;
}

}  // namespace tucvrp
