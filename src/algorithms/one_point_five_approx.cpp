#include "tucvrp/algorithms/one_point_five_approx.hpp"

#include "tucvrp/preprocessing.hpp"
#include "tucvrp/rng.hpp"

#include <cmath>
#include <stdexcept>

namespace tucvrp {

namespace {

int floor_div(int a, int b) {
    int q = a / b;
    int r = a % b;
    if (r != 0 && ((r > 0) != (b > 0))) {
        --q;
    }
    return q;
}

// Validate epsilon and snap it to the reciprocal grid used by the current scaffold.
OnePointFiveApproxParams validate_params(const OnePointFiveApproxParams& params) {
    if (params.epsilon <= 0.0 || params.epsilon >= 1.0) {
        throw std::invalid_argument("OnePointFiveApproxSolver requires epsilon in (0, 1)");
    }

    OnePointFiveApproxParams normalized = params;
    normalized.epsilon = 1.0 / std::ceil(1.0 / normalized.epsilon);
    return normalized;
}

}  // namespace

SolveResult OnePointFiveApproxSolver::solve(const Instance& instance) {
    return solve(instance, OnePointFiveApproxParams{});
}

// Apply the outer bounded-distance reduction before dispatching to the bounded-distance solver.
SolveResult OnePointFiveApproxSolver::solve(const Instance& instance, const OnePointFiveApproxParams& params) {
    SolveResult result;
    const OnePointFiveApproxParams normalized_params = validate_params(params);
    instance.validate();
    const int one_over_epsilon = static_cast<int>(std::round(1.0 / normalized_params.epsilon));
    const Instance normalized_instance = Preprocessor::make_binary_leaf_tree(instance);
    const RootedTreeData rooted_tree = RootedTreeBuilder::build(normalized_instance);
    const TreeDecomposition decomposition = DecompositionBuilder::make_trivial(rooted_tree);
    (void)decomposition;

    // Partition terminals by randomized logarithmic distance buckets.
    std::vector<Terminal> terminals = normalized_instance.terminals();
    std::vector<double> dists = normalized_instance.distances_from_depot();
    if (terminals.empty()) {
        return result;
    }
    std::sort(terminals.begin(), terminals.end(),
        [&dists](const Terminal& a, const Terminal& b) {
            return dists[a.vertex] < dists[b.vertex];
        });
    const int i_0 = Rng::uniform_int(0, one_over_epsilon - 1);

    int cur_exp_idx = INT_MIN;
    std::size_t i = 0;
    while (i < terminals.size()) {
        // Y_j captures one bounded-distance slab selected by the random offset.
        std::vector<Terminal> y_j;
        int current_bucket = static_cast<int>(std::floor(std::log(dists[terminals[i].vertex]) / std::log(one_over_epsilon)));
        cur_exp_idx = std::max(cur_exp_idx, floor_div(current_bucket - i_0, one_over_epsilon));
        int cur_exp = one_over_epsilon * cur_exp_idx + i_0;
        while (i < terminals.size() &&
               dists[terminals[i].vertex] < std::pow(one_over_epsilon, cur_exp + 1)) {
            y_j.push_back(terminals[i]);
            ++i;
        }
        if (!y_j.empty()) {
            const Instance subinstance = Instance::with_terminals(normalized_instance, y_j);
            const SolveResult subresult = solve_bounded_distance(subinstance, normalized_params);
            result.cost += subresult.cost;
            result.tours.insert(result.tours.end(), subresult.tours.begin(), subresult.tours.end());
        }

        // Z_j aggregates the intervening slabs that will be recursed on separately.
        std::vector<Terminal> z_j;
        for (int j = 1; j < one_over_epsilon; ++j) {
            cur_exp = one_over_epsilon * cur_exp_idx + i_0 + j;
            while (i < terminals.size() &&
                   dists[terminals[i].vertex] < std::pow(one_over_epsilon, cur_exp + 1)) {
                z_j.push_back(terminals[i]);
                ++i;
            }
        }
        if (!z_j.empty()) {
            const Instance subinstance = Instance::with_terminals(normalized_instance, z_j);
            const SolveResult subresult = solve_bounded_distance(subinstance, normalized_params);
            result.cost += subresult.cost;
            result.tours.insert(result.tours.end(), subresult.tours.begin(), subresult.tours.end());
        }

        ++cur_exp_idx;
    }
    
    return result;
}

SolveResult OnePointFiveApproxSolver::solve_bounded_distance(const Instance& instance,
                                                             const OnePointFiveApproxParams& params) {
    const RootedTreeData rooted_tree = RootedTreeBuilder::build(instance);
    const TreeDecomposition decomposition =
        DecompositionBuilder::decompose_bounded_instance(rooted_tree, params.epsilon);
    (void)decomposition;
    // TODO(phase 4): decompose each component further into blocks, clusters, and cells.
    // TODO(phase 4): local-solution machinery entry point.
    // TODO(phase 5): final DP and solution assembly entry point.

    throw std::logic_error(
        "OnePointFiveApproxSolver is scaffolded but the paper algorithm is not implemented yet");
}

}  // namespace tucvrp
