#include "tucvrp/algorithms/one_point_five_approx.hpp"

#include "tucvrp/preprocessing.hpp"
#include "tucvrp/rng.hpp"

#include <cmath>
#include <stdexcept>

namespace tucvrp {

namespace {

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

SolveResult OnePointFiveApproxSolver::solve(const Instance& instance, const OnePointFiveApproxParams& params) {
    SolveResult result;
    const OnePointFiveApproxParams normalized_params = validate_params(params);
    instance.validate();
    const int one_over_epsilon = static_cast<int>(std::round(1.0 / normalized_params.epsilon));
    const Instance normalized_instance = Preprocessor::make_binary_leaf_tree(instance);
    const RootedTreeData rooted_tree = RootedTreeBuilder::build(normalized_instance);
    const TreeDecomposition decomposition = DecompositionBuilder::make_trivial(rooted_tree);
    (void)decomposition;

    // reduce to bounded distance
    std::vector<Terminal> terminals = normalized_instance.terminals();
    const std::unordered_map<int, double> terminal_distances = normalized_instance.terminal_distances();
    std::sort(terminals.begin(), terminals.end(),
        [&terminal_distances](const Terminal& a, const Terminal& b) {
            return terminal_distances.at(a.vertex) < terminal_distances.at(b.vertex);
        });
    const int i_0 = Rng::uniform_int(0, one_over_epsilon - 1);

    const double smallest_terminal_distance = terminal_distances.at(terminals[0].vertex);
    int cur_exp_idx =
        (static_cast<int>(std::floor(std::log(smallest_terminal_distance) / std::log(one_over_epsilon))) - i_0) /
        one_over_epsilon;
    std::size_t i = 0;
    while (i < terminals.size()) {
        std::vector<Terminal> y_j;
        int cur_exp = one_over_epsilon * cur_exp_idx + i_0;
        while (i < terminals.size() &&
               terminal_distances.at(terminals[i].vertex) < std::pow(one_over_epsilon, cur_exp + 1)) {
            y_j.push_back(terminals[i]);
            ++i;
        }
        if (!y_j.empty()) {
            const Instance subinstance(normalized_instance, y_j);
            const SolveResult subresult = solve_bounded_distance(subinstance, normalized_params);
            result.cost += subresult.cost;
            result.tours.insert(result.tours.end(), subresult.tours.begin(), subresult.tours.end());
        }

        std::vector<Terminal> z_j;
        for (int j = 1; j < one_over_epsilon; ++j) {
            cur_exp = one_over_epsilon * cur_exp_idx + i_0 + j;
            while (i < terminals.size() &&
                   terminal_distances.at(terminals[i].vertex) < std::pow(one_over_epsilon, cur_exp + 1)) {
                z_j.push_back(terminals[i]);
                ++i;
            }
        }
        if (!z_j.empty()) {
            const Instance subinstance(normalized_instance, z_j);
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
    (void)instance;
    (void)params;
    // TODO(phase 2): binary-tree and structural normalization entry point.
    // TODO(phase 3): component decomposition entry point.
    // TODO(phase 4): local-solution machinery entry point.
    // TODO(phase 5): final DP and solution assembly entry point.

    throw std::logic_error(
        "OnePointFiveApproxSolver is scaffolded but the paper algorithm is not implemented yet");
}

}  // namespace tucvrp
