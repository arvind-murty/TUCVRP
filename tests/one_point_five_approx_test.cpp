#include "tucvrp/algorithms/one_point_five_approx.hpp"
#include "tucvrp/preprocessing.hpp"
#include "tucvrp/rng.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

using tucvrp::Instance;
using tucvrp::LocalSubtourType;
using tucvrp::OnePointFiveApproxParams;
using tucvrp::OnePointFiveApproxSolver;
using tucvrp::Rng;
using tucvrp::RootedTreeBuilder;
using tucvrp::TreeDecomposition;
using tucvrp::Block;
using tucvrp::Cell;
using tucvrp::Cluster;
using tucvrp::Component;
using tucvrp::BoundedDistanceContext;
using tucvrp::LocalPhaseState;
using tucvrp::SubtreeConfigurationTable;
using tucvrp::SubtreePhaseState;
using tucvrp::SubtreeConfigurationValue;

namespace {

double debug_alpha_value(double epsilon) {
    return std::pow(epsilon, 1.0 / epsilon + 1.0);
}

double debug_gamma_zero_value(double epsilon) {
    const double gamma = 12.0 / epsilon;
    return epsilon * debug_alpha_value(epsilon) / gamma;
}

bool debug_same_double(double a, double b) {
    return std::abs(a - b) <= 1e-9;
}

unsigned long long choose_count(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }
    k = std::min(k, n - k);
    unsigned long long result = 1;
    for (int i = 1; i <= k; ++i) {
        result = (result * static_cast<unsigned long long>(n - k + i)) / static_cast<unsigned long long>(i);
    }
    return result;
}

int count_q_c_parts(const BoundedDistanceContext& context, int component_id, double gamma_zero) {
    const Component& component = context.decomposition.components[component_id];
    const auto& rooted_tree = context.rooted_tree;

    int count = 0;

    for (const int v : component.vertices) {
        if (rooted_tree.is_terminal(v) && rooted_tree.demands[v] > gamma_zero) {
            ++count;
        }
    }

    for (const int block_id : component.block_ids) {
        const Block& block = context.decomposition.blocks[block_id];
        for (const int cluster_id : block.cluster_ids) {
            const Cluster& cluster = context.decomposition.clusters[cluster_id];
            for (const int cell_id : cluster.cell_ids) {
                const Cell& cell = context.decomposition.cells[cell_id];
                bool has_small_terminal = false;
                for (const int v : cell.vertices) {
                    if (rooted_tree.is_terminal(v) && rooted_tree.demands[v] <= gamma_zero) {
                        has_small_terminal = true;
                        break;
                    }
                }
                if (has_small_terminal) {
                    ++count;
                }
            }
        }
    }

    return count;
}

std::vector<double> compute_y_c_values(const BoundedDistanceContext& context, int component_id, double gamma_zero, double alpha) {
    const Component& component = context.decomposition.components[component_id];
    const auto& rooted_tree = context.rooted_tree;

    std::vector<double> part_demands;

    for (const int v : component.vertices) {
        if (rooted_tree.is_terminal(v) && rooted_tree.demands[v] > gamma_zero) {
            part_demands.push_back(rooted_tree.demands[v]);
        }
    }

    for (const int block_id : component.block_ids) {
        const Block& block = context.decomposition.blocks[block_id];
        for (const int cluster_id : block.cluster_ids) {
            const Cluster& cluster = context.decomposition.clusters[cluster_id];
            for (const int cell_id : cluster.cell_ids) {
                const Cell& cell = context.decomposition.cells[cell_id];
                double part_demand = 0.0;
                for (const int v : cell.vertices) {
                    if (rooted_tree.is_terminal(v) && rooted_tree.demands[v] <= gamma_zero) {
                        part_demand += rooted_tree.demands[v];
                    }
                }
                if (part_demand > 0.0) {
                    part_demands.push_back(part_demand);
                }
            }
        }
    }

    std::vector<double> y_values{alpha};
    const int part_count = static_cast<int>(part_demands.size());
    const int subset_count = 1 << part_count;
    for (int mask = 1; mask < subset_count; ++mask) {
        double demand = 0.0;
        for (int i = 0; i < part_count; ++i) {
            if (mask & (1 << i)) {
                demand += part_demands[i];
            }
        }
        if (demand > alpha + 1e-9 && demand <= 1.0 + 1e-9) {
            y_values.push_back(std::min(1.0, demand));
        }
    }

    std::sort(y_values.begin(), y_values.end());
    y_values.erase(
        std::unique(y_values.begin(), y_values.end(), [](double a, double b) { return debug_same_double(a, b); }),
        y_values.end());
    return y_values;
}

unsigned long long count_leaf_local_or_subtree_configurations(int q_size, int y_size) {
    unsigned long long total = 0;
    for (int l = 1; l <= q_size; ++l) {
        total += choose_count(y_size + l - 1, l);
    }
    return total;
}

BoundedDistanceContext build_count_only_bounded_context(const Instance& instance, const OnePointFiveApproxParams& params) {
    BoundedDistanceContext context;
    context.instance = tucvrp::Preprocessor::make_binary_leaf_tree(instance);
    context.rooted_tree = RootedTreeBuilder::build(context.instance);
    context.decomposition = tucvrp::DecompositionBuilder::decompose_bounded_instance(context.rooted_tree, params.epsilon);
    context.height_reduced =
        tucvrp::DecompositionBuilder::height_reduce_bounded_components(context.decomposition,
                                                                       context.rooted_tree,
                                                                       params.epsilon);
    return context;
}

}  // namespace

TEST_CASE("paper solver validates epsilon") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 1
1 3 1
2
2 0.4
3 0.5
)");
    const auto instance = Instance::parse(input);

    REQUIRE_THROWS_AS(OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 0.0}),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 1.0}),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = -0.5}),
                      std::invalid_argument);
}

TEST_CASE("paper solver returns a finite bounded-height-reduced value on valid input") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 1
1 3 1
2
2 0.4
3 0.5
)");
    const auto instance = Instance::parse(input);

    const auto result_default = OnePointFiveApproxSolver::solve(instance);
    const auto result = OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 0.25});
    REQUIRE(std::isfinite(result_default.cost));
    REQUIRE(std::isfinite(result.cost));
    REQUIRE_FALSE(result_default.tours.empty());
    REQUIRE_FALSE(result.tours.empty());
    REQUIRE(result_default.cost == Catch::Approx(result_default.tours[0].cost));
    REQUIRE(result.cost == Catch::Approx(result.tours[0].cost));
}

TEST_CASE("paper solver returns an empty solution when there are no terminals") {
    std::istringstream input(R"(
3 0
0 1 1
1 2 1
0
)");
    const auto instance = Instance::parse(input);

    const auto result = OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 0.25});
    REQUIRE(result.cost == Catch::Approx(0.0));
    REQUIRE(result.tours.empty());
}

TEST_CASE("bounded solver reconstructs one tour for a one-terminal bounded instance") {
    std::istringstream input(R"(
2 0
0 1 2
1
1 0.2
)");
    const auto instance = Instance::parse(input);

    const auto result = OnePointFiveApproxSolver::solve_bounded_distance(
        instance,
        OnePointFiveApproxParams{.epsilon = 0.9});

    REQUIRE(result.cost == Catch::Approx(4.0));
    REQUIRE(result.tours.size() == 1);
    REQUIRE(result.tours[0].terminals == std::vector<int>{1});
    REQUIRE(result.tours[0].walk == std::vector<int>({0, 1, 0}));
    REQUIRE(result.tours[0].demand == Catch::Approx(0.2));
    REQUIRE(result.tours[0].cost == Catch::Approx(4.0));
}

TEST_CASE("paper solver reconstructs a feasible partition of all terminals") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 1
1 3 1
2
2 0.4
3 0.5
)");
    const auto instance = Instance::parse(input);

    const auto result = OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 0.9});

    REQUIRE(result.cost == Catch::Approx(6.0));
    REQUIRE(result.tours.size() == 1);
    REQUIRE(result.tours[0].terminals == std::vector<int>({2, 3}));
    REQUIRE(result.tours[0].walk == std::vector<int>({0, 1, 2, 1, 3, 1, 0}));
    REQUIRE(result.tours[0].demand == Catch::Approx(0.9));
    REQUIRE(result.tours[0].cost == Catch::Approx(6.0));
}

TEST_CASE("paper solver projects normalized terminals back to the original instance") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 2
2 3 3
2
2 0.4
3 0.5
)");
    const auto instance = Instance::parse(input);

    const auto result = OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 0.9});

    REQUIRE(std::isfinite(result.cost));
    std::vector<int> covered_terminals;
    for (const auto& tour : result.tours) {
        covered_terminals.insert(covered_terminals.end(), tour.terminals.begin(), tour.terminals.end());
        for (const int vertex : tour.walk) {
            REQUIRE(vertex >= 0);
            REQUIRE(vertex < instance.vertex_count());
        }
    }
    std::sort(covered_terminals.begin(), covered_terminals.end());
    REQUIRE(covered_terminals == std::vector<int>({2, 3}));
}

TEST_CASE("paper solver handles terminals on logarithmic bucket boundaries") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 3
2 3 12
3
1 0.2
2 0.3
3 0.4
)");
    const auto instance = Instance::parse(input);

    Rng::seed(11);
    const auto result = OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 0.25});
    REQUIRE(std::isfinite(result.cost));
}

TEST_CASE("paper solver handles terminals at distances below one before bounded solver dispatch") {
    std::istringstream input(R"(
4 0
0 1 0.2
1 2 0.1
1 3 0.15
2
2 0.4
3 0.5
)");
    const auto instance = Instance::parse(input);

    Rng::seed(7);
    const auto result = OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 0.25});
    REQUIRE(std::isfinite(result.cost));
}

TEST_CASE("README sanity-check tree S1 has the documented local and subtree counts") {
    std::istringstream input(R"(
2 0
0 1 2
1
1 0.2
)");
    const auto instance = Instance::parse(input);

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {0, 1},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
        .cluster_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
        .cell_ids = {0},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
    });

    const auto local_table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});
    const auto subtree_table = OnePointFiveApproxSolver::compute_component_root_subtree_configurations(
        context,
        0,
        local_table,
        nullptr,
        OnePointFiveApproxParams{.epsilon = 0.25});

    REQUIRE(local_table.parts.size() == 1);
    REQUIRE(local_table.y_values.size() == 2);
    REQUIRE(local_table.values.size() == 2);
    REQUIRE(subtree_table.values.size() == 2);
}

TEST_CASE("README sanity-check tree S2 has the documented local and subtree counts") {
    std::istringstream input(R"(
3 0
0 1 1
0 2 1
2
1 0.2
2 0.3
)");
    const auto instance = Instance::parse(input);

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 2,
        .is_leaf = true,
        .is_big = true,
        .vertices = {0, 1, 2},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.5,
        .vertices = {0, 1, 2},
        .cluster_ids = {0, 1},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
        .cell_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 1,
        .block_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.3,
        .vertices = {0, 2},
        .cell_ids = {1},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 1,
        .cluster_id = 1,
        .root = 0,
        .exit = -1,
        .demand = 0.3,
        .vertices = {0, 2},
    });

    const auto local_table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});
    const auto subtree_table = OnePointFiveApproxSolver::compute_component_root_subtree_configurations(
        context,
        0,
        local_table,
        nullptr,
        OnePointFiveApproxParams{.epsilon = 0.25});

    REQUIRE(local_table.parts.size() == 2);
    REQUIRE(local_table.y_values.size() == 4);
    REQUIRE(local_table.values.size() == 14);
    REQUIRE(subtree_table.values.size() == 14);
}

TEST_CASE("README sanity-check tree S3 has the documented local and subtree counts") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 1
1 3 1
1
3 0.2
)");
    const auto instance = Instance::parse(input);

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {0, 1, 2, 3},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
        .cluster_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
        .cell_ids = {0},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
    });

    const auto local_table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});
    const auto subtree_table = OnePointFiveApproxSolver::compute_component_root_subtree_configurations(
        context,
        0,
        local_table,
        nullptr,
        OnePointFiveApproxParams{.epsilon = 0.25});

    REQUIRE(local_table.parts.size() == 1);
    REQUIRE(local_table.y_values.size() == 2);
    REQUIRE(local_table.values.size() == 2);
    REQUIRE(subtree_table.values.size() == 2);
}

TEST_CASE("README timeout tree T1 has the documented Q_c, Y_c, and count-only leaf state sizes") {
    std::istringstream input(R"(
10 0
0 1 3
0 2 7
1 3 6
2 4 3
4 5 8
0 6 8
1 7 3
7 8 2
0 9 3
5
1 0.177
2 0.207
3 0.410
8 0.401
9 0.205
)");
    const auto instance = Instance::parse(input);
    const auto params = OnePointFiveApproxParams{.epsilon = 0.25};
    const auto context = build_count_only_bounded_context(instance, params);

    REQUIRE(context.decomposition.components.size() == 1);
    REQUIRE(context.decomposition.components[0].is_leaf);
    REQUIRE(context.decomposition.components[0].exit == -1);

    const double alpha = debug_alpha_value(params.epsilon);
    const double gamma_zero = debug_gamma_zero_value(params.epsilon);
    const int q_size = count_q_c_parts(context, 0, gamma_zero);
    const auto y_values = compute_y_c_values(context, 0, gamma_zero, alpha);
    const auto count = count_leaf_local_or_subtree_configurations(q_size, static_cast<int>(y_values.size()));

    REQUIRE(q_size == 5);
    REQUIRE(y_values.size() == 26);
    REQUIRE(count == 169910ULL);
}

TEST_CASE("README timeout tree T2 has the documented Q_c, Y_c, and count-only leaf state sizes") {
    std::istringstream input(R"(
10 0
0 1 1
1 2 6
0 3 6
1 4 9
2 5 6
4 6 7
3 7 6
4 8 3
5 9 9
5
1 0.215
4 0.539
5 0.175
7 0.227
8 0.225
)");
    const auto instance = Instance::parse(input);
    const auto params = OnePointFiveApproxParams{.epsilon = 0.25};
    const auto context = build_count_only_bounded_context(instance, params);

    REQUIRE(context.decomposition.components.size() == 1);
    REQUIRE(context.decomposition.components[0].is_leaf);
    REQUIRE(context.decomposition.components[0].exit == -1);

    const double alpha = debug_alpha_value(params.epsilon);
    const double gamma_zero = debug_gamma_zero_value(params.epsilon);
    const int q_size = count_q_c_parts(context, 0, gamma_zero);
    const auto y_values = compute_y_c_values(context, 0, gamma_zero, alpha);
    const auto count = count_leaf_local_or_subtree_configurations(q_size, static_cast<int>(y_values.size()));

    REQUIRE(q_size == 5);
    REQUIRE(y_values.size() == 27);
    REQUIRE(count == 201375ULL);
}

TEST_CASE("README timeout tree T3 has the documented Q_c, Y_c, and count-only leaf state sizes") {
    std::istringstream input(R"(
10 0
0 1 2
0 2 8
0 3 5
3 4 4
4 5 3
1 6 9
4 7 8
6 8 3
7 9 2
5
2 0.417
3 0.127
4 0.238
5 0.275
8 0.524
)");
    const auto instance = Instance::parse(input);
    const auto params = OnePointFiveApproxParams{.epsilon = 0.25};
    const auto context = build_count_only_bounded_context(instance, params);

    REQUIRE(context.decomposition.components.size() == 1);
    REQUIRE(context.decomposition.components[0].is_leaf);
    REQUIRE(context.decomposition.components[0].exit == -1);

    const double alpha = debug_alpha_value(params.epsilon);
    const double gamma_zero = debug_gamma_zero_value(params.epsilon);
    const int q_size = count_q_c_parts(context, 0, gamma_zero);
    const auto y_values = compute_y_c_values(context, 0, gamma_zero, alpha);
    const auto count = count_leaf_local_or_subtree_configurations(q_size, static_cast<int>(y_values.size()));

    REQUIRE(q_size == 5);
    REQUIRE(y_values.size() == 22);
    REQUIRE(count == 80729ULL);
}

TEST_CASE("algorithm 2 computes local configurations for a leaf component") {
    std::istringstream input(R"(
2 0
0 1 2
1
1 0.2
)");
    const auto instance = Instance::parse(input);

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {0, 1},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
        .cluster_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
        .cell_ids = {0},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
    });

    const auto table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});

    REQUIRE(table.component_id == 0);
    REQUIRE(table.parts.size() == 1);
    REQUIRE(table.parts[0].terminals == std::vector<int>{1});
    REQUIRE(table.parts[0].demand == Catch::Approx(0.2));
    REQUIRE(table.y_values.size() == 2);
    REQUIRE(table.y_values[0] < 0.01);
    REQUIRE(table.y_values[1] == Catch::Approx(0.2));

    bool found_feasible = false;
    bool found_infeasible = false;
    for (const auto& value : table.values) {
        if (value.configuration.entries.size() != 1) {
            continue;
        }
        const auto& entry = value.configuration.entries[0];
        REQUIRE(entry.type == LocalSubtourType::Ending);
        if (entry.demand_bound == Catch::Approx(0.2)) {
            found_feasible = true;
            REQUIRE(value.cost == Catch::Approx(4.0));
        } else if (entry.demand_bound < 0.01) {
            found_infeasible = true;
            REQUIRE(std::isinf(value.cost));
        }
    }
    REQUIRE(found_feasible);
    REQUIRE(found_infeasible);
}

TEST_CASE("algorithm 2 distinguishes ending and passing local subtours") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(1, 2, 1.0);
    instance.add_edge(1, 3, 1.0);
    instance.add_terminal(3, 0.2);
    instance.validate();

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = 2,
        .terminal_count = 1,
        .is_leaf = false,
        .is_big = true,
        .vertices = {0, 1, 2, 3},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
        .cluster_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
        .cell_ids = {0},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
    });

    const auto table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});

    bool found_ending = false;
    bool found_passing = false;
    for (const auto& value : table.values) {
        if (value.configuration.entries.size() != 1) {
            continue;
        }
        const auto& entry = value.configuration.entries[0];
        if (entry.demand_bound != Catch::Approx(0.2)) {
            continue;
        }
        if (entry.type == LocalSubtourType::Ending) {
            found_ending = true;
            REQUIRE(value.cost == Catch::Approx(4.0));
        } else {
            found_passing = true;
            REQUIRE(value.cost == Catch::Approx(6.0));
        }
    }
    REQUIRE(found_ending);
    REQUIRE(found_passing);
}

TEST_CASE("algorithm 2 combines multiple parts into one or two ending subtours") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(1, 2, 1.0);
    instance.add_edge(1, 3, 1.0);
    instance.add_terminal(2, 0.2);
    instance.add_terminal(3, 0.2);
    instance.validate();

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 2,
        .is_leaf = true,
        .is_big = true,
        .vertices = {0, 1, 2, 3},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.4,
        .vertices = {0, 1, 2, 3},
        .cluster_ids = {0, 1},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1, 2},
        .cell_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 1,
        .block_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1, 3},
        .cell_ids = {1},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1, 2},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 1,
        .cluster_id = 1,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1, 3},
    });

    const auto table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});

    bool found_one_subtour = false;
    bool found_two_subtours = false;
    for (const auto& value : table.values) {
        if (value.configuration.entries.size() == 1 &&
            value.configuration.entries[0].type == LocalSubtourType::Ending &&
            value.configuration.entries[0].demand_bound == Catch::Approx(0.4)) {
            found_one_subtour = true;
            REQUIRE(value.cost == Catch::Approx(6.0));
        }
        if (value.configuration.entries.size() == 2 &&
            value.configuration.entries[0].type == LocalSubtourType::Ending &&
            value.configuration.entries[1].type == LocalSubtourType::Ending &&
            value.configuration.entries[0].demand_bound == Catch::Approx(0.2) &&
            value.configuration.entries[1].demand_bound == Catch::Approx(0.2)) {
            found_two_subtours = true;
            REQUIRE(value.cost == Catch::Approx(8.0));
        }
    }
    REQUIRE(found_one_subtour);
    REQUIRE(found_two_subtours);
}

TEST_CASE("algorithm 5 turns leaf-component local configurations into subtree configurations") {
    std::istringstream input(R"(
2 0
0 1 2
1
1 0.2
)");
    const auto instance = Instance::parse(input);

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {0, 1},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
        .cluster_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
        .cell_ids = {0},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
    });

    const auto local_table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});
    const auto subtree_table = OnePointFiveApproxSolver::compute_component_root_subtree_configurations(
        context,
        0,
        local_table,
        nullptr,
        OnePointFiveApproxParams{.epsilon = 0.25});

    REQUIRE(subtree_table.vertex == 0);
    bool found = false;
    for (const auto& value : subtree_table.values) {
        if (value.configuration.entries.size() == 1 &&
            value.configuration.entries[0].demand == Catch::Approx(0.2) &&
            value.configuration.entries[0].multiplicity == 1) {
            found = true;
            REQUIRE(value.cost == Catch::Approx(4.0));
        }
    }
    REQUIRE(found);
}

TEST_CASE("algorithm 5 produces the correct number of subtree configurations for a simple leaf component") {
    std::istringstream input(R"(
3 0
0 1 1
0 2 1
2
1 0.2
2 0.3
)");
    const auto instance = Instance::parse(input);

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 2,
        .is_leaf = true,
        .is_big = true,
        .vertices = {0, 1, 2},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.5,
        .vertices = {0, 1, 2},
        .cluster_ids = {0, 1},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
        .cell_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 1,
        .block_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.3,
        .vertices = {0, 2},
        .cell_ids = {1},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = -1,
        .demand = 0.2,
        .vertices = {0, 1},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 1,
        .cluster_id = 1,
        .root = 0,
        .exit = -1,
        .demand = 0.3,
        .vertices = {0, 2},
    });

    const auto local_table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});
    const auto subtree_table = OnePointFiveApproxSolver::compute_component_root_subtree_configurations(
        context,
        0,
        local_table,
        nullptr,
        OnePointFiveApproxParams{.epsilon = 0.25});

    // For this leaf component, the subtree configurations are exactly the multisets of size 1 or 2
    // drawn from Y_c = {alpha, 0.2, 0.3, 0.5}, so the theoretical count is C(4,1) + C(5,2) = 14.
    REQUIRE(local_table.parts.size() == 2);
    REQUIRE(local_table.y_values.size() == 4);
    REQUIRE(subtree_table.values.size() == 14);
}

TEST_CASE("algorithm 5 combines exit subtree configurations with internal-component local configurations") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(1, 2, 1.0);
    instance.add_edge(1, 3, 1.0);
    instance.add_terminal(3, 0.2);
    instance.validate();

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = 2,
        .terminal_count = 1,
        .is_leaf = false,
        .is_big = true,
        .vertices = {0, 1, 2, 3},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
        .cluster_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
        .cell_ids = {0},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
    });

    const auto local_table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});
    SubtreeConfigurationTable exit_table;
    exit_table.vertex = 2;
    exit_table.values.push_back(SubtreeConfigurationValue{
        .configuration = tucvrp::SubtreeConfiguration{
            .entries = {
                tucvrp::SubtreeConfigurationEntry{.demand = 0.3, .multiplicity = 1},
            },
        },
        .cost = 5.0,
    });

    const auto subtree_table = OnePointFiveApproxSolver::compute_component_root_subtree_configurations(
        context,
        0,
        local_table,
        &exit_table,
        OnePointFiveApproxParams{.epsilon = 0.25});

    bool found = false;
    for (const auto& value : subtree_table.values) {
        if (value.configuration.entries.size() == 1 &&
            value.configuration.entries[0].demand == Catch::Approx(0.5) &&
            value.configuration.entries[0].multiplicity == 1) {
            found = true;
            REQUIRE(value.cost == Catch::Approx(11.0));
        }
    }
    REQUIRE(found);
}

TEST_CASE("algorithm 5 produces the correct number of subtree configurations for a one-part internal component") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(1, 2, 1.0);
    instance.add_edge(1, 3, 1.0);
    instance.add_terminal(3, 0.2);
    instance.validate();

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = 2,
        .terminal_count = 1,
        .is_leaf = false,
        .is_big = true,
        .vertices = {0, 1, 2, 3},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
        .cluster_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
        .cell_ids = {0},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
    });

    const auto local_table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});
    SubtreeConfigurationTable exit_table;
    exit_table.vertex = 2;
    exit_table.values.push_back(SubtreeConfigurationValue{
        .configuration = tucvrp::SubtreeConfiguration{
            .entries = {
                tucvrp::SubtreeConfigurationEntry{.demand = 0.3, .multiplicity = 1},
            },
        },
        .cost = 5.0,
    });

    const auto subtree_table = OnePointFiveApproxSolver::compute_component_root_subtree_configurations(
        context,
        0,
        local_table,
        &exit_table,
        OnePointFiveApproxParams{.epsilon = 0.25});

    // The local configurations are [(alpha,E)], [(alpha,P)], [(0.2,E)], [(0.2,P)].
    // With one exit subtour of demand 0.3, these induce four distinct subtree configurations:
    // {0.3,alpha}, {alpha+0.3}, {0.3,0.2}, {0.5}.
    REQUIRE(local_table.parts.size() == 1);
    REQUIRE(local_table.y_values.size() == 2);
    REQUIRE(subtree_table.values.size() == 4);
}

TEST_CASE("algorithm 5 rejects missing exit subtree configurations for internal components") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(1, 2, 1.0);
    instance.add_edge(1, 3, 1.0);
    instance.add_terminal(3, 0.2);
    instance.validate();

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = 2,
        .terminal_count = 1,
        .is_leaf = false,
        .is_big = true,
        .vertices = {0, 1, 2, 3},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
        .cluster_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
        .cell_ids = {0},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 0,
        .exit = 2,
        .demand = 0.2,
        .vertices = {0, 1, 2, 3},
    });

    const auto local_table =
        OnePointFiveApproxSolver::compute_local_configurations(context, 0, OnePointFiveApproxParams{.epsilon = 0.25});

    REQUIRE_THROWS_AS(
        OnePointFiveApproxSolver::compute_component_root_subtree_configurations(
            context,
            0,
            local_table,
            nullptr,
            OnePointFiveApproxParams{.epsilon = 0.25}),
        std::invalid_argument);
}

TEST_CASE("algorithm 6 propagates one child table to a critical vertex") {
    std::istringstream input(R"(
2 0
0 1 2
1
1 0.2
)");
    const auto instance = Instance::parse(input);

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 0,
        .is_leaf = false,
        .is_big = false,
        .vertices = {0},
        .block_ids = {},
    });
    context.decomposition.components.push_back(Component{
        .id = 1,
        .root = 1,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {1},
        .block_ids = {},
    });
    context.height_reduced.critical_vertex_by_component = {0, 0};

    SubtreeConfigurationTable child_table;
    child_table.vertex = 1;
    child_table.values.push_back(SubtreeConfigurationValue{
        .configuration = tucvrp::SubtreeConfiguration{
            .entries = {
                tucvrp::SubtreeConfigurationEntry{.demand = 0.2, .multiplicity = 1},
            },
        },
        .cost = 5.0,
    });

    const auto table = OnePointFiveApproxSolver::compute_critical_vertex_subtree_configurations(
        context,
        0,
        {child_table},
        OnePointFiveApproxParams{.epsilon = 0.9});

    REQUIRE(table.vertex == 0);
    bool found = false;
    for (const auto& value : table.values) {
        if (value.configuration.entries.size() == 1 &&
            value.configuration.entries[0].demand == Catch::Approx(0.2) &&
            value.configuration.entries[0].multiplicity == 1) {
            found = true;
            REQUIRE(value.cost == Catch::Approx(9.0));
        }
    }
    REQUIRE(found);
}

TEST_CASE("algorithm 6 preserves the subtree-configuration count with one child") {
    std::istringstream input(R"(
2 0
0 1 2
1
1 0.2
)");
    const auto instance = Instance::parse(input);

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 0,
        .is_leaf = false,
        .is_big = false,
        .vertices = {0},
        .block_ids = {},
    });
    context.decomposition.components.push_back(Component{
        .id = 1,
        .root = 1,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {1},
        .block_ids = {},
    });
    context.height_reduced.critical_vertex_by_component = {0, 0};

    SubtreeConfigurationTable child_table;
    child_table.vertex = 1;
    child_table.values.push_back(SubtreeConfigurationValue{
        .configuration = tucvrp::SubtreeConfiguration{
            .entries = {
                tucvrp::SubtreeConfigurationEntry{.demand = 0.2, .multiplicity = 1},
            },
        },
        .cost = 5.0,
    });

    const auto table = OnePointFiveApproxSolver::compute_critical_vertex_subtree_configurations(
        context,
        0,
        {child_table},
        OnePointFiveApproxParams{.epsilon = 0.9});

    // Algorithm 6 guesses a set X and rounds child demands upward to values in X.
    // With a single child subtour of demand 0.2, the two distinct rounded outcomes are:
    // - {0.2} when X contains 0.2
    // - {1.0} when X is only the mandatory top sentinel {1.0}
    REQUIRE(table.values.size() == 2);
}

TEST_CASE("algorithm 6 combines subtree configurations from two child roots") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(0, 2, 2.0);
    instance.add_terminal(1, 0.2);
    instance.add_terminal(2, 0.3);
    instance.validate();

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 0,
        .is_leaf = false,
        .is_big = false,
        .vertices = {0},
        .block_ids = {},
    });
    context.decomposition.components.push_back(Component{
        .id = 1,
        .root = 1,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {1},
        .block_ids = {},
    });
    context.decomposition.components.push_back(Component{
        .id = 2,
        .root = 2,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {2},
        .block_ids = {},
    });
    context.height_reduced.critical_vertex_by_component = {0, 0, 0};

    SubtreeConfigurationTable first_child;
    first_child.vertex = 1;
    first_child.values.push_back(SubtreeConfigurationValue{
        .configuration = tucvrp::SubtreeConfiguration{
            .entries = {
                tucvrp::SubtreeConfigurationEntry{.demand = 0.2, .multiplicity = 1},
            },
        },
        .cost = 5.0,
    });

    SubtreeConfigurationTable second_child;
    second_child.vertex = 2;
    second_child.values.push_back(SubtreeConfigurationValue{
        .configuration = tucvrp::SubtreeConfiguration{
            .entries = {
                tucvrp::SubtreeConfigurationEntry{.demand = 0.3, .multiplicity = 1},
            },
        },
        .cost = 7.0,
    });

    const auto table = OnePointFiveApproxSolver::compute_critical_vertex_subtree_configurations(
        context,
        0,
        {first_child, second_child},
        OnePointFiveApproxParams{.epsilon = 0.9});

    bool found_split = false;
    bool found_combined = false;
    for (const auto& value : table.values) {
        if (value.configuration.entries.size() == 2 &&
            value.configuration.entries[0].demand == Catch::Approx(0.2) &&
            value.configuration.entries[0].multiplicity == 1 &&
            value.configuration.entries[1].demand == Catch::Approx(0.3) &&
            value.configuration.entries[1].multiplicity == 1) {
            found_split = true;
            REQUIRE(value.cost == Catch::Approx(18.0));
        }
        if (value.configuration.entries.size() == 1 &&
            value.configuration.entries[0].demand == Catch::Approx(0.5) &&
            value.configuration.entries[0].multiplicity == 1) {
            found_combined = true;
            REQUIRE(value.cost == Catch::Approx(18.0));
        }
    }
    REQUIRE(found_split);
    REQUIRE(found_combined);
}

TEST_CASE("algorithm 6 produces the correct number of subtree configurations for two one-tour children") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(0, 2, 2.0);
    instance.add_terminal(1, 0.2);
    instance.add_terminal(2, 0.3);
    instance.validate();

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 0,
        .is_leaf = false,
        .is_big = false,
        .vertices = {0},
        .block_ids = {},
    });
    context.decomposition.components.push_back(Component{
        .id = 1,
        .root = 1,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {1},
        .block_ids = {},
    });
    context.decomposition.components.push_back(Component{
        .id = 2,
        .root = 2,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {2},
        .block_ids = {},
    });
    context.height_reduced.critical_vertex_by_component = {0, 0, 0};

    SubtreeConfigurationTable first_child;
    first_child.vertex = 1;
    first_child.values.push_back(SubtreeConfigurationValue{
        .configuration = tucvrp::SubtreeConfiguration{
            .entries = {
                tucvrp::SubtreeConfigurationEntry{.demand = 0.2, .multiplicity = 1},
            },
        },
        .cost = 5.0,
    });

    SubtreeConfigurationTable second_child;
    second_child.vertex = 2;
    second_child.values.push_back(SubtreeConfigurationValue{
        .configuration = tucvrp::SubtreeConfiguration{
            .entries = {
                tucvrp::SubtreeConfigurationEntry{.demand = 0.3, .multiplicity = 1},
            },
        },
        .cost = 7.0,
    });

    const auto table = OnePointFiveApproxSolver::compute_critical_vertex_subtree_configurations(
        context,
        0,
        {first_child, second_child},
        OnePointFiveApproxParams{.epsilon = 0.9});

    // Algorithm 6 works over rounded subtree configurations rather than the exact child demands.
    // With child demands 0.2 and 0.3, the guessed-X rounding yields the following distinct
    // configurations at z:
    // - {0.2, 0.3}
    // - {0.2, 1.0}
    // - {0.3, 1.0}
    // - {1.0, 1.0}
    // - {0.5}
    // - {1.0}
    REQUIRE(table.values.size() == 6);
}

TEST_CASE("algorithm 6 rejects child tables that do not match the critical-vertex children") {
    std::istringstream input(R"(
2 0
0 1 2
1
1 0.2
)");
    const auto instance = Instance::parse(input);

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 0,
        .is_leaf = false,
        .is_big = false,
        .vertices = {0},
        .block_ids = {},
    });
    context.decomposition.components.push_back(Component{
        .id = 1,
        .root = 1,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {1},
        .block_ids = {},
    });
    context.height_reduced.critical_vertex_by_component = {0, 0};

    SubtreeConfigurationTable child_table;
    child_table.vertex = 7;

    REQUIRE_THROWS_AS(
        OnePointFiveApproxSolver::compute_critical_vertex_subtree_configurations(
            context,
            0,
            {child_table},
            OnePointFiveApproxParams{.epsilon = 0.9}),
        std::invalid_argument);
}

TEST_CASE("subtree phase runs bottom-up on a trivial root-critical-vertex instance") {
    std::istringstream input(R"(
2 0
0 1 2
1
1 0.2
)");
    const auto instance = Instance::parse(input);

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 0,
        .is_leaf = false,
        .is_big = false,
        .vertices = {0},
        .block_ids = {},
    });
    context.decomposition.components.push_back(Component{
        .id = 1,
        .root = 1,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {1},
        .block_ids = {0},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 1,
        .root = 1,
        .exit = -1,
        .demand = 0.2,
        .vertices = {1},
        .cluster_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 1,
        .exit = -1,
        .demand = 0.2,
        .vertices = {1},
        .cell_ids = {0},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 1,
        .exit = -1,
        .demand = 0.2,
        .vertices = {1},
    });
    context.height_reduced.original_parent_component = {-1, 0};
    context.height_reduced.critical_vertex_by_component = {0, 0};
    context.height_reduced.groups.push_back(
        tucvrp::HeightReducedComponentGroup{.id = 0, .class_index = 1, .critical_vertex = 0, .component_ids = {0, 1}});

    LocalPhaseState local_phase;
    local_phase.tables.resize(2);
    local_phase.tables[0] = tucvrp::LocalConfigurationTable{
        .component_id = 0,
        .alpha = 0.0,
        .parts = {},
        .y_values = {},
        .values = {tucvrp::LocalConfigurationValue{.configuration = {}, .cost = 0.0}},
    };
    local_phase.tables[1] =
        OnePointFiveApproxSolver::compute_local_configurations(context, 1, OnePointFiveApproxParams{.epsilon = 0.9});

    const auto subtree_phase =
        OnePointFiveApproxSolver::compute_subtree_phase(context, local_phase, OnePointFiveApproxParams{.epsilon = 0.9});

    REQUIRE(subtree_phase.component_root_tables[1].vertex == 1);
    REQUIRE(subtree_phase.critical_vertex_tables[0].vertex == 0);
    REQUIRE(OnePointFiveApproxSolver::bounded_height_reduced_opt_value(subtree_phase, 0) ==
            Catch::Approx(4.0));
}

TEST_CASE("subtree phase feeds a critical-vertex table into an internal component root table") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(1, 2, 1.0);
    instance.add_edge(1, 3, 1.0);
    instance.add_terminal(3, 0.2);
    instance.validate();

    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(instance);
    context.decomposition.depot = 0;
    context.decomposition.components.push_back(Component{
        .id = 0,
        .root = 0,
        .exit = -1,
        .terminal_count = 0,
        .is_leaf = false,
        .is_big = false,
        .vertices = {0},
        .block_ids = {},
    });
    context.decomposition.components.push_back(Component{
        .id = 1,
        .root = 1,
        .exit = 2,
        .terminal_count = 1,
        .is_leaf = false,
        .is_big = true,
        .vertices = {1, 2, 3},
        .block_ids = {0},
    });
    context.decomposition.components.push_back(Component{
        .id = 2,
        .root = 2,
        .exit = -1,
        .terminal_count = 0,
        .is_leaf = false,
        .is_big = false,
        .vertices = {2},
        .block_ids = {},
    });
    context.decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 1,
        .root = 1,
        .exit = 2,
        .demand = 0.2,
        .vertices = {1, 2, 3},
        .cluster_ids = {0},
    });
    context.decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = 1,
        .exit = 2,
        .demand = 0.2,
        .vertices = {1, 2, 3},
        .cell_ids = {0},
    });
    context.decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = 1,
        .exit = 2,
        .demand = 0.2,
        .vertices = {1, 2, 3},
    });
    context.height_reduced.original_parent_component = {-1, 0, 1};
    context.height_reduced.critical_vertex_by_component = {0, 0, 2};
    context.height_reduced.groups.push_back(
        tucvrp::HeightReducedComponentGroup{.id = 0, .class_index = 1, .critical_vertex = 2, .component_ids = {2}});
    context.height_reduced.groups.push_back(
        tucvrp::HeightReducedComponentGroup{.id = 1, .class_index = 1, .critical_vertex = 0, .component_ids = {0, 1}});

    LocalPhaseState local_phase;
    local_phase.tables.resize(3);
    local_phase.tables[0] = tucvrp::LocalConfigurationTable{
        .component_id = 0,
        .alpha = 0.0,
        .parts = {},
        .y_values = {},
        .values = {tucvrp::LocalConfigurationValue{.configuration = {}, .cost = 0.0}},
    };
    local_phase.tables[1] =
        OnePointFiveApproxSolver::compute_local_configurations(context, 1, OnePointFiveApproxParams{.epsilon = 0.9});
    local_phase.tables[2] = tucvrp::LocalConfigurationTable{
        .component_id = 2,
        .alpha = 0.0,
        .parts = {},
        .y_values = {},
        .values = {tucvrp::LocalConfigurationValue{.configuration = {}, .cost = 0.0}},
    };

    const auto subtree_phase =
        OnePointFiveApproxSolver::compute_subtree_phase(context, local_phase, OnePointFiveApproxParams{.epsilon = 0.9});

    REQUIRE(subtree_phase.critical_vertex_tables[2].vertex == 2);
    REQUIRE(subtree_phase.component_root_tables[1].vertex == 1);

    bool found = false;
    for (const auto& value : subtree_phase.component_root_tables[1].values) {
        if (value.configuration.entries.size() == 1 &&
            value.configuration.entries[0].demand == Catch::Approx(std::pow(0.9, 1.0 / 0.9 + 1.0)) &&
            value.configuration.entries[0].multiplicity == 1) {
            found = true;
            REQUIRE(value.cost == Catch::Approx(2.0));
        }
    }
    REQUIRE(found);
}

TEST_CASE("bounded height-reduced value picks the empty root configuration") {
    SubtreePhaseState subtree_phase;
    subtree_phase.critical_vertex_tables.resize(1);
    subtree_phase.critical_vertex_tables[0].vertex = 0;
    subtree_phase.critical_vertex_tables[0].values = {
        SubtreeConfigurationValue{
            .configuration = tucvrp::SubtreeConfiguration{},
            .cost = 7.0,
        },
        SubtreeConfigurationValue{
            .configuration =
                tucvrp::SubtreeConfiguration{
                    .entries = {tucvrp::SubtreeConfigurationEntry{.demand = 0.2, .multiplicity = 1}},
                },
            .cost = 3.0,
        },
    };

    REQUIRE(OnePointFiveApproxSolver::bounded_height_reduced_opt_value(subtree_phase, 0) ==
            Catch::Approx(3.0));
}
