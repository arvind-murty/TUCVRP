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
