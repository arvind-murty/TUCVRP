#include "tucvrp/instance.hpp"
#include "tucvrp/exact_solver.hpp"
#include "tucvrp/preprocessing.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <sstream>

using tucvrp::Instance;
using tucvrp::ExactSolver;
using tucvrp::Preprocessor;

namespace {

bool is_leaf_terminal_tree(const Instance& instance) {
    const auto parent = instance.parent_array();
    for (const auto& terminal : instance.terminals()) {
        for (const auto& edge : instance.neighbors(terminal.vertex)) {
            if (edge.to != parent[terminal.vertex]) {
                return false;
            }
        }
    }
    return true;
}

bool is_binary_rooted_tree(const Instance& instance) {
    const auto parent = instance.parent_array();
    for (int v : instance.vertices()) {
        int child_count = 0;
        for (const auto& edge : instance.neighbors(v)) {
            if (edge.to != parent[v]) {
                ++child_count;
            }
        }
        if (child_count > 2) {
            return false;
        }
    }
    return true;
}

}  // namespace

TEST_CASE("edge load lower bound matches subtree demand accounting") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 2
1 3 3
2
2 0.4
3 0.6
)");

    auto instance = Instance::parse(input);
    REQUIRE(Preprocessor::edge_load_lower_bound(instance) == Catch::Approx(12.0));
}

TEST_CASE("bounded distance statistics follow the paper definition") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 1
1 3 9
2
2 0.4
3 0.4
)");

    auto instance = Instance::parse(input);
    const auto stats = Preprocessor::bounded_distance_stats(instance, 0.5);
    REQUIRE(stats.min_distance == Catch::Approx(2.0));
    REQUIRE(stats.max_distance == Catch::Approx(10.0));
    REQUIRE_FALSE(stats.bounded);
}

TEST_CASE("make_binary_leaf_tree preserves demand and makes branching binary") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(0, 2, 1.0);
    instance.add_edge(0, 3, 1.0);
    instance.add_terminal(1, 0.2);
    instance.add_terminal(2, 0.3);
    instance.add_terminal(3, 0.4);
    instance.validate();

    auto normalized = Preprocessor::make_binary_leaf_tree(instance);

    REQUIRE(normalized.total_demand() == Catch::Approx(instance.total_demand()));
    REQUIRE(is_binary_rooted_tree(normalized));
    REQUIRE(is_leaf_terminal_tree(normalized));
    REQUIRE(normalized.depot() == 0);
}

TEST_CASE("make_binary_leaf_tree pushes internal terminals onto zero-cost leaf children") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(1, 2, 2.0);
    instance.add_edge(1, 3, 3.0);
    instance.add_terminal(1, 0.4);
    instance.add_terminal(2, 0.3);
    instance.validate();

    const auto normalized = Preprocessor::make_binary_leaf_tree(instance);
    const auto parent = normalized.parent_array();

    REQUIRE(normalized.total_demand() == Catch::Approx(instance.total_demand()));
    REQUIRE(normalized.terminal_count() == instance.terminal_count());

    for (const auto& terminal : normalized.terminals()) {
        std::size_t child_count = 0;
        for (const auto& edge : normalized.neighbors(terminal.vertex)) {
            if (edge.to != parent[terminal.vertex]) {
                ++child_count;
            }
        }
        REQUIRE(child_count == 0);
    }

    bool found_zero_cost_promoted_leaf = false;
    for (const auto& terminal : normalized.terminals()) {
        const int p = parent[terminal.vertex];
        if (p < 0) {
            continue;
        }
        for (const auto& edge : normalized.neighbors(terminal.vertex)) {
            if (edge.to == p && edge.weight == Catch::Approx(0.0)) {
                found_zero_cost_promoted_leaf = true;
            }
        }
    }
    REQUIRE(found_zero_cost_promoted_leaf);
}

TEST_CASE("make_binary_leaf_tree preserves the exact optimum on arbitrary trees") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(0, 2, 1.5);
    instance.add_edge(0, 3, 2.0);
    instance.add_edge(1, 4, 1.0);
    instance.add_edge(1, 5, 1.0);
    instance.add_edge(3, 6, 1.0);
    instance.add_terminal(0 + 1, 0.2);
    instance.add_terminal(2, 0.3);
    instance.add_terminal(3, 0.2);
    instance.add_terminal(4, 0.1);
    instance.add_terminal(6, 0.2);
    instance.validate();

    const auto exact = ExactSolver::solve(instance);
    const auto normalized = Preprocessor::make_binary_leaf_tree(instance);

    REQUIRE(normalized.total_demand() == Catch::Approx(instance.total_demand()));
    REQUIRE(is_leaf_terminal_tree(normalized));
    REQUIRE(is_binary_rooted_tree(normalized));
    REQUIRE(ExactSolver::solve(normalized).cost == Catch::Approx(exact.cost));
}
