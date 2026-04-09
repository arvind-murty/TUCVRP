#include "tucvrp/instance.hpp"
#include "tucvrp/exact_solver.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <sstream>

using tucvrp::ExactSolver;
using tucvrp::Instance;

TEST_CASE("instance parser validates a simple tree") {
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
    REQUIRE(instance.depot() == 0);
    REQUIRE(instance.vertex_count() == 4);
    REQUIRE(instance.edge_count() == 3);
    REQUIRE(instance.terminal_count() == 2);
    REQUIRE(instance.total_demand() == Catch::Approx(1.0));
    auto terminal_distances = instance.terminal_distances();
    REQUIRE(terminal_distances.size() == 2);
    REQUIRE(terminal_distances.at(2) == Catch::Approx(3.0));
    REQUIRE(terminal_distances.at(3) == Catch::Approx(4.0));
    const auto subtree_terminal_counts = instance.subtree_terminal_counts();
    REQUIRE(subtree_terminal_counts.size() == 4);
    REQUIRE(subtree_terminal_counts[0] == 2);
    REQUIRE(subtree_terminal_counts[1] == 2);
    REQUIRE(subtree_terminal_counts[2] == 1);
    REQUIRE(subtree_terminal_counts[3] == 1);
    REQUIRE(instance.tour_cost_for_terminals({2}) == Catch::Approx(6.0));
    REQUIRE(instance.tour_cost_for_terminals({2, 3}) == Catch::Approx(12.0));
}

TEST_CASE("instance copy constructor preserves tree and terminal state") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 2
1 3 3
2
2 0.4
3 0.6
)");

    const auto original = Instance::parse(input);
    const Instance copy(original);

    REQUIRE(copy.depot() == original.depot());
    REQUIRE(copy.vertex_count() == original.vertex_count());
    REQUIRE(copy.edge_count() == original.edge_count());
    REQUIRE(copy.terminal_count() == original.terminal_count());
    REQUIRE(copy.total_demand() == Catch::Approx(original.total_demand()));
    REQUIRE(copy.terminal_distances().at(2) == Catch::Approx(original.terminal_distances().at(2)));
    REQUIRE(copy.terminal_distances().at(3) == Catch::Approx(original.terminal_distances().at(3)));
    REQUIRE(copy.tour_cost_for_terminals({2, 3}) == Catch::Approx(original.tour_cost_for_terminals({2, 3})));
}

TEST_CASE("instance with_terminals replaces the terminal set on the same tree") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 2
1 3 3
2
2 0.4
3 0.6
)");

    const auto original = Instance::parse(input);
    const Instance subset = Instance::with_terminals(original, {{3, 0.6}});

    REQUIRE(subset.depot() == original.depot());
    REQUIRE(subset.vertex_count() == original.vertex_count());
    REQUIRE(subset.edge_count() == original.edge_count());
    REQUIRE(subset.terminal_count() == 1);
    REQUIRE(subset.total_demand() == Catch::Approx(0.6));
    REQUIRE(subset.is_terminal(3));
    REQUIRE_FALSE(subset.is_terminal(2));
    REQUIRE(subset.demand_of(3) == Catch::Approx(0.6));
    REQUIRE(subset.demand_of(2) == Catch::Approx(0.0));
    REQUIRE(subset.tour_cost_for_terminals({3}) == Catch::Approx(8.0));
}

TEST_CASE("exact solver partitions terminals into feasible tours") {
    std::istringstream input(R"(
5 0
0 1 1
1 2 1
1 3 1
0 4 5
3
2 0.5
3 0.5
4 0.6
)");

    auto instance = Instance::parse(input);
    auto result = ExactSolver::solve(instance);

    REQUIRE(result.tours.size() == 2);
    REQUIRE(result.cost == Catch::Approx(16.0));
}
