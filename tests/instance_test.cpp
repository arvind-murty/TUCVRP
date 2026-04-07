#include "tucvrp/instance.hpp"
#include "tucvrp/solver.hpp"

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
    REQUIRE(instance.tour_cost_for_terminals({2}) == Catch::Approx(6.0));
    REQUIRE(instance.tour_cost_for_terminals({2, 3}) == Catch::Approx(12.0));
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
