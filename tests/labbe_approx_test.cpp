#include "tucvrp/algorithms/labbe_approx.hpp"
#include "tucvrp/exact_solver.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <sstream>
#include <vector>

using tucvrp::ExactSolver;
using tucvrp::Instance;
using tucvrp::LabbeApproxSolver;

namespace {

std::vector<int> all_terminals_in_solution(const tucvrp::SolveResult& result) {
    std::vector<int> terminals;
    for (const auto& tour : result.tours) {
        terminals.insert(terminals.end(), tour.terminals.begin(), tour.terminals.end());
    }
    std::sort(terminals.begin(), terminals.end());
    return terminals;
}

void require_solution_is_feasible(const Instance& instance, const tucvrp::SolveResult& result) {
    REQUIRE(result.cost >= 0.0);

    std::vector<int> covered = all_terminals_in_solution(result);
    std::vector<int> expected;
    for (const auto& terminal : instance.terminals()) {
        expected.push_back(terminal.vertex);
    }
    std::sort(expected.begin(), expected.end());
    REQUIRE(covered == expected);

    double recomputed_cost = 0.0;
    for (const auto& tour : result.tours) {
        REQUIRE(tour.demand <= 1.0 + 1e-9);
        REQUIRE_FALSE(tour.walk.empty());
        REQUIRE(tour.walk.front() == instance.depot());
        REQUIRE(tour.walk.back() == instance.depot());
        REQUIRE(tour.cost == Catch::Approx(instance.tour_cost_for_terminals(tour.terminals)));
        REQUIRE(tour.walk == instance.tour_walk_for_terminals(tour.terminals));
        recomputed_cost += tour.cost;
    }
    REQUIRE(result.cost == Catch::Approx(recomputed_cost));
}

}  // namespace

TEST_CASE("labbe solver returns an empty solution when there are no terminals") {
    std::istringstream input(R"(
3 0
0 1 1
1 2 1
0
)");
    const auto instance = Instance::parse(input);

    const auto result = LabbeApproxSolver::solve(instance);
    REQUIRE(result.cost == Catch::Approx(0.0));
    REQUIRE(result.tours.empty());
}

TEST_CASE("labbe solver reconstructs one route for a single terminal") {
    std::istringstream input(R"(
2 0
0 1 2
1
1 0.2
)");
    const auto instance = Instance::parse(input);

    const auto result = LabbeApproxSolver::solve(instance);
    REQUIRE(result.tours.size() == 1);
    REQUIRE(result.tours[0].terminals == std::vector<int>({1}));
    REQUIRE(result.tours[0].demand == Catch::Approx(0.2));
    REQUIRE(result.tours[0].cost == Catch::Approx(4.0));
    REQUIRE(result.tours[0].walk == std::vector<int>({0, 1, 0}));
    require_solution_is_feasible(instance, result);
}

TEST_CASE("labbe solver merges a leaf aggregate into its parent when the combined demand fits") {
    std::istringstream input(R"(
3 0
0 1 1
1 2 1
2
1 0.2
2 0.3
)");
    const auto instance = Instance::parse(input);

    const auto result = LabbeApproxSolver::solve(instance);
    REQUIRE(result.tours.size() == 1);
    REQUIRE(result.tours[0].terminals == std::vector<int>({1, 2}));
    REQUIRE(result.tours[0].demand == Catch::Approx(0.5));
    REQUIRE(result.tours[0].cost == Catch::Approx(4.0));
    REQUIRE(result.tours[0].walk == std::vector<int>({0, 1, 2, 1, 0}));
    require_solution_is_feasible(instance, result);
}

TEST_CASE("labbe solver emits the leaf aggregate when it is heavier than its parent aggregate") {
    std::istringstream input(R"(
3 0
0 1 1
1 2 1
2
1 0.2
2 0.9
)");
    const auto instance = Instance::parse(input);

    const auto result = LabbeApproxSolver::solve(instance);
    REQUIRE(result.tours.size() == 2);
    REQUIRE(result.tours[0].terminals == std::vector<int>({2}));
    REQUIRE(result.tours[1].terminals == std::vector<int>({1}));
    require_solution_is_feasible(instance, result);
}

TEST_CASE("labbe solver emits the parent aggregate and replaces it by the leaf aggregate otherwise") {
    std::istringstream input(R"(
3 0
0 1 1
1 2 1
2
1 0.7
2 0.4
)");
    const auto instance = Instance::parse(input);

    const auto result = LabbeApproxSolver::solve(instance);
    REQUIRE(result.tours.size() == 2);
    REQUIRE(result.tours[0].terminals == std::vector<int>({1}));
    REQUIRE(result.tours[1].terminals == std::vector<int>({2}));
    require_solution_is_feasible(instance, result);
}

TEST_CASE("labbe solver handles multiple depot-child subtrees independently") {
    std::istringstream input(R"(
3 0
0 1 1
0 2 2
2
1 0.3
2 0.4
)");
    const auto instance = Instance::parse(input);

    const auto result = LabbeApproxSolver::solve(instance);
    REQUIRE(result.tours.size() == 2);
    REQUIRE(all_terminals_in_solution(result) == std::vector<int>({1, 2}));
    REQUIRE(result.cost == Catch::Approx(6.0));
    require_solution_is_feasible(instance, result);
}

TEST_CASE("labbe solver remains feasible on an instance with an internal terminal") {
    std::istringstream input(R"(
4 0
0 1 1
1 2 1
1 3 2
2
1 0.4
3 0.5
)");
    const auto instance = Instance::parse(input);

    const auto result = LabbeApproxSolver::solve(instance);
    require_solution_is_feasible(instance, result);
}

TEST_CASE("labbe solver cost stays within two times exact on small curated instances") {
    const std::vector<std::string> instances = {
        R"(
4 0
0 1 1
1 2 1
1 3 1
2
2 0.4
3 0.5
)",
        R"(
5 0
0 1 1
1 2 2
1 3 1
3 4 1
3
2 0.4
3 0.3
4 0.2
)",
        R"(
6 0
0 1 1
0 2 2
1 3 1
1 4 1
2 5 1
4
3 0.4
4 0.2
5 0.3
2 0.1
)",
    };

    for (const auto& text : instances) {
        std::istringstream input(text);
        const auto instance = Instance::parse(input);
        const auto exact = ExactSolver::solve(instance);
        const auto labbe = LabbeApproxSolver::solve(instance);
        require_solution_is_feasible(instance, labbe);
        REQUIRE(labbe.cost <= 2.0 * exact.cost + 1e-9);
    }
}
