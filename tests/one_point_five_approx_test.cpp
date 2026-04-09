#include "tucvrp/algorithms/one_point_five_approx.hpp"

#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include <stdexcept>

using tucvrp::Instance;
using tucvrp::OnePointFiveApproxParams;
using tucvrp::OnePointFiveApproxSolver;

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

TEST_CASE("paper solver stub throws not implemented on valid input") {
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

    REQUIRE_THROWS_WITH(OnePointFiveApproxSolver::solve(instance),
                        Catch::Matchers::ContainsSubstring("not implemented yet"));
    REQUIRE_THROWS_WITH(OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 0.25}),
                        Catch::Matchers::ContainsSubstring("not implemented yet"));
}
