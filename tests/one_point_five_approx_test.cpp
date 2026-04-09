#include "tucvrp/algorithms/one_point_five_approx.hpp"
#include "tucvrp/rng.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include <stdexcept>

using tucvrp::Instance;
using tucvrp::OnePointFiveApproxParams;
using tucvrp::OnePointFiveApproxSolver;
using tucvrp::Rng;

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
    REQUIRE_THROWS_WITH(OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 0.25}),
                        Catch::Matchers::ContainsSubstring("not implemented yet"));
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
    REQUIRE_THROWS_WITH(OnePointFiveApproxSolver::solve(instance, OnePointFiveApproxParams{.epsilon = 0.25}),
                        Catch::Matchers::ContainsSubstring("not implemented yet"));
}
