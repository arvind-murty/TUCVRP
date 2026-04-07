#include "tucvrp/instance.hpp"
#include "tucvrp/preprocessing.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <sstream>

using tucvrp::Instance;
using tucvrp::Preprocessor;

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

TEST_CASE("make_tree_binary preserves demand and does not increase branching above two") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(0, 2, 1.0);
    instance.add_edge(0, 3, 1.0);
    instance.add_terminal(1, 0.2);
    instance.add_terminal(2, 0.3);
    instance.add_terminal(3, 0.4);
    instance.validate();

    auto binarized = Preprocessor::make_tree_binary(instance);

    REQUIRE(binarized.total_demand() == Catch::Approx(instance.total_demand()));
    for (int v : binarized.vertices()) {
        REQUIRE(binarized.neighbors(v).size() <= 3);
    }
    REQUIRE(binarized.depot() == 0);
}
