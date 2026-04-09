#include "tucvrp/rooted_tree.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <sstream>

using tucvrp::Instance;
using tucvrp::RootedTreeBuilder;

TEST_CASE("rooted tree builder exposes parent children and terminal distance data") {
    std::istringstream input(R"(
5 0
0 1 1
1 2 2
1 3 3
0 4 4
3
2 0.4
3 0.5
4 0.2
)");

    const auto instance = Instance::parse(input);
    const auto rooted_tree = RootedTreeBuilder::build(instance);

    REQUIRE(rooted_tree.depot == 0);
    REQUIRE(rooted_tree.parent[0] == -1);
    REQUIRE(rooted_tree.parent[1] == 0);
    REQUIRE(rooted_tree.parent[2] == 1);
    REQUIRE(rooted_tree.parent[3] == 1);
    REQUIRE(rooted_tree.parent[4] == 0);
    REQUIRE(rooted_tree.distances_from_depot[2] == Catch::Approx(3.0));
    REQUIRE(rooted_tree.distances_from_depot[3] == Catch::Approx(4.0));
    REQUIRE(rooted_tree.distances_from_depot[4] == Catch::Approx(4.0));
    REQUIRE(rooted_tree.is_terminal(2));
    REQUIRE(rooted_tree.is_terminal(3));
    REQUIRE(rooted_tree.is_terminal(4));
    REQUIRE_FALSE(rooted_tree.is_terminal(1));
    REQUIRE(rooted_tree.children[0] == std::vector<int>{1, 4});
    REQUIRE(rooted_tree.children[1] == std::vector<int>{2, 3});
}
