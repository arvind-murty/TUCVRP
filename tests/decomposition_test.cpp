#include "tucvrp/decomposition.hpp"

#include <catch2/catch_test_macros.hpp>

#include <sstream>

using tucvrp::DecompositionBuilder;
using tucvrp::Instance;
using tucvrp::RootedTreeBuilder;

TEST_CASE("trivial decomposition contains one region at each level") {
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
    const auto rooted_tree = RootedTreeBuilder::build(instance);
    const auto decomposition = DecompositionBuilder::make_trivial(rooted_tree);

    REQUIRE(decomposition.depot == instance.depot());
    REQUIRE(decomposition.components.size() == 1);
    REQUIRE(decomposition.blocks.size() == 1);
    REQUIRE(decomposition.clusters.size() == 1);
    REQUIRE(decomposition.cells.size() == 1);

    REQUIRE(decomposition.components[0].root == instance.depot());
    REQUIRE(decomposition.components[0].block_ids == std::vector<int>{0});
    REQUIRE(decomposition.blocks[0].cluster_ids == std::vector<int>{0});
    REQUIRE(decomposition.clusters[0].cell_ids == std::vector<int>{0});
    REQUIRE(decomposition.cells[0].vertices == rooted_tree.vertices);
}
