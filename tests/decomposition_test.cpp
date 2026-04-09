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

TEST_CASE("bounded instance decomposition builds leaf and internal components") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(1, 2, 1.0);
    instance.add_edge(1, 3, 1.0);

    for (int i = 0; i < 8; ++i) {
        const int left_internal = 4 + 2 * i;
        const int left_leaf = left_internal + 1;
        instance.add_edge(2, left_internal, 0.0);
        instance.add_edge(left_internal, left_leaf, 0.0);
        instance.add_terminal(left_leaf, 0.01);
    }

    for (int i = 0; i < 9; ++i) {
        const int right_internal = 20 + 2 * i;
        const int right_leaf = right_internal + 1;
        instance.add_edge(3, right_internal, 0.0);
        instance.add_edge(right_internal, right_leaf, 0.0);
        instance.add_terminal(right_leaf, 0.01);
    }

    instance.validate();

    const auto rooted_tree = RootedTreeBuilder::build(instance);
    const auto decomposition = DecompositionBuilder::decompose_bounded_instance(rooted_tree, 0.9);

    REQUIRE(decomposition.components.size() == 2);
    REQUIRE(decomposition.components[0].root == 1);
    REQUIRE(decomposition.components[0].exit == -1);
    REQUIRE(decomposition.components[0].terminal_count == 17);
    REQUIRE(decomposition.components[0].is_leaf);
    REQUIRE(decomposition.components[0].is_big);

    REQUIRE(decomposition.components[1].root == 0);
    REQUIRE(decomposition.components[1].exit == 1);
    REQUIRE(decomposition.components[1].terminal_count == 0);
    REQUIRE_FALSE(decomposition.components[1].is_leaf);
    REQUIRE_FALSE(decomposition.components[1].is_big);
}
