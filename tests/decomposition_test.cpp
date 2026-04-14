#include "tucvrp/decomposition.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <stdexcept>
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
    REQUIRE(decomposition.blocks[0].component_id == 0);
    REQUIRE(decomposition.blocks[0].exit == -1);
    REQUIRE(decomposition.blocks[0].cluster_ids == std::vector<int>{0});
    REQUIRE(decomposition.clusters[0].cell_ids == std::vector<int>{0});
    REQUIRE(decomposition.cells[0].cluster_id == 0);
    REQUIRE(decomposition.cells[0].exit == -1);
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

    REQUIRE(decomposition.blocks.size() == 3);

    const auto& leaf_component = decomposition.components[0];
    REQUIRE(leaf_component.block_ids.size() == 2);
    REQUIRE(decomposition.blocks[leaf_component.block_ids[0]].root == 1);
    REQUIRE(decomposition.blocks[leaf_component.block_ids[0]].exit == -1);
    REQUIRE(decomposition.blocks[leaf_component.block_ids[1]].root == 1);
    REQUIRE(decomposition.blocks[leaf_component.block_ids[1]].exit == -1);

    const auto& internal_component = decomposition.components[1];
    REQUIRE(internal_component.block_ids.size() == 1);
    REQUIRE(decomposition.blocks[internal_component.block_ids[0]].root == 0);
    REQUIRE(decomposition.blocks[internal_component.block_ids[0]].exit == 1);
    REQUIRE(decomposition.blocks[internal_component.block_ids[0]].demand == 0.0);
}

TEST_CASE("block decomposition splits a component at big terminals") {
    std::istringstream input(R"(
3 0
0 1 1
1 2 1
1
2 0.2
)");

    const auto instance = Instance::parse(input);
    const auto rooted_tree = RootedTreeBuilder::build(instance);
    const auto decomposition = DecompositionBuilder::decompose_bounded_instance(rooted_tree, 0.9);

    REQUIRE(decomposition.components.size() == 1);
    REQUIRE(decomposition.blocks.size() == 1);

    const auto& component = decomposition.components[0];
    REQUIRE(component.block_ids.size() == 1);

    const auto& block = decomposition.blocks[component.block_ids[0]];
    REQUIRE(block.root == 0);
    REQUIRE(block.exit == 2);
    REQUIRE(block.vertices == std::vector<int>{0, 1, 2});
    REQUIRE(block.demand == 0.0);
}

TEST_CASE("height reduction groups same-class components under one critical vertex") {
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
    const auto reduced =
        DecompositionBuilder::height_reduce_bounded_components(decomposition, rooted_tree, 0.9);

    REQUIRE(reduced.groups.size() == 1);
    REQUIRE(reduced.groups[0].component_ids == std::vector<int>{0, 1});
    REQUIRE(reduced.groups[0].critical_vertex == 0);

    REQUIRE(reduced.original_parent_component[0] == 1);
    REQUIRE(reduced.original_parent_component[1] == -1);
    REQUIRE(reduced.class_index_by_component[0] == 1);
    REQUIRE(reduced.class_index_by_component[1] == 1);
    REQUIRE(reduced.group_id_by_component[0] == 0);
    REQUIRE(reduced.group_id_by_component[1] == 0);
    REQUIRE(reduced.critical_vertex_by_component[0] == 0);
    REQUIRE(reduced.critical_vertex_by_component[1] == 0);
    REQUIRE(reduced.attachment_length_by_component[0] == Catch::Approx(1.0));
    REQUIRE(reduced.attachment_length_by_component[1] == Catch::Approx(0.0));
}

TEST_CASE("height reduction separates components that lie in different distance classes") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(1, 2, 1.0);
    instance.add_edge(2, 3, 1.0);
    instance.add_edge(3, 4, 1.0);
    instance.add_edge(4, 5, 1.0);
    instance.add_terminal(5, 0.2);
    instance.validate();

    const auto rooted_tree = RootedTreeBuilder::build(instance);

    tucvrp::TreeDecomposition decomposition;
    decomposition.depot = 0;
    decomposition.components.push_back(tucvrp::Component{
        .id = 0,
        .root = 0,
        .exit = 4,
        .terminal_count = 0,
        .is_leaf = false,
        .is_big = false,
        .vertices = {0, 1, 2, 3, 4},
        .block_ids = {},
    });
    decomposition.components.push_back(tucvrp::Component{
        .id = 1,
        .root = 4,
        .exit = -1,
        .terminal_count = 1,
        .is_leaf = true,
        .is_big = true,
        .vertices = {4, 5},
        .block_ids = {},
    });

    const auto reduced =
        DecompositionBuilder::height_reduce_bounded_components(decomposition, rooted_tree, 0.5);

    REQUIRE(reduced.groups.size() == 2);
    REQUIRE(reduced.original_parent_component == std::vector<int>{-1, 0});
    REQUIRE(reduced.class_index_by_component[0] != reduced.class_index_by_component[1]);
    REQUIRE(reduced.group_id_by_component[0] != reduced.group_id_by_component[1]);
    REQUIRE(reduced.critical_vertex_by_component[0] == 0);
    REQUIRE(reduced.critical_vertex_by_component[1] == 4);
    REQUIRE(reduced.attachment_length_by_component[0] == Catch::Approx(0.0));
    REQUIRE(reduced.attachment_length_by_component[1] == Catch::Approx(0.0));
}

TEST_CASE("lifting a height-reduced solution preserves terminal partition and tour costs") {
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
    const tucvrp::SolveResult reduced_solution{
        .cost = 14.0,
        .tours = {
            tucvrp::Tour{.terminals = {2, 3}, .demand = 0.0, .cost = 14.0},
        },
    };

    const auto lifted =
        DecompositionBuilder::lift_solution_from_height_reduced_tree(reduced_solution, instance);

    REQUIRE(lifted.cost == Catch::Approx(14.0));
    REQUIRE(lifted.tours.size() == 1);
    REQUIRE(lifted.tours[0].terminals == std::vector<int>{2, 3});
    REQUIRE(lifted.tours[0].demand == Catch::Approx(0.9));
    REQUIRE(lifted.tours[0].cost == Catch::Approx(14.0));
}

TEST_CASE("lifting a height-reduced solution rejects invalid tours") {
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

    SECTION("non-terminal vertex") {
        const tucvrp::SolveResult reduced_solution{
            .cost = 2.0,
            .tours = {
                tucvrp::Tour{.terminals = {1}, .demand = 0.0, .cost = 2.0},
            },
        };

        REQUIRE_THROWS_AS(
            DecompositionBuilder::lift_solution_from_height_reduced_tree(reduced_solution, instance),
            std::invalid_argument);
    }

    SECTION("over-capacity terminal set") {
        const tucvrp::SolveResult reduced_solution{
            .cost = 12.0,
            .tours = {
                tucvrp::Tour{.terminals = {2, 2, 3}, .demand = 0.0, .cost = 12.0},
            },
        };

        REQUIRE_THROWS_AS(
            DecompositionBuilder::lift_solution_from_height_reduced_tree(reduced_solution, instance),
            std::invalid_argument);
    }
}

TEST_CASE("bounded decomposition satisfies basic component and block invariants") {
    Instance instance(0);
    instance.add_edge(0, 1, 1.0);
    instance.add_edge(1, 2, 1.0);
    instance.add_edge(1, 3, 1.0);
    instance.add_edge(2, 4, 0.0);
    instance.add_edge(2, 5, 0.0);
    instance.add_edge(3, 6, 0.0);
    instance.add_edge(3, 7, 0.0);
    instance.add_terminal(4, 0.11);
    instance.add_terminal(5, 0.12);
    instance.add_terminal(6, 0.13);
    instance.add_terminal(7, 0.14);
    instance.validate();

    const auto rooted_tree = RootedTreeBuilder::build(instance);
    const auto decomposition = DecompositionBuilder::decompose_bounded_instance(rooted_tree, 0.9);

    for (const auto& component : decomposition.components) {
        REQUIRE_FALSE(component.vertices.empty());
        REQUIRE(component.block_ids.size() >= 1);
        REQUIRE(component.vertices.front() == component.root);
        if (component.exit != -1) {
            REQUIRE(std::find(component.vertices.begin(), component.vertices.end(), component.exit) !=
                    component.vertices.end());
        }

        for (const int block_id : component.block_ids) {
            const auto& block = decomposition.blocks[block_id];
            REQUIRE(block.component_id == component.id);
            REQUIRE_FALSE(block.vertices.empty());
            REQUIRE(block.vertices.front() == block.root);
            REQUIRE(std::find(component.vertices.begin(), component.vertices.end(), block.root) !=
                    component.vertices.end());
            if (block.exit != -1) {
                REQUIRE(std::find(block.vertices.begin(), block.vertices.end(), block.exit) != block.vertices.end());
            }

            double expected_demand = 0.0;
            for (const int v : block.vertices) {
                if (v != block.root && v != block.exit && rooted_tree.is_terminal(v)) {
                    expected_demand += rooted_tree.demands[v];
                }
            }
            REQUIRE(block.demand == Catch::Approx(expected_demand));
        }
    }
}

TEST_CASE("cluster decomposition creates clusters consistent with their blocks") {
    std::istringstream input(R"(
3 0
0 1 1
1 2 1
1
2 0.2
)");

    const auto instance = Instance::parse(input);
    const auto rooted_tree = RootedTreeBuilder::build(instance);
    const auto decomposition = DecompositionBuilder::decompose_bounded_instance(rooted_tree, 0.9);

    REQUIRE(decomposition.blocks.size() == 1);
    REQUIRE(decomposition.clusters.size() >= 1);

    const auto& block = decomposition.blocks[0];
    REQUIRE_FALSE(block.cluster_ids.empty());

    for (const int cluster_id : block.cluster_ids) {
        const auto& cluster = decomposition.clusters[cluster_id];
        REQUIRE(cluster.block_id == block.id);
        REQUIRE_FALSE(cluster.vertices.empty());
        REQUIRE(cluster.vertices.front() == cluster.root);
        REQUIRE(std::find(block.vertices.begin(), block.vertices.end(), cluster.root) != block.vertices.end());
        if (cluster.exit != -1) {
            REQUIRE(std::find(cluster.vertices.begin(), cluster.vertices.end(), cluster.exit) != cluster.vertices.end());
        }

        double expected_demand = 0.0;
        for (const int v : cluster.vertices) {
            if (v != cluster.root && v != cluster.exit && rooted_tree.is_terminal(v)) {
                expected_demand += rooted_tree.demands[v];
            }
        }
        REQUIRE(cluster.demand == Catch::Approx(expected_demand));
    }
}

TEST_CASE("cluster decomposition can fall back to a whole-block cluster") {
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
    auto decomposition = DecompositionBuilder::make_trivial(rooted_tree);

    const auto original_cells = decomposition.cells.size();

    DecompositionBuilder::decompose_blocks_into_clusters(decomposition, rooted_tree, 0.25);

    REQUIRE(decomposition.blocks[0].cluster_ids.size() == 1);
    REQUIRE(decomposition.clusters.size() == 1);
    REQUIRE(decomposition.clusters[0].root == decomposition.blocks[0].root);
    REQUIRE(decomposition.clusters[0].exit == decomposition.blocks[0].exit);
    REQUIRE(decomposition.clusters[0].vertices == decomposition.blocks[0].vertices);
    REQUIRE(decomposition.clusters[0].demand == Catch::Approx(decomposition.blocks[0].demand));

    DecompositionBuilder::decompose_clusters_into_cells(decomposition, rooted_tree, 0.25);

    REQUIRE(decomposition.cells.size() == original_cells);
    REQUIRE(decomposition.clusters[0].cell_ids == std::vector<int>{0});
}

TEST_CASE("ending clusters decompose into a single cell") {
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
    auto decomposition = DecompositionBuilder::make_trivial(rooted_tree);

    DecompositionBuilder::decompose_clusters_into_cells(decomposition, rooted_tree, 0.25);

    REQUIRE(decomposition.clusters.size() == 1);
    REQUIRE(decomposition.cells.size() == 1);
    REQUIRE(decomposition.clusters[0].cell_ids == std::vector<int>{0});
    REQUIRE(decomposition.cells[0].cluster_id == 0);
    REQUIRE(decomposition.cells[0].root == decomposition.clusters[0].root);
    REQUIRE(decomposition.cells[0].exit == -1);
    REQUIRE(decomposition.cells[0].vertices == decomposition.clusters[0].vertices);
    REQUIRE(decomposition.cells[0].demand == Catch::Approx(decomposition.clusters[0].demand));
}

TEST_CASE("passing clusters decompose into cells by cutting the spine") {
    std::istringstream input(R"(
6 0
0 1 1
1 2 1
2 3 1
3 4 1
4 5 1
1
5 0.2
)");

    const auto instance = Instance::parse(input);
    const auto rooted_tree = RootedTreeBuilder::build(instance);
    const auto decomposition = DecompositionBuilder::decompose_bounded_instance(rooted_tree, 0.25);

    const auto* passing_cluster = static_cast<const tucvrp::Cluster*>(nullptr);
    for (const auto& cluster : decomposition.clusters) {
        if (cluster.exit == 5 && cluster.vertices.size() > 1) {
            passing_cluster = &cluster;
            break;
        }
    }
    REQUIRE(passing_cluster != nullptr);
    const auto& cluster = *passing_cluster;
    REQUIRE(cluster.exit == 5);
    REQUIRE(decomposition.cells.size() == 4);
    REQUIRE(cluster.cell_ids.size() == 4);

    for (const int cell_id : cluster.cell_ids) {
        const auto& cell = decomposition.cells[cell_id];
        REQUIRE(cell.cluster_id == cluster.id);
        REQUIRE_FALSE(cell.vertices.empty());
        REQUIRE(cell.exit != -1);
        REQUIRE(std::find(cell.vertices.begin(), cell.vertices.end(), cell.root) != cell.vertices.end());
        REQUIRE(std::find(cell.vertices.begin(), cell.vertices.end(), cell.exit) != cell.vertices.end());
        REQUIRE(cell.demand >= 0.0);
    }
}
