#include "tucvrp/decomposition.hpp"

namespace tucvrp {

TreeDecomposition DecompositionBuilder::make_trivial(const RootedTreeData& rooted_tree) {
    TreeDecomposition decomposition;
    decomposition.depot = rooted_tree.depot;

    decomposition.cells.push_back(Cell{
        .id = 0,
        .root = rooted_tree.depot,
        .vertices = rooted_tree.vertices,
    });

    decomposition.clusters.push_back(Cluster{
        .id = 0,
        .root = rooted_tree.depot,
        .vertices = rooted_tree.vertices,
        .cell_ids = {0},
    });

    decomposition.blocks.push_back(Block{
        .id = 0,
        .root = rooted_tree.depot,
        .vertices = rooted_tree.vertices,
        .cluster_ids = {0},
    });

    decomposition.components.push_back(Component{
        .id = 0,
        .root = rooted_tree.depot,
        .vertices = rooted_tree.vertices,
        .block_ids = {0},
    });

    return decomposition;
}

}  // namespace tucvrp
