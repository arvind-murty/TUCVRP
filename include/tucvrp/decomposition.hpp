#pragma once

#include "tucvrp/rooted_tree.hpp"

#include <vector>

namespace tucvrp {

// Lowest-level region used by the paper's decomposition hierarchy.
struct Cell {
    int id = -1;
    int root = -1;
    std::vector<int> vertices;
};

// A cluster groups cells inside one block.
struct Cluster {
    int id = -1;
    int root = -1;
    std::vector<int> vertices;
    std::vector<int> cell_ids;
};

// A block groups clusters inside one component.
struct Block {
    int id = -1;
    int root = -1;
    std::vector<int> vertices;
    std::vector<int> cluster_ids;
};

// A top-level component in the paper's recursive decomposition.
struct Component {
    int id = -1;
    int root = -1;
    std::vector<int> vertices;
    std::vector<int> block_ids;
};

// Full hierarchy of derived regions built over one rooted tree.
struct TreeDecomposition {
    int depot = -1;
    std::vector<Component> components;
    std::vector<Block> blocks;
    std::vector<Cluster> clusters;
    std::vector<Cell> cells;
};

class DecompositionBuilder {
  public:
    // Build the trivial decomposition where the whole tree is one component/block/cluster/cell.
    static TreeDecomposition make_trivial(const RootedTreeData& rooted_tree);
};

}  // namespace tucvrp
