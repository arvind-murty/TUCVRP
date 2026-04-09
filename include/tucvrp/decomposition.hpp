#pragma once

#include "tucvrp/rooted_tree.hpp"

#include <vector>

namespace tucvrp {

struct Cell {
    int id = -1;
    int root = -1;
    std::vector<int> vertices;
};

struct Cluster {
    int id = -1;
    int root = -1;
    std::vector<int> vertices;
    std::vector<int> cell_ids;
};

struct Block {
    int id = -1;
    int root = -1;
    std::vector<int> vertices;
    std::vector<int> cluster_ids;
};

struct Component {
    int id = -1;
    int root = -1;
    std::vector<int> vertices;
    std::vector<int> block_ids;
};

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
