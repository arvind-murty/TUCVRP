#pragma once

#include "tucvrp/solution.hpp"
#include "tucvrp/rooted_tree.hpp"

#include <vector>

namespace tucvrp {

// Lowest-level region used by the paper's decomposition hierarchy.
struct Cell {
    int id = -1;
    int cluster_id = -1;
    int root = -1;
    int exit = -1;
    double demand = 0.0;
    std::vector<int> vertices;
};

// A cluster groups cells inside one block.
struct Cluster {
    int id = -1;
    int block_id = -1;
    int root = -1;
    int exit = -1;
    double demand = 0.0;
    std::vector<int> vertices;
    std::vector<int> cell_ids;
};

// A block groups clusters inside one component.
struct Block {
    int id = -1;
    int component_id = -1;
    int root = -1;
    int exit = -1;
    double demand = 0.0;
    std::vector<int> vertices;
    std::vector<int> cluster_ids;
};

// A top-level component in the paper's recursive decomposition.
struct Component {
    int id = -1;
    int root = -1;
    int exit = -1;
    int terminal_count = 0;
    bool is_leaf = false;
    bool is_big = false;
    std::vector<int> vertices;
    std::vector<int> block_ids;
};

// One maximally connected set of components inside a single distance class C_i from Section 4.
struct HeightReducedComponentGroup {
    int id = -1;
    int class_index = -1;
    int critical_vertex = -1;
    std::vector<int> component_ids;
};

// Metadata for the Section 4 height reduction on the tree of components.
struct HeightReducedComponentTree {
    double d_min = 0.0;
    double d_tilde = 0.0;
    int h_epsilon = 0;
    std::vector<int> original_parent_component;
    std::vector<int> class_index_by_component;
    std::vector<int> group_id_by_component;
    std::vector<int> critical_vertex_by_component;
    std::vector<double> attachment_length_by_component;
    std::vector<HeightReducedComponentGroup> groups;
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
    // Build the component decomposition from Algorithm 5 using Gamma = 12 / epsilon and k = 1.
    static TreeDecomposition decompose_bounded_instance(const RootedTreeData& rooted_tree, double epsilon);
    // Apply the Section 4 height reduction to the current tree of components.
    static HeightReducedComponentTree height_reduce_bounded_components(const TreeDecomposition& decomposition,
                                                                      const RootedTreeData& rooted_tree,
                                                                      double epsilon);
    // Convert a solution on the height-reduced tree back to a feasible solution on the original bounded tree.
    static SolveResult lift_solution_from_height_reduced_tree(const SolveResult& reduced_solution,
                                                              const Instance& bounded_instance);
    // Split the current components into blocks using Section 4.1.
    static void decompose_components_into_blocks(TreeDecomposition& decomposition,
                                                 const RootedTreeData& rooted_tree,
                                                 double epsilon);
    // Placeholder for the Section 4.2 block-to-cluster decomposition.
    static void decompose_blocks_into_clusters(TreeDecomposition& decomposition,
                                               const RootedTreeData& rooted_tree,
                                               double epsilon);
    // Placeholder for the Section 4.3 cluster-to-cell decomposition.
    static void decompose_clusters_into_cells(TreeDecomposition& decomposition,
                                              const RootedTreeData& rooted_tree,
                                              double epsilon);
};

}  // namespace tucvrp
