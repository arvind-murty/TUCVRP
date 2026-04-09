#pragma once

#include "tucvrp/decomposition.hpp"

#include <vector>

namespace tucvrp::decomposition_detail {

double big_terminal_threshold(double epsilon, double gamma);

std::vector<int> collect_subtree_vertices(const RootedTreeData& rooted_tree, int root);
std::vector<int> compute_depths(const RootedTreeData& rooted_tree);

void append_component(TreeDecomposition& decomposition,
                      int root,
                      int exit,
                      int terminal_count,
                      bool is_leaf,
                      bool is_big,
                      std::vector<int> vertices);

void append_block(TreeDecomposition& decomposition,
                  int component_id,
                  int root,
                  int exit,
                  double demand,
                  std::vector<int> vertices);
void append_cluster(TreeDecomposition& decomposition,
                    int block_id,
                    int root,
                    int exit,
                    double demand,
                    std::vector<int> vertices);
void append_cell(TreeDecomposition& decomposition,
                 int cluster_id,
                 int root,
                 int exit,
                 double demand,
                 std::vector<int> vertices);

}  // namespace tucvrp::decomposition_detail
