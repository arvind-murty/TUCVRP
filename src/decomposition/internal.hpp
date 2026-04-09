#pragma once

#include "tucvrp/decomposition.hpp"

#include <vector>

namespace tucvrp::decomposition_detail {

double big_terminal_threshold(double epsilon, double gamma);

std::vector<int> collect_subtree_vertices(const RootedTreeData& rooted_tree, int root);
std::vector<int> compute_depths(const RootedTreeData& rooted_tree);

void analyze_backbone(const RootedTreeData& rooted_tree,
                      double gamma,
                      std::vector<int>& leaf_roots,
                      std::vector<int>& key_vertices);

int lowest_key_ancestor(int vertex, const std::vector<bool>& is_key, const RootedTreeData& rooted_tree);
int component_terminal_count(const RootedTreeData& rooted_tree, int root, int exit);
std::vector<int> collect_internal_component_vertices(const RootedTreeData& rooted_tree, int root, int exit);

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

std::vector<bool> component_membership(const RootedTreeData& rooted_tree, const Component& component);
std::vector<bool> spanning_subtree_for_block_decomposition(const RootedTreeData& rooted_tree,
                                                           const Component& component,
                                                           const std::vector<bool>& in_component,
                                                           double gamma0);
std::vector<bool> block_key_vertices(const RootedTreeData& rooted_tree,
                                     const Component& component,
                                     const std::vector<bool>& in_component,
                                     const std::vector<bool>& in_tu,
                                     double gamma0);
double block_demand(const RootedTreeData& rooted_tree,
                    const std::vector<int>& vertices,
                    int root,
                    int exit);

}  // namespace tucvrp::decomposition_detail
