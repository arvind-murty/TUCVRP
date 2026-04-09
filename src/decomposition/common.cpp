#include "internal.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tucvrp::decomposition_detail {

double big_terminal_threshold(double epsilon, double gamma) {
    // Definition 8: alpha = epsilon^(1/epsilon + 1) and Gamma_0 = epsilon * alpha / Gamma.
    // Any terminal with demand greater than Gamma_0 is considered big inside its component.
    const double alpha = std::pow(epsilon, 1.0 / epsilon + 1.0);
    return epsilon * alpha / gamma;
}

// Collect every vertex in the rooted subtree of `root`.
std::vector<int> collect_subtree_vertices(const RootedTreeData& rooted_tree, int root) {
    std::vector<int> vertices;
    std::vector<int> stack{root};
    while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        vertices.push_back(u);
        // Push children in reverse so the final DFS order is stable with respect to the
        // child order already stored in the rooted tree.
        for (auto it = rooted_tree.children[u].rbegin(); it != rooted_tree.children[u].rend(); ++it) {
            stack.push_back(*it);
        }
    }
    return vertices;
}

// Compute depth from the depot in the rooted tree so key vertices can be processed top-down.
std::vector<int> compute_depths(const RootedTreeData& rooted_tree) {
    std::vector<int> depth(rooted_tree.parent.size(), -1);
    std::vector<int> stack{rooted_tree.depot};
    depth[rooted_tree.depot] = 0;
    while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        for (const int child : rooted_tree.children[u]) {
            // Because the tree is already rooted, each child has a unique parent and a
            // uniquely determined depth.
            depth[child] = depth[u] + 1;
            stack.push_back(child);
        }
    }
    return depth;
}

// Append one component record with the next available id.
void append_component(TreeDecomposition& decomposition,
                      int root,
                      int exit,
                      int terminal_count,
                      bool is_leaf,
                      bool is_big,
                      std::vector<int> vertices) {
    decomposition.components.push_back(Component{
        .id = static_cast<int>(decomposition.components.size()),
        .root = root,
        .exit = exit,
        .terminal_count = terminal_count,
        .is_leaf = is_leaf,
        .is_big = is_big,
        .vertices = std::move(vertices),
        .block_ids = {},
    });
}

void append_block(TreeDecomposition& decomposition,
                  int component_id,
                  int root,
                  int exit,
                  double demand,
                  std::vector<int> vertices) {
    const int block_id = static_cast<int>(decomposition.blocks.size());
    decomposition.blocks.push_back(Block{
        .id = block_id,
        .component_id = component_id,
        .root = root,
        .exit = exit,
        .demand = demand,
        .vertices = std::move(vertices),
        .cluster_ids = {},
    });
    decomposition.components[component_id].block_ids.push_back(block_id);
}

void append_cluster(TreeDecomposition& decomposition,
                    int block_id,
                    int root,
                    int exit,
                    double demand,
                    std::vector<int> vertices) {
    const int cluster_id = static_cast<int>(decomposition.clusters.size());
    decomposition.clusters.push_back(Cluster{
        .id = cluster_id,
        .block_id = block_id,
        .root = root,
        .exit = exit,
        .demand = demand,
        .vertices = std::move(vertices),
        .cell_ids = {},
    });
    decomposition.blocks[block_id].cluster_ids.push_back(cluster_id);
}

}  // namespace tucvrp::decomposition_detail
