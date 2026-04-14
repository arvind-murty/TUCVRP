#include "internal.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace tucvrp {

namespace {

std::vector<bool> block_membership(const RootedTreeData& rooted_tree, const Block& block) {
    std::vector<bool> in_block(rooted_tree.parent.size(), false);
    for (const int v : block.vertices) {
        // Convert the explicit block vertex list into a constant-time membership map
        // for all subsequent local traversals inside this block.
        in_block[v] = true;
    }
    return in_block;
}

std::vector<int> block_children(const RootedTreeData& rooted_tree,
                                const std::vector<bool>& in_block,
                                int vertex) {
    std::vector<int> children;
    for (const int child : rooted_tree.children[vertex]) {
        // Only keep rooted children that still lie inside the current block.
        if (in_block[child]) {
            children.push_back(child);
        }
    }
    return children;
}

std::vector<int> block_postorder(const RootedTreeData& rooted_tree,
                                 const std::vector<bool>& in_block,
                                 int root) {
    std::vector<int> order;
    // Two-phase DFS to simulate recursive postorder without recursion.
    std::vector<std::pair<int, bool>> stack{{root, false}};
    while (!stack.empty()) {
        const auto [u, expanded] = stack.back();
        stack.pop_back();
        if (!expanded) {
            // Revisit u after every child in the block has been processed.
            stack.emplace_back(u, true);
            for (auto it = rooted_tree.children[u].rbegin(); it != rooted_tree.children[u].rend(); ++it) {
                if (in_block[*it]) {
                    stack.emplace_back(*it, false);
                }
            }
            continue;
        }
        order.push_back(u);
    }
    return order;
}

std::vector<double> block_subtree_demands(const RootedTreeData& rooted_tree,
                                          const Block& block,
                                          const std::vector<bool>& in_block) {
    std::vector<double> demand(rooted_tree.parent.size(), 0.0);
    for (const int v : block.vertices) {
        // Cluster demand only counts terminals strictly inside the block, so terminals at the
        // block root and block exit are excluded from every local subtree total.
        if (v != block.root && v != block.exit && rooted_tree.is_terminal(v)) {
            demand[v] = rooted_tree.demands[v];
        }
    }

    const std::vector<int> order = block_postorder(rooted_tree, in_block, block.root);
    for (const int u : order) {
        if (u == block.root) {
            continue;
        }
        const int p = rooted_tree.parent[u];
        if (p != -1 && in_block[p]) {
            // Push each child's already-computed local demand into its parent.
            demand[p] += demand[u];
        }
    }

    return demand;
}

std::vector<int> collect_block_subtree_vertices(const RootedTreeData& rooted_tree,
                                                const std::vector<bool>& in_block,
                                                int root) {
    std::vector<int> vertices;
    std::vector<int> stack{root};
    while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        vertices.push_back(u);
        for (auto it = rooted_tree.children[u].rbegin(); it != rooted_tree.children[u].rend(); ++it) {
            // This is an ordinary subtree walk, except that we ignore children leaving the block.
            if (in_block[*it]) {
                stack.push_back(*it);
            }
        }
    }
    return vertices;
}

double cluster_subgraph_demand(const std::vector<double>& subtree_demand, int root, int exit) {
    // Internal clusters are represented as subtree(root) minus subtree(exit),
    // with exit retained as the boundary vertex. The demand is therefore a subtraction.
    if (exit == -1) {
        return subtree_demand[root];
    }
    if (exit == root) {
        return 0.0;
    }
    return subtree_demand[root] - subtree_demand[exit];
}

std::vector<int> collect_internal_cluster_vertices(const RootedTreeData& rooted_tree,
                                                   const std::vector<bool>& in_block,
                                                   int root,
                                                   int exit) {
    if (root == exit) {
        // When the root and exit coincide, the cluster is just that boundary vertex.
        return {exit};
    }

    std::vector<int> vertices;
    std::vector<int> stack{root};
    while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        vertices.push_back(u);
        for (auto it = rooted_tree.children[u].rbegin(); it != rooted_tree.children[u].rend(); ++it) {
            // The exit subtree belongs to the next cluster toward the block root, so do not cross it.
            if (in_block[*it] && *it != exit) {
                stack.push_back(*it);
            }
        }
    }
    // Keep the exit vertex itself so later phases know where this cluster attaches upward.
    vertices.push_back(exit);
    return vertices;
}

int child_on_path_to_descendant(const RootedTreeData& rooted_tree,
                                const std::vector<bool>& in_block,
                                int ancestor,
                                int descendant) {
    int current = descendant;
    int child = -1;
    while (current != ancestor) {
        // Remember the last vertex before reaching the ancestor: that is the child of `ancestor`
        // lying on the ancestor-to-descendant path.
        child = current;
        current = rooted_tree.parent[current];
        if (current == -1 || !in_block[current]) {
            throw std::logic_error("invalid ancestor-descendant relation inside block");
        }
    }
    return child;
}

int deepest_big_cluster_root_on_path(const RootedTreeData& rooted_tree,
                                     const std::vector<double>& subtree_demand,
                                     int path_root,
                                     int exit,
                                     double gamma0) {
    for (int candidate = exit; candidate != path_root; candidate = rooted_tree.parent[candidate]) {
        // Scan upward from the lower boundary. The first candidate whose local demand is still at
        // least Gamma_0 is the deepest valid root for the next big internal cluster.
        if (cluster_subgraph_demand(subtree_demand, candidate, exit) >= gamma0) {
            return candidate;
        }
    }
    return path_root;
}

void analyze_block_backbone(const RootedTreeData& rooted_tree,
                            const Block& block,
                            const std::vector<bool>& in_block,
                            const std::vector<double>& subtree_demand,
                            double gamma0,
                            std::vector<int>& leaf_cluster_roots,
                            std::vector<int>& key_vertices) {
    leaf_cluster_roots.clear();
    key_vertices.clear();
    std::vector<bool> is_leaf_cluster_root(rooted_tree.parent.size(), false);
    std::vector<bool> on_backbone(rooted_tree.parent.size(), false);

    // Postorder traversal over the current block so every local child summary is available
    // before deciding whether the current vertex becomes a leaf-cluster root or a key vertex.
    std::vector<std::pair<int, bool>> stack{{block.root, false}};
    while (!stack.empty()) {
        const auto [u, expanded] = stack.back();
        stack.pop_back();
        if (!expanded) {
            stack.emplace_back(u, true);
            for (auto it = rooted_tree.children[u].rbegin(); it != rooted_tree.children[u].rend(); ++it) {
                if (in_block[*it]) {
                    stack.emplace_back(*it, false);
                }
            }
            continue;
        }

        const std::vector<int> children = block_children(rooted_tree, in_block, u);
        bool children_small = true;
        bool has_backbone_child = false;
        int backbone_child_count = 0;
        for (const int child : children) {
            // In Algorithm 3, a leaf-cluster root is a highest local subtree whose demand is at
            // least Gamma_0 while every child subtree has demand below Gamma_0.
            if (subtree_demand[child] >= gamma0) {
                children_small = false;
            }
            if (on_backbone[child]) {
                has_backbone_child = true;
                ++backbone_child_count;
            }
        }

        // Algorithm 3 leaf-cluster roots: a local subtree with demand at least Gamma_0 whose
        // child subtrees all have demand below Gamma_0.
        if (!children.empty() && subtree_demand[u] >= gamma0 && children_small) {
            is_leaf_cluster_root[u] = true;
            leaf_cluster_roots.push_back(u);
        }

        // The block backbone is the union of block-root-to-leaf-cluster-root paths.
        on_backbone[u] = is_leaf_cluster_root[u] || has_backbone_child;

        // Algorithm 4 key vertices for the cluster decomposition: the block root, every leaf-cluster
        // root, and every branching vertex of the block backbone.
        if (u == block.root || (on_backbone[u] && (is_leaf_cluster_root[u] || backbone_child_count >= 2))) {
            key_vertices.push_back(u);
        }
    }
}

int lowest_block_key_ancestor(int vertex,
                              const std::vector<bool>& is_key,
                              const std::vector<bool>& in_block,
                              const RootedTreeData& rooted_tree) {
    int current = rooted_tree.parent[vertex];
    while (current != -1 && in_block[current] && !is_key[current]) {
        current = rooted_tree.parent[current];
    }
    if (current == -1 || !in_block[current]) {
        throw std::logic_error("failed to find lowest block-key ancestor");
    }
    return current;
}

void create_leaf_clusters(TreeDecomposition& decomposition,
                          const RootedTreeData& rooted_tree,
                          const Block& block,
                          const std::vector<bool>& in_block,
                          const std::vector<double>& subtree_demand,
                          const std::vector<int>& leaf_cluster_roots) {
    for (const int root : leaf_cluster_roots) {
        int exit = -1;
        // If the block exit lies inside this leaf cluster, record it as the cluster exit.
        if (block.exit != -1) {
            for (int v = block.exit; v != -1 && in_block[v]; v = rooted_tree.parent[v]) {
                if (v == root) {
                    exit = block.exit;
                    break;
                }
                if (v == block.root) {
                    break;
                }
            }
        }

        decomposition_detail::append_cluster(
            decomposition,
            block.id,
            root,
            exit,
            cluster_subgraph_demand(subtree_demand, root, exit),
            // A leaf cluster is exactly the subtree rooted at `root` within the block.
            collect_block_subtree_vertices(rooted_tree, in_block, root));
    }
}

bool block_exit_is_covered(const TreeDecomposition& decomposition, const Block& block) {
    if (block.exit == -1) {
        return true;
    }
    for (const int cluster_id : block.cluster_ids) {
        if (decomposition.clusters[cluster_id].exit == block.exit) {
            return true;
        }
    }
    return false;
}

void ensure_block_exit_cluster(TreeDecomposition& decomposition, const Block& block) {
    if (block.exit == -1 || block_exit_is_covered(decomposition, block)) {
        return;
    }

    // Algorithm 3 explicitly adds a trivial leaf cluster for the block exit when no other
    // leaf cluster already contains it.
    decomposition_detail::append_cluster(
        decomposition, block.id, block.exit, block.exit, 0.0, {block.exit});
}

void create_internal_clusters(TreeDecomposition& decomposition,
                              const RootedTreeData& rooted_tree,
                              const Block& block,
                              const std::vector<bool>& in_block,
                              const std::vector<double>& subtree_demand,
                              const std::vector<int>& key_vertices,
                              double gamma0) {
    std::vector<bool> is_key(rooted_tree.parent.size(), false);
    for (const int v : key_vertices) {
        is_key[v] = true;
    }

    std::vector<int> non_root_keys;
    for (const int v : key_vertices) {
        if (v != block.root) {
            non_root_keys.push_back(v);
        }
    }

    const std::vector<int> depth = decomposition_detail::compute_depths(rooted_tree);
    // Process key vertices top-down so each path segment is split before any deeper segment nested below it.
    std::sort(non_root_keys.begin(), non_root_keys.end(), [&depth](int a, int b) { return depth[a] < depth[b]; });

    for (const int v2 : non_root_keys) {
        // v2 is the lower key vertex on the current backbone segment, and v1 is the closest
        // key ancestor above it inside the current block.
        const int v1 = lowest_block_key_ancestor(v2, is_key, in_block, rooted_tree);
        const int v1_child = child_on_path_to_descendant(rooted_tree, in_block, v1, v2);
        int x = v2;

        // Algorithm 4 repeatedly peels off the deepest internal cluster H(v) whose demand
        // is still at least Gamma_0.
        while (cluster_subgraph_demand(subtree_demand, v1_child, x) >= gamma0) {
            const int chosen = deepest_big_cluster_root_on_path(rooted_tree, subtree_demand, v1_child, x, gamma0);
            decomposition_detail::append_cluster(
                decomposition,
                block.id,
                chosen,
                x,
                cluster_subgraph_demand(subtree_demand, chosen, x),
                collect_internal_cluster_vertices(rooted_tree, in_block, chosen, x));
            // The remaining unexplained segment now stops at the newly created cluster root.
            x = chosen;
        }

        // The remaining top segment is the final small internal cluster with block-key root v1.
        std::vector<int> vertices{v1};
        const std::vector<int> remainder = collect_internal_cluster_vertices(rooted_tree, in_block, v1_child, x);
        vertices.insert(vertices.end(), remainder.begin(), remainder.end());
        decomposition_detail::append_cluster(
            decomposition,
            block.id,
            v1,
            x,
            cluster_subgraph_demand(subtree_demand, v1_child, x),
            std::move(vertices));
    }
}

void decompose_block_into_clusters(TreeDecomposition& decomposition,
                                   const RootedTreeData& rooted_tree,
                                   int block_id,
                                   double epsilon) {
    const Block& block = decomposition.blocks[block_id];
    const double gamma = 12.0 / epsilon;
    const double gamma0 = decomposition_detail::big_terminal_threshold(epsilon, gamma);
    const std::vector<bool> in_block = block_membership(rooted_tree, block);
    const std::vector<double> subtree_demand = block_subtree_demands(rooted_tree, block, in_block);

    // Algorithms 3 and 4 are driven by leaf-cluster roots and block key vertices.
    std::vector<int> leaf_cluster_roots;
    std::vector<int> key_vertices;
    analyze_block_backbone(rooted_tree, block, in_block, subtree_demand, gamma0, leaf_cluster_roots, key_vertices);

    // If no finer split is needed, the whole block itself is a single cluster.
    if (leaf_cluster_roots.empty() && key_vertices.size() == 1 && key_vertices[0] == block.root) {
        decomposition_detail::append_cluster(
            decomposition, block.id, block.root, block.exit, block.demand, block.vertices);
        return;
    }

    create_leaf_clusters(decomposition, rooted_tree, block, in_block, subtree_demand, leaf_cluster_roots);
    const bool needs_exit_singleton = (block.exit != -1 && !block_exit_is_covered(decomposition, block));
    ensure_block_exit_cluster(decomposition, block);
    if (needs_exit_singleton &&
        std::find(key_vertices.begin(), key_vertices.end(), block.exit) == key_vertices.end()) {
        // The singleton exit cluster behaves like an extra leaf cluster endpoint in Algorithm 4,
        // so promote the block exit to a key vertex before decomposing the remaining internal part.
        key_vertices.push_back(block.exit);
    }
    create_internal_clusters(decomposition, rooted_tree, block, in_block, subtree_demand, key_vertices, gamma0);

    // If Algorithms 3 and 4 do not create any finer split, the whole block itself is one cluster.
    if (decomposition.blocks[block_id].cluster_ids.empty()) {
        decomposition_detail::append_cluster(
            decomposition, block_id, block.root, block.exit, block.demand, block.vertices);
    }
}

}  // namespace

void DecompositionBuilder::decompose_blocks_into_clusters(TreeDecomposition& decomposition,
                                                          const RootedTreeData& rooted_tree,
                                                          double epsilon) {
    if (epsilon <= 0.0 || epsilon >= 1.0) {
        throw std::invalid_argument("decompose_blocks_into_clusters requires epsilon in (0, 1)");
    }

    decomposition.clusters.clear();
    for (auto& block : decomposition.blocks) {
        block.cluster_ids.clear();
    }

    for (int block_id = 0; block_id < static_cast<int>(decomposition.blocks.size()); ++block_id) {
        decompose_block_into_clusters(decomposition, rooted_tree, block_id, epsilon);
    }
}

}  // namespace tucvrp
