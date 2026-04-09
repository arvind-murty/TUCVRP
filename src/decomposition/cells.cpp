#include "internal.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tucvrp {

namespace {

using namespace decomposition_detail;

std::vector<bool> cluster_membership(const RootedTreeData& rooted_tree, const Cluster& cluster) {
    std::vector<bool> in_cluster(rooted_tree.parent.size(), false);
    for (const int v : cluster.vertices) {
        // Turn the explicit cluster vertex list into a constant-time membership test
        // for all local traversals inside this cluster.
        in_cluster[v] = true;
    }
    return in_cluster;
}

std::vector<int> spine_vertices(const RootedTreeData& rooted_tree,
                                const std::vector<bool>& in_cluster,
                                const Cluster& cluster) {
    std::vector<int> spine;
    if (cluster.exit == -1) {
        return spine;
    }

    // The spine is the unique path from the cluster root to the cluster exit.
    for (int v = cluster.exit; v != cluster.root; v = rooted_tree.parent[v]) {
        if (v == -1 || !in_cluster[v]) {
            throw std::logic_error("cluster exit is not below the cluster root");
        }
        spine.push_back(v);
    }
    spine.push_back(cluster.root);
    std::reverse(spine.begin(), spine.end());
    return spine;
}

double edge_weight_on_tree(const RootedTreeData& rooted_tree, int parent, int child) {
    // Rooted-tree distances are additive, so the local edge weight is the child distance
    // minus the parent distance.
    const double weight = rooted_tree.distances_from_depot[child] - rooted_tree.distances_from_depot[parent];
    if (weight < -1e-9) {
        throw std::logic_error("negative edge weight along rooted tree path");
    }
    return std::max(0.0, weight);
}

double cluster_spine_cost(const RootedTreeData& rooted_tree, const std::vector<int>& spine) {
    double cost = 0.0;
    for (std::size_t i = 1; i < spine.size(); ++i) {
        // Sum the edge costs along the root-to-exit spine.
        cost += edge_weight_on_tree(rooted_tree, spine[i - 1], spine[i]);
    }
    return cost;
}

std::set<std::pair<int, int>> spine_cut_edges(const RootedTreeData& rooted_tree,
                                              const std::vector<int>& spine,
                                              double epsilon) {
    std::set<std::pair<int, int>> cuts;
    if (spine.size() <= 1) {
        return cuts;
    }

    const int one_over_epsilon = static_cast<int>(std::round(1.0 / epsilon));
    const double spine_cost = cluster_spine_cost(rooted_tree, spine);
    if (one_over_epsilon <= 1 || spine_cost <= 1e-9) {
        return cuts;
    }

    std::vector<double> prefix(spine.size(), 0.0);
    for (std::size_t i = 1; i < spine.size(); ++i) {
        // prefix[i] is the cost of the spine prefix ending at spine[i].
        prefix[i] = prefix[i - 1] + edge_weight_on_tree(rooted_tree, spine[i - 1], spine[i]);
    }

    // Section 4.3: for each i in {1, ..., 1/epsilon - 1}, remove the unique spine edge
    // whose distance interval contains i * epsilon * l_x.
    for (int i = 1; i < one_over_epsilon; ++i) {
        const double threshold = i * epsilon * spine_cost;
        for (std::size_t j = 1; j < spine.size(); ++j) {
            // The threshold lies in exactly one edge interval of the prefix sums; removing that
            // edge makes the resulting cell spines have cost at most an epsilon fraction of l_x.
            if (prefix[j - 1] <= threshold && threshold < prefix[j]) {
                cuts.emplace(spine[j - 1], spine[j]);
                break;
            }
        }
    }

    return cuts;
}

std::vector<std::vector<int>> cell_components(const RootedTreeData& rooted_tree,
                                              const std::vector<bool>& in_cluster,
                                              const std::set<std::pair<int, int>>& cuts,
                                              const Cluster& cluster) {
    std::vector<std::vector<int>> components;
    std::vector<bool> visited(rooted_tree.parent.size(), false);
    // After removing cut edges from the spine, each connected component contains exactly one
    // highest spine vertex: either the cluster root, or the lower endpoint of a removed spine edge.
    std::vector<int> seeds{cluster.root};
    for (const auto& [_, lower] : cuts) {
        seeds.push_back(lower);
    }

    for (const int start : seeds) {
        if (visited[start]) {
            continue;
        }
        std::vector<int> component;
        std::vector<int> stack{start};
        visited[start] = true;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            component.push_back(u);

            if (rooted_tree.parent[u] != -1 && in_cluster[rooted_tree.parent[u]]) {
                const auto edge = std::make_pair(rooted_tree.parent[u], u);
                // The only forbidden upward move is across a removed spine edge.
                if (!cuts.contains(edge) && !visited[rooted_tree.parent[u]]) {
                    visited[rooted_tree.parent[u]] = true;
                    stack.push_back(rooted_tree.parent[u]);
                }
            }

            for (const int child : rooted_tree.children[u]) {
                if (!in_cluster[child]) {
                    continue;
                }
                const auto edge = std::make_pair(u, child);
                // Likewise, descend through every cluster edge except the removed spine edges.
                if (!cuts.contains(edge) && !visited[child]) {
                    visited[child] = true;
                    stack.push_back(child);
                }
            }
        }
        components.push_back(std::move(component));
    }

    return components;
}

std::pair<int, int> cell_root_and_exit(const std::vector<bool>& in_component,
                                       const std::vector<int>& spine) {
    int root = -1;
    int exit = -1;
    for (const int v : spine) {
        if (!in_component[v]) {
            continue;
        }
        // The first spine vertex inside this connected component is the cell root.
        if (root == -1) {
            root = v;
        }
        // The last spine vertex inside this connected component is the cell exit.
        exit = v;
    }

    if (root == -1) {
        throw std::logic_error("passing cell contains no spine vertex");
    }
    return {root, exit};
}

double cell_demand(const RootedTreeData& rooted_tree,
                   const std::vector<int>& vertices,
                   int root,
                   int exit) {
    double demand = 0.0;
    for (const int v : vertices) {
        // As for blocks and clusters, only terminals strictly inside the cell contribute to its demand.
        if (v != root && v != exit && rooted_tree.is_terminal(v)) {
            demand += rooted_tree.demands[v];
        }
    }
    return demand;
}

void decompose_ending_cluster_into_cell(TreeDecomposition& decomposition,
                                        const Cluster& cluster) {
    // Section 4.3: an ending cluster is already a single cell.
    append_cell(decomposition, cluster.id, cluster.root, -1, cluster.demand, cluster.vertices);
}

}  // namespace

void DecompositionBuilder::decompose_clusters_into_cells(TreeDecomposition& decomposition,
                                                         const RootedTreeData& rooted_tree,
                                                         double epsilon) {
    if (epsilon <= 0.0 || epsilon >= 1.0) {
        throw std::invalid_argument("decompose_clusters_into_cells requires epsilon in (0, 1)");
    }

    decomposition.cells.clear();
    for (auto& cluster : decomposition.clusters) {
        cluster.cell_ids.clear();
    }

    for (const Cluster& cluster : decomposition.clusters) {
        if (cluster.exit == -1) {
            decompose_ending_cluster_into_cell(decomposition, cluster);
            continue;
        }

        const std::vector<bool> in_cluster = cluster_membership(rooted_tree, cluster);
        const std::vector<int> spine = spine_vertices(rooted_tree, in_cluster, cluster);
        const std::set<std::pair<int, int>> cuts = spine_cut_edges(rooted_tree, spine, epsilon);
        const std::vector<std::vector<int>> components = cell_components(rooted_tree, in_cluster, cuts, cluster);

        // Removing the selected spine edges partitions the cluster into connected subgraphs,
        // and Section 4.3 defines each such subgraph to be a cell.
        for (const auto& component_vertices : components) {
            std::vector<bool> in_component(rooted_tree.parent.size(), false);
            for (const int v : component_vertices) {
                // Build membership for this one connected component so we can recover the
                // highest and lowest spine vertices lying inside it.
                in_component[v] = true;
            }

            const auto [root, exit] = cell_root_and_exit(in_component, spine);
            append_cell(
                decomposition,
                cluster.id,
                root,
                exit,
                cell_demand(rooted_tree, component_vertices, root, exit),
                component_vertices);
        }
    }
}

}  // namespace tucvrp
