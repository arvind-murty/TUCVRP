#include "internal.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace tucvrp {

namespace {

// In one postorder traversal, identify:
// 1. leaf-component roots,
// 2. the key vertices used by the internal-component decomposition.
void analyze_backbone(const RootedTreeData& rooted_tree,
                      double gamma,
                      std::vector<int>& leaf_roots,
                      std::vector<int>& key_vertices) {
    leaf_roots.clear();
    key_vertices.clear();
    std::vector<bool> is_leaf_root(rooted_tree.parent.size(), false);
    std::vector<bool> on_backbone(rooted_tree.parent.size(), false);

    // Use an explicit two-phase DFS stack to simulate recursive postorder:
    // first visit a node before its children, then revisit it after all children have been handled.
    std::vector<std::pair<int, bool>> stack{{rooted_tree.depot, false}};
    while (!stack.empty()) {
        const auto [u, expanded] = stack.back();
        stack.pop_back();

        if (!expanded) {
            // Schedule the postorder visit, then descend to the children.
            stack.emplace_back(u, true);
            for (auto it = rooted_tree.children[u].rbegin(); it != rooted_tree.children[u].rend(); ++it) {
                stack.emplace_back(*it, false);
            }
            continue;
        }

        // At this point every child summary is already known, so we can decide what role u plays.
        bool children_small = true;
        bool has_backbone_child = false;
        int backbone_child_count = 0;
        for (const int child : rooted_tree.children[u]) {
            // u is a leaf-component root only if its whole subtree is big but every child subtree is small.
            if (static_cast<double>(rooted_tree.subtree_terminal_counts[child]) >= gamma) {
                children_small = false;
            }
            // Any child already on the backbone means u is also on the backbone.
            if (on_backbone[child]) {
                has_backbone_child = true;
                ++backbone_child_count;
            }
        }

        if (!rooted_tree.children[u].empty() &&
            static_cast<double>(rooted_tree.subtree_terminal_counts[u]) >= gamma &&
            children_small) {
            is_leaf_root[u] = true;
            leaf_roots.push_back(u);
        }

        // The backbone is the union of depot-to-leaf-component-root paths.
        on_backbone[u] = is_leaf_root[u] || has_backbone_child;

        // Key vertices are the vertices where the internal-component decomposition can change direction:
        // the depot, every leaf-component root, and every backbone branching point.
        if (u == rooted_tree.depot || (on_backbone[u] && (is_leaf_root[u] || backbone_child_count >= 2))) {
            key_vertices.push_back(u);
        }
    }
}

// Walk upward until reaching the closest ancestor that is also a key vertex.
int lowest_key_ancestor(int vertex, const std::vector<bool>& is_key, const RootedTreeData& rooted_tree) {
    int current = rooted_tree.parent[vertex];
    while (current != -1 && !is_key[current]) {
        // Keep climbing until we hit the first ancestor that was marked as key.
        current = rooted_tree.parent[current];
    }
    if (current == -1) {
        throw std::logic_error("failed to find lowest key ancestor");
    }
    return current;
}

// Count terminals in the component rooted at `root` whose exit toward the depot is `exit`.
// For a leaf component, `exit == -1` and the whole subtree is included.
int component_terminal_count(const RootedTreeData& rooted_tree, int root, int exit) {
    // Internal components are represented as:
    //   subtree(root) minus subtree(exit), with exit kept as the boundary vertex.
    // Since subtree terminal counts are already available, the count is just a subtraction.
    if (exit == root || exit == -1) {
        return exit == -1 ? rooted_tree.subtree_terminal_counts[root] : 0;
    }
    return rooted_tree.subtree_terminal_counts[root] - rooted_tree.subtree_terminal_counts[exit];
}

// Collect the vertices in an internal component: the subtree of `root` with the exit subtree removed,
// plus the exit vertex itself so the attachment point is explicit in the representation.
std::vector<int> collect_internal_component_vertices(const RootedTreeData& rooted_tree, int root, int exit) {
    if (root == exit) {
        return {};
    }

    std::vector<int> vertices;
    std::vector<int> stack{root};
    while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        vertices.push_back(u);
        for (auto it = rooted_tree.children[u].rbegin(); it != rooted_tree.children[u].rend(); ++it) {
            // The exit subtree belongs to the next component closer to the depot, so omit it here.
            if (*it != exit) {
                stack.push_back(*it);
            }
        }
    }
    // Keep the exit vertex itself so later phases know where this component attaches upward.
    vertices.push_back(exit);
    return vertices;
}

}  // namespace

TreeDecomposition DecompositionBuilder::make_trivial(const RootedTreeData& rooted_tree) {
    TreeDecomposition decomposition;
    decomposition.depot = rooted_tree.depot;

    double total_demand = 0.0;
    for (const int terminal : rooted_tree.terminal_vertices) {
        total_demand += rooted_tree.demands[terminal];
    }

    decomposition.cells.push_back(Cell{
        .id = 0,
        .cluster_id = 0,
        .root = rooted_tree.depot,
        .exit = -1,
        .demand = total_demand,
        .vertices = rooted_tree.vertices,
    });

    decomposition.clusters.push_back(Cluster{
        .id = 0,
        .block_id = 0,
        .root = rooted_tree.depot,
        .exit = -1,
        .demand = total_demand,
        .vertices = rooted_tree.vertices,
        .cell_ids = {0},
    });

    decomposition.blocks.push_back(Block{
        .id = 0,
        .component_id = 0,
        .root = rooted_tree.depot,
        .exit = -1,
        .demand = total_demand,
        .vertices = rooted_tree.vertices,
        .cluster_ids = {0},
    });

    decomposition.components.push_back(Component{
        .id = 0,
        .root = rooted_tree.depot,
        .exit = -1,
        .terminal_count = static_cast<int>(rooted_tree.terminal_vertices.size()),
        .is_leaf = true,
        .is_big = true,
        .vertices = rooted_tree.vertices,
        .block_ids = {0},
    });

    return decomposition;
}

// Decompose a bounded-distance tree into components following Algorithm 5 / Lemma 9.
// In the current project we use the paper's threshold Gamma = 12 / epsilon with k = 1.
TreeDecomposition DecompositionBuilder::decompose_bounded_instance(const RootedTreeData& rooted_tree,
                                                                  double epsilon) {
    if (epsilon <= 0.0 || epsilon >= 1.0) {
        throw std::invalid_argument("decompose_bounded_instance requires epsilon in (0, 1)");
    }

    TreeDecomposition decomposition;
    decomposition.depot = rooted_tree.depot;

    const double gamma = 12.0 / epsilon;
    const std::vector<int> depth = decomposition_detail::compute_depths(rooted_tree);

    // Steps 1 and 2: a single postorder pass identifies leaf-component roots and key vertices.
    std::vector<int> leaf_roots;
    std::vector<int> key_vertices;
    analyze_backbone(rooted_tree, gamma, leaf_roots, key_vertices);
    for (const int v : leaf_roots) {
        // Each leaf-component root contributes one full big subtree component.
        decomposition_detail::append_component(
            decomposition,
            v,
            -1,
            rooted_tree.subtree_terminal_counts[v],
            true,
            true,
            decomposition_detail::collect_subtree_vertices(rooted_tree, v));
    }

    // If no leaf component exists, the whole bounded instance stays as one component.
    if (leaf_roots.empty()) {
        decomposition_detail::append_component(
            decomposition,
            rooted_tree.depot,
            -1,
            rooted_tree.subtree_terminal_counts[rooted_tree.depot],
            true,
            static_cast<double>(rooted_tree.subtree_terminal_counts[rooted_tree.depot]) >= gamma,
            rooted_tree.vertices);
    } else {
        std::vector<bool> is_key(rooted_tree.parent.size(), false);
        for (const int v : key_vertices) {
            is_key[v] = true;
        }

        // Process each non-root key vertex from top to bottom, decomposing the path segment between it
        // and its lowest key ancestor into maximal big internal components plus at most one final small one.
        std::vector<int> non_root_keys;
        for (const int v : key_vertices) {
            if (v != rooted_tree.depot) {
                non_root_keys.push_back(v);
            }
        }
        // Process key vertices top-down so each path segment [v1, v2] is handled before any deeper
        // segment that hangs below it.
        std::sort(non_root_keys.begin(), non_root_keys.end(), [&depth](int a, int b) { return depth[a] < depth[b]; });

        for (const int v2 : non_root_keys) {
            // v2 is the lower endpoint of a backbone segment, and v1 is the closest key ancestor above it.
            // The open path between them contains no other key vertices.
            const int v1 = lowest_key_ancestor(v2, is_key, rooted_tree);
            int x = v2;

            // While the segment between v1 and x still contains at least Gamma terminals, peel off one
            // maximal big internal component by choosing the lowest root whose component remains big.
            while (static_cast<double>(component_terminal_count(rooted_tree, v1, x)) >= gamma) {
                int chosen = v1;
                // Scan upward from x toward v1. The first candidate that still leaves a big component
                // is the lowest possible root, hence the maximal component we can peel off next.
                for (int candidate = x; candidate != v1; candidate = rooted_tree.parent[candidate]) {
                    if (static_cast<double>(component_terminal_count(rooted_tree, candidate, x)) >= gamma) {
                        chosen = candidate;
                        break;
                    }
                }

                const int terminal_count = component_terminal_count(rooted_tree, chosen, x);
                decomposition_detail::append_component(
                    decomposition,
                    chosen,
                    x,
                    terminal_count,
                    false,
                    static_cast<double>(terminal_count) >= gamma,
                    collect_internal_component_vertices(rooted_tree, chosen, x));
                // The remaining unexplained segment is now the path from v1 down to the new boundary `chosen`.
                x = chosen;
            }

            // Any leftover segment on the path contains fewer than Gamma terminals, so it becomes the
            // final small internal component between this pair of consecutive key vertices.
            if (v1 != x) {
                const int terminal_count = component_terminal_count(rooted_tree, v1, x);
                decomposition_detail::append_component(
                    decomposition,
                    v1,
                    x,
                    terminal_count,
                    false,
                    static_cast<double>(terminal_count) >= gamma,
                    collect_internal_component_vertices(rooted_tree, v1, x));
            }
        }
    }

    decompose_components_into_blocks(decomposition, rooted_tree, epsilon);
    decompose_blocks_into_clusters(decomposition, rooted_tree, epsilon);
    decompose_clusters_into_cells(decomposition, rooted_tree, epsilon);
    return decomposition;
}

}  // namespace tucvrp
