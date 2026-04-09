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

std::vector<bool> component_membership(const RootedTreeData& rooted_tree, const Component& component) {
    std::vector<bool> in_component(rooted_tree.parent.size(), false);
    for (const int v : component.vertices) {
        // Components already store their explicit vertex sets, so convert that list into a fast
        // membership bitmap for the block-decomposition helpers below.
        in_component[v] = true;
    }
    return in_component;
}

// Build the subtree T_U spanning the vertices in U from Section 4.1:
// the big terminals in the component, the component root, and possibly the component exit.
std::vector<bool> spanning_subtree_for_block_decomposition(const RootedTreeData& rooted_tree,
                                                           const Component& component,
                                                           const std::vector<bool>& in_component,
                                                           double gamma0) {
    std::vector<bool> in_u(rooted_tree.parent.size(), false);
    // By definition, the component root always belongs to U.
    in_u[component.root] = true;
    if (component.exit != -1) {
        // Internal components also include their exit vertex in U.
        in_u[component.exit] = true;
    }

    for (const int v : component.vertices) {
        // Big terminals are exactly the terminals of this component whose demand exceeds Gamma_0.
        if (v != component.root && v != component.exit && rooted_tree.is_terminal(v) &&
            rooted_tree.demands[v] > gamma0) {
            in_u[v] = true;
        }
    }

    std::vector<bool> in_tu(rooted_tree.parent.size(), false);
    for (const int start : component.vertices) {
        if (!in_u[start]) {
            continue;
        }
        // The subtree spanning U is the union of all root-to-u paths for u in U.
        // Since the input is already rooted, each such path is recovered by repeatedly following parent pointers.
        for (int v = start; v != component.root; v = rooted_tree.parent[v]) {
            if (v == -1 || !in_component[v]) {
                throw std::logic_error("component root is not an ancestor of a block key vertex");
            }
            in_tu[v] = true;
        }
        // Ensure the component root itself is marked once at the end of every path.
        in_tu[component.root] = true;
    }

    return in_tu;
}

std::vector<bool> block_key_vertices(const RootedTreeData& rooted_tree,
                                     const Component& component,
                                     const std::vector<bool>& in_component,
                                     const std::vector<bool>& in_tu,
                                     double gamma0) {
    std::vector<bool> is_block_key(rooted_tree.parent.size(), false);

    for (const int v : component.vertices) {
        if (!in_tu[v]) {
            // Vertices outside T_U do not participate in the block split.
            continue;
        }

        // A vertex belongs to U if it is the component root, the component exit,
        // or a big terminal in this component.
        bool in_u = (v == component.root || v == component.exit);
        if (rooted_tree.is_terminal(v) && v != component.root && v != component.exit &&
            rooted_tree.demands[v] > gamma0) {
            in_u = true;
        }

        int tu_child_count = 0;
        for (const int child : rooted_tree.children[v]) {
            if (in_component[child] && in_tu[child]) {
                // Only children that stay inside both the component and T_U count toward
                // the "two children in T_U" condition from Section 4.1.
                ++tu_child_count;
            }
        }

        // Section 4.1: a key vertex is either in U or has two children in T_U.
        if (in_u || tu_child_count >= 2) {
            is_block_key[v] = true;
        }
    }

    return is_block_key;
}

double block_demand(const RootedTreeData& rooted_tree,
                    const std::vector<int>& vertices,
                    int root,
                    int exit) {
    double demand = 0.0;
    for (const int v : vertices) {
        // By definition, only terminals strictly inside the block contribute to its demand.
        // So terminals at the root or the exit are excluded.
        if (v != root && v != exit && rooted_tree.is_terminal(v)) {
            demand += rooted_tree.demands[v];
        }
    }
    return demand;
}

}  // namespace tucvrp::decomposition_detail
