#include "tucvrp/algorithms/labbe_approx.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace tucvrp {
namespace {

constexpr double kTolerance = 1e-9;

struct Aggregate {
    double demand = 0.0;
    std::vector<int> terminals;
};

Tour make_tour(const Instance& instance, Aggregate aggregate) {
    std::sort(aggregate.terminals.begin(), aggregate.terminals.end());
    aggregate.terminals.erase(std::unique(aggregate.terminals.begin(), aggregate.terminals.end()), aggregate.terminals.end());

    Tour tour;
    tour.terminals = std::move(aggregate.terminals);
    tour.demand = aggregate.demand;
    tour.cost = instance.tour_cost_for_terminals(tour.terminals);
    tour.walk = instance.tour_walk_for_terminals(tour.terminals);
    return tour;
}

void append_terminals(std::vector<int>& dst, const std::vector<int>& src) {
    dst.insert(dst.end(), src.begin(), src.end());
}

// Collect the vertices in one depot-child subtree. With no fixed route cost in the current
// objective, these subtrees can be processed independently and their tours concatenated.
std::vector<int> subtree_vertices_from_child(const std::vector<std::vector<int>>& children, int child_root) {
    std::vector<int> subtree;
    std::vector<int> stack{child_root};
    while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        subtree.push_back(u);
        for (const int v : children[u]) {
            stack.push_back(v);
        }
    }
    return subtree;
}

// Pick the deepest active leaf in the current subtree, breaking ties by vertex id. This replaces
// the paper's unspecified "arbitrary leaf" choice with a deterministic rule for reproducible tests.
int choose_leaf(const std::vector<int>& subtree_vertices,
                const std::vector<bool>& active,
                const std::vector<int>& active_child_count,
                const std::vector<int>& depth,
                int subtree_root) {
    int best = -1;
    for (const int v : subtree_vertices) {
        if (!active[v] || active_child_count[v] != 0 || v == subtree_root) {
            continue;
        }
        if (best == -1 || depth[v] > depth[best] || (depth[v] == depth[best] && v < best)) {
            best = v;
        }
    }
    return best;
}

// Execute the adapted Labbé leaf-elimination heuristic on one depot-child subtree.
std::vector<Tour> solve_subtree(const Instance& instance,
                                const std::vector<int>& subtree_vertices,
                                const std::vector<int>& parent,
                                const std::vector<int>& depth,
                                int subtree_root) {
    std::vector<Tour> tours;

    const int vertex_count = instance.vertex_count();
    std::vector<bool> in_subtree(vertex_count, false);
    std::vector<bool> active(vertex_count, false);
    std::vector<int> active_child_count(vertex_count, 0);
    std::vector<Aggregate> aggregates(vertex_count);

    for (const int v : subtree_vertices) {
        in_subtree[v] = true;
        active[v] = true;
        aggregates[v].demand = instance.demand_of(v);
        if (instance.is_terminal(v)) {
            aggregates[v].terminals.push_back(v);
        }
    }

    for (const int v : subtree_vertices) {
        for (const auto& edge : instance.neighbors(v)) {
            if (edge.to < 0 || edge.to >= vertex_count || !in_subtree[edge.to]) {
                continue;
            }
            if (parent[edge.to] == v) {
                ++active_child_count[v];
            }
        }
    }

    while (true) {
        const int leaf = choose_leaf(subtree_vertices, active, active_child_count, depth, subtree_root);
        if (leaf == -1) {
            break;
        }

        const int parent_vertex = parent[leaf];
        if (parent_vertex < 0 || !active[parent_vertex]) {
            throw std::logic_error("labbe heuristic found a leaf whose parent is no longer active");
        }

        const double leaf_demand = aggregates[leaf].demand;
        const double parent_demand = aggregates[parent_vertex].demand;

        if (leaf_demand + parent_demand <= 1.0 + kTolerance) {
            // If the leaf aggregate and parent aggregate fit together, bubble the leaf aggregate up
            // and delete the leaf from the active tree.
            aggregates[parent_vertex].demand += leaf_demand;
            append_terminals(aggregates[parent_vertex].terminals, aggregates[leaf].terminals);
        } else if (leaf_demand > parent_demand + kTolerance) {
            // If they do not fit and the leaf aggregate is heavier, finalize the leaf as one route.
            if (leaf_demand > kTolerance) {
                tours.push_back(make_tour(instance, aggregates[leaf]));
            }
        } else {
            // Otherwise finalize the current parent aggregate, then replace the parent's residual
            // aggregate by the leaf aggregate and continue pushing that residual upward later.
            if (parent_demand > kTolerance) {
                tours.push_back(make_tour(instance, aggregates[parent_vertex]));
            }
            aggregates[parent_vertex] = aggregates[leaf];
        }

        active[leaf] = false;
        aggregates[leaf] = Aggregate{};
        --active_child_count[parent_vertex];
    }

    // After all non-root leaves are eliminated, the subtree root carries the final residual
    // aggregate for this depot-child subtree.
    if (active[subtree_root] && aggregates[subtree_root].demand > kTolerance) {
        tours.push_back(make_tour(instance, aggregates[subtree_root]));
    }

    return tours;
}

}  // namespace

SolveResult LabbeApproxSolver::solve(const Instance& instance) {
    instance.validate();

    SolveResult result;
    if (instance.terminal_count() == 0) {
        return result;
    }

    const int depot = instance.depot();
    const auto parent = instance.parent_array();

    std::vector<int> depth(instance.vertex_count(), 0);
    for (int v = 0; v < instance.vertex_count(); ++v) {
        if (v == depot || parent[v] < 0) {
            continue;
        }
        depth[v] = depth[parent[v]] + 1;
    }

    std::vector<std::vector<int>> children(instance.vertex_count());
    for (int v = 0; v < instance.vertex_count(); ++v) {
        if (v != depot && parent[v] >= 0) {
            children[parent[v]].push_back(v);
        }
    }

    // Solve each depot-child subtree independently. With no fixed cost per route, no route needs
    // to pass through the depot as an internal node.
    for (const int child : children[depot]) {
        const std::vector<int> subtree_vertices = subtree_vertices_from_child(children, child);
        std::vector<Tour> subtree_tours = solve_subtree(instance, subtree_vertices, parent, depth, child);
        for (Tour& tour : subtree_tours) {
            result.cost += tour.cost;
            result.tours.push_back(std::move(tour));
        }
    }

    return result;
}

}  // namespace tucvrp
