#include "tucvrp/preprocessing.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <stdexcept>

namespace tucvrp {

namespace {

constexpr double kEps = 1e-9;

// Return the next unused vertex id for introducing auxiliary nodes.
int next_available_vertex_id(const Instance& instance) {
    int next_id = 0;
    for (int v : instance.vertices()) {
        next_id = std::max(next_id, v + 1);
    }
    return next_id;
}

}

// Measure root-to-terminal distances and check whether they satisfy the bounded-distance regime.
BoundedDistanceStats Preprocessor::bounded_distance_stats(const Instance& instance, double epsilon) {
    if (epsilon <= 0.0 || epsilon >= 1.0) {
        throw std::invalid_argument("epsilon must be in (0, 1)");
    }
    instance.validate();

    const auto dist = instance.distance_from_depot();
    double min_distance = std::numeric_limits<double>::infinity();
    double max_distance = 0.0;

    for (const auto& terminal : instance.terminals()) {
        min_distance = std::min(min_distance, dist[terminal.vertex]);
        max_distance = std::max(max_distance, dist[terminal.vertex]);
    }

    if (!std::isfinite(min_distance)) {
        throw std::invalid_argument("instance must contain at least one terminal");
    }

    const double exponent = (1.0 / epsilon) - 1.0;
    const double ratio_limit = std::pow(1.0 / epsilon, exponent);
    return BoundedDistanceStats{
        .min_distance = min_distance,
        .max_distance = max_distance,
        .bounded = max_distance <= ratio_limit * min_distance + kEps,
    };
}

// Sum the mandatory load contribution of every edge based on subtree demand.
double Preprocessor::edge_load_lower_bound(const Instance& instance) {
    instance.validate();
    const auto parent = instance.parent_array();
    std::vector<int> order;
    order.reserve(instance.vertex_count());

    std::vector<std::pair<int, int>> stack;
    stack.emplace_back(instance.depot(), -1);
    while (!stack.empty()) {
        const auto [u, p] = stack.back();
        stack.pop_back();
        order.push_back(u);
        for (auto it = instance.neighbors(u).rbegin(); it != instance.neighbors(u).rend(); ++it) {
            const auto& edge = *it;
            if (edge.to != p) {
                stack.emplace_back(edge.to, u);
            }
        }
    }

    std::vector<double> subtree_demand(instance.vertex_count(), 0.0);
    for (const auto& terminal : instance.terminals()) {
        subtree_demand[terminal.vertex] += terminal.demand;
    }

    double lower_bound = 0.0;
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        int u = *it;
        if (u == instance.depot()) {
            continue;
        }
        int p = parent[u];
        double weight = -1.0;
        for (const auto& edge : instance.neighbors(u)) {
            if (edge.to == p) {
                weight = edge.weight;
                break;
            }
        }
        if (weight < 0.0) {
            throw std::logic_error("missing parent edge");
        }
        lower_bound += 2.0 * weight * std::ceil(subtree_demand[u] - kEps);
        subtree_demand[p] += subtree_demand[u];
    }
    return lower_bound;
}

// Normalize an instance so every terminal is a leaf and every vertex has at most two children.
Instance Preprocessor::make_binary_leaf_tree(const Instance& instance) {
    instance.validate();

    Instance out(instance.depot());
    int next_id = next_available_vertex_id(instance);
    std::vector<std::pair<int, int>> stack;
    stack.emplace_back(instance.depot(), -1);

    while (!stack.empty()) {
        const auto [u, parent] = stack.back();
        stack.pop_back();

        std::vector<Edge> children;
        for (const auto& edge : instance.neighbors(u)) {
            if (edge.to != parent) {
                children.push_back(edge);
            }
        }

        if (instance.is_terminal(u) && children.empty()) {
            out.add_terminal(u, instance.demand_of(u));
        }

        std::vector<Edge> outgoing = children;
        if (instance.is_terminal(u) && !children.empty()) {
            const int promoted_leaf = next_id++;
            outgoing.push_back(Edge{promoted_leaf, 0.0});
            out.add_terminal(promoted_leaf, instance.demand_of(u));
        }

        if (outgoing.size() <= 2) {
            for (const auto& edge : outgoing) {
                out.add_edge(u, edge.to, edge.weight);
            }
        } else {
            int chain_node = u;
            for (std::size_t i = 0; i + 2 < outgoing.size(); ++i) {
                out.add_edge(chain_node, outgoing[i].to, outgoing[i].weight);
                const int next_chain = next_id++;
                out.add_edge(chain_node, next_chain, 0.0);
                chain_node = next_chain;
            }
            out.add_edge(chain_node, outgoing[outgoing.size() - 2].to, outgoing[outgoing.size() - 2].weight);
            out.add_edge(chain_node, outgoing.back().to, outgoing.back().weight);
        }

        for (auto it = children.rbegin(); it != children.rend(); ++it) {
            stack.emplace_back(it->to, u);
        }
    }

    out.validate();
    return out;
}

}  // namespace tucvrp
