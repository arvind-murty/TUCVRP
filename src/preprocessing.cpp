#include "tucvrp/preprocessing.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace tucvrp {

namespace {

constexpr double kEps = 1e-9;

int next_available_vertex_id(const Instance& instance) {
    int next_id = 0;
    for (int v : instance.vertices()) {
        next_id = std::max(next_id, v + 1);
    }
    return next_id;
}

}

BoundedDistanceStats Preprocessor::bounded_distance_stats(const Instance& instance, double epsilon) {
    if (epsilon <= 0.0 || epsilon >= 1.0) {
        throw std::invalid_argument("epsilon must be in (0, 1)");
    }
    instance.validate();

    const auto dist = instance.distance_from_root();
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

double Preprocessor::edge_load_lower_bound(const Instance& instance) {
    instance.validate();
    const auto parent = instance.parent_array();
    std::vector<int> order;
    order.reserve(instance.vertex_count());

    std::function<void(int, int)> dfs = [&](int u, int p) {
        order.push_back(u);
        for (const auto& edge : instance.neighbors(u)) {
            if (edge.to != p) {
                dfs(edge.to, u);
            }
        }
    };
    dfs(instance.root(), -1);

    std::vector<double> subtree_demand(instance.vertex_count(), 0.0);
    for (const auto& terminal : instance.terminals()) {
        subtree_demand[terminal.vertex] += terminal.demand;
    }

    double lower_bound = 0.0;
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        int u = *it;
        if (u == instance.root()) {
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

Instance Preprocessor::binarize_tree(const Instance& instance) {
    instance.validate();

    Instance out(instance.root());
    int next_id = next_available_vertex_id(instance);
    const int copy_bound = std::max(next_id, instance.vertex_count());
    std::vector<bool> copied_terminal(copy_bound, false);

    std::function<void(int, int)> dfs = [&](int u, int parent) {
        if (instance.is_terminal(u) && !copied_terminal[u]) {
            out.add_terminal(u, instance.demand_of(u));
            copied_terminal[u] = true;
        }

        std::vector<Edge> children;
        for (const auto& edge : instance.neighbors(u)) {
            if (edge.to != parent) {
                children.push_back(edge);
            }
        }

        if (children.size() <= 2) {
            for (const auto& child : children) {
                out.add_edge(u, child.to, child.weight);
            }
        } else {
            int chain_node = u;
            for (std::size_t i = 0; i < children.size() - 2; ++i) {
                if (i == 0) {
                    out.add_edge(chain_node, children[i].to, children[i].weight);
                }
                int next_chain = next_id++;
                out.add_edge(chain_node, next_chain, 0.0);
                chain_node = next_chain;
            }
            out.add_edge(chain_node, children[children.size() - 2].to, children[children.size() - 2].weight);
            out.add_edge(chain_node, children.back().to, children.back().weight);
        }

        for (const auto& child : children) {
            dfs(child.to, u);
        }
    };

    dfs(instance.root(), -1);
    out.validate();
    return out;
}

}  // namespace tucvrp
