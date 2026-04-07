#include "tucvrp/instance.hpp"

#include <cmath>
#include <fstream>
#include <limits>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>

namespace tucvrp {

namespace {

constexpr double kEps = 1e-9;

}

Instance::Instance(int root) : root_(root), max_vertex_id_(root) {
    ensure_vertex(root);
}

void Instance::ensure_vertex(int vertex) {
    if (vertex < 0) {
        throw std::invalid_argument("vertex ids must be non-negative");
    }
    if (vertex >= static_cast<int>(adjacency_.size())) {
        adjacency_.resize(vertex + 1);
    }
    if (vertex > max_vertex_id_) {
        max_vertex_id_ = vertex;
    }
}

void Instance::set_root(int root) {
    ensure_vertex(root);
    root_ = root;
}

void Instance::add_edge(int u, int v, double weight) {
    if (weight < 0.0) {
        throw std::invalid_argument("edge weight must be non-negative");
    }
    ensure_vertex(u);
    ensure_vertex(v);
    adjacency_[u].push_back(Edge{v, weight});
    adjacency_[v].push_back(Edge{u, weight});
}

void Instance::add_terminal(int vertex, double demand) {
    if (demand <= 0.0 || demand > 1.0 + kEps) {
        throw std::invalid_argument("terminal demand must be in (0, 1]");
    }
    ensure_vertex(vertex);
    if (demand_by_vertex_.contains(vertex)) {
        throw std::invalid_argument("duplicate terminal");
    }
    terminals_.push_back(Terminal{vertex, demand});
    demand_by_vertex_[vertex] = demand;
}

int Instance::root() const { return root_; }

int Instance::vertex_count() const { return static_cast<int>(adjacency_.size()); }

int Instance::edge_count() const {
    int total = 0;
    for (const auto& edges : adjacency_) {
        total += static_cast<int>(edges.size());
    }
    return total / 2;
}

int Instance::terminal_count() const { return static_cast<int>(terminals_.size()); }

double Instance::total_demand() const {
    double total = 0.0;
    for (const auto& terminal : terminals_) {
        total += terminal.demand;
    }
    return total;
}

bool Instance::is_terminal(int vertex) const { return demand_by_vertex_.contains(vertex); }

double Instance::demand_of(int vertex) const {
    auto it = demand_by_vertex_.find(vertex);
    return it == demand_by_vertex_.end() ? 0.0 : it->second;
}

const std::vector<Edge>& Instance::neighbors(int vertex) const {
    if (vertex < 0 || vertex >= static_cast<int>(adjacency_.size())) {
        throw std::out_of_range("vertex out of range");
    }
    return adjacency_[vertex];
}

const std::vector<Terminal>& Instance::terminals() const { return terminals_; }

std::vector<int> Instance::vertices() const {
    std::vector<int> out;
    out.reserve(adjacency_.size());
    for (int v = 0; v < static_cast<int>(adjacency_.size()); ++v) {
        if (!adjacency_[v].empty() || v == root_ || demand_by_vertex_.contains(v)) {
            out.push_back(v);
        }
    }
    return out;
}

void Instance::validate() const {
    if (root_ < 0) {
        throw std::invalid_argument("root must be set");
    }
    const auto verts = vertices();
    if (verts.empty()) {
        throw std::invalid_argument("instance must contain at least one vertex");
    }
    if (edge_count() != static_cast<int>(verts.size()) - 1) {
        throw std::invalid_argument("graph must be a tree");
    }

    std::vector<bool> seen(adjacency_.size(), false);
    std::queue<int> queue;
    queue.push(root_);
    seen[root_] = true;

    while (!queue.empty()) {
        int u = queue.front();
        queue.pop();
        for (const auto& edge : adjacency_[u]) {
            if (!seen[edge.to]) {
                seen[edge.to] = true;
                queue.push(edge.to);
            }
        }
    }

    for (int v : verts) {
        if (!seen[v]) {
            throw std::invalid_argument("tree must be connected");
        }
    }

    for (const auto& terminal : terminals_) {
        if (terminal.vertex == root_) {
            throw std::invalid_argument("root cannot be a terminal");
        }
    }
}

std::vector<int> Instance::parent_array() const {
    validate();
    std::vector<int> parent(adjacency_.size(), -2);
    std::queue<int> queue;
    parent[root_] = -1;
    queue.push(root_);
    while (!queue.empty()) {
        int u = queue.front();
        queue.pop();
        for (const auto& edge : adjacency_[u]) {
            if (parent[edge.to] == -2) {
                parent[edge.to] = u;
                queue.push(edge.to);
            }
        }
    }
    return parent;
}

std::vector<double> Instance::distance_from_root() const {
    validate();
    std::vector<double> dist(adjacency_.size(), std::numeric_limits<double>::infinity());
    std::queue<int> queue;
    dist[root_] = 0.0;
    queue.push(root_);
    while (!queue.empty()) {
        int u = queue.front();
        queue.pop();
        for (const auto& edge : adjacency_[u]) {
            if (!std::isfinite(dist[edge.to])) {
                dist[edge.to] = dist[u] + edge.weight;
                queue.push(edge.to);
            }
        }
    }
    return dist;
}

double Instance::steiner_cost_for_terminal_subset(const std::vector<int>& terminal_vertices) const {
    validate();
    if (terminal_vertices.empty()) {
        return 0.0;
    }

    auto parent = parent_array();
    std::set<std::pair<int, int>> used_edges;
    for (int terminal : terminal_vertices) {
        if (!is_terminal(terminal)) {
            throw std::invalid_argument("subset contains a non-terminal");
        }
        int current = terminal;
        while (current != root_) {
            int p = parent[current];
            if (p < 0) {
                throw std::logic_error("failed to recover root path");
            }
            int a = std::min(current, p);
            int b = std::max(current, p);
            used_edges.emplace(a, b);
            current = p;
        }
    }

    double cost = 0.0;
    for (const auto& [a, b] : used_edges) {
        bool found = false;
        for (const auto& edge : adjacency_[a]) {
            if (edge.to == b) {
                cost += edge.weight;
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::logic_error("missing edge in steiner computation");
        }
    }
    return 2.0 * cost;
}

Instance Instance::parse(std::istream& input) {
    int n = 0;
    int root = -1;
    if (!(input >> n >> root)) {
        throw std::invalid_argument("failed to read header");
    }
    if (n <= 0) {
        throw std::invalid_argument("vertex count must be positive");
    }

    Instance instance(root);
    for (int i = 0; i < n; ++i) {
        instance.ensure_vertex(i);
    }

    for (int i = 0; i < n - 1; ++i) {
        int u = 0;
        int v = 0;
        double w = 0.0;
        if (!(input >> u >> v >> w)) {
            throw std::invalid_argument("failed to read edge");
        }
        instance.add_edge(u, v, w);
    }

    int terminal_count = 0;
    if (!(input >> terminal_count)) {
        throw std::invalid_argument("failed to read terminal count");
    }

    for (int i = 0; i < terminal_count; ++i) {
        int vertex = 0;
        double demand = 0.0;
        if (!(input >> vertex >> demand)) {
            throw std::invalid_argument("failed to read terminal");
        }
        instance.add_terminal(vertex, demand);
    }

    instance.validate();
    return instance;
}

Instance Instance::parse_file(const std::string& path) {
    std::ifstream stream(path);
    if (!stream) {
        throw std::runtime_error("failed to open instance file: " + path);
    }
    return parse(stream);
}

}  // namespace tucvrp
