#include "tucvrp/instance.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <sstream>
#include <stdexcept>

namespace tucvrp {

namespace {

constexpr double kEps = 1e-9;

}

// Create an instance whose depot is already known.
Instance::Instance(int depot) : depot_(depot), max_vertex_id_(depot) {
    ensure_vertex(depot);
}

// Copy an existing instance, preserving tree structure and terminal data.
Instance::Instance(const Instance& other, const std::vector<Terminal>& terminals)
    : depot_(other.depot_),
      max_vertex_id_(other.max_vertex_id_),
      adjacency_(other.adjacency_),
      terminals_(terminals),
      demand_by_vertex_(other.demand_by_vertex_) {}

// Expand storage so the given vertex id exists in the adjacency list.
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

// Update the depot vertex.
void Instance::set_depot(int depot) {
    ensure_vertex(depot);
    depot_ = depot;
}

// Insert one undirected weighted edge into the tree representation.
void Instance::add_edge(int u, int v, double weight) {
    if (weight < 0.0) {
        throw std::invalid_argument("edge weight must be non-negative");
    }
    ensure_vertex(u);
    ensure_vertex(v);
    adjacency_[u].push_back(Edge{v, weight});
    adjacency_[v].push_back(Edge{u, weight});
}

// Register one customer terminal and its unit-capacity-normalized demand.
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

// Return the depot vertex id.
int Instance::depot() const { return depot_; }

// Return the size of the internal adjacency array.
int Instance::vertex_count() const { return static_cast<int>(adjacency_.size()); }

// Count undirected edges by halving the adjacency-list total.
int Instance::edge_count() const {
    int total = 0;
    for (const auto& edges : adjacency_) {
        total += static_cast<int>(edges.size());
    }
    return total / 2;
}

// Return how many customer terminals are present.
int Instance::terminal_count() const { return static_cast<int>(terminals_.size()); }

// Sum all terminal demands.
double Instance::total_demand() const {
    double total = 0.0;
    for (const auto& terminal : terminals_) {
        total += terminal.demand;
    }
    return total;
}

// Check whether a vertex is marked as a terminal.
bool Instance::is_terminal(int vertex) const { return demand_by_vertex_.contains(vertex); }

// Look up the demand at a vertex, returning zero for non-terminals.
double Instance::demand_of(int vertex) const {
    auto it = demand_by_vertex_.find(vertex);
    return it == demand_by_vertex_.end() ? 0.0 : it->second;
}

// Return the adjacency list for a valid vertex id.
const std::vector<Edge>& Instance::neighbors(int vertex) const {
    if (vertex < 0 || vertex >= static_cast<int>(adjacency_.size())) {
        throw std::out_of_range("vertex out of range");
    }
    return adjacency_[vertex];
}

// Return the stored terminal list.
const std::vector<Terminal>& Instance::terminals() const { return terminals_; }

// Return the set of vertices that actually participate in the instance.
std::vector<int> Instance::vertices() const {
    std::vector<int> out;
    out.reserve(adjacency_.size());
    for (int v = 0; v < static_cast<int>(adjacency_.size()); ++v) {
        if (!adjacency_[v].empty() || v == depot_ || demand_by_vertex_.contains(v)) {
            out.push_back(v);
        }
    }
    return out;
}

// Verify that the graph is a connected tree rooted at a non-terminal depot.
void Instance::validate() const {
    if (depot_ < 0) {
        throw std::invalid_argument("depot must be set");
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
    queue.push(depot_);
    seen[depot_] = true;

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
        if (terminal.vertex == depot_) {
            throw std::invalid_argument("depot cannot be a terminal");
        }
    }
}

// Root the tree at the depot and record each vertex's parent.
std::vector<int> Instance::parent_array() const {
    validate();
    std::vector<int> parent(adjacency_.size(), -2);
    std::queue<int> queue;
    parent[depot_] = -1;
    queue.push(depot_);
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

// Compute weighted depot-to-vertex distances throughout the tree.
std::vector<double> Instance::distance_from_depot() const {
    validate();
    std::vector<double> dist(adjacency_.size(), std::numeric_limits<double>::infinity());
    std::queue<int> queue;
    dist[depot_] = 0.0;
    queue.push(depot_);
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

// Return depot-to-terminal distances
std::unordered_map<int, double> Instance::terminal_distances() const {
    const auto dist = distance_from_depot();
    std::unordered_map<int, double> terminal_dist;
    for (const auto& terminal : terminals_) {
        terminal_dist[terminal.vertex] = dist[terminal.vertex];
    }
    return terminal_dist;
}

// Compute the minimum round-trip tour cost to visit a given terminal subset from the depot.
double Instance::tour_cost_for_terminals(const std::vector<int>& terminal_vertices) const {
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
        while (current != depot_) {
            int p = parent[current];
            if (p < 0) {
                throw std::logic_error("failed to recover depot path");
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
            throw std::logic_error("missing edge in tour-cost computation");
        }
    }
    return 2.0 * cost;
}

// Parse the project's simple text instance format from an input stream.
Instance Instance::parse(std::istream& input) {
    int n = 0;
    int depot = -1;
    if (!(input >> n >> depot)) {
        throw std::invalid_argument("failed to read header");
    }
    if (n <= 0) {
        throw std::invalid_argument("vertex count must be positive");
    }

    Instance instance(depot);
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

// Open a file and parse an instance from it.
Instance Instance::parse_file(const std::string& path) {
    std::ifstream stream(path);
    if (!stream) {
        throw std::runtime_error("failed to open instance file: " + path);
    }
    return parse(stream);
}

std::ostream& operator<<(std::ostream& out, const Instance& instance) {
    instance.validate();

    const auto parent = instance.parent_array();
    const auto dist = instance.distance_from_depot();
    std::vector<std::vector<Edge>> children(instance.vertex_count());
    for (int u : instance.vertices()) {
        for (const auto& edge : instance.neighbors(u)) {
            if (edge.to != parent[u]) {
                children[u].push_back(edge);
            }
        }
    }

    out << "depot: " << instance.depot() << '\n';
    out << "vertices: " << instance.vertices().size() << ", edges: " << instance.edge_count() << ", terminals: "
        << instance.terminal_count() << '\n';
    out << "tree:\n";

    std::vector<std::pair<int, int>> stack;
    stack.emplace_back(instance.depot(), 0);
    while (!stack.empty()) {
        const auto [u, depth] = stack.back();
        stack.pop_back();

        out << std::string(depth * 2, ' ') << "- v" << u;
        if (u == instance.depot()) {
            out << " [depot]";
        }
        if (instance.is_terminal(u)) {
            out << " [term=" << std::fixed << std::setprecision(3) << instance.demand_of(u) << "]";
        }
        out << " [dep dst=" << std::fixed << std::setprecision(3) << dist[u] << "]";
        if (parent[u] >= 0) {
            double parent_weight = -1.0;
            for (const auto& edge : instance.neighbors(u)) {
                if (edge.to == parent[u]) {
                    parent_weight = edge.weight;
                    break;
                }
            }
            out << " [dst=" << std::fixed << std::setprecision(3) << parent_weight << "]";
        }
        out << '\n';

        for (auto it = children[u].rbegin(); it != children[u].rend(); ++it) {
            stack.emplace_back(it->to, depth + 1);
        }
    }

    return out;
}

}  // namespace tucvrp
