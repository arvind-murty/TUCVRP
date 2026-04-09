#pragma once

#include <istream>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

namespace tucvrp {

struct Edge {
    int to;
    double weight;
};

struct Terminal {
    int vertex;
    double demand;
};

class Instance {
  public:
    Instance() = default;
    explicit Instance(int depot);

    // Set the depot vertex for the instance.
    void set_depot(int depot);
    // Add an undirected weighted tree edge.
    void add_edge(int u, int v, double weight);
    // Register a customer vertex with demand in (0, 1].
    void add_terminal(int vertex, double demand);

    // Return the depot vertex id.
    [[nodiscard]] int depot() const;
    // Return the size of the backing vertex array.
    [[nodiscard]] int vertex_count() const;
    // Return the number of undirected edges in the tree.
    [[nodiscard]] int edge_count() const;
    // Return the number of terminals/customers.
    [[nodiscard]] int terminal_count() const;
    // Return the sum of all customer demands.
    [[nodiscard]] double total_demand() const;
    // Return whether a vertex is a customer terminal.
    [[nodiscard]] bool is_terminal(int vertex) const;
    // Return a terminal's demand, or 0 if it is not a terminal.
    [[nodiscard]] double demand_of(int vertex) const;
    // Return the adjacency list for one vertex.
    [[nodiscard]] const std::vector<Edge>& neighbors(int vertex) const;
    // Return all terminals in insertion order.
    [[nodiscard]] const std::vector<Terminal>& terminals() const;
    // Return the active vertices present in the instance.
    [[nodiscard]] std::vector<int> vertices() const;

    // Check that the graph is a connected rooted tree with valid terminals.
    void validate() const;

    // Return the rooted parent of every vertex after rooting the tree at the depot.
    [[nodiscard]] std::vector<int> parent_array() const;
    // Return the weighted distance from the depot to every vertex.
    [[nodiscard]] std::vector<double> distances_from_depot() const;
    // Return the depot distances for terminals
    [[nodiscard]] std::unordered_map<int, double> terminal_distances() const;
    // Return the number of terminals in the subtree rooted at each vertex id.
    [[nodiscard]] std::vector<int> subtree_terminal_counts() const;
    // Return the minimum round-trip tour cost to visit the given terminals from the depot.
    [[nodiscard]] double tour_cost_for_terminals(const std::vector<int>& terminal_vertices) const;

    // Parse an instance from the project's simple text format.
    static Instance parse(std::istream& input);
    // Parse an instance from a file on disk.
    static Instance parse_file(const std::string& path);
    // Copy the tree structure from another instance, but replace the terminal set entirely.
    static Instance with_terminals(const Instance& other, const std::vector<Terminal>& terminals);

  private:
    int depot_ = -1;
    int max_vertex_id_ = -1;
    std::vector<std::vector<Edge>> adjacency_;
    std::vector<Terminal> terminals_;
    std::unordered_map<int, double> demand_by_vertex_;

    // Grow internal storage so a vertex id can be used safely.
    void ensure_vertex(int vertex);
};

// Print a rooted view of the tree together with terminal demands and depot distances.
std::ostream& operator<<(std::ostream& out, const Instance& instance);

}  // namespace tucvrp
