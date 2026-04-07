#pragma once

#include <istream>
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
    explicit Instance(int root);

    void set_root(int root);
    void add_edge(int u, int v, double weight);
    void add_terminal(int vertex, double demand);

    [[nodiscard]] int root() const;
    [[nodiscard]] int vertex_count() const;
    [[nodiscard]] int edge_count() const;
    [[nodiscard]] int terminal_count() const;
    [[nodiscard]] double total_demand() const;
    [[nodiscard]] bool is_terminal(int vertex) const;
    [[nodiscard]] double demand_of(int vertex) const;
    [[nodiscard]] const std::vector<Edge>& neighbors(int vertex) const;
    [[nodiscard]] const std::vector<Terminal>& terminals() const;
    [[nodiscard]] std::vector<int> vertices() const;

    void validate() const;

    [[nodiscard]] std::vector<int> parent_array() const;
    [[nodiscard]] std::vector<double> distance_from_root() const;
    [[nodiscard]] double steiner_cost_for_terminal_subset(const std::vector<int>& terminal_vertices) const;

    static Instance parse(std::istream& input);
    static Instance parse_file(const std::string& path);

  private:
    int root_ = -1;
    int max_vertex_id_ = -1;
    std::vector<std::vector<Edge>> adjacency_;
    std::vector<Terminal> terminals_;
    std::unordered_map<int, double> demand_by_vertex_;

    void ensure_vertex(int vertex);
};

}  // namespace tucvrp
