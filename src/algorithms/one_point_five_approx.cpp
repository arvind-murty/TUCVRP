#include "tucvrp/algorithms/one_point_five_approx.hpp"

#include "tucvrp/preprocessing.hpp"
#include "tucvrp/rng.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace tucvrp {

namespace {

constexpr double kTolerance = 1e-9;

// Internal witness for one subtour while reconstructing the DP choices.
// `demand_bound` tracks the configuration value attached to this subtour in the current DP table,
// while `demand` stores the true demand of the terminals served by the subtour.
struct ReconstructedSubtour {
    std::vector<int> terminals;
    double demand = 0.0;
    double cost = 0.0;
    double demand_bound = 0.0;
};

// Local subtours additionally remember whether they are ending or passing inside one component.
struct LocalReconstructedSubtour {
    ReconstructedSubtour subtour;
    LocalSubtourType type = LocalSubtourType::Ending;
};

// Reconstruction proceeds top-down from one chosen optimum table entry. Caches avoid rebuilding
// the same subtree witness multiple times when a parent explores several combination candidates.
struct ReconstructionState {
    const BoundedDistanceContext& context;
    const LocalPhaseState& local_phase;
    const SubtreePhaseState& subtree_phase;
    OnePointFiveApproxParams params;
    std::vector<std::vector<std::optional<std::vector<ReconstructedSubtour>>>> component_cache;
    std::vector<std::vector<std::optional<std::vector<ReconstructedSubtour>>>> critical_cache;
};

int floor_div(int a, int b) {
    int q = a / b;
    int r = a % b;
    if (r != 0 && ((r > 0) != (b > 0))) {
        --q;
    }
    return q;
}

// Validate epsilon and snap it to the reciprocal grid used by the current implementation.
OnePointFiveApproxParams validate_params(const OnePointFiveApproxParams& params) {
    if (params.epsilon <= 0.0 || params.epsilon >= 1.0) {
        throw std::invalid_argument("OnePointFiveApproxSolver requires epsilon in (0, 1)");
    }

    OnePointFiveApproxParams normalized = params;
    normalized.epsilon = 1.0 / std::ceil(1.0 / normalized.epsilon);
    return normalized;
}

double alpha_value(double epsilon) {
    return std::pow(epsilon, 1.0 / epsilon + 1.0);
}

double gamma_zero_value(double epsilon) {
    const double gamma = 12.0 / epsilon;
    const double alpha = alpha_value(epsilon);
    return epsilon * alpha / gamma;
}

bool same_double(double a, double b) {
    return std::abs(a - b) <= kTolerance;
}

bool local_entry_less(const LocalConfigurationEntry& a, const LocalConfigurationEntry& b) {
    if (!same_double(a.demand_bound, b.demand_bound)) {
        return a.demand_bound < b.demand_bound;
    }
    return static_cast<int>(a.type) < static_cast<int>(b.type);
}

void normalize_configuration(LocalConfiguration& configuration) {
    std::sort(configuration.entries.begin(), configuration.entries.end(), local_entry_less);
}

bool same_configuration(const LocalConfiguration& a, const LocalConfiguration& b) {
    if (a.entries.size() != b.entries.size()) {
        return false;
    }
    for (std::size_t i = 0; i < a.entries.size(); ++i) {
        if (!same_double(a.entries[i].demand_bound, b.entries[i].demand_bound) ||
            a.entries[i].type != b.entries[i].type) {
            return false;
        }
    }
    return true;
}

bool subtree_entry_less(const SubtreeConfigurationEntry& a, const SubtreeConfigurationEntry& b) {
    if (!same_double(a.demand, b.demand)) {
        return a.demand < b.demand;
    }
    return a.multiplicity < b.multiplicity;
}

void normalize_subtree_configuration(SubtreeConfiguration& configuration) {
    std::sort(configuration.entries.begin(), configuration.entries.end(), subtree_entry_less);
    std::vector<SubtreeConfigurationEntry> merged;
    for (const auto& entry : configuration.entries) {
        if (entry.multiplicity == 0) {
            continue;
        }
        if (!merged.empty() && same_double(merged.back().demand, entry.demand)) {
            merged.back().multiplicity += entry.multiplicity;
        } else {
            merged.push_back(entry);
        }
    }
    configuration.entries = std::move(merged);
}

bool same_subtree_configuration(const SubtreeConfiguration& a, const SubtreeConfiguration& b) {
    if (a.entries.size() != b.entries.size()) {
        return false;
    }
    for (std::size_t i = 0; i < a.entries.size(); ++i) {
        if (!same_double(a.entries[i].demand, b.entries[i].demand) ||
            a.entries[i].multiplicity != b.entries[i].multiplicity) {
            return false;
        }
    }
    return true;
}

int subtree_configuration_size(const SubtreeConfiguration& configuration) {
    int total = 0;
    for (const auto& entry : configuration.entries) {
        total += entry.multiplicity;
    }
    return total;
}

double total_reconstructed_cost(const std::vector<ReconstructedSubtour>& subtours) {
    double total = 0.0;
    for (const auto& subtour : subtours) {
        total += subtour.cost;
    }
    return total;
}

void normalize_terminals(std::vector<int>& terminals) {
    std::sort(terminals.begin(), terminals.end());
}

SubtreeConfiguration subtree_configuration_from_reconstructed_tours(
    const std::vector<ReconstructedSubtour>& subtours) {
    SubtreeConfiguration configuration;
    configuration.entries.reserve(subtours.size());
    for (const auto& subtour : subtours) {
        configuration.entries.push_back(SubtreeConfigurationEntry{
            .demand = subtour.demand_bound,
            .multiplicity = 1,
        });
    }
    normalize_subtree_configuration(configuration);
    return configuration;
}

Tour to_public_tour(const ReconstructedSubtour& subtour) {
    Tour tour;
    tour.terminals = subtour.terminals;
    normalize_terminals(tour.terminals);
    tour.demand = subtour.demand;
    tour.cost = subtour.cost;
    return tour;
}

SolveResult to_public_solution(const std::vector<ReconstructedSubtour>& subtours) {
    SolveResult result;
    result.tours.reserve(subtours.size());
    for (const auto& subtour : subtours) {
        result.tours.push_back(to_public_tour(subtour));
        result.cost += subtour.cost;
    }
    return result;
}

int map_normalized_terminal_to_original(const Instance& original,
                                        const Instance& normalized,
                                        const std::vector<int>& normalized_parent,
                                        int terminal) {
    if (!normalized.is_terminal(terminal)) {
        throw std::invalid_argument("projection received a non-terminal from the normalized instance");
    }

    int current = terminal;
    while (current >= 0) {
        if (current < original.vertex_count() && original.is_terminal(current)) {
            return current;
        }

        const int parent = normalized_parent[current];
        if (parent < 0) {
            break;
        }

        double parent_weight = -1.0;
        for (const auto& edge : normalized.neighbors(current)) {
            if (edge.to == parent) {
                parent_weight = edge.weight;
                break;
            }
        }
        if (parent_weight < 0.0) {
            throw std::logic_error("failed to locate parent edge during terminal projection");
        }
        if (parent_weight > kTolerance) {
            break;
        }
        current = parent;
    }

    throw std::logic_error("failed to map normalized terminal back to an original terminal");
}

SolveResult project_normalized_solution_to_original(const Instance& original,
                                                    const Instance& normalized,
                                                    const SolveResult& normalized_solution) {
    const std::vector<int> normalized_parent = normalized.parent_array();

    SolveResult projected;
    for (const auto& normalized_tour : normalized_solution.tours) {
        Tour projected_tour;
        projected_tour.terminals.reserve(normalized_tour.terminals.size());
        for (const int terminal : normalized_tour.terminals) {
            projected_tour.terminals.push_back(
                map_normalized_terminal_to_original(original, normalized, normalized_parent, terminal));
        }
        normalize_terminals(projected_tour.terminals);
        projected_tour.terminals.erase(
            std::unique(projected_tour.terminals.begin(), projected_tour.terminals.end()),
            projected_tour.terminals.end());

        for (const int terminal : projected_tour.terminals) {
            projected_tour.demand += original.demand_of(terminal);
        }
        projected_tour.cost = original.tour_cost_for_terminals(projected_tour.terminals);
        projected_tour.walk = original.tour_walk_for_terminals(projected_tour.terminals);

        if (!same_double(projected_tour.cost, normalized_tour.cost)) {
            throw std::logic_error("normalized-tour projection changed the route cost");
        }

        projected.cost += projected_tour.cost;
        projected.tours.push_back(std::move(projected_tour));
    }

    return projected;
}

bool reconstruct_component_root_value(ReconstructionState& state,
                                      int component_id,
                                      int value_index,
                                      std::vector<ReconstructedSubtour>& out);
bool reconstruct_critical_vertex_value(ReconstructionState& state,
                                       int critical_vertex,
                                       int value_index,
                                       std::vector<ReconstructedSubtour>& out);

const SubtreeConfigurationTable* critical_table_for_vertex(const SubtreePhaseState& subtree_phase, int vertex) {
    if (vertex < 0 || vertex >= static_cast<int>(subtree_phase.critical_vertex_tables.size())) {
        return nullptr;
    }
    if (subtree_phase.critical_vertex_tables[vertex].vertex != vertex) {
        return nullptr;
    }
    return &subtree_phase.critical_vertex_tables[vertex];
}

int component_id_by_root_vertex(const TreeDecomposition& decomposition, int root_vertex) {
    for (const Component& component : decomposition.components) {
        if (component.root == root_vertex) {
            return component.id;
        }
    }
    return -1;
}

bool reconstruct_subtree_value_at_vertex(ReconstructionState& state,
                                         int vertex,
                                         const SubtreeConfiguration& configuration,
                                         double cost,
                                         std::vector<ReconstructedSubtour>& out) {
    if (const auto* critical_table = critical_table_for_vertex(state.subtree_phase, vertex)) {
        for (int i = 0; i < static_cast<int>(critical_table->values.size()); ++i) {
            if (same_subtree_configuration(critical_table->values[i].configuration, configuration) &&
                same_double(critical_table->values[i].cost, cost)) {
                return reconstruct_critical_vertex_value(state, vertex, i, out);
            }
        }
    }

    const int component_id = component_id_by_root_vertex(state.context.decomposition, vertex);
    if (component_id != -1 &&
        component_id < static_cast<int>(state.subtree_phase.component_root_tables.size()) &&
        state.subtree_phase.component_root_tables[component_id].vertex == vertex) {
        const auto& table = state.subtree_phase.component_root_tables[component_id];
        for (int i = 0; i < static_cast<int>(table.values.size()); ++i) {
            if (same_subtree_configuration(table.values[i].configuration, configuration) &&
                same_double(table.values[i].cost, cost)) {
                return reconstruct_component_root_value(state, component_id, i, out);
            }
        }
    }

    return false;
}

double beta_value(double epsilon) {
    return 0.25 * std::pow(epsilon, 4.0 / epsilon + 1.0);
}

double component_spine_cost(const BoundedDistanceContext& context, const Component& component) {
    if (component.exit == -1) {
        return 0.0;
    }
    return 2.0 * (context.rooted_tree.distances_from_depot[component.exit] -
                  context.rooted_tree.distances_from_depot[component.root]);
}

double edge_cost_in_height_reduced_tree(const BoundedDistanceContext& context,
                                        int critical_vertex,
                                        int child_root) {
    const double critical_distance = context.rooted_tree.distances_from_depot[critical_vertex];
    const double child_distance = context.rooted_tree.distances_from_depot[child_root];
    if (child_distance + kTolerance < critical_distance) {
        throw std::logic_error("height-reduced child root lies above its critical vertex");
    }
    return child_distance - critical_distance;
}

void record_subtree_value(SubtreeConfigurationTable& table,
                          SubtreeConfiguration configuration,
                          double cost) {
    normalize_subtree_configuration(configuration);
    for (auto& value : table.values) {
        if (same_subtree_configuration(value.configuration, configuration)) {
            value.cost = std::min(value.cost, cost);
            return;
        }
    }
    table.values.push_back(SubtreeConfigurationValue{
        .configuration = std::move(configuration),
        .cost = cost,
    });
}

std::vector<bool> component_membership(const RootedTreeData& rooted_tree, const Component& component) {
    std::vector<bool> in_component(rooted_tree.parent.size(), false);
    for (const int v : component.vertices) {
        in_component[v] = true;
    }
    return in_component;
}

// Compute the cost of the minimum subtour in component c spanning the component root together
// with the requested terminals, and the exit vertex when the subtour is passing.
double local_subtour_cost(const RootedTreeData& rooted_tree,
                          const Component& component,
                          const std::vector<bool>& in_component,
                          const std::vector<int>& terminals,
                          LocalSubtourType type) {
    std::vector<bool> used_edge_to_parent(rooted_tree.parent.size(), false);
    auto mark_path_to_root = [&](int vertex) {
        int current = vertex;
        while (current != component.root) {
            if (current == -1 || !in_component[current]) {
                throw std::logic_error("local subtour leaves the component");
            }
            const int parent = rooted_tree.parent[current];
            if (parent == -1 || !in_component[parent]) {
                throw std::logic_error("component root is not an ancestor of a local subtour vertex");
            }
            used_edge_to_parent[current] = true;
            current = parent;
        }
    };

    for (const int terminal : terminals) {
        mark_path_to_root(terminal);
    }
    if (type == LocalSubtourType::Passing) {
        if (component.exit == -1) {
            throw std::logic_error("passing local subtour requested in a leaf component");
        }
        mark_path_to_root(component.exit);
    }

    double cost = 0.0;
    for (int v = 0; v < static_cast<int>(used_edge_to_parent.size()); ++v) {
        if (!used_edge_to_parent[v]) {
            continue;
        }
        const int parent = rooted_tree.parent[v];
        cost += 2.0 * (rooted_tree.distances_from_depot[v] - rooted_tree.distances_from_depot[parent]);
    }
    return cost;
}

std::vector<LocalPart> component_parts(const BoundedDistanceContext& context,
                                       int component_id,
                                       double gamma_zero) {
    const Component& component = context.decomposition.components[component_id];
    const RootedTreeData& rooted_tree = context.rooted_tree;

    std::vector<LocalPart> parts;

    // Definition 26: every big terminal becomes a singleton part in Q_c.
    for (const int v : component.vertices) {
        if (rooted_tree.is_terminal(v) && rooted_tree.demands[v] > gamma_zero) {
            parts.push_back(LocalPart{
                .terminals = {v},
                .demand = rooted_tree.demands[v],
            });
        }
    }

    // Definition 26: every cell contributes one part consisting of all small terminals in that cell.
    for (const int block_id : component.block_ids) {
        const Block& block = context.decomposition.blocks[block_id];
        for (const int cluster_id : block.cluster_ids) {
            const Cluster& cluster = context.decomposition.clusters[cluster_id];
            for (const int cell_id : cluster.cell_ids) {
                const Cell& cell = context.decomposition.cells[cell_id];
                LocalPart part;
                for (const int v : cell.vertices) {
                    if (rooted_tree.is_terminal(v) && rooted_tree.demands[v] <= gamma_zero) {
                        part.terminals.push_back(v);
                        part.demand += rooted_tree.demands[v];
                    }
                }
                if (!part.terminals.empty()) {
                    parts.push_back(std::move(part));
                }
            }
        }
    }

    return parts;
}

std::vector<double> y_values_for_component(const std::vector<LocalPart>& parts, double alpha) {
    std::vector<double> values{alpha};
    const int part_count = static_cast<int>(parts.size());
    const int subset_count = 1 << part_count;
    for (int mask = 1; mask < subset_count; ++mask) {
        double demand = 0.0;
        for (int i = 0; i < part_count; ++i) {
            if (mask & (1 << i)) {
                demand += parts[i].demand;
            }
        }
        if (demand > alpha + kTolerance && demand <= 1.0 + kTolerance) {
            values.push_back(std::min(1.0, demand));
        }
    }

    std::sort(values.begin(), values.end());
    values.erase(
        std::unique(values.begin(), values.end(), [](double a, double b) { return same_double(a, b); }),
        values.end());
    return values;
}

std::vector<LocalConfiguration> enumerate_local_configurations(const std::vector<double>& y_values,
                                                               bool allow_passing,
                                                               int max_length) {
    std::vector<LocalConfiguration> configurations;
    if (max_length == 0) {
        configurations.push_back(LocalConfiguration{});
        return configurations;
    }

    std::vector<LocalConfigurationEntry> options;
    options.reserve(y_values.size() * (allow_passing ? 2 : 1));
    for (const double y : y_values) {
        options.push_back(LocalConfigurationEntry{.demand_bound = y, .type = LocalSubtourType::Ending});
        if (allow_passing) {
            options.push_back(LocalConfigurationEntry{.demand_bound = y, .type = LocalSubtourType::Passing});
        }
    }

    for (int length = 1; length <= max_length; ++length) {
        std::vector<int> digits(length, 0);
        bool done = false;
        while (!done) {
            LocalConfiguration current;
            current.entries.reserve(length);
            for (const int digit : digits) {
                current.entries.push_back(options[digit]);
            }

            normalize_configuration(current);
            bool duplicate = false;
            for (const LocalConfiguration& existing : configurations) {
                if (same_configuration(existing, current)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                configurations.push_back(std::move(current));
            }

            int position = length - 1;
            while (position >= 0) {
                ++digits[position];
                if (digits[position] < static_cast<int>(options.size())) {
                    break;
                }
                digits[position] = 0;
                --position;
            }
            done = (position < 0);
        }
    }
    return configurations;
}

double configuration_cost(const BoundedDistanceContext& context,
                          int component_id,
                          const std::vector<LocalPart>& parts,
                          const LocalConfiguration& configuration) {
    const Component& component = context.decomposition.components[component_id];
    const RootedTreeData& rooted_tree = context.rooted_tree;
    const std::vector<bool> in_component = component_membership(rooted_tree, component);
    const int part_count = static_cast<int>(parts.size());
    const int subtour_count = static_cast<int>(configuration.entries.size());

    if (part_count == 0) {
        return subtour_count == 0 ? 0.0 : std::numeric_limits<double>::infinity();
    }
    if (subtour_count == 0) {
        return std::numeric_limits<double>::infinity();
    }

    double best_cost = std::numeric_limits<double>::infinity();
    std::vector<double> assigned_demand(subtour_count, 0.0);
    std::vector<std::vector<int>> assigned_terminals(subtour_count);
    std::vector<int> assigned_part_count(subtour_count, 0);

    // Algorithm 2 line 3: brute-force every partition of Q_c into ℓ(A) labeled parts.
    std::vector<int> next_choice(part_count, 0);
    std::vector<int> chosen_bucket(part_count, -1);
    int part_index = 0;
    while (part_index >= 0) {
        if (part_index == part_count) {
            bool all_nonempty = true;
            for (int i = 0; i < subtour_count; ++i) {
                if (assigned_part_count[i] == 0) {
                    all_nonempty = false;
                    break;
                }
            }
            if (all_nonempty) {
                double total_cost = 0.0;
                for (int i = 0; i < subtour_count; ++i) {
                    total_cost += local_subtour_cost(rooted_tree,
                                                     component,
                                                     in_component,
                                                     assigned_terminals[i],
                                                     configuration.entries[i].type);
                }
                best_cost = std::min(best_cost, total_cost);
            }

            --part_index;
            if (part_index >= 0) {
                const int bucket = chosen_bucket[part_index];
                assigned_terminals[bucket].resize(assigned_terminals[bucket].size() -
                                                  parts[part_index].terminals.size());
                assigned_part_count[bucket] -= 1;
                assigned_demand[bucket] -= parts[part_index].demand;
            }
            continue;
        }

        bool advanced = false;
        while (next_choice[part_index] < subtour_count) {
            const int bucket = next_choice[part_index]++;
            const double next_demand = assigned_demand[bucket] + parts[part_index].demand;
            if (next_demand > configuration.entries[bucket].demand_bound + kTolerance) {
                continue;
            }

            chosen_bucket[part_index] = bucket;
            assigned_demand[bucket] = next_demand;
            assigned_part_count[bucket] += 1;
            assigned_terminals[bucket].insert(assigned_terminals[bucket].end(),
                                              parts[part_index].terminals.begin(),
                                              parts[part_index].terminals.end());

            ++part_index;
            if (part_index < part_count) {
                next_choice[part_index] = 0;
            }
            advanced = true;
            break;
        }

        if (!advanced) {
            next_choice[part_index] = 0;
            chosen_bucket[part_index] = -1;
            --part_index;
            if (part_index >= 0) {
                const int bucket = chosen_bucket[part_index];
                assigned_terminals[bucket].resize(assigned_terminals[bucket].size() -
                                                  parts[part_index].terminals.size());
                assigned_part_count[bucket] -= 1;
                assigned_demand[bucket] -= parts[part_index].demand;
            }
        }
    }
    return best_cost;
}

std::vector<SubtreeConfiguration> combine_exit_and_local_configurations(
    const SubtreeConfiguration& exit_configuration,
    const LocalConfiguration& local_configuration) {
    std::vector<int> passing_indices;
    for (int i = 0; i < static_cast<int>(local_configuration.entries.size()); ++i) {
        if (local_configuration.entries[i].type == LocalSubtourType::Passing) {
            passing_indices.push_back(i);
        }
    }

    if (passing_indices.empty()) {
        SubtreeConfiguration result;
        result.entries = exit_configuration.entries;
        for (const auto& entry : local_configuration.entries) {
            result.entries.push_back(SubtreeConfigurationEntry{
                .demand = entry.demand_bound,
                .multiplicity = 1,
            });
        }
        normalize_subtree_configuration(result);
        return {result};
    }

    if (subtree_configuration_size(exit_configuration) < static_cast<int>(passing_indices.size())) {
        return {};
    }

    std::vector<SubtreeConfiguration> out;
    std::vector<int> association_counts(exit_configuration.entries.size(), 0);
    std::vector<int> choice_by_passing(passing_indices.size(), -1);
    std::vector<int> next_choice(passing_indices.size(), 0);
    int position = 0;
    while (position >= 0) {
        if (position == static_cast<int>(passing_indices.size())) {
            SubtreeConfiguration result;

            for (int i = 0; i < static_cast<int>(passing_indices.size()); ++i) {
                const int local_index = passing_indices[i];
                const int exit_choice = choice_by_passing[i];
                result.entries.push_back(SubtreeConfigurationEntry{
                    .demand = local_configuration.entries[local_index].demand_bound +
                              exit_configuration.entries[exit_choice].demand,
                    .multiplicity = 1,
                });
            }
            for (int j = 0; j < static_cast<int>(exit_configuration.entries.size()); ++j) {
                const int remaining = exit_configuration.entries[j].multiplicity - association_counts[j];
                if (remaining > 0) {
                    result.entries.push_back(SubtreeConfigurationEntry{
                        .demand = exit_configuration.entries[j].demand,
                        .multiplicity = remaining,
                    });
                }
            }
            for (const auto& entry : local_configuration.entries) {
                if (entry.type == LocalSubtourType::Ending) {
                    result.entries.push_back(SubtreeConfigurationEntry{
                        .demand = entry.demand_bound,
                        .multiplicity = 1,
                    });
                }
            }

            normalize_subtree_configuration(result);
            bool duplicate = false;
            for (const auto& existing : out) {
                if (same_subtree_configuration(existing, result)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                out.push_back(std::move(result));
            }

            --position;
            if (position >= 0) {
                association_counts[choice_by_passing[position]] -= 1;
            }
            continue;
        }

        const int local_index = passing_indices[position];
        const double local_demand = local_configuration.entries[local_index].demand_bound;
        bool advanced = false;
        while (next_choice[position] < static_cast<int>(exit_configuration.entries.size())) {
            const int exit_choice = next_choice[position]++;
            if (association_counts[exit_choice] >= exit_configuration.entries[exit_choice].multiplicity) {
                continue;
            }
            if (local_demand + exit_configuration.entries[exit_choice].demand > 1.0 + kTolerance) {
                continue;
            }

            choice_by_passing[position] = exit_choice;
            association_counts[exit_choice] += 1;
            ++position;
            if (position < static_cast<int>(passing_indices.size())) {
                next_choice[position] = 0;
            }
            advanced = true;
            break;
        }

        if (!advanced) {
            next_choice[position] = 0;
            if (position == 0) {
                break;
            }
            --position;
            association_counts[choice_by_passing[position]] -= 1;
        }
    }
    return out;
}

std::vector<int> child_component_ids_of_critical_vertex(const BoundedDistanceContext& context,
                                                        int critical_vertex) {
    std::vector<int> children;
    for (const Component& component : context.decomposition.components) {
        if (context.height_reduced.critical_vertex_by_component[component.id] == critical_vertex &&
            component.root != critical_vertex) {
            children.push_back(component.id);
        }
    }
    std::sort(children.begin(), children.end(), [&](int a, int b) {
        const int root_a = context.decomposition.components[a].root;
        const int root_b = context.decomposition.components[b].root;
        return root_a < root_b;
    });
    return children;
}

std::vector<double> candidate_y_values_from_child_tables(const std::vector<SubtreeConfigurationTable>& child_tables) {
    std::vector<double> y_values;
    for (const auto& table : child_tables) {
        for (const auto& value : table.values) {
            for (const auto& entry : value.configuration.entries) {
                y_values.push_back(entry.demand);
            }
        }
    }
    // Keep 1.0 as a top sentinel so every child-subtour demand can be rounded upward to some x in X.
    y_values.push_back(1.0);
    std::sort(y_values.begin(), y_values.end());
    y_values.erase(
        std::unique(y_values.begin(), y_values.end(), [](double a, double b) { return same_double(a, b); }),
        y_values.end());
    return y_values;
}

std::vector<std::vector<double>> enumerate_x_sets(const std::vector<double>& y_values, int max_size) {
    std::vector<double> base_values;
    for (const double y : y_values) {
        if (!same_double(y, 1.0)) {
            base_values.push_back(y);
        }
    }

    const int capped_size = std::max(0, max_size - 1);
    std::vector<std::vector<double>> subsets{{}};
    for (const double y : base_values) {
        const int current_count = static_cast<int>(subsets.size());
        for (int i = 0; i < current_count; ++i) {
            if (static_cast<int>(subsets[i].size()) >= capped_size) {
                continue;
            }
            auto extended = subsets[i];
            extended.push_back(y);
            subsets.push_back(std::move(extended));
        }
    }
    for (auto& subset : subsets) {
        subset.push_back(1.0);
        std::sort(subset.begin(), subset.end());
    }
    std::sort(subsets.begin(), subsets.end());
    subsets.erase(std::unique(subsets.begin(), subsets.end()), subsets.end());
    subsets.erase(std::remove_if(subsets.begin(),
                                 subsets.end(),
                                 [](const std::vector<double>& subset) { return subset.empty(); }),
                  subsets.end());
    return subsets;
}

SubtreeConfiguration round_subtree_configuration_to_x(const SubtreeConfiguration& configuration,
                                                      const std::vector<double>& x_values) {
    if (x_values.empty()) {
        throw std::invalid_argument("X must be non-empty when rounding a subtree configuration");
    }

    SubtreeConfiguration rounded;
    for (const auto& entry : configuration.entries) {
        auto it = std::lower_bound(x_values.begin(), x_values.end(), entry.demand - kTolerance);
        if (it == x_values.end()) {
            throw std::logic_error("failed to round subtree demand to X");
        }
        rounded.entries.push_back(SubtreeConfigurationEntry{
            .demand = *it,
            .multiplicity = entry.multiplicity,
        });
    }
    normalize_subtree_configuration(rounded);
    return rounded;
}

std::vector<SubtreeConfiguration> combine_sum_lists(const SubtreeConfiguration& previous,
                                                    const SubtreeConfiguration& rounded_child) {
    if (previous.entries.empty()) {
        return {rounded_child};
    }

    std::vector<double> previous_demands;
    for (const auto& entry : previous.entries) {
        for (int i = 0; i < entry.multiplicity; ++i) {
            previous_demands.push_back(entry.demand);
        }
    }

    std::vector<double> child_demands;
    for (const auto& entry : rounded_child.entries) {
        for (int i = 0; i < entry.multiplicity; ++i) {
            child_demands.push_back(entry.demand);
        }
    }

    std::vector<SubtreeConfiguration> out;
    std::vector<bool> used_previous(previous_demands.size(), false);
    std::vector<int> pairing(child_demands.size(), -2);
    std::vector<int> next_choice(child_demands.size(), -1);
    int child_index = 0;
    while (child_index >= 0) {
        if (child_index == static_cast<int>(child_demands.size())) {
            SubtreeConfiguration result;
            for (std::size_t i = 0; i < previous_demands.size(); ++i) {
                if (!used_previous[i]) {
                    result.entries.push_back(SubtreeConfigurationEntry{
                        .demand = previous_demands[i],
                        .multiplicity = 1,
                    });
                }
            }
            for (std::size_t i = 0; i < child_demands.size(); ++i) {
                if (pairing[i] == -1) {
                    result.entries.push_back(SubtreeConfigurationEntry{
                        .demand = child_demands[i],
                        .multiplicity = 1,
                    });
                } else {
                    result.entries.push_back(SubtreeConfigurationEntry{
                        .demand = child_demands[i] + previous_demands[pairing[i]],
                        .multiplicity = 1,
                    });
                }
            }
            normalize_subtree_configuration(result);
            bool duplicate = false;
            for (const auto& existing : out) {
                if (same_subtree_configuration(existing, result)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                out.push_back(std::move(result));
            }

            --child_index;
            if (child_index >= 0 && pairing[child_index] >= 0) {
                used_previous[pairing[child_index]] = false;
            }
            continue;
        }

        bool advanced = false;
        while (next_choice[child_index] <= static_cast<int>(previous_demands.size())) {
            const int choice = next_choice[child_index]++;
            if (choice == -1) {
                pairing[child_index] = -1;
                ++child_index;
                if (child_index < static_cast<int>(child_demands.size())) {
                    next_choice[child_index] = -1;
                }
                advanced = true;
                break;
            }

            if (choice >= static_cast<int>(previous_demands.size())) {
                break;
            }
            if (used_previous[choice]) {
                continue;
            }
            if (previous_demands[choice] + child_demands[child_index] > 1.0 + kTolerance) {
                continue;
            }
            used_previous[choice] = true;
            pairing[child_index] = choice;
            ++child_index;
            if (child_index < static_cast<int>(child_demands.size())) {
                next_choice[child_index] = -1;
            }
            advanced = true;
            break;
        }

        if (!advanced) {
            next_choice[child_index] = -1;
            if (child_index == 0) {
                break;
            }
            --child_index;
            if (pairing[child_index] >= 0) {
                used_previous[pairing[child_index]] = false;
            }
        }
    }
    return out;
}

// Re-run the Algorithm 2 partition search until we find one partition that attains the chosen
// optimum f(c, A). The resulting local witness keeps the exact served terminals and the actual
// subtour cost inside the component.
bool reconstruct_local_configuration(const BoundedDistanceContext& context,
                                     int component_id,
                                     const std::vector<LocalPart>& parts,
                                     const LocalConfiguration& configuration,
                                     double target_cost,
                                     std::vector<LocalReconstructedSubtour>& out) {
    const Component& component = context.decomposition.components[component_id];
    const RootedTreeData& rooted_tree = context.rooted_tree;
    const std::vector<bool> in_component = component_membership(rooted_tree, component);
    const int part_count = static_cast<int>(parts.size());
    const int subtour_count = static_cast<int>(configuration.entries.size());

    if (part_count == 0) {
        if (subtour_count == 0 && same_double(target_cost, 0.0)) {
            out.clear();
            return true;
        }
        return false;
    }

    std::vector<double> assigned_demand(subtour_count, 0.0);
    std::vector<std::vector<int>> assigned_terminals(subtour_count);
    std::vector<int> assigned_part_count(subtour_count, 0);
    std::vector<int> next_choice(part_count, 0);
    std::vector<int> chosen_bucket(part_count, -1);

    int part_index = 0;
    while (part_index >= 0) {
        if (part_index == part_count) {
            bool all_nonempty = true;
            for (int i = 0; i < subtour_count; ++i) {
                if (assigned_part_count[i] == 0) {
                    all_nonempty = false;
                    break;
                }
            }

            if (all_nonempty) {
                std::vector<LocalReconstructedSubtour> candidate;
                double total_cost = 0.0;
                for (int i = 0; i < subtour_count; ++i) {
                    normalize_terminals(assigned_terminals[i]);
                    const double subtour_cost =
                        local_subtour_cost(rooted_tree,
                                           component,
                                           in_component,
                                           assigned_terminals[i],
                                           configuration.entries[i].type);
                    total_cost += subtour_cost;
                    candidate.push_back(LocalReconstructedSubtour{
                        .subtour =
                            ReconstructedSubtour{
                                .terminals = assigned_terminals[i],
                                .demand = assigned_demand[i],
                                .cost = subtour_cost,
                                .demand_bound = configuration.entries[i].demand_bound,
                            },
                        .type = configuration.entries[i].type,
                    });
                }

                if (same_double(total_cost, target_cost)) {
                    out = std::move(candidate);
                    return true;
                }
            }

            --part_index;
            if (part_index >= 0) {
                const int bucket = chosen_bucket[part_index];
                assigned_terminals[bucket].resize(assigned_terminals[bucket].size() -
                                                  parts[part_index].terminals.size());
                assigned_part_count[bucket] -= 1;
                assigned_demand[bucket] -= parts[part_index].demand;
            }
            continue;
        }

        bool advanced = false;
        while (next_choice[part_index] < subtour_count) {
            const int bucket = next_choice[part_index]++;
            const double next_demand = assigned_demand[bucket] + parts[part_index].demand;
            if (next_demand > configuration.entries[bucket].demand_bound + kTolerance) {
                continue;
            }

            chosen_bucket[part_index] = bucket;
            assigned_demand[bucket] = next_demand;
            assigned_part_count[bucket] += 1;
            assigned_terminals[bucket].insert(assigned_terminals[bucket].end(),
                                              parts[part_index].terminals.begin(),
                                              parts[part_index].terminals.end());
            ++part_index;
            if (part_index < part_count) {
                next_choice[part_index] = 0;
            }
            advanced = true;
            break;
        }

        if (!advanced) {
            next_choice[part_index] = 0;
            chosen_bucket[part_index] = -1;
            --part_index;
            if (part_index >= 0) {
                const int bucket = chosen_bucket[part_index];
                assigned_terminals[bucket].resize(assigned_terminals[bucket].size() -
                                                  parts[part_index].terminals.size());
                assigned_part_count[bucket] -= 1;
                assigned_demand[bucket] -= parts[part_index].demand;
            }
        }
    }

    return false;
}

ReconstructedSubtour merge_reconstructed_subtours(const ReconstructedSubtour& a,
                                                  const ReconstructedSubtour& b,
                                                  double merged_bound) {
    ReconstructedSubtour merged;
    merged.terminals = a.terminals;
    merged.terminals.insert(merged.terminals.end(), b.terminals.begin(), b.terminals.end());
    normalize_terminals(merged.terminals);
    merged.demand = a.demand + b.demand;
    merged.cost = a.cost + b.cost;
    merged.demand_bound = merged_bound;
    return merged;
}

bool combine_exit_and_local_witnesses(const std::vector<ReconstructedSubtour>& exit_subtours,
                                      const std::vector<LocalReconstructedSubtour>& local_subtours,
                                      double spine_cost,
                                      const SubtreeConfiguration& target_configuration,
                                      double target_cost,
                                      std::vector<ReconstructedSubtour>& out) {
    std::vector<int> passing_indices;
    for (int i = 0; i < static_cast<int>(local_subtours.size()); ++i) {
        if (local_subtours[i].type == LocalSubtourType::Passing) {
            passing_indices.push_back(i);
        }
    }

    std::vector<int> association_counts(exit_subtours.size(), 0);
    std::vector<int> choice_by_passing(passing_indices.size(), -1);
    std::vector<int> next_choice(passing_indices.size(), 0);
    int passing_position = 0;
    while (passing_position >= 0) {
        if (passing_position == static_cast<int>(passing_indices.size())) {
            std::vector<ReconstructedSubtour> candidate;
            for (int i = 0; i < static_cast<int>(passing_indices.size()); ++i) {
                const auto& local_subtour = local_subtours[passing_indices[i]].subtour;
                const auto& exit_subtour = exit_subtours[choice_by_passing[i]];
                candidate.push_back(merge_reconstructed_subtours(local_subtour,
                                                                 exit_subtour,
                                                                 local_subtour.demand_bound +
                                                                     exit_subtour.demand_bound));
            }
            for (int i = 0; i < static_cast<int>(exit_subtours.size()); ++i) {
                if (association_counts[i] == 0) {
                    ReconstructedSubtour lifted = exit_subtours[i];
                    lifted.cost += spine_cost;
                    candidate.push_back(std::move(lifted));
                }
            }
            for (const auto& local_subtour : local_subtours) {
                if (local_subtour.type == LocalSubtourType::Ending) {
                    candidate.push_back(local_subtour.subtour);
                }
            }

            if (same_subtree_configuration(subtree_configuration_from_reconstructed_tours(candidate),
                                           target_configuration) &&
                same_double(total_reconstructed_cost(candidate), target_cost)) {
                out = std::move(candidate);
                return true;
            }

            --passing_position;
            if (passing_position >= 0) {
                association_counts[choice_by_passing[passing_position]] -= 1;
            }
            continue;
        }

        const auto& local_subtour = local_subtours[passing_indices[passing_position]].subtour;
        bool advanced = false;
        while (next_choice[passing_position] < static_cast<int>(exit_subtours.size())) {
            const int exit_choice = next_choice[passing_position]++;
            if (association_counts[exit_choice] >= 1) {
                continue;
            }
            if (local_subtour.demand_bound + exit_subtours[exit_choice].demand_bound > 1.0 + kTolerance) {
                continue;
            }
            choice_by_passing[passing_position] = exit_choice;
            association_counts[exit_choice] += 1;
            ++passing_position;
            if (passing_position < static_cast<int>(passing_indices.size())) {
                next_choice[passing_position] = 0;
            }
            advanced = true;
            break;
        }

        if (!advanced) {
            next_choice[passing_position] = 0;
            if (passing_position == 0) {
                break;
            }
            --passing_position;
            association_counts[choice_by_passing[passing_position]] -= 1;
        }
    }

    return false;
}

std::vector<std::vector<ReconstructedSubtour>> combine_critical_vertex_witnesses(
    const std::vector<ReconstructedSubtour>& previous,
    const std::vector<ReconstructedSubtour>& child) {
    std::vector<std::vector<ReconstructedSubtour>> out;
    std::vector<bool> used_previous(previous.size(), false);
    std::vector<int> pairing(child.size(), -2);
    std::vector<int> next_choice(child.size(), -1);
    int child_index = 0;
    while (child_index >= 0) {
        if (child_index == static_cast<int>(child.size())) {
            std::vector<ReconstructedSubtour> result;
            for (int i = 0; i < static_cast<int>(previous.size()); ++i) {
                if (!used_previous[i]) {
                    result.push_back(previous[i]);
                }
            }
            for (int i = 0; i < static_cast<int>(child.size()); ++i) {
                if (pairing[i] == -1) {
                    result.push_back(child[i]);
                } else {
                    result.push_back(merge_reconstructed_subtours(child[i],
                                                                  previous[pairing[i]],
                                                                  child[i].demand_bound +
                                                                      previous[pairing[i]].demand_bound));
                }
            }
            out.push_back(std::move(result));

            --child_index;
            if (child_index >= 0 && pairing[child_index] >= 0) {
                used_previous[pairing[child_index]] = false;
            }
            continue;
        }

        bool advanced = false;
        while (next_choice[child_index] <= static_cast<int>(previous.size())) {
            const int choice = next_choice[child_index]++;
            if (choice == -1) {
                pairing[child_index] = -1;
                ++child_index;
                if (child_index < static_cast<int>(child.size())) {
                    next_choice[child_index] = -1;
                }
                advanced = true;
                break;
            }
            if (choice >= static_cast<int>(previous.size())) {
                break;
            }
            if (used_previous[choice]) {
                continue;
            }
            if (previous[choice].demand_bound + child[child_index].demand_bound > 1.0 + kTolerance) {
                continue;
            }
            used_previous[choice] = true;
            pairing[child_index] = choice;
            ++child_index;
            if (child_index < static_cast<int>(child.size())) {
                next_choice[child_index] = -1;
            }
            advanced = true;
            break;
        }

        if (!advanced) {
            next_choice[child_index] = -1;
            if (child_index == 0) {
                break;
            }
            --child_index;
            if (pairing[child_index] >= 0) {
                used_previous[pairing[child_index]] = false;
            }
        }
    }
    return out;
}

bool reconstruct_component_root_value(ReconstructionState& state,
                                      int component_id,
                                      int value_index,
                                      std::vector<ReconstructedSubtour>& out) {
    if (state.component_cache[component_id][value_index].has_value()) {
        out = *state.component_cache[component_id][value_index];
        return true;
    }

    const Component& component = state.context.decomposition.components[component_id];
    const auto& table = state.subtree_phase.component_root_tables[component_id];
    const auto& target_value = table.values[value_index];
    const auto& local_table = state.local_phase.tables[component_id];

    if (component.exit == -1) {
        for (const auto& local_value : local_table.values) {
            SubtreeConfiguration induced;
            for (const auto& entry : local_value.configuration.entries) {
                induced.entries.push_back(SubtreeConfigurationEntry{
                    .demand = entry.demand_bound,
                    .multiplicity = 1,
                });
            }
            normalize_subtree_configuration(induced);
            if (!same_subtree_configuration(induced, target_value.configuration) ||
                !same_double(local_value.cost, target_value.cost)) {
                continue;
            }

            std::vector<LocalReconstructedSubtour> local_subtours;
            if (!reconstruct_local_configuration(state.context,
                                                 component_id,
                                                 local_table.parts,
                                                 local_value.configuration,
                                                 local_value.cost,
                                                 local_subtours)) {
                continue;
            }

            std::vector<ReconstructedSubtour> reconstructed;
            reconstructed.reserve(local_subtours.size());
            for (const auto& local_subtour : local_subtours) {
                reconstructed.push_back(local_subtour.subtour);
            }
            state.component_cache[component_id][value_index] = reconstructed;
            out = reconstructed;
            return true;
        }
        return false;
    }

    const auto* critical_exit_table = critical_table_for_vertex(state.subtree_phase, component.exit);
    const int exit_component_id = component_id_by_root_vertex(state.context.decomposition, component.exit);
    const auto* component_exit_table =
        (exit_component_id != -1 &&
         exit_component_id < static_cast<int>(state.subtree_phase.component_root_tables.size()) &&
         state.subtree_phase.component_root_tables[exit_component_id].vertex == component.exit)
            ? &state.subtree_phase.component_root_tables[exit_component_id]
            : nullptr;

    const double spine_cost = component_spine_cost(state.context, component);
    auto try_exit_table = [&](const SubtreeConfigurationTable& exit_table) {
        for (int exit_index = 0; exit_index < static_cast<int>(exit_table.values.size()); ++exit_index) {
            const auto& exit_value = exit_table.values[exit_index];
            for (const auto& local_value : local_table.values) {
                const std::vector<SubtreeConfiguration> combinations =
                    combine_exit_and_local_configurations(exit_value.configuration, local_value.configuration);
                const int passing_count = static_cast<int>(std::count_if(
                    local_value.configuration.entries.begin(),
                    local_value.configuration.entries.end(),
                    [](const LocalConfigurationEntry& entry) {
                        return entry.type == LocalSubtourType::Passing;
                    }));
                const int exit_subtour_count = subtree_configuration_size(exit_value.configuration);
                const double combined_cost =
                    local_value.cost + exit_value.cost + spine_cost * (exit_subtour_count - passing_count);

                for (const auto& combination : combinations) {
                    if (!same_subtree_configuration(combination, target_value.configuration) ||
                        !same_double(combined_cost, target_value.cost)) {
                        continue;
                    }

                    std::vector<ReconstructedSubtour> exit_subtours;
                    if (!reconstruct_subtree_value_at_vertex(state,
                                                             component.exit,
                                                             exit_value.configuration,
                                                             exit_value.cost,
                                                             exit_subtours)) {
                        continue;
                    }

                    std::vector<LocalReconstructedSubtour> local_subtours;
                    if (!reconstruct_local_configuration(state.context,
                                                         component_id,
                                                         local_table.parts,
                                                         local_value.configuration,
                                                         local_value.cost,
                                                         local_subtours)) {
                        continue;
                    }

                    std::vector<ReconstructedSubtour> reconstructed;
                    if (!combine_exit_and_local_witnesses(exit_subtours,
                                                          local_subtours,
                                                          spine_cost,
                                                          target_value.configuration,
                                                          target_value.cost,
                                                          reconstructed)) {
                        continue;
                    }

                    state.component_cache[component_id][value_index] = reconstructed;
                    out = reconstructed;
                    return true;
                }
            }
        }
        return false;
    };

    if (critical_exit_table != nullptr && try_exit_table(*critical_exit_table)) {
        return true;
    }
    if (component_exit_table != nullptr && try_exit_table(*component_exit_table)) {
        return true;
    }
    return false;
}

bool reconstruct_critical_vertex_value(ReconstructionState& state,
                                       int critical_vertex,
                                       int value_index,
                                       std::vector<ReconstructedSubtour>& out) {
    if (state.critical_cache[critical_vertex][value_index].has_value()) {
        out = *state.critical_cache[critical_vertex][value_index];
        return true;
    }

    const auto& table = state.subtree_phase.critical_vertex_tables[critical_vertex];
    const auto& target_value = table.values[value_index];
    const std::vector<int> child_component_ids =
        child_component_ids_of_critical_vertex(state.context, critical_vertex);

    if (child_component_ids.empty()) {
        if (target_value.configuration.entries.empty() && same_double(target_value.cost, 0.0)) {
            state.critical_cache[critical_vertex][value_index] = std::vector<ReconstructedSubtour>{};
            out.clear();
            return true;
        }
        return false;
    }

    std::vector<SubtreeConfigurationTable> child_tables;
    child_tables.reserve(child_component_ids.size());
    for (const int child_component_id : child_component_ids) {
        child_tables.push_back(state.subtree_phase.component_root_tables[child_component_id]);
    }

    const double beta_for_max_size = beta_value(state.params.epsilon);
    const int max_x_size = std::max(1, static_cast<int>(std::ceil(1.0 / beta_for_max_size)));
    const std::vector<double> candidate_y = candidate_y_values_from_child_tables(child_tables);
    const std::vector<std::vector<double>> x_sets = enumerate_x_sets(candidate_y, max_x_size);

    for (const std::vector<double>& x_values : x_sets) {
        struct CriticalFrame {
            int child_index = 0;
            std::vector<ReconstructedSubtour> current;
            std::vector<std::vector<ReconstructedSubtour>> next_candidates;
            std::size_t next_candidate_index = 0;
            bool expanded = false;
        };

        std::vector<CriticalFrame> stack;
        stack.push_back(CriticalFrame{.child_index = 0, .current = {}});
        while (!stack.empty()) {
            CriticalFrame& frame = stack.back();
            if (frame.child_index == static_cast<int>(child_component_ids.size())) {
                if (same_subtree_configuration(subtree_configuration_from_reconstructed_tours(frame.current),
                                               target_value.configuration) &&
                    same_double(total_reconstructed_cost(frame.current), target_value.cost)) {
                    state.critical_cache[critical_vertex][value_index] = frame.current;
                    out = frame.current;
                    return true;
                }
                stack.pop_back();
                continue;
            }

            if (!frame.expanded) {
                const int child_component_id = child_component_ids[frame.child_index];
                const int child_root = state.context.decomposition.components[child_component_id].root;
                const double edge_cost =
                    edge_cost_in_height_reduced_tree(state.context, critical_vertex, child_root);
                const auto& child_table = state.subtree_phase.component_root_tables[child_component_id];

                frame.next_candidates.clear();
                for (int child_value_index = 0;
                     child_value_index < static_cast<int>(child_table.values.size());
                     ++child_value_index) {
                    std::vector<ReconstructedSubtour> child_subtours;
                    if (!reconstruct_component_root_value(state,
                                                          child_component_id,
                                                          child_value_index,
                                                          child_subtours)) {
                        continue;
                    }

                    for (auto& subtour : child_subtours) {
                        subtour.cost += 2.0 * edge_cost;
                        auto it = std::lower_bound(x_values.begin(),
                                                   x_values.end(),
                                                   subtour.demand_bound - kTolerance);
                        if (it == x_values.end()) {
                            throw std::logic_error("failed to round reconstructed child subtour");
                        }
                        subtour.demand_bound = *it;
                    }

                    const auto candidates = combine_critical_vertex_witnesses(frame.current, child_subtours);
                    frame.next_candidates.insert(frame.next_candidates.end(), candidates.begin(), candidates.end());
                }
                frame.next_candidate_index = 0;
                frame.expanded = true;
            }

            if (frame.next_candidate_index >= frame.next_candidates.size()) {
                stack.pop_back();
                continue;
            }

            auto next_current = frame.next_candidates[frame.next_candidate_index++];
            stack.push_back(CriticalFrame{
                .child_index = frame.child_index + 1,
                .current = std::move(next_current),
            });
        }
    }

    return false;
}

}  // namespace

SolveResult OnePointFiveApproxSolver::solve(const Instance& instance) {
    return solve(instance, OnePointFiveApproxParams{});
}

// Apply the outer bounded-distance reduction before dispatching to the bounded-distance solver.
SolveResult OnePointFiveApproxSolver::solve(const Instance& instance, const OnePointFiveApproxParams& params) {
    SolveResult result;
    const OnePointFiveApproxParams normalized_params = validate_params(params);
    instance.validate();
    const int one_over_epsilon = static_cast<int>(std::round(1.0 / normalized_params.epsilon));
    const Instance normalized_instance = Preprocessor::make_binary_leaf_tree(instance);

    // Partition terminals by randomized logarithmic distance buckets.
    std::vector<Terminal> terminals = normalized_instance.terminals();
    std::vector<double> dists = normalized_instance.distances_from_depot();
    if (terminals.empty()) {
        return result;
    }
    std::sort(terminals.begin(), terminals.end(),
        [&dists](const Terminal& a, const Terminal& b) {
            return dists[a.vertex] < dists[b.vertex];
        });
    const int i_0 = Rng::uniform_int(0, one_over_epsilon - 1);

    int cur_exp_idx = INT_MIN;
    std::size_t i = 0;
    while (i < terminals.size()) {
        // Y_j captures one bounded-distance slab selected by the random offset.
        std::vector<Terminal> y_j;
        int current_bucket = static_cast<int>(std::floor(std::log(dists[terminals[i].vertex]) / std::log(one_over_epsilon)));
        cur_exp_idx = std::max(cur_exp_idx, floor_div(current_bucket - i_0, one_over_epsilon));
        int cur_exp = one_over_epsilon * cur_exp_idx + i_0;
        while (i < terminals.size() &&
               dists[terminals[i].vertex] < std::pow(one_over_epsilon, cur_exp + 1)) {
            y_j.push_back(terminals[i]);
            ++i;
        }
        if (!y_j.empty()) {
            const Instance subinstance = Instance::with_terminals(normalized_instance, y_j);
            const SolveResult normalized_subresult = solve_bounded_distance(subinstance, normalized_params);
            const SolveResult subresult =
                project_normalized_solution_to_original(instance, normalized_instance, normalized_subresult);
            result.cost += subresult.cost;
            result.tours.insert(result.tours.end(), subresult.tours.begin(), subresult.tours.end());
        }

        // Z_j aggregates the intervening slabs that will be recursed on separately.
        std::vector<Terminal> z_j;
        for (int j = 1; j < one_over_epsilon; ++j) {
            cur_exp = one_over_epsilon * cur_exp_idx + i_0 + j;
            while (i < terminals.size() &&
                   dists[terminals[i].vertex] < std::pow(one_over_epsilon, cur_exp + 1)) {
                z_j.push_back(terminals[i]);
                ++i;
            }
        }
        if (!z_j.empty()) {
            const Instance subinstance = Instance::with_terminals(normalized_instance, z_j);
            const SolveResult normalized_subresult = solve_bounded_distance(subinstance, normalized_params);
            const SolveResult subresult =
                project_normalized_solution_to_original(instance, normalized_instance, normalized_subresult);
            result.cost += subresult.cost;
            result.tours.insert(result.tours.end(), subresult.tours.begin(), subresult.tours.end());
        }

        ++cur_exp_idx;
    }
    
    return result;
}

BoundedDistanceContext OnePointFiveApproxSolver::build_bounded_distance_context(
    const Instance& instance,
    const OnePointFiveApproxParams& params) {
    BoundedDistanceContext context;
    context.instance = instance;
    context.rooted_tree = RootedTreeBuilder::build(context.instance);
    context.decomposition =
        DecompositionBuilder::decompose_bounded_instance(context.rooted_tree, params.epsilon);
    context.height_reduced =
        DecompositionBuilder::height_reduce_bounded_components(context.decomposition,
                                                               context.rooted_tree,
                                                               params.epsilon);
    return context;
}

LocalConfigurationTable OnePointFiveApproxSolver::compute_local_configurations(
    const BoundedDistanceContext& context,
    int component_id,
    const OnePointFiveApproxParams& params) {
    if (component_id < 0 || component_id >= static_cast<int>(context.decomposition.components.size())) {
        throw std::out_of_range("component_id is outside the component decomposition");
    }

    const Component& component = context.decomposition.components[component_id];
    const double alpha = alpha_value(params.epsilon);
    const double gamma_zero = gamma_zero_value(params.epsilon);

    LocalConfigurationTable table;
    table.component_id = component_id;
    table.alpha = alpha;
    table.parts = component_parts(context, component_id, gamma_zero);
    table.y_values = y_values_for_component(table.parts, alpha);

    const bool allow_passing = (component.exit != -1);
    const std::vector<LocalConfiguration> configurations =
        enumerate_local_configurations(table.y_values, allow_passing, static_cast<int>(table.parts.size()));

    // Algorithm 2: for every local configuration A, try every feasible partition of Q_c into |A| groups
    // that respect the demand bounds y_i, and keep the minimum subtour cost.
    for (const LocalConfiguration& configuration : configurations) {
        table.values.push_back(LocalConfigurationValue{
            .configuration = configuration,
            .cost = configuration_cost(context, component_id, table.parts, configuration),
        });
    }

    return table;
}

LocalPhaseState OnePointFiveApproxSolver::compute_local_phase(const BoundedDistanceContext& context,
                                                              const OnePointFiveApproxParams& params) {
    LocalPhaseState state;
    state.tables.reserve(context.decomposition.components.size());
    for (int component_id = 0; component_id < static_cast<int>(context.decomposition.components.size()); ++component_id) {
        state.tables.push_back(compute_local_configurations(context, component_id, params));
    }
    return state;
}

SubtreeConfigurationTable OnePointFiveApproxSolver::compute_component_root_subtree_configurations(
    const BoundedDistanceContext& context,
    int component_id,
    const LocalConfigurationTable& local_table,
    const SubtreeConfigurationTable* exit_table,
    const OnePointFiveApproxParams&) {
    if (component_id < 0 || component_id >= static_cast<int>(context.decomposition.components.size())) {
        throw std::out_of_range("component_id is outside the component decomposition");
    }
    if (local_table.component_id != component_id) {
        throw std::invalid_argument("local_table does not correspond to the requested component");
    }

    const Component& component = context.decomposition.components[component_id];
    SubtreeConfigurationTable table;
    table.vertex = component.root;

    if (component.exit == -1) {
        // Leaf-component case from Appendix C.1: a local configuration induces a subtree configuration
        // directly by turning every local subtour into one subtree subtour rooted at r_c.
        for (const auto& local_value : local_table.values) {
            SubtreeConfiguration configuration;
            for (const auto& entry : local_value.configuration.entries) {
                configuration.entries.push_back(SubtreeConfigurationEntry{
                    .demand = entry.demand_bound,
                    .multiplicity = 1,
                });
            }
            record_subtree_value(table, std::move(configuration), local_value.cost);
        }
        return table;
    }

    if (exit_table == nullptr || exit_table->vertex != component.exit) {
        throw std::invalid_argument(
            "internal component root DP requires subtree configurations at the exit vertex");
    }

    const double spine_cost = component_spine_cost(context, component);
    for (const auto& exit_value : exit_table->values) {
        for (const auto& local_value : local_table.values) {
            const std::vector<SubtreeConfiguration> combinations =
                combine_exit_and_local_configurations(exit_value.configuration, local_value.configuration);
            const int passing_count = static_cast<int>(std::count_if(
                local_value.configuration.entries.begin(),
                local_value.configuration.entries.end(),
                [](const LocalConfigurationEntry& entry) {
                    return entry.type == LocalSubtourType::Passing;
                }));
            const int exit_subtour_count = subtree_configuration_size(exit_value.configuration);
            const double combined_cost =
                local_value.cost + exit_value.cost + spine_cost * (exit_subtour_count - passing_count);

            // Algorithm 5: each way to associate the passing local subtours with exit subtours yields
            // one resulting subtree configuration at r_c.
            for (const auto& configuration : combinations) {
                record_subtree_value(table, configuration, combined_cost);
            }
        }
    }

    return table;
}

SubtreeConfigurationTable OnePointFiveApproxSolver::compute_critical_vertex_subtree_configurations(
    const BoundedDistanceContext& context,
    int critical_vertex,
    const std::vector<SubtreeConfigurationTable>& child_tables,
    const OnePointFiveApproxParams& params) {
    const std::vector<int> child_component_ids = child_component_ids_of_critical_vertex(context, critical_vertex);
    if (child_component_ids.size() != child_tables.size()) {
        throw std::invalid_argument(
            "child_tables must correspond exactly to the child component roots of the critical vertex");
    }
    for (std::size_t i = 0; i < child_component_ids.size(); ++i) {
        const int expected_vertex = context.decomposition.components[child_component_ids[i]].root;
        if (child_tables[i].vertex != expected_vertex) {
            throw std::invalid_argument("child_tables are not ordered by the critical vertex's child roots");
        }
    }

    SubtreeConfigurationTable table;
    table.vertex = critical_vertex;
    if (child_tables.empty()) {
        record_subtree_value(table, SubtreeConfiguration{}, 0.0);
        return table;
    }

    const double beta = beta_value(params.epsilon);
    const int max_x_size = std::max(1, static_cast<int>(std::ceil(1.0 / beta)));
    const std::vector<double> candidate_y = candidate_y_values_from_child_tables(child_tables);
    const std::vector<std::vector<double>> x_sets = enumerate_x_sets(candidate_y, max_x_size);

    // Algorithm 6 processes the child component roots of z from left to right. For each guessed set X,
    // it rounds every child table to X and then performs a DP over prefixes of the children, where the
    // state is the current rounded sum list at z.
    for (const std::vector<double>& x_values : x_sets) {
        std::vector<SubtreeConfigurationValue> current_prefix = {
            SubtreeConfigurationValue{.configuration = SubtreeConfiguration{}, .cost = 0.0},
        };

        for (std::size_t child_index = 0; child_index < child_tables.size(); ++child_index) {
            const int component_id = child_component_ids[child_index];
            const int child_root = context.decomposition.components[component_id].root;
            const double edge_cost = edge_cost_in_height_reduced_tree(context, critical_vertex, child_root);

            std::vector<SubtreeConfigurationValue> next_prefix;
            for (const auto& prefix_value : current_prefix) {
                for (const auto& child_value : child_tables[child_index].values) {
                    const SubtreeConfiguration rounded_child =
                        round_subtree_configuration_to_x(child_value.configuration, x_values);
                    const std::vector<SubtreeConfiguration> combined =
                        combine_sum_lists(prefix_value.configuration, rounded_child);
                    const double total_cost =
                        prefix_value.cost +
                        child_value.cost +
                        2.0 * subtree_configuration_size(child_value.configuration) * edge_cost;

                    for (const auto& configuration : combined) {
                        bool updated = false;
                        for (auto& next_value : next_prefix) {
                            if (same_subtree_configuration(next_value.configuration, configuration)) {
                                next_value.cost = std::min(next_value.cost, total_cost);
                                updated = true;
                                break;
                            }
                        }
                        if (!updated) {
                            next_prefix.push_back(SubtreeConfigurationValue{
                                .configuration = configuration,
                                .cost = total_cost,
                            });
                        }
                    }
                }
            }
            current_prefix = std::move(next_prefix);
        }

        for (const auto& value : current_prefix) {
            record_subtree_value(table, value.configuration, value.cost);
        }
    }

    return table;
}

SubtreePhaseState OnePointFiveApproxSolver::compute_subtree_phase(const BoundedDistanceContext& context,
                                                                  const LocalPhaseState& local_phase,
                                                                  const OnePointFiveApproxParams& params) {
    if (local_phase.tables.size() != context.decomposition.components.size()) {
        throw std::invalid_argument("local_phase must contain one local table per component");
    }
    if (context.height_reduced.critical_vertex_by_component.size() != context.decomposition.components.size()) {
        throw std::invalid_argument(
            "compute_subtree_phase requires critical_vertex_by_component for every component");
    }
    if (context.height_reduced.original_parent_component.size() != context.decomposition.components.size()) {
        throw std::invalid_argument(
            "compute_subtree_phase requires original_parent_component for every component");
    }

    SubtreePhaseState state;
    state.component_root_tables.resize(context.decomposition.components.size());
    state.critical_vertex_tables.resize(context.rooted_tree.parent.size());
    for (auto& table : state.component_root_tables) {
        table.vertex = -1;
    }
    for (auto& table : state.critical_vertex_tables) {
        table.vertex = -1;
    }

    std::vector<bool> component_done(context.decomposition.components.size(), false);
    std::vector<bool> critical_done(context.rooted_tree.parent.size(), false);

    // Seed all leaf-component roots immediately: Algorithm 5 for a leaf component depends only on f(c, ·).
    for (const Component& component : context.decomposition.components) {
        if (component.exit == -1) {
            state.component_root_tables[component.id] =
                compute_component_root_subtree_configurations(context,
                                                              component.id,
                                                              local_phase.tables[component.id],
                                                              nullptr,
                                                              params);
            component_done[component.id] = true;
        }
    }

    int completed_components = 0;
    for (bool done : component_done) {
        if (done) {
            ++completed_components;
        }
    }
    int completed_critical_vertices = 0;

    while (completed_components < static_cast<int>(context.decomposition.components.size()) ||
           completed_critical_vertices < static_cast<int>(context.height_reduced.groups.size())) {
        bool made_progress = false;

        // First, compute every component-root table whose exit table is already available.
        for (const Component& component : context.decomposition.components) {
            if (component_done[component.id]) {
                continue;
            }

            const SubtreeConfigurationTable* exit_table = nullptr;
            if (component.exit != -1) {
                if (component.exit < static_cast<int>(critical_done.size()) && critical_done[component.exit]) {
                    exit_table = &state.critical_vertex_tables[component.exit];
                } else {
                    for (int i = 0; i < static_cast<int>(state.component_root_tables.size()); ++i) {
                        if (component_done[i] && state.component_root_tables[i].vertex == component.exit) {
                            exit_table = &state.component_root_tables[i];
                            break;
                        }
                    }
                }
                if (exit_table == nullptr) {
                    continue;
                }
            }

            state.component_root_tables[component.id] = compute_component_root_subtree_configurations(
                context,
                component.id,
                local_phase.tables[component.id],
                exit_table,
                params);
            component_done[component.id] = true;
            ++completed_components;
            made_progress = true;
        }

        // Then, compute every critical-vertex table whose child component-root tables are ready.
        for (const auto& group : context.height_reduced.groups) {
            const int critical_vertex = group.critical_vertex;
            if (critical_done[critical_vertex]) {
                continue;
            }

            const std::vector<int> child_component_ids =
                child_component_ids_of_critical_vertex(context, critical_vertex);
            bool children_ready = true;
            std::vector<SubtreeConfigurationTable> child_tables;
            child_tables.reserve(child_component_ids.size());
            for (const int child_component_id : child_component_ids) {
                if (!component_done[child_component_id]) {
                    children_ready = false;
                    break;
                }
                child_tables.push_back(state.component_root_tables[child_component_id]);
            }
            if (!children_ready) {
                continue;
            }

            state.critical_vertex_tables[critical_vertex] =
                compute_critical_vertex_subtree_configurations(context, critical_vertex, child_tables, params);
            critical_done[critical_vertex] = true;
            ++completed_critical_vertices;
            made_progress = true;
        }

        if (!made_progress) {
            throw std::logic_error("subtree phase could not make progress; height-reduced dependencies are inconsistent");
        }
    }

    return state;
}

double OnePointFiveApproxSolver::bounded_height_reduced_opt_value(const SubtreePhaseState& subtree_phase,
                                                                  int root_vertex) {
    if (root_vertex < 0 || root_vertex >= static_cast<int>(subtree_phase.critical_vertex_tables.size()) ||
        subtree_phase.critical_vertex_tables[root_vertex].vertex != root_vertex) {
        throw std::invalid_argument("root_vertex does not have a computed critical-vertex table");
    }

    double best = std::numeric_limits<double>::infinity();
    for (const auto& value : subtree_phase.critical_vertex_tables[root_vertex].values) {
        best = std::min(best, value.cost);
    }
    return best;
}

SolveResult OnePointFiveApproxSolver::solve_bounded_distance(const Instance& instance,
                                                             const OnePointFiveApproxParams& params) {
    const OnePointFiveApproxParams normalized_params = validate_params(params);
    if (instance.terminals().empty()) {
        return {};
    }
    const BoundedDistanceContext context = build_bounded_distance_context(instance, normalized_params);
    const LocalPhaseState local_phase = compute_local_phase(context, normalized_params);
    const SubtreePhaseState subtree_phase = compute_subtree_phase(context, local_phase, normalized_params);
    const int depot = context.instance.depot();
    const int root_component_id = component_id_by_root_vertex(context.decomposition, depot);
    const bool use_root_component =
        root_component_id != -1 &&
        (context.decomposition.components[root_component_id].terminal_count > 0 ||
         !context.decomposition.components[root_component_id].block_ids.empty());

    int best_index = -1;
    double best_cost = std::numeric_limits<double>::infinity();
    bool reconstruct_component = false;
    if (use_root_component) {
        const auto& root_table = subtree_phase.component_root_tables[root_component_id];
        for (int i = 0; i < static_cast<int>(root_table.values.size()); ++i) {
            if (root_table.values[i].cost < best_cost) {
                best_cost = root_table.values[i].cost;
                best_index = i;
            }
        }
        reconstruct_component = true;
    } else {
        const auto* root_table = critical_table_for_vertex(subtree_phase, depot);
        if (root_table == nullptr || root_table->values.empty()) {
            throw std::logic_error("bounded-distance solver did not produce a usable root table");
        }
        for (int i = 0; i < static_cast<int>(root_table->values.size()); ++i) {
            if (root_table->values[i].cost < best_cost) {
                best_cost = root_table->values[i].cost;
                best_index = i;
            }
        }
    }
    if (best_index == -1) {
        throw std::logic_error("failed to choose a bounded height-reduced root value");
    }

    ReconstructionState reconstruction{
        .context = context,
        .local_phase = local_phase,
        .subtree_phase = subtree_phase,
        .params = normalized_params,
    };
    reconstruction.component_cache.resize(subtree_phase.component_root_tables.size());
    for (std::size_t i = 0; i < subtree_phase.component_root_tables.size(); ++i) {
        reconstruction.component_cache[i].resize(subtree_phase.component_root_tables[i].values.size());
    }
    reconstruction.critical_cache.resize(subtree_phase.critical_vertex_tables.size());
    for (std::size_t i = 0; i < subtree_phase.critical_vertex_tables.size(); ++i) {
        reconstruction.critical_cache[i].resize(subtree_phase.critical_vertex_tables[i].values.size());
    }

    std::vector<ReconstructedSubtour> reduced_subtours;
    const bool reconstructed =
        reconstruct_component
            ? reconstruct_component_root_value(reconstruction, root_component_id, best_index, reduced_subtours)
            : reconstruct_critical_vertex_value(reconstruction, depot, best_index, reduced_subtours);
    if (!reconstructed) {
        throw std::logic_error("failed to reconstruct optimum bounded height-reduced witness");
    }

    const SolveResult reduced_solution = to_public_solution(reduced_subtours);
    if (!same_double(reduced_solution.cost, best_cost)) {
        throw std::logic_error("reconstructed bounded height-reduced witness has the wrong total cost");
    }
    return DecompositionBuilder::lift_solution_from_height_reduced_tree(reduced_solution, context.instance);
}

}  // namespace tucvrp
