#include "internal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

namespace tucvrp {

namespace {

// In the component tree, the parent of a non-root component is the unique component that contains
// both the component root r_c and its parent vertex in the original rooted tree. That parent
// component is exactly the one owning the edge that leaves r_c toward the depot.
std::vector<int> component_parents(const TreeDecomposition& decomposition, const RootedTreeData& rooted_tree) {
    std::vector<bool> seen_root(rooted_tree.parent.size(), false);
    for (const Component& component : decomposition.components) {
        if (component.root < 0 || component.root >= static_cast<int>(seen_root.size())) {
            throw std::logic_error("component root is outside the rooted tree");
        }
        if (seen_root[component.root]) {
            throw std::logic_error("height reduction requires component roots to be unique");
        }
        seen_root[component.root] = true;
    }

    std::vector<std::vector<bool>> in_component;
    in_component.reserve(decomposition.components.size());
    for (const Component& component : decomposition.components) {
        std::vector<bool> membership(rooted_tree.parent.size(), false);
        for (const int v : component.vertices) {
            membership[v] = true;
        }
        in_component.push_back(std::move(membership));
    }

    std::vector<int> parent_component(decomposition.components.size(), -1);
    for (const Component& component : decomposition.components) {
        const int parent_vertex = rooted_tree.parent[component.root];
        if (parent_vertex == -1) {
            continue;
        }

        int parent_id = -1;
        for (const Component& candidate : decomposition.components) {
            if (candidate.id == component.id) {
                continue;
            }
            // The upward edge leaving r_c must stay entirely inside the parent component.
            if (in_component[candidate.id][component.root] && in_component[candidate.id][parent_vertex]) {
                if (parent_id != -1) {
                    throw std::logic_error("component tree has multiple parents for one component");
                }
                parent_id = candidate.id;
            }
        }

        if (parent_id == -1) {
            throw std::logic_error("failed to find parent component for non-root component");
        }
        parent_component[component.id] = parent_id;
    }

    return parent_component;
}

double min_terminal_distance(const RootedTreeData& rooted_tree) {
    double d_min = std::numeric_limits<double>::infinity();
    for (const int terminal : rooted_tree.terminal_vertices) {
        d_min = std::min(d_min, rooted_tree.distances_from_depot[terminal]);
    }
    return d_min;
}

int class_index_for_distance(double distance, double d_tilde, int h_epsilon) {
    if (d_tilde <= 0.0 || h_epsilon <= 0) {
        throw std::logic_error("invalid height-reduction parameters");
    }
    // Lemma 21: C_i contains components whose root distances lie in
    // [ (i - 1) * D~, i * D~ ).
    const int index = static_cast<int>(std::floor(distance / d_tilde)) + 1;
    return std::clamp(index, 1, h_epsilon);
}

}  // namespace

HeightReducedComponentTree DecompositionBuilder::height_reduce_bounded_components(
    const TreeDecomposition& decomposition,
    const RootedTreeData& rooted_tree,
    double epsilon) {
    if (epsilon <= 0.0 || epsilon >= 1.0) {
        throw std::invalid_argument("height_reduce_bounded_components requires epsilon in (0, 1)");
    }

    HeightReducedComponentTree reduced;
    reduced.original_parent_component = component_parents(decomposition, rooted_tree);
    reduced.class_index_by_component.assign(decomposition.components.size(), -1);
    reduced.group_id_by_component.assign(decomposition.components.size(), -1);
    reduced.critical_vertex_by_component.assign(decomposition.components.size(), -1);
    reduced.attachment_length_by_component.assign(decomposition.components.size(), 0.0);

    if (decomposition.components.empty() || rooted_tree.terminal_vertices.empty()) {
        return reduced;
    }

    // Section 4 constants: alpha = epsilon^(1/epsilon + 1), D~ = alpha * epsilon * Dmin,
    // and H_epsilon = (1 / epsilon)^(2 / epsilon + 1).
    reduced.d_min = min_terminal_distance(rooted_tree);
    if (!std::isfinite(reduced.d_min) || reduced.d_min <= 0.0) {
        throw std::logic_error("height reduction expects positive terminal distances");
    }

    const double one_over_epsilon = 1.0 / epsilon;
    const double alpha = std::pow(epsilon, one_over_epsilon + 1.0);
    reduced.d_tilde = alpha * epsilon * reduced.d_min;
    reduced.h_epsilon = static_cast<int>(std::ceil(std::pow(one_over_epsilon, 2.0 * one_over_epsilon + 1.0)));

    for (const Component& component : decomposition.components) {
        reduced.class_index_by_component[component.id] = class_index_for_distance(
            rooted_tree.distances_from_depot[component.root], reduced.d_tilde, reduced.h_epsilon);
    }

    std::vector<std::vector<int>> component_children(decomposition.components.size());
    for (const Component& component : decomposition.components) {
        const int parent_id = reduced.original_parent_component[component.id];
        if (parent_id != -1) {
            component_children[parent_id].push_back(component.id);
        }
    }

    std::vector<bool> visited(decomposition.components.size(), false);
    for (const Component& start : decomposition.components) {
        if (visited[start.id]) {
            continue;
        }

        const int class_index = reduced.class_index_by_component[start.id];
        HeightReducedComponentGroup group;
        group.id = static_cast<int>(reduced.groups.size());
        group.class_index = class_index;

        std::queue<int> queue;
        queue.push(start.id);
        visited[start.id] = true;
        while (!queue.empty()) {
            const int component_id = queue.front();
            queue.pop();
            group.component_ids.push_back(component_id);

            const int parent_id = reduced.original_parent_component[component_id];
            if (parent_id != -1 && !visited[parent_id] &&
                reduced.class_index_by_component[parent_id] == class_index) {
                visited[parent_id] = true;
                queue.push(parent_id);
            }

            for (const int child_id : component_children[component_id]) {
                if (!visited[child_id] && reduced.class_index_by_component[child_id] == class_index) {
                    visited[child_id] = true;
                    queue.push(child_id);
                }
            }
        }

        // Definition 22: the critical vertex of a maximally connected set is the root vertex of the
        // component closest to the depot.
        int critical_component = group.component_ids.front();
        double best_distance = rooted_tree.distances_from_depot[decomposition.components[critical_component].root];
        for (const int component_id : group.component_ids) {
            const double distance = rooted_tree.distances_from_depot[decomposition.components[component_id].root];
            if (distance < best_distance - 1e-9 ||
                (std::abs(distance - best_distance) <= 1e-9 && component_id < critical_component)) {
                critical_component = component_id;
                best_distance = distance;
            }
        }
        group.critical_vertex = decomposition.components[critical_component].root;

        for (const int component_id : group.component_ids) {
            const int root = decomposition.components[component_id].root;
            reduced.group_id_by_component[component_id] = group.id;
            reduced.critical_vertex_by_component[component_id] = group.critical_vertex;
            // Algorithm 1 reattaches each component root r_c directly to the critical vertex z of its
            // maximally connected set, with edge weight equal to the original r_c-to-z distance.
            reduced.attachment_length_by_component[component_id] =
                rooted_tree.distances_from_depot[root] - rooted_tree.distances_from_depot[group.critical_vertex];
            if (reduced.attachment_length_by_component[component_id] < -1e-9) {
                throw std::logic_error("critical vertex is below a component root");
            }
            reduced.attachment_length_by_component[component_id] =
                std::max(0.0, reduced.attachment_length_by_component[component_id]);
        }

        std::sort(group.component_ids.begin(), group.component_ids.end());
        reduced.groups.push_back(std::move(group));
    }

    std::sort(reduced.groups.begin(),
              reduced.groups.end(),
              [](const HeightReducedComponentGroup& a, const HeightReducedComponentGroup& b) {
                  if (a.class_index != b.class_index) {
                      return a.class_index < b.class_index;
                  }
                  return a.critical_vertex < b.critical_vertex;
              });
    for (int i = 0; i < static_cast<int>(reduced.groups.size()); ++i) {
        reduced.groups[i].id = i;
        for (const int component_id : reduced.groups[i].component_ids) {
            reduced.group_id_by_component[component_id] = i;
        }
    }

    return reduced;
}

SolveResult DecompositionBuilder::lift_solution_from_height_reduced_tree(
    const SolveResult& reduced_solution,
    const Instance& bounded_instance) {
    bounded_instance.validate();

    SolveResult lifted;
    for (const Tour& reduced_tour : reduced_solution.tours) {
        Tour lifted_tour;
        lifted_tour.terminals = reduced_tour.terminals;

        for (const int terminal : lifted_tour.terminals) {
            if (!bounded_instance.is_terminal(terminal)) {
                throw std::invalid_argument(
                    "lift_solution_from_height_reduced_tree received a non-terminal vertex");
            }
            lifted_tour.demand += bounded_instance.demand_of(terminal);
        }

        if (lifted_tour.demand > 1.0 + 1e-9) {
            throw std::invalid_argument(
                "lift_solution_from_height_reduced_tree received an over-capacity tour");
        }

        const double original_min_cost = bounded_instance.tour_cost_for_terminals(lifted_tour.terminals);
        // Under the paper's construction, every height-reduced edge expands to an equal-length path
        // in the original bounded tree. Since the current Tour abstraction stores only the served
        // terminals and the route cost, we preserve the reduced-tree route cost and sanity-check that
        // a tour of no greater cost exists in the original tree.
        if (original_min_cost > reduced_tour.cost + 1e-9) {
            throw std::logic_error(
                "height-reduced tour cannot be realized in the original bounded tree at the same cost");
        }

        lifted_tour.walk = bounded_instance.tour_walk_for_terminals(lifted_tour.terminals);
        lifted_tour.cost = reduced_tour.cost;
        lifted.tours.push_back(std::move(lifted_tour));
        lifted.cost += reduced_tour.cost;
    }

    return lifted;
}

}  // namespace tucvrp
