#pragma once

#include "tucvrp/decomposition.hpp"
#include "tucvrp/instance.hpp"
#include "tucvrp/rooted_tree.hpp"
#include "tucvrp/solution.hpp"

#include <string>

namespace tucvrp {

// Tunable parameters for the paper's (1.5 + epsilon)-approximation.
struct OnePointFiveApproxParams {
    double epsilon = 0.1;
};

// Global bounded-distance state used by the decomposition, height reduction, and DP phases.
struct BoundedDistanceContext {
    Instance instance;
    RootedTreeData rooted_tree;
    TreeDecomposition decomposition;
    HeightReducedComponentTree height_reduced;
};

enum class LocalSubtourType {
    Passing,
    Ending,
};

// One part of Q_c from Definition 26: either one big terminal, or all small terminals in one cell.
struct LocalPart {
    std::vector<int> terminals;
    double demand = 0.0;
};

// One pair (y_i, b_i) in a local configuration from Definition 34.
struct LocalConfigurationEntry {
    double demand_bound = 0.0;
    LocalSubtourType type = LocalSubtourType::Ending;
};

// A local configuration A for one component.
struct LocalConfiguration {
    std::vector<LocalConfigurationEntry> entries;
};

// One computed value f(c, A).
struct LocalConfigurationValue {
    LocalConfiguration configuration;
    double cost = 0.0;
};

// The full local table for one component produced by Algorithm 2.
struct LocalConfigurationTable {
    int component_id = -1;
    double alpha = 0.0;
    std::vector<LocalPart> parts;
    std::vector<double> y_values;
    std::vector<LocalConfigurationValue> values;
};

// Output of the first DP phase from Section 9.1.
struct LocalPhaseState {
    std::vector<LocalConfigurationTable> tables;
};

// One pair (y_i, n_i) in a subtree configuration from Definition 37.
struct SubtreeConfigurationEntry {
    double demand = 0.0;
    int multiplicity = 0;
};

// A subtree configuration list A at a vertex v.
struct SubtreeConfiguration {
    std::vector<SubtreeConfigurationEntry> entries;
};

// One computed value g(v, A).
struct SubtreeConfigurationValue {
    SubtreeConfiguration configuration;
    double cost = 0.0;
};

// The table of subtree-configuration values at one vertex.
struct SubtreeConfigurationTable {
    int vertex = -1;
    std::vector<SubtreeConfigurationValue> values;
};

// Output of the second DP phase after propagating subtree tables through the height-reduced tree.
struct SubtreePhaseState {
    std::vector<SubtreeConfigurationTable> component_root_tables;
    std::vector<SubtreeConfigurationTable> critical_vertex_tables;
};

// Main entry point for the paper's approximation algorithm.
class OnePointFiveApproxSolver {
  public:
    // Solve using the paper algorithm with default approximation parameters.
    static SolveResult solve(const Instance& instance);
    // Solve using the paper algorithm with explicit approximation parameters.
    static SolveResult solve(const Instance& instance, const OnePointFiveApproxParams& params);
    // Solve the bounded-distance subproblem after the outer reduction step.
    static SolveResult solve_bounded_distance(const Instance& instance, const OnePointFiveApproxParams& params);
    // Compute the Algorithm 2 table of local configuration values for one component.
    static LocalConfigurationTable compute_local_configurations(const BoundedDistanceContext& context,
                                                               int component_id,
                                                               const OnePointFiveApproxParams& params);
    // Compute Algorithm 5 at the root of one component. For internal components, `exit_table`
    // must contain the already-computed subtree configurations at the exit vertex e_c.
    static SubtreeConfigurationTable compute_component_root_subtree_configurations(
        const BoundedDistanceContext& context,
        int component_id,
        const LocalConfigurationTable& local_table,
        const SubtreeConfigurationTable* exit_table,
        const OnePointFiveApproxParams& params);
    // Compute Algorithm 6 at one critical vertex z from the already-computed tables at its child roots.
    static SubtreeConfigurationTable compute_critical_vertex_subtree_configurations(
        const BoundedDistanceContext& context,
        int critical_vertex,
        const std::vector<SubtreeConfigurationTable>& child_tables,
        const OnePointFiveApproxParams& params);
    // Run the full bottom-up subtree DP over the height-reduced tree.
    static SubtreePhaseState compute_subtree_phase(const BoundedDistanceContext& context,
                                                   const LocalPhaseState& local_phase,
                                                   const OnePointFiveApproxParams& params);
    // Extract the optimum bounded height-reduced value from the final root table.
    static double bounded_height_reduced_opt_value(const SubtreePhaseState& subtree_phase, int root_vertex);

  private:
    // Build the bounded-distance context: decomposition plus the Section 4 height reduction.
    static BoundedDistanceContext build_bounded_distance_context(const Instance& instance,
                                                                const OnePointFiveApproxParams& params);
    // Section 9.1: compute all local-configuration tables, one per component.
    static LocalPhaseState compute_local_phase(const BoundedDistanceContext& context,
                                               const OnePointFiveApproxParams& params);
};

}  // namespace tucvrp
