#include "tucvrp/decomposition.hpp"

namespace tucvrp {

void DecompositionBuilder::decompose_clusters_into_cells(TreeDecomposition& decomposition,
                                                         const RootedTreeData& rooted_tree,
                                                         double epsilon) {
    (void)decomposition;
    (void)rooted_tree;
    (void)epsilon;
    // TODO(section 4.3): split each cluster into cells.
}

}  // namespace tucvrp
