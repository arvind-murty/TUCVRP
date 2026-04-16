# TUCVRP

C++20 implementation work for tree UCVRP, centered on the paper
["A Tight (1.5 + epsilon)-Approximation for Unsplittable Capacitated Vehicle Routing on Trees"](https://arxiv.org/pdf/2202.05691),
with the bounded-height reduction infrastructure from
["A PTAS for Unsplittable Capacitated Vehicle Routing on Trees"](https://arxiv.org/pdf/2111.03735).

The codebase currently contains:

- a general tree instance model and parser
- an exact exponential solver for small instances
- the paper solver pipeline:
  - preprocessing to a binary tree with leaf terminals
  - bounded-distance outer reduction
  - component / block / cluster / cell decomposition
  - height reduction on the component tree
  - local-configuration DP
  - subtree-configuration DP
  - tour reconstruction back to the original instance

## Build

This project uses `make` and a system-installed Catch2.

```bash
make
make test
```

The `Makefile` discovers Catch2 through `pkg-config`.

## Current Status

Implemented:

- rooted tree instance model with validation
- text parser for problem instances
- exact subset-DP solver for small instances
- edge-load lower bound
- project-wide RNG helper
- preprocessing:
  - binary-tree normalization
  - terminal-to-leaf normalization
  - removal of zero-distance terminals during normalization
- decomposition hierarchy:
  - bounded-instance component decomposition
  - block decomposition
  - cluster decomposition
  - cell decomposition
- height reduction of the component tree
- Section 9 dynamic-programming phases:
  - Algorithm 2 local configurations
  - Algorithm 5 subtree configurations at component roots
  - Algorithm 6 subtree configurations at critical vertices
- bottom-up bounded-instance DP driver
- reconstruction of tours and walks on the original input tree

Tested and instrumented:

- local-configuration count sanity checks against the paper definitions
- subtree-configuration count sanity checks on simple cases
- count-only checks for the larger timing-out local-configuration examples documented in
  [README.local-config-debug.md](/Users/arvindmurty/CS/CS583/TUCVRP/README.local-config-debug.md)

Practical caveat:

- the paper solver is implemented end-to-end and returns concrete tours, but it is still slow on some modest random instances, especially for smaller `epsilon`
- the current bottleneck is the local-configuration computation in Algorithm 2
- the project is in a good correctness/debugging state, but not yet in a performance-tuned state

## CLI

The CLI currently runs:

- the paper solver with `epsilon = 0.25`
- the edge-load lower bound
- the exact solver when the instance has at most `20` terminals

Example:

```bash
./bin/tucvrp data/examples/path_small.txt
```

The output includes:

- the pretty-printed rooted tree
- terminal count
- total demand
- edge-load lower bound
- exact optimum or a skip message
- paper solver cost
- exact tours and explicit walks
- paper-solver tours and explicit walks

## Instance Format

The CLI reads a simple text format:

```text
n depot
u0 v0 w0
u1 v1 w1
...
u(n-2) v(n-2) w(n-2)
m
terminal_0 demand_0
...
terminal_(m-1) demand_(m-1)
```

Example:

```text
4 0
0 1 1.0
1 2 2.0
1 3 3.0
2
2 0.4
3 0.6
```

## Example Instances

Small hand-checkable examples are in `data/examples/`:

- `path_small.txt`
- `binary_balanced.txt`
- `star_heavy_leaf.txt`
- `multi_tour_mixed.txt`
- `high_degree_for_binary.txt`
- `lower_bound_gap.txt`

Example usage:

```bash
./bin/tucvrp data/examples/multi_tour_mixed.txt
```

## Code Structure

Main public headers:

- [instance.hpp](/Users/arvindmurty/CS/CS583/TUCVRP/include/tucvrp/instance.hpp)
  tree model, parsing, distance utilities, tour-cost/walk helpers
- [solution.hpp](/Users/arvindmurty/CS/CS583/TUCVRP/include/tucvrp/solution.hpp)
  public `Tour` / `SolveResult` types
- [exact_solver.hpp](/Users/arvindmurty/CS/CS583/TUCVRP/include/tucvrp/exact_solver.hpp)
  exact exponential solver
- [preprocessing.hpp](/Users/arvindmurty/CS/CS583/TUCVRP/include/tucvrp/preprocessing.hpp)
  normalization and lower-bound utilities
- [rooted_tree.hpp](/Users/arvindmurty/CS/CS583/TUCVRP/include/tucvrp/rooted_tree.hpp)
  rooted-tree annotations derived from an instance
- [decomposition.hpp](/Users/arvindmurty/CS/CS583/TUCVRP/include/tucvrp/decomposition.hpp)
  component/block/cluster/cell and height-reduction data structures
- [one_point_five_approx.hpp](/Users/arvindmurty/CS/CS583/TUCVRP/include/tucvrp/algorithms/one_point_five_approx.hpp)
  paper solver API plus exposed DP table builders used by tests

Decomposition implementation is split by paper layer under `src/decomposition/`:

- `components.cpp`
- `blocks.cpp`
- `clusters.cpp`
- `cells.cpp`
- `height_reduction.cpp`

## Testing

Run the full suite with:

```bash
make test
```

The tests cover:

- instance parsing and utilities
- preprocessing invariants
- rooted-tree construction
- decomposition invariants
- height reduction
- RNG reproducibility
- exact solver behavior
- paper solver DP phases and reconstruction

Two useful test/debug files:

- [tests/one_point_five_approx_test.cpp](/Users/arvindmurty/CS/CS583/TUCVRP/tests/one_point_five_approx_test.cpp)
  detailed tests for local and subtree configurations, DP phase interaction, and tour reconstruction
- [README.local-config-debug.md](/Users/arvindmurty/CS/CS583/TUCVRP/README.local-config-debug.md)
  notes on local-configuration counting, including the reproduced timing-out trees

## What The Paper Solver Returns

`OnePointFiveApproxSolver::solve(...)` returns a `SolveResult` on the original input tree.

Each `Tour` contains:

- `terminals`
  the original terminal vertices assigned to that tour
- `demand`
  the total terminal demand covered by that tour
- `cost`
  the route cost
- `walk`
  an explicit depot-to-depot walk on the original tree

The paper solver internally works on normalized / bounded / height-reduced structures, but the public result is projected back to the original instance before being returned.

## Known Limitations

- The exact solver is intentionally limited to small instances and is skipped by the CLI above `20` terminals.
- The approximation solver is functionally implemented, but not yet performance-tuned.
- The current implementation is best treated as an instrumented research implementation, not a production solver.

## Next Steps

Highest-value remaining work:

1. optimize Algorithm 2 local-configuration evaluation
2. add more end-to-end randomized benchmark tooling
3. continue proof-alignment review of the DP state semantics and approximation guarantee
