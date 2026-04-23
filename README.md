# TUCVRP

This repository contains a C++20 implementation of routing algorithms for the
tree version of the Unsplittable Capacitated Vehicle Routing Problem (UCVRP).

The current codebase includes:

- an exact solver for small instances
- an implementation of the modern `(1.5 + ε)` approximation framework
- a simpler Labbé-style baseline heuristic adapted to the same model
- preprocessing, decomposition, height reduction, and dynamic-programming infrastructure
- scripts for comparing solver outputs on hand-built and random instances

The repository is best thought of as an implementation and experimentation
workspace for approximation algorithms on trees. It is in good shape for
inspection, testing, and solver comparison, but the `(1.5 + ε)` solver is still
slow on some modest instances and should not be treated as a production-tuned tool.

## Problem Model

The repository works with the following normalized model:

- the input is a weighted tree rooted at a depot
- each terminal has demand in `(0, 1]`
- each vehicle has capacity `1`
- terminal demands are unsplittable
- the objective is to minimize total route length

For a set of terminals, the cheapest one-tour route is the round trip over the
minimal subtree containing those terminals and the depot. The codebase uses that
fact throughout both the exact and approximate solvers.

## Implemented Solvers

### 1. Exact Solver

`ExactSolver` solves small instances exactly by subset dynamic programming.

- intended for instances with at most `20` terminals
- useful as a correctness oracle
- used heavily by tests and comparison scripts

Relevant files:

- `include/tucvrp/exact_solver.hpp`
- `src/exact_solver.cpp`

### 2. Labbé Baseline

`LabbeApproxSolver` is a deterministic adaptation of the older tree heuristic to
the current one-vehicle-type model used in this repository.

- works directly on the original rooted tree
- repeatedly eliminates leaves bottom-up
- produces explicit tours and walks on the original instance
- currently serves as the baseline for comparison against the `(1.5 + ε)` solver

Relevant files:

- `include/tucvrp/algorithms/labbe_approx.hpp`
- `src/algorithms/labbe_approx.cpp`
- `tests/labbe_approx_test.cpp`

### 3. `(1.5 + ε)` Paper Solver

`OnePointFiveApproxSolver` implements the main pipeline behind the
`(1.5 + ε)` approximation framework used in this project.

The implementation currently includes:

- preprocessing to a binary tree with leaf terminals
- bounded-distance outer reduction
- component / block / cluster / cell decomposition
- height reduction on the component tree
- local-configuration DP
- subtree-configuration DP
- reconstruction of explicit depot-to-depot tours on the original input tree

Relevant files:

- `include/tucvrp/algorithms/one_point_five_approx.hpp`
- `src/algorithms/one_point_five_approx.cpp`
- `src/decomposition/`
- `tests/one_point_five_approx_test.cpp`

## Repository Layout

### Public headers

- `include/tucvrp/instance.hpp`
  - tree instance model
  - parser
  - rooted-tree utilities
  - route cost and route walk helpers
- `include/tucvrp/solution.hpp`
  - `Tour` and `SolveResult`
- `include/tucvrp/preprocessing.hpp`
  - lower bound and normalization
- `include/tucvrp/rooted_tree.hpp`
  - cached rooted-tree annotations
- `include/tucvrp/decomposition.hpp`
  - decomposition and height-reduction data structures
- `include/tucvrp/exact_solver.hpp`
- `include/tucvrp/algorithms/labbe_approx.hpp`
- `include/tucvrp/algorithms/one_point_five_approx.hpp`

### Source tree

- `src/instance.cpp`
  - parsing, validation, rooted tree utilities, tour support computation
- `src/preprocessing.cpp`
  - binary/leaf normalization and edge-load lower bound
- `src/rooted_tree.cpp`
  - rooted-tree annotations
- `src/exact_solver.cpp`
  - exact subset DP
- `src/algorithms/labbe_approx.cpp`
  - baseline heuristic
- `src/algorithms/one_point_five_approx.cpp`
  - main approximation solver and DP logic
- `src/decomposition/`
  - decomposition implementation split by layer:
    - `components.cpp`
    - `blocks.cpp`
    - `clusters.cpp`
    - `cells.cpp`
    - `height_reduction.cpp`

### Tests

- `tests/instance_test.cpp`
- `tests/preprocessing_test.cpp`
- `tests/rooted_tree_test.cpp`
- `tests/decomposition_test.cpp`
- `tests/rng_test.cpp`
- `tests/labbe_approx_test.cpp`
- `tests/one_point_five_approx_test.cpp`
- `tests/compare_scripts_test.py`

### Scripts

- `scripts/compare_examples.py`
  - compare exact / Labbé / paper solver on checked-in example instances
- `scripts/compare_random.py`
  - generate random exact-solvable instances and compare the three solvers
- `scripts/compare_common.py`
  - shared helper code for script parsing and instance generation

### Data

- `data/examples/`
  - small hand-built example instances used for debugging and smoke testing

## Build and Test

### Requirements

- `clang++` with C++20 support
- `make`
- `pkg-config`
- Catch2 installed in a way visible to `pkg-config`
- Python 3 for the comparison scripts

The `Makefile` discovers Catch2 through `pkg-config`.

### Build

```bash
make
```

This builds the main CLI at:

```bash
bin/tucvrp
```

### Run all tests

```bash
make test
```

This runs:

- the full Catch2 C++ test suite
- the Python helper-script test

The current suite covers:

- instance parsing and validation
- route-cost and route-walk helpers
- preprocessing invariants
- rooted-tree annotations
- decomposition invariants
- height reduction
- exact-solver behavior
- Labbé baseline behavior
- local-configuration and subtree-configuration DP behavior
- bounded solver reconstruction
- comparison-script parsing/generation helpers

## Command-Line Usage

Run the CLI on one instance file:

```bash
./bin/tucvrp data/examples/path_small.txt
```

The CLI currently prints:

- the pretty-printed rooted tree
- terminal count
- total demand
- edge-load lower bound
- exact optimum, when `terminal_count <= 20`
- Labbé baseline cost and tours
- paper solver cost and tours

Each printed tour includes:

- served terminals
- total demand
- route cost
- explicit depot-to-depot walk

## Instance File Format

The input format is:

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

where:

- `n` is the number of vertices
- `depot` is the depot vertex id
- the next `n-1` lines are undirected weighted tree edges
- `m` is the number of terminals
- each terminal line gives `vertex demand`

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

The repository includes small hand-checkable examples:

- `data/examples/path_small.txt`
- `data/examples/binary_balanced.txt`
- `data/examples/star_heavy_leaf.txt`
- `data/examples/multi_tour_mixed.txt`
- `data/examples/high_degree_for_binary.txt`
- `data/examples/lower_bound_gap.txt`

These are useful for:

- inspecting solver output manually
- smoke-testing the CLI
- comparing the exact, Labbé, and paper solvers

## Comparison Scripts

### Compare on checked-in examples

```bash
python3 scripts/compare_examples.py
```

This:

- builds `bin/tucvrp` automatically if needed
- runs it on every file in `data/examples/`
- parses exact / lower-bound / Labbé / paper-solver output
- prints a compact comparison table

### Compare on random instances

```bash
python3 scripts/compare_random.py --instances 10 --vertices 10 --terminals 5 --timeout 20
```

Useful flags:

- `--seed`
- `--instances`
- `--vertices`
- `--terminals`
- `--timeout`

This script is meant for quick solver comparison on small exact-solvable trees.
It uses the existing CLI, not a separate benchmarking binary.

Important note:

- the script applies a per-instance timeout when invoking `bin/tucvrp`
- this is intentional
- the current `(1.5 + ε)` solver can still become slow on some modest random instances
- timeout rows are reported explicitly in the output instead of leaving background jobs running

## How To Work With The Codebase

### If you want to understand the basic data model

Start with:

- `include/tucvrp/instance.hpp`
- `src/instance.cpp`
- `include/tucvrp/solution.hpp`

### If you want to understand the exact baseline

Read:

- `include/tucvrp/exact_solver.hpp`
- `src/exact_solver.cpp`
- `tests/instance_test.cpp`

### If you want to understand the Labbé baseline

Read:

- `include/tucvrp/algorithms/labbe_approx.hpp`
- `src/algorithms/labbe_approx.cpp`
- `tests/labbe_approx_test.cpp`

### If you want to understand the paper solver

The cleanest order is:

1. `include/tucvrp/preprocessing.hpp`
2. `include/tucvrp/rooted_tree.hpp`
3. `include/tucvrp/decomposition.hpp`
4. `src/decomposition/`
5. `include/tucvrp/algorithms/one_point_five_approx.hpp`
6. `src/algorithms/one_point_five_approx.cpp`
7. `tests/one_point_five_approx_test.cpp`

### If you want to understand the performance bottleneck

Read:

- `README.local-config-debug.md`
- `tests/one_point_five_approx_test.cpp`

Those explain and test the local-configuration and subtree-configuration state counts,
including the larger timing-out examples.

## Current Status

What is implemented:

- exact solver
- Labbé baseline solver
- end-to-end `(1.5 + ε)` solver with explicit tour reconstruction
- decomposition hierarchy and height reduction
- comparison scripts and solver-side tests

What is still weak:

- the `(1.5 + ε)` solver is not performance-tuned
- Algorithm 2 local-configuration evaluation is the dominant runtime bottleneck
- the codebase is strongest as a correctness/debugging platform, not yet as a scalable experiment platform

## Known Limitations

- `ExactSolver` is only intended for small instances and is skipped by the CLI above `20` terminals.
- The paper solver can be slow on modest random instances, especially for smaller values of `ε`.
- The comparison scripts are lightweight wrappers around the CLI, not a dedicated benchmark framework.
- The paper draft in `paper/` is still a working draft, not a finalized manuscript.

## Related Documentation

- `README.local-config-debug.md`
  - notes on local-configuration counting and timing-out examples
- `paper/implementation_draft.tex`
  - current draft of the implementation paper

## Next Steps

The most useful next improvements are:

1. optimize Algorithm 2 local-configuration evaluation
2. extend solver comparison experiments
3. continue tightening the implementation paper and bibliography
4. keep auditing the `(1.5 + ε)` DP state semantics against the source papers
