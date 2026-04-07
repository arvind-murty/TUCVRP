# TUCVRP

Initial C++20 scaffold for implementing the paper
"A Tight (1.5 + epsilon)-Approximation for Unsplittable Capacitated Vehicle Routing on Trees"
([arXiv:2202.05691](https://arxiv.org/pdf/2202.05691)).

## Build

This project uses `make` and the system-installed Catch2 package.

```bash
make
make test
```

The current `Makefile` discovers Catch2 through `pkg-config`.

## Current Status

Implemented:

- rooted tree instance model with validation
- simple text parser for instances
- exact exponential solver for small instances, useful as a correctness oracle
- paper-relevant utilities:
- bounded-distance check from Definition 4
- edge-load lower bound
- conversion to a binary tree by expanding high-degree vertices

Scaffolded but not yet paper-complete:

- component decomposition from Lemma 6
- multi-level block/cluster/cell decomposition
- local theorem machinery from Section 5
- adaptive rounding and dynamic program from Sections 8 and 9

The code is structured so those pieces can be added incrementally without rewriting the basics.

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

- `path_small.txt`: a tiny path instance
- `binary_balanced.txt`: a balanced binary tree
- `star_heavy_leaf.txt`: a star with one far expensive terminal
- `multi_tour_mixed.txt`: a small instance that forces multiple tours
- `high_degree_for_binary.txt`: a high-degree tree useful for binary-tree preprocessing
- `lower_bound_gap.txt`: a small instance where the edge-load lower bound is strictly below the optimum

Example usage:

```bash
./bin/tucvrp data/examples/path_small.txt
```

## Next Steps

1. Implement the exact component partition from Lemma 6.
2. Add block/cluster/cell builders with invariants encoded in tests.
3. Represent local solutions explicitly and implement the simplification steps from Section 5.
4. Add the global structural rounding and the Section 9 dynamic program.
