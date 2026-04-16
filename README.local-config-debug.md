# Local Configuration Debug Notes

This note records the end-to-end investigation into the slow cases of
`OnePointFiveApproxSolver`, focusing on the local configuration dynamic program
from Definition 34 / Algorithm 2 of
["A Tight (1.5 + ε)-Approximation for Unsplittable Capacitated Vehicle Routing on Trees"](https://arxiv.org/pdf/2202.05691).

The goal of this note is narrow:
- identify whether the local-configuration state space is genuinely large on the
  timing-out instances, or whether the implementation is producing too many states
- document the exact trees used during the investigation
- record the theoretical counts implied by the paper and compare them to the
  implementation on small sanity-check instances

## Summary

The main findings are:

1. On small sanity-check trees, the implementation's local-configuration count
   matches the paper exactly.
2. On the reproduced timing-out random trees, the number of local configurations
   implied by the paper is already very large:
   - `80,729`
   - `169,910`
   - `201,375`
3. The timing problem is therefore not explained by the implementation inventing
   obviously bogus extra local configurations.
4. The hotspot is still `compute_local_configurations(...)`: even when the state
   count is theoretically correct, the current implementation evaluates all local
   configurations and then brute-forces all labeled partitions for each one.

## Paper Definitions Used

For a component `c`, the paper defines:

- `Q_c`:
  the partition of the terminals of component `c` such that each part is either
  a single big terminal or all small terminals in a cell in `c`
- `Y_c`:

  `Y_c = {α} ∪ { demand(Q~_c) : Q~_c ⊆ Q_c } ∩ (α, 1]`

  where
  `α = ε^(1/ε + 1)`

- A local configuration `(c, A)`:
  a list
  `A = ((y_1, b_1), ..., (y_l, b_l))`
  such that:
  - `l <= |Q_c|`
  - each `y_i ∈ Y_c`
  - each `b_i ∈ {passing, ending}`

For a leaf component, only `ending` is meaningful, because there is no exit
vertex.

## Counting Local Configurations

### Leaf component

If `c` is a leaf component, every entry in a local configuration is of the form:

- `(y_i, ending)`

So a local configuration is just a multiset of size `l` from the set `Y_c`, with
`1 <= l <= |Q_c|`.

If `|Y_c| = m` and `|Q_c| = q`, then the number of local configurations is:

`sum_{l=1..q} C(m + l - 1, l)`

This is the standard count of multisets of size `l` chosen from `m` values.

### Internal component

If `c` is internal, each entry also carries a label:

- `passing`
- `ending`

So each `y ∈ Y_c` yields two pair-types:

- `(y, passing)`
- `(y, ending)`

Ignoring semantic infeasibility for the moment, that gives `2|Y_c|` pair-types.
The same multiset formula then gives:

`sum_{l=1..q} C(2m + l - 1, l)`

where:
- `m = |Y_c|`
- `q = |Q_c|`

In the timing-out cases below, the component is a leaf component, so only the
leaf formula is relevant.

## Sanity-Check Trees

These were used to verify that the implementation is not obviously producing the
wrong number of local configurations on very small components.

### Tree S1: one terminal

```text
2 0
0 1 2
1
1 0.2
```

With `ε = 0.25`:

- `α = 0.25^5 = 0.0009765625`
- `Q_c = { {1} }`
- therefore `|Q_c| = 1`
- subset demands of `Q_c` are:
  - `0.2`
- so:
  - `Y_c = { α, 0.2 }`
  - `|Y_c| = 2`

Since this is a leaf component, the theoretical number of local configurations is:

- for `l = 1`: `C(2 + 1 - 1, 1) = C(2,1) = 2`

Total:

- theoretical count: `2`
- implementation count: `2`

Implementation output:

- `parts = 1`
- `y = 2`
- `values = 2`

### Tree S2: two terminals, one leaf component

```text
3 0
0 1 1
0 2 1
2
1 0.2
2 0.3
```

With `ε = 0.25`:

- `α = 0.0009765625`
- `Q_c = { {1}, {2} }`
- therefore `|Q_c| = 2`
- subset demands are:
  - `0.2`
  - `0.3`
  - `0.5`
- so:
  - `Y_c = { α, 0.2, 0.3, 0.5 }`
  - `|Y_c| = 4`

Leaf-component count:

- `l = 1`: `C(4 + 1 - 1, 1) = C(4,1) = 4`
- `l = 2`: `C(4 + 2 - 1, 2) = C(5,2) = 10`

Total:

- theoretical count: `14`
- implementation count: `14`

Implementation output:

- `parts = 2`
- `y = 4`
- `values = 14`

### Tree S3: one terminal deeper in the tree

```text
4 0
0 1 1
1 2 1
1 3 1
1
3 0.2
```

The local-state calculation is the same as Tree S1:

- `|Q_c| = 1`
- `|Y_c| = 2`
- theoretical count: `2`
- implementation count: `2`

Implementation output:

- `parts = 1`
- `y = 2`
- `values = 2`

## Timing-Out Trees

These trees were generated randomly during empirical end-to-end testing. The
approximation solver timed out externally while being evaluated at `ε = 0.25`.

The key observation is that all of them decompose into a single leaf component,
so the theoretical count of local configurations is exactly the leaf-component
formula above.

### Timeout Tree T1

```text
10 0
0 1 3
0 2 7
1 3 6
2 4 3
4 5 8
0 6 8
1 7 3
7 8 2
0 9 3
5
1 0.177
2 0.207
3 0.410
8 0.401
9 0.205
```

After preprocessing and decomposition at `ε = 0.25`:

- number of components: `1`
- that component is a leaf component
- `|Q_c| = 5`
- `|Y_c| = 26`

Theoretical local-configuration count:

`sum_{l=1..5} C(26 + l - 1, l)`

which is:

- `C(26,1) = 26`
- `C(27,2) = 351`
- `C(28,3) = 3,276`
- `C(29,4) = 23,751`
- `C(30,5) = 142,506`

Total:

- theoretical count: `169,910`

Implementation-side shape measurement:

- component `0`
- `exit = -1`
- `|Q_c| = 5`
- `|Y_c| = 26`

So the implementation is consistent with the paper-level shape.

### Timeout Tree T2

```text
10 0
0 1 1
1 2 6
0 3 6
1 4 9
2 5 6
4 6 7
3 7 6
4 8 3
5 9 9
5
1 0.215
4 0.539
5 0.175
7 0.227
8 0.225
```

Measured shape after preprocessing and decomposition at `ε = 0.25`:

- one leaf component
- `|Q_c| = 5`
- `|Y_c| = 27`

Theoretical count:

`sum_{l=1..5} C(27 + l - 1, l)`

which is:

- `C(27,1) = 27`
- `C(28,2) = 378`
- `C(29,3) = 3,654`
- `C(30,4) = 27,405`
- `C(31,5) = 169,911`

Total:

- theoretical count: `201,375`

Implementation-side shape measurement:

- component `0`
- `exit = -1`
- `|Q_c| = 5`
- `|Y_c| = 27`

Again, this matches the theoretical shape.

### Timeout Tree T3

```text
10 0
0 1 2
0 2 8
0 3 5
3 4 4
4 5 3
1 6 9
4 7 8
6 8 3
7 9 2
5
2 0.417
3 0.127
4 0.238
5 0.275
8 0.524
```

Measured shape after preprocessing and decomposition at `ε = 0.25`:

- one leaf component
- `|Q_c| = 5`
- `|Y_c| = 22`

Theoretical count:

`sum_{l=1..5} C(22 + l - 1, l)`

which is:

- `C(22,1) = 22`
- `C(23,2) = 253`
- `C(24,3) = 2,024`
- `C(25,4) = 12,650`
- `C(26,5) = 65,780`

Total:

- theoretical count: `80,729`

Implementation-side shape measurement:

- component `0`
- `exit = -1`
- `|Q_c| = 5`
- `|Y_c| = 22`

Again, the shape matches.

## How the Measurements Were Obtained

### Small sanity-check trees

For the sanity cases:

1. Construct a small instance by hand.
2. Run preprocessing and component decomposition.
3. Call `compute_local_configurations(...)`.
4. Compare:
   - the number of parts in `Q_c`
   - the number of values in `Y_c`
   - the number of local configurations produced
5. Cross-check with the formula from the paper.

### Timeout trees

For the timeout trees:

1. Generate random tree instances.
2. Run the full approximation solver externally with a timeout.
3. For instances that time out:
   - run preprocessing and decomposition only
   - compute `Q_c` using the paper definition:
     - one part per big terminal
     - one part per cell containing small terminals
   - compute `Y_c` as all subset sums in `(α, 1]` together with `α`
   - use the multiset formula for a leaf component
4. Compare that count to the observed implementation-side shape:
   - `|Q_c|`
   - `|Y_c|`

### Important note

For the timeout cases, I did **not** call `compute_local_configurations(...)` to count the values,
because that is exactly the slow routine being diagnosed.

Instead, I reconstructed only:

- `|Q_c|`
- `|Y_c|`
- the theoretical number of local configurations

using the paper definitions directly.

## Conclusions

The evidence from these checks is:

1. On simple trees, the implementation's local-configuration count matches the paper.
2. On the timing-out trees, the paper-level local-state count is already large.
3. Therefore, the local slowdown is **not** explained by the implementation producing clearly
   incorrect extra local configurations.

That does **not** prove the implementation is optimal. The current Algorithm 2 implementation
still evaluates all local configurations and then brute-forces all labeled partitions for each of
them, which is likely much slower than necessary in practice. But the state-space explosion in the
timing-out examples appears to be real, not obviously a counting bug.

## Current Hypothesis

The main performance problem is:

- `compute_local_configurations(...)`

and specifically the combination of:

- a large `|Y_c|`
- a large number of local configurations implied by the paper
- brute-force evaluation of all labeled partitions for each local configuration

So the next debugging/optimization step should focus on:

- pruning dominated local configurations
- reusing work across local configurations
- checking whether Algorithm 2 can be implemented more efficiently than the current literal
  brute-force loop
