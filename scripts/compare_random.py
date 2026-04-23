#!/usr/bin/env python3

from __future__ import annotations

import argparse
import random
import subprocess

from compare_common import generate_random_tree_instance, ratio, run_generated_instance


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate random small tree instances and compare exact, Labbé, and paper-solver costs. "
            "The script applies a per-instance timeout because the current paper solver can still "
            "be slow on some modest random instances."
        )
    )
    parser.add_argument("--seed", type=int, default=12345, help="Random seed for instance generation.")
    parser.add_argument(
        "--instances",
        type=int,
        default=10,
        help="Number of random instances to generate.",
    )
    parser.add_argument(
        "--vertices",
        type=int,
        default=10,
        help="Number of vertices in each generated tree.",
    )
    parser.add_argument(
        "--terminals",
        type=int,
        default=5,
        help="Number of terminals in each generated tree.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=20.0,
        help="Per-instance timeout in seconds when invoking bin/tucvrp.",
    )
    args = parser.parse_args()

    rng = random.Random(args.seed)

    print(
        f"{'instance':12} {'exact':>10} {'labbe':>10} {'paper':>10} "
        f"{'labbe/exact':>12} {'paper/exact':>12} {'status':>10}"
    )
    print("-" * 82)

    completed = 0
    labbe_ratios: list[float] = []
    paper_ratios: list[float] = []

    for index in range(args.instances):
        instance_name = f"random_{index:02d}.txt"
        instance_text = generate_random_tree_instance(rng, args.vertices, args.terminals)
        try:
            result = run_generated_instance(instance_text, instance_name, timeout_seconds=args.timeout)
        except subprocess.TimeoutExpired:
            print(f"{instance_name:12} {'-':>10} {'-':>10} {'-':>10} {'-':>12} {'-':>12} {'timeout':>10}")
            continue

        completed += 1
        if result.exact_cost is not None and result.exact_cost > 0.0:
            labbe_ratios.append(result.labbe_cost / result.exact_cost)
            paper_ratios.append(result.paper_cost / result.exact_cost)

        exact_text = "-" if result.exact_cost is None else f"{result.exact_cost:.3f}"
        print(
            f"{result.instance_name:12} "
            f"{exact_text:>10} "
            f"{result.labbe_cost:10.3f} "
            f"{result.paper_cost:10.3f} "
            f"{ratio(result.labbe_cost, result.exact_cost):>12} "
            f"{ratio(result.paper_cost, result.exact_cost):>12} "
            f"{'ok':>10}"
        )

    print()
    print(f"completed: {completed}/{args.instances}")
    if labbe_ratios:
        print(f"labbe ratio mean/max: {sum(labbe_ratios) / len(labbe_ratios):.4f} / {max(labbe_ratios):.4f}")
    if paper_ratios:
        print(f"paper ratio mean/max: {sum(paper_ratios) / len(paper_ratios):.4f} / {max(paper_ratios):.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
