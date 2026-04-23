#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from compare_common import REPO_ROOT, ratio, run_instance


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare exact, Labbé, and paper-solver costs on checked-in example instances."
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=20.0,
        help="Per-instance timeout in seconds when invoking bin/tucvrp.",
    )
    args = parser.parse_args()

    example_dir = REPO_ROOT / "data" / "examples"
    instance_paths = sorted(example_dir.glob("*.txt"))
    if not instance_paths:
        raise SystemExit(f"no example instances found in {example_dir}")

    print(
        f"{'instance':28} {'lb':>10} {'exact':>10} {'labbe':>10} {'paper':>10} "
        f"{'labbe/exact':>12} {'paper/exact':>12}"
    )
    print("-" * 98)

    for instance_path in instance_paths:
        result = run_instance(instance_path, timeout_seconds=args.timeout)
        exact_text = "-" if result.exact_cost is None else f"{result.exact_cost:.3f}"
        print(
            f"{result.instance_name:28} "
            f"{result.lower_bound:10.3f} "
            f"{exact_text:>10} "
            f"{result.labbe_cost:10.3f} "
            f"{result.paper_cost:10.3f} "
            f"{ratio(result.labbe_cost, result.exact_cost):>12} "
            f"{ratio(result.paper_cost, result.exact_cost):>12}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
