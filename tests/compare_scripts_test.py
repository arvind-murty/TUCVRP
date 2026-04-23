#!/usr/bin/env python3

"""Lightweight tests for the solver-comparison helper scripts."""

from __future__ import annotations

import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from compare_common import generate_random_tree_instance, parse_output, ratio  # noqa: E402


def test_parse_output_with_exact() -> None:
    output = """
edge-load lower bound: 12.000000
exact optimum: 14.000000
labbe baseline cost: 16.000000
paper solver cost: 15.000000
"""
    run = parse_output("sample.txt", output)
    assert run.instance_name == "sample.txt"
    assert run.lower_bound == 12.0
    assert run.exact_cost == 14.0
    assert run.labbe_cost == 16.0
    assert run.paper_cost == 15.0


def test_parse_output_with_skipped_exact() -> None:
    output = """
edge-load lower bound: 12.000000
exact optimum: skipped (more than 20 terminals)
labbe baseline cost: 16.000000
paper solver cost: 15.000000
"""
    run = parse_output("sample.txt", output)
    assert run.exact_cost is None
    assert ratio(run.labbe_cost, run.exact_cost) == "-"


def test_generate_random_tree_instance_has_expected_shape() -> None:
    import random

    text = generate_random_tree_instance(random.Random(7), vertex_count=8, terminal_count=4)
    lines = [line.strip() for line in text.splitlines() if line.strip()]

    assert lines[0] == "8 0"
    assert len(lines) == 1 + 7 + 1 + 4
    assert lines[8] == "4"

    edge_pattern = re.compile(r"^\d+ \d+ \d+$")
    demand_pattern = re.compile(r"^\d+ 0\.\d{3}$")

    for line in lines[1:8]:
        assert edge_pattern.match(line)
    for line in lines[9:]:
        assert demand_pattern.match(line)


def main() -> int:
    test_parse_output_with_exact()
    test_parse_output_with_skipped_exact()
    test_generate_random_tree_instance_has_expected_shape()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
