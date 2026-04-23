#!/usr/bin/env python3

"""Shared helpers for comparing exact, Labbé, and paper-solver outputs."""

from __future__ import annotations

import random
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


REPO_ROOT = Path(__file__).resolve().parent.parent
BIN = REPO_ROOT / "bin" / "tucvrp"

EXACT_PATTERN = re.compile(r"^exact optimum: (?P<value>.+)$", re.MULTILINE)
LABBE_PATTERN = re.compile(r"^labbe baseline cost: (?P<value>[0-9.]+)$", re.MULTILINE)
PAPER_PATTERN = re.compile(r"^paper solver cost: (?P<value>[0-9.]+)$", re.MULTILINE)
LB_PATTERN = re.compile(r"^edge-load lower bound: (?P<value>[0-9.]+)$", re.MULTILINE)


@dataclass
class SolverRun:
    instance_name: str
    lower_bound: float
    exact_cost: Optional[float]
    labbe_cost: float
    paper_cost: float
    raw_output: str


def ensure_binary() -> None:
    """Build the CLI lazily so the comparison scripts can be run from a clean checkout."""
    if BIN.exists():
        return
    subprocess.run(["make"], cwd=REPO_ROOT, check=True)


def parse_output(instance_name: str, output: str) -> SolverRun:
    """Parse the current CLI summary block into a structured record."""
    lb_match = LB_PATTERN.search(output)
    labbe_match = LABBE_PATTERN.search(output)
    paper_match = PAPER_PATTERN.search(output)
    exact_match = EXACT_PATTERN.search(output)

    if lb_match is None or labbe_match is None or paper_match is None or exact_match is None:
        raise ValueError(f"failed to parse solver output for {instance_name}")

    exact_text = exact_match.group("value").strip()
    exact_cost: Optional[float]
    if exact_text.startswith("skipped"):
        exact_cost = None
    else:
        exact_cost = float(exact_text)

    return SolverRun(
        instance_name=instance_name,
        lower_bound=float(lb_match.group("value")),
        exact_cost=exact_cost,
        labbe_cost=float(labbe_match.group("value")),
        paper_cost=float(paper_match.group("value")),
        raw_output=output,
    )


def run_instance(instance_path: Path, timeout_seconds: float) -> SolverRun:
    """Run the CLI on one existing instance file and parse the result."""
    ensure_binary()
    completed = subprocess.run(
        [str(BIN), str(instance_path)],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    return parse_output(instance_path.name, completed.stdout)


def ratio(numerator: float, denominator: Optional[float]) -> str:
    """Format a comparison ratio, leaving skipped exact runs blank."""
    if denominator is None or denominator == 0.0:
        return "-"
    return f"{numerator / denominator:.4f}"


def generate_random_tree_instance(
    rng: random.Random,
    vertex_count: int,
    terminal_count: int,
) -> str:
    """Generate a small exact-solvable random tree instance in the project's text format."""
    if vertex_count < 2:
        raise ValueError("vertex_count must be at least 2")
    if terminal_count < 1 or terminal_count >= vertex_count:
        raise ValueError("terminal_count must be in [1, vertex_count - 1]")

    lines = [f"{vertex_count} 0"]

    for vertex in range(1, vertex_count):
        parent = rng.randrange(0, vertex)
        weight = rng.randint(1, 9)
        lines.append(f"{parent} {vertex} {weight}")

    terminals = sorted(rng.sample(range(1, vertex_count), terminal_count))
    lines.append(str(terminal_count))
    for terminal in terminals:
        # Keep demands comfortably below 1 so instances remain feasible but nontrivial.
        demand = rng.uniform(0.10, 0.55)
        lines.append(f"{terminal} {demand:.3f}")

    return "\n".join(lines) + "\n"


def run_generated_instance(
    text: str,
    instance_name: str,
    timeout_seconds: float,
) -> SolverRun:
    """Run the CLI on one generated instance text without leaving files behind."""
    ensure_binary()
    with tempfile.TemporaryDirectory() as tmpdir:
        instance_path = Path(tmpdir) / instance_name
        instance_path.write_text(text, encoding="utf-8")
        completed = subprocess.run(
            [str(BIN), str(instance_path)],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
        return parse_output(instance_name, completed.stdout)
