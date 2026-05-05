# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# A/B comparator for TS1-cloud LHS stress metrics.
#
# Loads two metrics JSON files produced by the LHS stress harness
# (TS1_CLOUD_STRESS_METRICS_OUT) and prints a delta table covering convergence,
# step distribution, solver-cost aggregates, and a per-bucket histogram diff
# of warm-KPI step counts.

from __future__ import annotations

import argparse
import json
import math
import sys
from typing import Any, Dict, Tuple


# Regression tolerances used by the optional pytest entry. A regression
# is signaled when ANY tracked metric worsens by more than its tolerance.
DEFAULT_TOLERANCES: Dict[str, float] = {
    "total_wasted_work_rate": 0.005,   # +0.5 pp on the headline scalar
    "cold_failure_rate": 0.005,        # +0.5 percentage points allowed
    "warm_failure_rate": 0.002,        # +0.2 percentage points allowed
    "warm_median_steps": 2.0,          # +2 internal steps allowed
    "warm_p95_steps": 5.0,             # +5 internal steps allowed
    "warm_p99_steps": 10.0,            # +10 internal steps allowed
    "warm_high_step_rate": 0.02,       # +2 percentage points allowed
    "warm_kpi_aggregates.mean_function_calls": 50.0,
    "warm_kpi_aggregates.mean_jacobian_updates": 5.0,
    "warm_kpi_aggregates.mean_rejected_fraction": 0.02,
    "warm_kpi_aggregates.mean_ns_per_step": 5.0e5,  # 0.5 ms/step jitter
}


def _load(path: str) -> Dict[str, Any]:
    with open(path) as f:
        return json.load(f)


def _get(d: Dict[str, Any], dotted: str) -> Any:
    cur: Any = d
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return float("nan")
        cur = cur[part]
    return cur


def _delta(a: Any, b: Any) -> float:
    try:
        af = float(a)
        bf = float(b)
    except (TypeError, ValueError):
        return float("nan")
    if math.isnan(af) and math.isnan(bf):
        return 0.0
    if math.isnan(af) or math.isnan(bf):
        return float("nan")
    return bf - af


def _fmt(v: Any) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if math.isnan(f):
        return "nan"
    if abs(f) >= 1.0e6 or (0.0 < abs(f) < 1.0e-3):
        return f"{f:.3e}"
    return f"{f:.4g}"


_HEADLINE_METRICS = (
    "total_wasted_work_rate",
    "cold_failure_rate",
    "warm_failure_rate",
    "cold_median_steps",
    "cold_p95_steps",
    "warm_median_steps",
    "warm_p90_steps",
    "warm_p95_steps",
    "warm_p99_steps",
    "warm_max_steps",
    "warm_high_step_rate",
    "warm_kpi_aggregates.mean_function_calls",
    "warm_kpi_aggregates.mean_jacobian_updates",
    "warm_kpi_aggregates.mean_decompositions",
    "warm_kpi_aggregates.mean_solves",
    "warm_kpi_aggregates.mean_rejected_fraction",
    "warm_kpi_aggregates.mean_ns_per_step",
    "warm_kpi_aggregates.median_elapsed_ns",
    "warm_kpi_aggregates.p95_elapsed_ns",
)


def compute_delta(baseline: Dict[str, Any], candidate: Dict[str, Any]) -> Dict[str, Any]:
    """Compute structured A/B delta record."""
    deltas: Dict[str, Any] = {"metrics": {}}
    for key in _HEADLINE_METRICS:
        a = _get(baseline, key)
        b = _get(candidate, key)
        deltas["metrics"][key] = {
            "baseline": a,
            "candidate": b,
            "delta": _delta(a, b),
        }
    base_hist = baseline.get("warm_step_histogram", {})
    cand_hist = candidate.get("warm_step_histogram", {})
    keys = list(base_hist.keys()) + [k for k in cand_hist.keys() if k not in base_hist]
    deltas["histogram"] = {
        k: {
            "baseline": int(base_hist.get(k, 0)),
            "candidate": int(cand_hist.get(k, 0)),
            "delta": int(cand_hist.get(k, 0)) - int(base_hist.get(k, 0)),
        }
        for k in keys
    }
    deltas["context"] = {
        "baseline_label": baseline.get("label"),
        "candidate_label": candidate.get("label"),
        "baseline_samples": baseline.get("samples"),
        "candidate_samples": candidate.get("samples"),
        "seed_match": baseline.get("seed") == candidate.get("seed"),
    }
    return deltas


def find_regressions(
    deltas: Dict[str, Any],
    tolerances: Dict[str, float] = DEFAULT_TOLERANCES,
) -> Tuple[bool, list]:
    """Return (regressed, list_of_messages). Higher metric values are worse."""
    bad = []
    for key, tol in tolerances.items():
        rec = deltas["metrics"].get(key)
        if rec is None:
            continue
        d = rec["delta"]
        if isinstance(d, float) and math.isnan(d):
            continue
        if d > tol:
            bad.append(
                f"{key}: baseline={_fmt(rec['baseline'])} -> candidate={_fmt(rec['candidate'])} "
                f"delta=+{_fmt(d)} (tol +{_fmt(tol)})"
            )
    return (len(bad) > 0, bad)


def format_report(deltas: Dict[str, Any]) -> str:
    lines = []
    ctx = deltas["context"]
    lines.append(
        f"A/B context: baseline='{ctx['baseline_label']}' "
        f"({ctx['baseline_samples']} samples) vs candidate='{ctx['candidate_label']}' "
        f"({ctx['candidate_samples']} samples), seed_match={ctx['seed_match']}"
    )
    lines.append("")
    lines.append(f"{'metric':<48} {'baseline':>14} {'candidate':>14} {'delta':>14}")
    lines.append("-" * 92)
    for key, rec in deltas["metrics"].items():
        lines.append(
            f"{key:<48} {_fmt(rec['baseline']):>14} {_fmt(rec['candidate']):>14} {_fmt(rec['delta']):>14}"
        )
    lines.append("")
    lines.append("warm_step_histogram (warm KPI samples per bucket)")
    lines.append(f"{'bucket':<14} {'baseline':>10} {'candidate':>10} {'delta':>10}")
    lines.append("-" * 48)
    for key, rec in deltas["histogram"].items():
        sign = "+" if rec["delta"] > 0 else ""
        lines.append(f"{key:<14} {rec['baseline']:>10} {rec['candidate']:>10} {sign}{rec['delta']:>9}")
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", help="Path to baseline metrics JSON")
    parser.add_argument("candidate", help="Path to candidate metrics JSON")
    parser.add_argument(
        "--check-regression",
        action="store_true",
        help="Exit nonzero if any tracked metric regresses beyond tolerance.",
    )
    args = parser.parse_args(argv)

    deltas = compute_delta(_load(args.baseline), _load(args.candidate))
    print(format_report(deltas))

    if args.check_regression:
        regressed, messages = find_regressions(deltas)
        if regressed:
            print("\nRegressions detected:")
            for m in messages:
                print(f"  - {m}")
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
