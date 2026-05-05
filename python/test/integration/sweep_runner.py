# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# Parallel sweep driver for TS1-cloud solver-parameter optimization.
#
# Each variant is a dict of solver overrides passed to _create_micm().
# Variants run concurrently in worker processes; each worker writes its own
# metrics JSON, then the driver A/B-compares every variant against two
# reference files:
#   - committed pre-charge-fix baseline (recovery target)
#   - post-charge-fix reference (improvement target)
#
# Usage:
#   uv run python python/test/integration/sweep_runner.py \
#       --plan tolerance_1d \
#       --samples 200 \
#       --workers 12 \
#       --out-dir /tmp/ts1_sweeps
#
# Built-in plans:
#   "tolerance_1d"   1D sweeps over atol_aqueous, atol_gas_alg, atol_gas_diff,
#                    relative_tolerance, h_start, constraint_init_tolerance.
#   "smoke"          A 2-variant smoke-test plan for quick wiring validation.

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Dict, List, Optional

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# Import these eagerly in the parent so import errors surface fast.
from python.test.integration import ab_compare  # noqa: E402

REPO_BASELINE = os.path.join(
    os.path.dirname(__file__), "data", "ts1_cloud_baseline.json"
)
POST_CHARGE_FIX_REF = os.path.join(
    os.path.dirname(__file__), "data", "ts1_cloud_after_charge_fix.json"
)


# ------------------------- worker -------------------------

def _run_variant(
    name: str,
    overrides: Dict[str, Any],
    n_samples: int,
    seed: int,
    warm_solves: int,
    dt: float,
    out_path: str,
) -> Dict[str, Any]:
    """Run one variant in a worker process and emit its metrics file.

    The import of the harness is deferred to the worker so each process
    initializes its own MICM bindings without contention.
    """
    # Lazy import in worker.
    from python.test.integration import test_ts1_cloud_lhs_stress as harness

    t0 = time.perf_counter()
    summary = harness._run_stress(
        n_samples=n_samples,
        seed=seed,
        warm_solves=warm_solves,
        dt=dt,
        overrides=overrides,
    )
    elapsed = time.perf_counter() - t0

    payload = {k: v for k, v in summary.items() if k != "_per_sample"}
    payload["label"] = name
    payload["overrides"] = overrides
    payload["per_sample"] = summary["_per_sample"]
    payload["wallclock_seconds"] = elapsed

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2, default=harness._json_default)

    # Return a compact descriptor for the parent.
    return {
        "name": name,
        "out_path": out_path,
        "wallclock_seconds": elapsed,
        "warm_median_steps": summary["warm_median_steps"],
        "warm_p95_steps": summary["warm_p95_steps"],
        "warm_p99_steps": summary["warm_p99_steps"],
        "warm_max_steps": summary["warm_max_steps"],
        "warm_failure_rate": summary["warm_failure_rate"],
        "cold_failure_rate": summary["cold_failure_rate"],
    }


# ------------------------- plans -------------------------

def plan_smoke() -> List[Dict[str, Any]]:
    return [
        {"name": "smoke_default", "overrides": {}},
        {"name": "smoke_rtol_1e-4", "overrides": {"relative_tolerance": 1e-4}},
    ]


def plan_tolerance_1d() -> List[Dict[str, Any]]:
    """One-dim sweeps; 'default' establishes the per-sweep zero point."""
    variants: List[Dict[str, Any]] = [{"name": "default", "overrides": {}}]

    for v in (1e-10, 1e-9, 1e-7, 1e-6):
        variants.append({"name": f"atol_aq_{v:.0e}", "overrides": {"atol_aqueous": v}})

    for v in (1e-11, 1e-10, 1e-8, 1e-7):
        variants.append({"name": f"atol_gas_alg_{v:.0e}", "overrides": {"atol_gas_alg": v}})

    for v in (1e-5, 1e-4, 1e-2):
        variants.append({"name": f"atol_gas_diff_{v:.0e}", "overrides": {"atol_gas_diff": v}})

    for v in (1e-5, 1e-4, 1e-3):
        variants.append({"name": f"rtol_{v:.0e}", "overrides": {"relative_tolerance": v}})

    for v in (1e-3, 1e-1, 1.0):
        variants.append({"name": f"h_start_{v:.0e}", "overrides": {"h_start": v}})

    for v in (1e-10, 1e-9, 1e-7, 1e-6):
        variants.append({"name": f"cinit_tol_{v:.0e}", "overrides": {"constraint_init_tolerance": v}})

    return variants


def plan_tolerance_2d() -> List[Dict[str, Any]]:
    """2D refinement around cinit_tol=1e-9 (best 1D cold-robustness)."""
    variants: List[Dict[str, Any]] = [
        {"name": "default", "overrides": {}},
        {"name": "cinit_1e-9_solo", "overrides": {"constraint_init_tolerance": 1e-9}},
        {"name": "cinit_1e-10_solo", "overrides": {"constraint_init_tolerance": 1e-10}},
        # cinit=1e-9 paired with each promising knob
        {"name": "cinit_1e-9__atol_aq_1e-9",
         "overrides": {"constraint_init_tolerance": 1e-9, "atol_aqueous": 1e-9}},
        {"name": "cinit_1e-9__atol_aq_1e-10",
         "overrides": {"constraint_init_tolerance": 1e-9, "atol_aqueous": 1e-10}},
        {"name": "cinit_1e-9__rtol_1e-3",
         "overrides": {"constraint_init_tolerance": 1e-9, "relative_tolerance": 1e-3}},
        {"name": "cinit_1e-9__rtol_1e-4",
         "overrides": {"constraint_init_tolerance": 1e-9, "relative_tolerance": 1e-4}},
        {"name": "cinit_1e-9__hstart_1e-1",
         "overrides": {"constraint_init_tolerance": 1e-9, "h_start": 1e-1}},
        {"name": "cinit_1e-9__hstart_1.0",
         "overrides": {"constraint_init_tolerance": 1e-9, "h_start": 1.0}},
        # cinit=1e-10 paired
        {"name": "cinit_1e-10__rtol_1e-3",
         "overrides": {"constraint_init_tolerance": 1e-10, "relative_tolerance": 1e-3}},
        {"name": "cinit_1e-10__atol_aq_1e-9",
         "overrides": {"constraint_init_tolerance": 1e-10, "atol_aqueous": 1e-9}},
        # triples around the strongest cold-robust base
        {"name": "cinit_1e-9__atol_aq_1e-9__rtol_1e-3",
         "overrides": {"constraint_init_tolerance": 1e-9, "atol_aqueous": 1e-9,
                       "relative_tolerance": 1e-3}},
        {"name": "cinit_1e-9__atol_aq_1e-9__hstart_1e-1",
         "overrides": {"constraint_init_tolerance": 1e-9, "atol_aqueous": 1e-9,
                       "h_start": 1e-1}},
        {"name": "cinit_1e-9__atol_aq_1e-9__rtol_1e-3__hstart_1e-1",
         "overrides": {"constraint_init_tolerance": 1e-9, "atol_aqueous": 1e-9,
                       "relative_tolerance": 1e-3, "h_start": 1e-1}},
    ]
    return variants


def plan_min_halflife() -> List[Dict[str, Any]]:
    """Phase 4: per-reaction `min halflife [s]` sweep.

    The TS1-cloud config caps three rate-capped DISSOLVED_REACTIONs at
    `min halflife [s] = 1.0`:
       R1b: SO2OOH- + H+   -> SO4-- + 2 H+   (acid-catalyzed sulfate prod)
       R2:  HSO3- + O3(aq) -> SO4-- + H+     (O3 pathway, low rate)
       R3:  SO3-- + O3(aq) -> SO4--          (O3 pathway, very high k)

    Lower halflife = less aggressive cap (closer to true chemistry, but stiffer
    for the solver). Higher halflife = stronger cap (more numerically stable
    but more chemistry distortion). We sweep each reaction independently to
    find the smallest cap that doesn't degrade total_wasted_work or warm tail
    relative to the Regime B baseline.
    """
    variants: List[Dict[str, Any]] = [{"name": "default", "overrides": {}}]
    for tag in ("R1b", "R2", "R3"):
        for v in (0.1, 0.3, 3.0, 10.0):
            variants.append({
                "name": f"mhl_{tag}_{v:g}",
                "overrides": {"min_halflife": {tag: v}},
            })
    # All-three combos at uniform values (mirrors current default = 1.0).
    for v in (0.1, 0.3, 3.0, 10.0):
        variants.append({
            "name": f"mhl_all_{v:g}",
            "overrides": {"min_halflife": {"R1b": v, "R2": v, "R3": v}},
        })
    return variants


PLANS = {
    "smoke": plan_smoke,
    "tolerance_1d": plan_tolerance_1d,
    "tolerance_2d": plan_tolerance_2d,
    "min_halflife": plan_min_halflife,
}


# ------------------------- A/B reporting -------------------------

PHASE3_RECOVERY_GATES = {
    # Recover original (uncorrected mechanism) tail behavior.
    "warm_median_steps": 9.0 + 1.0,        # baseline 9
    "warm_p95_steps": 311.20 + 1e-6,        # baseline 311.2
    "warm_p99_steps": 675.92 + 1e-6,        # baseline 675.92
    "warm_failure_rate": 0.02269 + 1e-6,    # baseline 2.27%
    "cold_failure_rate": 0.244 + 1e-6,      # baseline 24.4%
}

PHASE3_IMPROVEMENT_GATES = {
    # Strict improvement over post-charge-fix snapshot: ΔP95 ≤ 0 etc.
    "warm_p95_steps": 0.0,
    "warm_p99_steps": 0.0,
    "warm_failure_rate": 1e-6,
    "cold_failure_rate": 1e-6,
}


def _eval_recovery(candidate: Dict[str, Any]) -> Dict[str, Any]:
    """Check candidate against absolute Phase 3 recovery gates."""
    failures = []
    for key, ceiling in PHASE3_RECOVERY_GATES.items():
        v = float(candidate.get(key, float("nan")))
        if v != v:  # NaN
            continue
        if v > ceiling:
            failures.append(f"{key}={v:.4g} > recovery_gate {ceiling:.4g}")
    return {"passed": len(failures) == 0, "violations": failures}


def _eval_improvement(deltas: Dict[str, Any]) -> Dict[str, Any]:
    failures = []
    for key, max_delta in PHASE3_IMPROVEMENT_GATES.items():
        rec = deltas["metrics"].get(key)
        if rec is None:
            continue
        d = rec["delta"]
        if isinstance(d, float) and d != d:
            continue
        if d > max_delta:
            failures.append(f"{key}: delta={d:+.4g} > improvement_gate {max_delta:.4g}")
    return {"passed": len(failures) == 0, "violations": failures}


def summarize_run(
    variant_name: str,
    candidate_path: str,
    pre_baseline: Dict[str, Any],
    post_baseline: Dict[str, Any],
) -> Dict[str, Any]:
    with open(candidate_path) as f:
        candidate = json.load(f)
    deltas_pre = ab_compare.compute_delta(pre_baseline, candidate)
    deltas_post = ab_compare.compute_delta(post_baseline, candidate)
    recovery = _eval_recovery(candidate)
    improvement = _eval_improvement(deltas_post)
    return {
        "name": variant_name,
        "candidate_path": candidate_path,
        "wallclock_seconds": candidate.get("wallclock_seconds"),
        "warm_median_steps": candidate.get("warm_median_steps"),
        "warm_p95_steps": candidate.get("warm_p95_steps"),
        "warm_p99_steps": candidate.get("warm_p99_steps"),
        "warm_max_steps": candidate.get("warm_max_steps"),
        "warm_failure_rate": candidate.get("warm_failure_rate"),
        "cold_failure_rate": candidate.get("cold_failure_rate"),
        "delta_vs_pre": {k: deltas_pre["metrics"][k]["delta"] for k in (
            "warm_median_steps", "warm_p95_steps", "warm_p99_steps", "warm_failure_rate", "cold_failure_rate"
        )},
        "delta_vs_post": {k: deltas_post["metrics"][k]["delta"] for k in (
            "warm_median_steps", "warm_p95_steps", "warm_p99_steps", "warm_failure_rate", "cold_failure_rate"
        )},
        "phase3_recovery_passed": recovery["passed"],
        "phase3_improvement_passed": improvement["passed"],
        "violations_recovery": recovery["violations"],
        "violations_improvement": improvement["violations"],
    }


def _print_table(rows: List[Dict[str, Any]]) -> None:
    cols = [
        ("name", 22),
        ("warm_median_steps", 8),
        ("warm_p95_steps", 9),
        ("warm_p99_steps", 9),
        ("warm_max_steps", 10),
        ("warm_failure_rate", 9),
        ("cold_failure_rate", 9),
    ]
    header = " ".join(f"{label:>{w}}" if i > 0 else f"{label:<{w}}"
                      for i, (label, w) in enumerate(cols))
    print(header)
    print("-" * len(header))
    for row in rows:
        cells = []
        for i, (label, w) in enumerate(cols):
            v = row.get(label, "")
            if isinstance(v, float):
                cell = f"{v:.4g}"
            else:
                cell = str(v)
            cells.append(f"{cell:>{w}}" if i > 0 else f"{cell:<{w}}")
        flag_r = "R" if row.get("phase3_recovery_passed") else "."
        flag_i = "I" if row.get("phase3_improvement_passed") else "."
        print(" ".join(cells) + f"  [{flag_r}{flag_i}]")
    print("\nLegend: [R] recovery vs original baseline, [I] improvement vs post-charge-fix")


# ------------------------- driver -------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", default="tolerance_1d", choices=PLANS.keys())
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260504)
    parser.add_argument("--warm-solves", type=int, default=3)
    parser.add_argument("--dt", type=float, default=10.0)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 4))
    parser.add_argument("--out-dir", default="/tmp/ts1_sweeps")
    parser.add_argument("--summary", default=None,
                        help="Optional path to write summary JSON for all variants")
    args = parser.parse_args(argv)

    plan = PLANS[args.plan]()
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Plan '{args.plan}': {len(plan)} variants, {args.samples} samples each, "
          f"{args.workers} workers")
    print(f"Output dir: {args.out_dir}")

    # Workers must use the 'spawn' start method to avoid CUDA/OpenMP fork hazards.
    ctx = mp.get_context("spawn")

    pending = []
    t0_all = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx) as pool:
        futures = {}
        for v in plan:
            out_path = os.path.join(args.out_dir, f"{v['name']}.json")
            fut = pool.submit(
                _run_variant,
                v["name"], v["overrides"],
                args.samples, args.seed, args.warm_solves, args.dt,
                out_path,
            )
            futures[fut] = v["name"]

        for fut in as_completed(futures):
            name = futures[fut]
            try:
                result = fut.result()
            except Exception as exc:  # pragma: no cover - surface worker errors
                print(f"  [FAIL] {name}: {exc!r}")
                continue
            print(f"  [done in {result['wallclock_seconds']:6.1f}s] {name}: "
                  f"med={result['warm_median_steps']:.0f} p95={result['warm_p95_steps']:.1f} "
                  f"p99={result['warm_p99_steps']:.1f} maxsteps={result['warm_max_steps']:.0f} "
                  f"warm_fail={result['warm_failure_rate']:.2%}")
            pending.append(result)

    elapsed_all = time.perf_counter() - t0_all
    print(f"\nAll variants finished in {elapsed_all:.1f}s wallclock.\n")

    # Load reference snapshots once and run A/B summaries.
    with open(REPO_BASELINE) as f:
        pre_baseline = json.load(f)
    with open(POST_CHARGE_FIX_REF) as f:
        post_baseline = json.load(f)

    summaries = []
    for r in pending:
        s = summarize_run(r["name"], r["out_path"], pre_baseline, post_baseline)
        summaries.append(s)

    # Sort: passing variants first, then by warm_p95_steps ascending.
    def _sort_key(s):
        return (
            not s["phase3_improvement_passed"],
            not s["phase3_recovery_passed"],
            float(s.get("warm_p95_steps") or float("inf")),
        )

    summaries.sort(key=_sort_key)
    _print_table(summaries)

    if args.summary:
        with open(args.summary, "w") as f:
            json.dump(summaries, f, indent=2)
        print(f"\nSummary written to {args.summary}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
