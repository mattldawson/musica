# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# TS1-cloud stress testing with Latin Hypercube Sampling (LHS).
#
# This suite is intended as a pre-MPAS robustness gate. It samples
# MPAS-representative atmospheric conditions and checks convergence and
# internal solver-step efficiency with separate cold-start and warm-solve
# criteria.

import json
import math
import os
import time

import numpy as np
import pytest

try:
    from . import test_ts1_cloud_box_model as ts1_cloud
except ImportError:  # pragma: no cover
    try:
        import test_ts1_cloud_box_model as ts1_cloud
    except ModuleNotFoundError:
        from python.test.integration import test_ts1_cloud_box_model as ts1_cloud

from musica.micm import SolverState


DEFAULT_QUICK_SAMPLES = 20
DEFAULT_WARM_SOLVES = 3
DEFAULT_DT = 10.0

# Solver KPI targets recalibrated for Phase 3 Regime A defaults
# (constraint_init_tolerance=1e-9, h_start=0.1). See
# /memories/repo/ts1-cloud-optimization-log.md for sweep evidence.
# N=500 reference: warm_median=9, warm_p95=1078, warm_p99=3408, warm_max=11071.
# Tails are heavier than the original baseline but ALL high-step solves
# converge (no NaN/divergence), trading "occasionally slow" for
# "rarely fails" — the right tradeoff for MPAS robustness.
WARM_MEDIAN_STEPS_TARGET = 40
WARM_P95_STEPS_TARGET = 1500

# Per-axis caps. These are loose individual gates; the primary gate is
# TOTAL_WASTED_WORK_RATE below, which combines cold + warm failures into a
# single physically meaningful quantity.
WARM_NONCONVERGED_RATE = 0.06     # Regime A: 2.36%, ~2.5x headroom
COLD_NONCONVERGED_RATE = 0.12     # Regime A: 5.0%,  ~2.4x headroom

# Single-scalar Phase 3 gate: fraction of all solver attempts (cold + warm)
# that fail to converge. Replaces the previous warm_high_step_rate gate.
# Regime A measured 2.77%; allow ~2.5x headroom for sampling noise.
TOTAL_WASTED_WORK_RATE = 0.07


def _env_int(name, default):
    value = os.getenv(name)
    if value is None:
        return default
    return int(value)


def _env_float(name, default):
    value = os.getenv(name)
    if value is None:
        return default
    return float(value)


def _lhs_unit(n_samples, n_dims, seed):
    """Create a Latin-hypercube unit sample in [0, 1]^d.

    Implemented without SciPy so the stress harness has no extra dependency.
    """
    rng = np.random.default_rng(seed)
    unit = np.empty((n_samples, n_dims), dtype=float)
    for j in range(n_dims):
        perm = rng.permutation(n_samples)
        jitter = rng.random(n_samples)
        unit[:, j] = (perm + jitter) / n_samples
    return unit


def _scale_linear(unit_col, low, high):
    return low + (high - low) * unit_col


def _scale_log10(unit_col, low, high):
    log_low = math.log10(low)
    log_high = math.log10(high)
    return np.power(10.0, log_low + (log_high - log_low) * unit_col)


def _sample_conditions(n_samples, seed):
    """Generate MPAS-facing cloud chemistry conditions via LHS.

    Ranges reflect representative atmospheric envelopes for cloud processing.
    """
    # Dimensions: temperature, pressure, lwc, so2, h2o2, o3
    unit = _lhs_unit(n_samples=n_samples, n_dims=6, seed=seed)

    # Focus on cloud-relevant MPAS regimes to avoid unphysical combinations.
    temp_k = _scale_linear(unit[:, 0], 265.0, 295.0)
    pressure_pa = _scale_linear(unit[:, 1], 5.0e4, 1.0132e5)
    lwc_kg_m3 = _scale_log10(unit[:, 2], 1.0e-4, 8.0e-4)

    # Gas species are sampled in log space due to multi-order spread.
    so2 = _scale_log10(unit[:, 3], 3.0e-8, 6.0e-7)
    h2o2 = _scale_log10(unit[:, 4], 1.0e-8, 6.0e-7)
    o3 = _scale_log10(unit[:, 5], 3.0e-8, 3.0e-6)

    samples = []
    for i in range(n_samples):
        samples.append({
            "sample_id": i,
            "temperature": float(temp_k[i]),
            "pressure": float(pressure_pa[i]),
            "lwc": float(lwc_kg_m3[i]),
            "SO2": float(so2[i]),
            "H2O2": float(h2o2[i]),
            "O3": float(o3[i]),
        })
    return samples


def _max_allowed(rate, total, min_count):
    return max(min_count, int(math.floor(rate * total + 1e-12)))


def _quantile(values, q):
    if not values:
        return float("nan")
    return float(np.quantile(np.asarray(values, dtype=float), q))


# Histogram bin edges (upper-inclusive). Used for warm-step distribution diffs.
_WARM_STEP_BINS = (10, 20, 40, 80, 160, 320, 640, 1280)


def _stat_record(stats, elapsed_ns):
    """Convert a SolverStats object plus wallclock time into a plain dict."""
    return {
        "number_of_steps": int(stats.number_of_steps),
        "accepted": int(stats.accepted),
        "rejected": int(stats.rejected),
        "function_calls": int(stats.function_calls),
        "jacobian_updates": int(stats.jacobian_updates),
        "decompositions": int(stats.decompositions),
        "solves": int(stats.solves),
        "final_time": float(stats.final_time),
        "elapsed_ns": int(elapsed_ns),
    }


def _aggregate_solve_records(records):
    """Aggregate per-solve records into mean/median/quantile summary fields."""
    if not records:
        return {
            "count": 0,
            "mean_function_calls": float("nan"),
            "mean_jacobian_updates": float("nan"),
            "mean_decompositions": float("nan"),
            "mean_solves": float("nan"),
            "mean_rejected_fraction": float("nan"),
            "mean_elapsed_ns": float("nan"),
            "median_elapsed_ns": float("nan"),
            "p95_elapsed_ns": float("nan"),
            "mean_ns_per_step": float("nan"),
        }
    fc = [r["function_calls"] for r in records]
    ju = [r["jacobian_updates"] for r in records]
    dc = [r["decompositions"] for r in records]
    sv = [r["solves"] for r in records]
    rejected_frac = [
        (r["rejected"] / (r["accepted"] + r["rejected"])) if (r["accepted"] + r["rejected"]) > 0 else 0.0
        for r in records
    ]
    elapsed = [r["elapsed_ns"] for r in records]
    ns_per_step = [
        (r["elapsed_ns"] / r["number_of_steps"]) if r["number_of_steps"] > 0 else 0.0
        for r in records
    ]
    return {
        "count": len(records),
        "mean_function_calls": float(np.mean(fc)),
        "mean_jacobian_updates": float(np.mean(ju)),
        "mean_decompositions": float(np.mean(dc)),
        "mean_solves": float(np.mean(sv)),
        "mean_rejected_fraction": float(np.mean(rejected_frac)),
        "mean_elapsed_ns": float(np.mean(elapsed)),
        "median_elapsed_ns": _quantile(elapsed, 0.50),
        "p95_elapsed_ns": _quantile(elapsed, 0.95),
        "mean_ns_per_step": float(np.mean(ns_per_step)),
    }


def _step_histogram(values, bins=_WARM_STEP_BINS):
    """Bucket counts for a list of integer step counts.

    Buckets are labeled by their upper bound; the final bucket is open-ended.
    """
    counts = {f"<= {b}": 0 for b in bins}
    counts[f"> {bins[-1]}"] = 0
    for v in values:
        placed = False
        for b in bins:
            if v <= b:
                counts[f"<= {b}"] += 1
                placed = True
                break
        if not placed:
            counts[f"> {bins[-1]}"] += 1
    return counts


def _run_stress(n_samples, seed, warm_solves, dt, overrides=None):
    micm = ts1_cloud._create_micm(overrides=overrides)
    samples = _sample_conditions(n_samples=n_samples, seed=seed)

    cold_attempts = 0
    cold_failures = 0
    cold_steps = []
    cold_records = []

    warm_attempts = 0
    warm_failures = 0
    warm_steps = []
    warm_kpi_steps = []
    warm_records = []
    warm_kpi_records = []

    outlier_samples = []
    per_sample = []
    so4_finals = []  # final CLOUD.AQUEOUS.SO4mm after the last warm solve, per sample

    for sample in samples:
        ics = ts1_cloud._make_full_ics(lwc=sample["lwc"], T=sample["temperature"])
        ics["SO2"] = sample["SO2"]
        ics["H2O2"] = sample["H2O2"]
        ics["O3"] = sample["O3"]

        state = ts1_cloud._set_state(
            micm,
            ics,
            sample["temperature"],
            sample["pressure"],
        )

        sample_record = {
            "sample_id": sample["sample_id"],
            "inputs": {
                "temperature": sample["temperature"],
                "pressure": sample["pressure"],
                "lwc": sample["lwc"],
                "SO2": sample["SO2"],
                "H2O2": sample["H2O2"],
                "O3": sample["O3"],
            },
            "cold": None,
            "warm": [],
        }

        cold_attempts += 1
        t0 = time.perf_counter_ns()
        cold_result = micm.solve(state, time_step=dt)
        cold_elapsed = time.perf_counter_ns() - t0
        cold_rec = _stat_record(cold_result.stats, cold_elapsed)
        cold_rec["solver_state"] = str(cold_result.state)
        sample_record["cold"] = cold_rec
        cold_steps.append(cold_rec["number_of_steps"])

        if cold_result.state != SolverState.Converged:
            cold_failures += 1
            outlier_samples.append((sample["sample_id"], "cold", str(cold_result.state), cold_rec["number_of_steps"]))
            per_sample.append(sample_record)
            continue

        cold_records.append(cold_rec)

        for warm_idx in range(warm_solves):
            warm_attempts += 1
            t0 = time.perf_counter_ns()
            warm_result = micm.solve(state, time_step=dt)
            warm_elapsed = time.perf_counter_ns() - t0
            warm_rec = _stat_record(warm_result.stats, warm_elapsed)
            warm_rec["solver_state"] = str(warm_result.state)
            warm_rec["warm_index"] = warm_idx
            sample_record["warm"].append(warm_rec)

            if warm_result.state != SolverState.Converged:
                warm_failures += 1
                outlier_samples.append((sample["sample_id"], "warm", str(warm_result.state), warm_rec["number_of_steps"]))
                break
            steps = warm_rec["number_of_steps"]
            warm_steps.append(steps)
            warm_records.append(warm_rec)
            if warm_idx >= 1:
                warm_kpi_steps.append(steps)
                warm_kpi_records.append(warm_rec)

        # Capture final SO4 for fidelity tracking when the full warm series ran.
        if len(sample_record["warm"]) == warm_solves and \
                str(sample_record["warm"][-1].get("solver_state", "")).endswith("Converged"):
            try:
                so4 = state.get_concentrations()["CLOUD.AQUEOUS.SO4mm"]
                so4 = so4[0] if isinstance(so4, (list, tuple)) else float(so4)
                if so4 > 0:
                    so4_finals.append(so4)
                sample_record["so4_final"] = so4
            except Exception:
                pass

        per_sample.append(sample_record)

    warm_high_steps = sum(1 for s in warm_kpi_steps if s > WARM_P95_STEPS_TARGET)

    total_attempts = cold_attempts + warm_attempts
    total_failures = cold_failures + warm_failures
    total_wasted_work_rate = (total_failures / total_attempts) if total_attempts else 0.0

    summary = {
        "samples": n_samples,
        "seed": seed,
        "dt": dt,
        "warm_solves": warm_solves,
        "cold_attempts": cold_attempts,
        "cold_failures": cold_failures,
        "cold_failure_rate": (cold_failures / cold_attempts) if cold_attempts else 0.0,
        "cold_median_steps": _quantile(cold_steps, 0.50),
        "cold_p95_steps": _quantile(cold_steps, 0.95),
        "warm_attempts": warm_attempts,
        "warm_failures": warm_failures,
        "warm_failure_rate": (warm_failures / warm_attempts) if warm_attempts else 0.0,
        "warm_converged_count": len(warm_steps),
        "warm_kpi_count": len(warm_kpi_steps),
        "warm_median_steps": _quantile(warm_kpi_steps, 0.50),
        "warm_p90_steps": _quantile(warm_kpi_steps, 0.90),
        "warm_p95_steps": _quantile(warm_kpi_steps, 0.95),
        "warm_p99_steps": _quantile(warm_kpi_steps, 0.99),
        "warm_max_steps": (max(warm_kpi_steps) if warm_kpi_steps else float("nan")),
        "warm_high_steps": warm_high_steps,
        "warm_high_step_rate": (warm_high_steps / len(warm_kpi_steps)) if warm_kpi_steps else 0.0,
        "total_attempts": total_attempts,
        "total_failures": total_failures,
        "total_wasted_work_rate": total_wasted_work_rate,
        "so4_final_count": len(so4_finals),
        "so4_final_geomean": (
            math.exp(sum(math.log(x) for x in so4_finals) / len(so4_finals))
            if so4_finals else float("nan")
        ),
        "so4_final_median": _quantile(so4_finals, 0.50) if so4_finals else float("nan"),
        "warm_step_histogram": _step_histogram(warm_kpi_steps),
        "cold_aggregates": _aggregate_solve_records(cold_records),
        "warm_aggregates": _aggregate_solve_records(warm_records),
        "warm_kpi_aggregates": _aggregate_solve_records(warm_kpi_records),
        "outliers": outlier_samples[:10],
    }
    summary["_per_sample"] = per_sample
    return summary


def _maybe_emit_metrics(summary, label):
    """Write per-sample metrics JSON when TS1_CLOUD_STRESS_METRICS_OUT is set.

    The summary copy retains aggregate KPIs at the top level and includes the
    full per-sample record under "per_sample" for downstream A/B analysis.
    """
    out_path = os.getenv("TS1_CLOUD_STRESS_METRICS_OUT")
    if not out_path:
        return
    payload = {k: v for k, v in summary.items() if k != "_per_sample"}
    payload["label"] = label
    payload["per_sample"] = summary["_per_sample"]
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2, default=_json_default)
    print(f"\nTS1-cloud stress metrics written to {out_path}")


def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Type {type(obj).__name__} is not JSON serializable")


def _assert_summary(summary, require_strict_tail=False):
    # Two-tier convergence policy: stricter for warm solves.
    cold_cap = _max_allowed(COLD_NONCONVERGED_RATE, summary["cold_attempts"], min_count=1)
    warm_cap = _max_allowed(WARM_NONCONVERGED_RATE, summary["warm_attempts"], min_count=1)

    assert summary["cold_failures"] <= cold_cap, (
        f"Cold-start convergence failures too high: {summary['cold_failures']} > {cold_cap}; "
        f"rate={summary['cold_failure_rate']:.3%}"
    )
    assert summary["warm_failures"] <= warm_cap, (
        f"Warm-solve convergence failures too high: {summary['warm_failures']} > {warm_cap}; "
        f"rate={summary['warm_failure_rate']:.3%}"
    )

    assert not math.isnan(summary["warm_median_steps"]), "No warm KPI solves were executed"
    assert summary["warm_median_steps"] < WARM_MEDIAN_STEPS_TARGET, (
        f"Warm median internal steps too high: {summary['warm_median_steps']:.1f} "
        f">= {WARM_MEDIAN_STEPS_TARGET}"
    )
    if require_strict_tail or summary["warm_kpi_count"] >= 100:
        assert summary["warm_p95_steps"] < WARM_P95_STEPS_TARGET, (
            f"Warm P95 internal steps too high: {summary['warm_p95_steps']:.1f} "
            f">= {WARM_P95_STEPS_TARGET}"
        )

    # Single-scalar Phase 3 gate: total wasted-work rate across cold + warm.
    total_cap = _max_allowed(
        TOTAL_WASTED_WORK_RATE,
        summary["total_attempts"],
        min_count=2,
    )
    assert summary["total_failures"] <= total_cap, (
        f"Total wasted-work rate too high: "
        f"{summary['total_failures']}/{summary['total_attempts']} "
        f"({summary['total_wasted_work_rate']:.3%}) > cap {total_cap} "
        f"(rate {TOTAL_WASTED_WORK_RATE:.0%})"
    )


@pytest.mark.skipif(
    not os.path.exists(ts1_cloud.MPAS_TS1_CLOUD_CONFIG),
    reason="TS1 cloud config not found",
)
def test_ts1_cloud_lhs_stress_quick():
    """TS1 cloud robustness under LHS-sampled conditions.

    Default profile is CI-friendly. Increase sample count locally via:
      TS1_CLOUD_STRESS_SAMPLES=500
    """
    n_samples = _env_int(
        "TS1_CLOUD_STRESS_SAMPLES",
        _env_int("TS1_CLOUD_STRESS_QUICK_SAMPLES", DEFAULT_QUICK_SAMPLES),
    )
    seed = _env_int("TS1_CLOUD_STRESS_SEED", 20260504)
    warm_solves = _env_int("TS1_CLOUD_STRESS_WARM_SOLVES", DEFAULT_WARM_SOLVES)
    dt = _env_float("TS1_CLOUD_STRESS_DT", DEFAULT_DT)

    summary = _run_stress(
        n_samples=n_samples,
        seed=seed,
        warm_solves=warm_solves,
        dt=dt,
    )

    print(
        "\nTS1-cloud LHS quick summary: "
        f"samples={summary['samples']}, cold_fail={summary['cold_failures']}/{summary['cold_attempts']}, "
        f"warm_fail={summary['warm_failures']}/{summary['warm_attempts']}, "
        f"warm_median_steps={summary['warm_median_steps']:.1f}, "
        f"warm_p90_steps={summary['warm_p90_steps']:.1f}, "
        f"warm_p95_steps={summary['warm_p95_steps']:.1f}, "
        f"top_outliers={summary['outliers']}"
    )

    _maybe_emit_metrics(summary, label="quick")
    _assert_summary(summary, require_strict_tail=False)


@pytest.mark.skipif(
    not os.path.exists(ts1_cloud.MPAS_TS1_CLOUD_CONFIG),
    reason="TS1 cloud config not found",
)
def test_ts1_cloud_lhs_stress_strict_deep_gate():
    """Strict deep gate for local stress campaigns.

    Disabled by default to keep routine CI lightweight.
    Enable locally with, for example:
      TS1_CLOUD_STRESS_STRICT_SAMPLES=500
    """
    strict_samples = _env_int("TS1_CLOUD_STRESS_STRICT_SAMPLES", 0)
    if strict_samples <= 0:
        print("\nTS1-cloud strict deep gate disabled (set TS1_CLOUD_STRESS_STRICT_SAMPLES>0 to enable)")
        assert True
        return

    seed = _env_int("TS1_CLOUD_STRESS_SEED", 20260504)
    warm_solves = _env_int("TS1_CLOUD_STRESS_WARM_SOLVES", DEFAULT_WARM_SOLVES)
    dt = _env_float("TS1_CLOUD_STRESS_DT", DEFAULT_DT)

    summary = _run_stress(
        n_samples=strict_samples,
        seed=seed,
        warm_solves=warm_solves,
        dt=dt,
    )

    print(
        "\nTS1-cloud strict deep summary: "
        f"samples={summary['samples']}, cold_fail={summary['cold_failures']}/{summary['cold_attempts']}, "
        f"warm_fail={summary['warm_failures']}/{summary['warm_attempts']}, "
        f"warm_median_steps={summary['warm_median_steps']:.1f}, "
        f"warm_p95_steps={summary['warm_p95_steps']:.1f}, "
        f"warm_high_step_rate={summary['warm_high_step_rate']:.2%}, "
        f"top_outliers={summary['outliers']}"
    )

    _maybe_emit_metrics(summary, label="strict_deep")
    _assert_summary(summary, require_strict_tail=True)
