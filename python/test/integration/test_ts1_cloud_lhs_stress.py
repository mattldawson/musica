# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# TS1-cloud stress testing with Latin Hypercube Sampling (LHS).
#
# This suite is intended as a pre-MPAS robustness gate. It samples
# MPAS-representative atmospheric conditions and checks convergence and
# internal solver-step efficiency with separate cold-start and warm-solve
# criteria.

import math
import os

import numpy as np
import pytest

try:
    import test_ts1_cloud_box_model as ts1_cloud
except ModuleNotFoundError:  # pragma: no cover
    from python.test.integration import test_ts1_cloud_box_model as ts1_cloud

from musica.micm import SolverState


DEFAULT_QUICK_SAMPLES = 20
DEFAULT_WARM_SOLVES = 3
DEFAULT_DT = 10.0

# Solver KPI targets requested for warm solves.
WARM_MEDIAN_STEPS_TARGET = 40
WARM_P95_STEPS_TARGET = 80

# Fraction caps used with small-sample guardrails.
WARM_NONCONVERGED_RATE = 0.10
COLD_NONCONVERGED_RATE = 0.25
WARM_HIGH_STEP_RATE = 0.15


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


def _run_stress(n_samples, seed, warm_solves, dt):
    micm = ts1_cloud._create_micm()
    samples = _sample_conditions(n_samples=n_samples, seed=seed)

    cold_attempts = 0
    cold_failures = 0
    cold_steps = []

    warm_attempts = 0
    warm_failures = 0
    warm_steps = []
    warm_kpi_steps = []

    outlier_samples = []

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

        cold_attempts += 1
        cold_result = micm.solve(state, time_step=dt)
        cold_steps.append(int(cold_result.stats.number_of_steps))

        if cold_result.state != SolverState.Converged:
            cold_failures += 1
            outlier_samples.append((sample["sample_id"], "cold", str(cold_result.state), int(cold_result.stats.number_of_steps)))
            continue

        for warm_idx in range(warm_solves):
            warm_attempts += 1
            warm_result = micm.solve(state, time_step=dt)
            if warm_result.state != SolverState.Converged:
                warm_failures += 1
                outlier_samples.append((sample["sample_id"], "warm", str(warm_result.state), int(warm_result.stats.number_of_steps)))
                break
            steps = int(warm_result.stats.number_of_steps)
            warm_steps.append(steps)
            if warm_idx >= 1:
                warm_kpi_steps.append(steps)

    warm_high_steps = sum(1 for s in warm_kpi_steps if s > WARM_P95_STEPS_TARGET)

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
        "warm_high_steps": warm_high_steps,
        "warm_high_step_rate": (warm_high_steps / len(warm_kpi_steps)) if warm_kpi_steps else 0.0,
        "outliers": outlier_samples[:10],
    }
    return summary


def _assert_summary(summary):
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
    if summary["warm_kpi_count"] >= 100:
        assert summary["warm_p95_steps"] < WARM_P95_STEPS_TARGET, (
            f"Warm P95 internal steps too high: {summary['warm_p95_steps']:.1f} "
            f">= {WARM_P95_STEPS_TARGET}"
        )

    warm_high_cap = _max_allowed(
        WARM_HIGH_STEP_RATE,
        summary["warm_kpi_count"],
        min_count=2,
    )
    assert summary["warm_high_steps"] <= warm_high_cap, (
        f"Too many warm high-step outliers (> {WARM_P95_STEPS_TARGET}): "
        f"{summary['warm_high_steps']} > {warm_high_cap}"
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

    _assert_summary(summary)
