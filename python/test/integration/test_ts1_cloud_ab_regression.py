# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# Regression gate that compares a candidate TS1-cloud LHS stress metrics file
# (produced via TS1_CLOUD_STRESS_METRICS_OUT) against the committed baseline.
#
# The candidate file is selected via the TS1_CLOUD_STRESS_CANDIDATE env var.
# When the env var is unset the test passes trivially, so the gate is a
# no-op in routine CI. To run the gate locally:
#
#   1. Generate a candidate run, e.g.
#        TS1_CLOUD_STRESS_STRICT_SAMPLES=500 \
#          TS1_CLOUD_STRESS_METRICS_OUT=/tmp/candidate.json \
#          uv run pytest -v -s \
#            python/test/integration/test_ts1_cloud_lhs_stress.py::test_ts1_cloud_lhs_stress_strict_deep_gate
#   2. Run the regression gate
#        TS1_CLOUD_STRESS_CANDIDATE=/tmp/candidate.json \
#          uv run pytest -v -s python/test/integration/test_ts1_cloud_ab_regression.py

import json
import os

import pytest

try:
    from . import ab_compare
except ImportError:  # pragma: no cover
    try:
        import ab_compare
    except ModuleNotFoundError:
        from python.test.integration import ab_compare


REPO_BASELINE = os.path.join(
    os.path.dirname(__file__), "data", "ts1_cloud_baseline.json"
)


def test_baseline_file_exists():
    """The committed baseline must be present and parseable."""
    assert os.path.exists(REPO_BASELINE), (
        f"Committed baseline missing: {REPO_BASELINE}"
    )
    with open(REPO_BASELINE) as f:
        payload = json.load(f)
    for required in (
        "samples",
        "warm_median_steps",
        "warm_p95_steps",
        "warm_step_histogram",
        "warm_kpi_aggregates",
    ):
        assert required in payload, f"Baseline missing key: {required}"


def test_ab_compare_self_baseline_is_zero_delta():
    """Comparing the baseline against itself must yield zero deltas."""
    with open(REPO_BASELINE) as f:
        payload = json.load(f)
    deltas = ab_compare.compute_delta(payload, payload)
    for key, rec in deltas["metrics"].items():
        d = rec["delta"]
        if isinstance(d, float) and d != d:  # NaN check without importing math
            continue
        assert abs(d) < 1.0e-9, f"Self-comparison should be zero, got {key}={d}"
    for key, rec in deltas["histogram"].items():
        assert rec["delta"] == 0, f"Self-histogram should be zero, got {key}={rec['delta']}"
    regressed, msgs = ab_compare.find_regressions(deltas)
    assert not regressed, f"Self-comparison flagged regressions: {msgs}"


def test_candidate_against_baseline():
    """If a candidate metrics file is provided, gate it against the baseline."""
    candidate_path = os.getenv("TS1_CLOUD_STRESS_CANDIDATE")
    if not candidate_path:
        print("\nA/B regression gate disabled (set TS1_CLOUD_STRESS_CANDIDATE to a metrics JSON path).")
        return

    assert os.path.exists(candidate_path), f"Candidate metrics file missing: {candidate_path}"

    with open(REPO_BASELINE) as f:
        baseline = json.load(f)
    with open(candidate_path) as f:
        candidate = json.load(f)

    deltas = ab_compare.compute_delta(baseline, candidate)
    print("\n" + ab_compare.format_report(deltas))

    regressed, messages = ab_compare.find_regressions(deltas)
    if regressed:
        pytest.fail("Candidate regresses vs committed baseline:\n  - " + "\n  - ".join(messages))
