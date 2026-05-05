# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# TS1-cloud Phase 5 acceptance gate: synthetic 3D column ensemble.
#
# The previous phases used Latin-Hypercube samples drawn from a 6-dim envelope
# of (T, P, LWC, SO2, H2O2, O3) — useful for stress-testing edge combinations
# but unweighted by realistic atmospheric structure. Phase 5 closes the loop:
# it builds an MPAS-A-style column ensemble where every cell carries a
# physically self-consistent T(z), P(z), gas profile (from the existing 13-
# level CSV) and only the cloudy levels (intersecting the prescribed-cloud
# pressure window) get a non-trivial LWC.
#
# Acceptance criteria (designed to fail loudly if a future change regresses):
#   - 0 substep non-convergences across the full ensemble
#   - Non-cloudy levels see no SO4 production (within rounding)
#   - Cloudy levels see strictly positive SO4 production
#   - Per-cell sulfur conservation within 1%
#
# We do NOT have actual MPAS-A 60km output checked into the repo, so this is
# a synthetic-but-realistic ensemble: the column comes straight from the MPAS
# `initial_conditions.csv` but horizontal variation is generated via small
# random perturbations on (T, LWC) representative of MPAS-A grid-cell spread.
# This keeps the gate runnable in CI while exercising the Phase-1..4 winning
# config across the full vertical envelope MPAS will hit.

import csv
import math
import os

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

REPO_ROOT = ts1_cloud._REPO_ROOT
MPAS_ICS_CSV = ts1_cloud.MPAS_TS1_CLOUD_ICS_CSV
PRESCRIBED_CLOUD_TXT = os.path.join(
    REPO_ROOT, "configs", "v1", "ts1_cloud", "prescribed_cloud.txt"
)


# US Standard Atmosphere (1976), simplified piecewise hydrostatic profile.
# Inputs in km; outputs (T [K], P [Pa]) good to within ~5% through 80 km — more
# than adequate for an acceptance gate that allows ±5K horizontal jitter.
def _us_std_atmos(height_km):
    z = height_km * 1000.0  # m
    T0 = 288.15
    P0 = 101325.0
    if z <= 11000.0:
        L = -0.0065  # K/m
        T = T0 + L * z
        P = P0 * (T / T0) ** (-9.80665 * 0.0289644 / (8.31446 * L))
    elif z <= 20000.0:
        T11 = 216.65
        P11 = 22632.06
        T = T11
        P = P11 * math.exp(-9.80665 * 0.0289644 * (z - 11000.0) / (8.31446 * T11))
    elif z <= 32000.0:
        T20, P20, L = 216.65, 5474.89, 0.001
        T = T20 + L * (z - 20000.0)
        P = P20 * (T / T20) ** (-9.80665 * 0.0289644 / (8.31446 * L))
    elif z <= 47000.0:
        T32, P32, L = 228.65, 868.02, 0.0028
        T = T32 + L * (z - 32000.0)
        P = P32 * (T / T32) ** (-9.80665 * 0.0289644 / (8.31446 * L))
    else:
        T = 270.65
        P = 110.91 * math.exp(-9.80665 * 0.0289644 * (z - 47000.0) / (8.31446 * T))
    return T, P


def _load_prescribed_cloud():
    """Read the (p_top, p_bot, lwc) triple from prescribed_cloud.txt."""
    with open(PRESCRIBED_CLOUD_TXT) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            return float(parts[0]), float(parts[1]), float(parts[2])
    raise RuntimeError("No usable line in prescribed_cloud.txt")


def _load_full_column_ics():
    """Return list of (height_km, gas_ics_dict) rows from the MPAS ICs CSV."""
    rows = []
    with open(MPAS_ICS_CSV) as f:
        reader = csv.reader(f)
        header = [c.strip() for c in next(reader)]
        for raw in reader:
            vals = [v.strip() for v in raw]
            height = float(vals[0])
            gas = {}
            for name, v in zip(header[1:], vals[1:]):
                try:
                    gas[name] = float(v)
                except ValueError:
                    pass
            rows.append((height, gas))
    return rows


def _build_3d_ensemble(n_horizontal, seed):
    """Build an N_horizontal × N_levels cell list for the acceptance run.

    For each (column, level):
      - Gas ICs come from the column's height row in the MPAS CSV.
      - T, P from US Standard Atmosphere with small horizontal jitter (±2K).
      - LWC: lognormal jitter (median = prescribed value, σ ≈ 0.3 in log10) for
        cloudy levels; 0 for clear-sky levels.
    A level is "cloudy" iff the standard-atmosphere pressure at that level
    falls within the prescribed_cloud pressure window.
    """
    p_top, p_bot, lwc_kg_kg = _load_prescribed_cloud()
    column_rows = _load_full_column_ics()
    rng = np.random.default_rng(seed)

    cells = []
    for col_id in range(n_horizontal):
        T_jitter_per_level = rng.normal(0.0, 2.0, size=len(column_rows))
        for lvl_idx, (height_km, gas_ics) in enumerate(column_rows):
            T_std, P_std = _us_std_atmos(height_km)
            T = T_std + float(T_jitter_per_level[lvl_idx])
            P = P_std
            cloudy = (P >= p_top) and (P <= p_bot)
            if cloudy:
                # MPAS LWC is given as kg(H2O)/kg(air); convert to kg/m³-air.
                rho_air = P / (287.05 * T)  # ideal gas, dry-air R
                lwc_median = lwc_kg_kg * rho_air
                lwc = float(lwc_median * np.exp(rng.normal(0.0, 0.3 * math.log(10.0))))
                lwc = max(lwc, 1.0e-5)  # don't go absurdly thin
            else:
                lwc = 0.0
            cells.append({
                "cell_id": (col_id, lvl_idx),
                "column_id": col_id,
                "level_index": lvl_idx,
                "height_km": height_km,
                "T": T,
                "P": P,
                "lwc": lwc,
                "cloudy": cloudy,
                "gas_ics": gas_ics,
            })
    return cells


def _build_cell_full_ics(cell):
    """Convert a 3D-ensemble cell into the MICM IC dict.

    Clear-sky (LWC=0) cells get a tiny placeholder LWC because the constraint
    solver still needs feasible aqueous concentrations to initialize. With
    LWC ~1e-12 kg/m³ the aqueous chemistry is effectively dormant.
    """
    lwc = cell["lwc"] if cell["cloudy"] else 1.0e-12
    ics = ts1_cloud._make_full_ics(lwc=lwc, T=cell["T"])
    # Override gas-phase species from this column's MPAS row.
    for name, val in cell["gas_ics"].items():
        ics[name] = val
    return ics


@pytest.mark.skipif(
    not os.path.exists(MPAS_ICS_CSV) or not os.path.exists(PRESCRIBED_CLOUD_TXT),
    reason="TS1-cloud reference column files not found",
)
class TestTS1Cloud3DAcceptance:
    """Phase 5 acceptance gate: synthetic 3D column ensemble."""

    def test_3d_ensemble_acceptance(self):
        """Run the MPAS-style 9 × 100s subcycled integration on every cell.

        Failure modes that this gate catches:
          • Any non-convergence across (n_horizontal × n_levels) cells.
          • Negative SO4 production in clear-sky cells (would imply the
            constraint solver is leaking S into aqueous phase even at LWC=0).
          • Non-positive SO4 production in cloudy cells.
          • >1% sulfur drift in any cell.
        """
        n_horizontal = int(os.getenv("TS1_CLOUD_3D_COLUMNS", "8"))
        seed = int(os.getenv("TS1_CLOUD_3D_SEED", "20260504"))

        cells = _build_3d_ensemble(n_horizontal, seed)
        # One MICM reused across all cells (mirrors how MPAS would run it).
        micm = ts1_cloud._create_micm()

        n_cloudy = sum(1 for c in cells if c["cloudy"])
        n_clear = len(cells) - n_cloudy
        print(
            f"\n  3D ensemble: {n_horizontal} columns × {len(cells)//n_horizontal} levels "
            f"= {len(cells)} cells ({n_cloudy} cloudy / {n_clear} clear)"
        )

        non_converged = []
        sulfur_drifts = []
        cloudy_so4_deltas = []
        clear_so4_deltas = []

        for cell in cells:
            ics = _build_cell_full_ics(cell)
            state = ts1_cloud._set_state(micm, ics, cell["T"], cell["P"])
            so4_init = ts1_cloud._get(state, "CLOUD.AQUEOUS.SO4mm")
            s_init = ts1_cloud._sulfur_total(state)

            ok = True
            for substep in range(ts1_cloud.N_SUBSTEPS):
                r = micm.solve(state, time_step=ts1_cloud.DT_SUBSTEP)
                if r.state != SolverState.Converged:
                    non_converged.append((cell["cell_id"], substep, str(r.state)))
                    ok = False
                    break
            if not ok:
                continue

            so4_final = ts1_cloud._get(state, "CLOUD.AQUEOUS.SO4mm")
            s_final = ts1_cloud._sulfur_total(state)
            drift = abs(s_final - s_init) / max(abs(s_init), 1e-30)
            sulfur_drifts.append(drift)
            delta = so4_final - so4_init
            (cloudy_so4_deltas if cell["cloudy"] else clear_so4_deltas).append(delta)

        max_drift = max(sulfur_drifts) if sulfur_drifts else 0.0
        print(
            f"  non_converged={len(non_converged)}  max_S_drift={max_drift:.2e}  "
            f"cloudy_dSO4 min={min(cloudy_so4_deltas, default=0):.2e} "
            f"max={max(cloudy_so4_deltas, default=0):.2e}  "
            f"clear_dSO4 max_abs={max((abs(d) for d in clear_so4_deltas), default=0):.2e}"
        )
        for entry in non_converged[:5]:
            print(f"    NC: {entry}")

        assert not non_converged, (
            f"{len(non_converged)} substeps failed across {len(cells)} cells; "
            f"first: {non_converged[:3]}"
        )
        assert max_drift < 0.01, (
            f"max sulfur drift {max_drift:.2%} exceeds 1% in some cell"
        )
        assert all(d > 0 for d in cloudy_so4_deltas), (
            "Every cloudy cell should produce SO4; some did not"
        )
        # Clear cells: aqueous chemistry is dormant (LWC=1e-12), so any SO4
        # change must be at the floating-point floor.
        assert all(abs(d) < 1e-12 for d in clear_so4_deltas), (
            f"Clear-sky cells produced non-trivial SO4 (max |Δ| = "
            f"{max((abs(d) for d in clear_so4_deltas), default=0):.2e})"
        )
