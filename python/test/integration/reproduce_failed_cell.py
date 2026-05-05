#!/usr/bin/env python
# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# CheMPAS-A ts1_cloud blowup reproduction.
#
# Loads a per-cell pre/post-solve dump produced by `diag_dump_failing_cells`
# in `mpas_chemistry_driver.F90` and replays it through the Python MICM API
# using the same MPAS config and the same tuned RosenbrockSolverParameters.
#
# Goal: confirm whether the explosion is reproducible outside MPAS, and if so,
# bisect the cause (single-cell vs batched, tuned vs default tolerances,
# Phase 2a aqueous defaults vs not).
#
# Usage:
#     python python/test/integration/reproduce_failed_cell.py \
#         --csv MPAS-Model/data/jw_480km_ts1_cloud/diag/failed_cell_00044_k03.csv \
#         [--config configs/v1/ts1_cloud/config.json] \
#         [--dt 100.0] [--n-substeps 9]
#
# Default config matches the box model (`configs/v1/ts1_cloud/config.json`).
# The MPAS deployed copy lives at `MPAS-Model/chemistry_data/ts1_cloud/config.json`
# and should be byte-identical; pass --config to use the MPAS copy explicitly.

import argparse
import csv
import math
import os
import sys

from musica.micm import MICM, SolverType, SolverState, RosenbrockSolverParameters


_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
DEFAULT_CONFIG = os.path.join(_REPO_ROOT, "configs", "v1", "ts1_cloud", "config.json")
MPAS_CONFIG = os.path.join(
    _REPO_ROOT, "MPAS-Model", "chemistry_data", "ts1_cloud", "config.json"
)


def _parse_dump(path):
    """Parse a diag_dump_failing_cells CSV.

    Returns (header: dict, species: dict[name -> (pre, post)], rates: dict[name -> (pre, post)]).
    """
    header = {}
    species = {}
    rates = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                # "# key=value"
                kv = line.lstrip("# ").split("=", 1)
                if len(kv) == 2:
                    header[kv[0].strip()] = kv[1].strip()
                continue
            if line.startswith("kind,"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 4:
                continue
            kind, name, pre, post = parts[0], parts[1], parts[2], parts[3]
            try:
                pre_f = float(pre)
            except ValueError:
                continue
            try:
                post_f = float(post)
            except ValueError:
                post_f = float("nan")  # rate_params have post="NA"
            if kind == "species":
                species[name] = (pre_f, post_f)
            elif kind == "rate_param":
                rates[name] = (pre_f, post_f)
    return header, species, rates


def _build_mpas_tuned_micm(config_path):
    """Recreate the exact MPAS-tuned solver setup."""
    micm = MICM(
        config_path=config_path,
        solver_type=SolverType.rosenbrock_dae4_standard_order,
    )
    tmp_state = micm.create_state()
    ordering = tmp_state.get_species_ordering()

    atol_aqueous = 1e-9
    atol_gas_alg = 1e-9
    atol_gas_diff = 1e-3

    abs_tols = [atol_gas_diff] * len(ordering)
    for name, idx in ordering.items():
        if "CLOUD.AQUEOUS." in name:
            abs_tols[idx] = atol_aqueous
        elif name in ("SO2", "H2O2", "O3"):
            abs_tols[idx] = atol_gas_alg

    micm.set_solver_parameters(RosenbrockSolverParameters(
        absolute_tolerances=abs_tols,
        h_start=0.1,
        constraint_init_max_iterations=100,
        constraint_init_tolerance=1e-9,
        max_number_of_steps=200000,
    ))
    return micm, ordering


def _build_default_micm(config_path):
    """Default (untuned) MICM — for comparison."""
    micm = MICM(
        config_path=config_path,
        solver_type=SolverType.rosenbrock_dae4_standard_order,
    )
    return micm, micm.create_state().get_species_ordering()


def _set_state_from_dump(micm, species_pre, rates_pre, header):
    state = micm.create_state()
    state.set_conditions(
        temperatures=float(header["T_K"]),
        pressures=float(header["P_Pa"]),
        air_densities=float(header["air_mol_per_m3"]),
    )
    sp_ord = state.get_species_ordering()
    rp_ord = state.get_user_defined_rate_parameters_ordering()
    sp_in = {n: v for n, (v, _) in species_pre.items() if n in sp_ord}
    rp_in = {n: v for n, (v, _) in rates_pre.items() if n in rp_ord}
    state.set_concentrations(sp_in)
    state.set_user_defined_rate_parameters(rp_in)
    return state, sp_ord, rp_ord


def _summarise(state, sp_ord, label, top_n=12):
    concs = state.get_concentrations()
    rows = []
    for name, idx in sp_ord.items():
        v = concs[name]
        v = v[0] if isinstance(v, list) else float(v)
        rows.append((name, v))
    # Detect blowup: NaN, Inf, or |v| > 1e6 (clearly nonphysical)
    bad = [(n, v) for (n, v) in rows
           if (not math.isfinite(v)) or abs(v) > 1.0e6]
    rows.sort(key=lambda r: abs(r[1]), reverse=True)
    print(f"\n[{label}] top {top_n} by |value|:")
    for n, v in rows[:top_n]:
        flag = "  <BAD>" if (not math.isfinite(v) or abs(v) > 1.0e6) else ""
        print(f"  {n:35s} = {v: .6e}{flag}")
    print(f"[{label}] # nonfinite/huge species = {len(bad)}")
    return bad


def _run_one(label, micm, header, species_pre, rates_pre, dt, n_substeps):
    state, sp_ord, rp_ord = _set_state_from_dump(micm, species_pre, rates_pre, header)
    print(f"\n=== {label} ===")
    print(f"  species set: {sum(1 for n in species_pre if n in sp_ord)} / {len(species_pre)} (ordering size {len(sp_ord)})")
    print(f"  rates   set: {sum(1 for n in rates_pre  if n in rp_ord)} / {len(rates_pre)} (ordering size {len(rp_ord)})")
    _summarise(state, sp_ord, f"{label} pre-solve", top_n=8)
    for step in range(n_substeps):
        result = micm.solve(state, time_step=dt)
        bad_now = _summarise(state, sp_ord,
                             f"{label} after substep {step+1} (state={result.state})",
                             top_n=8)
        if result.state != SolverState.Converged or bad_now:
            print(f"  ** stopped after substep {step+1}: solver_state={result.state}, bad={len(bad_now)}")
            return state, result, bad_now
    return state, result, []


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", required=True, help="failed_cell_*.csv produced by MPAS")
    p.add_argument("--config", default=DEFAULT_CONFIG,
                   help=f"MICM config.json (default: {DEFAULT_CONFIG})")
    p.add_argument("--mpas-config", action="store_true",
                   help=f"shortcut: use {MPAS_CONFIG}")
    p.add_argument("--dt", type=float, default=100.0,
                   help="substep dt seconds (default 100, MPAS subcycled)")
    p.add_argument("--n-substeps", type=int, default=9,
                   help="number of substeps to run (default 9 = MPAS 900s)")
    p.add_argument("--mode", choices=["tuned", "default", "both"], default="both")
    p.add_argument("--zero-aqueous", action="store_true",
                   help="Zero out all CLOUD.AQUEOUS.* species (except H2O) "
                        "before solving — tests whether Phase 2a defaults "
                        "are the trigger.")
    p.add_argument("--scale-aqueous-by-lwc", action="store_true",
                   help="Scale Phase 2a defaults by LWC fraction "
                        "(treat them as mol/L water rather than mol/m^3 cell).")
    p.add_argument("--box-model-seeds", action="store_true",
                   help="Replace aqueous pre-state with box-model-style "
                        "tiny algebraic seeds (Hp=1e-9, OHm=1e-14, "
                        "HSO3m=1e-12, SO3mm=1e-16). Tests proposed config fix.")
    args = p.parse_args()

    config_path = MPAS_CONFIG if args.mpas_config else args.config
    if not os.path.exists(config_path):
        sys.exit(f"config not found: {config_path}")

    header, species, rates = _parse_dump(args.csv)
    print(f"Loaded {os.path.basename(args.csv)}")
    print(f"  iCell={header.get('iCell')} k={header.get('k')} "
          f"lat={header.get('lat_deg')} lon={header.get('lon_deg')}")
    print(f"  T={header.get('T_K')} K  P={header.get('P_Pa')} Pa  "
          f"air={header.get('air_mol_per_m3')} mol/m^3")
    print(f"  species rows: {len(species)}  rate_param rows: {len(rates)}")
    print(f"  config: {config_path}")

    # LWC diagnostic: what's the implied [Hp] in mol/L water?
    lwc_mol_per_m3 = species.get("CLOUD.AQUEOUS.H2O", (0.0, 0.0))[0]
    if lwc_mol_per_m3 > 0:
        # mol H2O × 0.018 kg/mol / 997 kg/m^3 = m^3_water per m^3_cell
        lwc_vol_frac = lwc_mol_per_m3 * 0.018015 / 997.0
        Hp_mol_per_m3 = species.get("CLOUD.AQUEOUS.Hp", (0.0, 0.0))[0]
        if Hp_mol_per_m3 > 0:
            Hp_mol_per_L_water = Hp_mol_per_m3 / lwc_vol_frac / 1000.0
            print(f"  LWC = {lwc_mol_per_m3:.3e} mol/m^3 cell  "
                  f"= {lwc_vol_frac*1e6:.3f} cm^3_water/m^3_cell  "
                  f"({lwc_vol_frac*997.0*1000:.3f} g/m^3)")
            print(f"  [Hp] = {Hp_mol_per_m3:.3e} mol/m^3 cell  "
                  f"= {Hp_mol_per_L_water:.3e} mol/L water  "
                  f"(pH={-math.log10(max(Hp_mol_per_L_water, 1e-30)):.2f})")

    if args.zero_aqueous:
        for n in list(species.keys()):
            if n.startswith("CLOUD.AQUEOUS.") and n != "CLOUD.AQUEOUS.H2O":
                species[n] = (0.0, species[n][1])
        print("  ** zero-aqueous: cleared all CLOUD.AQUEOUS.* except H2O")
    if args.scale_aqueous_by_lwc and lwc_mol_per_m3 > 0:
        lwc_vol_frac = lwc_mol_per_m3 * 0.018015 / 997.0
        scale = lwc_vol_frac * 1000.0  # mol/L_water -> mol/m^3_cell
        for n in list(species.keys()):
            if n.startswith("CLOUD.AQUEOUS.") and n != "CLOUD.AQUEOUS.H2O":
                old = species[n][0]
                species[n] = (old * scale, species[n][1])
        print(f"  ** scale-aqueous-by-lwc: factor {scale:.3e}")
    if args.box_model_seeds:
        seeds = {
            "CLOUD.AQUEOUS.Hp":      1e-9,
            "CLOUD.AQUEOUS.OHm":     1e-14,
            "CLOUD.AQUEOUS.HSO3m":   1e-12,
            "CLOUD.AQUEOUS.SO3mm":   1e-16,
            "CLOUD.AQUEOUS.SO2":     1e-12,
            "CLOUD.AQUEOUS.H2O2":    1e-12,
            "CLOUD.AQUEOUS.O3":      1e-14,
            "CLOUD.AQUEOUS.SO2OOHm": 0.0,
        }
        for n, v in seeds.items():
            if n in species:
                species[n] = (v, species[n][1])
        print(f"  ** box-model-seeds: replaced {len(seeds)} aqueous species "
              f"with tiny algebraic seeds")

    # Show what MPAS produced post-solve (ground truth blowup)
    bad_mpas = [(n, post) for n, (_pre, post) in species.items()
                if (not math.isfinite(post)) or abs(post) > 1.0e6]
    bad_mpas.sort(key=lambda r: abs(r[1]) if math.isfinite(r[1]) else float('inf'), reverse=True)
    print(f"\n[MPAS ground-truth post-solve] # bad species = {len(bad_mpas)}")
    for n, v in bad_mpas[:10]:
        print(f"  {n:35s} = {v: .6e}")

    if args.mode in ("tuned", "both"):
        micm_t, _ = _build_mpas_tuned_micm(config_path)
        _run_one("TUNED", micm_t, header, species, rates, args.dt, args.n_substeps)

    if args.mode in ("default", "both"):
        micm_d, _ = _build_default_micm(config_path)
        _run_one("DEFAULT", micm_d, header, species, rates, args.dt, args.n_substeps)


if __name__ == "__main__":
    main()
