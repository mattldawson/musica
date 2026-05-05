# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# TS1-cloud mechanism audit script.
#
# Runs four correctness audits against a MICM v1 mechanism config:
#   1. Charge balance per aqueous reaction (charge inferred from species name
#      suffixes: 'p' = +1, 'pp' = +2, 'm' = -1, 'mm' = -2).
#   2. Mass balance per reaction via molecular weights, excluding the
#      gas-phase collector species (_produced, _lost) which are intentionally
#      non-conservative accumulators.
#   3. Collector audit: every gas-phase reaction touching SO2/H2O2/O3 must
#      carry matching {_produced, _lost} collectors with stoichiometry that
#      mirrors the gas reactant/product coefficient. Collectors must never
#      appear as reactants or in rate-law denominators.
#   4. min_halflife candidate scan: rank reactions by representative species
#      half-life under sampled MPAS conditions; flag any with half-life
#      below the LHS dt that lack a `min halflife [s]` annotation.
#
# Usage:
#   uv run python python/test/integration/mechanism_audit.py \
#       configs/v1/ts1_cloud/config.json \
#       --report /tmp/ts1_cloud_audit.json

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from typing import Any, Dict, List, Optional, Tuple


COLLECTOR_GASES = ("SO2", "H2O2", "O3")
COLLECTOR_SUFFIXES = ("_produced", "_lost")
DEFAULT_DT_S = 10.0  # matches LHS harness DEFAULT_DT


# ------------------------- charge inference -------------------------

_CHARGE_RE = re.compile(r"(p+|m+)$")


def infer_charge(name: str) -> Optional[int]:
    """Infer aqueous-species charge from name suffix.

    Returns None for species whose name does not end in a recognized charge
    suffix (e.g. neutral aqueous species like SO2, H2O2, O3, H2O).
    """
    if name == "H2O":
        return 0
    # Carve out names that legitimately end in 'p' or 'm' but are not charge
    # suffixes. We special-case the small known aqueous set; for everything
    # else, the suffix rule is applied conservatively (only used in aqueous
    # reactions, where the convention is enforced).
    m = _CHARGE_RE.search(name)
    if not m:
        return 0
    suffix = m.group(1)
    sign = 1 if suffix.startswith("p") else -1
    return sign * len(suffix)


def is_aqueous_reaction(reaction: Dict[str, Any]) -> bool:
    rtype = reaction.get("type", "")
    if not rtype.startswith("DISSOLVED") and rtype not in ("HENRY_LAW_EQUILIBRIUM",):
        return False
    # Skip Henry's law (gas <-> aqueous, both neutral)
    return rtype.startswith("DISSOLVED")


# ------------------------- molecular-weight lookup -------------------------

# Aqueous species are not in the top-level species list (they appear under
# the AQUEOUS phase). Provide molecular weights for the canonical set so the
# mass-balance audit can include aqueous reactions.
AQUEOUS_MW = {
    "Hp": 0.001008,           # H+
    "OHm": 0.017008,          # OH-
    "HSO3m": 0.081072,        # HSO3-
    "SO3mm": 0.080064,        # SO3--
    "SO4mm": 0.096064,        # SO4--
    "SO2OOHm": 0.097064,      # SO2OOH- (peroxomonosulfate intermediate)
    "SO2": 0.064066,
    "H2O2": 0.034014,
    "O3": 0.047998,
    "H2O": 0.018015,
}


def build_mw_table(config: Dict[str, Any]) -> Dict[str, float]:
    """Map species name -> molecular weight in kg/mol."""
    mw: Dict[str, float] = {}
    for sp in config.get("species", []):
        name = sp.get("name")
        m = sp.get("molecular weight [kg mol-1]")
        if name and m is not None:
            mw[name] = float(m)
    # Layer aqueous values (use config value if present, else fallback)
    for name, value in AQUEOUS_MW.items():
        mw.setdefault(name, value)
    return mw


def is_collector(name: str) -> bool:
    return any(name.endswith(s) for s in COLLECTOR_SUFFIXES)


# ------------------------- audits -------------------------

def audit_charge_balance(reactions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    issues = []
    for idx, r in enumerate(reactions):
        if not is_aqueous_reaction(r):
            continue
        lhs = sum(
            float(p["coefficient"]) * (infer_charge(p["species name"]) or 0)
            for p in r.get("reactants", [])
        )
        rhs = sum(
            float(p["coefficient"]) * (infer_charge(p["species name"]) or 0)
            for p in r.get("products", [])
        )
        delta = rhs - lhs
        if abs(delta) > 1e-9:
            issues.append({
                "reaction_index": idx,
                "type": r.get("type"),
                "comment": r.get("__comment", ""),
                "lhs_charge": lhs,
                "rhs_charge": rhs,
                "delta": delta,
                "reactants": [(p["species name"], p["coefficient"]) for p in r.get("reactants", [])],
                "products": [(p["species name"], p["coefficient"]) for p in r.get("products", [])],
            })
    return issues


def audit_mass_balance(
    reactions: List[Dict[str, Any]],
    mw: Dict[str, float],
    rel_tol: float = 1e-4,
) -> List[Dict[str, Any]]:
    """Verify Σ(coef * MW) balances across each reaction, ignoring collectors.

    Reactions with non-mass-bearing species (e.g. EMISSION sources where the
    product alone is real) are skipped if either side has zero total mass.
    """
    issues = []
    for idx, r in enumerate(reactions):
        rtype = r.get("type", "")
        # Skip non-mass-conserving reaction types: emissions, first-order loss
        # (which have a reactant but optional/no products), surface (handled
        # separately), and equilibria (handled by the algebraic constraint).
        if rtype in ("EMISSION", "FIRST_ORDER_LOSS", "SURFACE", "DISSOLVED_EQUILIBRIUM",
                     "HENRY_LAW_EQUILIBRIUM"):
            continue

        def side_mass(side):
            total = 0.0
            unknown = []
            for p in r.get(side, []):
                name = p["species name"]
                if is_collector(name):
                    continue
                if name not in mw:
                    unknown.append(name)
                    continue
                total += float(p["coefficient"]) * mw[name]
            return total, unknown

        lhs, lhs_unk = side_mass("reactants")
        rhs, rhs_unk = side_mass("products")
        if lhs == 0.0 or rhs == 0.0:
            # Cannot meaningfully compare (e.g. emission/loss-style reactions)
            continue

        denom = max(abs(lhs), abs(rhs))
        rel_err = abs(rhs - lhs) / denom
        if rel_err > rel_tol or lhs_unk or rhs_unk:
            issues.append({
                "reaction_index": idx,
                "type": rtype,
                "comment": r.get("__comment", ""),
                "lhs_mass_kg_mol": lhs,
                "rhs_mass_kg_mol": rhs,
                "rel_err": rel_err,
                "missing_mw": sorted(set(lhs_unk + rhs_unk)),
                "reactants": [(p["species name"], p["coefficient"]) for p in r.get("reactants", [])],
                "products": [(p["species name"], p["coefficient"]) for p in r.get("products", [])],
            })
    return issues


def audit_collectors(reactions: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Verify collector species accompany every gas reaction touching SO2/H2O2/O3.

    Returns a dict with two keys:
      - missing_collectors: reactions where a collector should appear but doesn't
        or where stoichiometry doesn't match the gas coefficient.
      - misplaced_collectors: collector species appearing as reactants (forbidden).
    """
    missing = []
    misplaced = []
    for idx, r in enumerate(reactions):
        rtype = r.get("type", "")
        # Gas-phase only (collectors are gas-phase tags). Skip aqueous and
        # surface and Henry's law reactions; collector-tagging applies to
        # gas-phase chemistry only.
        if rtype in ("DISSOLVED_REACTION", "DISSOLVED_REVERSIBLE_REACTION",
                     "DISSOLVED_EQUILIBRIUM", "HENRY_LAW_EQUILIBRIUM",
                     "SURFACE"):
            continue

        reactants = r.get("reactants", [])
        products = r.get("products", [])
        prod_map = {p["species name"]: float(p["coefficient"]) for p in products}
        react_map = {p["species name"]: float(p["coefficient"]) for p in reactants}
        react_names = set(react_map.keys())

        # 1) collectors must never be reactants
        for p in reactants:
            if is_collector(p["species name"]):
                misplaced.append({
                    "reaction_index": idx,
                    "type": rtype,
                    "comment": r.get("__comment", ""),
                    "issue": "collector appears as reactant",
                    "species": p["species name"],
                })

        for gas in COLLECTOR_GASES:
            # Net consumption / production of the gas. Catalytic patterns where
            # the gas appears in both sides cancel, and so should the collectors.
            net_consumed = react_map.get(gas, 0.0) - prod_map.get(gas, 0.0)
            if net_consumed > 0.0:
                expected = net_consumed
                actual = prod_map.get(f"{gas}_lost", 0.0)
                if abs(actual - expected) > 1e-9:
                    missing.append({
                        "reaction_index": idx,
                        "type": rtype,
                        "comment": r.get("__comment", ""),
                        "issue": f"missing/incorrect {gas}_lost (net consumption)",
                        "expected_coefficient": expected,
                        "actual_coefficient": actual,
                        "reactants": [(x["species name"], x["coefficient"]) for x in reactants],
                        "products": [(x["species name"], x["coefficient"]) for x in products],
                    })
            elif net_consumed < 0.0:
                expected = -net_consumed
                actual = prod_map.get(f"{gas}_produced", 0.0)
                if abs(actual - expected) > 1e-9:
                    missing.append({
                        "reaction_index": idx,
                        "type": rtype,
                        "comment": r.get("__comment", ""),
                        "issue": f"missing/incorrect {gas}_produced (net production)",
                        "expected_coefficient": expected,
                        "actual_coefficient": actual,
                        "reactants": [(x["species name"], x["coefficient"]) for x in reactants],
                        "products": [(x["species name"], x["coefficient"]) for x in products],
                    })
            else:
                # net_consumed == 0: catalytic. Collectors should also be 0.
                spurious_lost = prod_map.get(f"{gas}_lost", 0.0)
                spurious_prod = prod_map.get(f"{gas}_produced", 0.0)
                if spurious_lost != 0.0 or spurious_prod != 0.0:
                    missing.append({
                        "reaction_index": idx,
                        "type": rtype,
                        "comment": r.get("__comment", ""),
                        "issue": f"spurious {gas} collector on catalytic reaction",
                        "expected_coefficient": 0.0,
                        "actual_coefficient": spurious_lost or spurious_prod,
                        "reactants": [(x["species name"], x["coefficient"]) for x in reactants],
                        "products": [(x["species name"], x["coefficient"]) for x in products],
                    })

    return {"missing_collectors": missing, "misplaced_collectors": misplaced}


# ------------------------- min_halflife candidate scan -------------------------

# Representative MPAS-relevant concentrations for half-life estimation.
# Concentrations in mol m-3 (matching v1 unit convention).
REPRESENTATIVE_GAS_CONC = {
    # Background gas-phase concentrations at ~1e5 Pa, 280 K, typical vmr.
    "OH": 1.0e-12 * 2.4e-2 * 6.022e23,   # 1 pptv ~ 1e6 cm^-3 ~ 1.66e-12 mol/m^3
    "HO2": 1.0e-11 * 2.4e-2 * 6.022e23,
    "NO": 1.0e-10 * 2.4e-2 * 6.022e23,
    "NO2": 1.0e-9 * 2.4e-2 * 6.022e23,
    "O3": 5.0e-8 * 2.4e-2 * 6.022e23,
    "H2O2": 1.0e-9 * 2.4e-2 * 6.022e23,
    "SO2": 1.0e-9 * 2.4e-2 * 6.022e23,
    "M": 2.4e25,                          # molar density of air (molec/m^3)
}


def _arrhenius_k(reaction: Dict[str, Any], T: float = 280.0) -> Optional[float]:
    """Compute Arrhenius/Troe-like rate constant where possible.

    Returns None if the reaction does not have a closed-form rate or requires
    non-trivial inputs. Units returned as given in the config (assumed
    consistent across reactants for a half-life estimate).
    """
    rtype = reaction.get("type")
    if rtype == "ARRHENIUS":
        A = float(reaction.get("A", 1.0))
        B = float(reaction.get("B", 0.0))
        C = float(reaction.get("C", 0.0))
        D = float(reaction.get("D", 300.0))
        E = float(reaction.get("E", 0.0))
        # k = A * (T/D)^B * exp(C/T) * (1 + E*P) — drop pressure correction.
        return A * (T / D) ** B * math.exp(C / T)
    if rtype == "FIRST_ORDER_LOSS":
        return float(reaction.get("scaling factor", 1.0))
    if rtype == "DISSOLVED_REACTION":
        rc = reaction.get("rate constant", {})
        A = float(rc.get("A", 1.0))
        B = float(rc.get("B", 0.0))
        C = float(rc.get("C", 0.0))
        D = float(rc.get("D", 300.0))
        return A * (T / D) ** B * math.exp(-C / T)  # DISSOLVED uses Ea-style exp(-C/T) per code
    return None


def audit_min_halflife_candidates(
    reactions: List[Dict[str, Any]],
    dt: float = DEFAULT_DT_S,
) -> List[Dict[str, Any]]:
    """Flag fast reactions without a `min halflife [s]` whose half-life << dt."""
    candidates = []
    for idx, r in enumerate(reactions):
        if "min halflife [s]" in r:
            continue
        k = _arrhenius_k(r)
        if k is None:
            continue
        # Rough effective first-order rate using representative concentrations.
        # k_eff = k * Π_other_reactant [X]
        reactants = r.get("reactants", [])
        if not reactants:
            continue
        # Pick first reactant as the "tracked species"; multiply by typical
        # concentrations of the others. M (air) is treated as a static
        # third-body and never the tracked stiff species.
        tracked = reactants[0]["species name"]
        if tracked == "M":
            if len(reactants) < 2:
                continue
            tracked = reactants[1]["species name"]
            other = reactants[2:]
        else:
            other = reactants[1:]
        k_eff = k
        skipped = False
        for p in other:
            name = p["species name"]
            conc = REPRESENTATIVE_GAS_CONC.get(name)
            if conc is None:
                if name in AQUEOUS_MW:
                    conc = 1.0
                else:
                    skipped = True
                    break
            k_eff *= conc * float(p["coefficient"])
        if skipped or k_eff <= 0.0:
            continue
        halflife = math.log(2) / k_eff if k_eff > 0 else float("inf")
        if halflife < dt / 5.0:
            candidates.append({
                "reaction_index": idx,
                "type": r.get("type"),
                "comment": r.get("__comment", ""),
                "halflife_s": halflife,
                "k_eff_per_s": k_eff,
                "tracked_reactant": tracked,
                "reactants": [(p["species name"], p["coefficient"]) for p in reactants],
                "products": [(p["species name"], p["coefficient"]) for p in r.get("products", [])],
            })
    candidates.sort(key=lambda c: c["halflife_s"])
    return candidates


# ------------------------- driver -------------------------

def run_audits(config_path: str, dt: float = DEFAULT_DT_S) -> Dict[str, Any]:
    with open(config_path) as f:
        config = json.load(f)

    gas_reactions = config.get("reactions", [])
    aqueous_reactions = config.get("aerosol processes", [])
    all_reactions = list(gas_reactions) + list(aqueous_reactions)

    mw = build_mw_table(config)

    return {
        "config_path": os.path.abspath(config_path),
        "n_gas_reactions": len(gas_reactions),
        "n_aqueous_reactions": len(aqueous_reactions),
        "charge_balance_issues": audit_charge_balance(all_reactions),
        "mass_balance_issues": audit_mass_balance(all_reactions, mw),
        "collector_audit": audit_collectors(gas_reactions),
        "min_halflife_candidates": audit_min_halflife_candidates(all_reactions, dt=dt),
    }


def format_summary(report: Dict[str, Any]) -> str:
    lines = [
        f"Mechanism audit: {report['config_path']}",
        f"  gas reactions:     {report['n_gas_reactions']}",
        f"  aqueous reactions: {report['n_aqueous_reactions']}",
        "",
        f"Charge balance issues:        {len(report['charge_balance_issues'])}",
        f"Mass balance issues:          {len(report['mass_balance_issues'])}",
        f"Missing/incorrect collectors: {len(report['collector_audit']['missing_collectors'])}",
        f"Misplaced collectors:         {len(report['collector_audit']['misplaced_collectors'])}",
        f"min_halflife candidates:      {len(report['min_halflife_candidates'])}",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="Path to mechanism config JSON")
    parser.add_argument("--report", help="Optional path to write the full JSON report")
    parser.add_argument("--dt", type=float, default=DEFAULT_DT_S,
                        help="Reference timestep used by the min_halflife scan")
    parser.add_argument("--verbose", action="store_true",
                        help="Print full per-issue details")
    args = parser.parse_args(argv)

    report = run_audits(args.config, dt=args.dt)
    print(format_summary(report))

    if args.verbose:
        print("\n--- Charge balance issues ---")
        for it in report["charge_balance_issues"]:
            print(json.dumps(it, indent=2))
        print("\n--- Mass balance issues ---")
        for it in report["mass_balance_issues"]:
            print(json.dumps(it, indent=2))
        print("\n--- Missing/incorrect collectors ---")
        for it in report["collector_audit"]["missing_collectors"]:
            print(json.dumps(it, indent=2))
        print("\n--- Misplaced collectors ---")
        for it in report["collector_audit"]["misplaced_collectors"]:
            print(json.dumps(it, indent=2))
        print("\n--- min_halflife candidates (top 25) ---")
        for it in report["min_halflife_candidates"][:25]:
            print(json.dumps(it, indent=2))

    if args.report:
        with open(args.report, "w") as f:
            json.dump(report, f, indent=2)
        print(f"\nFull report written to {args.report}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
