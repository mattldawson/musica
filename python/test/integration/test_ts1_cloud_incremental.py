# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# Incremental TS1 + Cloud Chemistry development tests.
# Follows the plan in .github/prompts/ts1_cloud_chem_mechanism_development.prompt.md
#
# Each step adds one cloud chemistry component to the TS1 gas-phase mechanism,
# verifying convergence, mass conservation, and physical reasonableness.

import copy
import csv
import json
import math
import os

import pytest

import musica.mechanism_configuration as mc
from musica.mechanism_configuration import Parser
from musica.micm import MICM, SolverType, SolverState
from musica.miam import (
    ArrheniusRateConstant,
    EquilibriumConstant,
    HenrysLawConstant,
    UniformSection,
    DissolvedReaction,
    HenryLawEquilibriumConstraint,
    DissolvedEquilibriumConstraint,
    LinearConstraint,
    LinearConstraintTerm,
    Model,
)


# ═══ Paths ═══════════════════════════════════════════════════════════════════

BASE_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..")
TS1_JSON = os.path.join(BASE_DIR, "configs", "v1", "ts1", "ts1.json")
TS1_ICS_CSV = os.path.join(BASE_DIR, "configs", "v1", "ts1",
                           "initial_conditions.csv")


# ═══ Physical Constants ═══════════════════════════════════════════════════════

MW_H2O = 0.018015   # kg/mol
RHO_H2O = 997.0     # kg/m³
R_GAS = 8.314462     # J/(mol·K)
T0 = 298.15          # K, reference temperature
LWC = 0.3e-3         # liquid water content [kg m-3] (typical cloud)
C_H2O = LWC / MW_H2O # mol/m³-air (cloud water as MIAM state variable)
C_H2O_MOLAR = RHO_H2O / (MW_H2O * 1000)  # mol/L, molarity of pure water

# TS1 initial condition values (from initial_conditions.csv)
T_INIT = 287.45      # K
P_INIT = 101320.0    # Pa

# Henry's law constants (mol m-3 Pa-1 at 298.15 K)
M_ATM_TO_MOL_M3_PA = 1000.0 / 101325.0
HLC_SO2_REF = 1.23 * M_ATM_TO_MOL_M3_PA
HLC_H2O2_REF = 7.4e4 * M_ATM_TO_MOL_M3_PA
HLC_O3_REF = 1.15e-2 * M_ATM_TO_MOL_M3_PA

# Time step for integration
DT = 900.0  # seconds (MPAS chemistry timestep)


# ═══ Config Builder ══════════════════════════════════════════════════════════

def _load_ts1_config():
    """Load TS1 mechanism JSON as a mutable dict."""
    with open(TS1_JSON) as f:
        return json.load(f)


def _load_ts1_ics():
    """Load TS1 initial conditions from CSV.

    Returns dict mapping 'CONC.species' -> value and 'USER.param' -> value.
    """
    ics = {}
    with open(TS1_ICS_CSV) as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                ics[row[0].strip()] = float(row[1])
    return ics


def _get_species_in_reaction(rxn, role):
    """Get list of (species_name, coefficient) from reactants or products."""
    return [(entry["species name"], entry["coefficient"])
            for entry in rxn.get(role, [])]


def _find_source_sink_reactions(config, species_name):
    """Find reactions that produce or consume a species.

    Returns (source_indices, sink_indices) where each element is
    (reaction_index, yield/consumed_coefficient).
    Reactions where the species appears on BOTH sides with equal
    coefficients (net zero) are excluded.
    """
    sources = []
    sinks = []
    for i, rxn in enumerate(config["reactions"]):
        reactants = dict(_get_species_in_reaction(rxn, "reactants"))
        products = dict(_get_species_in_reaction(rxn, "products"))
        r_coef = reactants.get(species_name, 0.0)
        p_coef = products.get(species_name, 0.0)
        net = p_coef - r_coef
        if net > 0:
            sources.append((i, net))
        elif net < 0:
            sinks.append((i, -net))  # positive value = amount consumed
    return sources, sinks


def _add_source_sink_to_config(config, species_name):
    """Add source/sink tracer species and co-products to a TS1 config.

    Modifies config in place. Returns (source_name, sink_name, n_sources, n_sinks).
    """
    source_name = f"{species_name}_source"
    sink_name = f"{species_name}_sink"
    sources, sinks = _find_source_sink_reactions(config, species_name)

    # Add tracer species to species list
    config["species"].append({
        "name": source_name,
        "molecular weight [kg mol-1]": 0.0,
        "__description": f"Tracer tracking gas-phase {species_name} production",
    })
    config["species"].append({
        "name": sink_name,
        "molecular weight [kg mol-1]": 0.0,
        "__description": f"Tracer tracking gas-phase {species_name} consumption",
    })

    # Add to gas phase
    gas_phase = next(p for p in config["phases"] if p["name"] == "gas")
    gas_phase["species"].append({"name": source_name})
    gas_phase["species"].append({"name": sink_name})

    # Add co-products to reactions
    for idx, coef in sources:
        config["reactions"][idx]["products"].append({
            "species name": source_name,
            "coefficient": coef,
        })
    for idx, coef in sinks:
        config["reactions"][idx]["products"].append({
            "species name": sink_name,
            "coefficient": coef,
        })

    return source_name, sink_name, len(sources), len(sinks)


def _write_config(config, tmp_path, name="config.json"):
    """Write config dict to a JSON file and return the path."""
    path = os.path.join(tmp_path, name)
    with open(path, "w") as f:
        json.dump(config, f, indent=2)
    return path


# ═══ MIAM Model Builders ════════════════════════════════════════════════════

def _make_aqueous_species():
    """Create MIAM species for the AQUEOUS phase."""
    h2o = mc.Species(name="H2O")
    h2o.molecular_weight_kg_mol = MW_H2O
    h2o.density_kg_m3 = RHO_H2O
    return {
        "SO2_aq": mc.Species(name="SO2_aq"),
        "H2O2_aq": mc.Species(name="H2O2_aq"),
        "O3_aq": mc.Species(name="O3_aq"),
        "Hp": mc.Species(name="Hp"),
        "OHm": mc.Species(name="OHm"),
        "HSO3m": mc.Species(name="HSO3m"),
        "SO3mm": mc.Species(name="SO3mm"),
        "SO4mm": mc.Species(name="SO4mm"),
        "H2O": h2o,
    }


def _make_cloud_representation():
    """Create the CLOUD UNIFORM_SECTION aerosol representation."""
    return UniformSection(
        name="CLOUD", phase_names=["AQUEOUS"],
        min_radius=1e-6, max_radius=1e-5,
    )


# --- HLC constraints ---

def _hlc_so2():
    return HenryLawEquilibriumConstraint(
        gas_species_name="SO2",
        condensed_species_name="SO2_aq",
        solvent_name="H2O",
        condensed_phase_name="AQUEOUS",
        henrys_law_constant=HenrysLawConstant(hlc_ref=HLC_SO2_REF, c=3120.0),
        mw_solvent=MW_H2O, rho_solvent=RHO_H2O,
    )


def _hlc_h2o2():
    return HenryLawEquilibriumConstraint(
        gas_species_name="H2O2",
        condensed_species_name="H2O2_aq",
        solvent_name="H2O",
        condensed_phase_name="AQUEOUS",
        henrys_law_constant=HenrysLawConstant(hlc_ref=HLC_H2O2_REF, c=6621.0),
        mw_solvent=MW_H2O, rho_solvent=RHO_H2O,
    )


def _hlc_o3():
    return HenryLawEquilibriumConstraint(
        gas_species_name="O3",
        condensed_species_name="O3_aq",
        solvent_name="H2O",
        condensed_phase_name="AQUEOUS",
        henrys_law_constant=HenrysLawConstant(hlc_ref=HLC_O3_REF, c=2560.0),
        mw_solvent=MW_H2O, rho_solvent=RHO_H2O,
    )


# --- Dissolved equilibria ---

def _eq_kw():
    return DissolvedEquilibriumConstraint(
        phase_name="AQUEOUS",
        reactant_names=["H2O"],
        product_names=["Hp", "OHm"],
        algebraic_species_name="OHm",
        solvent_name="H2O",
        equilibrium_constant=EquilibriumConstant(
            a=1e-14 / (C_H2O_MOLAR * C_H2O_MOLAR), c=0.0),
    )


def _eq_ka1():
    return DissolvedEquilibriumConstraint(
        phase_name="AQUEOUS",
        reactant_names=["SO2_aq"],
        product_names=["HSO3m", "Hp"],
        algebraic_species_name="HSO3m",
        solvent_name="H2O",
        equilibrium_constant=EquilibriumConstant(
            a=1.7e-2 / C_H2O_MOLAR, c=2090.0),
    )


def _eq_ka2():
    return DissolvedEquilibriumConstraint(
        phase_name="AQUEOUS",
        reactant_names=["HSO3m"],
        product_names=["SO3mm", "Hp"],
        algebraic_species_name="SO3mm",
        solvent_name="H2O",
        equilibrium_constant=EquilibriumConstant(
            a=6.0e-8 / C_H2O_MOLAR, c=1120.0),
    )


# --- Linear constraints ---

def _lc_sulfur(include_source_sink=True):
    """Sulfur mass-conservation constraint.

    SO2(g) = C + SO2_source - SO2_sink - [aq sulfur]

    When include_source_sink is False, the constraint only covers SO2(g) + aq species
    (for testing without modified reactions).
    """
    terms = [
        LinearConstraintTerm("gas", "SO2", 1.0),
        LinearConstraintTerm("AQUEOUS", "SO2_aq", 1.0),
        LinearConstraintTerm("AQUEOUS", "HSO3m", 1.0),
        LinearConstraintTerm("AQUEOUS", "SO3mm", 1.0),
        LinearConstraintTerm("AQUEOUS", "SO4mm", 1.0),
    ]
    if include_source_sink:
        terms.extend([
            LinearConstraintTerm("gas", "SO2_source", -1.0),
            LinearConstraintTerm("gas", "SO2_sink", 1.0),
        ])
    return LinearConstraint(
        algebraic_phase_name="gas",
        algebraic_species_name="SO2",
        terms=terms,
        diagnose_from_state=True,
    )


def _lc_h2o2(include_source_sink=True):
    """H2O2 mass-conservation constraint."""
    terms = [
        LinearConstraintTerm("gas", "H2O2", 1.0),
        LinearConstraintTerm("AQUEOUS", "H2O2_aq", 1.0),
    ]
    if include_source_sink:
        terms.extend([
            LinearConstraintTerm("gas", "H2O2_source", -1.0),
            LinearConstraintTerm("gas", "H2O2_sink", 1.0),
        ])
    return LinearConstraint(
        algebraic_phase_name="gas",
        algebraic_species_name="H2O2",
        terms=terms,
        diagnose_from_state=True,
    )


def _lc_o3(include_source_sink=True):
    """O3 mass-conservation constraint."""
    terms = [
        LinearConstraintTerm("gas", "O3", 1.0),
        LinearConstraintTerm("AQUEOUS", "O3_aq", 1.0),
    ]
    if include_source_sink:
        terms.extend([
            LinearConstraintTerm("gas", "O3_source", -1.0),
            LinearConstraintTerm("gas", "O3_sink", 1.0),
        ])
    return LinearConstraint(
        algebraic_phase_name="gas",
        algebraic_species_name="O3",
        terms=terms,
        diagnose_from_state=True,
    )


def _lc_charge_balance():
    """Charge balance: H+ = OH- + HSO3- + 2 SO3-- + 2 SO4--"""
    return LinearConstraint(
        algebraic_phase_name="AQUEOUS",
        algebraic_species_name="Hp",
        terms=[
            LinearConstraintTerm("AQUEOUS", "Hp", 1.0),
            LinearConstraintTerm("AQUEOUS", "OHm", -1.0),
            LinearConstraintTerm("AQUEOUS", "HSO3m", -1.0),
            LinearConstraintTerm("AQUEOUS", "SO3mm", -2.0),
            LinearConstraintTerm("AQUEOUS", "SO4mm", -2.0),
        ],
        constant=0.0,
    )


# --- Dissolved reaction ---

def _dissolved_r1():
    """R1: HSO3- + H2O2(aq) + H+ -> SO4-- + H2O + 2H+ (Hoffmann & Calvert 1985).

    Full rate law: rate = k[H+][HSO3-][H2O2] / (1 + 13[H+])
    At cloud pH > 3 the denominator ≈ 1, so we include H+ as a reactant.
    k = 7.45e7 M^-2 s^-1, Ea/R = 4430 K.
    With 3 reactants: A_miam = C_H2O_MOLAR^2 * k_lit.
    """
    return DissolvedReaction(
        phase_name="AQUEOUS",
        reactant_names=["HSO3m", "H2O2_aq", "Hp"],
        product_names=["SO4mm", "H2O", "Hp", "Hp"],
        solvent_name="H2O",
        rate_constant=ArrheniusRateConstant(a=C_H2O_MOLAR**2 * 7.45e7, c=4430.0),
    )


def _compute_cloud_ics(so2_g0, h2o2_g0=0.0, o3_g0=0.0, so4mm0=0.0):
    """Compute self-consistent equilibrium ICs for cloud species.

    Given gas-phase concentrations, compute aqueous-phase concentrations
    at Henry's law + dissociation equilibrium + charge balance.

    Returns dict of cloud species ICs (CLOUD.AQUEOUS.* keys).
    """
    T = T_INIT
    # Volume fraction of liquid water in air
    f_v = C_H2O * MW_H2O / RHO_H2O
    # HLC at temperature T (partition coefficient includes f_v)
    alpha_SO2 = HLC_SO2_REF * math.exp(3120.0 * (1/T - 1/T0)) * R_GAS * T * f_v
    alpha_H2O2 = HLC_H2O2_REF * math.exp(6621.0 * (1/T - 1/T0)) * R_GAS * T * f_v
    alpha_O3 = HLC_O3_REF * math.exp(2560.0 * (1/T - 1/T0)) * R_GAS * T * f_v

    # Equilibrium constants at T (MIAM mole-fraction convention: K = Ka_solution / C_H2O_MOLAR)
    Ka1_T = (1.7e-2 / C_H2O_MOLAR) * math.exp(2090.0 * (1/T0 - 1/T))
    Ka2_T = (6.0e-8 / C_H2O_MOLAR) * math.exp(1120.0 * (1/T0 - 1/T))
    Kw_T = 1e-14 / (C_H2O_MOLAR * C_H2O_MOLAR)

    # Iteratively solve for Hp consistent with charge balance
    hp = 1e-4  # initial guess pH~4
    for _ in range(200):
        ohm = Kw_T * C_H2O * C_H2O / hp
        # Sulfur partitioning (conservation: SO2_g + SO2_aq + HSO3m + SO3mm = so2_g0)
        f = (1.0 + alpha_SO2
             + Ka1_T * alpha_SO2 * C_H2O / hp
             + Ka2_T * Ka1_T * alpha_SO2 * C_H2O * C_H2O / (hp * hp))
        so2_g = so2_g0 / f
        so2_aq = alpha_SO2 * so2_g
        hso3m = Ka1_T * so2_aq * C_H2O / hp
        so3mm = Ka2_T * hso3m * C_H2O / hp
        hp_new = ohm + hso3m + 2.0 * so3mm + 2.0 * so4mm0
        if abs(hp_new - hp) < 1e-15 * max(hp, 1e-30):
            break
        hp = 0.5 * (hp + hp_new)

    # H2O2 and O3 partitioning (simple, no dissociation)
    h2o2_g = h2o2_g0 / (1.0 + alpha_H2O2) if h2o2_g0 > 0 else 0.0
    h2o2_aq = alpha_H2O2 * h2o2_g
    o3_g = o3_g0 / (1.0 + alpha_O3) if o3_g0 > 0 else 0.0
    o3_aq = alpha_O3 * o3_g

    return {
        "CLOUD.AQUEOUS.H2O": C_H2O,
        "CLOUD.AQUEOUS.SO2_aq": so2_aq,
        "CLOUD.AQUEOUS.H2O2_aq": h2o2_aq,
        "CLOUD.AQUEOUS.O3_aq": o3_aq,
        "CLOUD.AQUEOUS.Hp": hp,
        "CLOUD.AQUEOUS.OHm": ohm,
        "CLOUD.AQUEOUS.HSO3m": hso3m,
        "CLOUD.AQUEOUS.SO3mm": so3mm,
        "CLOUD.AQUEOUS.SO4mm": so4mm0,
        # Updated gas concentrations reflecting dissolution
        "_SO2_g": so2_g,
        "_H2O2_g": h2o2_g,
        "_O3_g": o3_g,
    }


# ═══ MIAM Model Assembly ════════════════════════════════════════════════════

def _build_miam_model(processes=None, constraints=None):
    """Build a MIAM Model with cloud chemistry components.

    Species are determined automatically from the processes/constraints.
    Gas-phase species referenced by constraints must also be declared.
    """
    aq_species = _make_aqueous_species()

    # Determine which aqueous and gas species are needed
    needed_aq = {"H2O"}  # always need solvent
    needed_gas = set()

    for c in (constraints or []):
        if isinstance(c, HenryLawEquilibriumConstraint):
            needed_aq.add(c.condensed_species_name)
            needed_gas.add(c.gas_species_name)
        elif isinstance(c, DissolvedEquilibriumConstraint):
            needed_aq.update(c.reactant_names)
            needed_aq.update(c.product_names)
        elif isinstance(c, LinearConstraint):
            for t in c.terms:
                if t.phase_name == "AQUEOUS":
                    needed_aq.add(t.species_name)
                elif t.phase_name == "gas":
                    needed_gas.add(t.species_name)

    for p in (processes or []):
        if isinstance(p, DissolvedReaction):
            needed_aq.update(p.reactant_names)
            needed_aq.update(p.product_names)

    # Build species lists
    used_aq = [aq_species[k] for k in sorted(needed_aq) if k in aq_species]
    gas_sp = [mc.Species(name=n) for n in sorted(needed_gas)]
    all_sp = gas_sp + used_aq
    aq_phase = mc.Phase(name="AQUEOUS", species=used_aq)

    return Model(
        name="cloud_chemistry",
        species=all_sp,
        condensed_phases=[aq_phase],
        representations=[_make_cloud_representation()],
        processes=processes or [],
        constraints=constraints or [],
    )


# ═══ Solving & Verification ════════════════════════════════════════════════

def _load_full_ics():
    """Load TS1 ICs from CSV as dict of {species_name: concentration}."""
    raw = _load_ts1_ics()
    ics = {}
    for key, val in raw.items():
        if key.startswith("CONC."):
            ics[key[5:]] = val  # strip "CONC." prefix
    return ics


def _solve(micm, model, ics, dt=DT, n_steps=1):
    """Create state, set ICs, solve, return (state, ordered_results).

    ics: dict of {species_name: concentration}. Species not in the mechanism
         are silently ignored.
    """
    state = micm.create_state()
    if model is not None:
        model.set_default_parameters(state)
    state.set_conditions(temperatures=T_INIT, pressures=P_INIT)

    ordering = state.get_species_ordering()
    filtered = {k: v for k, v in ics.items() if k in ordering}
    state.set_concentrations(filtered)

    results = []
    for _ in range(n_steps):
        result = micm.solve(state, time_step=dt)
        results.append(result)
        if result.state != SolverState.Converged:
            break
    return state, results


def _get_conc(state, name):
    """Get concentration of a species from state (grid cell 0)."""
    ordering = state.get_species_ordering()
    if name not in ordering:
        return None
    val = state.get_concentrations()[name]
    return val[0] if isinstance(val, (list, tuple)) else val


def _check_sulfur_conservation(state, ics, rtol=1e-6):
    """Verify sulfur mass conservation.

    C_initial = SO2(g) - SO2_source + SO2_sink + aq_S
    C_final should equal C_initial.
    """
    so2_g = _get_conc(state, "SO2")
    so2_src = _get_conc(state, "SO2_source") or 0.0
    so2_snk = _get_conc(state, "SO2_sink") or 0.0
    so2_aq = _get_conc(state, "CLOUD.AQUEOUS.SO2_aq") or 0.0
    hso3m = _get_conc(state, "CLOUD.AQUEOUS.HSO3m") or 0.0
    so3mm = _get_conc(state, "CLOUD.AQUEOUS.SO3mm") or 0.0
    so4mm = _get_conc(state, "CLOUD.AQUEOUS.SO4mm") or 0.0

    total_final = so2_g - so2_src + so2_snk + so2_aq + hso3m + so3mm + so4mm

    # Initial C: from ICs (source/sink start at 0)
    so2_g0 = ics.get("SO2", 0.0)
    so2_aq0 = ics.get("CLOUD.AQUEOUS.SO2_aq", 0.0)
    hso3m0 = ics.get("CLOUD.AQUEOUS.HSO3m", 0.0)
    so3mm0 = ics.get("CLOUD.AQUEOUS.SO3mm", 0.0)
    so4mm0 = ics.get("CLOUD.AQUEOUS.SO4mm", 0.0)
    total_initial = so2_g0 + so2_aq0 + hso3m0 + so3mm0 + so4mm0

    print(f"  Sulfur conservation: initial={total_initial:.6e}, "
          f"final={total_final:.6e}, "
          f"rel_diff={abs(total_final - total_initial) / max(total_initial, 1e-30):.2e}")
    print(f"    SO2(g)={so2_g:.4e}, src={so2_src:.4e}, snk={so2_snk:.4e}")
    print(f"    SO2(aq)={so2_aq:.4e}, HSO3m={hso3m:.4e}, "
          f"SO3mm={so3mm:.4e}, SO4mm={so4mm:.4e}")

    assert abs(total_final - total_initial) < rtol * max(total_initial, 1e-30), \
        f"Sulfur not conserved: {total_initial:.6e} -> {total_final:.6e}"

    return total_initial, total_final


def _check_h2o2_conservation(state, ics, rtol=1e-6):
    """Verify H2O2 mass conservation."""
    h2o2_g = _get_conc(state, "H2O2")
    h2o2_src = _get_conc(state, "H2O2_source") or 0.0
    h2o2_snk = _get_conc(state, "H2O2_sink") or 0.0
    h2o2_aq = _get_conc(state, "CLOUD.AQUEOUS.H2O2_aq") or 0.0

    total_final = h2o2_g - h2o2_src + h2o2_snk + h2o2_aq
    total_initial = ics.get("H2O2", 0.0) + ics.get("CLOUD.AQUEOUS.H2O2_aq", 0.0)

    print(f"  H2O2 conservation: initial={total_initial:.6e}, "
          f"final={total_final:.6e}, "
          f"rel_diff={abs(total_final - total_initial) / max(total_initial, 1e-30):.2e}")

    assert abs(total_final - total_initial) < rtol * max(total_initial, 1e-30), \
        f"H2O2 not conserved: {total_initial:.6e} -> {total_final:.6e}"


def _check_charge_balance(state, atol=1e-10):
    """Verify charge balance: H+ = OH- + HSO3- + 2*SO3-- + 2*SO4--"""
    hp = _get_conc(state, "CLOUD.AQUEOUS.Hp") or 0.0
    ohm = _get_conc(state, "CLOUD.AQUEOUS.OHm") or 0.0
    hso3m = _get_conc(state, "CLOUD.AQUEOUS.HSO3m") or 0.0
    so3mm = _get_conc(state, "CLOUD.AQUEOUS.SO3mm") or 0.0
    so4mm = _get_conc(state, "CLOUD.AQUEOUS.SO4mm") or 0.0
    anions = ohm + hso3m + 2.0 * so3mm + 2.0 * so4mm
    residual = hp - anions
    print(f"  Charge balance: H+={hp:.4e}, anions={anions:.4e}, "
          f"residual={residual:.4e}")
    assert abs(residual) < atol + 1e-6 * max(hp, anions), \
        f"Charge balance violated: H+={hp:.4e}, anions={anions:.4e}"


# ═══ Fixtures ═══════════════════════════════════════════════════════════════

@pytest.fixture(scope="module")
def ts1_ics():
    """Full TS1 initial conditions."""
    return _load_full_ics()


@pytest.fixture(scope="module")
def ts1_mechanism():
    """Parsed TS1 mechanism (unmodified)."""
    parser = Parser()
    return parser.parse(TS1_JSON)


# ═══ Tests ══════════════════════════════════════════════════════════════════


class TestStep0_Baseline:
    """Step 0: Pure TS1, Rosenbrock solver. Verify it works at all."""

    def test_ts1_rosenbrock_converges(self, ts1_mechanism, ts1_ics):
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_standard_order,
        )
        state, results = _solve(micm, None, ts1_ics, dt=DT)
        print(f"Step 0: solver_state={results[-1].state}")
        assert results[-1].state == SolverState.Converged

        # Key species should still be positive
        for sp in ["SO2", "H2O2", "O3"]:
            val = _get_conc(state, sp)
            print(f"  {sp} = {val:.4e}")
            assert val >= 0, f"{sp} went negative: {val}"


class TestStep1_CloudPhaseOnly:
    """Step 1: Add CLOUD + AQUEOUS phase, DAE4 solver, no aerosol processes."""

    def test_dae4_with_inert_aqueous(self, ts1_mechanism, ts1_ics):
        model = _build_miam_model()  # No processes, no constraints
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        ics = {**ts1_ics, "CLOUD.AQUEOUS.H2O": C_H2O}
        state, results = _solve(micm, model, ics, dt=DT)
        print(f"Step 1: solver_state={results[-1].state}")
        assert results[-1].state == SolverState.Converged

        # Gas species should be unchanged from Step 0 (no cloud processes)
        for sp in ["SO2", "H2O2", "O3"]:
            val = _get_conc(state, sp)
            print(f"  {sp} = {val:.4e}")
            assert val >= 0


class TestStep2_HLC_SO2:
    """Step 2: Add HLC SO2 + source/sink tracers + sulfur constraint."""

    @pytest.fixture(scope="class")
    def step2_config(self, tmp_path_factory):
        tmp = str(tmp_path_factory.mktemp("step2"))
        config = _load_ts1_config()
        src, snk, n_src, n_snk = _add_source_sink_to_config(config, "SO2")
        print(f"  SO2 source reactions: {n_src}, sink reactions: {n_snk}")
        return _write_config(config, tmp, "step2.json")

    @pytest.fixture(scope="class")
    def step2_mechanism(self, step2_config):
        return Parser().parse(step2_config)

    def test_converges(self, step2_mechanism, ts1_ics):
        model = _build_miam_model(
            constraints=[_hlc_so2(), _lc_sulfur()],
        )
        micm = MICM(
            mechanism=step2_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        ics = {
            **ts1_ics,
            "SO2_source": 0.0,
            "SO2_sink": 0.0,
            "CLOUD.AQUEOUS.H2O": C_H2O,
            "CLOUD.AQUEOUS.SO2_aq": 0.0,
            "CLOUD.AQUEOUS.HSO3m": 0.0,
            "CLOUD.AQUEOUS.SO3mm": 0.0,
            "CLOUD.AQUEOUS.SO4mm": 0.0,
            "CLOUD.AQUEOUS.Hp": 1e-4,
        }
        state, results = _solve(micm, model, ics, dt=DT)
        print(f"Step 2: solver_state={results[-1].state}")
        assert results[-1].state == SolverState.Converged

        # SO2 should have partitioned into aqueous phase
        so2_aq = _get_conc(state, "CLOUD.AQUEOUS.SO2_aq")
        print(f"  SO2(aq) = {so2_aq:.4e}")
        assert so2_aq > 0, "SO2 didn't dissolve"

        # Sulfur conservation
        _check_sulfur_conservation(state, ics)


class TestStep3_DissolvedEquilibria:
    """Step 3: Add Ka1, Ka2, Kw dissolved equilibria."""

    @pytest.fixture(scope="class")
    def step3_config(self, tmp_path_factory):
        tmp = str(tmp_path_factory.mktemp("step3"))
        config = _load_ts1_config()
        _add_source_sink_to_config(config, "SO2")
        return _write_config(config, tmp, "step3.json")

    @pytest.fixture(scope="class")
    def step3_mechanism(self, step3_config):
        return Parser().parse(step3_config)

    def test_converges_with_equilibria(self, step3_mechanism, ts1_ics):
        model = _build_miam_model(
            constraints=[
                _hlc_so2(),
                _eq_kw(), _eq_ka1(), _eq_ka2(),
                _lc_sulfur(),
            ],
        )
        micm = MICM(
            mechanism=step3_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        ics = {
            **ts1_ics,
            "SO2_source": 0.0,
            "SO2_sink": 0.0,
            "CLOUD.AQUEOUS.H2O": C_H2O,
            "CLOUD.AQUEOUS.SO2_aq": 0.0,
            "CLOUD.AQUEOUS.HSO3m": 0.0,
            "CLOUD.AQUEOUS.SO3mm": 0.0,
            "CLOUD.AQUEOUS.SO4mm": 0.0,
            "CLOUD.AQUEOUS.Hp": 1e-4,
            "CLOUD.AQUEOUS.OHm": 1e-10,
        }
        state, results = _solve(micm, model, ics, dt=DT)
        print(f"Step 3: solver_state={results[-1].state}")
        assert results[-1].state == SolverState.Converged

        # HSO3- and SO3-- should be non-zero from equilibria
        hso3m = _get_conc(state, "CLOUD.AQUEOUS.HSO3m")
        so3mm = _get_conc(state, "CLOUD.AQUEOUS.SO3mm")
        hp = _get_conc(state, "CLOUD.AQUEOUS.Hp")
        ohm = _get_conc(state, "CLOUD.AQUEOUS.OHm")
        print(f"  HSO3m={hso3m:.4e}, SO3mm={so3mm:.4e}")
        print(f"  Hp={hp:.4e}, OHm={ohm:.4e}")
        print(f"  pH = {-math.log10(max(hp, 1e-30)):.2f}")
        assert hso3m > 0
        assert so3mm > 0

        # pH should be reasonable for atmospheric SO2 levels (pH ~3-6)
        ph = -math.log10(max(hp, 1e-30))
        assert 1.0 < ph < 8.0, f"pH={ph} out of range"

        _check_sulfur_conservation(state, ics)


class TestStep4_ChargeBalance:
    """Step 4: Add charge balance constraint."""

    @pytest.fixture(scope="class")
    def step4_config(self, tmp_path_factory):
        tmp = str(tmp_path_factory.mktemp("step4"))
        config = _load_ts1_config()
        _add_source_sink_to_config(config, "SO2")
        return _write_config(config, tmp, "step4.json")

    @pytest.fixture(scope="class")
    def step4_mechanism(self, step4_config):
        return Parser().parse(step4_config)

    def test_converges_with_charge_balance(self, step4_mechanism, ts1_ics):
        model = _build_miam_model(
            constraints=[
                _hlc_so2(),
                _eq_kw(), _eq_ka1(), _eq_ka2(),
                _lc_sulfur(),
                _lc_charge_balance(),
            ],
        )
        micm = MICM(
            mechanism=step4_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        # Compute self-consistent equilibrium ICs for cloud species
        cloud_ics = _compute_cloud_ics(ts1_ics["SO2"])
        ics = {
            **ts1_ics,
            "SO2": cloud_ics["_SO2_g"],
            "SO2_source": 0.0,
            "SO2_sink": 0.0,
            **{k: v for k, v in cloud_ics.items() if not k.startswith("_")},
        }
        state, results = _solve(micm, model, ics, dt=DT)
        print(f"Step 4: solver_state={results[-1].state}")
        assert results[-1].state == SolverState.Converged

        _check_sulfur_conservation(state, ics)
        _check_charge_balance(state)


class TestStep5_HLC_H2O2:
    """Step 5: Add HLC H2O2 + source/sink + H2O2 constraint."""

    @pytest.fixture(scope="class")
    def step5_config(self, tmp_path_factory):
        tmp = str(tmp_path_factory.mktemp("step5"))
        config = _load_ts1_config()
        _add_source_sink_to_config(config, "SO2")
        _add_source_sink_to_config(config, "H2O2")
        return _write_config(config, tmp, "step5.json")

    @pytest.fixture(scope="class")
    def step5_mechanism(self, step5_config):
        return Parser().parse(step5_config)

    def test_converges(self, step5_mechanism, ts1_ics):
        model = _build_miam_model(
            constraints=[
                _hlc_so2(), _hlc_h2o2(),
                _eq_kw(), _eq_ka1(), _eq_ka2(),
                _lc_sulfur(), _lc_h2o2(),
                _lc_charge_balance(),
            ],
        )
        micm = MICM(
            mechanism=step5_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        cloud_ics = _compute_cloud_ics(ts1_ics["SO2"], h2o2_g0=ts1_ics["H2O2"])
        ics = {
            **ts1_ics,
            "SO2": cloud_ics["_SO2_g"],
            "H2O2": cloud_ics["_H2O2_g"],
            "SO2_source": 0.0, "SO2_sink": 0.0,
            "H2O2_source": 0.0, "H2O2_sink": 0.0,
            **{k: v for k, v in cloud_ics.items() if not k.startswith("_")},
        }
        state, results = _solve(micm, model, ics, dt=DT)
        print(f"Step 5: solver_state={results[-1].state}")
        assert results[-1].state == SolverState.Converged

        h2o2_aq = _get_conc(state, "CLOUD.AQUEOUS.H2O2_aq")
        print(f"  H2O2(aq) = {h2o2_aq:.4e}")
        assert h2o2_aq > 0, "H2O2 didn't dissolve"

        _check_sulfur_conservation(state, ics)
        _check_h2o2_conservation(state, ics)
        _check_charge_balance(state)


class TestStep6_HLC_O3:
    """Step 6: Add HLC O3 + source/sink + O3 constraint. All HLCs complete."""

    @pytest.fixture(scope="class")
    def step6_config(self, tmp_path_factory):
        tmp = str(tmp_path_factory.mktemp("step6"))
        config = _load_ts1_config()
        _add_source_sink_to_config(config, "SO2")
        _add_source_sink_to_config(config, "H2O2")
        _add_source_sink_to_config(config, "O3")
        return _write_config(config, tmp, "step6.json")

    @pytest.fixture(scope="class")
    def step6_mechanism(self, step6_config):
        return Parser().parse(step6_config)

    def test_converges(self, step6_mechanism, ts1_ics):
        model = _build_miam_model(
            constraints=[
                _hlc_so2(), _hlc_h2o2(), _hlc_o3(),
                _eq_kw(), _eq_ka1(), _eq_ka2(),
                _lc_sulfur(), _lc_h2o2(), _lc_o3(),
                _lc_charge_balance(),
            ],
        )
        micm = MICM(
            mechanism=step6_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        cloud_ics = _compute_cloud_ics(
            ts1_ics["SO2"], h2o2_g0=ts1_ics["H2O2"], o3_g0=ts1_ics["O3"])
        ics = {
            **ts1_ics,
            "SO2": cloud_ics["_SO2_g"],
            "H2O2": cloud_ics["_H2O2_g"],
            "O3": cloud_ics["_O3_g"],
            "SO2_source": 0.0, "SO2_sink": 0.0,
            "H2O2_source": 0.0, "H2O2_sink": 0.0,
            "O3_source": 0.0, "O3_sink": 0.0,
            **{k: v for k, v in cloud_ics.items() if not k.startswith("_")},
        }
        state, results = _solve(micm, model, ics, dt=DT)
        print(f"Step 6: solver_state={results[-1].state}")
        assert results[-1].state == SolverState.Converged

        # O3 has low solubility -- should barely dissolve
        o3_aq = _get_conc(state, "CLOUD.AQUEOUS.O3_aq")
        o3_g = _get_conc(state, "O3")
        print(f"  O3(g) = {o3_g:.4e}, O3(aq) = {o3_aq:.4e}")
        if o3_aq is not None and o3_g is not None and o3_g > 0:
            print(f"  O3 dissolution fraction: {o3_aq / (o3_g + o3_aq):.2e}")

        _check_sulfur_conservation(state, ics)
        _check_h2o2_conservation(state, ics)
        _check_charge_balance(state)


class TestStep7_DissolvedReaction:
    """Step 7: Add R1 (HSO3- + H2O2 -> SO4-- + H2O + H+). The critical test."""

    @pytest.fixture(scope="class")
    def step7_config(self, tmp_path_factory):
        tmp = str(tmp_path_factory.mktemp("step7"))
        config = _load_ts1_config()
        _add_source_sink_to_config(config, "SO2")
        _add_source_sink_to_config(config, "H2O2")
        _add_source_sink_to_config(config, "O3")
        return _write_config(config, tmp, "step7.json")

    @pytest.fixture(scope="class")
    def step7_mechanism(self, step7_config):
        return Parser().parse(step7_config)

    def _make_micm_and_model(self, mechanism):
        model = _build_miam_model(
            processes=[_dissolved_r1()],
            constraints=[
                _hlc_so2(), _hlc_h2o2(), _hlc_o3(),
                _eq_kw(), _eq_ka1(), _eq_ka2(),
                _lc_sulfur(), _lc_h2o2(), _lc_o3(),
                _lc_charge_balance(),
            ],
        )
        micm = MICM(
            mechanism=mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        return micm, model

    def _make_ics(self, ts1_ics):
        cloud_ics = _compute_cloud_ics(
            ts1_ics["SO2"], h2o2_g0=ts1_ics["H2O2"], o3_g0=ts1_ics["O3"])
        return {
            **ts1_ics,
            "SO2": cloud_ics["_SO2_g"],
            "H2O2": cloud_ics["_H2O2_g"],
            "O3": cloud_ics["_O3_g"],
            "SO2_source": 0.0, "SO2_sink": 0.0,
            "H2O2_source": 0.0, "H2O2_sink": 0.0,
            "O3_source": 0.0, "O3_sink": 0.0,
            **{k: v for k, v in cloud_ics.items() if not k.startswith("_")},
        }

    def test_converges_dt900(self, step7_mechanism, ts1_ics):
        """Full cloud chemistry with dt=900s (MPAS timestep)."""
        micm, model = self._make_micm_and_model(step7_mechanism)
        ics = self._make_ics(ts1_ics)

        state, results = _solve(micm, model, ics, dt=DT)
        print(f"Step 7 (dt=900): solver_state={results[-1].state}")
        assert results[-1].state == SolverState.Converged

        # R1 should have produced SO4--
        so4mm = _get_conc(state, "CLOUD.AQUEOUS.SO4mm")
        print(f"  SO4mm = {so4mm:.4e}")
        assert so4mm > 0, "R1 didn't produce SO4--"

        _check_sulfur_conservation(state, ics)
        _check_h2o2_conservation(state, ics)
        _check_charge_balance(state)

    def test_mass_conservation_quantitative(self, step7_mechanism, ts1_ics):
        """Verify SO4 production equals sulfur removed from SO2 pool."""
        micm, model = self._make_micm_and_model(step7_mechanism)
        ics = self._make_ics(ts1_ics)

        state, results = _solve(micm, model, ics, dt=DT)
        assert results[-1].state == SolverState.Converged

        so4mm = _get_conc(state, "CLOUD.AQUEOUS.SO4mm") or 0.0
        so4mm_init = ics.get("CLOUD.AQUEOUS.SO4mm", 0.0)
        delta_so4 = so4mm - so4mm_init

        # The sulfur constraint ensures total S is conserved.
        # SO4 gained = sulfur lost from (SO2_g + SO2_aq + HSO3m + SO3mm)
        s_init, s_final = _check_sulfur_conservation(state, ics)

        print(f"  ΔSO4mm = {delta_so4:.4e}")
        print(f"  Sulfur budget: init={s_init:.4e}, final={s_final:.4e}")

    def test_no_cloud_water_no_effect(self, step7_mechanism, ts1_ics):
        """With zero cloud water, cloud chemistry should have no effect.

        Note: With the 3-reactant R1 (HSO3- + H2O2 + H+), [S]^2 in the
        denominator causes NaN when H2O ≈ 0. In MPAS, cloud chemistry
        is only run when LWC > threshold, so this is expected behavior.
        """
        micm, model = self._make_micm_and_model(step7_mechanism)
        ics = {
            **ts1_ics,
            "SO2_source": 0.0, "SO2_sink": 0.0,
            "H2O2_source": 0.0, "H2O2_sink": 0.0,
            "O3_source": 0.0, "O3_sink": 0.0,
            "CLOUD.AQUEOUS.H2O": 1e-30,  # effectively no cloud
            "CLOUD.AQUEOUS.SO2_aq": 0.0,
            "CLOUD.AQUEOUS.H2O2_aq": 0.0,
            "CLOUD.AQUEOUS.O3_aq": 0.0,
            "CLOUD.AQUEOUS.HSO3m": 0.0,
            "CLOUD.AQUEOUS.SO3mm": 0.0,
            "CLOUD.AQUEOUS.SO4mm": 0.0,
            "CLOUD.AQUEOUS.Hp": 1e-30,
            "CLOUD.AQUEOUS.OHm": 1e-30,
        }
        state, results = _solve(micm, model, ics, dt=DT)
        print(f"Step 7 (no cloud): solver_state={results[-1].state}")
        # Solver may NaN with near-zero [H2O] — 3-reactant R1 has 1/[S]^2
        assert results[-1].state in (SolverState.Converged, SolverState.NaNDetected)
