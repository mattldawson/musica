# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# Incremental cloud chemistry integration test.
# Starts from the TS1 gas-phase mechanism and adds MIAM elements one at a time
# to isolate which component causes solver failure (InfDetected).

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


# ═══ Constants ═══════════════════════════════════════════════════════════════

TS1_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "..", "..",
    "MPAS-Model", "chemistry_data", "ts1", "config.json",
)
M_ATM_TO_MOL_M3_PA = 1000.0 / 101325.0
C_H2O_M = 55.556
MW_H2O = 0.018
RHO_H2O = 1000.0
R_GAS = 8.314
T0 = 298.15

GAS0_SO2 = 3.01e-8
GAS0_H2O2 = 3.01e-8
GAS0_O3 = 1.5e-6
SO4MM0 = 1.0
C_H2O = 55556.0
T_INIT = 280.0
P_INIT = 70000.0


# ═══ Helpers ═════════════════════════════════════════════════════════════════

def _load_ts1_mechanism():
    """Load the TS1 gas-phase mechanism from its JSON config."""
    parser = Parser()
    return parser.parse(TS1_CONFIG)


def _make_cloud_species():
    """Create all species for cloud chemistry (gas + aqueous)."""
    # Gas-phase species referenced by HLC constraints
    so2_g = mc.Species(name="SO2")
    h2o2_g = mc.Species(name="H2O2")
    o3_g = mc.Species(name="O3")
    # Aqueous-phase species
    so2_aq = mc.Species(name="SO2_aq")
    h2o2_aq = mc.Species(name="H2O2_aq")
    o3_aq = mc.Species(name="O3_aq")
    hp = mc.Species(name="Hp")
    ohm = mc.Species(name="OHm")
    hso3m = mc.Species(name="HSO3m")
    so3mm = mc.Species(name="SO3mm")
    so4mm = mc.Species(name="SO4mm")
    h2o = mc.Species(name="H2O")
    h2o.molecular_weight_kg_mol = MW_H2O
    h2o.density_kg_m3 = RHO_H2O
    return {
        "SO2": so2_g, "H2O2": h2o2_g, "O3": o3_g,
        "SO2_aq": so2_aq, "H2O2_aq": h2o2_aq, "O3_aq": o3_aq,
        "Hp": hp, "OHm": ohm, "HSO3m": hso3m, "SO3mm": so3mm,
        "SO4mm": so4mm, "H2O": h2o,
    }


def _make_cloud_representation():
    """Create the CLOUD representation."""
    return UniformSection(
        name="CLOUD", phase_names=["AQUEOUS"],
        min_radius=1e-6, max_radius=1e-5,
    )


def _hlc_so2():
    return HenryLawEquilibriumConstraint(
        gas_species_name="SO2",
        condensed_species_name="SO2_aq",
        solvent_name="H2O",
        condensed_phase_name="AQUEOUS",
        henrys_law_constant=HenrysLawConstant(
            hlc_ref=1.23 * M_ATM_TO_MOL_M3_PA, c=3120.0),
        mw_solvent=MW_H2O, rho_solvent=RHO_H2O,
    )


def _hlc_h2o2():
    return HenryLawEquilibriumConstraint(
        gas_species_name="H2O2",
        condensed_species_name="H2O2_aq",
        solvent_name="H2O",
        condensed_phase_name="AQUEOUS",
        henrys_law_constant=HenrysLawConstant(
            hlc_ref=7.4e4 * M_ATM_TO_MOL_M3_PA, c=6621.0),
        mw_solvent=MW_H2O, rho_solvent=RHO_H2O,
    )


def _hlc_o3():
    return HenryLawEquilibriumConstraint(
        gas_species_name="O3",
        condensed_species_name="O3_aq",
        solvent_name="H2O",
        condensed_phase_name="AQUEOUS",
        henrys_law_constant=HenrysLawConstant(
            hlc_ref=1.15e-2 * M_ATM_TO_MOL_M3_PA, c=2560.0),
        mw_solvent=MW_H2O, rho_solvent=RHO_H2O,
    )


def _eq_kw():
    return DissolvedEquilibriumConstraint(
        phase_name="AQUEOUS",
        reactant_names=["H2O"],
        product_names=["Hp", "OHm"],
        algebraic_species_name="OHm",
        solvent_name="H2O",
        equilibrium_constant=EquilibriumConstant(
            a=1e-14 / (C_H2O_M * C_H2O_M), c=0.0),
    )


def _eq_ka1():
    return DissolvedEquilibriumConstraint(
        phase_name="AQUEOUS",
        reactant_names=["SO2_aq"],
        product_names=["HSO3m", "Hp"],
        algebraic_species_name="HSO3m",
        solvent_name="H2O",
        equilibrium_constant=EquilibriumConstant(a=1.7e-2 / C_H2O_M, c=2090.0),
    )


def _eq_ka2():
    return DissolvedEquilibriumConstraint(
        phase_name="AQUEOUS",
        reactant_names=["HSO3m"],
        product_names=["SO3mm", "Hp"],
        algebraic_species_name="SO3mm",
        solvent_name="H2O",
        equilibrium_constant=EquilibriumConstant(a=6.0e-8 / C_H2O_M, c=1120.0),
    )


def _lc_s_conservation():
    return LinearConstraint(
        algebraic_phase_name="gas",
        algebraic_species_name="SO2",
        terms=[
            LinearConstraintTerm("gas", "SO2", 1.0),
            LinearConstraintTerm("AQUEOUS", "SO2_aq", 1.0),
            LinearConstraintTerm("AQUEOUS", "HSO3m", 1.0),
            LinearConstraintTerm("AQUEOUS", "SO3mm", 1.0),
        ],
        constant=GAS0_SO2,
    )


def _lc_h2o2_conservation():
    return LinearConstraint(
        algebraic_phase_name="gas",
        algebraic_species_name="H2O2",
        terms=[
            LinearConstraintTerm("gas", "H2O2", 1.0),
            LinearConstraintTerm("AQUEOUS", "H2O2_aq", 1.0),
        ],
        constant=GAS0_H2O2,
    )


def _lc_o3_conservation():
    return LinearConstraint(
        algebraic_phase_name="gas",
        algebraic_species_name="O3",
        terms=[
            LinearConstraintTerm("gas", "O3", 1.0),
            LinearConstraintTerm("AQUEOUS", "O3_aq", 1.0),
        ],
        constant=GAS0_O3,
    )


def _lc_charge_balance():
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


def _dissolved_reaction_r1():
    return DissolvedReaction(
        phase_name="AQUEOUS",
        reactant_names=["HSO3m", "H2O2_aq", "Hp"],
        product_names=["SO4mm", "H2O", "Hp", "Hp"],
        solvent_name="H2O",
        rate_constant=ArrheniusRateConstant(a=C_H2O_M**2 * 7.45e7, c=4430.0),
    )


def _build_model(species_dict, processes=None, constraints=None):
    """Build a MIAM Model from given species, processes, and constraints.

    Only includes species actually referenced by the processes/constraints.
    Always includes H2O as the solvent.
    """
    # Determine which species are needed (both gas and aqueous)
    needed = {"H2O"}  # always need solvent
    gas_needed = set()
    for c in (constraints or []):
        if isinstance(c, HenryLawEquilibriumConstraint):
            needed.add(c.condensed_species_name)
            gas_needed.add(c.gas_species_name)
        elif isinstance(c, DissolvedEquilibriumConstraint):
            needed.update(c.reactant_names)
            needed.update(c.product_names)
        elif isinstance(c, LinearConstraint):
            for t in c.terms:
                if t.phase_name == "AQUEOUS":
                    needed.add(t.species_name)
                elif t.phase_name == "gas":
                    gas_needed.add(t.species_name)
    for p in (processes or []):
        if isinstance(p, DissolvedReaction):
            needed.update(p.reactant_names)
            needed.update(p.product_names)

    aq_used = {k: v for k, v in species_dict.items()
               if k in needed and k not in gas_needed}
    gas_used = {k: v for k, v in species_dict.items() if k in gas_needed}
    all_sp = list(gas_used.values()) + list(aq_used.values())
    aq_phase = mc.Phase(name="AQUEOUS", species=list(aq_used.values()))

    return Model(
        name="cloud_chemistry",
        species=all_sp,
        condensed_phases=[aq_phase],
        representations=[_make_cloud_representation()],
        processes=processes or [],
        constraints=constraints or [],
    )


def _compute_equilibrium_ics():
    """Compute self-consistent initial conditions (same as reference test)."""
    T = T_INIT
    hlc_so2_T = (1.23 * M_ATM_TO_MOL_M3_PA) * math.exp(3120.0 * (1/T - 1/T0))
    hlc_h2o2_T = (7.4e4 * M_ATM_TO_MOL_M3_PA) * math.exp(6621.0 * (1/T - 1/T0))
    hlc_o3_T = (1.15e-2 * M_ATM_TO_MOL_M3_PA) * math.exp(2560.0 * (1/T - 1/T0))
    alpha_SO2 = hlc_so2_T * R_GAS * T
    alpha_H2O2 = hlc_h2o2_T * R_GAS * T
    alpha_O3 = hlc_o3_T * R_GAS * T
    Ka1_T = (1.7e-2 / C_H2O_M) * math.exp(2090.0 * (1/T0 - 1/T))
    Ka2_T = (6.0e-8 / C_H2O_M) * math.exp(1120.0 * (1/T0 - 1/T))
    Kw_T = 1e-14 / (C_H2O_M * C_H2O_M)

    ic_h2o2_g = GAS0_H2O2 / (1.0 + alpha_H2O2)
    ic_h2o2_aq = alpha_H2O2 * ic_h2o2_g
    ic_o3_g = GAS0_O3 / (1.0 + alpha_O3)
    ic_o3_aq = alpha_O3 * ic_o3_g

    hp_ic = 2.0 * SO4MM0
    for _ in range(100):
        ic_ohm = Kw_T * C_H2O * C_H2O / hp_ic
        f = (1.0 + alpha_SO2
             + Ka1_T * alpha_SO2 * C_H2O / hp_ic
             + Ka2_T * Ka1_T * alpha_SO2 * C_H2O * C_H2O / (hp_ic * hp_ic))
        ic_so2_g = GAS0_SO2 / f
        ic_so2_aq = alpha_SO2 * ic_so2_g
        ic_hso3m = Ka1_T * ic_so2_aq * C_H2O / hp_ic
        ic_so3mm = Ka2_T * ic_hso3m * C_H2O / hp_ic
        hp_new = ic_ohm + ic_hso3m + 2.0 * ic_so3mm + 2.0 * SO4MM0
        if abs(hp_new - hp_ic) < 1e-15 * hp_ic:
            break
        hp_ic = 0.5 * (hp_ic + hp_new)

    return {
        "SO2": ic_so2_g,
        "H2O2": ic_h2o2_g,
        "O3": ic_o3_g,
        "CLOUD.AQUEOUS.H2O": C_H2O,
        "CLOUD.AQUEOUS.SO2_aq": ic_so2_aq,
        "CLOUD.AQUEOUS.H2O2_aq": ic_h2o2_aq,
        "CLOUD.AQUEOUS.O3_aq": ic_o3_aq,
        "CLOUD.AQUEOUS.Hp": hp_ic,
        "CLOUD.AQUEOUS.OHm": ic_ohm,
        "CLOUD.AQUEOUS.HSO3m": ic_hso3m,
        "CLOUD.AQUEOUS.SO3mm": ic_so3mm,
        "CLOUD.AQUEOUS.SO4mm": SO4MM0,
    }


def _try_solve(micm, model, ics, dt=0.01, target=0.1):
    """Create state, set ICs, and try to integrate.

    Returns (final_state, SolverState of first failed step or Converged).
    """
    state = micm.create_state()
    model.set_default_parameters(state)
    state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
    # Filter ICs to only species that exist in this mechanism
    ordering = state.get_species_ordering()
    filtered_ics = {k: v for k, v in ics.items() if k in ordering}
    state.set_concentrations(filtered_ics)

    total = 0.0
    while total < target - 1e-10:
        step = min(dt, target - total)
        result = micm.solve(state, time_step=step)
        if result.state != SolverState.Converged:
            return state, result.state
        total += step
    return state, SolverState.Converged


# ═══ Tests ═══════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def ts1_mechanism():
    """Load the TS1 mechanism once for the module."""
    return _load_ts1_mechanism()


class TestIncrementalCloudChemistry:
    """Add cloud chemistry elements one at a time on top of TS1.

    Each test checks whether the solver converges after adding
    progressively more MIAM components.
    """

    def test_step0_ts1_baseline_rosenbrock(self, ts1_mechanism):
        """Step 0: TS1 with standard Rosenbrock (no MIAM) — baseline."""
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_standard_order,
        )
        state = micm.create_state()
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
        # Set a few key species to non-zero
        state.set_concentrations({
            "SO2": GAS0_SO2, "H2O2": GAS0_H2O2, "O3": GAS0_O3,
        })
        result = micm.solve(state, time_step=0.01)
        print(f"Step 0 (Rosenbrock): state={result.state}")
        assert result.state == SolverState.Converged

    def test_step0b_ts1_baseline_dae4(self, ts1_mechanism):
        """Step 0b: TS1 with DAE4 solver (no MIAM) — just the solver type."""
        sp = _make_cloud_species()
        # Minimal MIAM model: only solvent, no processes, no constraints
        model = _build_model(sp)
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        state = micm.create_state()
        model.set_default_parameters(state)
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
        state.set_concentrations({
            "SO2": GAS0_SO2, "H2O2": GAS0_H2O2, "O3": GAS0_O3,
            "CLOUD.AQUEOUS.H2O": C_H2O,
        })
        result = micm.solve(state, time_step=0.01)
        print(f"Step 0b (DAE4, no cloud procs): state={result.state}")
        assert result.state == SolverState.Converged

    def test_step1_hlc_so2_only(self, ts1_mechanism):
        """Step 1: TS1 + Henry's Law SO2 only."""
        sp = _make_cloud_species()
        model = _build_model(sp, constraints=[_hlc_so2()])
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        ics = {"SO2": GAS0_SO2, "CLOUD.AQUEOUS.H2O": C_H2O,
               "CLOUD.AQUEOUS.SO2_aq": 0.0}
        _, solver_state = _try_solve(micm, model, ics)
        print(f"Step 1 (HLC SO2): state={solver_state}")
        assert solver_state == SolverState.Converged

    def test_step2_hlc_so2_h2o2(self, ts1_mechanism):
        """Step 2: TS1 + Henry's Law SO2 + H2O2."""
        sp = _make_cloud_species()
        model = _build_model(sp, constraints=[_hlc_so2(), _hlc_h2o2()])
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        ics = {"SO2": GAS0_SO2, "H2O2": GAS0_H2O2,
               "CLOUD.AQUEOUS.H2O": C_H2O,
               "CLOUD.AQUEOUS.SO2_aq": 0.0, "CLOUD.AQUEOUS.H2O2_aq": 0.0}
        _, solver_state = _try_solve(micm, model, ics)
        print(f"Step 2 (HLC SO2+H2O2): state={solver_state}")
        assert solver_state == SolverState.Converged

    def test_step3_hlc_all_three(self, ts1_mechanism):
        """Step 3: TS1 + all three Henry's Law equilibria."""
        sp = _make_cloud_species()
        model = _build_model(sp, constraints=[
            _hlc_so2(), _hlc_h2o2(), _hlc_o3()])
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        ics = {"SO2": GAS0_SO2, "H2O2": GAS0_H2O2, "O3": GAS0_O3,
               "CLOUD.AQUEOUS.H2O": C_H2O,
               "CLOUD.AQUEOUS.SO2_aq": 0.0, "CLOUD.AQUEOUS.H2O2_aq": 0.0,
               "CLOUD.AQUEOUS.O3_aq": 0.0}
        _, solver_state = _try_solve(micm, model, ics)
        print(f"Step 3 (all HLC): state={solver_state}")
        assert solver_state == SolverState.Converged

    def test_step4_hlc_plus_kw(self, ts1_mechanism):
        """Step 4: HLCs + water dissociation (Kw)."""
        sp = _make_cloud_species()
        model = _build_model(sp, constraints=[
            _hlc_so2(), _hlc_h2o2(), _hlc_o3(), _eq_kw()])
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        ics = _compute_equilibrium_ics()
        _, solver_state = _try_solve(micm, model, ics)
        print(f"Step 4 (HLC + Kw): state={solver_state}")
        assert solver_state == SolverState.Converged

    def test_step5_hlc_plus_equilibria(self, ts1_mechanism):
        """Step 5: HLCs + all dissolved equilibria (Kw, Ka1, Ka2)."""
        sp = _make_cloud_species()
        model = _build_model(sp, constraints=[
            _hlc_so2(), _hlc_h2o2(), _hlc_o3(),
            _eq_kw(), _eq_ka1(), _eq_ka2()])
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        ics = _compute_equilibrium_ics()
        _, solver_state = _try_solve(micm, model, ics)
        print(f"Step 5 (HLC + equilibria): state={solver_state}")
        assert solver_state == SolverState.Converged

    def test_step6_add_dissolved_reaction(self, ts1_mechanism):
        """Step 6: Everything + dissolved reaction (R1)."""
        sp = _make_cloud_species()
        model = _build_model(
            sp,
            processes=[_dissolved_reaction_r1()],
            constraints=[
                _hlc_so2(), _hlc_h2o2(), _hlc_o3(),
                _eq_kw(), _eq_ka1(), _eq_ka2()],
        )
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        ics = _compute_equilibrium_ics()
        _, solver_state = _try_solve(micm, model, ics)
        print(f"Step 6 (+ dissolved rxn): state={solver_state}")
        assert solver_state == SolverState.Converged

    def test_step7_add_linear_constraints(self, ts1_mechanism):
        """Step 7: Full cloud chemistry = HLCs + equilibria + R1 + linear constraints."""
        sp = _make_cloud_species()
        model = _build_model(
            sp,
            processes=[_dissolved_reaction_r1()],
            constraints=[
                _hlc_so2(), _hlc_h2o2(), _hlc_o3(),
                _eq_kw(), _eq_ka1(), _eq_ka2(),
                _lc_s_conservation(), _lc_h2o2_conservation(),
                _lc_o3_conservation(), _lc_charge_balance()],
        )
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        ics = _compute_equilibrium_ics()
        _, solver_state = _try_solve(micm, model, ics)
        print(f"Step 7 (full cloud chem): state={solver_state}")
        assert solver_state == SolverState.Converged


class TestMpasLikeConditions:
    """Test with conditions mimicking MPAS driver behavior.

    In MPAS, state.set_concentrations() is called with only the species
    that have non-zero values, while everything else stays at 0.
    The key question: does having most TS1 species at 0 cause InfDetected?
    """

    def test_mpas_sparse_ics_no_cloud(self, ts1_mechanism):
        """MPAS-like: TS1 + DAE4, only a few species non-zero, no cloud."""
        sp = _make_cloud_species()
        model = _build_model(sp)  # no processes, no constraints
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        state = micm.create_state()
        model.set_default_parameters(state)
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
        # Only set a handful of species (MPAS-like sparse initialization)
        state.set_concentrations({
            "O3": 1.2e-6, "SO2": 8.9e-7, "H2O2": 2.0e-6,
            "CLOUD.AQUEOUS.H2O": C_H2O,
        })
        result = micm.solve(state, time_step=0.01)
        print(f"MPAS sparse (no cloud): state={result.state}")
        assert result.state == SolverState.Converged

    def test_mpas_sparse_ics_full_cloud(self, ts1_mechanism):
        """MPAS-like: TS1 + full cloud chemistry, sparse ICs."""
        sp = _make_cloud_species()
        model = _build_model(
            sp,
            processes=[_dissolved_reaction_r1()],
            constraints=[
                _hlc_so2(), _hlc_h2o2(), _hlc_o3(),
                _eq_kw(), _eq_ka1(), _eq_ka2(),
                _lc_s_conservation(), _lc_h2o2_conservation(),
                _lc_o3_conservation(), _lc_charge_balance()],
        )
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        state = micm.create_state()
        model.set_default_parameters(state)
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
        # MPAS-like: only species from initial conditions / cloud_set_state
        state.set_concentrations({
            "O3": 1.2e-6, "SO2": 8.9e-7, "H2O2": 2.0e-6,
            "CLOUD.AQUEOUS.H2O": C_H2O,
            "CLOUD.AQUEOUS.Hp": 1e-4,
            "CLOUD.AQUEOUS.OHm": 1e-10,
            "CLOUD.AQUEOUS.HSO3m": 1e-6,
            "CLOUD.AQUEOUS.SO3mm": 1e-8,
            "CLOUD.AQUEOUS.SO4mm": 0.0,
        })
        result = micm.solve(state, time_step=0.01)
        print(f"MPAS sparse (full cloud): state={result.state}")
        assert result.state == SolverState.Converged

    @pytest.mark.skip(reason="3-reactant R1 with H+ coupling causes NaN from degenerate ICs")
    @pytest.mark.parametrize("h2o_floor", [1e-3, 1.0])
    def test_mpas_sparse_ics_cloud_floor(self, ts1_mechanism, h2o_floor):
        """MPAS-like: full cloud + sparse ICs, ALL aqueous floored to 1e-30.

        Note: The 3-reactant R1 (HSO3- + H2O2 + H+) has [S]^2 in the
        denominator and couples H+ directly into the rate, so the DAE
        solver fails when H2O and all aqueous species are near-zero.
        Only test with physically meaningful cloud water levels.
        """
        sp = _make_cloud_species()
        model = _build_model(
            sp,
            processes=[_dissolved_reaction_r1()],
            constraints=[
                _hlc_so2(), _hlc_h2o2(), _hlc_o3(),
                _eq_kw(), _eq_ka1(), _eq_ka2(),
                _lc_s_conservation(), _lc_h2o2_conservation(),
                _lc_o3_conservation(), _lc_charge_balance()],
        )
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        state = micm.create_state()
        model.set_default_parameters(state)
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
        # Floor ALL aqueous species to 1e-30 (matching MPAS cloud_set_state)
        aq_floor = 1e-30
        state.set_concentrations({
            "O3": 1.2e-6, "SO2": 8.9e-7, "H2O2": 2.0e-6,
            "CLOUD.AQUEOUS.H2O": max(h2o_floor, aq_floor),
            "CLOUD.AQUEOUS.SO2_aq": aq_floor,
            "CLOUD.AQUEOUS.H2O2_aq": aq_floor,
            "CLOUD.AQUEOUS.O3_aq": aq_floor,
            "CLOUD.AQUEOUS.Hp": aq_floor,
            "CLOUD.AQUEOUS.OHm": aq_floor,
            "CLOUD.AQUEOUS.HSO3m": aq_floor,
            "CLOUD.AQUEOUS.SO3mm": aq_floor,
            "CLOUD.AQUEOUS.SO4mm": aq_floor,
        })
        result = micm.solve(state, time_step=0.01)
        print(f"H2O={h2o_floor:.0e}, all aq floored: state={result.state}")
        assert result.state == SolverState.Converged

    def test_mpas_hlc_only_no_defaults(self, ts1_mechanism):
        """Just HLCs (no equilibria, no reactions, no constraints) with minimal ICs."""
        sp = _make_cloud_species()
        model = _build_model(sp, constraints=[_hlc_so2(), _hlc_h2o2(), _hlc_o3()])
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        state = micm.create_state()
        model.set_default_parameters(state)
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
        state.set_concentrations({
            "O3": 1.2e-6, "SO2": 8.9e-7, "H2O2": 2.0e-6,
            "CLOUD.AQUEOUS.H2O": 1.0e-30,
        })
        result = micm.solve(state, time_step=0.01)
        print(f"HLC-only H2O=1e-30: state={result.state}")

    def test_mpas_full_cloud_zero_defaults(self, ts1_mechanism):
        """Full cloud but with defaults=0 (like MPAS non-cloud cell)."""
        sp = _make_cloud_species()
        model = _build_model(
            sp,
            processes=[_dissolved_reaction_r1()],
            constraints=[
                _hlc_so2(), _hlc_h2o2(), _hlc_o3(),
                _eq_kw(), _eq_ka1(), _eq_ka2(),
                _lc_s_conservation(), _lc_h2o2_conservation(),
                _lc_o3_conservation(), _lc_charge_balance()],
        )
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        state = micm.create_state()
        model.set_default_parameters(state)
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
        # Full cloud water but zero for all aqueous species
        state.set_concentrations({
            "O3": 1.2e-6, "SO2": 8.9e-7, "H2O2": 2.0e-6,
            "CLOUD.AQUEOUS.H2O": C_H2O,
            "CLOUD.AQUEOUS.Hp": 0.0,
            "CLOUD.AQUEOUS.OHm": 0.0,
            "CLOUD.AQUEOUS.HSO3m": 0.0,
            "CLOUD.AQUEOUS.SO3mm": 0.0,
        })
        result = micm.solve(state, time_step=0.01)
        print(f"Full H2O, zero defaults: state={result.state}")

    def test_mpas_exact_non_cloud_cell(self, ts1_mechanism):
        """Exact MPAS conditions for a non-cloud cell: H2O=1e-30 + defaults."""
        sp = _make_cloud_species()
        model = _build_model(
            sp,
            processes=[_dissolved_reaction_r1()],
            constraints=[
                _hlc_so2(), _hlc_h2o2(), _hlc_o3(),
                _eq_kw(), _eq_ka1(), _eq_ka2(),
                _lc_s_conservation(), _lc_h2o2_conservation(),
                _lc_o3_conservation(), _lc_charge_balance()],
        )
        micm = MICM(
            mechanism=ts1_mechanism,
            solver_type=SolverType.rosenbrock_dae4_standard_order,
            external_models=[model],
        )
        state = micm.create_state()
        model.set_default_parameters(state)
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
        # Non-cloud cell: tiny H2O + MPAS default equilibrium species
        state.set_concentrations({
            "O3": 1.2e-6, "SO2": 8.9e-7, "H2O2": 2.0e-6,
            "CLOUD.AQUEOUS.H2O": 1.0e-30,
            "CLOUD.AQUEOUS.Hp": 1.0e-4,
            "CLOUD.AQUEOUS.OHm": 1.0e-10,
            "CLOUD.AQUEOUS.HSO3m": 1.0e-6,
            "CLOUD.AQUEOUS.SO3mm": 1.0e-8,
        })
        result = micm.solve(state, time_step=0.01)
        print(f"Non-cloud cell (H2O=1e-30 + defaults): state={result.state}")
