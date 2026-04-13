# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# Integration test for MIAM aerosol model created from a config file.
# Validates the config → parser → builder → solver path.

import math
import pytest

from musica.micm import MICM, SolverType, SolverState
from musica.utils import find_config_path


# ═══ Constants (matching test_miam_cloud_chemistry.py) ═══════════════════════

M_ATM_TO_MOL_M3_PA = 1000.0 / 101325.0
C_H2O_M = 55.556         # mol/L
MW_H2O = 0.018            # kg/mol
RHO_H2O = 1000.0          # kg/m3
R_GAS = 8.314             # J/(mol·K)
T0 = 298.15               # K

# Initial conditions
GAS0_SO2 = 3.01e-8        # mol/m3 (~ 1 ppb)
GAS0_H2O2 = 3.01e-8       # mol/m3
GAS0_O3 = 1.5e-6          # mol/m3
SO4MM0 = 1.0              # mol/m3 (background sulfate)
C_H2O = 55556.0           # mol/m3 (liquid water in droplet)
T_INIT = 280.0            # K
P_INIT = 70000.0          # Pa


def _config_path():
    return find_config_path("v1", "cam_cloud_chemistry", "config.json")


def _compute_equilibrium_ics():
    """Compute self-consistent initial conditions for the config-driven model.

    The config uses unified species names (SO2 in both gas and AQUEOUS),
    so state variable names are e.g. CLOUD.AQUEOUS.SO2 (not CLOUD.AQUEOUS.SO2_aq).
    """
    T = T_INIT

    # Temperature-adjusted Henry's Law constants
    # Config HLC_ref values are already in SI [mol m-3 Pa-1]
    # (the config stores values that were multiplied by M_ATM_TO_MOL_M3_PA)
    hlc_so2_T = 1.214e-2 * math.exp(3120.0 * (1.0 / T - 1.0 / T0))
    hlc_h2o2_T = 7.306e2 * math.exp(6621.0 * (1.0 / T - 1.0 / T0))
    hlc_o3_T = 1.135e-4 * math.exp(2560.0 * (1.0 / T - 1.0 / T0))
    alpha_SO2 = hlc_so2_T * R_GAS * T
    alpha_H2O2 = hlc_h2o2_T * R_GAS * T
    alpha_O3 = hlc_o3_T * R_GAS * T

    # Temperature-adjusted equilibrium constants
    Ka1_T = 3.06e-4 * math.exp(2090.0 * (1.0 / T0 - 1.0 / T))
    Ka2_T = 1.08e-9 * math.exp(1120.0 * (1.0 / T0 - 1.0 / T))
    Kw_T = 3.24e-18

    # H2O2 and O3: simple HLC split
    ic_h2o2_g = GAS0_H2O2 / (1.0 + alpha_H2O2)
    ic_h2o2_aq = alpha_H2O2 * ic_h2o2_g
    ic_o3_g = GAS0_O3 / (1.0 + alpha_O3)
    ic_o3_aq = alpha_O3 * ic_o3_g

    # Iterate on [H+] for SO2 equilibria + charge balance
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
        "CLOUD.AQUEOUS.SO2": ic_so2_aq,
        "CLOUD.AQUEOUS.H2O2": ic_h2o2_aq,
        "CLOUD.AQUEOUS.O3": ic_o3_aq,
        "CLOUD.AQUEOUS.Hp": hp_ic,
        "CLOUD.AQUEOUS.OHm": ic_ohm,
        "CLOUD.AQUEOUS.HSO3m": ic_hso3m,
        "CLOUD.AQUEOUS.SO3mm": ic_so3mm,
        "CLOUD.AQUEOUS.SO4mm": SO4MM0,
    }


def _integrate(micm, state, target_time, dt_init=0.01):
    """Adaptive time stepping loop."""
    total_time = 0.0
    dt = dt_init
    while total_time < target_time - 1e-10:
        step = min(dt, target_time - total_time)
        result = micm.solve(state, time_step=step)
        assert result.state == SolverState.Converged, \
            f"Solver failed at t={total_time:.4f}s"
        total_time += step
        if total_time > 0.1 and dt < 0.1:
            dt = 0.1
        if total_time > 1.0 and dt < 1.0:
            dt = 1.0


def _get_conc(state, name):
    """Get scalar concentration from state."""
    return state.get_concentrations()[name][0]


# ═══ Tests ═══════════════════════════════════════════════════════════════════

class TestConfigDrivenSolverCreation:
    """Test creating a MIAM-enabled solver from a config file."""

    def test_create_solver_dae4(self):
        micm = MICM(
            config_path=_config_path(),
            solver_type=SolverType.rosenbrock_dae4_standard_order,
        )
        assert micm is not None
        assert micm.solver_type() == SolverType.rosenbrock_dae4_standard_order

    def test_create_solver_dae6(self):
        micm = MICM(
            config_path=_config_path(),
            solver_type=SolverType.rosenbrock_dae6_standard_order,
        )
        assert micm is not None

    def test_create_state(self):
        micm = MICM(
            config_path=_config_path(),
            solver_type=SolverType.rosenbrock_dae4_standard_order,
        )
        state = micm.create_state()
        assert state is not None

    def test_state_has_cloud_species(self):
        """The state should contain both gas and condensed-phase species."""
        micm = MICM(
            config_path=_config_path(),
            solver_type=SolverType.rosenbrock_dae4_standard_order,
        )
        state = micm.create_state()
        concs = state.get_concentrations()
        # Gas phase species
        assert "SO2" in concs
        assert "H2O2" in concs
        assert "O3" in concs
        # Condensed-phase species (prefixed with representation.phase.)
        assert "CLOUD.AQUEOUS.SO2" in concs
        assert "CLOUD.AQUEOUS.H2O2" in concs
        assert "CLOUD.AQUEOUS.O3" in concs
        assert "CLOUD.AQUEOUS.Hp" in concs
        assert "CLOUD.AQUEOUS.OHm" in concs
        assert "CLOUD.AQUEOUS.HSO3m" in concs
        assert "CLOUD.AQUEOUS.SO3mm" in concs
        assert "CLOUD.AQUEOUS.SO4mm" in concs
        assert "CLOUD.AQUEOUS.H2O" in concs


class TestConfigDrivenSolve:
    """Test solving with a config-driven MIAM solver."""

    def test_solve_converges(self):
        """Solver should converge and produce SO4."""
        micm = MICM(
            config_path=_config_path(),
            solver_type=SolverType.rosenbrock_dae4_standard_order,
        )
        state = micm.create_state()

        # Set conditions
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)

        # Set self-consistent initial conditions
        ics = _compute_equilibrium_ics()
        state.set_concentrations(ics)

        # Integrate for 10 seconds
        _integrate(micm, state, target_time=10.0)

        # SO4 should increase from kinetics
        so4_f = _get_conc(state, "CLOUD.AQUEOUS.SO4mm")
        assert so4_f >= SO4MM0, "SO4 should only increase from kinetics"

    def test_mass_conservation_h2o2(self):
        """Total H2O2 (gas + aq) should decrease as it's consumed."""
        micm = MICM(
            config_path=_config_path(),
            solver_type=SolverType.rosenbrock_dae4_standard_order,
        )
        state = micm.create_state()
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
        ics = _compute_equilibrium_ics()
        state.set_concentrations(ics)

        initial_total_h2o2 = ics["H2O2"] + ics["CLOUD.AQUEOUS.H2O2"]

        _integrate(micm, state, target_time=10.0)

        # H2O2 is consumed by the dissolved reaction (HSO3- + H2O2 → products)
        final_h2o2_g = _get_conc(state, "H2O2")
        final_h2o2_aq = _get_conc(state, "CLOUD.AQUEOUS.H2O2")
        final_total_h2o2 = final_h2o2_g + final_h2o2_aq

        # With diagnose-from-state, the constraint is re-diagnosed each step.
        # The total should be <= initial (H2O2 consumed by reaction).
        assert final_total_h2o2 <= initial_total_h2o2 + 1e-15

    def test_charge_balance(self):
        """Charge balance should be maintained at all times."""
        micm = MICM(
            config_path=_config_path(),
            solver_type=SolverType.rosenbrock_dae4_standard_order,
        )
        state = micm.create_state()
        state.set_conditions(temperatures=T_INIT, pressures=P_INIT)
        ics = _compute_equilibrium_ics()
        state.set_concentrations(ics)

        _integrate(micm, state, target_time=10.0)

        # Check: H+ - OH- - HSO3- - 2*SO3-- - 2*SO4-- ≈ 0
        hp = _get_conc(state, "CLOUD.AQUEOUS.Hp")
        ohm = _get_conc(state, "CLOUD.AQUEOUS.OHm")
        hso3m = _get_conc(state, "CLOUD.AQUEOUS.HSO3m")
        so3mm = _get_conc(state, "CLOUD.AQUEOUS.SO3mm")
        so4mm = _get_conc(state, "CLOUD.AQUEOUS.SO4mm")

        charge_imbalance = hp - ohm - hso3m - 2.0 * so3mm - 2.0 * so4mm
        assert abs(charge_imbalance) < 1e-6, \
            f"Charge balance violated: imbalance = {charge_imbalance}"


class TestConfigDrivenMechanismRoundTrip:
    """Test that a programmatic mechanism exported to config produces same results."""

    def test_export_and_reload(self, tmp_path):
        """Export a gas-phase-only mechanism and reload it."""
        import musica.mechanism_configuration as mc

        so2 = mc.Species(name="SO2")
        h2o2 = mc.Species(name="H2O2")
        o3 = mc.Species(name="O3")
        gas = mc.Phase(name="gas", species=[so2, h2o2, o3])
        mechanism = mc.Mechanism(
            species=[so2, h2o2, o3], phases=[gas], reactions=[]
        )

        # Export to JSON
        config_file = tmp_path / "test_mechanism.json"
        mechanism.export(str(config_file))

        # Reload from config
        micm_from_config = MICM(
            config_path=str(config_file),
            solver_type=SolverType.rosenbrock_standard_order,
        )
        assert micm_from_config is not None

        # Create from mechanism directly
        micm_from_mech = MICM(
            mechanism=mechanism,
            solver_type=SolverType.rosenbrock_standard_order,
        )

        # Both should create identical state structure
        state_config = micm_from_config.create_state()
        state_mech = micm_from_mech.create_state()
        concs_config = state_config.get_concentrations()
        concs_mech = state_mech.get_concentrations()
        assert set(concs_config.keys()) == set(concs_mech.keys())
