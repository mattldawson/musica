# Fix: Cloud Water Concentration Units (H2O = mol/m³ of AIR, not solution)

## Problem

Several test files and Tutorial 14 set `CLOUD.AQUEOUS.H2O` to **55556 mol/m³** (the molarity of pure liquid water). This is wrong — MICM/MIAM state variables are all in **mol/m³ of air**, so the H2O concentration should reflect the cloud liquid water content (LWC) per cubic meter of air:

```
[H2O]_air = LWC [kg/m³] / MW_H2O [kg/mol]
```

A typical cloud has LWC ≈ 0.3 g/m³, giving `[H2O] ≈ 0.017 mol/m³-air`.

Using 55556 implies a volume fraction `f_v = [H2O] × MW / ρ ≈ 1.0` — the entire cubic meter of "air" is liquid water, which is nonphysical.

## Impact

The `[solvent]` appears in three MIAM formulas:

| Component | Formula | Effect of wrong H2O |
|---|---|---|
| **HenryLawEquilibriumConstraint** | `[X_aq] = HLC·R·T·f_v·[X_g]` where `f_v = [H2O]·MW/ρ` | Dissolution ~3.3M× too high; nearly all gas dissolves |
| **DissolvedEquilibriumConstraint** | `Keq·∏[R]/[S]^(nr-1) = ∏[P]/[S]^(np-1)` | Equilibrium ratios shift by powers of 3.3M |
| **DissolvedReaction** | `rate = k/[S]^(nr-1) · ∏[Ri]` | R1 rate ~3.3M× too slow (divides by huge solvent) |

The MPAS driver (`mpas_chemistry_cloud.F90`) correctly computes `cloud_conc = prescribed_lwc * rho_d / cloud_water_mw`, so the deployed code is fine. Only tests are affected.

## Separately correct: `C_H2O_M = 55.556 mol/L` in rate constants

The value `C_H2O_M * 7.45e7` used in `ArrheniusRateConstant(a=...)` is a **unit conversion factor** (molarity of pure water converts M⁻¹s⁻¹ literature values into MIAM's internal units). This is correct and should NOT be changed.

## Files to fix

All of these set `C_H2O = 55556.0` as the `CLOUD.AQUEOUS.H2O` state variable:

| File | Line | Notes |
|---|---|---|
| `python/test/integration/test_miam_cloud_chemistry.py` | 40 | `C_H2O = 55556.0` used as IC |
| `python/test/integration/test_miam_from_config.py` | 28 | Same |
| `python/test/integration/test_incremental_cloud.py` | 46 | Same |
| `python/test/integration/test_ts1_cloud_incremental.py` | 49 | Same |
| `fortran/test/integration/test_miam_cloud_chemistry.F90` | 133 | `C_H2O = 55556.0_real64` |
| `tutorials/14. cam_cloud_chemistry.ipynb` | cell ~98 | Same |

### What to change in each file

1. Replace `C_H2O = 55556.0` (the IC value) with a realistic cloud value:
   ```python
   LWC = 0.3e-3        # kg/m³ of air (typical cumulus)
   MW_H2O = 0.018015   # kg/mol
   C_H2O = LWC / MW_H2O  # ≈ 0.0167 mol/m³-air
   ```

2. Keep `C_H2O_M = 55.556` (the rate constant conversion factor) unchanged.

3. Re-derive all equilibrium IC helper functions. The IC iteration uses `C_H2O` as `[S]` in Ka1/Ka2/Kw expressions — the formulas are correct, only the initial value changes. With realistic water, expect:
   - Most SO2 stays gaseous (~96%, only ~4% as aqueous HSO3⁻)
   - Most H2O2 dissolves (~68%, α ≈ 2.16 — H2O2 is extremely soluble)
   - O3 mostly stays gaseous (very low solubility)
   - pH ≈ 5.4 (not 4.0 as with the old value)
   - R1 is extremely fast (τ ~ 0.01s for total sulfur)

4. Update test assertions to match the new equilibrium. The `test_standalone_cloud_config.py` file already has the correct values as a reference.

## Verification

`python/test/integration/test_standalone_cloud_config.py` tests the standalone `configs/v1/cam_cloud_chemistry/config.json` with realistic H2O and confirms:
- Solver converges at dt=0.01s and dt=1.0s
- SO4²⁻ is produced by R1
- Charge balance holds
- f_v = 55556 gives unphysical volume fraction > 1
