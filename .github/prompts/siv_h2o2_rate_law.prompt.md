# S(IV) + H₂O₂ Oxidation Rate Law — Correct Implementation

## Summary

The aqueous-phase oxidation of S(IV) by H₂O₂ has a **non-standard rate law**
with an acid-catalyzed mechanism. The original MIAM R1 used a simple bimolecular
rate (HSO₃⁻ + H₂O₂), which was ~10⁴× too fast at cloud pH because it omitted the
[H⁺]-dependence. This document describes the correct rate law, the current fix,
and the long-term plan.

## Literature Rate Law

From Hoffmann & Calvert (1985), Martin & Damschen (1981), Seinfeld & Pandis Ch. 7:

```
-d[S(IV)]/dt = k · [H⁺] · [HSO₃⁻] · [H₂O₂(aq)] / (1 + K · [H⁺])
```

Parameters at 298 K:
- **k** = 7.45 × 10⁷ M⁻² s⁻¹  (third-order in [H⁺][HSO₃⁻][H₂O₂])
- **K** = 13 M⁻¹  (saturation constant for the intermediate)
- **Eₐ/R** = 4430 K

### CAM Implementation (mo_setsox / cloud_aqueous_chemistry.F90)

CAM uses an equivalent form with [SO₂(aq)] instead of [HSO₃⁻]:

```fortran
k_siv_h2o2 = 8.e4 * EXP(-3650.*(1./T - 1./298.)) / (.1 + h_plus_conc)
dso4_dt = k_siv_h2o2 * h2o2_aq * so2_aq * molar_to_mixing_ratio
```

Note: CAM computes pH first (bisection on charge balance), then uses it as a
constant in the oxidation step. The `/ (.1 + [H⁺])` encodes the full mechanism.

### Mechanism (Hoffmann 1986)

The rate law arises from a two-step mechanism:

```
Step 1:  HSO₃⁻ + H₂O₂(aq) ⇌ SO₂OOH⁻ + H₂O     (fast pre-equilibrium, Keq)
Step 2:  SO₂OOH⁻ + H⁺ → H₂SO₄(aq) → 2H⁺ + SO₄²⁻  (rate-determining)
```

Applying steady-state to the intermediate SO₂OOH⁻ gives:

```
rate = k₂ · Keq · [HSO₃⁻] · [H₂O₂] · [H⁺] / (1 + Keq · [H₂O₂])
```

Under dilute cloud conditions ([H₂O₂] ≪ 1/Keq), this simplifies to the
standard form above (with k = k₂·Keq and K = 13 from the [H⁺] saturation).

## The Bug

The original MIAM R1 was defined as a bimolecular dissolved reaction:

```
Reactants: HSO₃⁻, H₂O₂     →  rate = A · exp(-C/T') / [S] · [HSO₃⁻] · [H₂O₂]
```

with A = C_H2O_M × 7.45e7 (≈ 4.14e9), C = 4430.

This is missing the [H⁺] factor in the numerator. At cloud pH ≈ 4:
- [H⁺] ≈ 10⁻⁴ M
- Missing factor ≈ [H⁺] / (1 + K·[H⁺]) ≈ 10⁻⁴
- **Rate was ~10,000× too fast**, causing extreme stiffness (τ ≈ 0.2 ms)

## Current Fix (April 2026)

Include H⁺ as a **third reactant**. Since MIAM's DissolvedReaction computes:

```
rate = k(T) / [S]^(n_r - 1) · ∏[reactants]
```

With 3 reactants:
```
rate = k(T) / [S]² · [HSO₃⁻] · [H₂O₂] · [H⁺]
```

This correctly captures the [H⁺] numerator. The (1 + K·[H⁺]) denominator is
dropped — at cloud pH > 3, this term is < 1.3% (pH 3: 1.3%, pH 4: 0.13%).

### Updated Stoichiometry

```
Reactants:  HSO₃⁻ + H₂O₂(aq) + H⁺
Products:   SO₄²⁻ + H₂O + H⁺ + H⁺    (repeat for coefficient 2)
Net:        HSO₃⁻ + H₂O₂ → SO₄²⁻ + H₂O + H⁺  (unchanged)
```

Note: MIAM's builder encodes stoichiometric coefficients by **repeating species
names** in the product list. The JSON `"coefficient"` field is ignored by
`miam_config_convert.cpp`.

### Updated Rate Constant

For 3 reactants in MIAM units (mol/m³-air):

```
A_miam = C_H2O_M² × k_lit = 55.34² × 7.45e7 ≈ 2.282e11
C_miam = 4430.0 K   (unchanged)
```

### Physical Effect

At pH 4.27, T = 287.45 K, the corrected rate is:
- Old timescale:  τ ≈ 0.175 ms  (caused NaN/stiffness failure)
- New timescale:  τ ≈ minutes    (compatible with 900 s timestep)

## Long-Term Plan: Unfold Into Intermediate Reactions

To capture the full rate law including the (1 + K·[H⁺]) denominator at low pH,
unfold R1 into two elementary reactions with the SO₂OOH⁻ intermediate:

```json
{
  "__comment": "R1a: HSO3- + H2O2 <-> SO2OOH- + H2O (fast pre-equilibrium)",
  "type": "DISSOLVED_REVERSIBLE_REACTION",
  "reactants": ["HSO3m", "H2O2"],
  "products": ["SO2OOHm", "H2O"],
  "equilibrium_constant": { "A": Keq, "C": ... }
}

{
  "__comment": "R1b: SO2OOH- + H+ -> SO4-- + 2H+ (rate-determining)",
  "type": "DISSOLVED_REACTION",
  "reactants": ["SO2OOHm", "Hp"],
  "products": ["SO4mm", "Hp", "Hp"],
  "rate_constant": { "A": k2_miam, "C": 4430.0 }
}
```

### Required Changes for Intermediate Approach

1. **Add SO₂OOH⁻ species** to the mechanism (peroxysulfurous acid anion)
2. **Determine Keq and k₂ separately** — literature gives their product (k = k₂·Keq)
   and K = 13 M⁻¹ for the saturation. Need to find individual values from
   Hoffmann (1986) or McArdle & Hoffmann (1983).
3. The fast pre-equilibrium (R1a) could use `DissolvedReversibleReaction` with
   a large forward rate or be handled as a `DissolvedEquilibriumConstraint` if
   the solver can handle the stiffness.

### When to Do This

The intermediate approach is only needed if:
- Cloud pH drops below ~3 (volcanic/industrial plumes)
- SO₂ concentrations are very high (> 10 ppb)
- Sub-second accuracy is needed in the oxidation kinetics

For typical tropospheric cloud chemistry (pH 3–6, SO₂ < 10 ppb), the current
3-reactant approach is accurate to better than 1.3%.

## References

- Hoffmann, M.R. and Calvert, J.G. (1985). Chemical transformation modules
  for Eulerian acid deposition models, Vol. 2, EPA/600/3-85/017.
- Martin, L.R. and Damschen, D.E. (1981). Aqueous oxidation of sulfur dioxide
  by hydrogen peroxide at low pH. Atmos. Environ., 15, 1615–1621.
- Seinfeld, J.H. and Pandis, S.N. (2016). Atmospheric Chemistry and Physics,
  3rd ed., Chapter 7.
- McArdle, J.V. and Hoffmann, M.R. (1983). Kinetics and mechanism of the
  oxidation of aquated sulfur dioxide by hydrogen peroxide at low pH.
  J. Phys. Chem., 87, 5425–5429.
- CAM source: `cloud_aqueous_chemistry.F90` in NCAR/CAM-ACOM-dev
  (commit bb9cb1f, Matt Dawson's refactoring of mo_setsox.F90)
- MIAM design notes: `.github/design/notes_on_cam_cloud_chem.md` in NCAR/miam
