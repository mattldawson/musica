# Phase 7: MIAM in Mechanism Configuration — Detailed Plan

## Executive Summary

Phase 7 enables MPAS-A (and any Fortran host) to use MIAM aerosol chemistry
entirely through configuration files. The key architectural insight is that the
existing Fortran API — `micm_t(config_path, solver_type, error)` — should
**transparently** create a MIAM-enabled DAE solver when the mechanism
configuration file contains aerosol model definitions. No new Fortran types for
MIAM representations, processes, or constraints are needed. The config file drives
everything.

### Current State (Phases 0–3 of MIAM Integration Plan)

| Layer | MIAM Status |
|-------|-------------|
| MIAM C++ library | ✅ Header-only, fetched via CMake |
| MUSICA C++ wrapper | ✅ `CreateMicmWithMiam()` in `miam_builder.cpp` |
| MUSICA Python bindings | ✅ All types bound via pybind11 |
| Mechanism Configuration | ❌ No MIAM parsing |
| MUSICA C interface | ❌ No MIAM entry point |
| MUSICA Fortran API | ❌ Missing DAE solver enums; no MIAM path |

### Target State (after Phase 7)

| Layer | MIAM Status |
|-------|-------------|
| Mechanism Configuration | ✅ Parses aerosol models from v1 JSON/YAML |
| MUSICA C interface | ✅ `CreateMicm()` auto-detects MIAM and uses DAE solver |
| MUSICA Fortran API | ✅ DAE enums added; `micm_t(path, solver)` works with MIAM configs |
| End-to-end test | ✅ Fortran CAM Cloud Chemistry matches Python reference |

---

## Design Principle: Config-Driven, Not API-Driven

The MIAM Python API exposes typed classes (`UniformSection`,
`DissolvedReaction`, `HenryLawEquilibriumConstraint`, etc.) for programmatic
construction. Fortran does **not** need this. Fortran hosts read a config file
path and get back a solver. The config file is the API.

This means:
- **No** Fortran `miam_types` module
- **No** Fortran builder pattern for representations/processes/constraints
- **No** new C structs for MIAM types passed across the FFI boundary
- The `ReadConfiguration()` C++ function parses both gas-phase reactions
  **and** aerosol model definitions from a single config file
- `CreateMicm()` detects the presence of aerosol config and internally calls
  `CreateMicmWithMiam()` with the appropriate DAE solver

From the Fortran user's perspective:

```fortran
! Before Phase 7: gas-phase only
type(micm_t) :: solver
solver = micm_t("config.json", RosenbrockStandardOrder, error)

! After Phase 7: same call, but config.json now includes aerosol model
! → solver automatically uses DAE and includes MIAM
solver = micm_t("config.json", RosenbrockDAE4StandardOrder, error)
```

MIAM contributions appear as additional state variables and parameters in the
solver state, just as they do through the Python API. The host model sets/gets
them by name (e.g., `CLOUD.AQUEOUS.SO2_aq`, `CLOUD.MODE.AQUEOUS.H2O`).

---

## Proposed Mechanism Configuration v1 Schema for MIAM

The v1 format currently has four top-level keys: `version`, `name`, `species`,
`phases`, `reactions`. We add three new top-level keys: `condensed phases`,
`aerosol representations`, and `aerosol processes`. Constraints are inferred from
the aerosol process definitions (Henry's Law equilibrium constraints from
Henry's Law phase transfers; dissolved equilibrium constraints from reversible
reactions with `"equilibrium": true`; linear constraints from explicit
conservation/balance entries).

### Annotated Example: CAM Cloud Chemistry

```yaml
version: "1.0.0"
name: "TS1 + CAM Cloud Chemistry"

# ── Gas-Phase Species (unchanged from current v1) ──────────────────────
species:
  - name: SO2
    molecular weight [kg mol-1]: 0.06407
  - name: H2O2
    molecular weight [kg mol-1]: 0.03401
  - name: O3
    molecular weight [kg mol-1]: 0.04800
  # ... other TS1 gas species ...

# ── Gas Phase (unchanged) ──────────────────────────────────────────────
phases:
  - name: gas
    species:
      - name: SO2
        diffusion coefficient [m2 s-1]: 1.28e-5
      - name: H2O2
        diffusion coefficient [m2 s-1]: 1.46e-5
      - name: O3
        diffusion coefficient [m2 s-1]: 1.48e-5
      # ... other species ...

# ── Gas-Phase Reactions (unchanged) ────────────────────────────────────
reactions:
  - type: ARRHENIUS
    # ... existing TS1 gas-phase reactions ...

# ══════════════════════════════════════════════════════════════════════
# NEW: Aerosol Model Configuration
# ══════════════════════════════════════════════════════════════════════

# ── Condensed Phases ───────────────────────────────────────────────────
# Define condensed-phase species and their properties.
# These are separate from gas-phase species — they exist in aerosol/cloud
# droplets and are tracked by MIAM representations.
condensed phases:
  - name: AQUEOUS
    species:
      - name: SO2_aq
        molecular weight [kg mol-1]: 0.06407
      - name: H2O2_aq
        molecular weight [kg mol-1]: 0.03401
      - name: O3_aq
        molecular weight [kg mol-1]: 0.04800
      - name: Hp
        molecular weight [kg mol-1]: 0.00101
      - name: OHm
        molecular weight [kg mol-1]: 0.01701
      - name: HSO3m
        molecular weight [kg mol-1]: 0.08107
      - name: SO3mm
        molecular weight [kg mol-1]: 0.08007
      - name: SO4mm
        molecular weight [kg mol-1]: 0.09607
      - name: H2O
        molecular weight [kg mol-1]: 0.01801
        density [kg m-3]: 997.0

# ── Aerosol Representations ───────────────────────────────────────────
# How aerosol particles are represented (size distribution).
# Each representation contains one or more condensed phases.
aerosol representations:
  - type: UNIFORM_SECTION
    name: CLOUD
    phases:
      - AQUEOUS
    minimum radius [m]: 1.0e-6
    maximum radius [m]: 1.0e-5

  # Other representation types (not used in CAM Cloud, shown for schema):
  # - type: SINGLE_MOMENT_MODE
  #   name: ACCUMULATION
  #   phases:
  #     - SULFATE
  #   geometric mean radius [m]: 5.0e-8
  #   geometric standard deviation: 1.6
  #
  # - type: TWO_MOMENT_MODE
  #   name: AITKEN
  #   phases:
  #     - SULFATE
  #   geometric standard deviation: 1.8

# ── Aerosol Processes ─────────────────────────────────────────────────
# Processes acting on condensed-phase species. Each process type maps
# to a MIAM process and optionally a MIAM constraint.
aerosol processes:
  # --- Henry's Law Equilibrium (Constraint) ---
  # Gas ↔ condensed steady-state equilibrium treated as an algebraic
  # constraint in the DAE system. The condensed species concentration
  # is determined by: [aq] = HLC · R · T · liquid_water_fraction · [gas].
  # Use this when gas-liquid equilibration is fast relative to the
  # chemistry timestep (typical for cloud droplets).
  - type: HENRY_LAW_EQUILIBRIUM
    gas phase: gas
    gas-phase species: SO2
    condensed phase: AQUEOUS
    condensed-phase species: SO2_aq
    solvent: H2O
    Henry's law constant:
      HLC_ref [mol m-3 Pa-1]: 1.23
      C [K]: 3120.0

  - type: HENRY_LAW_EQUILIBRIUM
    gas phase: gas
    gas-phase species: H2O2
    condensed phase: AQUEOUS
    condensed-phase species: H2O2_aq
    solvent: H2O
    Henry's law constant:
      HLC_ref [mol m-3 Pa-1]: 7.4e4
      C [K]: 6621.0

  - type: HENRY_LAW_EQUILIBRIUM
    gas phase: gas
    gas-phase species: O3
    condensed phase: AQUEOUS
    condensed-phase species: O3_aq
    solvent: H2O
    Henry's law constant:
      HLC_ref [mol m-3 Pa-1]: 1.15e-2
      C [K]: 2560.0

  # --- Henry's Law Phase Transfer (Kinetic Process) ---
  # Gas ↔ condensed mass transfer treated as an ODE with explicit
  # diffusion and accommodation kinetics. Use this when the transfer
  # timescale is comparable to the chemistry timestep (e.g., coarse
  # aerosol or slowly-accommodating species).
  # NOT used in the CAM Cloud Chemistry example (shown for schema):
  #
  # - type: HENRY_LAW_PHASE_TRANSFER
  #   gas phase: gas
  #   gas-phase species: HNO3
  #   condensed phase: AQUEOUS
  #   condensed-phase species: HNO3_aq
  #   solvent: H2O
  #   Henry's law constant:
  #     HLC_ref [mol m-3 Pa-1]: 2.1e5
  #     C [K]: 8700.0
  #   diffusion coefficient [m2 s-1]: 1.32e-5
  #   accommodation coefficient: 0.05

  # --- Dissolved Reactions ---
  # Irreversible kinetic reaction in the condensed phase.
  - type: DISSOLVED_REACTION
    condensed phase: AQUEOUS
    solvent: H2O
    reactants:
      - species name: HSO3m
        coefficient: 1
      - species name: H2O2_aq
        coefficient: 1
    products:
      - species name: SO4mm
        coefficient: 1
      - species name: H2O
        coefficient: 1
      - species name: Hp
        coefficient: 1
    rate constant:
      type: ARRHENIUS
      A: 4.13e10       # C_H2O_M × 7.45e7 = 55556 × 7.45e7
      C: 4430.0         # [K]

  # --- Dissolved Reversible Reactions (Kinetic) ---
  # Reversible reaction in the condensed phase treated as a pair of
  # forward/reverse ODEs (NOT an algebraic constraint). Exactly 2 of
  # the 3 rate parameters must be provided:
  #   - forward rate constant
  #   - reverse rate constant
  #   - equilibrium constant
  # When the equilibrium constant is one of the two, the missing
  # rate constant is diagnosed: k_fwd = K_eq · k_rev (or vice versa).
  # Use this when the equilibration timescale is comparable to the
  # chemistry timestep; use DISSOLVED_EQUILIBRIUM when it is fast.
  # NOT used in the CAM Cloud Chemistry example (shown for schema):
  #
  # - type: DISSOLVED_REVERSIBLE_REACTION
  #   condensed phase: AQUEOUS
  #   solvent: H2O
  #   reactants:
  #     - species name: CO2_aq
  #       coefficient: 1
  #     - species name: H2O
  #       coefficient: 1
  #   products:
  #     - species name: Hp
  #       coefficient: 1
  #     - species name: HCO3m
  #       coefficient: 1
  #   forward rate constant:          # provide exactly 2 of 3
  #     type: ARRHENIUS
  #     A: 4.3e-2
  #     C: 913.0
  #   equilibrium constant:           # missing reverse is diagnosed
  #     A: 7.74e-9
  #     C: 1000.0

  # --- Dissolved Equilibria ---
  # Fast reversible reactions treated as algebraic constraints (not ODEs).
  # Each generates a DissolvedEquilibriumConstraint in the DAE system.
  - type: DISSOLVED_EQUILIBRIUM
    condensed phase: AQUEOUS
    solvent: H2O
    reactants:
      - species name: H2O
        coefficient: 2
    products:
      - species name: Hp
        coefficient: 1
      - species name: OHm
        coefficient: 1
    algebraic species: OHm
    equilibrium constant:
      A: 3.24e-18       # Kw / C_H2O_M^2 = 1e-14 / (55556^2)
      C: 0.0

  - type: DISSOLVED_EQUILIBRIUM
    condensed phase: AQUEOUS
    solvent: H2O
    reactants:
      - species name: SO2_aq
        coefficient: 1
    products:
      - species name: Hp
        coefficient: 1
      - species name: HSO3m
        coefficient: 1
    algebraic species: HSO3m
    equilibrium constant:
      A: 3.06e-4        # Ka1 / C_H2O_M = 1.7e-2 / 55556
      C: 2090.0

  - type: DISSOLVED_EQUILIBRIUM
    condensed phase: AQUEOUS
    solvent: H2O
    reactants:
      - species name: HSO3m
        coefficient: 1
    products:
      - species name: Hp
        coefficient: 1
      - species name: SO3mm
        coefficient: 1
    algebraic species: SO3mm
    equilibrium constant:
      A: 1.08e-9        # Ka2 / C_H2O_M = 6.0e-8 / 55556
      C: 1120.0

  # --- Linear Constraints ---
  # Mass conservation and charge balance. Each defines a linear
  # relationship among species concentrations that must hold at all times.
  # The "algebraic species" is the variable solved algebraically to
  # enforce the constraint.
  - type: LINEAR_CONSTRAINT
    name: sulfur conservation
    algebraic phase: gas
    algebraic species: SO2
    terms:
      - phase: gas
        species: SO2
        coefficient: 1.0
      - phase: AQUEOUS
        species: SO2_aq
        coefficient: 1.0
      - phase: AQUEOUS
        species: HSO3m
        coefficient: 1.0
      - phase: AQUEOUS
        species: SO3mm
        coefficient: 1.0
    constant [mol m-3]: 3.01e-8

  - type: LINEAR_CONSTRAINT
    name: H2O2 conservation
    algebraic phase: gas
    algebraic species: H2O2
    terms:
      - phase: gas
        species: H2O2
        coefficient: 1.0
      - phase: AQUEOUS
        species: H2O2_aq
        coefficient: 1.0
    constant [mol m-3]: 3.01e-8

  - type: LINEAR_CONSTRAINT
    name: O3 conservation
    algebraic phase: gas
    algebraic species: O3
    terms:
      - phase: gas
        species: O3
        coefficient: 1.0
      - phase: AQUEOUS
        species: O3_aq
        coefficient: 1.0
    constant [mol m-3]: 1.5e-6

  - type: LINEAR_CONSTRAINT
    name: charge balance
    algebraic phase: AQUEOUS
    algebraic species: Hp
    terms:
      - phase: AQUEOUS
        species: Hp
        coefficient: 1.0
      - phase: AQUEOUS
        species: OHm
        coefficient: -1.0
      - phase: AQUEOUS
        species: HSO3m
        coefficient: -1.0
      - phase: AQUEOUS
        species: SO3mm
        coefficient: -2.0
      - phase: AQUEOUS
        species: SO4mm
        coefficient: -2.0
    constant [mol m-3]: 0.0
```

### Schema Summary

| Top-Level Key | Status | Content |
|---|---|---|
| `version` | Existing | `"1.0.0"` (no version bump needed — additive extension) |
| `name` | Existing | Mechanism name |
| `species` | Existing | Gas-phase species with properties |
| `phases` | Existing | Gas-phase definition with diffusion coefficients |
| `reactions` | Existing | Gas-phase reactions (Arrhenius, Troe, Photolysis, etc.) |
| `condensed phases` | **New** | Condensed-phase species and properties |
| `aerosol representations` | **New** | Particle size representations (UNIFORM_SECTION, SINGLE_MOMENT_MODE, TWO_MOMENT_MODE) |
| `aerosol processes` | **New** | Aerosol processes and constraints (HENRY_LAW_EQUILIBRIUM, HENRY_LAW_PHASE_TRANSFER, DISSOLVED_REACTION, DISSOLVED_REVERSIBLE_REACTION, DISSOLVED_EQUILIBRIUM, LINEAR_CONSTRAINT) |

### Design Decisions

1. **Constraints and processes share `aerosol processes`.** Both kinetic
   processes (ODE) and algebraic constraints live in the same section because
   they participate in the same coupled system. For Henry's Law partitioning,
   users choose per species between the equilibrium form
   (`HENRY_LAW_EQUILIBRIUM` — algebraic constraint, no kinetics) and the
   kinetic form (`HENRY_LAW_PHASE_TRANSFER` — ODE with diffusion and
   accommodation coefficients). The equilibrium form is appropriate when
   gas-liquid equilibration is fast; the kinetic form when transfer timescales
   matter.

2. **Condensed-phase species are separate from gas-phase species.** The `species`
   array remains gas-phase only. Condensed-phase species live in
   `condensed phases`. This matches the physical separation (gas vs.
   aqueous/solid) and avoids ambiguity about which phase a species belongs to.

3. **Representation ↔ phase binding is by name.** Each representation lists which
   condensed phases it contains. The MIAM builder resolves the binding.

4. **Rate constant objects are inline.** Rather than referencing named rate
   constants, each process embeds its rate constant. This is consistent with how
   gas-phase reactions already embed their parameters.

5. **No version bump.** These are additive keys. Configs without `condensed phases` /
   `aerosol representations` / `aerosol processes` parse identically to before.
   The v1 parser simply checks whether these optional sections exist.

---

## Implementation Plan

### Step 7.1: Add DAE Solver Types to Fortran (`fortran/micm/micm.F90`)

**Effort:** Small  
**Files:** `fortran/micm/micm.F90`

Add the four missing DAE solver enumerators to match the C++ `MICMSolver` enum:

```fortran
enum, bind(c)
  enumerator :: UndefinedSolver            = 0
  enumerator :: Rosenbrock                 = 1
  enumerator :: RosenbrockStandardOrder    = 2
  enumerator :: BackwardEuler              = 3
  enumerator :: BackwardEulerStandardOrder = 4
  enumerator :: CudaRosenbrock             = 5
  enumerator :: RosenbrockDAE4             = 6   ! NEW
  enumerator :: RosenbrockDAE4StandardOrder = 7  ! NEW
  enumerator :: RosenbrockDAE6             = 8   ! NEW
  enumerator :: RosenbrockDAE6StandardOrder = 9  ! NEW
end enum
```

No other Fortran changes. The existing `micm_t` constructor passes solver_type
as an integer to the C `CreateMicm()` function, which already supports these
values.

**Verification:** Extend the existing Fortran MICM unit test to create a solver
with `RosenbrockDAE4StandardOrder` using a gas-phase-only config. Confirm it
solves without error.

---

### Step 7.2: Extend Mechanism Configuration Parser for MIAM

**Effort:** Large (core of Phase 7)  
**Files:**
- `build-miam/_deps/mechanism_configuration-src/` (or wherever the library lives)
  — extend v1 parser
- Alternatively, if Mechanism Configuration is an external dependency:
  submit changes upstream or fork

**Sub-tasks:**

**7.2a — Data Model Extension**

Add types to the Mechanism Configuration library's internal data model:

```
CondensedPhase        { name, species[] }
AerosolRepresentation { type, name, phases[], min_radius, max_radius,
                        geometric_mean_radius, geometric_standard_deviation }
AerosolProcess        { type, ... (varies by type) }
```

These live alongside the existing `Species`, `Phase`, `Reaction` types. The
top-level `Mechanism` struct gains:

```cpp
struct Mechanism {
  Version version;
  std::string name;
  std::vector<Species> species;           // existing
  std::vector<Phase> phases;              // existing
  std::vector<Reaction> reactions;        // existing
  std::vector<CondensedPhase> condensed_phases;          // NEW
  std::vector<AerosolRepresentation> aerosol_representations;  // NEW
  std::vector<AerosolProcess> aerosol_processes;               // NEW
};
```

**7.2b — v1 Parser Extension**

Extend the v1 YAML/JSON parser to read the three new sections. Each aerosol
process type maps to a parser function:

| Config `type` | Parser Function | MIAM Type |
|---|---|---|
| `HENRY_LAW_EQUILIBRIUM` | `ParseHenryLawEquilibrium()` | Constraint (algebraic) |
| `HENRY_LAW_PHASE_TRANSFER` | `ParseHenryLawPhaseTransfer()` | Process (kinetic ODE) |
| `DISSOLVED_REACTION` | `ParseDissolvedReaction()` | Process (irreversible) |
| `DISSOLVED_REVERSIBLE_REACTION` | `ParseDissolvedReversibleReaction()` | Process (forward + reverse ODEs) |
| `DISSOLVED_EQUILIBRIUM` | `ParseDissolvedEquilibrium()` | Constraint (algebraic) |
| `LINEAR_CONSTRAINT` | `ParseLinearConstraint()` | Constraint |

Unknown keys starting with `__` are ignored (existing convention for comments).

**7.2c — Conversion to MIAM Config Types**

Add a conversion function (analogous to existing `ParserV1()` → `Chemistry`):

```cpp
miam_config::ModelConfig ConvertToMiamConfig(const Mechanism& mechanism);
```

This maps parsed `CondensedPhase` → `miam_config::PhaseDef`,
`AerosolRepresentation` → `miam_config::Representation`,
and `AerosolProcess` → `miam_config::Process` / `miam_config::Constraint`.

**7.2d — Unit Tests**

- Parse the CAM Cloud Chemistry YAML and verify all fields round-trip correctly.
- Parse a config with no aerosol sections and confirm it produces an empty
  `ModelConfig` (backward compatibility).
- Parse configs with each representation type (UNIFORM_SECTION,
  SINGLE_MOMENT_MODE, TWO_MOMENT_MODE).
- Error handling: missing required fields, unknown types, invalid references.

---

### Step 7.3: Update `ReadConfiguration()` and `CreateMicm()` for MIAM

**Effort:** Medium  
**Files:**
- `src/micm/parse.hpp` / `src/micm/parse.cpp`
- `src/micm/micm_c_interface.cpp`
- `include/musica/micm/micm_c_interface.hpp`

**7.3a — Extend `ReadConfiguration()` Return Type**

Currently returns `Chemistry` (gas-phase only). Change to return a richer type:

```cpp
struct MechanismConfig {
  Chemistry chemistry;                        // gas-phase (existing)
  std::optional<miam_config::ModelConfig> miam_config;  // aerosol (new)
};

MechanismConfig ReadConfiguration(const std::string& config_path);
```

When the parsed `Mechanism` has non-empty `condensed_phases` /
`aerosol_representations` / `aerosol_processes`, populate `miam_config`.
Otherwise leave it as `std::nullopt`.

**7.3b — Update `CreateMicm()` to Auto-Detect MIAM**

```cpp
MICM* CreateMicm(const char* config_path, MICMSolver solver_type, Error* error)
{
  auto config = ReadConfiguration(config_path);
  if (config.miam_config.has_value()) {
    return CreateMicmWithMiam(
        config.chemistry, solver_type, *config.miam_config, error);
  } else {
    MICM* micm = new MICM(config.chemistry, solver_type);
    NoError(error);
    return micm;
  }
}
```

This is the key change that makes MIAM transparent to Fortran. The same
`CreateMicm()` C function (already bound from Fortran) now handles both
gas-only and gas+aerosol configs.

**7.3c — Similarly Update `CreateMicmFromConfigString()`**

Same pattern for the string-based variant, if used.

**7.3d — Update `CreateMicmFromChemistryMechanism()`**

This function takes a pre-parsed `Chemistry` object. It cannot support MIAM
through the current signature (no MIAM config). Options:
- Leave as-is (gas-only path for programmatic API users).
- Add a new overload that accepts both `Chemistry` and `miam_config::ModelConfig`.

For Fortran/config-driven use, this doesn't matter — they go through
`CreateMicm(config_path)`.

---

### Step 7.4: Create CAM Cloud Chemistry Config File

**Effort:** Small  
**Files:** `configs/v1/cam_cloud_chemistry/config.json` (or `.yaml`)

Write the complete CAM Cloud Chemistry configuration in v1 format, following the
schema above. Combine:
- Gas-phase species and reactions from TS1
- Aerosol model (CLOUD representation, AQUEOUS phase, all processes/constraints)

This config is both a test fixture and a reference for users.

---

### Step 7.5: Python Integration Test from Config

**Effort:** Medium  
**Files:** `python/test/integration/test_miam_from_config.py`

The existing `test_miam_cloud_chemistry.py` builds the MIAM model
programmatically. Add a new test that:

1. Parses the CAM Cloud Chemistry config file from Step 7.4
2. Creates the solver via `MICM(config_path, solver_type)` (not programmatic API)
3. Runs the same 10-second integration
4. Compares results against the programmatic test (should be identical)

This validates the full config → parser → builder → solver path before we
involve Fortran.

---

### Step 7.6: Fortran Unit Test

**Effort:** Medium  
**Files:** `fortran/test/test_miam_cloud_chemistry.F90`

Create a Fortran test that:

1. Creates `micm_t(cam_cloud_config_path, RosenbrockDAE4StandardOrder, error)`
2. Gets a state: `state = solver%get_state(1, error)`
3. Sets initial conditions (T, P, gas concentrations, liquid water)
4. Calls `solver%solve(time_step, state, ...)`
5. Verifies SO4²⁻ production and mass conservation
6. Compares key values against the Python reference

The Fortran test should be added to CMake's test suite.

---

### Step 7.7: MPAS Integration (Phase 8 Preview)

**Not in scope for Phase 7**, but the design enables it:

The MPAS chemistry driver already calls `micm_t(config_path, solver_type, error)`.
When Phase 8 uses a TS1+CAM Cloud config:

1. The solver automatically includes MIAM (detected from config)
2. The driver queries state variable names — MIAM variables appear alongside
   gas-phase variables (e.g., `CLOUD.AQUEOUS.SO2_aq`, `CLOUD.AQUEOUS.Hp`)
3. The driver maps MPAS cloud liquid water content → `CLOUD.MODE.AQUEOUS.H2O`
4. MIAM aerosol species are advected as tracers (same as gas-phase species)

No special MIAM code in the Fortran driver — it's all config-driven.

---

## Dependency Graph

```
Step 7.1 ──────────────────────────────────────────┐
  (DAE enums in Fortran)                           │
                                                   │
Step 7.2 ─────────────────┐                        │
  (MechConfig parser)     │                        │
                          ▼                        │
                   Step 7.3                        │
                   (CreateMicm auto-detect)        │
                          │                        │
                          ▼                        │
                   Step 7.4                        │
                   (CAM Cloud config file)         │
                          │                        │
                     ┌────┴────┐                   │
                     ▼         ▼                   │
              Step 7.5      Step 7.6 ◄─────────────┘
              (Python test) (Fortran test)
```

Steps 7.1 and 7.2 can proceed in parallel. Step 7.3 depends on 7.2. Steps 7.5
and 7.6 depend on 7.3 + 7.4, and 7.6 additionally depends on 7.1.

---

## Risk Areas and Mitigations

### 1. Mechanism Configuration Is an External Dependency

The Mechanism Configuration library is fetched via CMake (`FetchContent`). If
it's maintained in a separate repo, the schema extension (Step 7.2) requires
either:
- A PR to that upstream repo (preferred — keeps the library canonical)
- A fork or local patch

**Mitigation:** Coordinate with the Mechanism Configuration maintainers early.
The extension is additive (no breaking changes to existing parsing).

### 2. Rate Constant Functions vs. Named Forms

MIAM processes use `std::function<double(const Conditions&)>` for rate
constants. The config schema uses named forms (`ARRHENIUS` type with A/C
parameters). This covers the CAM Cloud Chemistry case and most atmospheric
chemistry use cases. User-defined rate constant functions cannot be specified
in config and would require programmatic construction.

**Mitigation:** Document the supported rate constant forms. For CAM Cloud
Chemistry (the Phase 8 target), Arrhenius + temperature-dependent forms are
sufficient.

### 3. Constraint Initialization Parameters

DAE solvers need `constraint_init_max_iterations` and
`constraint_init_tolerance` for Newton-based constraint initialization. These
are solver parameters, not mechanism parameters.

**Mitigation:** Expose these through the existing
`RosenbrockSolverParameters` C/Fortran interface (the C struct already has
these fields: `constraint_init_max_iterations` and
`constraint_init_tolerance`). The Fortran API's
`set_rosenbrock_solver_parameters()` needs to pass these through. Default
values (100 iterations, 1e-8 tolerance) work for CAM Cloud Chemistry.

### 4. `CreateMicm()` Solver Type Validation

If the config contains aerosol processes but the caller passes a non-DAE
solver type (e.g., `RosenbrockStandardOrder`), the builder will fail because
MIAM constraints require a DAE solver.

**Mitigation:** `CreateMicm()` should either:
- (a) Auto-upgrade to DAE4 if MIAM config is present and a non-DAE solver was
  requested (with a log warning), or
- (b) Return an error with a clear message: "Config contains aerosol model;
  use a DAE solver type (RosenbrockDAE4StandardOrder recommended)."

Option (b) is safer and more explicit.

---

## Verification Strategy

### Fortran Unit Test (Step 7.6)

The CAM Cloud Chemistry test case is well-defined:
- 3 gas species, 9 aqueous species, 1 kinetic reaction, 9 constraints
- 10-second integration at T=280K, P=70000Pa
- Expected: SO₂ oxidation to SO₄²⁻, mass conservation of S/H2O2/O3

Compare against the Python reference in `test_miam_cloud_chemistry.py`.

### MPAS Notebook (Phase 8, future)

`verification/phase08_cloud_chemistry.ipynb` will prescribe a cloud layer in the
JW test and verify sulfate production in cloudy cells.

### Regression

All existing tests (gas-phase configs, Chapman, TS1) must continue to pass.
`CreateMicm()` with configs that have no `condensed phases` section must behave
identically to before.

---

## Open Questions for Review

1. **Should `HENRY_LAW_EQUILIBRIUM` entries require the solvent's molecular
   weight and density, or should those be inferred from the condensed phase
   species definition?** The MIAM `HenryLawEquilibriumConstraint` builder
   needs `Mw_solvent` and `rho_solvent`. The proposed schema puts these on
   the solvent species in `condensed phases` (via `molecular weight` and
   `density` properties). The parser would look them up by solvent name.

2. **Should `CreateMicm()` auto-upgrade non-DAE solver types when MIAM is
   present, or return an error?** Proposed: return an error (explicit is
   better than implicit).

3. **Should `LINEAR_CONSTRAINT` `constant` values be fixed in config, or
   should the host model be able to set them at runtime?** For CAM Cloud
   Chemistry, the constants are the total mass of S, H2O2, O3 in each cell —
   these may need to be computed from initial conditions rather than
   hardcoded. This could be handled by having the host set them as solver
   parameters, similar to how emissions rates are set.

4. **Should the Mechanism Configuration extension live in the mechanism_configuration
   repo (upstream) or in MUSICA?** The conversion layer (Step 7.2c,
   `ConvertToMiamConfig`) can live in MUSICA regardless, but the parser types
   (Step 7.2a–b) should ideally be upstream.
