---
description: "CheMPAS-A development plan: integrating MUSICA chemistry (MICM, TUV-x, MIAM) into MPAS-A with runtime-configurable species and aerosol representations. Use when working on MPAS chemistry integration, CheMPAS-A, or the MIAM Fortran bindings."
---

# CheMPAS-A: MUSICA Chemistry in MPAS-A

## Repositories

- **MUSICA**: development branch (this repo)
- **MPAS-A**: https://github.com/mattldawson/MPAS-Model

## Goal

Integrate TS1 gas-phase chemistry (MICM), CAM Cloud aqueous chemistry (MIAM), and
TUV-x photolysis into MPAS-A. Docker-based local development on the 480-km
quasi-uniform mesh (2,562 cells). Start with Chapman chemistry (hardcoded species),
then add **runtime species configuration** so MPAS tracers are dynamically allocated
from the MICM mechanism at init-time (following the proven CAM-MPAS
`atm_allocate_scalars` pattern). After MIAM Fortran bindings and Mechanism
Configuration updates, add **runtime aerosol representation configuration** so
aerosol fields are driven by the MIAM config. One build runs any mechanism.

---

## Phase 0: Local Build & Test Infrastructure

1. **Create Dockerfile** — Ubuntu with gfortran, OpenMPI, NetCDF, PNetCDF, PIO2, CMake,
   and MUSICA-Fortran pre-installed via `pkg-config`. Multi-stage build: base image with
   deps, dev image with source mounted.
2. **Create `devcontainer.json`** — VS Code dev container with `docker-compose.yml`.
3. **Mesh download script** — Fetch the 480-km mesh (`x1.2562.tar.gz`, 1.5 MB) and static
   file (`x1.2562_static.tar.gz`, 1.0 MB) from `www2.mmm.ucar.edu` into a gitignored
   `data/` directory. *(parallel with 1)*
4. **Build vanilla MPAS-A** — Verify `make CORE=atmosphere gfortran` and
   `make CORE=init_atmosphere gfortran`. Run the Jablonowski-Williamson baroclinic wave
   idealized test case on the 480-km mesh with 1–4 MPI ranks. *(depends on 1, 3)*
5. **GFS-initialized real-data run** — Download sample GFS analysis data, run
   `init_atmosphere_model` → 6-hour forecast on 480-km mesh. *(parallel with 4)*
6. **Test harness** — Shell scripts + CTest: build verification, idealized/real-data
   execution, `ncdump`-based output validation. *(depends on 4, 5)*
7. **CI pipeline** — GitHub Actions: build Docker → compile → test on every PR. 480-km
   mesh, 1–2 MPI ranks. *(depends on 6)*

**Verification:**
- JW baroclinic wave 24 simulated hours in < 5 min wall clock (2 MPI ranks)
- GFS 6-hour forecast completes on 480-km mesh
- CI pipeline passes green on a test PR

---

## Phase 1: Passive Tracer Infrastructure

8. **Add passive tracers to `Registry.xml`** — 2–3 test tracers advected by the dynamical
   core. Validates the tracer advection pathway before chemistry.
9. **Initialize tracers** — Gaussian blob + uniform background.
10. **Validate tracer advection** — Mass conservation, no negatives with monotonic limiter.
    Add output stream for tracer diagnostics.
11. **Tracer indexing unit test** — Standalone Fortran test for name ↔ index mapping
    (reused in Phase 3 species mapping).

**Verification:**
- Passive tracers conserve mass to machine precision over 48-hour JW test
- `ncdump -h` on output shows expected tracer variables
- Tracer indexing unit test passes

---

## Phase 2: Chapman Chemistry (Simplest MICM + TUV-x)

12. **Build with MUSICA** — `make CORE=atmosphere MUSICA=true gfortran`. Docker image
    already has MUSICA installed.
13. **Create `mpas_chemistry_driver.F90`** — New module in
    `src/core_atmosphere/chemistry/`:
    - `chemistry_init()` → creates `micm_t` (Chapman config) + `tuvx_t` (Chapman TUV-x)
    - `chemistry_timestep()` → per-column: TUV-x photolysis → map state →
      `micm%solve()` → update tracers
    - `chemistry_finalize()` → cleanup
14. **Register Chapman species** — O, O1D, O3 as advected tracers in `Registry.xml`.
    M and O2 diagnosed from air density.
15. **Wire into physics driver** — Operator-split after radiation. Namelist options:
    `config_chemistry_enabled`, `config_chemistry_config_path`.
16. **Package Chapman configs** — Include mechanism JSON + TUV-x config.
17. **Integration test** — JW + Chapman on 480-km mesh; verify O3 diurnal cycle; compare
    single-column against standalone Chapman tutorial (tutorial 10).

**Verification:**
- MPAS + Chapman compiles and runs without crashes on 480-km mesh
- O3 shows diurnal cycle (production dayside, none nightside)
- Single-column values match tutorial 10 within solver tolerance
- Mass conservation for O + O1D + O3
- Real-data + Chapman runs for 6 hours

---

## Phase 3: Runtime Species Configuration from MICM

*Key architectural phase. Instead of editing `Registry.xml` per mechanism, the
chemistry driver reads the MICM mechanism JSON at init-time and dynamically
allocates tracers. MPAS supports this via `atm_allocate_scalars` (used when MPAS
is a CAM dycore).*

18. **Implement `atm_allocate_scalars_musica()`** — At init-time:
    - Read MICM mechanism config (JSON)
    - Query MICM for species names/properties (`micm%get_species_ordering()`)
    - Allocate `scalars` var_array with correct tracer count
    - Build name → index mapping arrays
    - Classify species as "advected" (long-lived) vs. "diagnostic" (short-lived
      radicals) based on species properties or a config file
19. **Refactor chemistry driver** — Remove hardcoded Chapman species references. Replace
    with generic loops over MICM species indices. Driver becomes **mechanism-agnostic**.
20. **Mechanism-switching test** — Swap `config_chemistry_config_path` from Chapman →
    analytical test mechanism (`configs/v0/analytical/`) without recompilation. Both run.

**Verification:**
- Chapman produces identical results to Phase 2
- Chapman → analytical switch works without recompiling
- Species names in output NetCDF match the mechanism definition
- No regressions in existing tests

---

## Phase 4: Emissions & Deposition Stubs

*Provides minimum-viable emissions, dry deposition, and wet deposition so that TS1
(Phase 6) produces reasonable-ish results. These are simple parameterizations that
feed MICM's existing `EMISSION` and `FIRST_ORDER_LOSS` reaction types — not
detailed inventories. Detailed emissions from files and refined deposition handling
are deferred to future work.*

21. **Stub emissions module** — `mpas_chemistry_emissions.F90`. For each MICM `EMISSION`
    reaction, compute a surface emission rate and pass it to MICM as a user-defined
    parameter before each `solve()` call. Rates come from a simple config file
    (species name → constant flux in molecules cm⁻² s⁻¹). Inject into lowest model
    level only, scaled by grid-cell area and layer thickness. Enough for NOx, CO, VOCs
    at order-of-magnitude levels.
22. **Stub dry deposition module** — `mpas_chemistry_dry_deposition.F90`. For each MICM
    `FIRST_ORDER_LOSS` reaction flagged as dry deposition, compute a first-order loss
    rate = v_d / Δz at the lowest level, where v_d (deposition velocity, cm s⁻¹) is
    read per-species from a config file. Simple land/ocean distinction using MPAS
    land-use fraction if available, otherwise a single global value per species.
23. **Stub wet deposition module** — `mpas_chemistry_wet_deposition.F90`. For each MICM
    `FIRST_ORDER_LOSS` reaction flagged as wet deposition, apply a scavenging loss rate
    = α × P, where α (scavenging coefficient, per mm hr⁻¹) is per-species from config
    and P is the MPAS precipitation rate. Applied in levels where MPAS cloud water > 0.
    Soluble species only (flagged in config).
24. **Wire into chemistry driver & test** — Call emissions before `micm%solve()`;
    apply deposition loss rates as MICM parameters. Test with Chapman (verify it still
    works — Chapman has no emission/loss reactions so stubs are no-ops) and with a
    trivial test mechanism that has one emitted and one deposited species. Verify mass
    budget: total source = emissions − dry dep − wet dep ± chemistry.

**Verification:**
- Chapman results unchanged (stubs are no-ops for mechanisms without EMISSION/FIRST_ORDER_LOSS)
- Test mechanism with one emitted + one deposited species reaches steady state
- Emission and deposition rates in diagnostics output match config values
- Mass budget closes to solver tolerance
- No hardcoded species names — stubs are mechanism-agnostic via MICM reaction metadata

---

## Phase 5: TUV-x TS1/TSMLT Photolysis

25. **Upgrade TUV-x config** — Switch to `configs/tuvx/ts1_tsmlt.json` for the full
    photolysis rate set needed by TS1.
26. **Couple MPAS vertical coordinate to TUV-x** — Map terrain-following height →
    TUV-x altitude grid. Update O2/O3/air profiles from MPAS state each chemistry step.
27. **Implement alias mapping** — Route TUV-x reaction labels (`jo2_a`, `jo3->O1D`, …)
    to MICM photolysis parameters using the aliasing config.
28. **Regression test** — Verify Chapman photolysis rates unchanged after TUV-x upgrade.

**Verification:**
- Photolysis rates for Chapman reactions match previous phase
- Alias mapping routes TUV-x output correctly
- No numerical issues from MPAS → TUV-x profile updates

---

## Phase 6: TS1 Mechanism (~70 species, ~200 reactions)

29. **Switch config to TS1** — Point `config_chemistry_config_path` at
    `configs/v1/ts1/ts1.json`. Thanks to Phase 3, this should "just work" — runtime
    allocation handles the new species set automatically.
30. **Initial conditions** — Default from `configs/v1/ts1/initial_conditions.csv`. For
    real-data runs, read from WACCM/CAM-chem output.
31. **Populate emissions/deposition configs for TS1** — Fill in the Phase 4 stub config
    files with order-of-magnitude rates for TS1 species (NOx surface emissions, O3/NO2
    dry deposition velocities, soluble species wet scavenging coefficients). Values from
    literature or CAM-chem defaults — good enough for plausible concentrations, not
    publication-quality.
32. **Performance profiling** — 480-km mesh (2,562 × 55 levels ≈ 141K points). Target:
    < 2 sec/chemistry step. GPTL timers.
33. **Integration test** — 24-hour forecast; O3, NO, NO2, CO diurnal patterns; mass budgets.

**Verification:**
- TS1 runs with NO code changes to the chemistry driver (config-only switch)
- Stable 48-hour run; no negative species
- O3 and NO2 show expected diurnal patterns
- < 2 sec/step on workstation
- NOy/Ox mass budgets close

---

## Phase 7: MIAM Fortran Bindings & Mechanism Configuration in MUSICA

*Blocker for Phases 8–9. Work happens in this MUSICA repo.*

34. **Design MIAM Fortran API** — Mirror C++/Python pattern. Key types: `miam_model_t`,
    representations, processes, constraints. Design review before implementation.
35. **Implement MIAM C interface** — `src/miam/miam_c_interface.cpp` following
    `src/micm/micm.cpp` pattern. Functions: `CreateMiamModel`, `DeleteMiamModel`,
    `AddMiamToMicm`.
36. **Implement MIAM Fortran bindings** — `fortran/miam/miam.F90` using `bind(C)`,
    following `fortran/micm/micm.F90` as template.
37. **Add DAE solver types to Fortran API** — `RosenbrockDAE4`, `RosenbrockDAE6` in
    `SolverType` enum. `micm_t` constructor accepts external models.
38. **Extend Mechanism Configuration for MIAM** — If not already done by this point
    (see `docs/MIAM_INTEGRATION_PLAN.md` Phase 4):
    - Extend MechanismConfiguration C++ library to parse MIAM aerosol model definitions
      (representations, processes, constraints) from JSON/YAML
    - Add Fortran bindings for the extended parser so that config-driven MIAM setup
      works from Fortran (not just programmatic API)
    - Wire parsed MIAM configs into solver construction through the C/Fortran interface
    - This enables MPAS to set up MIAM from a config file path, just like it does for
      MICM — no programmatic Fortran construction of MIAM objects needed
39. **Fortran unit test** — CAM Cloud Chemistry config: create MIAM model from config
    file, attach to MICM with TS1, solve one timestep. Compare against Python tutorial 14
    output.

**Verification:**
- Fortran test produces same results as Python tutorial 14 (CAM Cloud Chemistry)
- Config-driven MIAM setup from Fortran works (JSON path → MIAM model → attached to MICM)
- `pkg-config musica-fortran` includes MIAM headers/libs
- MPAS `MUSICA=true` still builds with the updated MUSICA library

---

## Phase 8: Runtime Aerosol Representation Configuration from MIAM

*Mirrors Phase 3 but for aerosol/condensed-phase fields. MIAM defines representations
(e.g., `UniformSection`) and condensed-phase species. These are dynamically allocated
at runtime.*

40. **Query MIAM for representations and species** — At init-time:
    - Enumerate condensed-phase representations (sections/modes) and properties
      (radius range, phase name)
    - Enumerate aqueous/condensed species within each representation
    - Allocate MPAS fields for aerosol number concentrations and mass mixing ratios
      per representation
41. **Generic cloud coupling** — Map MPAS cloud fraction, LWC, droplet radius to MIAM
    representations based on representation metadata, not hardcoded.
42. **Representation ↔ MPAS field mapping** — Generic: if MIAM config defines 2 cloud
    sections instead of 1, MPAS fields adjust automatically.
43. **Test with minimal aerosol config** — 1 representation, 1 condensed species.
    Validates allocation machinery.
44. **Switching test** — Change MIAM config (1 section → 2 sections) without
    recompilation.

**Verification:**
- Aerosol fields in output NetCDF match the MIAM representation definition
- Switching MIAM configs changes the field set without recompiling
- Trivial aerosol config runs without crashes
- No gas-phase regressions from earlier phases

---

## Phase 9: CAM Cloud Chemistry in MPAS

45. **Extend chemistry driver for MIAM** — Create MIAM model at init from config file
    (using Phase 7 config-driven setup); attach via `AddExternalModel`; switch to
    `RosenbrockDAE4StandardOrder`.
46. **Couple cloud properties** — Via generic Phase 8 mapping. Skip aqueous chemistry
    in clear-sky cells.
47. **DAE constraint initialization** — Max iterations, tolerance for Henry's Law,
    dissociation equilibria, charge balance, mass conservation.
48. **Integration test** — TS1 + CAM Cloud Chemistry on 480-km mesh with real-data init;
    sulfate production in cloudy regions.

**Verification:**
- DAE solver converges in cloudy grid cells
- Sulfate production rates match tutorial 14 for equivalent conditions
- No crashes in clear-sky cells (zero LWC handled gracefully)
- 48-hour real-data run completes without solver failures

---

## Phase 10: Hardening & Validation

49. **End-to-end regression test** — Full config (TS1 + CAM Cloud + TUV-x TS1/TSMLT);
    store reference output.
50. **Multi-resolution testing** — Verify on 240-km mesh (10,242 cells).
51. **Documentation** — Build, configure, run, add new mechanisms. Docker-first.
52. **Performance scaling** — 240-km and 120-km mesh benchmarks.

**Verification:**
- Regression test detects bitwise differences in chemistry output
- 240-km mesh run completes with chemistry
- New team member builds and runs CheMPAS-A in < 1 hour (using Docker)

---

## Key Files

### MUSICA (existing)

| File | Role |
|------|------|
| `fortran/micm/micm.F90` | MICM Fortran API (`micm_t`, `solve()`); template for MIAM bindings |
| `fortran/tuvx/tuvx.F90` | TUV-x Fortran API (`tuvx_t`, `run()`) |
| `fortran/util.F90` | Shared types: `string_t`, `error_t`, `mappings_t` |
| `configs/v0/chapman/` | Chapman mechanism (Phase 2) |
| `configs/v0/analytical/` | Analytical test mechanism (Phase 3 switching test) |
| `configs/v1/ts1/ts1.json` | TS1 mechanism (Phase 6) |
| `configs/tuvx/ts1_tsmlt.json` | TUV-x TS1/TSMLT photolysis (Phase 5) |
| `src/miam/miam_builder.cpp` | MIAM C++ ↔ MICM integration reference |
| `docs/MIAM_INTEGRATION_PLAN.md` | Existing MIAM integration plan (Phases 0–5) |

### MUSICA (new — Phase 7)

| File | Role |
|------|------|
| `fortran/miam/miam.F90` | MIAM Fortran bindings |
| `src/miam/miam_c_interface.cpp` | MIAM C interface |
| `include/musica/miam/` | C interface headers |

### MPAS-A (new)

| File | Role |
|------|------|
| `docker/Dockerfile.chempas` | Build environment + MUSICA |
| `.devcontainer/devcontainer.json` | VS Code dev container |
| `src/core_atmosphere/chemistry/mpas_chemistry_driver.F90` | Central chemistry coupling |
| `src/core_atmosphere/chemistry/mpas_chemistry_species_mapping.F90` | Tracer ↔ MICM mapping |
| `src/core_atmosphere/chemistry/mpas_chemistry_allocate_scalars.F90` | Runtime species allocation (Phase 3) |
| `src/core_atmosphere/chemistry/mpas_chemistry_emissions.F90` | Stub emissions (Phase 4) |
| `src/core_atmosphere/chemistry/mpas_chemistry_dry_deposition.F90` | Stub dry deposition (Phase 4) |
| `src/core_atmosphere/chemistry/mpas_chemistry_wet_deposition.F90` | Stub wet deposition (Phase 4) |
| `src/core_atmosphere/chemistry/mpas_chemistry_aerosol_mapping.F90` | Aerosol ↔ MIAM mapping (Phase 8) |
| `test/` | Test harness, reference data |
| `scripts/download_mesh.sh` | Mesh downloader |

---

## Decisions

- **Repos**: MUSICA work on this dev branch; MPAS on `mattldawson/MPAS-Model`
- **Mesh**: 480-km (2,562 cells) for dev/CI; 240-km (10,242 cells) for validation
- **Build**: `gfortran` (GNU toolchain); Docker handles all dependencies
- **I/O**: SMIOL (no external PIO dep) for 480-km; PIO2 optional for larger meshes
- **Solver**: `RosenbrockStandardOrder` (Phases 2–6) →
  `RosenbrockDAE4StandardOrder` (Phase 9+)
- **Chemistry coupling**: operator-split, called after physics
- **Runtime species allocation** (Phase 3): CAM-MPAS `atm_allocate_scalars` pattern
- **Emissions/deposition** (Phase 4): stub parameterizations feeding MICM's
  `EMISSION`/`FIRST_ORDER_LOSS` reaction types; detailed inventories deferred
- **Runtime aerosol allocation** (Phase 8): query MIAM model at init
- **MIAM Fortran bindings + Mechanism Configuration** (Phase 7) block Phases 8–9
- **Short-lived radicals**: diagnosed locally, not advected (~50% tracer cost reduction)

## Further Considerations

1. **Chemistry timestep**: configurable `config_chemistry_dt` (default: physics timestep,
   ~15 min at 480-km). Separate from dynamics timestep.
2. **Species advection cost**: TS1 has ~70 species; only advect long-lived, diagnose
   short-lived radicals locally. Phase 3's classification system handles this.
3. **SMIOL vs PIO**: SMIOL eliminates an external dep for the 480-km mesh. At 120-km+,
   PIO2 parallel I/O becomes important.
4. **Mechanism Configuration for MIAM**: The existing MIAM integration plan
   (`docs/MIAM_INTEGRATION_PLAN.md`, Phase 4) covers extending the C++/Python config
   parser. Phase 7 of this plan adds the Fortran bindings for that parser so MPAS can
   use config-driven MIAM setup.
5. **Emissions refinement**: Phase 4 stubs use constant or simple config-driven rates.
   Future work: read CAMS/EDGAR inventories, diurnal/seasonal scaling, plume rise for
   elevated sources. The stub interface (MICM `EMISSION` parameters) remains the same —
   only the rate computation changes.
