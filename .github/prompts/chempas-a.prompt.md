---
description: "CheMPAS-A development plan: integrating MUSICA chemistry (MICM, TUV-x, MIAM) into MPAS-A with runtime-configurable species and aerosol representations. Use when working on MPAS chemistry integration, CheMPAS-A, or the MIAM Fortran bindings."
---

# CheMPAS-A: MUSICA Chemistry in MPAS-A

## Repositories

- **MUSICA**: development branch (this repo)
- **MPAS-A**: https://github.com/mattldawson/MPAS-Model

## Goal

Integrate TS1 gas-phase chemistry (MICM), CAM Cloud aqueous chemistry (MIAM), and
TUV-x photolysis into MPAS-A. Container-based local development on the 480-km
quasi-uniform mesh (2,562 cells). All phases use the Jablonowski-Williamson
baroclinic wave idealized test case (no GFS real-data dependency). Start with
Chapman chemistry (hardcoded species), then add **runtime species configuration** so
MPAS tracers are dynamically allocated from the MICM mechanism at init-time
(following the proven CAM-MPAS `atm_allocate_scalars` pattern). After MIAM Fortran
bindings and Mechanism Configuration updates, add **runtime aerosol representation
configuration** so aerosol fields are driven by the MIAM config. One build runs any
mechanism.

**Test strategy:** The JW case provides realistic T/P/density profiles, multi-level
atmosphere, multi-cell lat/lon coverage, and working tracer advection. For later
phases that need cloud water (MIAM) or aerosol optical depth (TUV-x), we prescribe
synthetic fields in the JW initialization (e.g., a cloud layer at 700–850 hPa with
prescribed LWC). This avoids a GFS dependency while still exercising all code paths.

---

## Phase 0: Local Build & Test Infrastructure

1. **Create Containerfile** — Fedora with gfortran, OpenMPI, NetCDF, PNetCDF, CMake,
   and MUSICA-Fortran pre-installed via `pkg-config`. Multi-stage build:
   `deps` (libraries), `build` (compiled MPAS), `run` (lean runtime), `dev` (interactive).
2. **Create `devcontainer.json`** — VS Code dev container with `docker-compose.yml`.
3. **Mesh download script** — Fetch the 480-km mesh (`x1.2562.tar.gz`, 1.5 MB) and static
   file (`x1.2562_static.tar.gz`, 1.0 MB) from `www2.mmm.ucar.edu` into a gitignored
   `data/` directory. *(parallel with 1)*
4. **Build vanilla MPAS-A** — Verify `make gnu CORE=init_atmosphere` and
   `make gnu CORE=atmosphere` inside the container. Run the Jablonowski-Williamson
   baroclinic wave idealized test case on the 480-km mesh with 2 MPI ranks.
   *(depends on 1, 3)*
5. **Test harness** — Shell scripts: build verification, JW execution,
   `ncdump`-based output validation. *(depends on 4)*
6. **CI pipeline** — GitHub Actions: build container → compile → test on every PR.
   480-km mesh, 2 MPI ranks. *(depends on 5)*

### Verify

See `verification/README.md` for setup, then open
[`verification/phase00_jw_baseline.ipynb`](verification/phase00_jw_baseline.ipynb).
Plots surface pressure, potential temperature, zonal/meridional wind maps and
profiles, time evolution statistics, and runs automated pass/fail checks.

---

## Phase 1: Passive Tracer Infrastructure

8. **Add passive tracers to `Registry.xml`** — 2–3 test tracers advected by the dynamical
   core. Validates the tracer advection pathway before chemistry.
9. **Initialize tracers** — Gaussian blob + uniform background.
10. **Validate tracer advection** — Mass conservation, no negatives with monotonic limiter.
    Add output stream for tracer diagnostics.
11. **Tracer indexing unit test** — Standalone Fortran test for name ↔ index mapping
    (reused in Phase 3 species mapping).

### Verify

[`verification/phase01_tracers.ipynb`](verification/phase01_tracers.ipynb) —
Plots tracer mass conservation timeseries, spatial distribution maps, and
negative-value checks.

---

## Phase 2: Chapman Chemistry (Simplest MICM + TUV-x)

12. **Build with MUSICA** — `make gnu CORE=atmosphere MUSICA=true`. Container image
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

### Verify

[`verification/phase02_chapman.ipynb`](verification/phase02_chapman.ipynb) —
Plots O3 diurnal cycle (map + timeseries), odd-oxygen (O+O1D+O3) mass
conservation, and single-column comparison against tutorial 10.

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

### Verify

[`verification/phase03_runtime_species.ipynb`](verification/phase03_runtime_species.ipynb) —
Compares species sets between Chapman and analytical mechanisms,
and verifies bitwise-identical output when rerunning Chapman.

---

## Phase 4: Emissions & Deposition Stubs

*Provides minimum-viable emissions, dry deposition, and wet deposition so that TS1
(Phase 6) produces reasonable-ish results. These are simple parameterizations that
feed MICM's existing `EMISSION` and `FIRST_ORDER_LOSS` reaction types — not
detailed inventories. Detailed emissions from files and refined deposition handling
are deferred to future work.*

__Notes from the human__:
* Remember that we'll be adding the TS1 mechanism soon, so instead of hardcoding emissions
  and deposition for specific species, it might make sense to create a handful of generic
  2/3-D profiles for emissions and deposition that can be configured for specific species
  in the config files for each mechanism.
* Because we're using an idealized MPAS configuration, I don't think we have alnd-use types
  and things like this available, and the goal here is not really to create a realistic
  configuration for CheMPAS-A, but just to verify our implementation of each Phase, For
  the emissions and deposition, we want to be able to see downstream effects of gases
  being emitted and deposited as air parcels move around the domain.
* You can add some extra species to the Chapman mechanism for verifying the emissions
  and deposition logic are working as expected. They dont have to correspond to real-
  world species. foo, bar, etc. are fine.

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

### Verify

[`verification/phase04_emissions_deposition.ipynb`](verification/phase04_emissions_deposition.ipynb) —
Plots emitted species accumulation, deposited species decay,
and mass budget closure timeseries.

---

## Phase 4b: Mid-Development Review
*At this point, let's step back and look at what's been implemented, and come up with a plan
to address any issues with what we've done so far. We can review the plan as a team and
determine whether and how to implement any changes.

Things to consider:
- Now that species are configured at run-time, we should have no references in code to
specific chemical species. Please review all the code we've added and ensure this is the
case. If there are any hard-coded species, describe how that part of the code can be refactored to be mechanism-agnostic.
- The mapping between MPAS, TUV-x and MICM grid cells is very complex and error prone. Please
develop a plan for adding a Jupyter notebook as part of this mini-phase that can convince us
that the mapping is precisely correct in a quantitative way.
- Review our development thus far for any debug code, anything that looks like it could be
a bug, and anything that just needs to be cleaned up. Develop as part of this mini-phase a plan for how to address any issues you find.

---

## Phase 5: TUV-x TS1/TSMLT Photolysis
*Remember: No species or mechanism names anywhere in the Fortran code.*

25. **Upgrade TUV-x config** — Switch to `configs/tuvx/ts1_tsmlt.json` for the full
    photolysis rate set needed by TS1.
26. **Couple MPAS vertical coordinate to TUV-x** — Map terrain-following height →
    TUV-x altitude grid. Update O2/O3/air profiles from MPAS state each chemistry step.
27. **Implement alias mapping** — Route TUV-x reaction labels (`jo2_a`, `jo3->O1D`, …)
    to MICM photolysis parameters using the aliasing config.
28. **Regression test** — Verify Chapman photolysis rates unchanged after TUV-x upgrade.

### Verify

[`verification/phase05_tuvx_photolysis.ipynb`](verification/phase05_tuvx_photolysis.ipynb) —
Plots photolysis rate vertical profiles and regression comparison
against Phase 2 reference output.

---

## Phase 6: TS1 Mechanism (~70 species, ~200 reactions)
*Remember: No species or mechanism names anywhere in the Fortran code.*

__Notes from the human__:
* Let's make sure that we have enough precurors in place (either by initial conditions or
emissions) that almost all species have non-zero concentrations after a few chemistry steps.
I think there are 4 or 5 species that don't acutally participate in anything, so they can
be zero.
* There's also surface reactions on aerosols, but we don't have aerosols yet, so
let's add an aerosol stub that let's us provide an initial vertical aerosol profile in a csv
file, just like we set up initial conditions for gas-phase species. This will let us verify that surface reactions are working as expected.
* Mapping of photolysis reactions from TUV-x to MICM will be important to verify carefully.
Let's make sure that we have a plan for quantitatively verifying this mapping in the Jupyter notebook for this phase.
* Before starting on step 26. let's reorganize the chemistry data folder. We essentially want one folder
under `chemistry_data/` for each mechanism (analytical, chapman, ts1, etc.) that contains all the config files needed for that mechanism (MICM mechanism JSON, TUV-x photolysis config, emissions/deposition config, etc.). This will make it easier to manage the growing number of config files and ensure that when we switch mechanisms we have all the right configs in place. And we want to make sure we don't have any hard-coded
paths to these folders anywhere in the Fortran code. The path to the mechanism folder should be passed in as an argument (e.g., `config_chemistry_config_path`) and then all config files for that mechanism should be read relative to that path.

29. **Switch config to TS1** — Point `config_chemistry_config_path` at
    `configs/v1/ts1/ts1.json`. Thanks to Phase 3, this should "just work" — runtime
    allocation handles the new species set automatically.
30. **Initial conditions** — Default from `configs/v1/ts1/initial_conditions.csv`.
    Uniform background profiles adequate for JW test.
31. **Populate emissions/deposition configs for TS1** — Fill in the Phase 4 stub config
    files with order-of-magnitude rates for TS1 species (NOx surface emissions, O3/NO2
    dry deposition velocities, soluble species wet scavenging coefficients). Values from
    literature or CAM-chem defaults — good enough for plausible concentrations, not
    publication-quality.
32. **Performance profiling** — 480-km mesh (2,562 × 55 levels ≈ 141K points). Target:
    < 2 sec/chemistry step. GPTL timers.
33. **Integration test** — 24-hour forecast; O3, NO, NO2, CO diurnal patterns; mass budgets.

### Verify

[`verification/phase06_ts1.ipynb`](verification/phase06_ts1.ipynb) —
Plots O3/NO/NO2/CO diurnal cycles, spatial O3 maps, and checks
all ~70 species for negative values.

---

## Phase 7: MIAM Fortran Bindings & Mechanism Configuration in MUSICA

*Blocker for Phases 8–9. Work happens in this MUSICA repo.*

34. **Develop Detailed Plan** - Work with the human to create a detailed plan for this
    phase. Include in the plan the proposed updated Mechanism Configuration v1 schema. Put
    the plan in a markdown doc alongside this one. Below are notes and things to include:
35. **MUSICA Fortran API** - We shouldn't actually have to update the MUSICA Fortran
    API to allow MIAM to be used from Fortran. The Fortran API for all components requires
    that most of the configuration be in config files. As long as Mechanism Configuration
    is update to include MIAM configuration, all we should have to do is update the MICM
    constructor to build with MIAM. API users will see MIAM contributions just as additional
    state variables and parameters. Convince yourself that this is true, and ask the human if you have any doubts.
36. **Add DAE solver types to Fortran API** — `RosenbrockDAE4`, `RosenbrockDAE6` in
    `SolverType` enum. `micm_t` constructor accepts external models.
37. **Extend Mechanism Configuration for MIAM** — If not already done by this point
    (see `docs/MIAM_INTEGRATION_PLAN.md` Phase 4):
    - Extend MechanismConfiguration C++ library to parse MIAM aerosol model definitions
      (representations, processes, constraints) from JSON/YAML
    - Add Fortran bindings for the extended parser so that config-driven MIAM setup
      works from Fortran (not just programmatic API)
    - Wire parsed MIAM configs into solver construction through the C/Fortran interface
    - This enables MPAS to set up MIAM from a config file path, just like it does for
      MICM — no programmatic Fortran construction of MIAM objects needed
38. **Fortran unit test** — CAM Cloud Chemistry config: create MIAM model from config
    file, attach to MICM with TS1, solve one timestep. Compare against Python tutorial 14
    output.

### Verify

[`verification/phase07_miam_bindings.ipynb`](verification/phase07_miam_bindings.ipynb) —
Compares MIAM Fortran unit test output against Python tutorial 14
(CAM Cloud Chemistry) reference values.

---

## Phase 8: CAM Cloud Chemistry in MPAS
*Remember: No species or mechanism names anywhere in the Fortran code.*

39. **Extend chemistry driver for MIAM** — Create a MICM solver that includes
    MIAM cloud chemistry (the CAM cloud chemistry configuration) and a
    `RosenbrockDAE4StandardOrder` solver.
40. **Couple cloud properties** — I think MPAS might just have cloud liquid water content
    available. If so, we can just look for a MICM solver variable named:
    `CLOUD.MODE.AQUEOUS.H2O` and when it's present set it using the cloud liquid water,
    and if it's not there, assume there is no cloud chemistry. Other cloud species can be
    treated like tracers to advect just like gas-phase species. Better yet, make this mapping
    from cloud liquid water part of the configuration, so there's not hard-coded references
    species species in the code. 
41. **DAE constraint initialization** — Max iterations, tolerance for Henry's Law,
    dissociation equilibria, charge balance, mass conservation.
42. **Integration test** — TS1 + CAM Cloud Chemistry on 480-km JW mesh with prescribed
    cloud layer (700–850 hPa, LWC=0.3 g/m³); sulfate production in cloudy regions.

### Verify

[`verification/phase08_cloud_chemistry.ipynb`](verification/phase08_cloud_chemistry.ipynb) —
Plots sulfate production timeseries in cloud layer, compares cloudy
vs clear-sky cells, and checks DAE solver convergence.

---

## Phase 9: Hardening & Validation

43. **End-to-end regression test** — Full config (TS1 + CAM Cloud + TUV-x TS1/TSMLT);
    store reference output.
44. **Multi-resolution testing** — Verify on 240-km mesh (10,242 cells).
45. **Documentation** — Build, configure, run, add new mechanisms. Container-first.
46. **Performance scaling** — 240-km and 120-km mesh benchmarks.

### Verify

[`verification/phase09_regression.ipynb`](verification/phase09_regression.ipynb) —
Runs bitwise regression against stored reference, plots per-species
differences, and compares 480-km vs 240-km mesh statistics.

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

### MPAS-A (new)

| File | Role |
|------|------|
| `docker/Containerfile` | Build environment + MUSICA |
| `.devcontainer/devcontainer.json` | VS Code dev container |
| `src/core_atmosphere/chemistry/mpas_chemistry_driver.F90` | Central chemistry coupling |
| `src/core_atmosphere/chemistry/mpas_chemistry_species_mapping.F90` | Tracer ↔ MICM mapping |
| `src/core_atmosphere/chemistry/mpas_chemistry_allocate_scalars.F90` | Runtime species allocation (Phase 3) |
| `src/core_atmosphere/chemistry/mpas_chemistry_emissions.F90` | Stub emissions (Phase 4) |
| `src/core_atmosphere/chemistry/mpas_chemistry_dry_deposition.F90` | Stub dry deposition (Phase 4) |
| `src/core_atmosphere/chemistry/mpas_chemistry_wet_deposition.F90` | Stub wet deposition (Phase 4) |
| `test/` | Test harness, reference data |
| `scripts/download_mesh.sh` | Mesh downloader |

---

## Decisions

- **Repos**: MUSICA work on this dev branch; MPAS on `mattldawson/MPAS-Model`
- **Mesh**: 480-km (2,562 cells) for dev/CI; 240-km (10,242 cells) for validation
- **Build**: `gfortran` (GNU toolchain); containers handle all dependencies
- **I/O**: SMIOL (no external PIO dep) for 480-km; PIO2 optional for larger meshes
- **Solver**: `RosenbrockStandardOrder` (Phases 2–6) →
  `RosenbrockDAE4StandardOrder` (Phase 9+)
- **Chemistry coupling**: operator-split, called after physics
- **Runtime species allocation** (Phase 3): CAM-MPAS `atm_allocate_scalars` pattern
- **Emissions/deposition** (Phase 4): stub parameterizations feeding MICM's
  `EMISSION`/`FIRST_ORDER_LOSS` reaction types; detailed inventories deferred
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
