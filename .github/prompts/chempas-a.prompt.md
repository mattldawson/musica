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

### Build / Run / Verify

```bash
# Pre-requisite: podman (or docker)
# 1. Pre-checkout physics externals and lookup tables (one-time, on host)
podman run --rm -v "$(pwd):/mpas:Z" -w /mpas/src/core_atmosphere/physics \
    chempas-deps ./../tools/manage_externals/checkout_externals --externals ./../Externals.cfg
podman run --rm -v "$(pwd):/mpas:Z" -w /mpas/src/core_atmosphere/physics/physics_wrf/files \
    chempas-deps bash -c 'curl -sSf -o MPAS-Data.tar.gz \
    https://codeload.github.com/MPAS-Dev/MPAS-Data/tar.gz/v8.2 && \
    tar xzf MPAS-Data.tar.gz && mv MPAS-Data-8.2/atmosphere/physics_wrf/files/* . && \
    rm -rf MPAS-Data.tar.gz MPAS-Data-8.2'

# 2. Build all container stages
podman build -f docker/Containerfile --target deps  -t chempas-deps  .
podman build -f docker/Containerfile --target build -t chempas-build .
podman build -f docker/Containerfile --target run   -t chempas-run   .

# 3. Download mesh data
bash scripts/download_data.sh data

# 4. Run JW baroclinic wave test (2 MPI ranks, 1-day simulation)
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2

# Expected: "JW baroclinic wave test PASSED", output.nc ~ 33 MB
# Verify output contains expected fields:
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'ncdump -h /mpas/data/jw_480km/output.nc | grep -E "surface_pressure|theta|uReconstructZonal"'
```

**Done when:** JW test prints `PASSED`, output.nc exists with pressure/theta/wind fields.

---

## Phase 1: Passive Tracer Infrastructure

8. **Add passive tracers to `Registry.xml`** — 2–3 test tracers advected by the dynamical
   core. Validates the tracer advection pathway before chemistry.
9. **Initialize tracers** — Gaussian blob + uniform background.
10. **Validate tracer advection** — Mass conservation, no negatives with monotonic limiter.
    Add output stream for tracer diagnostics.
11. **Tracer indexing unit test** — Standalone Fortran test for name ↔ index mapping
    (reused in Phase 3 species mapping).

### Build / Run / Verify

```bash
# Rebuild with tracer changes
podman build -f docker/Containerfile --target build -t chempas-build .

# Run JW test (48 hours, 2 MPI ranks) with passive tracers
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2

# Verify tracers exist in output
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'ncdump -h /mpas/data/jw_480km/output.nc | grep -E "tracer[0-9]+"'

# Verify mass conservation (total tracer mass at t=0 vs t=final)
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'python3 -c "
import netCDF4 as nc
ds = nc.Dataset(\"/mpas/data/jw_480km/output.nc\")
for v in ds.variables:
    if v.startswith(\"tracer\"):
        print(f\"{v}: t=0 sum={ds[v][0].sum():.6e}, t=-1 sum={ds[v][-1].sum():.6e}\")
"'
```

**Done when:** Tracers appear in output, mass conserved to machine precision over 48 hours.

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

### Build / Run / Verify

```bash
# Rebuild with MUSICA=true
podman build -f docker/Containerfile --target build -t chempas-build \
    --build-arg MUSICA=true .

# Run JW + Chapman (24 hours)
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2  # config_chemistry_enabled=true, config=chapman

# Verify O3 diurnal cycle
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'python3 -c "
import netCDF4 as nc
ds = nc.Dataset(\"/mpas/data/jw_480km/output.nc\")
o3 = ds[\"o3\"]  # shape: (time, cell, level)
for t in range(ds.dimensions[\"Time\"].size):
    print(f\"t={t}: O3 min={o3[t].min():.3e} max={o3[t].max():.3e} mean={o3[t].mean():.3e}\")
"'

# Compare single-column against tutorial 10 reference
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'python3 scripts/compare_chapman_tutorial10.py'
```

**Done when:** O3 shows diurnal variation, single-column matches tutorial 10 within solver tolerance, O+O1D+O3 mass conserved.

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

### Build / Run / Verify

```bash
# No rebuild needed — same binary, different config

# Run with Chapman config (should produce identical output to Phase 2)
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2  # config_chemistry_config_path=configs/chapman

# Switch to analytical test mechanism — NO RECOMPILATION
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2  # config_chemistry_config_path=configs/analytical

# Verify output has mechanism-specific species names
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'ncdump -h /mpas/data/jw_480km/output.nc | grep -E "float .*(nCells, nVertLevels)"'
```

**Done when:** Chapman gives identical results to Phase 2; analytical mechanism runs without recompilation; species names in NetCDF match mechanism JSON.

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

### Build / Run / Verify

```bash
# Rebuild with emissions/deposition modules
podman build -f docker/Containerfile --target build -t chempas-build .

# Run Chapman — stubs should be no-ops (no EMISSION/FIRST_ORDER_LOSS in Chapman)
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2  # config_chemistry_config_path=configs/chapman

# Run test mechanism with 1 emitted + 1 deposited species
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2  # config_chemistry_config_path=configs/emis_dep_test

# Verify mass budget
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'python3 scripts/check_mass_budget.py /mpas/data/jw_480km/output.nc'
```

**Done when:** Chapman results unchanged from Phase 2/3; test mechanism reaches steady state; mass budget (emissions − deposition ± chemistry) closes to solver tolerance.

---

## Phase 5: TUV-x TS1/TSMLT Photolysis

25. **Upgrade TUV-x config** — Switch to `configs/tuvx/ts1_tsmlt.json` for the full
    photolysis rate set needed by TS1.
26. **Couple MPAS vertical coordinate to TUV-x** — Map terrain-following height →
    TUV-x altitude grid. Update O2/O3/air profiles from MPAS state each chemistry step.
27. **Implement alias mapping** — Route TUV-x reaction labels (`jo2_a`, `jo3->O1D`, …)
    to MICM photolysis parameters using the aliasing config.
28. **Regression test** — Verify Chapman photolysis rates unchanged after TUV-x upgrade.

### Build / Run / Verify

```bash
# Rebuild with updated TUV-x config
podman build -f docker/Containerfile --target build -t chempas-build .

# Run JW + Chapman with TS1/TSMLT TUV-x
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2

# Compare photolysis rates against Phase 2 reference
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'python3 scripts/compare_photolysis_rates.py \
    /mpas/data/jw_480km/output.nc /mpas/data/reference/phase2_output.nc'
```

**Done when:** Photolysis rates for Chapman reactions match Phase 2 within tolerance; alias mapping verified via diagnostics.

---

## Phase 6: TS1 Mechanism (~70 species, ~200 reactions)

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

### Build / Run / Verify

```bash
# Config-only switch to TS1 — NO recompilation needed
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2  # config_chemistry_config_path=configs/ts1

# Verify no negative species and basic diurnal patterns
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'python3 -c "
import netCDF4 as nc
ds = nc.Dataset(\"/mpas/data/jw_480km/output.nc\")
for sp in [\"O3\", \"NO\", \"NO2\", \"CO\"]:
    if sp in ds.variables:
        vals = ds[sp][:]
        print(f\"{sp}: min={vals.min():.3e} max={vals.max():.3e} any_negative={bool((vals<0).any())}\")
"'

# Performance timing (GPTL output)
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'grep chemistry /mpas/data/jw_480km/timing.* || echo "Add GPTL timers"'
```

**Done when:** TS1 runs with config-only switch (no code changes); 48-hour run stable; no negative species; O3/NO2 diurnal patterns visible; < 2 sec/step.

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

### Build / Run / Verify

```bash
# Build updated MUSICA with MIAM Fortran bindings
cd /path/to/musica && mkdir -p build && cd build
cmake .. -DMUSICA_BUILD_FORTRAN_INTERFACE=ON -DMUSICA_BUILD_CARMA=OFF
make -j$(nproc) && ctest --output-on-failure

# Specifically run MIAM Fortran test
ctest -R miam_fortran --output-on-failure

# Verify pkg-config exposes MIAM
pkg-config --libs musica-fortran | grep -i miam

# Rebuild MPAS container with updated MUSICA
cd /path/to/MPAS-Model
podman build -f docker/Containerfile --target build -t chempas-build .
```

**Done when:** MIAM Fortran unit test matches tutorial 14 output; `pkg-config musica-fortran` includes MIAM; MPAS still builds with `MUSICA=true`.

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

### Build / Run / Verify

```bash
# Rebuild MPAS with aerosol mapping support
podman build -f docker/Containerfile --target build -t chempas-build .

# Run with minimal aerosol config (1 representation, 1 species)
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2  # config includes minimal MIAM

# Verify aerosol fields in output
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'ncdump -h /mpas/data/jw_480km/output.nc | grep -E "aerosol|section|representation"'

# Switch to 2-section config — NO recompilation
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2  # config includes 2-section MIAM

# Verify new fields appeared
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'ncdump -h /mpas/data/jw_480km/output.nc | grep -c section'
```

**Done when:** Aerosol fields in NetCDF match MIAM config; switching 1→2 sections works without recompilation; gas-phase results unchanged.

---

## Phase 9: CAM Cloud Chemistry in MPAS

45. **Extend chemistry driver for MIAM** — Create MIAM model at init from config file
    (using Phase 7 config-driven setup); attach via `AddExternalModel`; switch to
    `RosenbrockDAE4StandardOrder`.
46. **Couple cloud properties** — Via generic Phase 8 mapping. Skip aqueous chemistry
    in clear-sky cells.
47. **DAE constraint initialization** — Max iterations, tolerance for Henry's Law,
    dissociation equilibria, charge balance, mass conservation.
48. **Integration test** — TS1 + CAM Cloud Chemistry on 480-km JW mesh with prescribed
    cloud layer (700–850 hPa, LWC=0.3 g/m³); sulfate production in cloudy regions.

### Build / Run / Verify

```bash
# Rebuild with MIAM integration
podman build -f docker/Containerfile --target build -t chempas-build .

# Run TS1 + CAM Cloud Chemistry (JW with prescribed cloud layer)
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2  # config includes TS1 + MIAM cloud chem

# Verify sulfate production in cloudy cells
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'python3 -c "
import netCDF4 as nc
ds = nc.Dataset(\"/mpas/data/jw_480km/output.nc\")
if \"SO4\" in ds.variables:
    so4 = ds[\"SO4\"][:]
    print(f\"SO4: t=0 max={so4[0].max():.3e}, t=-1 max={so4[-1].max():.3e}\")
    print(f\"SO4 increased: {so4[-1].max() > so4[0].max()}\")
"'

# Verify DAE solver convergence (no failures in log)
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'grep -i "solver.*fail\|DAE.*error" /mpas/data/jw_480km/log.atmosphere.* || echo "No solver failures"'
```

**Done when:** DAE solver converges in cloudy cells; sulfate production visible; no crashes in clear-sky cells; 48-hour run completes.

---

## Phase 10: Hardening & Validation

49. **End-to-end regression test** — Full config (TS1 + CAM Cloud + TUV-x TS1/TSMLT);
    store reference output.
50. **Multi-resolution testing** — Verify on 240-km mesh (10,242 cells).
51. **Documentation** — Build, configure, run, add new mechanisms. Container-first.
52. **Performance scaling** — 240-km and 120-km mesh benchmarks.

### Build / Run / Verify

```bash
# Full-stack build
podman build -f docker/Containerfile --target build -t chempas-build .

# Generate regression reference (TS1 + CAM Cloud + TUV-x TS1/TSMLT)
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2  # full config
cp data/jw_480km/output.nc data/reference/phase10_reference.nc

# Regression test: compare against reference
podman run --rm -v "$(pwd)/data:/mpas/data:Z" chempas-build \
    bash -c 'python3 scripts/regression_test.py \
    /mpas/data/jw_480km/output.nc /mpas/data/reference/phase10_reference.nc'

# 240-km mesh test
bash scripts/download_data.sh data 240km
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 4  # 240-km, 4 MPI ranks
```

**Done when:** Regression test detects bitwise differences; 240-km run completes; full build-to-run < 1 hour for a new user.

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
| `docker/Containerfile` | Build environment + MUSICA |
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
- **Build**: `gfortran` (GNU toolchain); containers handle all dependencies
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
