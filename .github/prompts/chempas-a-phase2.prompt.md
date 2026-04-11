---
description: "CheMPAS-A Phase 2 detailed plan: Chapman Chemistry integration (MICM + TUV-x) into MPAS-A. Use when implementing Phase 2 of the CheMPAS-A development plan."
---

# CheMPAS-A Phase 2: Chapman Chemistry (MICM + TUV-x)

## Reference

- **Overall plan**: [chempas-a.prompt.md](chempas-a.prompt.md)
- **CAM-SIMA integration**: https://github.com/ESCOMP/atmospheric_physics/tree/development/schemes/musica
- **Key reference files**:
  - `musica_ccpp.F90` — top-level register→init→run→final lifecycle
  - `musica_ccpp_micm.F90` — MICM init/run with two-state residual pattern
  - `musica_ccpp_micm_util.F90` — strided copy between host model and MICM state
  - `musica_ccpp_tuvx.F90` — TUV-x run with per-column loop
  - `musica_ccpp_tuvx_height_grid.F90` — height grid conversion logic
  - `musica_ccpp_tuvx_gas_species.F90` — gas species profile setup

---

## 1. Architecture Overview

### 1.1 Data Flow (per chemistry timestep)
```
MPAS state (T, P, ρ, scalars)
    │
    ├─► TUV-x (per column)
    │     Input: height grid, T, ρ_air, O2, O3, surface albedo, SZA, cloud frac
    │     Output: photolysis_rate_constants(nVertLevels, n_photo_rxns)
    │
    └─► MICM (all grid cells at once, vectorized)
          Input: T, P, ρ_air, species concentrations (mol/m³), rate_parameters
          Output: updated species concentrations
    │
    ▼
Updated MPAS scalars (O, O1D, O3)
```

### 1.2 Operator Splitting
Chemistry is operator-split with dynamics, called from `physics_driver()` after
radiation (which provides solar zenith angle and surface albedo). The call sequence:

```
atm_do_timestep
├── physics_driver()
│   ├── MPAS_to_physics()
│   ├── [radiation, PBL, convection, ...]
│   ├── chemistry_timestep()         ← NEW: TUV-x → MICM
│   └── physics_to_MPAS()           ← tendencies already applied inside chemistry_timestep
├── atm_timestep()                   ← dynamics + scalar transport
```

Chemistry directly updates `scalars` at time level 1 (current state) rather than
accumulating mass-weighted tendencies into `tend_scalars`. This is because:
1. MICM solves the full ODE system (not a tendency)
2. The updated scalars are then transported by the dynamical core in `atm_timestep()`
3. This matches the microphysics pattern which also directly updates state

### 1.3 Species Strategy for Phase 2

| Species | MICM Type     | MPAS Treatment | Notes |
|---------|---------------|----------------|-------|
| O3      | Reactive      | Advected scalar | Primary verification target |
| O       | Reactive      | Advected scalar | Short-lived radical, but advected for completeness |
| O1D     | Reactive      | Advected scalar | Very short-lived, effectively diagnostic |
| O2      | Constant      | Diagnosed | Fixed fraction of dry air: 0.2095 by volume |
| M       | Third body    | Diagnosed | Total air number density from ρ_air |
| N2      | Constant      | Diagnosed | Fixed fraction: 0.7808 by volume |

For Phase 2 (hardcoded Chapman), O, O1D, O3 are registered as advected scalars in
`Registry.xml`. O2 and N2 are computed from air density at each grid cell. M is
implicit in MICM's third-body handling.

---

## 2. Vertical Grid Alignment

### 2.1 MPAS Vertical Grid (surface-to-top)
```
k index:    1      2      3     ...    nVertLevels    nVertLevelsP1
            │      │      │              │              │
zgrid:      ├──────┼──────┼── ... ──────┼──────────────┤
           sfc                                        top
            ▲ interface heights (nVertLevelsP1 values)

midpoints:     *      *      *    ...    *
              k=1    k=2    k=3        k=nVertLevels
              zmid(k) = 0.5*(zgrid(k) + zgrid(k+1))
```

Both MPAS and TUV-x use surface-to-top ordering. Unlike CAM-SIMA (which is top-to-
surface and requires inversion), MPAS grids are already in the correct orientation
for TUV-x.

### 2.2 TUV-x Grid Mapping

TUV-x calculates photolysis rate constants at its grid **interfaces**. To get rates
at MPAS layer midpoints, we set TUV-x interfaces at MPAS midpoints:

```
MPAS (surface-to-top)                    TUV-x (surface-to-top)
                                         
zgrid(nVertLevelsP1) ──── top            interface(nVL+2) ──── top = zgrid(nVL+1)  [exo layer]
                                         midpoint(nVL+1)  ──── avg(zmid(nVL), zgrid(nVL+1))
zmid(nVertLevels) ════ mid nVL           interface(nVL+1) ──── zmid(nVL)
                                         midpoint(nVL)    ──── zgrid(nVL)
zmid(nVertLevels-1) ═══ mid nVL-1       interface(nVL)   ──── zmid(nVL-1)
     ⋮                                       ⋮
zmid(2) ════════════ mid 2               interface(3)     ──── zmid(2)
                                         midpoint(2)      ──── zgrid(2)
zmid(1) ════════════ mid 1               interface(2)     ──── zmid(1)
                                         midpoint(1)      ──── avg(zgrid(1), zmid(1)) [half layer]
zgrid(1) ─────────── surface             interface(1)     ──── zgrid(1) = surface
```

**Key**: TUV-x has `nVertLevels + 1` sections (= `nVertLevels + 2` interfaces).
The extra layers are:
- Bottom: half-layer between surface and lowest MPAS midpoint
- Top: layer above highest MPAS midpoint to model top (holds exo-atmosphere values)

This is the same approach as CAM-SIMA (`musica_ccpp_tuvx_height_grid.F90`), but
**without the inversion** since MPAS is already surface-to-top.

### 2.3 Grid Setup Code (simplified)

```fortran
! MPAS midpoint heights in km (surface-to-top, no inversion needed)
do k = 1, nVertLevels
  zmid_km(k) = 0.5e-3 * (zgrid(k,iCell) + zgrid(k+1,iCell))
end do
zint_km(1:nVertLevelsP1) = 1.0e-3 * zgrid(1:nVertLevelsP1,iCell)

! TUV-x interfaces = MPAS midpoints (already in correct order)
tuvx_interfaces(1) = zint_km(1)                        ! surface
tuvx_interfaces(2:nVertLevels+1) = zmid_km(1:nVertLevels)  ! MPAS midpoints
tuvx_interfaces(nVertLevels+2) = zint_km(nVertLevelsP1) ! model top

! TUV-x midpoints
tuvx_midpoints(1) = 0.5 * (tuvx_interfaces(1) + tuvx_interfaces(2))  ! half layer
do k = 2, nVertLevels
  tuvx_midpoints(k) = zint_km(k)  ! MPAS interfaces become TUV-x midpoints
end do
tuvx_midpoints(nVertLevels+1) = 0.5 * (tuvx_interfaces(nVertLevels+1) &
                                       + tuvx_interfaces(nVertLevels+2))
```

### 2.4 Photolysis Rate Mapping Back to MPAS

TUV-x returns `photolysis_rate_constants(nVertLevels+2, n_photo_rxns)` on its
interface grid. Extract rates at TUV-x interfaces 2 through nVertLevels+1 to get
rates at MPAS midpoints (layers 1:nVertLevels).

---

## 3. MICM State Copy Strategy

### 3.1 State Layout

MICM uses flat 1-D arrays with configurable strides:
```
state%conditions(i_grid_cell)        ! T, P, air_density per grid cell
state%concentrations(flat_index)     ! species in mol/m³
state%rate_parameters(flat_index)    ! photolysis rates + any user-defined rates
```

Strides for concentrations: `state%species_strides%grid_cell` and
`state%species_strides%variable`.

For a single column with `nVertLevels` grid cells:
```
flat_index = (i_cell - 1) * species_strides%grid_cell
           + (i_species - 1) * species_strides%variable + 1
```

### 3.2 Unit Conversions

**MPAS → MICM (before solve)**:
```
concentration [mol/m³] = mmr [kg/kg] × ρ_air [kg/m³] / M_w [kg/mol]
```

**MICM → MPAS (after solve)**:
```
mmr [kg/kg] = concentration [mol/m³] × M_w [kg/mol] / ρ_air [kg/m³]
```

Where `M_w` is molar mass and `ρ_air` is **dry** air density (= `zz(k,i) * rho_zz(k,i)`).

### 3.3 Batched Multi-Column Solving

Following the CAM-SIMA pattern in `musica_ccpp_micm.F90`, we batch **all** grid
cells into MICM at once: `n_grid_cells = nCellsSolve × nVertLevels`. MICM's
internal vectorization groups grid cells into fixed-size sets (set at MUSICA build
time). If the total count is not evenly divisible, we use a second `state_t` for
the residual — exactly as CAM-SIMA does with `state_1` / `state_2`.

This avoids per-column solve overhead and lets MICM's SIMD-grouped matrix layout
work optimally from the start.

```fortran
! Initialization: allocate states for the full mesh
n_grid_cells = nCellsSolve * nVertLevels
max_cells = micm%get_maximum_number_of_grid_cells()
size_1 = min(n_grid_cells, max_cells)
size_2 = mod(n_grid_cells - size_1, max_cells)
state_1 => micm%get_state(size_1, error)
if (size_2 > 0) state_2 => micm%get_state(size_2, error)
```

### 3.4 Copy Pattern

Following the CAM-SIMA `musica_ccpp_micm_util.F90` pattern, but adapted for MPAS
array layouts. Key differences from CAM-SIMA:
- MPAS scalars: `scalars(index, nVertLevels, nCells)` (scalar index is leading dim)
- CAM-SIMA constituents: `constituents(column, layer, constituent)` (column leading)
- Grid cells are flattened: MPAS cell `iCell`, level `k` → grid cell index
  `(iCell-1)*nVertLevels + k`

```fortran
! Copy MPAS state → MICM state (all cells batched)
i_cell = 0
do iCell = 1, nCellsSolve
  do k = 1, nVertLevels
    i_cell = i_cell + 1
    state%conditions(i_cell)%temperature = temperature(k, iCell)
    state%conditions(i_cell)%pressure = pressure(k, iCell)
    state%conditions(i_cell)%air_density = rho_dry(k, iCell)

    ! Species: mmr → mol/m³
    do s = 1, n_species
      flat_idx = (i_cell-1) * state%species_strides%grid_cell &
               + (micm_species_index(s)-1) * state%species_strides%variable + 1
      state%concentrations(flat_idx) = scalars(mpas_species_index(s), k, iCell) &
                                     * rho_dry(k, iCell) / molar_mass(s)
    end do

    ! Rate parameters (photolysis rates from TUV-x)
    do r = 1, n_rate_params
      flat_idx = (i_cell-1) * state%rate_parameters_strides%grid_cell &
               + (r-1) * state%rate_parameters_strides%variable + 1
      state%rate_parameters(flat_idx) = photo_rates(k, iCell, r)
    end do
  end do
end do
```

The solve loop iterates over states in batches of `state_1_size`, exactly matching
the `musica_ccpp_micm.F90` pattern:

```fortran
state_1_size = state_1%number_of_grid_cells
do i_state = 1, ceiling(real(n_grid_cells) / state_1_size)
  state_size = min(n_grid_cells - (i_state-1)*state_1_size, state_1_size)
  if (state_size == state_1_size) then
    state => state_1
  else
    state => state_2
  end if
  offset = (i_state - 1) * state_1_size
  call update_micm_state(state, offset, ...)
  call micm%solve(dt, state, solver_state, solver_stats, error)
  call state%update_references(error)  ! C++ may swap pointers
  call extract_from_micm_state(state, offset, ...)
end do
```

---

## 4. Solar Geometry

### 4.1 Solar Zenith Angle

MPAS physics radiation already computes `coszs` (cosine of solar zenith angle) per
cell. This is available in the `diag_physics` pool. We convert:
```fortran
sza_radians = acos(max(min(coszs(iCell), 1.0), -1.0))
```

If `physics_suite = 'none'` (our JW test case), radiation is not active and `coszs`
is not computed. We must compute SZA ourselves from:
- `latCell(iCell)`, `lonCell(iCell)` — cell center coordinates (radians)
- Solar declination from day of year
- Hour angle from simulation time

### 4.2 Earth-Sun Distance

Compute from day of year using standard orbital formula, or use a constant
(1 AU) for the JW idealized case.

### 4.3 Surface Albedo

For JW idealized case (no land surface model): use a constant (e.g., 0.1 over
ocean). When real physics is active, use `sfc_albedo` from the `diag_physics` pool.

---

## 5. Files to Create

### 5.1 `src/core_atmosphere/chemistry/mpas_chemistry_driver.F90`
Central chemistry coupling module. Contains:
- `chemistry_init(domain)` — Initialize MICM + TUV-x from config paths
- `chemistry_timestep(block, mesh, state, diag, dt)` — Per-block chemistry
- `chemistry_finalize()` — Cleanup

### 5.2 `src/core_atmosphere/chemistry/mpas_chemistry_tuvx.F90`
TUV-x interface module. Contains:
- `tuvx_setup(config_path, nVertLevels)` — Create TUV-x instance with height/wavelength grids and profile/radiator handles
- `tuvx_run_column(iCell, ...)` — Set height grid, T, species profiles, compute photolysis rates for one column
- `tuvx_cleanup()` — Deallocate TUV-x resources

### 5.3 `src/core_atmosphere/chemistry/mpas_chemistry_micm.F90`
MICM interface module. Contains:
- `micm_setup(config_path, n_grid_cells)` — Create MICM solver, allocate state_1/state_2 for batched solving (CAM-SIMA two-state pattern), extract species/rate ordering
- `micm_solve_batch(dt)` — Solve all grid cells in state_1/state_2 batches
- `micm_cleanup()` — Deallocate states and solver

### 5.4 `src/core_atmosphere/chemistry/mpas_chemistry_state.F90`
State copy utilities. Contains:
- `update_micm_from_mpas(state, scalars, ...)` — Copy MPAS → MICM with unit conversion
- `update_mpas_from_micm(state, scalars, ...)` — Copy MICM → MPAS with unit conversion

### 5.5 `src/core_atmosphere/chemistry/mpas_chemistry_utils.F90`
Utility functions. Contains:
- `compute_solar_zenith_angle(lat, lon, julday, ut_hours)` — SZA calculation
- `compute_earth_sun_distance(julday)` — Earth-Sun distance in AU

### 5.6 Config files
- **Chapman MICM config**: use existing `configs/v0/chapman/` (reactions.json + species.json).
  MICM photolysis rate parameter names: `jO2`, `jO3->O`, `jO3->O1D`.
- **TUV-x config**: create a minimal Chapman-only TUV-x config at
  `configs/tuvx/chapman.json`, extracting just the Chapman photolysis reaction
  configurations from `configs/tuvx/ts1_tsmlt.json`. This keeps TUV-x focused on
  only the 3 reactions Chapman needs, and lets us combine O2 channels directly:
  - `jO2` — combine the `jo2_a` + `jo2_b` cross-section/quantum-yield configs
    into a single reaction (O2 + hv → 2O). The TS1/TSMLT config splits these
    into separate Lyman-α/Schumann-Runge continuum (`jo2_a`) and
    Schumann-Runge bands/Herzberg continuum (`jo2_b`) channels; combining
    them gives us the total jO2 that MICM expects as a single rate parameter.
  - `jO3->O1D` — from `jo3_a` (O3 + hv → O2 + O(1D))
  - `jO3->O` — from `jo3_b` (O3 + hv → O2 + O(3P))
  The wavelength grid, O2 absorption band settings, and cross-section data files
  are all reused from `configs/tuvx/data/`. Only the photolysis reaction list and
  top-level config structure need to be written.
- **Photolysis-to-MICM mapping config**: with 3 TUV-x reactions matching the 3 MICM
  photolysis rate parameters by name, the mapping is 1:1:
  ```json
  [
    {"tuvx_name": "jO2",       "micm_name": "PHOTO.jO2"},
    {"tuvx_name": "jO3->O1D",  "micm_name": "PHOTO.jO3->O1D"},
    {"tuvx_name": "jO3->O",    "micm_name": "PHOTO.jO3->O"}
  ]
  ```
  Uses the `index_mappings_t` approach from CAM-SIMA (`musica_ccpp_tuvx.F90` init)
  with `MUSICA_INDEX_MAPPINGS_MAP_ANY` semantics.

---

## 6. Files to Modify

### 6.1 `src/core_atmosphere/Registry.xml`
Add Chapman species to the scalars var_array:
```xml
<var name="o3"  array_group="chapman" units="kg kg^{-1}"
     description="Ozone mass mixing ratio"/>
<var name="o"   array_group="chapman" units="kg kg^{-1}"
     description="Atomic oxygen mass mixing ratio"/>
<var name="o1d" array_group="chapman" units="kg kg^{-1}"
     description="Excited atomic oxygen O(1D) mass mixing ratio"/>
```

Also add to `scalars_tend` and `lbc_scalars`. Add namelist options:
```xml
<nml_option name="config_chemistry_enabled" type="logical" default_value=".false."/>
<nml_option name="config_chemistry_dt" type="real" default_value="0."/>
<nml_option name="config_micm_config_path" type="character" default_value=""/>
<nml_option name="config_tuvx_config_path" type="character" default_value=""/>
<nml_option name="config_tuvx_micm_mapping_path" type="character" default_value=""/>
```

### 6.2 `src/core_atmosphere/physics/mpas_atmphys_driver.F90`
Add call to `chemistry_timestep()` inside `physics_driver()`, gated by
`config_chemistry_enabled`. Place after radiation (so SZA is available when
physics is active).

### 6.3 `src/core_atmosphere/mpas_atm_core.F`
Add `chemistry_init()` call in `atm_core_init()` and `chemistry_finalize()` in
`atm_core_finalize()`.

### 6.4 `src/core_atmosphere/CMakeLists.txt` (if using CMake) or `Makefile`
Add chemistry source files to the build.

### 6.5 `src/core_init_atmosphere/mpas_init_atm_cases.F`
Initialize O3 with a realistic vertical profile (e.g., stratospheric peak at ~25 km,
tropospheric background ~30 ppbv). Initialize O and O1D to zero (they're produced
photochemically). This gives MICM something meaningful to work with.

### 6.6 `docker/Containerfile`
Ensure MUSICA is built with TUV-x support (`-DMUSICA_ENABLE_TUVX=ON`). Verify the
existing container config includes this.

### 6.7 `scripts/run_jw_test.sh`
Add Chapman species to output stream. Set chemistry namelist options.

---

## 7. Implementation Steps

### Step 2.1: Verify MUSICA Build Flags in Container
**Goal**: Confirm the container builds MUSICA with both MICM and TUV-x enabled.
- Check Containerfile for `-DMUSICA_ENABLE_TUVX=ON` and `-DMUSICA_ENABLE_MICM=ON`
- Verify `pkg-config --libs musica-fortran` includes TUV-x
- If needed, update Containerfile and rebuild deps image

### Step 2.2: Add Chapman Species to Registry.xml
**Goal**: Register O3, O, O1D as advected scalars.
- Add species to `scalars` var_array (array_group="chapman")
- Add matching entries to `scalars_tend` and `lbc_scalars`
- Add to `Registry.xml` in both `core_atmosphere` and `core_init_atmosphere`
- Add chemistry namelist options

### Step 2.3: Initialize Chapman Species in JW Case
**Goal**: Provide realistic initial O3 profile.
- In `mpas_init_atm_cases.F`, JW case section:
  - O3: approximate stratospheric profile (e.g., from US Standard Atmosphere or
    a simple analytic formula: mixing ratio peaks ~8 ppmv at 35 km, ~30 ppbv at surface)
  - O, O1D: initialize to zero (photochemically produced)

### Step 2.4: Create Chemistry Utility Module
**Goal**: Solar geometry calculations.
- `mpas_chemistry_utils.F90`:
  - `compute_solar_zenith_angle()` using MPAS cell lat/lon and simulation time
  - `compute_earth_sun_distance()` from Julian day

### Step 2.5: Create MICM Interface Module
**Goal**: Wrap MICM Fortran API for MPAS use with batched solving.
- `mpas_chemistry_micm.F90`:
  - `micm_setup()`: Create `micm_t` from config path, compute `n_grid_cells = nCellsSolve * nVertLevels`, allocate `state_1` and (if needed) `state_2` using the CAM-SIMA two-state residual pattern (`max_cells = micm%get_maximum_number_of_grid_cells()`), extract `species_ordering` and `rate_parameters_ordering`, build index mapping arrays
  - `micm_solve_batch(dt)`: Loop over states — for each batch, call `update_micm_state()` with offset, `micm%solve()`, `state%update_references()`, `extract_from_micm_state()` with offset
  - `micm_cleanup()`: Deallocate states and solver

### Step 2.6: Create State Copy Module
**Goal**: Efficient MPAS ↔ MICM data transfer with unit conversion.
- `mpas_chemistry_state.F90`:
  - `update_micm_from_mpas()`: Copy T, P, ρ to conditions; convert species mmr→mol/m³
    using strided access; copy rate_parameters
  - `update_mpas_from_micm()`: Convert species mol/m³→mmr; write back to scalars array
  - Molar masses: O3=0.048, O=0.016, O1D=0.016 kg/mol
  - Air density: `rho_dry = zz(k,iCell) * rho_zz(k,iCell)`

### Step 2.7: Create TUV-x Interface Module
**Goal**: Wrap TUV-x Fortran API for MPAS use.
- `mpas_chemistry_tuvx.F90`:
  - `tuvx_setup()`: Create height grid (nVertLevels+1 sections), wavelength grid
    (from TUV-x config), temperature/albedo/ET flux profiles, dry air/O2/O3 profiles,
    cloud/aerosol radiators. Construct `tuvx_t`. Get photolysis rate constant ordering.
    Build mapping from TUV-x photo rates → MICM rate parameters.
  - `tuvx_run_column()`: For one column:
    1. Compute height midpoints/interfaces in km from `zgrid`
    2. Set height grid values (no inversion needed — MPAS is already surface-to-top)
    3. Set temperature profile from MPAS temperature
    4. Set gas species profiles (dry air, O2, O3) — convert mmr to mol/cm³
    5. Set surface albedo (constant for JW case)
    6. Call `tuvx%run(sza, earth_sun_dist, photo_rates, heating_rates)`
    7. Map photo rates to MICM rate parameter ordering
    8. Return rate_parameters(nVertLevels, n_rate_params)
  - `tuvx_cleanup()`: Deallocate

### Step 2.8: Create Chemistry Driver Module
**Goal**: Orchestrate TUV-x → MICM workflow.
- `mpas_chemistry_driver.F90`:
  - `chemistry_init(domain)`:
    1. Read namelist options
    2. Call `micm_setup(config_micm_config_path, nVertLevels)`
    3. Call `tuvx_setup(config_tuvx_config_path, nVertLevels)`
    4. Build species index mapping: MPAS scalar index ↔ MICM species index
       (match by name: state%species_ordering maps MICM names → MICM indices;
        MPAS `index_o3`, `index_o`, `index_o1d` from pool dimension queries)
  - `chemistry_timestep(block, mesh, state, diag, configs, dt)`:
    1. Get pointers: scalars, zgrid, zz, rho_zz, theta_m, exner, pressure
    2. Loop over cells (1 to nCellsSolve):
       a. Compute T(k), P(k), rho_dry(k) for all levels
       b. Compute SZA and earth-sun distance
       c. Call `tuvx_run_column()` → photolysis rates
       d. Call `update_micm_from_mpas()` — load state with T,P,ρ,species,rates
       e. Call `micm%solve(dt, state, ...)` or `micm_solve_column()`
       f. Call `state%update_references()` (C++ may swap pointers)
       g. Call `update_mpas_from_micm()` — write updated species back to scalars
    3. End cell loop
  - `chemistry_finalize()`: Call micm_cleanup(), tuvx_cleanup()

### Step 2.9: Wire Chemistry into MPAS
**Goal**: Connect chemistry to the MPAS physics/dynamics cycle.
- In `mpas_atm_core.F`:
  - `atm_core_init()`: Call `chemistry_init(domain)` after physics init
  - `atm_core_finalize()`: Call `chemistry_finalize()`
- In `mpas_atmphys_driver.F`:
  - Add call to `chemistry_timestep()` after radiation, gated by
    `config_chemistry_enabled`
  - **Alternative** (simpler for JW with `physics_suite='none'`): Call
    `chemistry_timestep()` directly from `atm_do_timestep()` before `atm_timestep()`,
    bypassing the physics framework entirely. This avoids modifying the physics driver
    and works even when no physics suite is configured.

### Step 2.10: Update Build System
**Goal**: Compile chemistry modules.
- Add `chemistry/` subdirectory to the atmosphere core Makefile
- Add source files with correct dependencies
- Ensure MUSICA Fortran includes and libraries are found via `pkg-config`

### Step 2.11: Prepare Chapman + TUV-x Config Files
**Goal**: Package mechanism and photolysis configs.
- Verify `configs/v0/chapman/` has correct species.json and reactions.json
- Create `configs/tuvx/chapman.json` — a minimal TUV-x config with only 3
  photolysis reactions, extracted from `configs/tuvx/ts1_tsmlt.json`:
  - `jO2`: combine `jo2_a` + `jo2_b` cross-section and quantum yield entries
    into a single reaction. Both use `data/cross_sections/O2_1.nc` with
    `apply O2 bands: true`; merge the quantum yield override bands so the
    single reaction produces the total O2 photolysis rate.
  - `jO3->O1D`: copy `jo3_a` config (O3 cross sections from `O3_1-4.nc`,
    quantum yield type `O3+hv->O2+O(1D)`)
  - `jO3->O`: copy `jo3_b` config (same cross sections,
    quantum yield type `O3+hv->O2+O(3P)`)
  - Include wavelength grid, height grid placeholder, standard profiles
    (temperature, surface albedo, ET flux, air, O2, O3) — same structure as
    TS1/TSMLT but with only 3 reactions in the photolysis block
  - All data files come from existing `configs/tuvx/data/`
- Create 1:1 photolysis-to-MICM mapping config:
  `jO2` → `PHOTO.jO2`, `jO3->O1D` → `PHOTO.jO3->O1D`, `jO3->O` → `PHOTO.jO3->O`
- Verify the mapping by dumping both orderings at init time and logging the
  resulting index pairs

### Step 2.12: Update JW Test Script
**Goal**: Enable chemistry in test runs.
- Set `config_chemistry_enabled = .true.`
- Set config paths for MICM, TUV-x, mapping
- Add O3, O, O1D to output stream
- Set `config_chemistry_dt` (= dynamics dt for simplicity, or a multiple)

### Step 2.13: Build and Debug
**Goal**: Get a clean compilation and first successful run.
- Rebuild container with updated source
- Fix compilation errors iteratively
- Run JW + Chapman on 480-km mesh
- Debug runtime errors (MICM solver convergence, TUV-x grid dimension mismatches, etc.)

### Step 2.14: Create Phase 2 Verification Notebook
**Goal**: Validate Chapman chemistry results.
- `verification/phase02_chapman.ipynb`:

**Section 1: O3 Diurnal Cycle**
- Map of O3 column at local noon vs. local midnight
- Timeseries of globally-averaged O3 showing diurnal oscillation
- Verify: O3 increases on dayside (jO2 → O + O → O3), decreases at night
  (only loss reactions active)

**Section 2: Odd Oxygen Conservation**
- Compute Ox = O + O1D + O3 (in moles) at each timestep
- Verify: total Ox is conserved within transport + chemistry splitting error
  (Chapman has no external sources/sinks of Ox)
- Tolerance: < 0.1% change over test period

**Section 3: Photolysis Rate Profiles**
- If we add diagnostic output of photolysis rates:
  - Vertical profile of jO2, jO3 at noon (should decrease with optical depth)
  - jO3→O1D should be significant only in UV, decrease rapidly below ~20 km

**Section 4: Single-Column Comparison**
- Extract one column (e.g., equatorial, subsolar point) from MPAS output
- Run the same column through standalone MUSICA Python solver (tutorial 10)
- Compare O3 timeseries — should match within solver tolerance
- This validates the coupling (unit conversions, grid mapping, rate parameter passing)

**Section 5: Species Non-negativity**
- Verify O3, O, O1D ≥ 0 everywhere at all times
- MICM solver should maintain non-negativity, but verify

**Section 6: Automated Pass/Fail Checks**
1. O3 diurnal amplitude > 0 at subsolar point
2. Ox conservation: max |ΔOx/Ox₀| < 0.1%
3. All species non-negative
4. O3 stratospheric peak exists (max O3 above 20 km > 10× surface value)
5. Nightside photolysis rates = 0 (SZA > 110°)
6. Single-column regression within 1% of Python reference

---

## 8. Technical Notes

### 8.1 Grid Alignment (from project lead)
Both MPAS and TUV-x vertical grids start at the surface and move up. This is unlike
CAM-SIMA, which uses a top-to-surface grid and must invert before passing to TUV-x.
The CAM-SIMA height grid code (`musica_ccpp_tuvx_height_grid.F90`) reverses array
elements — we skip that reversal but keep the same TUV-x grid structure (extra half-
layer at bottom, exo-layer at top).

### 8.2 TUV-x Interfaces ↔ MPAS Midpoints (from project lead)
TUV-x calculates photolysis rate constants at grid **interfaces**. We set up the
TUV-x grid so that TUV-x interfaces are at MPAS vertical grid **midpoints** (with a
half-size layer near the surface and an exo-layer at the top). This ensures the
computed photolysis rates are representative of the conditions at MPAS layer centers.

### 8.3 Efficient Data Copy (from project lead)
Copying data into and out of MICM/TUV-x needs to be very efficient. MICM's internal
matrices have a vector-grouped structure, and the state arrays are flat with
configurable strides. The reference copy pattern is in `musica_ccpp_micm_util.F90`:
- Uses pointer remapping to collapse 2D host arrays → 1D
- Accesses MICM state arrays with `species_strides%grid_cell` and
  `species_strides%variable` for stride-aware indexing
- Converts between mass mixing ratio (host) and molar concentration (MICM)

For Phase 2 (single column, ~26 levels, 3 species), raw performance is not critical.
But the copy routines should be designed with the Phase 6 multi-column batched solve
in mind — avoid unnecessary copies, use stride-aware access from the start.

### 8.4 MICM Two-State Pattern
MICM pre-allocates states for a fixed "maximum" number of grid cells (set at build
time for SIMD vectorization). When the total number of grid cells exceeds this,
`musica_ccpp_micm.F90` uses two states: `state_1` for full batches and `state_2` for
the residual. For Phase 2 (one column = nVertLevels cells), a single state suffices.

### 8.5 MPAS Single Precision
MPAS-A is built with `SINGLE_PRECISION`. The MUSICA Fortran API uses `real64`
(double precision). We need explicit precision conversion at the MPAS ↔ MUSICA
boundary. This is a potential source of tolerance issues in verification.

### 8.6 O3 Initial Condition Profile
A simple analytic O3 profile for the JW case:
```fortran
! US Standard Atmosphere O3 approximation (mmr in kg/kg)
z_km = zmid_km(k)
if (z_km < 20.0) then
  o3_vmr = 30.0e-9 + (z_km / 20.0) * 70.0e-9  ! 30-100 ppbv linear increase
else if (z_km < 35.0) then
  o3_vmr = 100.0e-9 + ((z_km - 20.0) / 15.0) * 7900.0e-9  ! ramp to 8 ppmv peak
else
  o3_vmr = 8000.0e-9 * exp(-(z_km - 35.0) / 8.0)  ! exponential decay above
end if
o3_mmr = o3_vmr * 0.048 / 0.029  ! vmr → mmr (M_O3/M_air)
```

### 8.7 Config Format: v0 vs v1
The `configs/v0/chapman/` directory uses the camp-data format (species.json +
reactions.json, referenced from config.json). The `configs/v1/chapman/` has both
config.json and config.yaml. Either format works with the current MICM Fortran
API — use v1 if available as it's the newer format.

---

## 9. Risk Register

| Risk | Mitigation |
|------|-----------|
| MUSICA container not built with TUV-x | Check build flags in Step 2.1; rebuild if needed |
| Combining O2 photolysis channels in TUV-x config | Verify combined jO2 against sum of jo2_a + jo2_b from full TS1 config |
| Precision mismatch (MPAS float32 ↔ MUSICA float64) | Explicit conversion; wider tolerances in verification |
| MICM solver divergence with bad initial conditions | Start with reasonable O3 profile; check for NaN/Inf after solve |
| SZA calculation wrong without active radiation | Implement standalone SZA computation in Step 2.4 |
| Photolysis rate mapping names don't match | Dump MICM rate_parameters_ordering and TUV-x photolysis ordering; verify names match |
| Large Makefile changes for MPAS atmosphere core | Keep chemistry in separate subdirectory; minimal Makefile edits |
| JW test has no cloud water for TUV-x cloud optics | Set cloud optics to zero (clear sky); fine for Chapman |

---

## 10. Dependencies and Prerequisites

- ✅ Phase 0: Container infrastructure (working)
- ✅ Phase 1: Passive tracer infrastructure (working)
- ✅ MUSICA libraries built in container (MICM + TUV-x)
- ✅ Chapman mechanism configs exist in `configs/v0/chapman/` and `configs/v1/chapman/`
- ✅ TUV-x cross-section data available in `configs/tuvx/data/`
- ✅ TUV-x reference: `configs/tuvx/ts1_tsmlt.json` has reaction configs to extract from
- ⬜ Chapman TUV-x config `configs/tuvx/chapman.json` (needs creation from TS1/TSMLT)
- ⬜ Photolysis-to-MICM rate parameter mapping config (needs creation)

---

## 11. Estimated Complexity

| Component | New Lines | Difficulty | Notes |
|-----------|-----------|-----------|-------|
| Registry.xml changes | ~30 | Low | Follow Phase 1 pattern |
| JW O3 initialization | ~30 | Low | Analytic profile |
| mpas_chemistry_utils.F90 | ~80 | Medium | SZA/earth-sun distance |
| mpas_chemistry_micm.F90 | ~150 | Medium | MICM wrapper |
| mpas_chemistry_state.F90 | ~120 | Medium | Copy routines |
| mpas_chemistry_tuvx.F90 | ~250 | High | Most complex: grid setup, profile management |
| mpas_chemistry_driver.F90 | ~200 | Medium | Orchestration |
| Build system changes | ~20 | Low | Makefile additions |
| Config files | ~50 | Medium | TUV-x config, mapping |
| Test script updates | ~15 | Low | Namelist + output |
| Verification notebook | ~300 | Medium | 6 sections, automated checks |
| **Total** | **~1250** | | |
