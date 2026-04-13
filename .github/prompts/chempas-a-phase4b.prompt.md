# Phase 4b: Mid-Development Review — Findings & Plan

## Summary

This document captures all issues found during the Phase 4b code review and proposes
remediation plans. Three areas were reviewed: (1) hardcoded species references,
(2) MPAS/TUV-x/MICM grid cell mapping correctness, and (3) debug code, bugs, and cleanup.

**Goal:** After this phase, adding a new mechanism requires **zero Fortran code changes**.
The workflow is: (1) add a `chemistry_data/<mechanism>/` folder with MICM config +
optional TUV-x config + optional IC profile file, (2) run a script to generate
Registry.xml entries, (3) rebuild and run.

**Invariant:** No species names appear anywhere in Fortran code — not as string
literals, not as variable names, not as named constants. All species identity is
resolved at runtime from configuration files.

---

## 1. Hardcoded Species References

### 1a. `mpas_chemistry_species.F90` — O3 name check + TUV-x profile discovery

**Location:** `chem_species_init()`, line ~160

```fortran
if (trim(species_name) == 'O3') then
   tuvx_o3_mpas_idx = mpas_idx
end if
```

**Issue:** Hard-codes the string `'O3'` to identify which species to pass to TUV-x.
More broadly, the entire TUV-x integration assumes specific species (O2, O3) by name.

**Fix:** Replace with a fully generic TUV-x profile discovery system. Each MICM
species that maps to a TUV-x profile gets these properties:

```json
{
  "name": "O3",
  "molecular weight [kg mol-1]": 0.0479982,
  "__is_tuvx_profile": "O3",
  "__tuvx_scale_height [km]": 4.5,
  "__tuvx_column_density_method": "arithmetic"
}
```

During `chem_species_init()`, build a dynamic array of TUV-x profile descriptors:

```fortran
type :: tuvx_profile_info
  character(len=64) :: profile_name     ! TUV-x profile name (e.g., "O3")
  integer           :: mpas_idx         ! MPAS scalar index (-1 if constant)
  integer           :: micm_idx         ! MICM species index
  real(real64)      :: molecular_weight  ! from MICM property
  real(real64)      :: scale_height      ! from __tuvx_scale_height [km]
  real(real64)      :: constant_vmr      ! from __mpas_constant_vmr (0 if advected)
  logical           :: use_arithmetic_mean ! from __tuvx_column_density_method
end type

type(tuvx_profile_info), allocatable :: tuvx_profiles(:)
```

The species module iterates all MICM species, and for any with `__is_tuvx_profile`,
populates a `tuvx_profile_info` entry. No species name strings appear in Fortran.
The variable `tuvx_o3_mpas_idx` is eliminated — replaced by the generic array.

The TUV-x module receives this array and uses it to create/set profiles generically.
See 1b for the TUV-x side.

---

### 1b. `mpas_chemistry_tuvx.F90` — Fully generic TUV-x gas profiles

**Location:** Throughout `tuvx_setup()` and `tuvx_run_column()`

Currently hardcoded:
- `MW_O3 = 0.048`, `MW_O2 = 0.032`, `VMR_O2 = 0.2095`
- `SCALE_HEIGHT_O3 = 4.5`, `SCALE_HEIGHT_O2 = 8.01`
- String literals `"O2"`, `"O3"` for profile creation/retrieval
- Separate `o2_profile`, `o3_profile` module pointers
- Separate calls to `set_gas_profile` for O2 and O3 with different parameters

**Fix:** Replace ALL of the above with a generic loop over `tuvx_profiles(:)` (the
array built by species discovery in 1a). Every species-specific value comes from
the `tuvx_profile_info` struct — molecular weight, scale height, constant VMR,
column density method, and the TUV-x profile name string.

**`tuvx_setup()`:**
- Replace separate `o2_profile`, `o3_profile` pointers with an array of
  `profile_t` pointers, one per `tuvx_profiles` entry.
- Create each gas profile using `profile_name` from the config, not a string literal.
- After constructing the TUV-x solver, re-obtain handles via the same names.

**`tuvx_run_column()`:**
- Replace the two separate `set_gas_profile` calls (one for O2, one for O3) with a
  single loop over `tuvx_profiles`:

```fortran
do i = 1, size(tuvx_profiles)
  if (tuvx_profiles(i)%constant_vmr > 0) then
    ! Constant species: derive MMR = constant_vmr * species_mw / mw_air
    mmr_col(:) = tuvx_profiles(i)%constant_vmr &
               * tuvx_profiles(i)%molecular_weight / MW_AIR
  else
    ! Advected species: extract from MPAS scalars
    mmr_col(:) = scalars(tuvx_profiles(i)%mpas_idx, 1:nVertLevels, iCell)
  end if
  call set_gas_profile(gas_profiles(i), nVertLevels, rho_col, mmr_col, &
                       tuvx_profiles(i)%molecular_weight, &
                       height_deltas, tuvx_profiles(i)%scale_height, &
                       tuvx_profiles(i)%use_arithmetic_mean, errmsg, errcode)
end do
```

The named constants `MW_O3`, `MW_O2`, `VMR_O2`, `SCALE_HEIGHT_O3`, `SCALE_HEIGHT_O2`
are eliminated entirely — these values now come from the MICM species properties at
runtime. Only `MW_AIR`, `AVOGADRO`, and `SCALE_HEIGHT_AIR` remain as universal
physical constants (they describe air, not any individual species).

**`tuvx_run_column()` signature change:**

The current signature accepts `o3_mmr(:)` as a dedicated parameter — an O3-specific
interface. Replace with a generic interface that receives the full MPAS scalars column:

```fortran
subroutine tuvx_run_column(nVertLevels, zgrid, temperature, rho_dry, &
                            scalars_col, sza, earth_sun_dist,        &
                            photo_rates, errmsg, errcode)
  real (kind=RKIND), intent(in) :: scalars_col(:,:)  ! (nScalars, nVertLevels)
```

Inside `tuvx_run_column`, the generic loop over `tuvx_profiles` extracts MMR data
from `scalars_col(profile%mpas_idx, :)` for advected species, or computes from
`constant_vmr` for constant species. The caller (`run_tuvx_photolysis`) passes
`scalars(:, :, iCell)` without knowing which species are TUV-x profiles.

The TUV-x config already defines which profiles exist (the radiators reference
`"O2"`, `"O3"` etc.). The MICM config's `__is_tuvx_profile` property provides the
binding between MICM species and TUV-x profiles. The Fortran code never mentions
any species by name.

---

### 1c. `mpas_chemistry_tuvx.F90` — `is_o3` flag in `set_gas_profile()`

**Location:** `set_gas_profile()`, line ~474

```fortran
logical, intent(in) :: is_o3
```

Used to select arithmetic mean (O3) vs geometric mean (well-mixed) for column density
calculation.

**Fix:** Rename parameter to `use_arithmetic_mean` (logical). The value comes from the
`tuvx_profile_info` struct, which reads it from the MICM species property
`"__tuvx_column_density_method": "arithmetic"` or `"geometric"` (default:
`"geometric"`). No species name needed to select the method.

---

### 1d. Registry.xml — Per-mechanism species packages (AUTOMATE)

**Location:** Lines 453-455, 1724-1756

```xml
<package name="chem_chapman_in" description="..."/>
<var name="o3" ... packages="chem_chapman_in;chem_emis_dep_in"/>
```

Each mechanism requires hand-written package and `<var>` entries.

**Fix:** Write `scripts/generate_registry.py` — a Python script that:
1. Reads a MICM mechanism config (v0 or v1 format)
2. Identifies advected species (those with `"__is_advected": true` property)
3. Generates the Registry.xml `<package>`, `<var>`, and `<var_struct>` snippets
4. Generates the `mpas_atm_core_interface.F` package activation block

This makes adding a new mechanism a one-command operation. Details in Section 4.

---

### 1e. `mpas_atm_core_interface.F` — Mechanism name switch (ELIMINATE)

**Location:** Lines 232-240

```fortran
if (trim(config_chemistry_mechanism) == 'chapman') then
   chem_chapman_inActive = .true.
...
```

**Fix:** Replace `config_chemistry_mechanism` with `config_chemistry_config_path` — a
single namelist option pointing to the top-level mechanism folder (e.g.,
`"chemistry_data/chapman_emis_dep"`). The chemistry driver derives all sub-paths from
this one root:
- MICM config: `<root>/config.json` or `<root>/micm/config.json` (auto-detect)
- TUV-x config: `<root>/tuvx/config.json` (if present → TUV-x enabled)
- TUV-x mapping: `<root>/tuvx_micm_mapping.json` (if present)

The package activation in `mpas_atm_core_interface.F` uses a naming convention: the
mechanism folder basename maps to `chem_<basename>_inActive`. The code iterates all
registered chemistry packages and activates the one matching the basename. No
if/else-if chain; no code change to add a mechanism.

This also eliminates three separate namelist options (`config_micm_config_path`,
`config_tuvx_config_path`, `config_tuvx_micm_mapping_path`) in favor of one path.

---

### 1f. `mpas_init_atm_cases.F` — IC initialization (CONFIG FILE)

**Location:** Lines 1295-1330

Hard-codes `index_o3`, `index_o`, `index_o1d` for initial conditions with an analytic
O3 profile formula baked into Fortran.

**Fix:** Replace all hardcoded IC code with a generic config-file-driven approach.
Each mechanism provides an optional initial conditions file:

```
chemistry_data/<mechanism>/initial_conditions.csv
```

Format: a simple CSV with one column per species and one row per MPAS vertical level
(k=1=top to k=nVertLevels=bottom). Column headers are species names matching the
MICM config. Values are mass mixing ratios [kg kg⁻¹]. Species not listed start at 0
(appropriate for short-lived intermediates like O, O1D).

Example `initial_conditions.csv` for chapman_emis_dep:
```csv
O3
0.0
0.0
...
2.85e-8
```

(For the JW test, the O3 profile is the US Standard Atmosphere discretized onto the
26 MPAS levels. This is computed once by a helper script and stored in the CSV.)

The `mpas_init_atm_cases.F` code becomes fully generic:
```fortran
! Read IC profile file (path derived from config_chemistry_config_path)
! For each column header, look up index_<name> from Registry
! Apply the profile to every column
! Species with no entry in the file default to 0
```

No species names appear in the Fortran code. The `generate_registry.py` script does
NOT generate IC code — the init_atmosphere core reads the CSV at runtime.

A helper script `scripts/generate_initial_conditions.py` can create the CSV from
analytic profiles (US Standard O3, etc.) for convenience.

---

### 1g. `DEFAULT_SURFACE_ALBEDO = 0.1` in `mpas_chemistry_tuvx.F90`

**Fix:** Add `config_chemistry_surface_albedo` namelist parameter. Low priority but
included in this phase for completeness.

---

## 2. MPAS / TUV-x / MICM Grid Cell Mapping Verification

### Current Mapping Architecture

Three distinct index spaces are used:

1. **MPAS**: `scalars(index_species, k, iCell)` — species-first, then vertical level
   (k=1=top, k=nVertLevels=bottom/surface), then horizontal cell.

2. **MICM**: Flat 1D arrays with stride-based indexing.
   `flat_idx = (i_cell - 1) * gc_stride + (var_idx - 1) * var_stride + 1`
   where `i_cell` counts sequentially: outer loop iCell=1..nCellsSolve,
   inner loop k=1..nVertLevels, so `i_cell = (iCell-1)*nVertLevels + k`.

3. **TUV-x**: Per-column, with nVertLevels+1 layers (nVertLevels+2 interfaces).
   - MPAS interface heights → TUV-x interfaces (with half-layer at surface +
     exo-layer at top)
   - MPAS midpoints → TUV-x midpoints
   - TUV-x layer k (1-based) maps to MPAS level k

### Specific Mapping Concerns

- **MPAS→MICM cell ordering**: The flattening `i_cell = (iCell-1)*nVertLevels + k`
  must be consistent across `update_micm_from_mpas`, `update_mpas_from_micm`,
  `emissions_set_rates`, and `deposition_set_rates`. All use the same nested loop
  pattern, so this should be consistent — but worth verifying explicitly.

- **MPAS→TUV-x height mapping**: MPAS `zgrid` is on interfaces (k=1..nVertLevels+1).
  TUV-x has nVertLevels+2 interfaces building in a half-layer at surface and
  exo-layer at top. The mapping in `tuvx_run_column` splits each MPAS layer in two
  (midpoint → interface), which is subtle and error-prone.

- **TUV-x→MICM photolysis mapping**: After TUV-x produces rates on its layers,
  `photo_rates(k, :)` for k=1..nVertLevels is copied from TUV-x layer k → MICM grid
  cell. The layer→level correspondence must align with the height grid mapping.

- **Emission/Deposition surface level**: Both modules use `k == nVertLevels` as the
  surface. This must match the MPAS convention that `zgrid(nVertLevels+1)` is the
  surface interface.

### Verification Notebook: `verification/phase4b_grid_mapping.ipynb`

**Goal:** Quantitatively verify that the MPAS→MICM→TUV-x grid cell mapping is correct
by running a single-timestep simulation and examining the internal state.

**Approach:**
1. **Add diagnostic output to a single-timestep run**: Add temporary write statements
   or use a special diagnostic mode to dump:
   - MPAS `zgrid(1:nVertLevels+1, iCell)` for a representative cell
   - TUV-x `tuvx_interfaces` and `tuvx_midpoints` for that column
   - MICM `concentrations` and `rate_parameters` flat arrays for that column
   - MICM strides (gc_stride, var_stride)
   - MPAS→MICM index mappings (advected_mpas_idx, advected_micm_idx)
   - Photolysis rates: TUV-x `photo_out` vs MICM rate_parameters

2. **Notebook cells**:
   - **Cell 1**: Run a 1-timestep chapman simulation, collect diagnostic output
   - **Cell 2**: Parse and display the MPAS vertical grid for one cell: interface
     heights, midpoint heights, layer thicknesses. Verify monotonically increasing
     from surface to top.
   - **Cell 3**: Reconstruct the TUV-x grid from MPAS data using the same algorithm
     as `tuvx_run_column`, display interfaces/midpoints/deltas. Verify
     nVertLevels+2 interfaces, nVertLevels+1 layers, all positive deltas.
   - **Cell 4**: Verify MICM flat index mapping. For each species and each vertical
     level, compute the expected flat index and verify it matches what the code
     produces. Check roundtrip: MPAS scalar → MICM concentration → MPAS scalar.
   - **Cell 5**: Verify emission/deposition rates. Check that emission rate is
     non-zero only at k=nVertLevels (surface). Check that the volumetric rate
     formula `flux * 1e4 / (N_A * dz)` gives the expected value.
   - **Cell 6**: Verify TUV-x→MICM photolysis rate mapping. For a daytime column,
     check that TUV-x layer k rates match the MICM rate parameter for that grid
     cell. For a nighttime column (SZA > 110°), verify all rates are zero.
   - **Cell 7**: Summary — pass/fail assertions for each mapping component.

**Implementation approach:** Rather than instrumenting the Fortran code, the notebook
can reconstruct the mappings from Python using the same algorithms and known MICM
strides, then compare against the output file. This avoids modifying production code.

---

## 3. Debug Code, Bugs, and Cleanup

### 3a. **BUG: Temperature not passed to TUV-x** (CRITICAL)

**Location:** `mpas_chemistry_driver.F90`, `run_tuvx_photolysis()`, lines 423-425

```fortran
call tuvx_run_column(nVertLevels, zgrid_col, rho_col, &
                      rho_col, o3_col, sza_r64, earth_sun_r64, &
                      photo_col, errmsg, errcode)
```

`rho_col` (dry air density) is passed for **both** the `temperature` and `rho_dry`
arguments. The `temp_2d` array is computed in `chemistry_timestep` but never passed
to `run_tuvx_photolysis`.

**Impact:** TUV-x receives density values (~1.2 kg/m³ at surface, ~0.03 at top)
instead of temperature values (~290 K at surface, ~220 at tropopause). This affects
the temperature profile used in TUV-x's cross-section calculations and Rayleigh
scattering. Photolysis rates will be incorrect, though the Chapman mechanism's simple
cross-sections may partially mask the error.

**Fix:**
1. Add `temp_2d` as a parameter to `run_tuvx_photolysis`
2. Extract `temp_col(1:nVertLevels) = temp_2d(1:nVertLevels, iCell)` per column
3. Pass `temp_col` as the temperature argument to `tuvx_run_column`

---

### 3b. All hardcoded array sizes → allocatable / automatic

All fixed-size arrays must be replaced with dynamically-sized arrays. No array limits
anywhere. Specifically:

**`mpas_chemistry_driver.F90` — `run_tuvx_photolysis()`:**
- `o3_col(200)` → `o3_col(nVertLevels)` (automatic)
- `zgrid_col(201)` → `zgrid_col(nVertLevels+1)` (automatic)
- `rho_col(200)` → `rho_col(nVertLevels)` (automatic)
- `photo_col(200, 10)` → `photo_col(nVertLevels, n_photo_rxns_local)` (automatic)

**`mpas_chemistry_driver.F90` — `chemistry_timestep()`:**
- Remove dead `o3_col(200)`, `zgrid_col(201)`, `rho_col(200)`, `photo_col(200,10)`
  declarations (these are only used in `run_tuvx_photolysis`).

**`mpas_chemistry_species.F90` — `chem_species_init()`:**
- `tmp_adv_mpas(100)` etc → `allocatable`, sized to `n_species` from MICM
- `tmp_con_micm(100)` etc → same
- Remove error check "Too many advected species (max 100)"

**`mpas_chemistry_emissions.F90` — `emissions_init()`:**
- `tmp_rp_idx(100)`, `tmp_flux(100)` → `allocatable`, sized to `n_species`
- Remove error check "Too many emission species (max 100)"

**`mpas_chemistry_deposition.F90` — `deposition_init()`:**
- `tmp_rp_idx(100)`, `tmp_vel(100)` → `allocatable`, sized to `n_species`
- Remove error check "Too many deposited species (max 100)"

---

### 3c. `n_deposited` unnecessarily public

**Location:** `mpas_chemistry_deposition.F90`, line 24

**Fix:** Remove from `public` statement.

---

### 3d. Debug code check — CLEAN

No `write(*,*)`, `print *`, `TODO`, `FIXME`, `HACK`, `DEBUG`, or `TEMP` markers found
in any chemistry module. All logging uses the proper `mpas_log_write` facility.

---

### 3e. Unused variable in `chemistry_timestep`

**Location:** `mpas_chemistry_driver.F90`, lines 241-242

**Fix:** Remove the unused `o3_col`, `zgrid_col`, `rho_col` declarations from
`chemistry_timestep` (covered by 3b).

---

### 3f. Duplicate Avogadro constant definitions

Both `mpas_chemistry_tuvx.F90` and `mpas_chemistry_emissions.F90` define their own
`AVOGADRO` parameter with the same value. `mpas_chemistry_species.F90` defines
`MW_AIR`.

**Fix:** Consolidate universal physical constants into `mpas_chemistry_utils.F90`:
`AVOGADRO`, `MW_AIR`, `SCALE_HEIGHT_AIR`. Import from there everywhere.

Species-specific values (`MW_O3`, `MW_O2`, `VMR_O2`, `SCALE_HEIGHT_O3`,
`SCALE_HEIGHT_O2`) are **eliminated as named constants** — they now come from MICM
species properties at runtime via the `tuvx_profiles` array (see 1a/1b).

---

## 4. Registry Generation Script: `scripts/generate_registry.py`

### Purpose

Eliminate all manual Registry.xml / Fortran edits when adding a mechanism. The script
reads a MICM mechanism config and produces:
1. Registry.xml snippet (package + var entries)
2. `mpas_atm_core_interface.F` package activation snippet

### Species Classification via `__is_advected`

Add a new species property to MICM configs:

```json
{
  "name": "O3",
  "molecular weight [kg mol-1]": 0.0479982,
  "__is_advected": true,
  "__is_tuvx_profile": "O3",
  "__tuvx_scale_height [km]": 4.5,
  "__tuvx_column_density_method": "arithmetic"
}
```

Classification rules:
- `"__is_advected": true` → MPAS tracer scalar, gets a Registry `<var>` entry,
  gets an `index_<name>` dimension, is transported by MPAS dynamics.
- `"__mpas_constant_vmr": <value>` → Constant species diagnosed from air composition.
  No Registry entry. No MPAS scalar. Set in MICM from VMR × air density each step.
- `"is third body": true` → Third body (M). No Registry entry. Handled by MICM.
- Otherwise (no `__is_advected`, no `__mpas_constant_vmr`, no `is third body`) →
  the script warns and skips.

### Full Set of MICM Species Properties for Runtime Configuration

| Property | Used by | Example |
|---|---|---|
| `"molecular weight [kg mol-1]"` | MICM (existing), TUV-x profile MW | `0.0479982` |
| `"__is_advected"` | Registry script, species module | `true` |
| `"__mpas_constant_vmr"` | Species module (existing) | `0.2095` |
| `"__is_tuvx_profile"` | Species module → TUV-x | `"O3"` |
| `"__tuvx_scale_height [km]"` | TUV-x exo-layer | `4.5` |
| `"__tuvx_column_density_method"` | TUV-x `set_gas_profile` | `"arithmetic"` |
| `"__mpas_surface_emission_flux [molec cm-2 s-1]"` | Emissions (existing) | `1.0e10` |
| `"__mpas_surface_deposition_velocity [cm s-1]"` | Deposition (existing) | `0.5` |

No species names appear in Fortran code. All species identity resolution happens
via these properties at runtime.

### Script Interface

```
python scripts/generate_registry.py <mechanism_config_path> [--name <mechanism_name>]
```

- `<mechanism_config_path>`: Path to top-level mechanism folder (e.g.,
  `chemistry_data/chapman_emis_dep`). Auto-detects v0 vs v1 format.
- `--name`: Override mechanism name (default: folder basename).

### Output Files

The script writes two files to `<mechanism_config_path>/generated/`:

**`registry_snippet.xml`:**
```xml
<!-- Auto-generated for mechanism: chapman_emis_dep -->
<package name="chem_chapman_emis_dep_in"
         description="Chemistry species for chapman_emis_dep mechanism"/>

<!-- Add these inside the scalars var_array in Registry.xml -->
<var name="o3" array_group="chem_chapman_emis_dep" units="kg kg^{-1}"
     description="O3 mass mixing ratio"
     packages="chem_chapman_emis_dep_in"/>
<var name="o" array_group="chem_chapman_emis_dep" units="kg kg^{-1}"
     description="O mass mixing ratio"
     packages="chem_chapman_emis_dep_in"/>
...
```

**`core_interface_snippet.F`:**
```fortran
! Auto-generated package activation for: chapman_emis_dep
nullify(chem_chapman_emis_dep_inActive)
call mpas_pool_get_package(packages, 'chem_chapman_emis_dep_inActive', &
                           chem_chapman_emis_dep_inActive)
if (associated(chem_chapman_emis_dep_inActive)) then
   chem_chapman_emis_dep_inActive = .false.
end if
```

### Mechanism Config Changes

Update all existing mechanism configs to add `__is_advected` and TUV-x profile
properties where applicable:

**`chemistry_data/chapman_emis_dep/config.json`** (v1):
```json
{"name": "O3", ..., "__is_advected": true,
 "__is_tuvx_profile": "O3", "__tuvx_scale_height [km]": 4.5,
 "__tuvx_column_density_method": "arithmetic"}
{"name": "O",  ..., "__is_advected": true}
{"name": "O1D",..., "__is_advected": true}
{"name": "foo",..., "__is_advected": true}
{"name": "bar",..., "__is_advected": true}
{"name": "O2", ..., "__mpas_constant_vmr": 0.2095,
 "__is_tuvx_profile": "O2", "__tuvx_scale_height [km]": 8.01,
 "__tuvx_column_density_method": "geometric"}
{"name": "N2", ..., "__mpas_constant_vmr": 0.7808}
```
N2 and M keep their existing properties. O2 gains TUV-x profile metadata.
O3 gains TUV-x profile metadata. No species names appear in Fortran.

**`chemistry_data/chapman/micm/species.json`** (v0): Same pattern.
**`chemistry_data/analytical/micm/species.json`** (v0): A, B, C get `__is_advected`
(no TUV-x profiles needed for the analytical mechanism).

### Initial Conditions Config File

Each mechanism provides an optional `initial_conditions.csv` (see 1f):

```
chemistry_data/chapman_emis_dep/initial_conditions.csv
chemistry_data/chapman/initial_conditions.csv
chemistry_data/analytical/initial_conditions.csv
```

The init_atmosphere code reads this file generically at runtime. A helper script
`scripts/generate_initial_conditions.py` can produce the CSV from analytic profiles.

### Package Activation Refactor

In `mpas_atm_core_interface.F`, replace the if/else-if chain with a convention-based
approach:

```fortran
! Derive mechanism name from config path basename
call mpas_pool_get_config(configs, 'config_chemistry_config_path', config_chem_path)
! Extract basename: "chemistry_data/chapman_emis_dep" → "chapman_emis_dep"
mech_name = basename(config_chem_path)
pkg_name = 'chem_' // trim(mech_name) // '_inActive'

! Deactivate all chemistry packages, then activate the matching one
! (iterate registered packages or use known list from Registry)
nullify(pkg_ptr)
call mpas_pool_get_package(packages, trim(pkg_name), pkg_ptr)
if (associated(pkg_ptr)) then
   pkg_ptr = .true.
end if
```

This means the Fortran code never mentions specific mechanisms by name.

### Namelist Simplification

Replace the four config paths with one:

```xml
<nml_option name="config_chemistry_config_path" type="character" default_value=""
     description="Path to chemistry mechanism configuration folder"
     possible_values="Valid directory path (e.g., chemistry_data/chapman)"/>
```

Remove: `config_micm_config_path`, `config_tuvx_config_path`,
`config_tuvx_micm_mapping_path`, `config_chemistry_mechanism`.

The driver derives sub-paths:
```fortran
micm_config  = trim(config_path) // '/config.json'   ! or /micm/config.json
tuvx_config  = trim(config_path) // '/tuvx/config.json'
tuvx_mapping = trim(config_path) // '/tuvx_micm_mapping.json'
```

### `run_jw_test.sh` Simplification

Replace the mechanism-specific config path logic with:
```bash
CHEM_CONFIG="chemistry_data/${MECHANISM}"
```

And in the namelist:
```
config_chemistry_config_path = '${CHEM_CONFIG}'
```

Remove the `case` validation for mechanism names — any folder that exists is valid.

---

## 5. Verification Notebook Update

After all code changes, re-run all existing Phase 0-4 notebooks to confirm no
regressions. The Phase 4b grid mapping notebook (Section 2) serves as the new
verification artifact.

---

## Implementation Plan

All items will be implemented in this phase. Ordered by dependency.

### Step 1 — Bug fix: Temperature not passed to TUV-x
- [ ] **3a**: Fix `run_tuvx_photolysis` to pass `temp_2d` and use it as temperature

### Step 2 — Consolidate universal physical constants
- [ ] **3f**: Move `AVOGADRO`, `MW_AIR`, `SCALE_HEIGHT_AIR` into
  `mpas_chemistry_utils.F90`. Update all importers.
  Species-specific constants (`MW_O3`, `MW_O2`, `VMR_O2`, scale heights) are
  removed entirely — they become runtime values from MICM config.

### Step 3 — Dynamic arrays everywhere
- [ ] **3b**: Replace all hardcoded-size arrays with allocatable/automatic arrays
  in driver, species, emissions, deposition modules
- [ ] **3e**: Remove dead declarations from `chemistry_timestep`
- [ ] **3c**: Remove `n_deposited` from public

### Step 4 — Add config properties to all mechanism configs
- [ ] Add `__is_advected` to all advected species
- [ ] Add `__is_tuvx_profile`, `__tuvx_scale_height [km]`,
  `__tuvx_column_density_method` to O2/O3 in chapman and chapman_emis_dep
- [ ] Update `chemistry_data/analytical/micm/species.json` (A, B, C: `__is_advected`)

### Step 5 — Generic TUV-x profile discovery (species module)
- [ ] **1a**: Build `tuvx_profiles(:)` array from `__is_tuvx_profile` property.
  Eliminate `tuvx_o3_mpas_idx` and all species name checks.

### Step 6 — Generic TUV-x gas profile handling
- [ ] **1b**: Replace `o2_profile`/`o3_profile` pointers with dynamic array.
  Loop over `tuvx_profiles` in `tuvx_setup()` and `tuvx_run_column()`.
  Eliminate `MW_O3`, `MW_O2`, `VMR_O2`, scale height constants.
- [ ] **1c**: Rename `is_o3` to `use_arithmetic_mean` in `set_gas_profile`

### Step 7 — Write `scripts/generate_registry.py`
- [ ] Implement v0 and v1 config parsing
- [ ] Generate Registry.xml snippet
- [ ] Generate `mpas_atm_core_interface.F` snippet
- [ ] Verify generated output matches existing hand-written entries

### Step 8 — Namelist and path refactor
- [ ] Replace four config paths with single `config_chemistry_config_path` in
  Registry.xml
- [ ] Update driver to derive sub-paths from root
- [ ] Refactor `mpas_atm_core_interface.F` to use convention-based package activation
- [ ] Update `run_jw_test.sh` to use single config path

### Step 9 — Initial conditions from config file
- [ ] **1f**: Write generic IC reader in `mpas_init_atm_cases.F` that reads
  `initial_conditions.csv` from `config_chemistry_config_path`
- [ ] Create IC CSV files for each mechanism
- [ ] Remove all hardcoded species IC code from `mpas_init_atm_cases.F`

### Step 10 — Surface albedo namelist
- [ ] **1g**: Add `config_chemistry_surface_albedo` to Registry.xml and
  `mpas_chemistry_tuvx.F90`

### Step 11 — Grid mapping verification notebook
- [ ] **Section 2**: Create `verification/phase4b_grid_mapping.ipynb`

### Step 12 — Build, run, and regression test
- [ ] Clean build in container
- [ ] Run all three mechanisms (chapman, analytical, chapman_emis_dep)
- [ ] Re-run Phase 0-4 verification notebooks — all must pass
- [ ] Run Phase 4b grid mapping notebook — all assertions pass

### Step 13 — Commit
- [ ] Commit and push all changes to `develop-something-ambitious`
