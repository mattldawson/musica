! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Integration test: create a MIAM-enabled solver from the CAM Cloud Chemistry
! config file and verify that it produces physically reasonable results.
!
program test_miam_cloud_chemistry

  use, intrinsic :: iso_c_binding
  use iso_fortran_env, only: real64
  use musica_util, only: assert, error_t, string_t
  use musica_micm, only: micm_t, solver_stats_t
  use musica_micm, only: RosenbrockDAE4StandardOrder
  use musica_state, only: conditions_t, state_t

#include "musica/error.hpp"

#define ASSERT( expr ) call assert( expr, __FILE__, __LINE__ )
#define ASSERT_EQ( a, b ) call assert( a == b, __FILE__, __LINE__ )
#define ASSERT_GT( a, b ) call assert( a > b, __FILE__, __LINE__ )

  implicit none

  write(*,*) "=== MIAM Cloud Chemistry Fortran Integration Tests ==="

  write(*,*) "Test 1: Create solver from config..."
  call test_create_solver()

  write(*,*) "Test 2: State has cloud species..."
  call test_state_species()

  write(*,*) "Test 3: Solve and check SO4 production..."
  call test_solve_so4_production()

  write(*,*) "=== All MIAM Cloud Chemistry tests passed ==="

contains

  ! ─── Test 1: solver creation ───────────────────────────────────────────────
  subroutine test_create_solver()
    type(micm_t), pointer :: micm
    type(error_t)         :: error

    micm => micm_t("configs/v1/cam_cloud_chemistry/config.json", &
                    RosenbrockDAE4StandardOrder, error)
    ASSERT( error%is_success() )
    ASSERT( associated(micm) )

    deallocate(micm)
    write(*,*) "  PASSED"
  end subroutine test_create_solver

  ! ─── Test 2: state contains both gas and cloud species ─────────────────────
  subroutine test_state_species()
    type(micm_t), pointer  :: micm
    type(state_t), pointer :: state
    type(error_t)          :: error
    integer                :: idx

    micm => micm_t("configs/v1/cam_cloud_chemistry/config.json", &
                    RosenbrockDAE4StandardOrder, error)
    ASSERT( error%is_success() )

    state => micm%get_state(1, error)
    ASSERT( error%is_success() )

    ! Gas-phase species
    idx = state%species_ordering%index("SO2", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    idx = state%species_ordering%index("H2O2", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    idx = state%species_ordering%index("O3", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    ! Condensed-phase species (representation.phase.species)
    idx = state%species_ordering%index("CLOUD.AQUEOUS.SO2", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    idx = state%species_ordering%index("CLOUD.AQUEOUS.H2O2", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    idx = state%species_ordering%index("CLOUD.AQUEOUS.O3", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    idx = state%species_ordering%index("CLOUD.AQUEOUS.Hp", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    idx = state%species_ordering%index("CLOUD.AQUEOUS.OHm", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    idx = state%species_ordering%index("CLOUD.AQUEOUS.HSO3m", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    idx = state%species_ordering%index("CLOUD.AQUEOUS.SO3mm", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    idx = state%species_ordering%index("CLOUD.AQUEOUS.SO4mm", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    idx = state%species_ordering%index("CLOUD.AQUEOUS.H2O", error)
    ASSERT( error%is_success() )
    ASSERT_GT( idx, 0 )

    deallocate(state)
    deallocate(micm)
    write(*,*) "  PASSED"
  end subroutine test_state_species

  ! ─── Test 3: solve produces SO4 ───────────────────────────────────────────
  subroutine test_solve_so4_production()
    type(micm_t), pointer  :: micm
    type(state_t), pointer :: state
    type(string_t)         :: solver_state
    type(solver_stats_t)   :: solver_stats
    type(error_t)          :: error
    real(real64), parameter :: R_GAS = 8.31446261815324_real64
    real(real64), parameter :: T0 = 298.15_real64
    real(real64), parameter :: T_INIT = 280.0_real64
    real(real64), parameter :: P_INIT = 70000.0_real64
    real(real64), parameter :: C_H2O = 55556.0_real64
    real(real64), parameter :: GAS0_SO2 = 3.01e-8_real64
    real(real64), parameter :: GAS0_H2O2 = 3.01e-8_real64
    real(real64), parameter :: GAS0_O3  = 1.5e-6_real64
    real(real64), parameter :: SO4MM0 = 1.0_real64

    ! Henry's law constants at T0 [mol m-3 Pa-1]
    real(real64), parameter :: HLC_SO2  = 1.214e-2_real64
    real(real64), parameter :: HLC_H2O2 = 7.306e2_real64
    real(real64), parameter :: HLC_O3   = 1.135e-4_real64
    ! dH/R for temperature correction
    real(real64), parameter :: DHR_SO2  = 3120.0_real64
    real(real64), parameter :: DHR_H2O2 = 6621.0_real64
    real(real64), parameter :: DHR_O3   = 2560.0_real64
    ! Equilibrium constants at T0
    real(real64), parameter :: KA1_T0 = 3.06e-4_real64
    real(real64), parameter :: DHR_KA1 = 2090.0_real64
    real(real64), parameter :: KA2_T0 = 1.08e-9_real64
    real(real64), parameter :: DHR_KA2 = 1120.0_real64
    real(real64), parameter :: KW_T0  = 3.24e-18_real64

    real(real64) :: hlc_so2_T, hlc_h2o2_T, hlc_o3_T
    real(real64) :: alpha_SO2, alpha_H2O2, alpha_O3
    real(real64) :: Ka1_T, Ka2_T, Kw_T
    real(real64) :: hp, ohm, so2_g, so2_aq, hso3m, so3mm
    real(real64) :: h2o2_g, h2o2_aq, o3_g, o3_aq
    real(real64) :: f, hp_new
    real(real64) :: so4_initial, so4_final
    real(real64) :: time_step, total_time, dt
    integer :: i, iter
    integer :: idx_SO2, idx_H2O2, idx_O3
    integer :: idx_aq_SO2, idx_aq_H2O2, idx_aq_O3
    integer :: idx_aq_Hp, idx_aq_OHm, idx_aq_HSO3m
    integer :: idx_aq_SO3mm, idx_aq_SO4mm, idx_aq_H2O

    micm => micm_t("configs/v1/cam_cloud_chemistry/config.json", &
                    RosenbrockDAE4StandardOrder, error)
    ASSERT( error%is_success() )

    state => micm%get_state(1, error)
    ASSERT( error%is_success() )

    ! Look up species indices
    idx_SO2      = state%species_ordering%index("SO2", error)
    idx_H2O2     = state%species_ordering%index("H2O2", error)
    idx_O3       = state%species_ordering%index("O3", error)
    idx_aq_SO2   = state%species_ordering%index("CLOUD.AQUEOUS.SO2", error)
    idx_aq_H2O2  = state%species_ordering%index("CLOUD.AQUEOUS.H2O2", error)
    idx_aq_O3    = state%species_ordering%index("CLOUD.AQUEOUS.O3", error)
    idx_aq_Hp    = state%species_ordering%index("CLOUD.AQUEOUS.Hp", error)
    idx_aq_OHm   = state%species_ordering%index("CLOUD.AQUEOUS.OHm", error)
    idx_aq_HSO3m = state%species_ordering%index("CLOUD.AQUEOUS.HSO3m", error)
    idx_aq_SO3mm = state%species_ordering%index("CLOUD.AQUEOUS.SO3mm", error)
    idx_aq_SO4mm = state%species_ordering%index("CLOUD.AQUEOUS.SO4mm", error)
    idx_aq_H2O   = state%species_ordering%index("CLOUD.AQUEOUS.H2O", error)

    ! ── Compute self-consistent equilibrium initial conditions ──
    ! Temperature-adjusted HLC
    hlc_so2_T  = HLC_SO2  * exp(DHR_SO2  * (1.0_real64/T_INIT - 1.0_real64/T0))
    hlc_h2o2_T = HLC_H2O2 * exp(DHR_H2O2 * (1.0_real64/T_INIT - 1.0_real64/T0))
    hlc_o3_T   = HLC_O3   * exp(DHR_O3   * (1.0_real64/T_INIT - 1.0_real64/T0))

    alpha_SO2  = hlc_so2_T  * R_GAS * T_INIT
    alpha_H2O2 = hlc_h2o2_T * R_GAS * T_INIT
    alpha_O3   = hlc_o3_T   * R_GAS * T_INIT

    ! Temperature-adjusted equilibrium constants
    Ka1_T = KA1_T0 * exp(DHR_KA1 * (1.0_real64/T0 - 1.0_real64/T_INIT))
    Ka2_T = KA2_T0 * exp(DHR_KA2 * (1.0_real64/T0 - 1.0_real64/T_INIT))
    Kw_T  = KW_T0

    ! Simple HLC partitioning for H2O2 and O3
    h2o2_g  = GAS0_H2O2 / (1.0_real64 + alpha_H2O2)
    h2o2_aq = alpha_H2O2 * h2o2_g
    o3_g    = GAS0_O3 / (1.0_real64 + alpha_O3)
    o3_aq   = alpha_O3 * o3_g

    ! Iterate on [H+] for SO2 equilibria + charge balance
    hp = 2.0_real64 * SO4MM0
    do iter = 1, 100
      ohm = Kw_T * C_H2O * C_H2O / hp
      f = 1.0_real64 + alpha_SO2 &
          + Ka1_T * alpha_SO2 * C_H2O / hp &
          + Ka2_T * Ka1_T * alpha_SO2 * C_H2O * C_H2O / (hp * hp)
      so2_g   = GAS0_SO2 / f
      so2_aq  = alpha_SO2 * so2_g
      hso3m   = Ka1_T * so2_aq * C_H2O / hp
      so3mm   = Ka2_T * hso3m  * C_H2O / hp
      hp_new  = ohm + hso3m + 2.0_real64 * so3mm + 2.0_real64 * SO4MM0
      if (abs(hp_new - hp) < 1.0e-15_real64 * hp) exit
      hp = 0.5_real64 * (hp + hp_new)
    end do

    ! Set conditions
    state%conditions(1)%temperature  = T_INIT
    state%conditions(1)%pressure     = P_INIT
    state%conditions(1)%air_density  = P_INIT / (R_GAS * T_INIT)

    ! Set species concentrations
    associate( vs => state%species_strides%variable )
      state%concentrations(1 + (idx_SO2      - 1) * vs) = so2_g
      state%concentrations(1 + (idx_H2O2     - 1) * vs) = h2o2_g
      state%concentrations(1 + (idx_O3       - 1) * vs) = o3_g
      state%concentrations(1 + (idx_aq_SO2   - 1) * vs) = so2_aq
      state%concentrations(1 + (idx_aq_H2O2  - 1) * vs) = h2o2_aq
      state%concentrations(1 + (idx_aq_O3    - 1) * vs) = o3_aq
      state%concentrations(1 + (idx_aq_Hp    - 1) * vs) = hp
      state%concentrations(1 + (idx_aq_OHm   - 1) * vs) = ohm
      state%concentrations(1 + (idx_aq_HSO3m - 1) * vs) = hso3m
      state%concentrations(1 + (idx_aq_SO3mm - 1) * vs) = so3mm
      state%concentrations(1 + (idx_aq_SO4mm - 1) * vs) = SO4MM0
      state%concentrations(1 + (idx_aq_H2O   - 1) * vs) = C_H2O
    end associate

    ! Record initial SO4
    associate( vs => state%species_strides%variable )
      so4_initial = state%concentrations(1 + (idx_aq_SO4mm - 1) * vs)
    end associate

    ! ── Adaptive time-stepping integration for 10 s ──
    total_time = 0.0_real64
    dt = 0.01_real64
    do while (total_time < 10.0_real64 - 1.0e-10_real64)
      time_step = min(dt, 10.0_real64 - total_time)
      call micm%solve(time_step, state, solver_state, solver_stats, error)
      ASSERT( error%is_success() )
      ASSERT_EQ( solver_state%get_char_array(), "Converged" )
      total_time = total_time + time_step
      if (total_time > 0.1_real64 .and. dt < 0.1_real64) dt = 0.1_real64
      if (total_time > 1.0_real64 .and. dt < 1.0_real64) dt = 1.0_real64
    end do

    ! ── Verify SO4 increased (kinetic production) ──
    associate( vs => state%species_strides%variable )
      so4_final = state%concentrations(1 + (idx_aq_SO4mm - 1) * vs)
    end associate

    write(*,*) "  SO4 initial:", so4_initial, " final:", so4_final
    ASSERT_GT( so4_final, so4_initial )

    deallocate(state)
    deallocate(micm)
    write(*,*) "  PASSED"
  end subroutine test_solve_so4_production

end program test_miam_cloud_chemistry
