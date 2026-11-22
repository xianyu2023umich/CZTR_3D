module cz3d_app_mod
  ! The 3D convection-zone run on ARBOR (Yin-Yang shell; a radial column in
  ! one block or in nb_r-cell radial blocks, optionally multirate): CZTR
  ! step 1 in 3D. Perturbations rho1, v = (v_r,
  ! v_theta, v_phi), s1 about the fixed 1D background (rho0, s0)(r):
  !
  !   drho1/dt = -div(rho v)
  !   dv/dt    = -(v.grad) v - (grad p1 + g rho1 e_r)/rho      (with the spherical curvature terms)
  !   ds1/dt   = -v.grad s1 - v_r ds0/dr
  !
  ! The background force (dp0/dr + rho0 g) is left out: the CZ background is
  ! taken as hydrostatic (well balanced, as CZTR step 1, where it is off in
  ! every CZ cell). One step: 4-stage RK U(k) = U(0) + dt R/(5 - k); then the
  ! radial energy transport of every column (mixing length + grey,
  ! cz3d_transport_mod); then the Rempel-form diffusion with minmod slopes
  ! (switch h) in all three directions. The order and the radial operators are
  ! CZTR's: radial derivatives use CZTR's stencil and metric by radial index
  ! (cz3d_bg_mod), so a horizontally uniform state evolves as in CZTR.
  ! Walls at the bottom and the top (wall_bottom, wall_top: 'face' or 'cell'),
  ! the mass flux r^2 rho v_r odd about the wall.
  ! This file holds the module's data and the time step (the hot path); the
  ! submodules hold setup and state transfer (src/cz3d_setup_smod.f90),
  ! diagnostics (src/cz3d_diag_smod.f90) and the synthetic observables
  ! (src/cz3d_synth_smod.f90), sharing the data by host association.
  use arbor
  use cz3d_params_mod
  use pfss_mod, only: pfss_legendre, pfss_radial
  use cz3d_bg_mod, only: k_join, n_r, n_c, bg_view, r_c, r_f, drdxi_c, rho0, s0, p0, lnp0, g_c, z_c, ds0dr, rho1_init, s1_init, &
      v_init, tab, area_f, xi, T0, ddr_c, cs0, vol_c, r_top_cz
  use cz3d_transport_mod, only: transport_column, grey_atm_column, cool_f, cool_c, grey_i0, tr_view, col_bottom, col_top
  use cztr3d_atm_mod, only: atm_losses, states_interior, atm_active, atm_setup, atm_block, atm_floors, atm_dt, atm_tr_heights, S_bgf, bgforce, atm_switch, trac_f
  use eos_types, only: eos_state
  use eos_thermo, only: eos_get_from_rho_T
  use const_mod, only: Lsun
  use eos_table, only: eos_eval, eos_s_from_rho_T, eos_field, eos_field_n, eos_fields_n, n_field, F_LNP, F_LNT, F_GAMMA1, F_GAMMA3M1, F_LNCP, F_LNCV
  use eos_math, only: expm1
  use vmath_mod, only: vexp, vexpm1, vlog
  use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD
  implicit none
  private

  integer, parameter, public :: rho1_ = 1, vr_ = 2, vt_ = 3, vp_ = 4, s1_ = 5, nVarCz = 5
  ! use_B: the magnetic field (B_r, B_theta, B_phi) [G] and the GLM scalar psi
  integer, parameter, public :: br_ = 6, bt_ = 7, bp_ = 8, psi_ = 9, nVarB = 9
  integer, public :: nv = nVarCz                  ! state variables of this run (cz_n_var)
  integer :: ie_ = nVarCz + 1                     ! the diffusion's energy flux in F1_FV..Div_CV
  real(dp), parameter :: four_pi = 4.0_dp*pi

  type, extends(arbor_app_t), public :: cz3d_app_t
  contains
    procedure :: init_block => cz_init_block
    procedure :: stable_dt => cz_stable_dt
    procedure :: advance => cz_advance
    procedure :: after_step => cz_after_step
    procedure :: write_plot => cz_write_plot
    procedure :: restart_save => cz_restart_save
    procedure :: restart_load => cz_restart_load
    procedure :: flag_block => cz_flag_block
  end type cz3d_app_t

  ! per block slot: the columns' relaxing convective flux and radiative flux
  ! at the radial faces 0..nb(1) (kept between steps). The convective flux
  ! is face variable 1 of ARBOR (FaceState_FVDB(0:nb(1), j, k, 1, 1, ib):
  ! saved in restart files, moved with blocks), or qconv_FB (qconv_face false)
  real(dp), allocatable :: qconv_FB(:,:,:,:), qrad_FB(:,:,:,:)
  ! energy flux of the numerical diffusion at the radial faces of each
  ! column, from the last step (diagnostic)
  real(dp), allocatable :: qdif_FB(:,:,:,:)
  real(dp), pointer :: f_mlt(:) => null()        ! subgrid factor of the mixing-length flux, faces 0:n_r (the block's level)
  real(dp), pointer :: f_mlt_g(:,:) => null()     ! (0:n_r, level 0:max(refine_h, patch_level))
  integer :: lev_view = 0                         ! the bound block's horizontal level
  real(dp), allocatable :: rho0_w(:), rq_w(:)      ! rho0, rq of the whole column, cells -1:n_r+2 (wall ghosts)
  real(dp), allocatable :: s0_w(:), r_w(:)         ! s0, r_c of the whole column (outflow ghosts)
  ! import: the export files' index (import_open)
  integer, allocatable :: imp_key(:,:), imp_file(:)
  integer(i8), allocatable :: imp_pos(:)
  integer :: imp_nb(3) = 0, imp_nv = 0
  real(dp) :: dtheta = 0.0_dp                     ! angular cell size (Yin-Yang: the same in theta and phi)
  ! rq(i): the physical density perturbation is rq rho1 (rsst_form =
  ! 'hotta2015': 1/xi^2, rho1 being the scaled xi^2 rho1; else 1). The
  ! inertia, mass and momentum fluxes use rho0 + rho1, the EOS and buoyancy
  ! rho0 + rq rho1.
  real(dp), pointer :: rq(:) => null(), rq_g(:) => null()
  ! local time stepping (dt1 > 0): shell i advances tsc(i) dt per step
  ! (tsc >= 1, from the largest stable rate rate_r of the shell, which
  ! stable_dt collects); tloc(i) the shell's own clock (in restart files)
  real(dp), pointer, public :: tsc(:) => null()
  real(dp), pointer :: tsc_g(:) => null()
  real(dp), allocatable, public :: tloc(:)
  real(dp), pointer :: rate_r(:) => null(), rate_r_g(:) => null()
  integer :: n_adv = 0                            ! steps advanced in this run (the dt1 ramp)
  real(dp) :: t_damp0 = -1.0_dp                   ! the run's start time (the friction schedule in model time)
  real(dp), public :: dt1_now = 0.0_dp            ! the current dt1 (with the ramp)

  ! si_radial (see cz3d_params_mod): S per cell, the stored RK tendencies of
  ! the blocks in the implicit range, the diffusion's v_r before/change; the
  ! range: global cells si_i0..si_i1 (radial blocks si_kb0..si_kb1); the local
  ! columns: si_cb(q, c) the local block of column c's q-th radial block in
  ! the range
  real(dp), allocatable :: siS_B(:,:,:,:), siR_B(:,:,:,:,:), siV_B(:,:,:,:), siX_B(:,:,:,:)
  ! b_diff_form = 'resistive': the cells' eta [cm^2/s] (with one ghost layer)
  real(dp), allocatable :: EtaB_GB(:,:,:,:,:)
  real(dp), allocatable :: siSa_B(:,:,:,:)   ! si_alfven: S of v_theta, v_phi (the radial Alfven speed)
  integer :: si_i0 = 0, si_i1 = 0, si_kb0 = 0, si_kb1 = -1, si_ncol = 0
  logical :: si_cols = .false.                   ! the columns distributed (cz_columns_to_ranks)
  integer, allocatable :: si_cb(:,:)

  ! one-block scratch
  real(dp), allocatable :: rho_G(:,:,:), p1_G(:,:,:), cs_G(:,:,:), W_G(:,:,:), D_C(:,:,:), R_CV(:,:,:,:)
  real(dp), allocatable :: Dr_C(:,:,:), F1_FV(:,:,:,:), F2_FV(:,:,:,:), F3_FV(:,:,:,:), Div_CV(:,:,:,:)
  real(dp), allocatable :: lim_G(:,:,:), T_G(:,:,:), hP_G(:,:,:)
  real(dp), allocatable :: C1_F(:,:,:), C2_F(:,:,:), C3_F(:,:,:)   ! the diffusion's face speeds
  ! use_B scratch: E = v x B (or B itself) with ghosts, a curl, J = curl B, div B, the Alfven speed
  real(dp), allocatable :: B_GV(:,:,:,:), C_CV(:,:,:,:), J_CV(:,:,:,:), divB_C(:,:,:), vA_G(:,:,:)
  logical :: b_gv_zeroed = .false.                ! (rhs_mhd: B_GV's unused edge ghosts set to 0 once)
  real(dp), allocatable :: dB_C(:,:,:,:)          ! b_diff_form = 'resistive': the cells' central differences of B (diffusion: bcell_derivs)

  ! diagnostics
  integer, public :: i_prof = 0
  integer :: u_log = 0
  logical, public :: restarted = .false.      ! continuing a run: append to log.txt
  real(dp), public :: t_next_prof = 0.0_dp, t_prof0 = 0.0_dp   ! profile i at t_prof0 + i dt_prof
  real(dp) :: wall0 = 0.0_dp
  real(dp) :: dt_last = 0.0_dp
  real(dp) :: dt_now = 0.0_dp                     ! the step being taken (the GLM speed 'cfl')
  ! angular-momentum correction: the reference L_z (its start value plus the
  ! Coriolis impulse since), the impulse accumulated since the last
  ! correction, the next correction time
  real(dp) :: lz_ref = 0.0_dp, lz_imp = 0.0_dp, t_next_lzc = 0.0_dp
  ! synthetic EUV (euv_write): the response table (uniform log10 T rows) and the schedule
  real(dp), allocatable :: euv_K(:,:), euv_lgT(:)
  integer, allocatable :: euv_ch(:)
  real(dp) :: t_next_euv = -huge(1.0_dp)
  ! synthetic continuum (phot_write): the wavelengths [A], their channel and weight, the channel ids
  real(dp), allocatable :: ph_lam(:), ph_wgt(:)
  integer, allocatable :: ph_chn(:), ph_id(:)
  logical :: lz_ref_set = .false.
  ! wall-clock timers per stage of a step (tic/toc; timers_report at the end of the run)
  integer, parameter :: T_ADV = 1, T_GHOST = 2, T_RHS = 3, T_SI = 4, T_TRANS = 5, T_ATM = 6, T_DIFF = 7, T_DT = 8, &
      T_AFTER = 9, T_D_ETA = 10, T_D_EOS = 11, T_D_FA = 12, T_D_FB = 13, T_D_FR = 14, T_D_FE = 15, T_D_DIV = 16, &
      T_D_UPD = 17, T_D_IMP = 18, T_D_PRE = 19, T_D_SIP = 20, T_R_EOS = 21, T_R_MASS = 22, T_R_ADV = 23, &
      T_R_PRS = 24, T_R_MCURL = 25, T_R_MLOOP = 26, T_R_UPD = 27, T_R_SIC = 28, T_G_RK = 29, T_G_ATM = 30, &
      T_G_DIF = 31, T_G_JOIN = 32, T_G_TR = 33, T_G_ETA = 34, T_R_OVL = 35, T_A_OVL = 36, T_S_SOL = 37, n_tim = 37
  character(len=24), parameter :: tim_name(n_tim) = [character(len=24) :: 'advance (total)', &
      '  ghost fill + join', '  rhs + RK update', '  radial implicit solves', '  energy transport', &
      '  atmosphere', '  diffusion', 'stable_dt (local)', 'after_step (outputs)', &
      '    eta (+ its ghosts)', '    EOS (derived)', '    faces A + limiter', '    faces B', &
      '    walls, radial E flux', '    energy flux (E)', '    flux divergence', '    cell update', &
      '    radial B implicit', '    weights, B derivatives', '    si_diff_precond', &
      '    EOS (derived, p1)', '    mass flux div', '    advection (4 vars)', '    pressure, curvature', &
      '    MHD: curls, div B', '    MHD: Lorentz, psi', '    RK update', '    si_coeff', &
      '    fill: 4 RK stages', '    fill: atmosphere', '    fill: diffusion', '    join_bc (all)', &
      '    fill: transport', '    fill: eta', '    during the fill: EOS', '    during the fill: losses', '    the column solves']
  real(dp) :: tim_acc(n_tim) = 0.0_dp, tim_t0(n_tim) = 0.0_dp
  ! the dt report (cz_log): per constraint (1 CFL sum, 2 wave cap, 3 losses)
  ! this step's smallest dt on this rank and its cell: (c, 1) dt, 2 global
  ! cell, 3 z [km], 4 |v_r|, 5 the radial wave speed in the CFL, 6 the
  ! horizontal fast speed, 7 the dominant CFL term (1 radial, 2 theta, 3 phi),
  ! 8 dr [km]; reset after each step (cz_after_step)
  real(dp) :: dtl(3, 8) = huge(1.0_dp)
  real(dp) :: dt_run = 0.0_dp            ! model time advanced in this run (floors_report)
  real(dp) :: eta_cap_n(2) = 0.0_dp      ! b_eta_radial = 'explicit': cells capped, cells (summed over steps)

  ! ---- energy budget per radial shell (dt_budget > 0) ----
  ! E = rho e + rho v^2/2 + rho Phi summed over each shell (cells in the
  ! Yin-Yang overlap weighted 1/2, margin cells 0: an integral over the
  ! sphere). Per operator, its energy change accumulated since the run's
  ! start, to first order per step: d(rho e) = rho T ds + h d rho (T, h of
  ! the state before the operator; h = h0(r) + T s1 + p1/rho, h0 from the EOS
  ! core), the kinetic change exactly, Phi d rho:
  !   1 hydro (the RK step), 2 radiation (incl. L_bottom), 3 cooling layer
  !   (exact), 4 numerical diffusion; 5 the RSST source
  !   S = (1 - 1/xi^2)(h + v^2/2 + Phi) div(rho v), over a step
  !   -(xi^2 - 1)(h + v^2/2 + Phi) drho with the RK step's drho (part of 1:
  !   energy the RSST system creates, 0 without RSST; summed over the domain,
  !   1 - 5 is the error of the non-conservative hydro scheme).
  logical :: budget_on = .false.
  ! options tested in cell and face loops, as flags (set in cz_setup: a string
  ! comparison per cell cost ~10% of a step)
  logical :: rsst_h2012 = .false., rsst_h2015 = .false.
  real(dp), pointer :: h0_r(:) => null(), phi_r(:) => null(), h0_r_g(:) => null(), phi_r_g(:) => null()
  real(dp), allocatable :: wyy_B(:,:,:), acc_IV(:,:), E0_I(:)
  ! the current block's radial offset (cells below it) and walls (cz_bind)
  integer :: r_off = 0, grey_i0_g = 1
  logical :: at_bottom = .true., at_top = .true.
  ! the bound block lies above the join (use_atmosphere; the join is a radial
  ! block face, so a block is all CZ or all atmosphere): CZ-only work skipped
  logical :: in_atm = .false.
  ! per bound block (cz_bind): atm_switch(z) of its radial cells and sin theta
  ! of its theta cells, ghosts included (cell_cap's ingredients, once per bind)
  real(dp), allocatable :: sw_r(:), sn_t(:)
  ! the radial size in the Alfven cap's dx_min (cz_bind): drdxi_c, or huge
  ! where va_cap_horizontal leaves it out (the semi-implicit radial range)
  real(dp), allocatable :: drc_r(:)
  logical :: glm_cfl = .false.                    ! glm_ch == 'cfl' (a flag for the cell loops)
  ! derived's outputs (need: a sum of them) and cell sets
  integer, parameter :: D_P1 = 1, D_CS = 2, D_T = 4
  integer, parameter :: C_ALL = 1, C_FACES1 = 2, C_INTERIOR = 3, C_SHELL = 4
  real(dp), allocatable :: p1_IB(:,:,:,:)        ! the RK stages' interior p1 per block, computed during the ghost fill
  logical :: eos_fresh = .false.                  ! si_coeff's derived(.true.) is this state's: rhs reuses rho_G, p1_G
  ! freeze_cz: each block's state at its first step (cells r_off + i <= k_join restored after every step)
  real(dp), allocatable :: frz_B(:,:,:,:,:)
  real(dp), allocatable :: fzb_B(:,:,:,:,:)            ! test_freeze_b: the stored B, psi
  logical :: join_wall_lo = .false., join_wall_hi = .false.   ! (cz_bind: the join as a face wall for the mass flux, ghost mode)
  logical, allocatable :: frz_have(:)
  real(dp), allocatable :: work_CB(:,:,:,:)       ! W per cell from the step's first RK stage
  real(dp), allocatable :: workh_CB(:,:,:,:)      ! its horizontal pressure part -v_h.grad_h p1
  logical :: want_work = .false.
  integer :: u_bsp = 0, u_bkh = 0
  real(dp) :: t_bud0 = 0.0_dp, t_next_bud = 0.0_dp
  integer :: u_bud = 0, u_bsum = 0
  ! 6: work of the background force the momentum equation leaves out,
  !    v_r (dp0/dr + rho0 g) dV dt (our system's energy gains it relative to
  !    the full equations); 7: the fluctuating part of 5, with h' + v^2/2
  !    instead of h + v^2/2 + Phi (5 - 7: the mean part, h0 + Phi)
  ! accumulators 8-13 split the hydro operator into kinetic and thermal
  ! energy (budget_split.dat): 8 dE_kin, 9 dE_th (rho T ds + (h + Phi) drho),
  ! 10 the pressure and buoyancy work W dt (W = -v.grad p1 - g rho1 v_r, the
  ! scheme's stencils, old state), 11 the RSST source's kinetic part, 12 and 13
  ! the resolved radial energy flux at the cell centres times dt (total
  ! rho v_r (h + v^2/2 + Phi), kinetic rho v_r v^2/2), over the sphere
  ! 14-18: the horizontal kinetic energy rho (v_theta^2 + v_phi^2)/2 per shell
  ! (budget_kh.dat): 14 its change by the hydro step, 15 the horizontal
  ! pressure work -v_h.grad_h p1 dt, 16 the radial flux rho v_r v_h^2/2 at the
  ! cell centres times dt (over the sphere), 17 its change by the diffusion,
  ! 18 the total kinetic energy's change by the diffusion
  integer, parameter :: nAcc = 18
  real(dp), allocatable :: bgf_r(:)

  public :: timers_report, euv_write, cz_smooth_b, cz_columns_to_ranks, cz_import_cz, cz_match_base, cz_seed_pfss, max_mach_report, cz_n_var, cz_field_from_hydro, export_state_write, import_open, import_max_level, import_criterion, import_meta, cz_setup, cz_transport_all, cz_log, cz_profile, wall_bc, state_checksum, budget_start, dump_shells_write

  ! The routines in the submodules (src/cz3d_setup_smod.f90, cz3d_diag_smod.f90,
  ! cz3d_synth_smod.f90) called from outside their submodule: their interfaces
  ! the parent's routines the submodules call: public, so the compiler keeps
  ! their symbols (a private procedure inlined at all its calls in the parent
  ! is dropped, and the submodules' references fail to link)
  public :: tic, toc, cz_bind, cz_unbind, si_in, ddr_block, curl_block, div_block, lorentz_f, cell_cap, yy_core, angular_momentum

  interface
    ! --- cz3d_setup_smod
    module integer function cz_n_var()
    end function cz_n_var
    module subroutine cz_flag_block(self, iBlock, flags)
      class(cz3d_app_t), intent(inout) :: self
      integer, intent(in) :: iBlock
      integer, intent(out) :: flags(3)
    end subroutine cz_flag_block
    module subroutine cz_setup()
    end subroutine cz_setup
    module subroutine cz_init_block(self, iBlock)
      class(cz3d_app_t), intent(inout) :: self
      integer, intent(in) :: iBlock
    end subroutine cz_init_block
    module subroutine cz_field_from_hydro()
    end subroutine cz_field_from_hydro
    module subroutine cz_seed_pfss()
    end subroutine cz_seed_pfss
    module subroutine export_state_write(t, iStep)
      real(dp), intent(in) :: t
      integer, intent(in) :: iStep
    end subroutine export_state_write
    module subroutine import_open()
    end subroutine import_open
    module integer function import_max_level() result(l)
    end function import_max_level
    module subroutine import_criterion(key, x_lo, x_hi, flags)
      integer, intent(in) :: key(KEY_LEN)
      real(dp), intent(in) :: x_lo(3), x_hi(3)
      integer, intent(out) :: flags(3)
    end subroutine import_criterion
    module subroutine cz_import_cz()
    end subroutine cz_import_cz
    module subroutine cz_match_base()
    end subroutine cz_match_base
    module subroutine import_meta(t, iStep, found)
      real(dp), intent(inout) :: t
      integer, intent(inout) :: iStep
      logical, intent(out) :: found
    end subroutine import_meta
    module subroutine cz_columns_to_ranks(tree)
      type(arbor_tree_t), intent(in) :: tree
    end subroutine cz_columns_to_ranks
    module subroutine cz_smooth_b()
    end subroutine cz_smooth_b
    module subroutine cz_restart_save(self, Data_I)
      class(cz3d_app_t), intent(inout) :: self
      real(dp), allocatable, intent(out) :: Data_I(:)
    end subroutine cz_restart_save
    module subroutine cz_restart_load(self, Data_I)
      class(cz3d_app_t), intent(inout) :: self
      real(dp), intent(in) :: Data_I(:)
    end subroutine cz_restart_load
    ! --- cz3d_diag_smod
    module subroutine budget_start(t)
      real(dp), intent(in) :: t
    end subroutine budget_start
    module subroutine budget_hydro(ib, dt)
      integer, intent(in) :: ib
      real(dp), intent(in) :: dt
    end subroutine budget_hydro
    module subroutine cz_after_step(self, iStep, t)
      class(cz3d_app_t), intent(inout) :: self
      integer, intent(in) :: iStep
      real(dp), intent(in) :: t
    end subroutine cz_after_step
    module subroutine cz_log(iStep, t)
      integer, intent(in) :: iStep
      real(dp), intent(in) :: t
    end subroutine cz_log
    module subroutine cz_profile(iStep, t)
      integer, intent(in) :: iStep
      real(dp), intent(in) :: t
    end subroutine cz_profile
    module subroutine dump_shells_write()
    end subroutine dump_shells_write
    module function state_checksum() result(c)
      integer(i8) :: c
    end function state_checksum
    module subroutine max_mach_report(iStep, t)
      integer, intent(in) :: iStep
      real(dp), intent(in) :: t
    end subroutine max_mach_report
    module subroutine cz_write_plot(self, dom, iStep, t)
      class(cz3d_app_t), intent(inout) :: self
      type(arbor_domain_t), intent(in) :: dom
      integer, intent(in) :: iStep
      real(dp), intent(in) :: t
    end subroutine cz_write_plot
    module subroutine timers_report(wall, n_step)
      real(dp), intent(in) :: wall
      integer, intent(in) :: n_step
    end subroutine timers_report
    ! --- cz3d_synth_smod
    module subroutine euv_write(iStep, t)
      integer, intent(in) :: iStep
      real(dp), intent(in) :: t
    end subroutine euv_write
  end interface

contains

  ! ================================================================ block binding (setup: cz3d_setup_smod)

  ! The radial view of block ib: the 1D arrays (background, metric, step
  ! factors, subgrid factor, cooling) indexed by the block's local cell 1..nb(1),
  ! and its walls. cz_unbind: the whole column (r_off = 0, both walls).
  subroutine cz_bind(ib)
    integer, intent(in) :: ib
    integer :: i
    lev_view = Key_IB(KL+2, ib)
    call column_views(Key_IB(KI+1, ib)*nb(1), nb(1))
    at_bottom = arbor_physical_side(ib, 1, 1)
    at_top = arbor_physical_side(ib, 1, 2)
    col_bottom = at_bottom
    col_top = at_top
    if (use_atmosphere) col_top = Key_IB(KI+1, ib)*nb(1) + nb(1) == k_join   ! (grey and MLT end at the join)
    in_atm = use_atmosphere .and. Key_IB(KI+1, ib)*nb(1) >= k_join
    if (.not. allocated(sw_r)) allocate(sw_r(1-ng:nb(1)+ng), sn_t(1-ng:nb(2)+ng), drc_r(1-ng:nb(1)+ng))
    drc_r = drdxi_c(1-ng:nb(1)+ng)
    if (va_cap_horizontal .and. si_radial .and. si_i1 > 0) then
      do i = 1 - ng, nb(1) + ng
        if (r_off + i >= si_i0 + 2 .and. r_off + i <= si_i1 - 2) drc_r(i) = huge(1.0_dp)
      end do
    end if
    sw_r = 0.0_dp
    if (use_atmosphere) sw_r = atm_switch(z_c(1-ng:nb(1)+ng))
    sn_t = sin(Metric_DB(2, ib)%X_I(1-ng:nb(2)+ng))
    join_wall_lo = freeze_cz_mode == 'ghost' .and. freeze_side == 'cz' .and. Key_IB(KI+1, ib)*nb(1) == k_join
    join_wall_hi = freeze_cz_mode == 'ghost' .and. freeze_side == 'atmosphere' .and. Key_IB(KI+1, ib)*nb(1) + nb(1) == k_join
  end subroutine cz_bind

  subroutine cz_unbind()
    lev_view = refine_h
    call column_views(0, n_r)
    at_bottom = .true.
    at_top = .true.
    col_bottom = .true.
    col_top = .true.
    in_atm = .false.
    join_wall_lo = .false.
    join_wall_hi = .false.
  end subroutine cz_unbind

  ! freeze_side: is block ib on the frozen side (the join is a block face)?
  logical function block_frozen(ib)
    integer, intent(in) :: ib
    integer :: lo
    lo = Key_IB(KI+1, ib)*nb(1)
    select case (freeze_side)
    case ('cz');         block_frozen = lo + nb(1) <= k_join
    case ('atmosphere'); block_frozen = lo >= k_join
    case default;        block_frozen = .false.
    end select
  end function block_frozen

  ! ghost mode: the frozen side's blocks are not computed
  logical function skip_frozen(ib)
    integer, intent(in) :: ib
    skip_frozen = freeze_cz_mode == 'ghost' .and. block_frozen(ib)
  end function skip_frozen

  ! freeze_cz_mode = 'ghost': after a ghost fill, the evolving side's radial
  ! ghost cells at the join (see freeze_cz_mode): v odd, B_theta, B_phi, psi
  ! even; for freeze_side = 'atmosphere' also rho1 (relative to rho0) and s1
  ! even (a wall); B_r (and for 'cz' rho1, s1) keep the frozen side's values
  ! from the fill
  subroutine join_bc()
    integer :: jb, ib, g, d, s, m, o, v
    if (freeze_side == 'none' .or. freeze_cz_mode /= 'ghost') return
    m = nb(1)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      o = Key_IB(KI+1, ib)*m
      if (freeze_side == 'cz' .and. o == k_join) then
        do g = 1, 2
          d = 1 - g;  s = g                           ! ghost 0 <- cell 1, ghost -1 <- cell 2
          State_GVB(d, :, :, vr_:vp_, ib) = -State_GVB(s, :, :, vr_:vp_, ib)
          if (use_B) State_GVB(d, :, :, bt_:psi_, ib) = State_GVB(s, :, :, bt_:psi_, ib)
        end do
      else if (freeze_side == 'atmosphere' .and. o + m == k_join) then
        do g = 1, 2
          d = m + g;  s = m + 1 - g                   ! ghost m+1 <- cell m, ghost m+2 <- cell m-1
          State_GVB(d, :, :, vr_:vp_, ib) = -State_GVB(s, :, :, vr_:vp_, ib)
          if (use_B) State_GVB(d, :, :, bt_:psi_, ib) = State_GVB(s, :, :, bt_:psi_, ib)
          State_GVB(d, :, :, s1_, ib) = State_GVB(s, :, :, s1_, ib)
          State_GVB(d, :, :, rho1_, ib) = State_GVB(s, :, :, rho1_, ib) &
              *(rq_w(o+s)*rho0_w(o+d))/(rq_w(o+d)*rho0_w(o+s))
        end do
      end if
    end do
  end subroutine join_bc

  subroutine column_views(o, m)
    integer, intent(in) :: o, m
    call bg_view(o, m)
    call tr_view(o, m)
    if (associated(tsc_g)) then
      tsc(1:m) => tsc_g(o+1:o+m)
      rate_r(1:m) => rate_r_g(o+1:o+m)
    end if
    if (associated(rq_g)) rq(-1:m+2) => rq_g(o-1:o+m+2)
    if (associated(h0_r_g)) then
      h0_r(-1:m+2) => h0_r_g(o-1:o+m+2)
      phi_r(-1:m+2) => phi_r_g(o-1:o+m+2)
    end if
    if (associated(f_mlt_g)) f_mlt(0:m) => f_mlt_g(o:o+m, min(lev_view, ubound(f_mlt_g, 2)))
    r_off = o
    grey_i0 = grey_i0_g - o
  end subroutine column_views

  ! the stage timers (see tim_name)
  subroutine tic(k)
    integer, intent(in) :: k
    tim_t0(k) = arbor_wtime()
  end subroutine tic
  subroutine toc(k)
    integer, intent(in) :: k
    tim_acc(k) = tim_acc(k) + (arbor_wtime() - tim_t0(k))
  end subroutine toc

  ! ================================================================ boundaries

  ! Walls in r: v_r odd, everything else even; 'cell': v_r = 0 in the
  ! boundary cell. (theta, phi: Yin-Yang patch edges, filled by ARBOR.)
  subroutine wall_bc(A, iVarMin, iVarMax, iBlock, iDir, iSide, nGhost)
    real(dp), intent(inout) :: A(1-ng:, 1-ng:, 1-ng:, :, :)
    integer, intent(in) :: iVarMin, iVarMax, iBlock, iDir, iSide, nGhost
    integer :: v, wall
    character(len=8) :: wt
    character(len=12) :: bc
    logical :: odd
    wt = merge(wall_bottom, wall_top, iSide == 1)
    if (iDir == 1 .and. iSide == 2 .and. wt == 'outflow') then
      call outflow_top(Key_IB(KI+1, iBlock)*nb(1))
      return
    end if
    wall = merge(WALL_FACE, WALL_CELL, wt == 'face')
    bc = merge(b_bc_bottom, b_bc_top, iSide == 1)
    do v = iVarMin, iVarMax
      ! odd about the wall: v_r; B_r at a conductor; B_theta, B_phi for a radial field
      odd = iDir == 1 .and. (v == vr_ .or. (v == br_ .and. bc == 'conductor') &
          .or. ((v == bt_ .or. v == bp_) .and. bc == 'radial'))
      if (v > nv) cycle
      if (odd) then
        if (wall == WALL_CELL) call arbor_bc_set_cell(A, v, v, iBlock, iDir, iSide, 0.0_dp)
        call arbor_bc_mirror(A, v, v, iBlock, iDir, iSide, nGhost, -1, wall)
      else
        call arbor_bc_mirror(A, v, v, iBlock, iDir, iSide, nGhost, +1, wall)
      end if
      ! the density: its relative perturbation rq rho1/rho0 mirrored (the
      ! background changes steeply across the wall; mirroring rho1 itself
      ! gives negative ghost densities above an evacuated top cell)
      if (v == rho1_ .and. iDir == 1 .and. allocated(rho0_w)) call relative_rho1(Key_IB(KI+1, iBlock)*nb(1))
    end do
  contains
    ! wall_top = 'outflow': rho ~ r^-2, v_r >= 0 copied, v_h and s copied;
    ! B_r ~ r^-2, B_h and psi copied
    subroutine outflow_top(off)
      integer, intent(in) :: off
      integer :: g, m, dst, v
      real(dp) :: q2
      m = nb(1)
      do g = 1, nGhost
        dst = m + g
        q2 = (r_w(off+m)/r_w(off+dst))**2
        do v = iVarMin, min(iVarMax, nv)
          if (v == rho1_) then
            A(dst, :, :, v, iBlock) = ((rho0_w(off+m) + rq_w(off+m)*A(m, :, :, v, iBlock))*q2 - rho0_w(off+dst)) &
                /rq_w(off+dst)
          else if (v == vr_) then
            A(dst, :, :, v, iBlock) = max(A(m, :, :, v, iBlock), 0.0_dp)
          else if (v == s1_) then
            A(dst, :, :, v, iBlock) = A(m, :, :, v, iBlock) + s0_w(off+m) - s0_w(off+dst)
          else if (use_B .and. v == br_) then
            A(dst, :, :, v, iBlock) = A(m, :, :, v, iBlock)*q2
          else
            A(dst, :, :, v, iBlock) = A(m, :, :, v, iBlock)
          end if
        end do
      end do
    end subroutine outflow_top

    subroutine relative_rho1(off)
      integer, intent(in) :: off
      integer :: g, dst, src
      do g = 1, nGhost
        if (iSide == 1) then
          dst = 1 - g
          src = merge(g, 1 + g, wall == WALL_FACE)
        else
          dst = nb(1) + g
          src = merge(nb(1) + 1 - g, nb(1) - g, wall == WALL_FACE)
        end if
        A(dst, :, :, rho1_, iBlock) = A(dst, :, :, rho1_, iBlock) &
            *(rq_w(off+src)*rho0_w(off+dst))/(rq_w(off+dst)*rho0_w(off+src))
      end do
    end subroutine relative_rho1
  end subroutine wall_bc

  ! The radial mass flux W = r^2 rho v_r odd about both walls (mirroring
  ! rho1 and v_r alone would leak: r^2 rho0 differs across the wall)
  subroutine odd_wall_flux(W)
    real(dp), intent(inout) :: W(1-ng:, 1-ng:, 1-ng:)
    integer :: m
    m = nb(1)
    if (at_bottom) then
      if (wall_bottom == 'face') then
        W(0, :, :) = -W(1, :, :);  W(-1, :, :) = -W(2, :, :)
      else
        W(1, :, :) = 0.0_dp;  W(0, :, :) = -W(2, :, :);  W(-1, :, :) = -W(3, :, :)
      end if
    end if
    if (join_wall_lo) then                           ! (the join, ghost mode: a face wall for the evolving side)
      W(0, :, :) = -W(1, :, :);  W(-1, :, :) = -W(2, :, :)
    end if
    if (join_wall_hi) then
      W(m+1, :, :) = -W(m, :, :);  W(m+2, :, :) = -W(m-1, :, :)
    end if
    if (at_top .and. wall_top /= 'outflow') then
      if (wall_top == 'face') then
        W(m+1, :, :) = -W(m, :, :);  W(m+2, :, :) = -W(m-1, :, :)
      else
        W(m, :, :) = 0.0_dp;  W(m+1, :, :) = -W(m-1, :, :);  W(m+2, :, :) = -W(m-2, :, :)
      end if
    end if
  end subroutine odd_wall_flux

  ! ================================================================ equation of state

  ! The equation of state on block ib: rho_G (the inertia density rho0 + rho1)
  ! on every cell, and the outputs in need (a sum of D_P1, D_CS, D_T) on the
  ! cells of cells:
  !   D_P1  p1_G = p0 expm1(ln p - ln p0)
  !   D_CS  cs_G = sqrt(Gamma1 p/rho)/xi (+ the capped v_A with the field)
  !   D_T   T_G
  !   C_ALL      every cell, both ghost layers, edges and corners included
  !   C_FACES1   the interior and one ghost layer across each face
  !   C_INTERIOR the interior cells
  !   C_SHELL    the ghost cells only (C_ALL less the interior)
  ! Only the table fields the outputs need are interpolated (ln p for D_P1 or
  ! D_CS, Gamma1 for D_CS, ln T for D_T), a radial segment of a column at a
  ! time (eos_hhe's column evaluators: the same values as the scalar calls,
  ! whatever the segment).
  subroutine derived(ib, need, cells)
    integer, intent(in) :: ib, need, cells
    real(dp) :: rp(1-ng:nb(1)+ng), s(1-ng:nb(1)+ng), q(1-ng:nb(1)+ng, 3), d(1-ng:nb(1)+ng), &
        e1(1-ng:nb(1)+ng), e2(1-ng:nb(1)+ng), e3(1-ng:nb(1)+ng)
    integer :: i, j, k, nf, ip, iT, ig, ids(3), seg(2, 2), ns, m
    logical :: wp1, wcs, wT, jin, kin
    wp1 = iand(need, D_P1) /= 0
    wcs = iand(need, D_CS) /= 0
    wT = iand(need, D_T) /= 0
    ! the table fields and their columns in q
    nf = 0;  ip = 0;  iT = 0;  ig = 0
    if (wp1 .or. wcs) then
      nf = nf + 1;  ids(nf) = F_LNP;  ip = nf
    end if
    if (wT) then
      nf = nf + 1;  ids(nf) = F_LNT;  iT = nf
    end if
    if (wcs) then
      nf = nf + 1;  ids(nf) = F_GAMMA1;  ig = nf
    end if
    do k = 1 - ng, nb(3) + ng
      do j = 1 - ng, nb(2) + ng
        do i = 1 - ng, nb(1) + ng
          rho_G(i, j, k) = rho0(i) + State_GVB(i, j, k, rho1_, ib)
          rp(i) = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)                     ! physical density
          s(i) = s0(i) + State_GVB(i, j, k, s1_, ib)
        end do
        if (nf == 0) cycle
        ! the column's radial segments
        jin = j >= 1 .and. j <= nb(2)
        kin = k >= 1 .and. k <= nb(3)
        ns = 1
        select case (cells)
        case (C_ALL)
          seg(:, 1) = [1 - ng, nb(1) + ng]
        case (C_FACES1)
          if (jin .and. kin) then                     ! interior column: its radial face ghosts too
            seg(:, 1) = [0, nb(1) + 1]
          else if ((jin .and. (k == 0 .or. k == nb(3) + 1)) .or. (kin .and. (j == 0 .or. j == nb(2) + 1))) then
            seg(:, 1) = [1, nb(1)]                    ! a horizontal face ghost column
          else
            cycle
          end if
        case (C_INTERIOR)
          if (.not. (jin .and. kin)) cycle
          seg(:, 1) = [1, nb(1)]
        case (C_SHELL)
          if (jin .and. kin) then                     ! interior column: its radial ghosts
            ns = 2
            seg(:, 1) = [1 - ng, 0]
            seg(:, 2) = [nb(1) + 1, nb(1) + ng]
          else
            seg(:, 1) = [1 - ng, nb(1) + ng]
          end if
        case default
          call arbor_abort('derived: unknown cells')
        end select
        do m = 1, ns
          call eval_seg(seg(1, m), seg(2, m))
        end do
      end do
    end do
  contains
    subroutine eval_seg(i0, i1)
      integer, intent(in) :: i0, i1
      integer :: n
      n = i1 - i0 + 1
      call eos_fields_n(tab, ids(1:nf), n, rp(i0:i1), s(i0:i1), q(i0:i1, 1:nf))
      if (wp1) then
        d(i0:i1) = q(i0:i1, ip) - lnp0(i0:i1)
        call vexpm1(n, d(i0:i1), e1(i0:i1))
        p1_G(i0:i1, j, k) = p0(i0:i1)*e1(i0:i1)
      end if
      if (wcs) then
        call vexp(n, q(i0:i1, ip), e2(i0:i1))
        do i = i0, i1
          cs_G(i, j, k) = sqrt(q(i, ig)*e2(i)/rp(i))/xi(i)                         ! reduced (RSST)
        end do
        ! with B: the fast speed's bound c_red + v_A (v_A capped, see va_eff; va_cap the same arithmetic)
        if (use_B) then
!GCC$ ivdep
!DIR$ IVDEP
          do i = i0, i1
            cs_G(i, j, k) = cs_G(i, j, k) + va_cap(sum(State_GVB(i, j, k, br_:bp_, ib)**2), rp(i), cs_G(i, j, k), &
                norm2(State_GVB(i, j, k, vr_:vp_, ib)), min(drc_r(i), r_c(i)*Metric_DB(2, ib)%Dx_I(j), &
                r_c(i)*sn_t(j)*Metric_DB(3, ib)%Dx_I(k)), sw_r(i))
          end do
        end if
      end if
      if (wT) then
        call vexp(n, q(i0:i1, iT), e3(i0:i1))
        T_G(i0:i1, j, k) = e3(i0:i1)
      end if
    end subroutine eval_seg
  end subroutine derived

  ! ================================================================ radial semi-implicit

  ! the variables the column solves act on: v_r, and v_theta, v_phi with si_alfven
  pure logical function si_var(v)
    integer, intent(in) :: v
    si_var = v == vr_ .or. (si_alfven .and. (v == vt_ .or. v == vp_))
  end function si_var

  logical function si_in(ib)
    integer, intent(in) :: ib
    si_in = si_radial .and. Key_IB(KI+1, ib) >= si_kb0 .and. Key_IB(KI+1, ib) <= si_kb1
  end function si_in

  ! S of block ib's cells in the implicit range for the step dt (0 elsewhere
  ! and in the two cells at each end of the range):
  ! S = max(0, c^2 dt^2/(si_cfl - C_f)^2 - dr^2/4), c the fast speed bound of
  ! derived (c_s + v_A), C_f = |v_r| dt/dr (at most si_cfl/2 in the formula);
  ! si_alfven: siSa_B the same with c the radial Alfven speed (Boris reduced)
  subroutine si_coeff(ib, dt, cells)
    integer, intent(in) :: ib
    real(dp), intent(in) :: dt
    integer, intent(in) :: cells                    ! derived's cells: C_ALL (rhs reuses rho_G, p1_G) or C_INTERIOR
    integer :: i, j, k, ilo, ihi
    real(dp) :: dr, cf, rp
    if (.not. si_cols) call arbor_abort('cz3d: si_radial: the columns were not distributed (cz_columns_to_ranks)')
    siS_B(:, :, :, ib) = 0.0_dp
    if (si_alfven) siSa_B(:, :, :, ib) = 0.0_dp
    if (.not. si_in(ib)) return
    call derived(ib, D_P1 + D_CS, cells)        ! (T is not needed here nor by rhs)
    eos_fresh = cells == C_ALL                          ! (rhs, next on this block and state, reuses rho_G and p1_G)
    ! the block's cells inside the range less its two end cells at each end
    ilo = max(1, si_i0 + 2 - r_off)
    ihi = min(nb(1), si_i1 - 2 - r_off)
    do k = 1, nb(3)
      do j = 1, nb(2)
!GCC$ ivdep
!DIR$ IVDEP
        do i = ilo, ihi
          dr = drdxi_c(i)
          cf = min(abs(State_GVB(i, j, k, vr_, ib))*dt/dr, 0.5_dp*si_cfl)
          siS_B(i, j, k, ib) = max(0.0_dp, (cs_G(i, j, k)*dt/(si_cfl - cf))**2 - 0.25_dp*dr**2)
        end do
        if (.not. si_alfven) cycle
!GCC$ ivdep
!DIR$ IVDEP
        do i = ilo, ihi
          dr = drdxi_c(i)
          cf = min(abs(State_GVB(i, j, k, vr_, ib))*dt/dr, 0.5_dp*si_cfl)
          rp = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
          ! (cell_cs, cell_cap, va_eff through va_cap: the same arithmetic)
          siSa_B(i, j, k, ib) = max(0.0_dp, (va_cap(State_GVB(i, j, k, br_, ib)**2, rp, &
              sqrt(5.0_dp/3.0_dp*max(p0(i) + p1_G(i, j, k), 1.0e-30_dp)/rp)/xi(i), &
              norm2(State_GVB(i, j, k, vr_:vp_, ib)), min(drc_r(i), r_c(i)*Metric_DB(2, ib)%Dx_I(j), &
              r_c(i)*sn_t(j)*Metric_DB(3, ib)%Dx_I(k)), sw_r(i))*dt/(si_cfl - cf))**2 - 0.25_dp*dr**2)
        end do
      end do
    end do
  end subroutine si_coeff

  ! X <- x with (rho - (1/r^2) d/dr r^2 rho S d/dr) x = rho X along every
  ! local column's range (conservative 3-point form on the cell widths; face
  ! S and rho the means of the two cells; no flux through the range's end
  ! faces). X(i, j, k, ib) on the blocks of the range, S = SS (siS_B, or
  ! siSa_B for v_theta, v_phi); Thomas algorithm (the rows are diagonally
  ! dominant). X2 (optional): a second right-hand side with the same matrix
  ! (v_theta and v_phi), eliminated with the same multipliers.
  ! All (j, k) columns of a column stack at once, (j, k) the inner index:
  ! the recursion runs in r, each step a vector operation over the columns
  ! (the same arithmetic per column); a column whose S is 0 everywhere keeps
  ! its X (the identity, not solved)
  subroutine si_column_solve(X, SS, X2)
    real(dp), intent(inout) :: X(:, :, :, :)
    real(dp), intent(in) :: SS(:, :, :, :)
    real(dp), intent(inout), optional :: X2(:, :, :, :)
    integer :: c, q, i, j, k, m, nk, n, g, ib, njk, jk
    real(dp), allocatable :: a(:,:), b(:,:), cc(:,:), d(:,:), d2(:,:), wf(:,:), rp(:,:), sv(:,:), w(:), smax(:)
    real(dp) :: vol
    logical :: two
    two = present(X2)
    nk = si_kb1 - si_kb0 + 1
    m = nk*nb(1)
    njk = nb(2)*nb(3)
    allocate(a(njk, m), b(njk, m), cc(njk, m), d(njk, m), wf(njk, 0:m), rp(njk, m), sv(njk, m), w(njk), smax(njk))
    if (two) allocate(d2(njk, m))
    call cz_unbind()                                   ! (the global radial arrays)
    do c = 1, si_ncol
      ! the column stack's values, (j, k) inner
      do q = 1, nk
        ib = si_cb(q, c)
        do k = 1, nb(3)
          do j = 1, nb(2)
            jk = j + (k - 1)*nb(2)
            do i = 1, nb(1)
              n = (q - 1)*nb(1) + i
              g = si_i0 - 1 + n
              rp(jk, n) = rho0(g) + rq(g)*State_GVB(i, j, k, rho1_, ib)
              sv(jk, n) = SS(i, j, k, ib)
              d(jk, n) = X(i, j, k, ib)
              if (two) d2(jk, n) = X2(i, j, k, ib)
            end do
          end do
        end do
      end do
      smax = maxval(sv, dim=2)
      if (all(smax <= 0.0_dp)) cycle
      wf(:, 0) = 0.0_dp;  wf(:, m) = 0.0_dp
      do n = 1, m - 1
        g = si_i0 - 1 + n
        wf(:, n) = r_f(g)**2*0.25_dp*(rp(:, n) + rp(:, n+1))*(sv(:, n) + sv(:, n+1))/(r_c(g+1) - r_c(g))
      end do
      do n = 1, m
        g = si_i0 - 1 + n
        vol = r_c(g)**2*(r_f(g) - r_f(g-1))
        a(:, n) = -wf(:, n-1)
        cc(:, n) = -wf(:, n)
        b(:, n) = vol*rp(:, n) + wf(:, n-1) + wf(:, n)
        d(:, n) = vol*rp(:, n)*d(:, n)
        if (two) d2(:, n) = vol*rp(:, n)*d2(:, n)
      end do
      do n = 2, m
        w = a(:, n)/b(:, n-1)
        b(:, n) = b(:, n) - w*cc(:, n-1)
        d(:, n) = d(:, n) - w*d(:, n-1)
        if (two) d2(:, n) = d2(:, n) - w*d2(:, n-1)
      end do
      d(:, m) = d(:, m)/b(:, m)
      if (two) d2(:, m) = d2(:, m)/b(:, m)
      do n = m - 1, 1, -1
        d(:, n) = (d(:, n) - cc(:, n)*d(:, n+1))/b(:, n)
        if (two) d2(:, n) = (d2(:, n) - cc(:, n)*d2(:, n+1))/b(:, n)
      end do
      do q = 1, nk
        ib = si_cb(q, c)
        do k = 1, nb(3)
          do j = 1, nb(2)
            jk = j + (k - 1)*nb(2)
            if (smax(jk) <= 0.0_dp) cycle              ! (S = 0 in the whole column: X as it was)
            do i = 1, nb(1)
              n = (q - 1)*nb(1) + i
              X(i, j, k, ib) = d(jk, n)
              if (two) X2(i, j, k, ib) = d2(jk, n)
            end do
          end do
        end do
      end do
    end do
  end subroutine si_column_solve

  ! b_diff_form = 'resistive': each cell's eta = b_diff_factor Phi (|v| + c)/2 dh
  ! (see b_diff_form; b_eta_dr > 0: b_eta_dr (|v| [+ c_s]) dr/2), then its ghosts (one layer; at the radial walls the
  ! boundary cell's value). Needs the ghosts of State (two layers) filled.
  subroutine b_eta_fill(iBlocks)
    integer, intent(in) :: iBlocks(:)
    integer :: jb, ib, i, j, k, c, d, di, dj, dk, s
    real(dp) :: phi, w, du, jq, dh
    do jb = 1, size(iBlocks)
      ib = iBlocks(jb)
      call cz_bind(ib)
      if (b_eta_dr > 0.0_dp) then                        ! (CZTR_2D's eta = b_eta_dr (|v| [+ c_s]) dr/2)
        if (b_eta_cs) call derived(ib, D_P1, C_ALL)
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, nb(1)
              w = sqrt(sum(State_GVB(i, j, k, vr_:vp_, ib)**2))
              if (b_eta_cs) w = w + cell_cs(i, j, k, rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib))
              EtaB_GB(i, j, k, 1, ib) = b_eta_dr*w*0.5_dp*(r_f(i) - r_f(i-1))
            end do
          end do
        end do
        cycle
      end if
      call derived(ib, D_CS, C_ALL)
      associate (th => Metric_DB(2, ib)%X_I, dth_ => Metric_DB(2, ib)%Dx_I, dph_ => Metric_DB(3, ib)%Dx_I)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            phi = 0.0_dp
            do d = 1, 3
              di = merge(1, 0, d == 1);  dj = merge(1, 0, d == 2);  dk = merge(1, 0, d == 3)
              do s = -1, 0                               ! (the faces below and above the cell)
                do c = br_, bp_
                  associate (q1 => State_GVB(i+(s-1)*di, j+(s-1)*dj, k+(s-1)*dk, c, ib), &
                             q2 => State_GVB(i+s*di, j+s*dj, k+s*dk, c, ib), &
                             q3 => State_GVB(i+(s+1)*di, j+(s+1)*dj, k+(s+1)*dk, c, ib), &
                             q4 => State_GVB(i+(s+2)*di, j+(s+2)*dj, k+(s+2)*dk, c, ib))
                    du = q3 - q2
                    if (du == 0.0_dp) cycle
                    jq = jump4b(q1, q2, q3, q4)
                    w = max(0.0_dp, min(1.0_dp, jq/du))
                  end associate
                  phi = max(phi, w)
                end do
              end do
            end do
            dh = r_c(i)*min(dth_(j), sin(th(j))*dph_(k))
            EtaB_GB(i, j, k, 1, ib) = b_diff_factor*phi*0.5_dp*(sqrt(sum(State_GVB(i, j, k, vr_:vp_, ib)**2)) &
                + cs_G(i, j, k))*dh
          end do
        end do
      end do
      end associate
    end do
    ! b_eta_radial = 'explicit': eta capped at the explicit limit (see cz3d_params_mod)
    if (b_eta_explicit) then
      do jb = 1, size(iBlocks)
        ib = iBlocks(jb)
        call cz_bind(ib)
        associate (th => Metric_DB(2, ib)%X_I, dth_ => Metric_DB(2, ib)%Dx_I, dph_ => Metric_DB(3, ib)%Dx_I)
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, nb(1)
              w = 0.4_dp/(tsc(i)*dt_now*(1.0_dp/min(r_c(i+1) - r_c(i), r_c(i) - r_c(i-1))**2 &
                  + 1.0_dp/(r_c(i)*dth_(j))**2 + 1.0_dp/(r_c(i)*sin(th(j))*dph_(k))**2))
              eta_cap_n(2) = eta_cap_n(2) + 1.0_dp
              if (EtaB_GB(i, j, k, 1, ib) > w) then
                EtaB_GB(i, j, k, 1, ib) = w
                eta_cap_n(1) = eta_cap_n(1) + 1.0_dp
              end if
            end do
          end do
        end do
        end associate
      end do
    end if
    call cz_unbind()
    call tic(T_G_ETA)
    call arbor_fill_ghosts(EtaB_GB, n_ghost=1, faces_only=.true.)
    call toc(T_G_ETA)
    do jb = 1, size(iBlocks)
      ib = iBlocks(jb)
      if (Key_IB(KI+1, ib) == 0) EtaB_GB(0, :, :, 1, ib) = EtaB_GB(1, :, :, 1, ib)
      if ((Key_IB(KI+1, ib) + 1)*nb(1) == n_r) EtaB_GB(nb(1)+1, :, :, 1, ib) = EtaB_GB(nb(1), :, :, 1, ib)
    end do
  end subroutine b_eta_fill

  ! b_diff_form = 'resistive': the radial normal term of B_theta and B_phi
  ! implicit, along every local column (the whole column):
  !   x - x* = tsc dt (1/V) [A_+ eta_+ (r_+ x_+ - r x)/(r_f dr_+) - A_- eta_- (r x - r_- x_-)/(r_f dr_-)]
  ! (the radial flux -eta (r_+ x_+ - r x)/(r_f dr) of the explicit form, A =
  ! r_f^2, V = the cell's volume per steradian); at the column's end faces the ghost value is
  ! explicit (the fill before the diffusion). x* the state after the
  ! explicit fluxes. With diffusion_heating the field energy removed from a
  ! cell goes into its s1 (with the atmosphere: the CZ cells, up to k_join).
  subroutine b_radial_implicit(iBlocks, dt)
    integer, intent(in) :: iBlocks(:)
    real(dp), intent(in) :: dt
    integer :: c, q, i, j, k, m, nk, n, g, ib, v, jb
    real(dp), allocatable :: a(:), b(:), cc(:), d(:), et(:), kf(:), xg(:), de(:,:,:,:)
    real(dp) :: vol, w, rf, rpw(nb(1)), sw(nb(1)), ltw(nb(1)), tw(nb(1))
    logical :: heat
    nk = n_r/nb(1)
    m = n_r
    heat = diffusion_heating
    allocate(a(m), b(m), cc(m), d(m), et(0:m+1), kf(0:m), xg(2))
    if (heat) then
      allocate(de(nb(1), nb(2), nb(3), size(State_GVB, 5)))
      de = 0.0_dp
    end if
    call cz_unbind()                                   ! (the global radial arrays)
    do c = 1, si_ncol
      do k = 1, nb(3)
        do j = 1, nb(2)
          do q = 1, nk
            ib = si_cb(q, c)
            do i = 1, nb(1)
              et((q - 1)*nb(1) + i) = EtaB_GB(i, j, k, 1, ib)
            end do
          end do
          et(0) = EtaB_GB(0, j, k, 1, si_cb(1, c))
          et(m+1) = EtaB_GB(nb(1)+1, j, k, 1, si_cb(nk, c))
          if (maxval(et) <= 0.0_dp) cycle
          ! kf(n): A eta/(r_f dr) of the face above cell n (0: the bottom face)
          do n = 0, m
            rf = 0.5_dp*(r_c(n) + r_c(n+1))
            kf(n) = r_f(n)**2*0.5_dp*(et(n) + et(n+1))/(rf*(r_c(n+1) - r_c(n)))
          end do
          do v = bt_, bp_
            do q = 1, nk
              ib = si_cb(q, c)
              do i = 1, nb(1)
                d((q - 1)*nb(1) + i) = State_GVB(i, j, k, v, ib)
              end do
            end do
            xg(1) = State_GVB(0, j, k, v, si_cb(1, c))
            xg(2) = State_GVB(nb(1)+1, j, k, v, si_cb(nk, c))
            do n = 1, m
              vol = vol_c(n)/(tsc_g(n)*dt)
              a(n) = -kf(n-1)*r_c(n-1)
              cc(n) = -kf(n)*r_c(n+1)
              b(n) = vol + (kf(n-1) + kf(n))*r_c(n)
              if (heat .and. (n <= k_join .or. .not. use_atmosphere)) then
                ib = si_cb((n - 1)/nb(1) + 1, c)
                de(n - ((n - 1)/nb(1))*nb(1), j, k, ib) = de(n - ((n - 1)/nb(1))*nb(1), j, k, ib) - d(n)**2
              end if
              d(n) = vol*d(n)
            end do
            d(1) = d(1) - a(1)*xg(1);  a(1) = 0.0_dp
            d(m) = d(m) - cc(m)*xg(2);  cc(m) = 0.0_dp
            do n = 2, m
              w = a(n)/b(n-1)
              b(n) = b(n) - w*cc(n-1)
              d(n) = d(n) - w*d(n-1)
            end do
            d(m) = d(m)/b(m)
            do n = m - 1, 1, -1
              d(n) = (d(n) - cc(n)*d(n+1))/b(n)
            end do
            do q = 1, nk
              ib = si_cb(q, c)
              do i = 1, nb(1)
                n = (q - 1)*nb(1) + i
                State_GVB(i, j, k, v, ib) = d(n)
                if (heat .and. (n <= k_join .or. .not. use_atmosphere)) de(i, j, k, ib) = de(i, j, k, ib) + d(n)**2
              end do
            end do
          end do
        end do
      end do
    end do
    if (heat) then
      ! s1 += -dE_B/(rho T), dE_B = d(B^2)/(8 pi)
      do jb = 1, size(iBlocks)
        ib = iBlocks(jb)
        if (use_atmosphere .and. Key_IB(KI+1, ib)*nb(1) >= k_join) cycle
        call cz_bind(ib)
        ! T of the cells (the column EOS evaluator, ln T only: the values derived gives)
        n = merge(min(nb(1), k_join - r_off), nb(1), use_atmosphere)
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, n
              rpw(i) = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
              sw(i) = s0(i) + State_GVB(i, j, k, s1_, ib)
            end do
            call eos_field_n(tab, F_LNT, n, rpw(1:n), sw(1:n), ltw(1:n))
            call vexp(n, ltw(1:n), tw(1:n))
            do i = 1, n
              State_GVB(i, j, k, s1_, ib) = State_GVB(i, j, k, s1_, ib) &
                  - de(i, j, k, ib)/(2.0_dp*four_pi*(rho0(i) + State_GVB(i, j, k, rho1_, ib))*tw(i))
            end do
          end do
        end do
      end do
      call cz_unbind()
    end if
  end subroutine b_radial_implicit

  ! ================================================================ time step

  ! dt < cfl / sum_d (|v_d| + cs/xi)/l_d, l_r = CZTR's dr/dxi
  real(dp) function cz_stable_dt(self, iBlock) result(dt)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iBlock
    real(dp) :: cs, cs0, crad, ratew, rate, rr, sth, rc, rp(nb(1)), s(nb(1)), q(nb(1), 2), e(nb(1))
    real(dp) :: cs0v(nb(1)), csv(nb(1)), cradv(nb(1)), rcv(nb(1)), insi(nb(1)), dxm, vm, crb
    integer :: i, j, k
    call tic(T_DT)
    call cz_bind(iBlock)
    rate = 0.0_dp
    ratew = 0.0_dp
    insi = 0.0_dp
    do i = 1, nb(1)
      if (si_radial .and. r_off + i >= si_i0 + 2 .and. r_off + i <= si_i1 - 2) insi(i) = 1.0_dp
    end do
    do k = 1, nb(3)
      do j = 1, nb(2)
        sth = sn_t(j)
        do i = 1, nb(1)
          rp(i) = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, iBlock)
          s(i) = s0(i) + State_GVB(i, j, k, s1_, iBlock)
        end do
        call eos_fields_n(tab, [F_LNP, F_GAMMA1], nb(1), rp, s, q)
        call vexp(nb(1), q(:, 1), e)
        ! the speeds and the rate of each cell (a loop that vectorizes: va_cap
        ! is va_eff(cell_cap(...)) from scalars)
        ! si_radial: in the implicit range (not its two end cells) the radial
        ! wave speed leaves the CFL (with a field the radial Alfven speed stays,
        ! unless si_alfven); wave_cfl_max: a looser radial wave CFL there
        ! (wave_cfl_cs_only: with c_s alone)
!GCC$ ivdep
!DIR$ IVDEP
        do i = 1, nb(1)
          cs0v(i) = sqrt(q(i, 2)*e(i)/rp(i))/xi(i)
          csv(i) = cs0v(i)
          if (use_B) then
            dxm = min(drc_r(i), r_c(i)*Metric_DB(2, iBlock)%Dx_I(j), r_c(i)*sn_t(j)*Metric_DB(3, iBlock)%Dx_I(k))
            vm = norm2(State_GVB(i, j, k, vr_:vp_, iBlock))
            csv(i) = cs0v(i) + va_cap(sum(State_GVB(i, j, k, br_:bp_, iBlock)**2), rp(i), cs0v(i), vm, dxm, sw_r(i))
            crb = 0.0_dp
            if (.not. si_alfven) crb = va_cap(State_GVB(i, j, k, br_, iBlock)**2, rp(i), cs0v(i), vm, dxm, sw_r(i))
            cradv(i) = merge(crb, csv(i), insi(i) > 0.0_dp)
          else
            cradv(i) = merge(0.0_dp, csv(i), insi(i) > 0.0_dp)
          end if
          rr = r_c(i)
          rcv(i) = (abs(State_GVB(i, j, k, vr_, iBlock)) + cradv(i))/drdxi_c(i) &
              + (abs(State_GVB(i, j, k, vt_, iBlock)) + csv(i))/(rr*Metric_DB(2, iBlock)%Dx_I(j)) &
              + (abs(State_GVB(i, j, k, vp_, iBlock)) + csv(i))/(rr*sth*Metric_DB(3, iBlock)%Dx_I(k))
        end do
        ! the maxima and the dt report's minima
        do i = 1, nb(1)
          cs = csv(i);  cs0 = cs0v(i);  crad = cradv(i);  rc = rcv(i);  rr = r_c(i)
          if (insi(i) > 0.0_dp .and. wave_cfl_max > 0.0_dp) ratew = max(ratew, merge(cs0, cs, wave_cfl_cs_only)/drdxi_c(i))
          if (cfl/rc < dtl(1, 1)) then
            block
              real(dp) :: t3(3)
              t3 = [(abs(State_GVB(i, j, k, vr_, iBlock)) + crad)/drdxi_c(i), &
                  (abs(State_GVB(i, j, k, vt_, iBlock)) + cs)/(rr*Metric_DB(2, iBlock)%Dx_I(j)), &
                  (abs(State_GVB(i, j, k, vp_, iBlock)) + cs)/(rr*sth*Metric_DB(3, iBlock)%Dx_I(k))]
              dtl(1, :) = [cfl/rc, real(r_off + i, dp), z_c(i)/1.0e5_dp, abs(State_GVB(i, j, k, vr_, iBlock)), crad, cs, &
                  real(maxloc(t3, 1), dp), drdxi_c(i)/1.0e5_dp]
            end block
          end if
          if (insi(i) > 0.0_dp .and. wave_cfl_max > 0.0_dp) then
            if (wave_cfl_max*drdxi_c(i)/merge(cs0, cs, wave_cfl_cs_only) < dtl(2, 1)) &
                dtl(2, :) = [wave_cfl_max*drdxi_c(i)/merge(cs0, cs, wave_cfl_cs_only), real(r_off + i, dp), &
                z_c(i)/1.0e5_dp, abs(State_GVB(i, j, k, vr_, iBlock)), merge(cs0, cs, wave_cfl_cs_only), cs, 2.0_dp, &
                drdxi_c(i)/1.0e5_dp]
          end if
          rate = max(rate, rc)
          rate_r(i) = max(rate_r(i), rc)
        end do
      end do
    end do
    dt = min(cfl/rate, dt_max)
    if (ratew > 0.0_dp) dt = min(dt, wave_cfl_max/ratew)
    if (use_atmosphere) then
      block
        real(dp) :: da
        integer :: ia
        da = atm_dt(iBlock, rate_r(1:nb(1)), ia)
        if (da < dtl(3, 1) .and. ia > 0) dtl(3, :) = [da, real(r_off + ia, dp), z_c(ia)/1.0e5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            3.0_dp, drdxi_c(ia)/1.0e5_dp]
        dt = min(dt, da)
      end block
    end if
    call cz_unbind()
    call toc(T_DT)
  end function cz_stable_dt

  subroutine cz_advance(self, iBlocks, t, dt)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iBlocks(:)
    real(dp), intent(in) :: t, dt
    integer :: stage, jb, ib, i, j, k, v
    logical :: overlap, atm_ovl
    real(dp) :: cdt, wdamp, fdamp
    real(dp), parameter :: c_stage(4) = [0.0_dp, 0.25_dp, 1.0_dp/3.0_dp, 0.5_dp]
    call tic(T_ADV)
    dt_run = dt_run + dt
    call local_steps(dt)
    dt_now = dt
    if (test_freeze_b .and. use_B .and. .not. allocated(fzb_B)) then
      allocate(fzb_B(nb(1), nb(2), nb(3), br_:psi_, size(State_GVB, 5)))
      do jb = 1, size(iBlocks)
        fzb_B(:, :, :, :, iBlocks(jb)) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:psi_, iBlocks(jb))
      end do
    end if
    if (test_b_diff_only) then                         ! (a test: the diffusion's B fluxes only)
      call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.not. (use_B .and. b_diff_form == 'resistive'), time=t + dt)
      call join_bc()
      if (allocated(EtaB_GB)) call b_eta_fill(iBlocks)
      do jb = 1, size(iBlocks)
        call cz_bind(iBlocks(jb))
        call diffusion(iBlocks(jb), dt)
      end do
      call cz_unbind()
      if (allocated(EtaB_GB) .and. .not. b_eta_explicit) call b_radial_implicit(iBlocks, dt)
      call toc(T_ADV)
      return
    end if
    if (freeze_side /= 'none') then
      do jb = 1, size(iBlocks)
        ib = iBlocks(jb)
        if (frz_have(ib)) cycle
        if (.not. block_frozen(ib)) cycle
        if (freeze_cz_mode == 'line_tied') then       ! a static CZ: v = 0 in its cells, from now on
          call cz_bind(ib)
          do i = 1, min(nb(1), k_join - r_off)
            State_GVB(i, 1:nb(2), 1:nb(3), vr_:vp_, ib) = 0.0_dp
          end do
        end if
        frz_B(:, :, :, :, ib) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), 1:nv, ib)
        frz_have(ib) = .true.
      end do
      call cz_unbind()
    end if
    if (lz_correct_steps > 0 .or. lz_correct_dt > 0.0_dp) call lz_impulse(dt)
    ! Coriolis split: half a rotation first (also of StateOld, the RK stages' base)
    if (omega_rot /= 0.0_dp .and. coriolis == 'split') then
      do jb = 1, size(iBlocks)
        call cz_bind(iBlocks(jb))
        call coriolis_rotate(iBlocks(jb), 0.5_dp*dt, .true.)
      end do
      call cz_unbind()
    end if
    overlap = overlap_fills .and. mr_levels == 0     ! (multirate: begin interpolates sources in time)
    if (overlap) then
      if (.not. allocated(p1_IB)) allocate(p1_IB(nb(1), nb(2), nb(3), size(State_GVB, 5)))
      if (size(p1_IB, 4) /= size(State_GVB, 5)) then
        deallocate(p1_IB)
        allocate(p1_IB(nb(1), nb(2), nb(3), size(State_GVB, 5)))
      end if
    end if
    do stage = 1, 4
      call tic(T_GHOST)
      call tic(T_G_RK)
      if (overlap) then
        call arbor_fill_ghosts_begin(State_GVB, n_ghost=2, faces_only=.true., time=t + c_stage(stage)*dt)
      else
        call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.true., time=t + c_stage(stage)*dt)
      end if
      call toc(T_G_RK)
      call toc(T_GHOST)
      if (overlap) then
        ! while the ghost messages travel: the interior cells' equation of
        ! state (local, no ghosts), stored per block, and at stage 1 the
        ! semi-implicit coefficients (interior cells only)
        call tic(T_R_OVL)
        do jb = 1, size(iBlocks)
          ib = iBlocks(jb)
          if (skip_frozen(ib)) cycle
          call cz_bind(ib)
          call tic(T_R_SIC)
          if (si_radial .and. stage == 1) call si_coeff(ib, dt, C_INTERIOR)   ! (its derived: the interior's p1 too)
          call toc(T_R_SIC)
          if (.not. (si_radial .and. stage == 1 .and. si_in(ib))) call derived(ib, D_P1, C_INTERIOR)
          p1_IB(:, :, :, ib) = p1_G(1:nb(1), 1:nb(2), 1:nb(3))
          call mpi_progress()
        end do
        call cz_unbind()
        call toc(T_R_OVL)
        call tic(T_GHOST)
        call tic(T_G_RK)
        call arbor_fill_ghosts_end(State_GVB)
        call toc(T_G_RK)
        call toc(T_GHOST)
      end if
      call tic(T_GHOST)
      call tic(T_G_JOIN)
      call join_bc()
      call toc(T_G_JOIN)
      call toc(T_GHOST)
      call tic(T_RHS)
      do jb = 1, size(iBlocks)
        ib = iBlocks(jb)
        if (skip_frozen(ib)) cycle                    ! (the frozen side: not computed, restored below)
        call cz_bind(ib)
        want_work = budget_on .and. stage == 1
        if (overlap) then
          ! the ghost shell's p1 now, the interior's from the wait
          call tic(T_R_EOS)
          call derived(ib, D_P1, C_SHELL)
          p1_G(1:nb(1), 1:nb(2), 1:nb(3)) = p1_IB(:, :, :, ib)
          eos_fresh = .true.
          call toc(T_R_EOS)
        else
          call tic(T_R_SIC)
          if (si_radial .and. stage == 1) call si_coeff(ib, dt, C_ALL)
          call toc(T_R_SIC)
        end if
        call rhs(ib)
        call tic(T_R_UPD)
        if (allocated(fzb_B)) R_CV(:, :, :, br_:psi_) = 0.0_dp     ! (test_freeze_b)
        want_work = .false.
        ! the RK update, with each shell's step factor tsc (local time stepping):
        ! one pass in memory order (the same arithmetic as scaling R_CV first)
        cdt = dt/real(5 - stage, dp)
        if (si_in(ib)) then
          ! the solved velocities (v_r; with si_alfven v_theta, v_phi) and rho1
          ! (the solves' operator uses the stage's starting density) after the
          ! column solves below, the other variables now (as outside the range)
          do v = 1, nv
            if (si_var(v) .or. v == rho1_) then
              siR_B(:, :, :, v, ib) = R_CV(:, :, :, v)
            else
              do k = 1, nb(3)
                do j = 1, nb(2)
                  do i = 1, nb(1)
                    State_GVB(i, j, k, v, ib) = StateOld_GVB(i, j, k, v, ib) + cdt*R_CV(i, j, k, v)
                  end do
                end do
              end do
            end if
          end do
        else if (dt1_now > 0.0_dp) then
          do v = 1, nv
            do k = 1, nb(3)
              do j = 1, nb(2)
                do i = 1, nb(1)
                  State_GVB(i, j, k, v, ib) = StateOld_GVB(i, j, k, v, ib) + cdt*(tsc(i)*R_CV(i, j, k, v))
                end do
              end do
            end do
          end do
        else
          do v = 1, nv
            do k = 1, nb(3)
              do j = 1, nb(2)
                do i = 1, nb(1)
                  State_GVB(i, j, k, v, ib) = StateOld_GVB(i, j, k, v, ib) + cdt*R_CV(i, j, k, v)
                end do
              end do
            end do
          end do
        end if
        call toc(T_R_UPD)
      end do
      call cz_unbind()
      call toc(T_RHS)
      ! si_radial: the v_r tendencies of the implicit range through the column
      ! operator (si_alfven: v_theta, v_phi through theirs), then those blocks' RK update
      call tic(T_SI)
      if (si_radial) then
        call tic(T_S_SOL)
        call si_column_solve(siR_B(:, :, :, vr_, :), siS_B)
        if (si_alfven) then
          call si_column_solve(siR_B(:, :, :, vt_, :), siSa_B, siR_B(:, :, :, vp_, :))   ! (the same matrix)
        end if
        call toc(T_S_SOL)
        cdt = dt/real(5 - stage, dp)
        do jb = 1, size(iBlocks)
          ib = iBlocks(jb)
          if (.not. si_in(ib)) cycle
          do v = 1, nv
            if (.not. (si_var(v) .or. v == rho1_)) cycle
            do k = 1, nb(3)
              do j = 1, nb(2)
                do i = 1, nb(1)
                  State_GVB(i, j, k, v, ib) = StateOld_GVB(i, j, k, v, ib) + cdt*siR_B(i, j, k, v, ib)
                end do
              end do
            end do
          end do
        end do
      end if
      call toc(T_SI)
    end do
    if (budget_on) then
      do jb = 1, size(iBlocks)
        call budget_hydro(iBlocks(jb), dt)
      end do
    end if
    call tic(T_TRANS)
    call cz_transport_all(iBlocks, dt, .false., t + dt)
    call toc(T_TRANS)
    if (use_atmosphere) then
      ! conduction (the ghosts in every direction), losses, heating
      call tic(T_GHOST)
      call tic(T_G_ATM)
      ! (not psi: the atmosphere physics does not read it, and the diffusion's
      ! fill refreshes every variable before anything else does)
      atm_ovl = overlap .and. si_cols .and. (use_losses .or. use_heating)
      if (atm_ovl) then
        call arbor_fill_ghosts_begin(State_GVB, iVarMax=merge(psi_ - 1, nv, use_B), n_ghost=2, faces_only=.true., &
            time=t + dt)
        call toc(T_G_ATM)
        call toc(T_GHOST)
        ! while the messages travel: the losses and heating (per cell, the
        ! state before the conduction; the radial ghosts they read are this
        ! rank's, filled by begin)
        call tic(T_A_OVL)
        do jb = 1, size(iBlocks)
          ib = iBlocks(jb)
          if (skip_frozen(ib)) cycle
          call cz_bind(ib)
          if (.not. atm_active(ib)) cycle
          call states_interior(ib)
          call atm_losses(ib, at_top)
          call mpi_progress()
        end do
        call cz_unbind()
        call toc(T_A_OVL)
        call tic(T_GHOST)
        call tic(T_G_ATM)
        call arbor_fill_ghosts_end(State_GVB)
      else
        call arbor_fill_ghosts(State_GVB, iVarMax=merge(psi_ - 1, nv, use_B), n_ghost=2, faces_only=.true., time=t + dt)
      end if
      call toc(T_G_ATM)
      call tic(T_G_JOIN)
      call join_bc()
      call toc(T_G_JOIN)
      call toc(T_GHOST)
      call tic(T_ATM)
      do jb = 1, size(iBlocks)
        if (skip_frozen(iBlocks(jb))) cycle
        call cz_bind(iBlocks(jb))
        call atm_block(iBlocks(jb), dt, at_bottom, at_top, shell_steps(), losses_ready=atm_ovl)
      end do
      call cz_unbind()
      call toc(T_ATM)
    end if
    if (use_diffusion) then
      ! (b_diff_form = 'resistive': the edge ghosts too, for the tangential derivatives at block faces)
      call tic(T_GHOST)
      call tic(T_G_DIF)
      call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.not. (use_B .and. b_diff_form == 'resistive'), time=t + dt)
      call toc(T_G_DIF)
      call tic(T_G_JOIN)
      call join_bc()
      call toc(T_G_JOIN)
      call toc(T_GHOST)
      call tic(T_DIFF)
      call tic(T_D_SIP)
      if (si_radial .and. si_diff_precond) then             ! (v_r before the diffusion)
        do jb = 1, size(iBlocks)
          if (si_in(iBlocks(jb))) siV_B(:, :, :, iBlocks(jb)) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), vr_, iBlocks(jb))
        end do
      end if
      call toc(T_D_SIP)
      call tic(T_D_ETA)
      if (allocated(EtaB_GB)) call b_eta_fill(iBlocks)       ! (b_diff_form = 'resistive')
      call toc(T_D_ETA)
      do jb = 1, size(iBlocks)
        if (skip_frozen(iBlocks(jb))) cycle
        call cz_bind(iBlocks(jb))
        call diffusion(iBlocks(jb), dt)
      end do
      call cz_unbind()
      call tic(T_D_IMP)
      if (allocated(EtaB_GB) .and. .not. b_eta_explicit) call b_radial_implicit(iBlocks, dt)
      call toc(T_D_IMP)
      ! si_radial: the diffusion's v_r change through the same column operator
      call tic(T_D_SIP)
      if (si_radial .and. si_diff_precond) then
        do jb = 1, size(iBlocks)
          ib = iBlocks(jb)
          if (si_in(ib)) siX_B(:, :, :, ib) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), vr_, ib) - siV_B(:, :, :, ib)
        end do
        call si_column_solve(siX_B, siS_B)
        do jb = 1, size(iBlocks)
          ib = iBlocks(jb)
          if (si_in(ib)) State_GVB(1:nb(1), 1:nb(2), 1:nb(3), vr_, ib) = siV_B(:, :, :, ib) + siX_B(:, :, :, ib)
        end do
      end if
      call toc(T_D_SIP)
      call toc(T_DIFF)
    end if
    ! the atmosphere's floors (rho_floor; T_floor where atm_switch > 0.5), once
    ! per step after all the operators (the floors report counts their action)
    if (use_atmosphere) then
      do jb = 1, size(iBlocks)
        if (skip_frozen(iBlocks(jb))) cycle
        call cz_bind(iBlocks(jb))
        call atm_floors(iBlocks(jb))
      end do
      call cz_unbind()
    end if
    ! freeze_cz (modes line_tied, all): the CZ cells back to their stored state
    if (freeze_side == 'cz' .and. freeze_cz_mode /= 'ghost') then
      do jb = 1, size(iBlocks)
        call cz_bind(iBlocks(jb))
        do i = 1, min(nb(1), k_join - r_off)
          if (freeze_cz_mode == 'line_tied' .and. use_B .and. r_off + i > k_join - freeze_float_cells) then
            ! (the top cells: B_theta, B_phi and psi float)
            State_GVB(i, 1:nb(2), 1:nb(3), 1:br_, iBlocks(jb)) = frz_B(i, :, :, 1:br_, iBlocks(jb))
          else
            State_GVB(i, 1:nb(2), 1:nb(3), 1:nv, iBlocks(jb)) = frz_B(i, :, :, :, iBlocks(jb))
          end if
        end do
      end do
      call cz_unbind()
    end if
    ! friction on v above the join (v_damp_n > 0), for the run's first steps (n_adv: local_steps)
    if (v_damp_n > 0.0_dp .and. use_atmosphere) then
      if (t_damp0 < 0.0_dp) t_damp0 = t                ! (the run's start: the friction schedule's zero)
      if (v_damp_time_s > 0.0_dp) then                 ! the schedule in model time
        if (t + dt - t_damp0 <= v_damp_time_s) then
          cdt = 1.0_dp
        else if (t + dt - t_damp0 <= v_damp_time_s + v_damp_ramp_s) then
          cdt = 1.0_dp - (t + dt - t_damp0 - v_damp_time_s)/max(v_damp_ramp_s, tiny(1.0_dp))
        else
          cdt = 0.0_dp
        end if
      else if (n_adv <= v_damp_steps) then
        cdt = 1.0_dp
      else if (n_adv <= v_damp_steps + v_damp_ramp) then
        cdt = 1.0_dp - real(n_adv - v_damp_steps, dp)/real(v_damp_ramp + 1, dp)
      else
        cdt = 0.0_dp
      end if
      if (cdt > 0.0_dp) then
        cdt = 1.0_dp - cdt/(v_damp_n + 1.0_dp)
        do jb = 1, size(iBlocks)
          if (skip_frozen(iBlocks(jb))) cycle
          call cz_bind(iBlocks(jb))
          do i = 1, nb(1)
            wdamp = damp_weight(z_c(i))
            if (wdamp <= 0.0_dp) cycle
            fdamp = 1.0_dp - wdamp*(1.0_dp - cdt)     ! (cdt: the full factor, 1 - w(t)/(v_damp_n + 1))
            ! (v_r towards the init file's profile: the wind keeps flowing; 0 in the CZ)
            State_GVB(i, 1:nb(2), 1:nb(3), vr_, iBlocks(jb)) = v_init(i) &
                + fdamp*(State_GVB(i, 1:nb(2), 1:nb(3), vr_, iBlocks(jb)) - v_init(i))
            State_GVB(i, 1:nb(2), 1:nb(3), vt_:vp_, iBlocks(jb)) = fdamp*State_GVB(i, 1:nb(2), 1:nb(3), vt_:vp_, iBlocks(jb))
          end do
        end do
        call cz_unbind()
      end if
    end if
    ! a cell-centred wall keeps v_r = 0 in its cell (as the next fill would),
    ! and the field components odd about it (B_r for a conductor, B_h for a radial field)
    do jb = 1, size(iBlocks)
      call cz_bind(iBlocks(jb))
      if (wall_bottom == 'cell' .and. at_bottom) State_GVB(1, :, :, vr_, iBlocks(jb)) = 0.0_dp
      if (wall_top == 'cell' .and. at_top) State_GVB(nb(1), :, :, vr_, iBlocks(jb)) = 0.0_dp
      if (use_B) then
        call glm_damp(iBlocks(jb), dt)
        if (wall_bottom == 'cell' .and. at_bottom) call zero_odd_b(1, b_bc_bottom)
        if (wall_top == 'cell' .and. at_top) call zero_odd_b(nb(1), b_bc_top)
      end if
    end do
    call cz_unbind()
    ! ... and the second half after the step
    if (omega_rot /= 0.0_dp .and. coriolis == 'split') then
      do jb = 1, size(iBlocks)
        call cz_bind(iBlocks(jb))
        call coriolis_rotate(iBlocks(jb), 0.5_dp*dt, .false.)
      end do
      call cz_unbind()
    end if
    ! ghost mode: the frozen side's blocks exactly as stored
    if (freeze_side /= 'none' .and. freeze_cz_mode == 'ghost') then
      do jb = 1, size(iBlocks)
        ib = iBlocks(jb)
        if (block_frozen(ib)) State_GVB(1:nb(1), 1:nb(2), 1:nb(3), 1:nv, ib) = frz_B(:, :, :, :, ib)
      end do
    end if
    if (allocated(fzb_B)) then                        ! test_freeze_b: the field back
      do jb = 1, size(iBlocks)
        State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:psi_, iBlocks(jb)) = fzb_B(:, :, :, :, iBlocks(jb))
      end do
    end if
    dt_last = dt
    if (UseMultirate) then
      tloc = max(tloc, t + dt)                          ! (every level ends at the level-1 step's end)
    else
      tloc = tloc + tsc_g*dt
    end if
    call toc(T_ADV)
  contains
    ! the friction's weight at height z [cm] (see v_damp_below_km)
    real(dp) function damp_weight(z) result(w)
      real(dp), intent(in) :: z
      real(dp) :: zb, za, tw
      zb = -v_damp_below_km*1.0e5_dp
      za = v_damp_above_km*1.0e5_dp
      tw = max(v_damp_taper_km, 1.0e-6_dp)*1.0e5_dp
      w = 1.0_dp
      if (z < 0.0_dp .or. (z <= 0.0_dp .and. v_damp_below_km <= 0.0_dp)) then
        if (v_damp_below_km <= 0.0_dp) then
          w = 0.0_dp
        else if (z < zb) then
          w = cos(0.5_dp*pi*min(1.0_dp, (zb - z)/tw))**2
        end if
      else if (v_damp_above_km >= 0.0_dp .and. z > za) then
        w = cos(0.5_dp*pi*min(1.0_dp, (z - za)/tw))**2
      end if
    end function damp_weight
    subroutine zero_odd_b(i, bc)
      integer, intent(in) :: i
      character(*), intent(in) :: bc
      if (bc == 'conductor') then
        State_GVB(i, :, :, br_, iBlocks(jb)) = 0.0_dp
      else
        State_GVB(i, :, :, bt_:bp_, iBlocks(jb)) = 0.0_dp
      end if
    end subroutine zero_odd_b
  end subroutine cz_advance

  ! Local time stepping: tsc(i) = max(1, min(dt_i, dt1)/dt) with dt_i the
  ! stable step of shell i (cfl over its largest rate on all ranks), then the
  ! smallest tsc within +-2 shells (a smooth step that stays stable). Uses
  ! and resets rate_r, filled by this step's stable_dt calls. dt1_ramp_steps
  ! > 0: dt1 scaled by max(0, 1 - (n - dt1_hold_steps)/dt1_ramp_steps), n this
  ! run's steps (dt1 held for the first dt1_hold_steps).
  subroutine local_steps(dt)
    real(dp), intent(in) :: dt
    real(dp) :: a(1:n_r)
    integer :: i
    dt1_now = dt1
    if (dt1_ramp_steps > 0) dt1_now = dt1*min(1.0_dp, max(0.0_dp, &
        1.0_dp - real(n_adv - dt1_hold_steps, dp)/real(dt1_ramp_steps, dp)))
    n_adv = n_adv + 1
    if (dt1_now <= 0.0_dp) then
      tsc_g = 1.0_dp
      rate_r_g = 0.0_dp
      return
    end if
    a = rate_r_g
    call MPI_Allreduce(MPI_IN_PLACE, a, n_r, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD)
    do i = 1, n_r
      a(i) = max(1.0_dp, min(cfl/max(a(i), tiny(1.0_dp)), dt1_now)/dt)
    end do
    do i = 1, n_r
      tsc_g(i) = minval(a(max(1, i-2):min(n_r, i+2)))
    end do
    rate_r_g = 0.0_dp
  end subroutine local_steps

  ! ================================================================ hydrodynamics

  ! Radial derivative (CZTR's stencil and metric) of a one-block array at
  ! the interior cells
  subroutine ddr_block(F_G, Deriv_C)
    real(dp), intent(in) :: F_G(1-ng:, 1-ng:, 1-ng:)
    real(dp), intent(out) :: Deriv_C(:, :, :)
    integer :: i, j, k
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, nb(1)
          Deriv_C(i, j, k) = (F_G(i-2, j, k) - 8.0_dp*F_G(i-1, j, k) + 8.0_dp*F_G(i+1, j, k) - F_G(i+2, j, k)) &
              /(12.0_dp*drdxi_c(i))
        end do
      end do
    end do
  end subroutine ddr_block

  ! Right-hand side R_CV of block ib (ghost cells filled)
  subroutine rhs(ib)
    integer, intent(in) :: ib
    integer :: i, j, k, c
    real(dp) :: rr, sth, cth, rho, vr, vt, vp
    call tic(T_R_EOS)
    if (.not. eos_fresh) call derived(ib, D_P1, C_ALL)
    eos_fresh = .false.
    call toc(T_R_EOS)
    call tic(T_R_MASS)
    associate (th => Metric_DB(2, ib)%X_I)
      ! ---- mass: -div(rho v)
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          do i = 1 - ng, nb(1) + ng
            W_G(i, j, k) = r_c(i)**2*rho_G(i, j, k)*State_GVB(i, j, k, vr_, ib)
          end do
        end do
      end do
      call odd_wall_flux(W_G)
      call ddr_block(W_G, Dr_C)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            R_CV(i, j, k, rho1_) = -Dr_C(i, j, k)/r_c(i)**2
          end do
        end do
      end do
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          sth = sin(th(j))
          do i = 1 - ng, nb(1) + ng
            W_G(i, j, k) = sth*rho_G(i, j, k)*State_GVB(i, j, k, vt_, ib)
          end do
        end do
      end do
      call arbor_deriv_cell(W_G, 2, ib, D_C)
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          do i = 1 - ng, nb(1) + ng
            W_G(i, j, k) = rho_G(i, j, k)*State_GVB(i, j, k, vp_, ib)
          end do
        end do
      end do
      call arbor_deriv_cell(W_G, 3, ib, Dr_C)
      do k = 1, nb(3)
        do j = 1, nb(2)
          sth = sin(th(j))
          do i = 1, nb(1)
            R_CV(i, j, k, rho1_) = R_CV(i, j, k, rho1_) - (D_C(i, j, k) + Dr_C(i, j, k))/(r_c(i)*sth)
            if (rsst_h2012) R_CV(i, j, k, rho1_) = R_CV(i, j, k, rho1_)/xi(i)**2
          end do
        end do
      end do

      call toc(T_R_MASS)
      ! ---- momentum and entropy: advection -(v.grad) q for q = v_r, v_theta, v_phi, s1
      call tic(T_R_ADV)
      do c = vr_, s1_
        call ddr_block(State_GVB(:, :, :, c, ib), Dr_C)
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, nb(1)
              R_CV(i, j, k, c) = -State_GVB(i, j, k, vr_, ib)*Dr_C(i, j, k)
            end do
          end do
        end do
        call arbor_deriv_cell(State_GVB, c, 2, ib, D_C)
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, nb(1)
              R_CV(i, j, k, c) = R_CV(i, j, k, c) - State_GVB(i, j, k, vt_, ib)*D_C(i, j, k)/r_c(i)
            end do
          end do
        end do
        call arbor_deriv_cell(State_GVB, c, 3, ib, D_C)
        do k = 1, nb(3)
          do j = 1, nb(2)
            sth = sin(th(j))
            do i = 1, nb(1)
              R_CV(i, j, k, c) = R_CV(i, j, k, c) - State_GVB(i, j, k, vp_, ib)*D_C(i, j, k)/(r_c(i)*sth)
            end do
          end do
        end do
      end do

      call toc(T_R_ADV)
      ! ---- pressure gradient, gravity on rho1, curvature terms, ds0/dr
      call tic(T_R_PRS)
      call ddr_block(p1_G, Dr_C)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            rho = rho_G(i, j, k)
            R_CV(i, j, k, vr_) = R_CV(i, j, k, vr_) - (Dr_C(i, j, k) + g_c(i)*rq(i)*State_GVB(i, j, k, rho1_, ib))/rho
            if (use_atmosphere) R_CV(i, j, k, vr_) = R_CV(i, j, k, vr_) - S_bgf(r_off + i)*bgforce(r_off + i)/rho
            if (want_work) work_CB(i, j, k, ib) = -State_GVB(i, j, k, vr_, ib) &
                *(Dr_C(i, j, k) + g_c(i)*rq(i)*State_GVB(i, j, k, rho1_, ib))
          end do
        end do
      end do
      call arbor_deriv_cell(p1_G, 2, ib, D_C)
      call arbor_deriv_cell(p1_G, 3, ib, Dr_C)
      do k = 1, nb(3)
        do j = 1, nb(2)
          sth = sin(th(j))
          cth = cos(th(j))
          do i = 1, nb(1)
            rho = rho_G(i, j, k)
            rr = r_c(i)
            vr = State_GVB(i, j, k, vr_, ib)
            vt = State_GVB(i, j, k, vt_, ib)
            vp = State_GVB(i, j, k, vp_, ib)
            R_CV(i, j, k, vr_) = R_CV(i, j, k, vr_) + (vt**2 + vp**2)/rr
            if (want_work) then
              workh_CB(i, j, k, ib) = -vt*D_C(i, j, k)/rr - vp*Dr_C(i, j, k)/(rr*sth)
              work_CB(i, j, k, ib) = work_CB(i, j, k, ib) + workh_CB(i, j, k, ib)
            end if
            R_CV(i, j, k, vt_) = R_CV(i, j, k, vt_) - D_C(i, j, k)/(rr*rho) - (vr*vt - vp**2*cth/sth)/rr
            R_CV(i, j, k, vp_) = R_CV(i, j, k, vp_) - Dr_C(i, j, k)/(rr*sth*rho) - (vr*vp + vt*vp*cth/sth)/rr
            R_CV(i, j, k, s1_) = R_CV(i, j, k, s1_) - vr*ds0dr(i)
          end do
        end do
      end do
    end associate
    call toc(T_R_PRS)
    if (use_B) call rhs_mhd(ib)
    if (omega_rot /= 0.0_dp .and. coriolis == 'rk') call rhs_coriolis(ib)
  end subroutine rhs

  ! The exact solution of dv/dt = -2 Omega x v over tau (each shell's own
  ! tau tsc(i)): v rotated about the axis by phi = -2 omega_rot tau (Rodrigues,
  ! the axis unit vector in the cell's spherical components as rhs_coriolis):
  !   v' = v_par + cos(phi) v_perp + sin(phi) (n x v)
  ! with_old: StateOld too (before the RK stages, which build on it)
  subroutine coriolis_rotate(ib, tau, with_old)
    integer, intent(in) :: ib
    real(dp), intent(in) :: tau
    logical, intent(in) :: with_old
    integer :: i, j, k
    real(dp) :: e(3,3), n(3), v(3), vpar(3), phi
    real(dp), parameter :: zg(3) = [0.0_dp, 0.0_dp, 1.0_dp]
    do k = 1, nb(3)
      do j = 1, nb(2)
        e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
        if (Key_IB(KP, ib) == 2) then
          n = matmul(transpose(e), yy_cart(zg))
        else
          n = matmul(transpose(e), zg)
        end if
        do i = 1, nb(1)
          phi = -2.0_dp*omega_rot*tau*tsc(i)
          v = State_GVB(i, j, k, vr_:vp_, ib)
          vpar = dot_product(n, v)*n
          State_GVB(i, j, k, vr_:vp_, ib) = vpar + cos(phi)*(v - vpar) &
              + sin(phi)*[n(2)*v(3) - n(3)*v(2), n(3)*v(1) - n(1)*v(3), n(1)*v(2) - n(2)*v(1)]
          if (with_old) StateOld_GVB(i, j, k, vr_:vp_, ib) = State_GVB(i, j, k, vr_:vp_, ib)
        end do
      end do
    end do
  end subroutine coriolis_rotate

  ! Coriolis acceleration -2 Omega x v, Omega = omega_rot along the global z
  ! axis: in the Yin frame z itself, in the Yang frame its local y (yy_cart
  ! maps (0,0,1) to (0,1,0)); its spherical components from sph_basis, the
  ! cross product in the orthonormal (r, theta, phi) basis
  subroutine rhs_coriolis(ib)
    integer, intent(in) :: ib
    integer :: i, j, k
    real(dp) :: e(3,3), om(3), v(3)
    real(dp), parameter :: zg(3) = [0.0_dp, 0.0_dp, 1.0_dp]
    do k = 1, nb(3)
      do j = 1, nb(2)
        e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
        if (Key_IB(KP, ib) == 2) then
          om = omega_rot*matmul(transpose(e), yy_cart(zg))
        else
          om = omega_rot*matmul(transpose(e), zg)
        end if
        do i = 1, nb(1)
          v = State_GVB(i, j, k, vr_:vp_, ib)
          R_CV(i, j, k, vr_:vp_) = R_CV(i, j, k, vr_:vp_) &
              - 2.0_dp*[om(2)*v(3) - om(3)*v(2), om(3)*v(1) - om(1)*v(3), om(1)*v(2) - om(2)*v(1)]
        end do
      end do
    end do
  end subroutine rhs_coriolis

  ! ================================================================ magnetic field

  ! The MHD terms of block ib (ghosts filled; rho_G from derived):
  !   R_v   += f (curl B) x B / (4 pi rho)   (f: the Alfven-speed cap, see lorentz_f)
  !   R_B    = curl(v x B) - grad psi
  !   R_psi  = -c_h^2 div B,  c_h = c_s0/xi + v_A
  ! Curl and divergence in spherical components: radial derivatives with
  ! CZTR's stencil (ddr_block), angular ones with ARBOR's (arbor_deriv_cell).
  subroutine rhs_mhd(ib)
    integer, intent(in) :: ib
    integer :: i, j, k
    real(dp) :: b(3), v(3), jj(3), ch, rho, rp, b2, vm, dxm, dxg, csl, cap, f, s_inv
    logical :: capon, glmc
    integer :: i0, i1
    call tic(T_R_MCURL)
    ! J = curl B (B straight from the state: a contiguous section, no copy)
    call curl_block(State_GVB(:, :, :, br_:bp_, ib), ib, J_CV)
    ! E = v x B where the curl's derivatives read it: the interior and the ng
    ! ghost layers across each face (no edge or corner ghosts: no derivative
    ! reads them; zeroed at the first use, other routines' scratch values
    ! after that), then its curl
    if (.not. b_gv_zeroed) then
      B_GV = 0.0_dp
      b_gv_zeroed = .true.
    end if
    do k = 1 - ng, nb(3) + ng
      do j = 1 - ng, nb(2) + ng
        if ((j < 1 .or. j > nb(2)) .and. (k < 1 .or. k > nb(3))) cycle      ! (an edge column)
        i0 = merge(1 - ng, 1, j >= 1 .and. j <= nb(2) .and. k >= 1 .and. k <= nb(3))
        i1 = merge(nb(1) + ng, nb(1), i0 == 1 - ng)
        do i = i0, i1
          B_GV(i, j, k, 1) = State_GVB(i, j, k, vt_, ib)*State_GVB(i, j, k, bp_, ib) &
              - State_GVB(i, j, k, vp_, ib)*State_GVB(i, j, k, bt_, ib)
          B_GV(i, j, k, 2) = State_GVB(i, j, k, vp_, ib)*State_GVB(i, j, k, br_, ib) &
              - State_GVB(i, j, k, vr_, ib)*State_GVB(i, j, k, bp_, ib)
          B_GV(i, j, k, 3) = State_GVB(i, j, k, vr_, ib)*State_GVB(i, j, k, bt_, ib) &
              - State_GVB(i, j, k, vt_, ib)*State_GVB(i, j, k, br_, ib)
        end do
      end do
    end do
    call curl_block(B_GV, ib, C_CV)
    ! div B and grad psi
    call div_block(State_GVB(:, :, :, br_:bp_, ib), ib, divB_C)
    call ddr_block(State_GVB(:, :, :, psi_, ib), Dr_C)
    R_CV(:, :, :, br_) = C_CV(:, :, :, 1) - Dr_C
    call arbor_deriv_cell(State_GVB, psi_, 2, ib, D_C)
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, nb(1)
          R_CV(i, j, k, bt_) = C_CV(i, j, k, 2) - D_C(i, j, k)/r_c(i)
        end do
      end do
    end do
    call arbor_deriv_cell(State_GVB, psi_, 3, ib, D_C)
    call toc(T_R_MCURL)
    call tic(T_R_MLOOP)
    ! The Lorentz force with its cap factor and the GLM source. cell_cs,
    ! cell_cap, lorentz_f, va_eff and ch_cell written out (the same expressions
    ! in the same order): gfortran does not inline them here, and the calls
    ! kept this loop scalar
    capon = va_dt_target > 0.0_dp
    glmc = glm_cfl .and. dt_now > 0.0_dp
    associate (dth => Metric_DB(2, ib)%Dx_I, dph => Metric_DB(3, ib)%Dx_I)
      do k = 1, nb(3)
        do j = 1, nb(2)
          ! (each iteration writes only its own cell: no dependence between them)
!GCC$ ivdep
!DIR$ IVDEP
          do i = 1, nb(1)
            R_CV(i, j, k, bp_) = C_CV(i, j, k, 3) - D_C(i, j, k)/(r_c(i)*sn_t(j))
            b = State_GVB(i, j, k, br_:bp_, ib)
            jj = J_CV(i, j, k, :)
            rho = rho_G(i, j, k)
            rp = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            b2 = sum(b**2)
            vm = norm2(State_GVB(i, j, k, vr_:vp_, ib))
            dxm = min(drc_r(i), r_c(i)*dth(j), r_c(i)*sn_t(j)*dph(k))
            dxg = min(drdxi_c(i), r_c(i)*dth(j), r_c(i)*sn_t(j)*dph(k))      ! (the GLM speed's cap: the old one)
            ! the force: cap(cell_cs), lorentz_f
            csl = sqrt(5.0_dp/3.0_dp*max(p0(i) + p1_G(i, j, k), 1.0e-30_dp)/rp)/xi(i)
            cap = merge(max(cfl*dxm/va_dt_target - (csl + vm), max(1.0e5_dp, va_floor_cs*csl*sw_r(i))), va_max, capon)
            f = merge(1.0_dp/sqrt(1.0_dp + (b2/(four_pi*rp*cap**2))**2), 1.0_dp, cap > 0.0_dp)
            R_CV(i, j, k, vr_) = R_CV(i, j, k, vr_) + f*(jj(2)*b(3) - jj(3)*b(2))/(four_pi*rho)
            R_CV(i, j, k, vt_) = R_CV(i, j, k, vt_) + f*(jj(3)*b(1) - jj(1)*b(3))/(four_pi*rho)
            R_CV(i, j, k, vp_) = R_CV(i, j, k, vp_) + f*(jj(1)*b(2) - jj(2)*b(1))/(four_pi*rho)
            ! ch_cell: c_s0/xi + va_eff(cap(c_s0/xi)), 'cfl' the step's bound
            csl = cs0(i)/xi(i)
            cap = merge(max(cfl*dxg/va_dt_target - (csl + vm), max(1.0e5_dp, va_floor_cs*csl*sw_r(i))), va_max, capon)
            f = merge(1.0_dp/sqrt(1.0_dp + (b2/(four_pi*rp*cap**2))**2), 1.0_dp, cap > 0.0_dp)
            ch = csl + sqrt(f*b2/(four_pi*rp))
            s_inv = 1.0_dp/drdxi_c(i) + 1.0_dp/(r_c(i)*dth(j)) + 1.0_dp/(r_c(i)*sn_t(j)*dph(k))
            ch = merge(max(ch, glm_cfl_frac*cfl/(s_inv*(tsc(i)*dt_now))), ch, glmc)
            R_CV(i, j, k, psi_) = -ch**2*divB_C(i, j, k)
          end do
        end do
      end do
    end associate
    call toc(T_R_MLOOP)
  end subroutine rhs_mhd

  ! Curl of the spherical-component vector A_GV (cells with ghosts) at the
  ! interior cells of block ib:
  !   C_r     = [d_theta(sin theta A_phi) - d_phi A_theta] / (r sin theta)
  !   C_theta = [d_phi A_r / sin theta - d_r(r A_phi)] / r
  !   C_phi   = [d_r(r A_theta) - d_theta A_r] / r
  subroutine curl_block(A_GV, ib, C)
    real(dp), intent(in) :: A_GV(1-ng:, 1-ng:, 1-ng:, :)
    integer, intent(in) :: ib
    real(dp), intent(out) :: C(:, :, :, :)
    integer :: i, j, k
    ! (the components go to the derivatives directly; the scaled ones through
    ! W_G, filled with i innermost)
    associate (th => Metric_DB(2, ib)%X_I)
      ! C_r
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          W_G(:, j, k) = sin(th(j))*A_GV(:, j, k, 3)
        end do
      end do
      call arbor_deriv_cell(W_G, 2, ib, D_C)
      call arbor_deriv_cell(A_GV(:, :, :, 2), 3, ib, Dr_C)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            C(i, j, k, 1) = (D_C(i, j, k) - Dr_C(i, j, k))/(r_c(i)*sin(th(j)))
          end do
        end do
      end do
      ! C_theta
      call arbor_deriv_cell(A_GV(:, :, :, 1), 3, ib, D_C)
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          do i = 1 - ng, nb(1) + ng
            W_G(i, j, k) = r_c(i)*A_GV(i, j, k, 3)
          end do
        end do
      end do
      call ddr_block(W_G, Dr_C)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            C(i, j, k, 2) = (D_C(i, j, k)/sin(th(j)) - Dr_C(i, j, k))/r_c(i)
          end do
        end do
      end do
      ! C_phi
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          do i = 1 - ng, nb(1) + ng
            W_G(i, j, k) = r_c(i)*A_GV(i, j, k, 2)
          end do
        end do
      end do
      call ddr_block(W_G, Dr_C)
      call arbor_deriv_cell(A_GV(:, :, :, 1), 2, ib, D_C)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            C(i, j, k, 3) = (Dr_C(i, j, k) - D_C(i, j, k))/r_c(i)
          end do
        end do
      end do
    end associate
  end subroutine curl_block

  ! Divergence of the spherical-component vector A_GV at the interior cells:
  !   (1/r^2) d_r(r^2 A_r) + [d_theta(sin theta A_theta) + d_phi A_phi] / (r sin theta)
  subroutine div_block(A_GV, ib, Dv)
    real(dp), intent(in) :: A_GV(1-ng:, 1-ng:, 1-ng:, :)
    integer, intent(in) :: ib
    real(dp), intent(out) :: Dv(:, :, :)
    integer :: i, j, k
    associate (th => Metric_DB(2, ib)%X_I)
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          do i = 1 - ng, nb(1) + ng
            W_G(i, j, k) = r_c(i)**2*A_GV(i, j, k, 1)
          end do
        end do
      end do
      call ddr_block(W_G, Dr_C)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            Dv(i, j, k) = Dr_C(i, j, k)/r_c(i)**2
          end do
        end do
      end do
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          W_G(:, j, k) = sin(th(j))*A_GV(:, j, k, 2)
        end do
      end do
      call arbor_deriv_cell(W_G, 2, ib, D_C)
      call arbor_deriv_cell(A_GV(:, :, :, 3), 3, ib, Dr_C)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            Dv(i, j, k) = Dv(i, j, k) + (D_C(i, j, k) + Dr_C(i, j, k))/(r_c(i)*sin(th(j)))
          end do
        end do
      end do
    end associate
  end subroutine div_block

  ! The GLM cleaning speed of a cell: c_s0/xi + v_A ('fast'; v_A capped, see
  ! va_eff), or the largest its step allows ('cfl', see glm_ch), at least the former
  real(dp) function ch_cell(i, j, k, ib) result(ch)
    integer, intent(in) :: i, j, k, ib
    real(dp) :: s_inv, dti
    ! (the cap with the radial size too, as before va_cap_horizontal: psi
    ! propagates explicitly in r; va_cap is va_eff(cell_cap(...)) from scalars)
    ch = cs0(i)/xi(i) + va_cap(sum(State_GVB(i, j, k, br_:bp_, ib)**2), rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib), &
        cs0(i)/xi(i), norm2(State_GVB(i, j, k, vr_:vp_, ib)), &
        min(drdxi_c(i), r_c(i)*Metric_DB(2, ib)%Dx_I(j), r_c(i)*sn_t(j)*Metric_DB(3, ib)%Dx_I(k)), sw_r(i))
    ! ('cfl': branch-free, so the cell loops calling it vectorize)
    s_inv = 1.0_dp/drdxi_c(i) + 1.0_dp/(r_c(i)*Metric_DB(2, ib)%Dx_I(j)) &
        + 1.0_dp/(r_c(i)*sn_t(j)*Metric_DB(3, ib)%Dx_I(k))
    dti = tsc(i)*dt_now
    ch = merge(max(ch, glm_cfl_frac*cfl/(s_inv*dti)), ch, glm_cfl .and. dt_now > 0.0_dp)
  end function ch_cell

  ! The bound block's step factors (local time stepping) with those of the
  ! cells below and above it (cells 0:nb(1)+1; at a wall the end cell's own)
  function shell_steps() result(sc)
    real(dp) :: sc(0:nb(1)+1)
    integer :: m
    m = nb(1)
    sc(1:m) = tsc(1:m)
    sc(0) = merge(tsc(1), tsc_g(max(r_off, 1)), at_bottom)
    sc(m+1) = merge(tsc(m), tsc_g(min(r_off + m + 1, n_r)), at_top)
  end function shell_steps

  ! The Alfven-speed cap (va_max > 0, after Rempel 2009/2017): the Lorentz
  ! force's factor f = 1/sqrt(1 + (v_A/va_max)^4) for |B|^2 = b2 and the
  ! physical density rho; 1 without the cap
  ! (cap: the cell's cap from cell_cap instead of va_max)
  pure real(dp) function lorentz_f(b2, rho, cap) result(f)
    real(dp), intent(in) :: b2, rho
    real(dp), intent(in), optional :: cap
    real(dp) :: c
    c = va_max
    if (present(cap)) c = cap
    ! (branch-free, so the cell loops calling it vectorize; c = 0: 1)
    f = merge(1.0_dp/sqrt(1.0_dp + (b2/(four_pi*rho*c**2))**2), 1.0_dp, c > 0.0_dp)
  end function lorentz_f

  ! The Alfven speed the scheme sees: sqrt(f) v_A (-> va_max for v_A >> va_max)
  pure real(dp) function va_eff(b2, rho, cap) result(va)
    real(dp), intent(in) :: b2, rho
    real(dp), intent(in), optional :: cap
    va = sqrt(lorentz_f(b2, rho, cap)*b2/(four_pi*rho))
  end function va_eff

  ! The Alfven-speed cap of cell (i, j, k) of the bound block ib (cs: its signal
  ! speed c_s/xi): va_max, or with va_dt_target > 0 the speed that keeps the
  ! cell's step at the target, max(cfl dx_min/va_dt_target - (cs + |v|),
  ! va_floor_cs cs S(z), 1 km/s)
  real(dp) function cell_cap(i, j, k, ib, cs) result(c)
    integer, intent(in) :: i, j, k, ib
    real(dp), intent(in) :: cs
    real(dp) :: dxm, fl
    ! (branch-free, so the cell loops calling it vectorize: sw_r = 0 without
    ! the atmosphere; va_dt_target = 0: va_max. sw_r, sn_t: cz_bind's, of block ib)
    dxm = min(drc_r(i), r_c(i)*Metric_DB(2, ib)%Dx_I(j), r_c(i)*sn_t(j)*Metric_DB(3, ib)%Dx_I(k))
    fl = max(1.0e5_dp, va_floor_cs*cs*sw_r(i))
    c = merge(max(cfl*dxm/va_dt_target - (cs + norm2(State_GVB(i, j, k, vr_:vp_, ib))), fl), va_max, &
        va_dt_target > 0.0_dp)
  end function cell_cap

  ! va_eff(b2, rho, cell_cap(i, j, k, ib, cs)) from scalars only (the same
  ! arithmetic): |v| = vm, dx_min = dxm, the atmosphere switch sw; small and
  ! pure so the cell loops inline it and vectorize (cell_cap is not inlined)
  pure real(dp) function va_cap(b2, rho, cs, vm, dxm, sw) result(va)
    real(dp), intent(in) :: b2, rho, cs, vm, dxm, sw
    real(dp) :: c, f
    c = merge(max(cfl*dxm/va_dt_target - (cs + vm), max(1.0e5_dp, va_floor_cs*cs*sw)), va_max, va_dt_target > 0.0_dp)
    f = merge(1.0_dp/sqrt(1.0_dp + (b2/(four_pi*rho*c**2))**2), 1.0_dp, c > 0.0_dp)
    va = sqrt(f*b2/(four_pi*rho))
  end function va_cap

  ! The cell's signal speed from the pressure p0 + p1 (Gamma1 ~ 5/3; for the
  ! cap only), reduced by xi
  real(dp) function cell_cs(i, j, k, rho) result(cs)
    integer, intent(in) :: i, j, k
    real(dp), intent(in) :: rho
    cs = sqrt(5.0_dp/3.0_dp*max(p0(i) + p1_G(i, j, k), 1.0e-30_dp)/rho)/xi(i)
  end function cell_cs

  ! The Coriolis impulse of this step on L_z, each shell with its own step
  ! (local stepping): sum of -2 Omega rho s v_s dV tsc(i) dt; also sets the
  ! reference L_z at the first call (the run's start or restart)
  subroutine lz_impulse(dt)
    real(dp), intent(in) :: dt
    integer :: jb, ib, i, j, k
    real(dp) :: e(3,3), x(3), vc(3), w, th, ph, th2, ph2, imp, lz
    if (.not. lz_ref_set) then
      call angular_momentum(lz)
      lz_ref = lz;  lz_imp = 0.0_dp;  lz_ref_set = .true.
      t_next_lzc = 0.0_dp
    end if
    if (omega_rot == 0.0_dp) return
    imp = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          th = Metric_DB(2, ib)%X_I(j);  ph = Metric_DB(3, ib)%X_I(k)
          if (.not. yy_core(th, ph)) cycle
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          call yy_angles(th, ph, th2, ph2)
          if (yy_core(th2, ph2)) w = 0.5_dp*w
          e = sph_basis(th, ph)
          do i = 1, nb(1)
            x = r_c(i)*e(:, 1)
            vc = matmul(e, State_GVB(i, j, k, vr_:vp_, ib))
            if (Key_IB(KP, ib) == 2) then
              x = yy_cart(x);  vc = yy_cart(vc)
            end if
            imp = imp + w*Metric_DB(1, ib)%Vol_I(i)*(rho0(i) + State_GVB(i, j, k, rho1_, ib)) &
                *(x(1)*vc(1) + x(2)*vc(2))*tsc(i)
          end do
        end do
      end do
    end do
    call cz_unbind()
    lz_imp = lz_imp - 2.0_dp*omega_rot*dt*arbor_sum(imp)
  end subroutine lz_impulse

  ! the Yin-Yang core (theta in [pi/4, 3pi/4], |phi| <= 3pi/4)
  pure logical function yy_core(t, p)
    real(dp), intent(in) :: t, p
    yy_core = t >= 0.25_dp*pi .and. t <= 0.75_dp*pi .and. abs(p) <= 0.75_dp*pi
  end function yy_core

  ! Angular momentum of the flow about the global z axis, L_z = int rho (x
  ! cross v)_z dV with rho = rho0 + rho1 and the global Cartesian position and
  ! velocity (Yang mapped by yy_cart); Yin-Yang overlap weights as the budget's
  subroutine angular_momentum(lz, tcor)
    real(dp), intent(out) :: lz
    real(dp), intent(out), optional :: tcor
    integer :: jb, ib, i, j, k
    real(dp) :: e(3,3), x(3), vc(3), w, th, ph, th2, ph2, svs
    lz = 0.0_dp;  svs = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          th = Metric_DB(2, ib)%X_I(j);  ph = Metric_DB(3, ib)%X_I(k)
          if (.not. core(th, ph)) cycle
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          call yy_angles(th, ph, th2, ph2)
          if (core(th2, ph2)) w = 0.5_dp*w
          e = sph_basis(th, ph)
          do i = 1, nb(1)
            x = r_c(i)*e(:, 1)
            vc = matmul(e, State_GVB(i, j, k, vr_:vp_, ib))
            if (Key_IB(KP, ib) == 2) then
              x = yy_cart(x);  vc = yy_cart(vc)
            end if
            lz = lz + w*Metric_DB(1, ib)%Vol_I(i)*(rho0(i) + State_GVB(i, j, k, rho1_, ib))*(x(1)*vc(2) - x(2)*vc(1))
            ! s v_s = x vx + y vy (the distance from the axis times the velocity away from it)
            svs = svs + w*Metric_DB(1, ib)%Vol_I(i)*(rho0(i) + State_GVB(i, j, k, rho1_, ib))*(x(1)*vc(1) + x(2)*vc(2))
          end do
        end do
      end do
    end do
    call cz_unbind()
    lz = arbor_sum(lz)
    ! the Coriolis torque on the flow: -2 Omega int rho s v_s dV (= -Omega dI_z/dt)
    if (present(tcor)) tcor = -2.0_dp*omega_rot*arbor_sum(svs)
  contains
    logical function core(t, p)
      real(dp), intent(in) :: t, p
      core = t >= 0.25_dp*pi .and. t <= 0.75_dp*pi .and. abs(p) <= 0.75_dp*pi
    end function core
  end subroutine angular_momentum

  ! GLM: psi damped over the step by exp(-dt alpha_glm c_h / dx_min) (each
  ! shell's own step with local stepping)
  subroutine glm_damp(ib, dt)
    integer, intent(in) :: ib
    real(dp), intent(in) :: dt
    integer :: i, j, k
    real(dp) :: ch, dxm, sth
    do k = 1, nb(3)
      do j = 1, nb(2)
        sth = sin(Metric_DB(2, ib)%X_I(j))
        do i = 1, nb(1)
          ch = ch_cell(i, j, k, ib)
          dxm = min(drdxi_c(i), r_c(i)*Metric_DB(2, ib)%Dx_I(j), r_c(i)*sth*Metric_DB(3, ib)%Dx_I(k))
          State_GVB(i, j, k, psi_, ib) = State_GVB(i, j, k, psi_, ib)*exp(-tsc(i)*dt*alpha_glm*ch/dxm)
        end do
      end do
    end do
  end subroutine glm_damp

  ! ================================================================ radial transport

  ! Mixing length + grey in every column of the given blocks (one fill of
  ! the radial ghosts first). first: set the convective fluxes to their
  ! targets (call once at t = 0 with dt = 0, as CZTR).
  subroutine cz_transport_all(iBlocks, dt, first, time)
    integer, intent(in) :: iBlocks(:)
    real(dp), intent(in) :: dt
    logical, intent(in) :: first
    real(dp), intent(in), optional :: time          ! the blocks' time (multirate: for the ghosts of other levels)
    real(dp) :: rho(-1:nb(1)+2), T(-1:nb(1)+2), P(-1:nb(1)+2), cp(-1:nb(1)+2), cv(-1:nb(1)+2), Q(-1:nb(1)+2)
    real(dp) :: s1(-1:nb(1)+2), ss(-1:nb(1)+2), q6(-1:nb(1)+2, 6), lnrho(-1:nb(1)+2), Qrad(1:nb(1)), w
    real(dp) :: qc(0:nb(1)), tsx(0:nb(1)+1)
    integer :: jb, ib, i, j, k, m
    if (.not. (use_radiation .or. use_convection)) return
    m = nb(1)
    ! the columns use the interior cells and the radial ghosts only: with whole
    ! columns per block the wall conditions alone give them (no exchange);
    ! with radial blocks the ghosts between blocks need a fill
    if (m < n_r) then
      call tic(T_G_TR)
      if (si_cols .and. use_diffusion .and. .not. first) then
        ! the columns are whole on their rank (si_radial): the radial ghosts
        ! alone, copied between the column's blocks (what the fill's radial
        ! portals copy; the diffusion's fill later refreshes every ghost)
        call radial_ghosts_local()
      else if (present(time)) then
        call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.true., time=time)
      else
        call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.true.)
      end if
      call join_bc()
      call toc(T_G_TR)
    end if
    do jb = 1, size(iBlocks)
      ib = iBlocks(jb)
      if (use_atmosphere .and. Key_IB(KI+1, ib)*nb(1) >= k_join) then
        ! above the join: no grey, no MLT (their fluxes 0, for the profiles)
        if (qconv_face) then
          FaceState_FVDB(:, :, :, 1, 1, ib) = 0.0_dp
        else
          qconv_FB(:, :, :, ib) = 0.0_dp
        end if
        qrad_FB(:, :, :, ib) = 0.0_dp
        cycle
      end if
      call cz_bind(ib)
      if (m == n_r) then
        call wall_bc(State_GVB, 1, nv, ib, 1, 1, 2)
        call wall_bc(State_GVB, 1, nv, ib, 1, 2, 2)
      end if
      ! the step factors of the column part and of its neighbours in the
      ! blocks below and above (the faces between radial blocks take the
      ! smaller of their two cells', as inside a column)
      tsx(1:m) = tsc(1:m)
      tsx(0) = merge(tsc(1), tsc_g(max(r_off, 1)), at_bottom)
      tsx(m+1) = merge(tsc(m), tsc_g(min(r_off + m + 1, n_r)), at_top)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = -1, m + 2
            rho(i) = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)        ! physical density
            s1(i) = State_GVB(i, j, k, s1_, ib)
            ss(i) = s0(i) + s1(i)
          end do
          call eos_fields_n(tab, [F_LNP, F_LNT, F_GAMMA1, F_GAMMA3M1, F_LNCP, F_LNCV], m + 4, rho, ss, q6)
          call vlog(m + 4, rho, lnrho)
          call vexp(m + 4, q6(:, 1), P)
          call vexp(m + 4, q6(:, 2), T)
          call vexp(m + 4, q6(:, 5), cp)
          call vexp(m + 4, q6(:, 6), cv)
          do i = -1, m + 2
            Q(i) = q6(i, 4)*rho(i)*T(i)*cp(i)/(P(i)*q6(i, 3))
          end do
          if (qconv_face) then
            qc = FaceState_FVDB(0:m, j, k, 1, 1, ib)
          else
            qc = qconv_FB(:, j, k, ib)
          end if
          call transport_column(dt, rho, T, P, cp, cv, Q, lnrho, q6(:, 2), q6(:, 1), q6(:, 5), State_GVB(:, j, k, vr_, ib), &
              s1, qc, &
              f_mlt, alpha_mlt, cfl, L_bottom_Lsun*Lsun, grey_tau_above, use_radiation, use_convection, &
              use_hyperbolic, first, qrad_FB(:, j, k, ib), Qrad, wall_top == 'face' .or. use_atmosphere, tsx)
          if (budget_on .and. .not. first) then
            ! radiation (with L_bottom) and the exact cooling of each cell
            w = wyy_B(j, k, ib)
            do i = 1, m
              acc_IV(i, 3) = acc_IV(i, 3) - w*Metric_DB(1, ib)%Vol_I(i)*cool_c(i)*dt
              acc_IV(i, 2) = acc_IV(i, 2) + w*Metric_DB(1, ib)%Vol_I(i)*(rho(i)*T(i)*(s1(i) - State_GVB(i, j, k, s1_, ib)) &
                  + cool_c(i)*dt)
            end do
          end if
          State_GVB(1:m, j, k, s1_, ib) = s1(1:m)
          if (qconv_face) then
            FaceState_FVDB(0:m, j, k, 1, 1, ib) = qc
          else
            qconv_FB(:, j, k, ib) = qc
          end if
        end do
      end do
    end do
    call cz_unbind()
    if (grey_atm_blocks > 0 .and. use_radiation) call grey_atm_all(dt)
  end subroutine cz_transport_all

  ! Let MPI move the ghost messages in flight (Open MPI's shared memory and
  ! most networks send large messages by rendezvous: the data moves only
  ! while the ranks are inside MPI calls); a cheap probe does it
  subroutine mpi_progress()
    use mpi_f08, only: MPI_Iprobe, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_STATUS_IGNORE
    logical :: flag
    call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, flag, MPI_STATUS_IGNORE)
  end subroutine mpi_progress

  ! The two radial ghost layers of the transport's blocks (the CZ blocks;
  ! all blocks without the atmosphere) at their interior (j, k), copied from
  ! the blocks below and above in the same column (si_cb: whole columns on
  ! this rank), all variables; at the domain's bottom (and top) the walls
  ! (wall_bc, the fill's physical boundary). The same values the fill's
  ! radial portals and boundary give these cells.
  subroutine radial_ghosts_local()
    integer :: c, q, nk, ib, ibd, ibu, kb
    nk = si_kb1 - si_kb0 + 1
    if (si_kb0 /= 0 .or. (si_kb1 + 1)*nb(1) /= n_r) &
        call arbor_abort('radial_ghosts_local: the si_radial range must be the whole column')
    do c = 1, si_ncol
      do q = 1, nk
        ib = si_cb(q, c)
        kb = Key_IB(KI+1, ib)
        if (use_atmosphere .and. kb*nb(1) >= k_join) cycle          ! (above the join: no transport)
        if (q > 1) then
          ibd = si_cb(q - 1, c)
          State_GVB(1-ng:0, 1:nb(2), 1:nb(3), :, ib) = State_GVB(nb(1)-ng+1:nb(1), 1:nb(2), 1:nb(3), :, ibd)
        else
          call wall_bc(State_GVB, 1, nv, ib, 1, 1, ng)
        end if
        if (q < nk) then
          ibu = si_cb(q + 1, c)
          State_GVB(nb(1)+1:nb(1)+ng, 1:nb(2), 1:nb(3), :, ib) = State_GVB(1:ng, 1:nb(2), 1:nb(3), :, ibu)
        else
          call wall_bc(State_GVB, 1, nv, ib, 1, 2, ng)
        end if
      end do
    end do
  end subroutine radial_ghosts_local

  ! grey_atm_blocks: the grey solve in the first grey_atm_blocks radial blocks
  ! above the join, in radial order, each column part's bottom flux the
  ! radiative flux leaving the block below at their common face (qrad_FB
  ! of this step; the columns are whole on their rank: si_cb). Its heating
  ! goes into s1, its face fluxes into qrad_FB (profiles, the next block).
  subroutine grey_atm_all(dt)
    real(dp), intent(in) :: dt
    real(dp) :: rho(nb(1)), T(nb(1)), cv(nb(1)), lnrho(nb(1)), ss(nb(1)), q2(nb(1), 2), sc(nb(1)), Qrad(nb(1)), &
        qr(0:nb(1))
    integer :: kb0, kb, c, ib, ibd, i, j, k, m
    if (.not. si_cols) call arbor_abort('grey_atm_blocks: the columns are not distributed (si_radial)')
    m = nb(1)
    kb0 = k_join/m                                     ! the first radial block above the join
    if (si_kb0 > kb0 - 1 .or. si_kb1 < kb0 + grey_atm_blocks - 1) &
        call arbor_abort('grey_atm_blocks: the blocks must lie in the si_radial range (si_z_km)')
    do kb = kb0, min(kb0 + grey_atm_blocks - 1, n_r/m - 1)
      do c = 1, si_ncol
        ib = si_cb(kb - si_kb0 + 1, c)
        ibd = si_cb(kb - si_kb0, c)                    ! the block below in the column
        call cz_bind(ib)
        sc = tsc(1:m)
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, m
              rho(i) = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
              ss(i) = s0(i) + State_GVB(i, j, k, s1_, ib)
            end do
            call eos_fields_n(tab, [F_LNT, F_LNCV], m, rho, ss, q2)
            call vlog(m, rho, lnrho)
            call vexp(m, q2(:, 1), T)
            call vexp(m, q2(:, 2), cv)
            call grey_atm_column(dt, sc, qrad_FB(m, j, k, ibd), rho, T, cv, lnrho, q2(:, 1), grey_tau_above, qr, Qrad)
            do i = 1, m
              State_GVB(i, j, k, s1_, ib) = State_GVB(i, j, k, s1_, ib) + dt*sc(i)*Qrad(i)/(rho(i)*T(i))
            end do
            qrad_FB(:, j, k, ib) = qr
          end do
        end do
      end do
    end do
    call cz_unbind()
  end subroutine grey_atm_all

  ! ================================================================ diffusion

  ! Rempel-form diffusion (as CZTR's diffusion_rempel), in all three
  ! directions: at each face
  !   F_rho  = -cf J(rho1)                       (donor-side positivity limiter)
  !   F_rhov = v_u F_rho - rho_m cf J(v)         (each component)
  !   F_rhos = s_u F_rho - rho_m cf J(s1)
  ! J the minmod-limited, h-switched jump, cf = (c_L + c_R)/4 with c = |v_d| + cs,
  ! u the donor cell, rho_m the smaller density. The limiter factor of a
  ! ghost cell is 1 (the limiter never acts in the CZ; where it would, the
  ! fluxes of neighbouring blocks could differ).
  subroutine diffusion(ib, dt)
    integer, intent(in) :: ib
    real(dp), intent(in) :: dt
    integer :: i, j, k, c
    real(dp) :: wr(0:nb(1)), wh(1:nb(1))
    real(dp) :: vol, outflow, rho_new, a1, a2, a3, ho, s1o, ke_o, dE_act, dti, kh_o, eb_o
    ! sin theta at the cell centres (snt) and at the theta faces between them (snf: face j between j and j+1)
    real(dp) :: snt(0:nb(2)+1), snf(0:nb(2))
    logical :: blend_r, blend_h, b_resist, b_j4
    call tic(T_D_PRE)
    ! the top-layer blend weights at the radial faces (wr) and the cells (wh: horizontal faces)
    wr = 0.0_dp;  wh = 0.0_dp
    if (diff_top_km > 0.0_dp) then
      do i = 0, nb(1)
        wr(i) = top_weight(r_top_cz - r_f(i))
        if (i >= 1) wh(i) = top_weight(r_top_cz - r_c(i))
      end do
    end if
    ! the top-layer blend only where its weights are nonzero (adding 0 x the
    ! blend elsewhere changes nothing): radial faces (wr), horizontal faces (wh)
    blend_r = diff_top_km > 0.0_dp .and. any(wr /= 0.0_dp)
    blend_h = diff_top_km > 0.0_dp .and. any(wh /= 0.0_dp)
    b_resist = use_B .and. b_diff_form == 'resistive'
    b_j4 = b_diffusion_h < 0.0_dp .and. b_diff_limiter == 'minmod'
    if (b_resist) then
      associate (th => Metric_DB(2, ib)%X_I)
        do j = 0, nb(2) + 1
          snt(j) = sin(th(j))
        end do
        do j = 0, nb(2)
          snf(j) = sin(0.5_dp*(th(j) + th(j+1)))
        end do
      end associate
      if (.not. allocated(dB_C)) allocate(dB_C(0:nb(1)+1, 0:nb(2)+1, 0:nb(3)+1, 6))
      call bcell_derivs()
    end if
    call toc(T_D_PRE)
    call tic(T_D_EOS)
    ! (above the join T and p1 are not used: no energy flux, no heating there; the budget needs them)
    call derived(ib, merge(D_CS, D_P1 + D_CS + D_T, in_atm .and. .not. budget_on), C_FACES1)
    call toc(T_D_EOS)
    call tic(T_D_FA)
    associate (m1 => Metric_DB(1, ib), m2 => Metric_DB(2, ib), m3 => Metric_DB(3, ib))
      ! face speeds and mass fluxes (unlimited), then the limiter
      call faces_a(1, 1, 0, 0, F1_FV, C1_F)
      call faces_a(2, 0, 1, 0, F2_FV, C2_F)
      call faces_a(3, 0, 0, 1, F3_FV, C3_F)
      lim_G = 1.0_dp
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            vol = m1%Vol_I(i)*m2%Vol_I(j)*m3%Vol_I(k)
            a1 = m2%Tang_ID(j, 1)*m3%Tang_ID(k, 1)
            a2 = m1%Tang_ID(i, 2)*m3%Tang_ID(k, 2)
            a3 = m1%Tang_ID(i, 3)*m2%Tang_ID(j, 3)
            outflow = tsc(i)*dt*(a1*(m1%Area_F(i)*max(F1_FV(i, j, k, 1), 0.0_dp) + m1%Area_F(i-1)*max(-F1_FV(i-1, j, k, 1), 0.0_dp)) &
                + a2*(m2%Area_F(j)*max(F2_FV(i, j, k, 1), 0.0_dp) + m2%Area_F(j-1)*max(-F2_FV(i, j-1, k, 1), 0.0_dp)) &
                + a3*(m3%Area_F(k)*max(F3_FV(i, j, k, 1), 0.0_dp) + m3%Area_F(k-1)*max(-F3_FV(i, j, k-1, 1), 0.0_dp)))/vol
            if (outflow > 0.5_dp*rho_G(i, j, k)) lim_G(i, j, k) = 0.5_dp*rho_G(i, j, k)/outflow
          end do
        end do
      end do
      call toc(T_D_FA)
      ! all fluxes with the donor values
      call tic(T_D_FB)
      call faces_b(1, 0, 0, F1_FV, C1_F)
      call faces_b(0, 1, 0, F2_FV, C2_F)
      call faces_b(0, 0, 1, F3_FV, C3_F)
    end associate
    call toc(T_D_FB)
    call tic(T_D_FR)
    call wall_fluxes(1, nVarCz)
    call radial_energy_flux(ib)
    call toc(T_D_FR)
    call tic(T_D_FE)
    if (diffusion_heating .and. .not. in_atm) then
      ! the energy each face's fluxes carry (variable 6), with h + Phi per cell
      ! (CZ blocks; above the join none, see below)
      ! (on the cells derived's faces1 evaluates: the interior and one layer across each face)
      do k = 0, nb(3) + 1
        do j = 0, nb(2) + 1
          if ((j == 0 .or. j == nb(2) + 1) .and. (k == 0 .or. k == nb(3) + 1)) cycle
          if (j >= 1 .and. j <= nb(2) .and. k >= 1 .and. k <= nb(3)) then
            do i = 0, nb(1) + 1
              hP_G(i, j, k) = h0_r(i) + T_G(i, j, k)*State_GVB(i, j, k, s1_, ib) + p1_G(i, j, k)/rho_G(i, j, k) + phi_r(i)
            end do
          else
            do i = 1, nb(1)
              hP_G(i, j, k) = h0_r(i) + T_G(i, j, k)*State_GVB(i, j, k, s1_, ib) + p1_G(i, j, k)/rho_G(i, j, k) + phi_r(i)
            end do
          end if
        end do
      end do
      call faces_e(1, 0, 0, F1_FV)
      call faces_e(0, 1, 0, F2_FV)
      call faces_e(0, 0, 1, F3_FV)
      ! with the atmosphere: only in the CZ (cells up to the join; no energy
      ! flux above the join face; its enthalpy is linearized about the
      ! background, which fails in the atmosphere: CZTR_2D). The join face
      ! keeps its energy flux.
      if (use_atmosphere) then
        do i = max(0, k_join - r_off + 1), nb(1)
          F1_FV(i, :, :, ie_) = 0.0_dp
        end do
        do i = max(1, k_join - r_off + 1), nb(1)
          F2_FV(i, :, :, ie_) = 0.0_dp;  F3_FV(i, :, :, ie_) = 0.0_dp
        end do
      end if
    else
      F1_FV(:, :, :, ie_) = 0.0_dp;  F2_FV(:, :, :, ie_) = 0.0_dp;  F3_FV(:, :, :, ie_) = 0.0_dp
    end if
    call wall_fluxes(ie_, ie_)
    call toc(T_D_FE)
    call tic(T_D_DIV)
    call arbor_flux_divergence(ib, F1_FV, F2_FV, F3_FV, Div_CV)
    call toc(T_D_DIV)
    call tic(T_D_UPD)
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, nb(1)
          dti = tsc(i)*dt                                   ! (the shell's own step)
          if (test_b_diff_only) then
            if (use_B) State_GVB(i, j, k, br_:psi_, ib) = State_GVB(i, j, k, br_:psi_, ib) - dti*Div_CV(i, j, k, br_:psi_)
            cycle
          end if
          rho_new = rho_G(i, j, k) - dti*Div_CV(i, j, k, 1)
          s1o = State_GVB(i, j, k, s1_, ib)
          if (rsst_h2015) then
            ! rho dq/dt = -div F_q for q = v, s1; energy rho T s1 + rho v^2/2
            ke_o = 0.5_dp*rho_G(i, j, k)*sum(State_GVB(i, j, k, vr_:vp_, ib)**2)
            State_GVB(i, j, k, rho1_, ib) = rho_new - rho0(i)
            do c = vr_, s1_
              State_GVB(i, j, k, c, ib) = State_GVB(i, j, k, c, ib) - dti*Div_CV(i, j, k, c)/rho_G(i, j, k)
            end do
            if (diffusion_heating .and. (.not. use_atmosphere .or. r_off + i <= k_join)) then
              dE_act = T_G(i, j, k)*(rho_new*State_GVB(i, j, k, s1_, ib) - rho_G(i, j, k)*s1o) &
                  + 0.5_dp*rho_new*sum(State_GVB(i, j, k, vr_:vp_, ib)**2) - ke_o
              State_GVB(i, j, k, s1_, ib) = State_GVB(i, j, k, s1_, ib) &
                  + (-dti*Div_CV(i, j, k, ie_) - dE_act)/(rho_new*T_G(i, j, k))
            end if
            cycle
          end if
          if (.not. in_atm .or. budget_on) ho = h0_r(i) + T_G(i, j, k)*s1o + p1_G(i, j, k)/rho_G(i, j, k) + phi_r(i)   ! h + Phi
          ke_o = 0.5_dp*rho_G(i, j, k)*sum(State_GVB(i, j, k, vr_:vp_, ib)**2)
          kh_o = rho_G(i, j, k)*sum(State_GVB(i, j, k, vt_:vp_, ib)**2)
          State_GVB(i, j, k, rho1_, ib) = rho_new - rho0(i)
          do c = vr_, vp_
            State_GVB(i, j, k, c, ib) = (rho_G(i, j, k)*State_GVB(i, j, k, c, ib) - dti*Div_CV(i, j, k, c))/rho_new
          end do
          State_GVB(i, j, k, s1_, ib) = (rho_G(i, j, k)*(s0(i) + State_GVB(i, j, k, s1_, ib)) - dti*Div_CV(i, j, k, s1_)) &
              /rho_new - s0(i)
          if (use_B) then
            eb_o = sum(State_GVB(i, j, k, br_:bp_, ib)**2)/(2.0_dp*four_pi)
            do c = br_, psi_
              State_GVB(i, j, k, c, ib) = State_GVB(i, j, k, c, ib) - dti*Div_CV(i, j, k, c)
            end do
          end if
          if (diffusion_heating .and. (.not. use_atmosphere .or. r_off + i <= k_join)) then
            ! the cell's energy change to first order, corrected to -dt div F_E by its entropy
            dE_act = rho_G(i, j, k)*T_G(i, j, k)*(State_GVB(i, j, k, s1_, ib) - s1o) + ho*(rho_new - rho_G(i, j, k)) &
                + 0.5_dp*rho_new*sum(State_GVB(i, j, k, vr_:vp_, ib)**2) - ke_o
            if (use_B) dE_act = dE_act + sum(State_GVB(i, j, k, br_:bp_, ib)**2)/(2.0_dp*four_pi) - eb_o
            State_GVB(i, j, k, s1_, ib) = State_GVB(i, j, k, s1_, ib) &
                + (-dti*Div_CV(i, j, k, ie_) - dE_act)/(rho_G(i, j, k)*T_G(i, j, k))
          end if
          if (budget_on) then
            acc_IV(i, 17) = acc_IV(i, 17) + wyy_B(j, k, ib)*Metric_DB(1, ib)%Vol_I(i)*0.5_dp &
                *(rho_new*sum(State_GVB(i, j, k, vt_:vp_, ib)**2) - kh_o)
            acc_IV(i, 18) = acc_IV(i, 18) + wyy_B(j, k, ib)*Metric_DB(1, ib)%Vol_I(i) &
                *(0.5_dp*rho_new*sum(State_GVB(i, j, k, vr_:vp_, ib)**2) - ke_o)
          end if
          if (budget_on) acc_IV(i, 4) = acc_IV(i, 4) + wyy_B(j, k, ib)*Metric_DB(1, ib)%Vol_I(i) &
              *(rho_G(i, j, k)*T_G(i, j, k)*(State_GVB(i, j, k, s1_, ib) - s1o) + ho*(rho_new - rho_G(i, j, k)) &
              + 0.5_dp*rho_new*sum(State_GVB(i, j, k, vr_:vp_, ib)**2) - ke_o)
        end do
      end do
    end do

    call toc(T_D_UPD)
  contains

    ! 1 above the depth diff_top_km, a cos^2 ramp to 0 at 2 diff_top_km
    pure real(dp) function top_weight(depth) result(w)
      real(dp), intent(in) :: depth
      real(dp) :: x
      x = depth/(diff_top_km*1.0e5_dp)
      if (depth < 0.0_dp .and. use_atmosphere) then
        ! above the photosphere: down to 0 at diff_top_above_km
        w = cos(0.5_dp*pi*min(1.0_dp, -depth/(diff_top_above_km*1.0e5_dp)))**2
      else if (x <= 1.0_dp) then
        w = 1.0_dp
      else if (x >= 2.0_dp) then
        w = 0.0_dp
      else
        w = cos(0.5_dp*pi*(x - 1.0_dp))**2
      end if
    end function top_weight

    ! The radial fluxes at the walls from the symmetry, not the ghosts (exact
    ! wall conservation whatever the ghosts): a cell-centred wall's outer face
    ! carries the mirror of the face below it (a variable even about the wall:
    ! the opposite flux; v_r, odd: the same), so the wall cell keeps half its
    ! budget; a face wall passes no flux of the even variables. Variables
    ! c1..c2 (the hydro variables and the energy; B keeps its ghost fluxes).
    subroutine wall_fluxes(c1, c2)
      integer, intent(in) :: c1, c2
      integer :: c, m
      real(dp) :: par
      m = nb(1)
      do c = c1, c2
        par = merge(-1.0_dp, 1.0_dp, c == vr_)
        if (at_top .and. wall_top /= 'outflow') then
          if (wall_top == 'cell') then
            F1_FV(m, :, :, c) = -par*F1_FV(m - 1, :, :, c)
          else if (c /= vr_) then
            F1_FV(m, :, :, c) = 0.0_dp
          end if
        end if
        if (at_bottom) then
          if (wall_bottom == 'cell') then
            F1_FV(0, :, :, c) = -par*F1_FV(1, :, :, c)
          else if (c /= vr_) then
            F1_FV(0, :, :, c) = 0.0_dp
          end if
        end if
      end do
    end subroutine wall_fluxes

    ! The face passes of one direction: faces between cells (i, j, k) and
    ! (i+di, j+dj, k+dk), loops with i innermost that vectorize. Each face's
    ! arithmetic is the per-face formula, unchanged.

    ! A: face speed cf = (c_L + c_R)/4, c = |v_d| + cs, and the unlimited mass
    ! flux F_rho = -cf J(rho1)
    subroutine faces_a(d, di, dj, dk, F, C)
      integer, intent(in) :: d, di, dj, dk
      real(dp), intent(inout) :: F(1-di:nb(1), 1-dj:nb(2), 1-dk:nb(3), ie_)
      real(dp), intent(out) :: C(1-di:nb(1), 1-dj:nb(2), 1-dk:nb(3))
      integer :: i, j, k, vd
      vd = vr_ + d - 1
      do k = 1 - dk, nb(3)
        do j = 1 - dj, nb(2)
          do i = 1 - di, nb(1)
            C(i, j, k) = 0.25_dp*(abs(State_GVB(i, j, k, vd, ib)) + cs_G(i, j, k) &
                + abs(State_GVB(i+di, j+dj, k+dk, vd, ib)) + cs_G(i+di, j+dj, k+dk))
          end do
        end do
      end do
      ! si_radial: the radial face speed capped at si_diff_cap dr/dt where
      ! either cell is implicit (the explicit diffusion stays stable)
      if (si_radial .and. di == 1) then
        do i = 0, nb(1)
          if (r_off + i + 1 < si_i0 + 2 .or. r_off + i > si_i1 - 2) cycle
          C(i, :, :) = min(C(i, :, :), si_diff_cap*(r_c(i+1) - r_c(i))/dt)
        end do
      end if
      do k = 1 - dk, nb(3)
        do j = 1 - dj, nb(2)
          do i = 1 - di, nb(1)
            F(i, j, k, rho1_) = -C(i, j, k)*jump4(State_GVB(i-di, j-dj, k-dk, rho1_, ib), State_GVB(i, j, k, rho1_, ib), &
                State_GVB(i+di, j+dj, k+dk, rho1_, ib), State_GVB(i+2*di, j+2*dj, k+2*dk, rho1_, ib))
          end do
        end do
      end do
      ! the top layer: blended toward the full jump (diff_top_km)
      if (merge(blend_r, blend_h, di == 1)) then
        do k = 1 - dk, nb(3)
          do j = 1 - dj, nb(2)
            do i = 1 - di, nb(1)
              F(i, j, k, rho1_) = F(i, j, k, rho1_) + merge(wr(max(i, 0)), wh(max(i, 1)), di == 1) &
                  *(-C(i, j, k)*(State_GVB(i+di, j+dj, k+dk, rho1_, ib) - State_GVB(i, j, k, rho1_, ib)) - F(i, j, k, rho1_))
            end do
          end do
        end do
      end if
    end subroutine faces_a

    ! B: the mass flux limited by its donor cell (F_rho >= 0: the left cell),
    ! then F_q = q_u F_rho - rho_m cf J(q) for q = v, s0 + s1 (u the donor,
    ! rho_m the smaller density), -rho_m cf J(q) for hotta2015, and
    ! -b_diff_factor cf J(q) for the field and psi. (The limiter factor is
    ! positive: the limited flux keeps the sign, so it still names the donor.)
    subroutine faces_b(di, dj, dk, F, C)
      integer, intent(in) :: di, dj, dk
      real(dp), intent(inout) :: F(1-di:nb(1), 1-dj:nb(2), 1-dk:nb(3), ie_)
      real(dp), intent(in) :: C(1-di:nb(1), 1-dj:nb(2), 1-dk:nb(3))
      integer :: i, j, k, v
      real(dp) :: rho_m, jq
      logical :: up
      do k = 1 - dk, nb(3)
        do j = 1 - dj, nb(2)
          do i = 1 - di, nb(1)
            F(i, j, k, rho1_) = merge(lim_G(i, j, k), lim_G(i+di, j+dj, k+dk), F(i, j, k, rho1_) >= 0.0_dp)*F(i, j, k, rho1_)
          end do
        end do
      end do
      do v = vr_, s1_
        do k = 1 - dk, nb(3)
          do j = 1 - dj, nb(2)
            do i = 1 - di, nb(1)
              up = F(i, j, k, rho1_) >= 0.0_dp
              rho_m = min(rho_G(i, j, k), rho_G(i+di, j+dj, k+dk))
              jq = jump4(State_GVB(i-di, j-dj, k-dk, v, ib), State_GVB(i, j, k, v, ib), &
                  State_GVB(i+di, j+dj, k+dk, v, ib), State_GVB(i+2*di, j+2*dj, k+2*dk, v, ib))
              if (rsst_h2015) then
                ! each variable diffused on its own (Hotta et al. 2015 App. B): the
                ! diffused mass is the scaled xi^2 rho1 and must not carry v and s
                F(i, j, k, v) = -rho_m*C(i, j, k)*jq
              else if (v == s1_) then
                F(i, j, k, v) = merge(s0(i) + State_GVB(i, j, k, v, ib), s0(i+di) + State_GVB(i+di, j+dj, k+dk, v, ib), up) &
                    *F(i, j, k, rho1_) - rho_m*C(i, j, k)*jq
              else
                F(i, j, k, v) = merge(State_GVB(i, j, k, v, ib), State_GVB(i+di, j+dj, k+dk, v, ib), up)*F(i, j, k, rho1_) &
                    - rho_m*C(i, j, k)*jq
              end if
            end do
          end do
        end do
      end do
      ! the top layer: the diffusive part blended toward the full jump (diff_top_km)
      if (merge(blend_r, blend_h, di == 1)) then
        do v = vr_, s1_
          do k = 1 - dk, nb(3)
            do j = 1 - dj, nb(2)
              do i = 1 - di, nb(1)
                rho_m = min(rho_G(i, j, k), rho_G(i+di, j+dj, k+dk))
                jq = jump4(State_GVB(i-di, j-dj, k-dk, v, ib), State_GVB(i, j, k, v, ib), &
                    State_GVB(i+di, j+dj, k+dk, v, ib), State_GVB(i+2*di, j+2*dj, k+2*dk, v, ib))
                F(i, j, k, v) = F(i, j, k, v) - rho_m*C(i, j, k)*merge(wr(max(i, 0)), wh(max(i, 1)), di == 1) &
                    *((State_GVB(i+di, j+dj, k+dk, v, ib) - State_GVB(i, j, k, v, ib)) - jq)
              end do
            end do
          end do
        end do
      end if
      ! the magnetic field and psi: -c J(q) (no density weight); b_diff_form =
      ! 'resistive': B's fluxes -eta J (see bfaces), none of the normal component
      do v = nVarCz + 1, nv
        if (b_resist .and. v <= bp_) then
          ! -eta J_dm, eta the face's mean of the cells' eta; on a radial
          ! face only the mixed term (the normal one: b_radial_implicit),
          ! both with b_eta_radial = 'explicit'; none of the normal component
          if (v - br_ + 1 == di + 2*dj + 3*dk) then
            F(:, :, :, v) = 0.0_dp
            cycle
          end if
          call bfaces(di + 2*dj + 3*dk, v - br_ + 1, F(:, :, :, v))
        else if (b_j4) then
          ! (jump4b with b_diffusion_h < 0 and the minmod limiter is jump4: a loop that vectorizes)
          do k = 1 - dk, nb(3)
            do j = 1 - dj, nb(2)
              do i = 1 - di, nb(1)
                jq = jump4(State_GVB(i-di, j-dj, k-dk, v, ib), State_GVB(i, j, k, v, ib), &
                    State_GVB(i+di, j+dj, k+dk, v, ib), State_GVB(i+2*di, j+2*dj, k+2*dk, v, ib))
                F(i, j, k, v) = -b_diff_factor*C(i, j, k)*jq
              end do
            end do
          end do
        else
          do k = 1 - dk, nb(3)
            do j = 1 - dj, nb(2)
              do i = 1 - di, nb(1)
                jq = jump4b(State_GVB(i-di, j-dj, k-dk, v, ib), State_GVB(i, j, k, v, ib), &
                    State_GVB(i+di, j+dj, k+dk, v, ib), State_GVB(i+2*di, j+2*dj, k+2*dk, v, ib))
                F(i, j, k, v) = -b_diff_factor*C(i, j, k)*jq
              end do
            end do
          end do
        end if
      end do
    end subroutine faces_b

    ! b_diff_form = 'resistive': the face fluxes -eta_f (qn + qt)/ds of B_m
    ! across the faces normal to d (m /= d), i.e. -eta J with
    ! ds_d (d_d B_m - d_m B_d) in spherical form, ds_d/(h_d h_m) [d_d(h_m B_m)
    ! - d_m(h_d B_d)], h = (1, r, r sin theta): qn the normal derivative's term
    ! (a compact jump across the face), qt the tangential one (the mean of the
    ! two cells' central differences, dB_C from bcell_derivs); ds = ds_d the
    ! distance of the two centres, eta_f the mean of the two cells' eta. On a
    ! radial face qn only with b_eta_radial = 'explicit' (implicit:
    ! b_radial_implicit). Fv: the faces of one direction (as F in faces_b).
    ! (State_GVB and EtaB_GB indexed directly: an associate name for an array
    ! section would start at 1, not at the ghosts' 1 - ng)
    subroutine bfaces(d, m, Fv)
      integer, intent(in) :: d, m
      real(dp), intent(out) :: Fv(1-merge(1, 0, d == 1):nb(1), 1-merge(1, 0, d == 2):nb(2), 1-merge(1, 0, d == 3):nb(3))
      integer :: i, j, k
      real(dp) :: dx, rf, sf, a, qn, qt, ds
      associate (th => Metric_DB(2, ib)%X_I, ph => Metric_DB(3, ib)%X_I)
      select case (10*d + m)
      case (12)                                           ! radial faces, B_theta: d_r(r B_theta) - d_theta B_r
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 0, nb(1)
              dx = r_c(i+1) - r_c(i)
              rf = 0.5_dp*(r_c(i) + r_c(i+1))
              a = dB_C(i, j, k, 1) + dB_C(i+1, j, k, 1)
              qn = merge((r_c(i+1)*State_GVB(i+1, j, k, bt_, ib) - r_c(i)*State_GVB(i, j, k, bt_, ib))/rf, 0.0_dp, b_eta_explicit)
              qt = -0.5_dp*dx*a/rf
              Fv(i, j, k) = -0.5_dp*(EtaB_GB(i, j, k, 1, ib) + EtaB_GB(i+1, j, k, 1, ib))*(qn + qt)/dx
            end do
          end do
        end do
      case (13)                                           ! radial faces, B_phi: d_r(r B_phi) - d_phi B_r/sin
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 0, nb(1)
              dx = r_c(i+1) - r_c(i)
              rf = 0.5_dp*(r_c(i) + r_c(i+1))
              a = dB_C(i, j, k, 2) + dB_C(i+1, j, k, 2)
              qn = merge((r_c(i+1)*State_GVB(i+1, j, k, bp_, ib) - r_c(i)*State_GVB(i, j, k, bp_, ib))/rf, 0.0_dp, b_eta_explicit)
              qt = -0.5_dp*dx*a/(rf*snt(j))
              Fv(i, j, k) = -0.5_dp*(EtaB_GB(i, j, k, 1, ib) + EtaB_GB(i+1, j, k, 1, ib))*(qn + qt)/dx
            end do
          end do
        end do
      case (21)                                           ! theta faces, B_r: d_theta B_r - d_r(r B_theta)
        do k = 1, nb(3)
          do j = 0, nb(2)
            dx = th(j+1) - th(j)
            do i = 1, nb(1)
              ds = r_c(i)*dx
              a = dB_C(i, j, k, 3) + dB_C(i, j+1, k, 3)
              qn = State_GVB(i, j+1, k, br_, ib) - State_GVB(i, j, k, br_, ib)
              qt = -0.5_dp*dx*a
              Fv(i, j, k) = -0.5_dp*(EtaB_GB(i, j, k, 1, ib) + EtaB_GB(i, j+1, k, 1, ib))*(qn + qt)/ds
            end do
          end do
        end do
      case (23)                                           ! theta faces, B_phi: d_theta(sin B_phi)/sin - d_phi B_theta/sin
        do k = 1, nb(3)
          do j = 0, nb(2)
            dx = th(j+1) - th(j)
            sf = snf(j)
            do i = 1, nb(1)
              ds = r_c(i)*dx
              a = dB_C(i, j, k, 4) + dB_C(i, j+1, k, 4)
              qn = (snt(j+1)*State_GVB(i, j+1, k, bp_, ib) - snt(j)*State_GVB(i, j, k, bp_, ib))/sf
              qt = -0.5_dp*dx*a/sf
              Fv(i, j, k) = -0.5_dp*(EtaB_GB(i, j, k, 1, ib) + EtaB_GB(i, j+1, k, 1, ib))*(qn + qt)/ds
            end do
          end do
        end do
      case (31)                                           ! phi faces, B_r: d_phi B_r - sin d_r(r B_phi)
        do k = 0, nb(3)
          dx = ph(k+1) - ph(k)
          do j = 1, nb(2)
            do i = 1, nb(1)
              ds = r_c(i)*snt(j)*dx
              a = dB_C(i, j, k, 5) + dB_C(i, j, k+1, 5)
              qn = State_GVB(i, j, k+1, br_, ib) - State_GVB(i, j, k, br_, ib)
              qt = -0.5_dp*dx*snt(j)*a
              Fv(i, j, k) = -0.5_dp*(EtaB_GB(i, j, k, 1, ib) + EtaB_GB(i, j, k+1, 1, ib))*(qn + qt)/ds
            end do
          end do
        end do
      case (32)                                           ! phi faces, B_theta: d_phi B_theta - d_theta(sin B_phi)
        do k = 0, nb(3)
          dx = ph(k+1) - ph(k)
          do j = 1, nb(2)
            do i = 1, nb(1)
              ds = r_c(i)*snt(j)*dx
              a = dB_C(i, j, k, 6) + dB_C(i, j, k+1, 6)
              qn = State_GVB(i, j, k+1, bt_, ib) - State_GVB(i, j, k, bt_, ib)
              qt = -0.5_dp*dx*a
              Fv(i, j, k) = -0.5_dp*(EtaB_GB(i, j, k, 1, ib) + EtaB_GB(i, j, k+1, 1, ib))*(qn + qt)/ds
            end do
          end do
        end do
      end select
      end associate
    end subroutine bfaces

    ! the cells' central differences the face currents use (dB_C(:, :, :, n)),
    ! each once per block, where the faces need them (one cell beyond the
    ! block in the faces' direction): 1 d_theta B_r, 2 d_phi B_r (cells i =
    ! 0..nb+1; radial faces), 3 d_r(r B_theta), 4 d_phi B_theta (j = 0..nb+1;
    ! theta faces), 5 d_r(r B_phi), 6 d_theta(sin theta B_phi) (k = 0..nb+1; phi faces)
    subroutine bcell_derivs()
      integer :: i, j, k
      associate (th => Metric_DB(2, ib)%X_I, ph => Metric_DB(3, ib)%X_I)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 0, nb(1) + 1
            dB_C(i, j, k, 1) = (State_GVB(i, j+1, k, br_, ib) - State_GVB(i, j-1, k, br_, ib))/(th(j+1) - th(j-1))
            dB_C(i, j, k, 2) = (State_GVB(i, j, k+1, br_, ib) - State_GVB(i, j, k-1, br_, ib))/(ph(k+1) - ph(k-1))
          end do
        end do
      end do
      do k = 1, nb(3)
        do j = 0, nb(2) + 1
          do i = 1, nb(1)
            dB_C(i, j, k, 3) = (r_c(i+1)*State_GVB(i+1, j, k, bt_, ib) - r_c(i-1)*State_GVB(i-1, j, k, bt_, ib))/(r_c(i+1) - r_c(i-1))
            dB_C(i, j, k, 4) = (State_GVB(i, j, k+1, bt_, ib) - State_GVB(i, j, k-1, bt_, ib))/(ph(k+1) - ph(k-1))
          end do
        end do
      end do
      do k = 0, nb(3) + 1
        do j = 1, nb(2)
          do i = 1, nb(1)
            dB_C(i, j, k, 5) = (r_c(i+1)*State_GVB(i+1, j, k, bp_, ib) - r_c(i-1)*State_GVB(i-1, j, k, bp_, ib))/(r_c(i+1) - r_c(i-1))
            dB_C(i, j, k, 6) = (snt(j+1)*State_GVB(i, j+1, k, bp_, ib) - snt(j-1)*State_GVB(i, j-1, k, bp_, ib))/(th(j+1) - th(j-1))
          end do
        end do
      end do
      end associate
    end subroutine bcell_derivs

    ! E: the energy the fluxes carry, with d(rho e) = T d(rho s) + (h - T s) d rho
    ! and d(rho v^2/2) = v.d(rho v) - (v^2/2) d rho:
    !   F_E = T (F_rhos - s F_rho) + (h + Phi) F_rho + v.F_rhov - (v^2/2) F_rho
    ! with the face means of T, s, h + Phi, v (hotta2015: T F_s1 + v.F_rhov;
    ! with B also the magnetic energy B^2/(8 pi) the field's fluxes carry)
    subroutine faces_e(di, dj, dk, F)
      integer, intent(in) :: di, dj, dk
      real(dp), intent(inout) :: F(1-di:nb(1), 1-dj:nb(2), 1-dk:nb(3), ie_)
      integer :: i, j, k
      real(dp) :: Tf, sf, hf, vf(3), fe
      do k = 1 - dk, nb(3)
        do j = 1 - dj, nb(2)
          do i = 1 - di, nb(1)
            Tf = 0.5_dp*(T_G(i, j, k) + T_G(i+di, j+dj, k+dk))
            sf = 0.5_dp*(s0(i) + State_GVB(i, j, k, s1_, ib) + s0(i+di) + State_GVB(i+di, j+dj, k+dk, s1_, ib))
            hf = 0.5_dp*(hP_G(i, j, k) + hP_G(i+di, j+dj, k+dk))
            vf = 0.5_dp*(State_GVB(i, j, k, vr_:vp_, ib) + State_GVB(i+di, j+dj, k+dk, vr_:vp_, ib))
            if (rsst_h2015) then
              fe = Tf*F(i, j, k, s1_) + sum(vf*F(i, j, k, vr_:vp_))
            else
              fe = Tf*(F(i, j, k, s1_) - sf*F(i, j, k, rho1_)) + hf*F(i, j, k, rho1_) + sum(vf*F(i, j, k, vr_:vp_)) &
                  - 0.5_dp*sum(vf**2)*F(i, j, k, rho1_)
              if (use_B) fe = fe + sum(0.5_dp*(State_GVB(i, j, k, br_:bp_, ib) &
                  + State_GVB(i+di, j+dj, k+dk, br_:bp_, ib))*F(i, j, k, br_:bp_))/four_pi
            end if
            F(i, j, k, ie_) = fe
          end do
        end do
      end do
    end subroutine faces_e

  end subroutine diffusion

  ! ================================================================ the energy flux of the diffusion; limiters (energy budget: cz3d_diag_smod)

  ! Energy flux of the diffusion through the radial faces (diagnostic):
  ! with d(rho e) = T d(rho s) + (h - T s) d rho and d(rho v^2/2) = v.d(rho v)
  ! - (v^2/2) d rho, the fluxes F_rho, F_rhov, F_rhos carry
  !   F_E = T (F_rhos - s F_rho) + h1 F_rho + v.F_rhov - (v^2/2) F_rho
  ! with face means of T, s, v and of the enthalpy perturbation h1 = T s1 +
  ! p1/rho (the background enthalpy's share, h0 <F_rho>, as the mean mass
  ! flux term of the enthalpy flux, is left out). The first term is
  ! solar_model's (commented-out) diffusive entropy flux in energy units.
  subroutine radial_energy_flux(ib)
    integer, intent(in) :: ib
    integer :: i, j, k
    real(dp) :: Tf, sf, h1f, vf(3), Fr
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 0, nb(1)
          Tf = 0.5_dp*(T_G(i, j, k) + T_G(i+1, j, k))
          sf = 0.5_dp*(s0(i) + State_GVB(i, j, k, s1_, ib) + s0(i+1) + State_GVB(i+1, j, k, s1_, ib))
          h1f = 0.5_dp*(T_G(i, j, k)*State_GVB(i, j, k, s1_, ib) + p1_G(i, j, k)/rho_G(i, j, k) &
              + T_G(i+1, j, k)*State_GVB(i+1, j, k, s1_, ib) + p1_G(i+1, j, k)/rho_G(i+1, j, k))
          vf = 0.5_dp*(State_GVB(i, j, k, vr_:vp_, ib) + State_GVB(i+1, j, k, vr_:vp_, ib))
          Fr = F1_FV(i, j, k, rho1_)
          if (rsst_h2015) then
            ! separate diffusion: energy rho T s1 + rho v^2/2, F_s the flux of rho s1
            qdif_FB(i, j, k, ib) = Tf*F1_FV(i, j, k, s1_) + sum(vf*F1_FV(i, j, k, vr_:vp_))
            cycle
          end if
          qdif_FB(i, j, k, ib) = Tf*(F1_FV(i, j, k, s1_) - sf*Fr) + h1f*Fr + sum(vf*F1_FV(i, j, k, vr_:vp_)) &
              - 0.5_dp*sum(vf**2)*Fr
        end do
      end do
    end do
  end subroutine radial_energy_flux

  ! Minmod-limited, h-switched jump across the face between q(2) and q(3)
  ! (CZTR's limited_jump)
  pure real(dp) function jump(q) result(jj)
    real(dp), intent(in) :: q(4)
    real(dp) :: du, url, phi
    du = q(3) - q(2)
    if (.not. use_minmod) then
      jj = du
      return
    end if
    url = du - 0.5_dp*(minmod(q(2) - q(1), du) + minmod(du, q(4) - q(3)))
    phi = 0.0_dp
    if (du /= 0.0_dp) phi = max(0.0_dp, 1.0_dp + diffusion_h*(url/du - 1.0_dp))
    jj = phi*url
  end function jump

  ! jump with four scalar arguments, branch-free (merge: both sides are
  ! computed, the division by du = 0 discarded), the same values as jump
  elemental real(dp) function jump4(q1, q2, q3, q4) result(jj)
    real(dp), intent(in) :: q1, q2, q3, q4
    real(dp) :: du, url, phi
    du = q3 - q2
    url = du - 0.5_dp*(minmod2(q2 - q1, du) + minmod2(du, q4 - q3))
    phi = merge(max(0.0_dp, 1.0_dp + diffusion_h*(url/du - 1.0_dp)), 0.0_dp, du /= 0.0_dp)
    jj = merge(phi*url, du, use_minmod)
  end function jump4

  ! jump4 for the field and psi: h = b_diffusion_h (< 0: diffusion_h), the
  ! slopes minmod or MC (b_diff_limiter); the face mismatch u_r - u_l =
  ! du - (slope_2 + slope_3)/2, no flux where it opposes du (no anti-diffusion)
  elemental real(dp) function jump4b(q1, q2, q3, q4) result(jj)
    real(dp), intent(in) :: q1, q2, q3, q4
    real(dp) :: du, url, phi, h, s2, s3
    if (b_diffusion_h < 0.0_dp .and. b_diff_limiter == 'minmod') then
      jj = jump4(q1, q2, q3, q4)
      return
    end if
    h = merge(diffusion_h, b_diffusion_h, b_diffusion_h < 0.0_dp)
    du = q3 - q2
    if (b_diff_limiter == 'mc') then
      s2 = minmod2(0.5_dp*(q3 - q1), minmod2(2.0_dp*(q3 - q2), 2.0_dp*(q2 - q1)))
      s3 = minmod2(0.5_dp*(q4 - q2), minmod2(2.0_dp*(q4 - q3), 2.0_dp*(q3 - q2)))
      url = du - 0.5_dp*(s2 + s3)
    else
      url = du - 0.5_dp*(minmod2(q2 - q1, du) + minmod2(du, q4 - q3))
    end if
    phi = 0.0_dp
    if (du /= 0.0_dp) then
      if (url/du > 0.0_dp) phi = max(0.0_dp, 1.0_dp + h*(url/du - 1.0_dp))
    end if
    jj = merge(phi*url, du, use_minmod)
  end function jump4b

  elemental real(dp) function minmod2(a, b)
    real(dp), intent(in) :: a, b
    minmod2 = merge(0.0_dp, merge(a, b, abs(a) <= abs(b)), a*b <= 0.0_dp)
  end function minmod2

  pure real(dp) function minmod(a, b)
    real(dp), intent(in) :: a, b
    if (a*b <= 0.0_dp) then
      minmod = 0.0_dp
    else if (abs(a) <= abs(b)) then
      minmod = a
    else
      minmod = b
    end if
  end function minmod

  ! diagnostics: cz3d_diag_smod; synthetic EUV and continuum: cz3d_synth_smod

end module cz3d_app_mod
