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
  use arbor
  use cz3d_params_mod
  use pfss_mod, only: pfss_legendre, pfss_radial
  use cz3d_bg_mod, only: k_join, n_r, n_c, bg_view, r_c, r_f, drdxi_c, rho0, s0, p0, lnp0, g_c, z_c, ds0dr, rho1_init, s1_init, &
      v_init, tab, area_f, xi, T0, ddr_c, cs0, vol_c, r_top_cz
  use cz3d_transport_mod, only: transport_column, cool_f, cool_c, grey_i0, tr_view, col_bottom, col_top
  use cztr3d_atm_mod, only: atm_setup, atm_block, atm_floors, atm_dt, atm_tr_heights, S_bgf, bgforce, atm_switch
  use eos_types, only: eos_state
  use eos_thermo, only: eos_get_from_rho_T
  use const_mod, only: Lsun
  use eos_table, only: eos_eval, eos_field, eos_field_n, eos_fields_n, n_field, F_LNP, F_LNT, F_GAMMA1, F_GAMMA3M1, F_LNCP, F_LNCV
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
  real(dp), public :: dt1_now = 0.0_dp            ! the current dt1 (with the ramp)

  ! one-block scratch
  real(dp), allocatable :: rho_G(:,:,:), p1_G(:,:,:), cs_G(:,:,:), W_G(:,:,:), D_C(:,:,:), R_CV(:,:,:,:)
  real(dp), allocatable :: Dr_C(:,:,:), F1_FV(:,:,:,:), F2_FV(:,:,:,:), F3_FV(:,:,:,:), Div_CV(:,:,:,:)
  real(dp), allocatable :: lim_G(:,:,:), T_G(:,:,:), hP_G(:,:,:)
  real(dp), allocatable :: C1_F(:,:,:), C2_F(:,:,:), C3_F(:,:,:)   ! the diffusion's face speeds
  ! use_B scratch: E = v x B (or B itself) with ghosts, a curl, J = curl B, div B, the Alfven speed
  real(dp), allocatable :: B_GV(:,:,:,:), C_CV(:,:,:,:), J_CV(:,:,:,:), divB_C(:,:,:), vA_G(:,:,:)

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
  logical :: lz_ref_set = .false.

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
  ! freeze_cz: each block's state at its first step (cells r_off + i <= k_join restored after every step)
  real(dp), allocatable :: frz_B(:,:,:,:,:)
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

  public :: cz_seed_pfss, max_mach_report, cz_n_var, cz_field_from_hydro, export_state_write, import_open, import_max_level, import_criterion, import_meta, cz_setup, cz_transport_all, cz_log, cz_profile, wall_bc, state_checksum, budget_start, dump_shells_write

contains

  ! ================================================================ setup

  ! The number of state variables (5 hydro; 9 with the magnetic field)
  integer function cz_n_var()
    cz_n_var = merge(nVarB, nVarCz, use_B)
  end function cz_n_var

  ! After arbor_block_sync: per-block arrays, the subgrid factor, scratch
  ! Refinement to refine_h levels in theta and phi everywhere and to
  ! patch_level near the patch centre (arbor_regrid)
  subroutine cz_flag_block(self, iBlock, flags)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iBlock
    integer, intent(out) :: flags(3)
    real(dp) :: x0(3), x(3), e(3,3), la, lo
    integer :: j, k
    flags = 0
    if (Key_IB(KL+2, iBlock) < refine_h) flags(2:3) = 1
    if (Key_IB(KL+2, iBlock) >= patch_level) return
    la = patch_lat_deg*pi/180.0_dp;  lo = patch_lon_deg*pi/180.0_dp
    x0 = [cos(la)*cos(lo), cos(la)*sin(lo), sin(la)]
    do k = 1, nb(3)
      do j = 1, nb(2)
        e = sph_basis(Metric_DB(2, iBlock)%X_I(j), Metric_DB(3, iBlock)%X_I(k))
        x = e(:, 1)
        if (Key_IB(KP, iBlock) == 2) x = yy_cart(x)
        if (dot_product(x, x0) >= cos(patch_radius_deg*pi/180.0_dp)) then
          flags(2:3) = 1
          return
        end if
      end do
    end do
  end subroutine cz_flag_block

  ! The radial view of block ib: the 1D arrays (background, metric, step
  ! factors, subgrid factor, cooling) indexed by the block's local cell 1..nb(1),
  ! and its walls. cz_unbind: the whole column (r_off = 0, both walls).
  subroutine cz_bind(ib)
    integer, intent(in) :: ib
    lev_view = Key_IB(KL+2, ib)
    call column_views(Key_IB(KI+1, ib)*nb(1), nb(1))
    at_bottom = arbor_physical_side(ib, 1, 1)
    at_top = arbor_physical_side(ib, 1, 2)
    col_bottom = at_bottom
    col_top = at_top
    if (use_atmosphere) col_top = Key_IB(KI+1, ib)*nb(1) + nb(1) == k_join   ! (grey and MLT end at the join)
  end subroutine cz_bind

  subroutine cz_unbind()
    lev_view = refine_h
    call column_views(0, n_r)
    at_bottom = .true.
    at_top = .true.
    col_bottom = .true.
    col_top = .true.
  end subroutine cz_unbind

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

  ! After arbor_block_sync or a regrid (blocks may have changed): the
  ! per-block arrays (zeroed: the convective flux restarts at its targets),
  ! the subgrid factor, the scratch arrays (once)
  subroutine cz_setup()
    integer :: i, lev
    real(dp) :: Hp, R
    if (mod(n_r, nb(1)) /= 0) call arbor_abort('cz3d: the radial block size must divide the radial cells')
    if (nb(1) < n_r .and. dt_budget > 0.0_dp) call arbor_abort('cz3d: the energy budget needs whole columns per block')
    if (mr_levels > 1 .and. (dt1 > 0.0_dp .or. dt_budget > 0.0_dp .or. lz_correct_steps > 0 .or. lz_correct_dt > 0.0_dp)) &
        call arbor_abort('cz3d: multirate needs dt1 = 0 and no energy budget or L_z correction')
    call cz_unbind()
    rsst_h2012 = rsst_form == 'hotta2012'
    grey_i0_g = 1                                 ! first cell of the implicit grey solve (whole column)
    if (grey_top_km > 0.0_dp) then
      do i = 1, n_r
        if (-z_c(i) < grey_top_km*1.0e5_dp) then
          grey_i0_g = i
          exit
        end if
      end do
    end if
    if (use_atmosphere .and. mod(k_join, nb(1)) /= 0) call arbor_abort('cztr3d: the join cell must end a radial block '// &
        '(k_join a multiple of the radial block size)')
    if (nb(1) < merge(k_join, n_r, use_atmosphere) .and. grey_i0_g <= merge(k_join, n_r, use_atmosphere) - nb(1)) &
        call arbor_abort('cz3d: with radial blocks the implicit grey layer '// &
        '(grey_top_km) must lie within the top block')
    grey_i0 = grey_i0_g
    rsst_h2015 = rsst_form == 'hotta2015'
    if (use_atmosphere) call atm_setup()
    if (allocated(frz_B)) deallocate(frz_B, frz_have)
    if (freeze_cz) then
      allocate(frz_B(nb(1), nb(2), nb(3), cz_n_var(), nBlockMax), frz_have(nBlockMax))
      frz_have = .false.
    end if
    if (allocated(qrad_FB)) deallocate(qrad_FB, qdif_FB)
    if (allocated(qconv_FB)) deallocate(qconv_FB)
    if (associated(f_mlt_g)) deallocate(f_mlt_g)
    if (.not. qconv_face) then
      allocate(qconv_FB(0:nb(1), nb(2), nb(3), nBlockMax))
      qconv_FB = 0.0_dp
    end if
    allocate(qrad_FB(0:nb(1), nb(2), nb(3), nBlockMax))
    allocate(qdif_FB(0:nb(1), nb(2), nb(3), nBlockMax))
    qrad_FB = 0.0_dp
    qdif_FB = 0.0_dp
    ! the subgrid factor per horizontal level (the cell size halves per level)
    allocate(f_mlt_g(0:n_r, 0:max(refine_h, patch_level)))
    f_mlt_g = 1.0_dp
    do lev = 0, ubound(f_mlt_g, 2)
      dtheta = 0.5_dp*pi/real(n_theta - 2*merge(yy_margin, ng + 1, yy_margin >= 0), dp)/2.0_dp**lev
      if (.not. mlt_switch) cycle
      do i = 0, n_r
        ! background pressure scale height at the face per horizontal cell size
        Hp = 0.5_dp*(p0(i)/(rho0(i)*g_c(i)) + p0(i+1)/(rho0(i+1)*g_c(i+1)))
        R = Hp/(r_f(i)*dtheta)
        f_mlt_g(i, lev) = 1.0_dp/(1.0_dp + (R/mlt_R0)**mlt_p)
      end do
    end do
    lev_view = refine_h
    f_mlt(0:n_r) => f_mlt_g(0:n_r, refine_h)
    if (.not. allocated(tloc)) then
      allocate(tloc(1:n_r))
      tloc = 0.0_dp
    end if
    if (allocated(rho_G)) return
    allocate(tsc_g(1:n_r), rate_r_g(1:n_r))
    tsc => tsc_g
    rate_r => rate_r_g
    tsc = 1.0_dp
    rate_r = 0.0_dp
    allocate(rq_g(-1:n_r+2))
    rq => rq_g
    rq = 1.0_dp
    if (rsst_h2015) rq = 1.0_dp/xi**2
    call background_h_phi()
    allocate(rho0_w(-1:n_r+2), rq_w(-1:n_r+2), s0_w(-1:n_r+2), r_w(-1:n_r+2))
    rho0_w = rho0(-1:n_r+2)
    s0_w = s0(-1:n_r+2)
    r_w = r_c(-1:n_r+2)
    rq_w = rq(-1:n_r+2)
    allocate(rho_G(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng))
    allocate(p1_G, cs_G, W_G, lim_G, T_G, hP_G, mold=rho_G)
    nv = cz_n_var()
    ie_ = nv + 1
    allocate(D_C(nb(1), nb(2), nb(3)), R_CV(nb(1), nb(2), nb(3), nv))
    allocate(Dr_C, mold=D_C)
    ! (variable ie_ = nv + 1: the diffusion's energy flux)
    allocate(F1_FV(0:nb(1), nb(2), nb(3), ie_), F2_FV(nb(1), 0:nb(2), nb(3), ie_), &
        F3_FV(nb(1), nb(2), 0:nb(3), ie_), Div_CV(nb(1), nb(2), nb(3), ie_))
    allocate(C1_F(0:nb(1), nb(2), nb(3)), C2_F(nb(1), 0:nb(2), nb(3)), C3_F(nb(1), nb(2), 0:nb(3)))
    if (use_B) allocate(B_GV(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng, 3), C_CV(nb(1), nb(2), nb(3), 3), &
        J_CV(nb(1), nb(2), nb(3), 3), divB_C(nb(1), nb(2), nb(3)), vA_G(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng))
  end subroutine cz_setup

  ! Initial state: the 1D state in every column, plus an optional random s1
  ! (from each cell's global index, so independent of the number of ranks)
  subroutine cz_init_block(self, iBlock)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iBlock
    integer :: i, j, k
    call cz_bind(iBlock)
    State_GVB(:, :, :, :, iBlock) = 0.0_dp
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, nb(1)
          State_GVB(i, j, k, rho1_, iBlock) = rho1_init(i)/rq(i)
          State_GVB(i, j, k, vr_, iBlock) = v_init(i)
          State_GVB(i, j, k, s1_, iBlock) = s1_init(i)
          if (ds1_random > 0.0_dp .and. (.not. use_atmosphere .or. r_off + i <= k_join)) &
              State_GVB(i, j, k, s1_, iBlock) = s1_init(i) + ds1_random &
              *urand(Key_IB(KP, iBlock), r_off + i, Key_IB(KI+2, iBlock)*nb(2) + j, Key_IB(KI+3, iBlock)*nb(3) + k)
        end do
      end do
    end do
    if (len_trim(import_dir) > 0) call import_block(iBlock)
    if (use_B .and. b_seed == 'toroidal') call seed_toroidal(iBlock)
    if (use_B .and. b_seed == 'dipole') call seed_dipole(iBlock)
    call cz_unbind()
  end subroutine cz_init_block

  ! Toroidal seed field about the global z axis, added to the block's field:
  ! B = b_seed_G w(r) sin(theta_g) e_phi_g = b_seed_G w(r) (-y, x, 0) for the
  ! unit position (x, y, z) in the global (Yin) frame; divergence free. w = 1
  ! below b_seed_r1_rsun, 0 above b_seed_r2_rsun, sin^2 between.
  subroutine seed_toroidal(iBlock)
    use const_mod, only: Rsun
    integer, intent(in) :: iBlock
    real(dp) :: e(3,3), x(3), bg(3), bl(3), w, r1, r2
    integer :: i, j, k
    r1 = b_seed_r1_rsun*Rsun
    r2 = b_seed_r2_rsun*Rsun
    do k = 1, nb(3)
      do j = 1, nb(2)
        e = sph_basis(Metric_DB(2, iBlock)%X_I(j), Metric_DB(3, iBlock)%X_I(k))
        x = e(:, 1)
        if (Key_IB(KP, iBlock) == 2) x = yy_cart(x)          ! the position in the global frame
        bg = [-x(2), x(1), 0.0_dp]
        if (Key_IB(KP, iBlock) == 2) bg = yy_cart(bg)        ! back to this patch's frame
        bl = matmul(transpose(e), bg)                         ! (B_r, B_theta, B_phi) per unit amplitude
        do i = 1, nb(1)
          if (r_c(i) <= r1) then
            w = 1.0_dp
          else if (r_c(i) >= r2) then
            w = 0.0_dp
          else
            w = cos(0.5_dp*pi*(r_c(i) - r1)/(r2 - r1))**2
          end if
          State_GVB(i, j, k, br_:bp_, iBlock) = State_GVB(i, j, k, br_:bp_, iBlock) + b_seed_G*w*bl
        end do
      end do
    end do
  end subroutine seed_toroidal

  ! Dipole seed (b_seed = 'dipole'): the l = 1 potential field about the
  ! global z axis with a source surface at r_ss = b_seed_rss_rsun, Phi = A
  ! cos(theta) [(r_ss/r)^2 - r/r_ss]: B_r = A cos(theta) (2 r_ss^2/r^3 +
  ! 1/r_ss), B_theta = A sin(theta) (r_ss^2/r^3 - 1/r_ss) (0 at r_ss); above r_ss
  ! radial, B_r(r_ss, theta) (r_ss/r)^2; A from b_seed_G = B_r at the pole at
  ! r_top_cz. Added to the block's field (theta: the global colatitude)
  subroutine seed_dipole(iBlock)
    use const_mod, only: Rsun
    integer, intent(in) :: iBlock
    real(dp) :: e(3,3), x(3), bg(3), bl(3), rss, A, ct, st, r, br, bth, rr, thg(3)
    integer :: i, j, k
    rss = b_seed_rss_rsun*Rsun
    A = b_seed_G/(2.0_dp*rss**2/r_top_cz**3 + 1.0_dp/rss)
    do k = 1, nb(3)
      do j = 1, nb(2)
        e = sph_basis(Metric_DB(2, iBlock)%X_I(j), Metric_DB(3, iBlock)%X_I(k))
        x = e(:, 1)
        if (Key_IB(KP, iBlock) == 2) x = yy_cart(x)          ! the position in the global frame
        ct = x(3)
        st = sqrt(max(1.0_dp - ct*ct, 0.0_dp))
        thg = [0.0_dp, 0.0_dp, -1.0_dp]                       ! e_theta of the global frame (at the poles: any)
        if (st > 1.0e-12_dp) thg = [ct*x(1)/st, ct*x(2)/st, -st]
        do i = 1, nb(1)
          r = r_c(i)
          if (r <= rss) then
            br = A*ct*(2.0_dp*rss**2/r**3 + 1.0_dp/rss)
            bth = A*st*(rss**2/r**3 - 1.0_dp/rss)
          else
            rr = rss/r
            br = A*ct*3.0_dp/rss*rr**2
            bth = 0.0_dp
          end if
          bg = br*x + bth*thg                                 ! global Cartesian
          if (Key_IB(KP, iBlock) == 2) bg = yy_cart(bg)      ! to this patch's frame
          bl = matmul(transpose(e), bg)
          State_GVB(i, j, k, br_:bp_, iBlock) = State_GVB(i, j, k, br_:bp_, iBlock) + bl
        end do
      end do
    end do
  end subroutine seed_dipole

  ! A hydrodynamic restart read into an MHD run (arbor_restart_read with
  ! n_var_file < nv): B and psi came back zero; add the seed field to every
  ! block (init_block is not called on a restart)
  subroutine cz_field_from_hydro()
    integer :: jb
    character(len=200) :: msg
    if (.not. use_B) return
    do jb = 1, nBlock
      call cz_bind(iBlockList_I(jb))
      State_GVB(:, :, :, br_:psi_, iBlockList_I(jb)) = 0.0_dp
      if (b_seed == 'toroidal') call seed_toroidal(iBlockList_I(jb))
      if (b_seed == 'dipole') call seed_dipole(iBlockList_I(jb))
    end do
    call cz_unbind()
    write(msg, '(A,A,A,F8.2,A)') '  magnetic field added to a hydrodynamic restart: seed ', trim(b_seed), ' (', b_seed_G, ' G)'
    call arbor_print(msg)
  end subroutine cz_field_from_hydro

  ! b_seed = 'pfss': the field above the join (cells > k_join) replaced by the
  ! PFSS field of the join cell's B_r (pfss_mod; R0 = that cell's centre, the
  ! source surface at b_seed_rss_rsun, radial above it). The coefficients
  ! come from a midpoint quadrature over the join shell (solid-angle weights,
  ! 1/2 where both Yin-Yang cores cover a point, as dump_shells_write), refined
  ! pfss_iter times on the residual, summed over the ranks; the monopole is
  ! dropped. Each block column then sums the harmonics once (O(L^2)) and each
  ! cell its l series (O(L)), with l cut where (R0/r)^l < 1e-10. psi = 0
  ! above the join. pfss_check: the CZ field set to the analytic dipole
  ! first, and the result compared with it.
  subroutine cz_seed_pfss()
    use const_mod, only: Rsun
    real(dp), allocatable :: c(:,:), s(:,:), dc(:,:), ds(:,:), p(:,:), dpt(:,:), Sl(:), Tl(:), Ul(:), &
        fr(:,:), gh(:,:), frs(:), ghs(:), cmp(:), smp(:), Bsave(:,:,:,:), phi_B(:,:,:,:)
    real(dp) :: R0, rss, dth, x(3), e(3,3), eth(3), eph(3), bg(3), ct, st, ph, w, b, syn, a, &
        st_(5), br, bt, bp, x0, err(3), wt, t0
    integer :: lmx, lcut, it, jb, ib, i, j, k, l, m, k_ss
    logical :: want_phi
    character(len=200) :: msg
    if (pfss_check) then
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        State_GVB(:, :, :, br_:psi_, ib) = 0.0_dp
        call seed_dipole(ib)
      end do
      call cz_unbind()
    end if
    call cz_unbind()
    t0 = arbor_wtime()
    R0 = r_c(k_join)                       ! (unbound view: the global grid)
    rss = b_seed_rss_rsun*Rsun
    k_ss = 0
    if (b_seed == 'dpf') then                          ! the source surface at the nearest radial-block face
      k_ss = k_join + nb(1)
      do i = k_join + nb(1), n_r, nb(1)
        if (abs(r_f(i) - rss) < abs(r_f(k_ss) - rss)) k_ss = i
      end do
      rss = r_f(k_ss)
    end if
    dth = huge(1.0_dp)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      dth = min(dth, minval(Metric_DB(2, ib)%Dx_I(1:nb(2))))
    end do
    dth = arbor_min(dth)
    lmx = pfss_lmax
    if (lmx <= 0) lmx = nint(pi/dth)
    allocate(c(0:lmx, 0:lmx), s(0:lmx, 0:lmx), dc(0:lmx, 0:lmx), ds(0:lmx, 0:lmx), p(0:lmx, 0:lmx), dpt(0:lmx, 0:lmx), &
        Sl(0:lmx), Tl(0:lmx), Ul(0:lmx), fr(0:lmx, nb(1)), gh(0:lmx, nb(1)), frs(0:lmx), ghs(0:lmx), cmp(0:lmx), smp(0:lmx))
    c = 0.0_dp
    s = 0.0_dp
    ! --- the fit: pass it synthesises the current coefficients at the join
    ! shell; the residual's quadrature corrects them (not after the last pass)
    do it = 0, pfss_iter
      dc = 0.0_dp
      ds = 0.0_dp
      st_ = 0.0_dp                          ! sum w, sum w B, sum w B^2, sum w res^2, max |res|
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        if (Key_IB(KI+1, ib)*nb(1) + nb(1) /= k_join) cycle
        do k = 1, nb(3)
          do j = 1, nb(2)
            w = shell_weight(ib, j, k)
            if (w <= 0.0_dp) cycle
            call global_angles(ib, j, k, x, ct, st, ph)
            call pfss_legendre(lmx, ct, st, p)
            do m = 0, lmx
              cmp(m) = cos(m*ph)
              smp(m) = sin(m*ph)
            end do
            b = State_GVB(nb(1), j, k, br_, ib)
            syn = 0.0_dp
            do m = 0, lmx
              do l = m, lmx
                syn = syn + p(l, m)*(c(l, m)*cmp(m) + s(l, m)*smp(m))
              end do
            end do
            st_ = [st_(1) + w, st_(2) + w*b, st_(3) + w*b*b, st_(4) + w*(b - syn)**2, max(st_(5), abs(b - syn))]
            if (it == pfss_iter) cycle
            do m = 0, lmx
              a = merge(1.0_dp, 2.0_dp, m == 0)*w*(b - syn)
              dc(m:lmx, m) = dc(m:lmx, m) + a*cmp(m)*p(m:lmx, m)
              ds(m:lmx, m) = ds(m:lmx, m) + a*smp(m)*p(m:lmx, m)
            end do
          end do
        end do
      end do
      call arbor_sum_inplace(st_(1:4))
      st_(5) = arbor_max(st_(5))
      write(msg, '(A,I0,A,I0,A,F9.5,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A)') '  pfss: l_max ', lmx, ', pass ', it, &
          ': sum of weights/4pi ', st_(1)/(4.0_dp*pi), ', <B_r> ', st_(2)/st_(1), ' G, rms B_r ', sqrt(st_(3)/st_(1)), &
          ' G, rms residual ', sqrt(st_(4)/st_(1)), ', max ', st_(5), ' G'
      call arbor_print(msg)
      if (it == pfss_iter) exit
      call arbor_sum_inplace(dc)
      call arbor_sum_inplace(ds)
      c = c + dc
      s = s + ds
    end do
    c(0, 0) = 0.0_dp                          ! (no monopole)
    ! --- the field above the join
    want_phi = dpf_check .or. b_seed == 'dpf'
    if (want_phi) then
      allocate(phi_B(nb(1), nb(2), nb(3), nBlockMax))
      phi_B = 0.0_dp
    end if
    err = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (r_off + nb(1) < k_join .or. (r_off + nb(1) == k_join .and. .not. want_phi)) cycle
      x0 = R0/min(r_f(0), rss)
      lcut = lmx
      if (x0 < 1.0_dp - 1.0e-12_dp) lcut = min(lmx, max(1, ceiling(log(1.0e-10_dp)/log(x0))))
      do i = 1, nb(1)
        call pfss_radial(lcut, min(r_c(i), rss), R0, rss, fr(0:lcut, i), gh(0:lcut, i))
      end do
      call pfss_radial(lcut, rss, R0, rss, frs(0:lcut), ghs(0:lcut))
      if (pfss_check) Bsave = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:bp_, ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          call global_angles(ib, j, k, x, ct, st, ph)
          call pfss_legendre(lcut, ct, st, p(0:lcut, 0:lcut), dpt(0:lcut, 0:lcut))
          Sl(0:lcut) = 0.0_dp
          Tl(0:lcut) = 0.0_dp
          Ul(0:lcut) = 0.0_dp
          do m = 0, lcut
            cmp(m) = cos(m*ph)
            smp(m) = sin(m*ph)
            do l = m, lcut
              a = c(l, m)*cmp(m) + s(l, m)*smp(m)
              Sl(l) = Sl(l) + p(l, m)*a
              Tl(l) = Tl(l) + dpt(l, m)*a
              Ul(l) = Ul(l) + p(l, m)*m*(s(l, m)*cmp(m) - c(l, m)*smp(m))
            end do
          end do
          Ul(0:lcut) = Ul(0:lcut)/max(st, 1.0e-12_dp)
          eth = [ct*cos(ph), ct*sin(ph), -st]
          eph = [-sin(ph), cos(ph), 0.0_dp]
          e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
          do i = 1, nb(1)
            if (want_phi) then                         ! the potential: Phi = -r sum gh_l S_l
              phi_B(i, j, k, ib) = 0.0_dp
              if (r_c(i) <= rss) phi_B(i, j, k, ib) = -r_c(i)*sum(gh(0:lcut, i)*Sl(0:lcut))
              if (r_off + i <= k_join) cycle           ! (the CZ cells: Phi only)
            end if
            if (r_c(i) <= rss) then
              br = sum(fr(0:lcut, i)*Sl(0:lcut))
              bt = sum(gh(0:lcut, i)*Tl(0:lcut))
              bp = sum(gh(0:lcut, i)*Ul(0:lcut))
            else
              br = sum(frs(0:lcut)*Sl(0:lcut))*(rss/r_c(i))**2
              bt = 0.0_dp
              bp = 0.0_dp
            end if
            bg = br*x + bt*eth + bp*eph                 ! global Cartesian
            if (Key_IB(KP, ib) == 2) bg = yy_cart(bg)  ! to this patch's frame
            State_GVB(i, j, k, br_:bp_, ib) = matmul(transpose(e), bg)
            State_GVB(i, j, k, psi_, ib) = 0.0_dp
            if (pfss_check) then
              wt = sum((State_GVB(i, j, k, br_:bp_, ib) - Bsave(i, j, k, :))**2)
              err = [err(1) + wt, err(2) + sum(Bsave(i, j, k, :)**2), max(err(3), sqrt(wt))]
            end if
          end do
        end do
      end do
    end do
    call cz_unbind()
    if (pfss_check) then
      call arbor_sum_inplace(err(1:2))
      err(3) = arbor_max(err(3))
      write(msg, '(A,ES10.3,A,ES10.3,A)') '  pfss_check: |B - B_dipole| above the join: rms/rms|B| ', &
          sqrt(err(1)/max(err(2), tiny(1.0_dp))), ', max ', err(3)/b_seed_G, ' b_seed_G'
      call arbor_print(msg)
    end if
    write(msg, '(A,ES10.3,A,F8.2,A)') '  pfss: field above the join set (source surface ', rss/Rsun, ' R_sun), ', &
        arbor_wtime() - t0, ' s'
    call arbor_print(msg)
    if (dpf_check) call dpf_check_run(phi_B, rss)
    if (b_seed == 'dpf') then
      call dpf_solve(phi_B, k_ss)
      if (mf_iter > 0) call mf_relax(r_f(k_ss))
      if (pfss_check) then                           ! the solved field against the analytic dipole (below 0.95 r_ss)
        err = 0.0_dp
        do jb = 1, nBlock
          ib = iBlockList_I(jb)
          call cz_bind(ib)
          if (r_off + nb(1) <= k_join) cycle
          Bsave = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:bp_, ib)
          State_GVB(:, :, :, br_:bp_, ib) = 0.0_dp
          call seed_dipole(ib)
          do k = 1, nb(3)
            do j = 1, nb(2)
              do i = 1, nb(1)
                if (r_c(i) > 0.95_dp*rss) cycle
                wt = sum((Bsave(i, j, k, :) - State_GVB(i, j, k, br_:bp_, ib))**2)
                err = [err(1) + wt, err(2) + sum(State_GVB(i, j, k, br_:bp_, ib)**2), max(err(3), sqrt(wt))]
              end do
            end do
          end do
          State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:bp_, ib) = Bsave
        end do
        call cz_unbind()
        call arbor_sum_inplace(err(1:2))
        err(3) = arbor_max(err(3))
        write(msg, '(A,ES10.3,A,ES10.3,A)') '  pfss_check (dpf): |B - B_dipole| below 0.95 r_ss: rms/rms|B| ', &
            sqrt(err(1)/max(err(2), tiny(1.0_dp))), ', max ', err(3)/b_seed_G, ' b_seed_G'
        call arbor_print(msg)
      end if
    end if
  contains
    ! the cell's solid angle, halved where both cores cover the point, 0 outside the core
    real(dp) function shell_weight(ib, j, k) result(w)
      integer, intent(in) :: ib, j, k
      real(dp) :: th, ph_, th2, ph2
      th = Metric_DB(2, ib)%X_I(j)
      ph_ = Metric_DB(3, ib)%X_I(k)
      w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
      if (.not. core(th, ph_)) then
        w = 0.0_dp
      else
        call yy_angles(th, ph_, th2, ph2)
        if (core(th2, ph2)) w = 0.5_dp*w
      end if
    end function shell_weight
    logical function core(t, p_)
      real(dp), intent(in) :: t, p_
      core = t >= 0.25_dp*pi .and. t <= 0.75_dp*pi .and. abs(p_) <= 0.75_dp*pi
    end function core
    ! the cell centre's global unit vector, cos/sin of the global colatitude and the longitude
    subroutine global_angles(ib, j, k, x, ct, st, ph)
      integer, intent(in) :: ib, j, k
      real(dp), intent(out) :: x(3), ct, st, ph
      real(dp) :: e(3,3)
      e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
      x = e(:, 1)
      if (Key_IB(KP, ib) == 2) x = yy_cart(x)
      ct = max(-1.0_dp, min(1.0_dp, x(3)))
      st = sqrt(max(1.0_dp - ct*ct, 0.0_dp))
      ph = atan2(x(2), x(1))
    end subroutine global_angles
  end subroutine cz_seed_pfss

  ! dpf_check: the curl of the analytic PFSS field (now in the state) and of
  ! B = -grad Phi from the code's own stencils (Phi = the PFSS potential,
  ! phi_B), both with filled ghosts, measured by the code's curl as
  ! q = |curl B| r dtheta/|B| (B^2-weighted rms, and the max) in atmosphere
  ! cells whose stencils stay above the join and below 0.95 r_ss, split into
  ! patch-interior and patch-edge cells (within 2.5 cells of the Yin-Yang
  ! patch boundary: their stencils reach interpolated ghosts) and three
  ! heights (z < 5 Mm, 5-50 Mm, > 50 Mm). Also rms(div B dx)/rms|B| of both.
  ! The state's field is restored afterwards (the run keeps the analytic PFSS).
  subroutine dpf_check_run(phi_B, rss)
    real(dp), intent(in) :: phi_B(:,:,:,:), rss
    real(dp), allocatable :: Ban(:,:,:,:,:)
    real(dp) :: acc(4, 2, 3, 2), th, ph, dth, dph, q, b2, z, c2, d2, rj
    integer :: jb, ib, i, j, k, m, ireg, ih
    character(len=200) :: msg
    character(len=8), parameter :: fname(2) = ['PFSS    ', '-grad Ph']
    character(len=8), parameter :: rname(2) = ['interior', 'edge    ']
    character(len=10), parameter :: hname(3) = ['z<5Mm     ', '5-50Mm    ', 'z>50Mm    ']
    call cz_unbind()
    rj = r_f(k_join)                                     ! the join face (unbound view)
    allocate(Ban(nb(1), nb(2), nb(3), 3, nBlockMax))
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      Ban(:, :, :, :, ib) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:bp_, ib)
    end do
    acc = 0.0_dp            ! (sum B^2, sum B^2 q^2, max q, sum (div B dx)^2) x region x height x field
    do m = 1, 2
      if (m == 2) then
        ! Phi into psi (above the join and the top CZ block), then B = -grad Phi above the join
        do jb = 1, nBlock
          ib = iBlockList_I(jb)
          State_GVB(1:nb(1), 1:nb(2), 1:nb(3), psi_, ib) = phi_B(:, :, :, ib)
        end do
        call arbor_fill_ghosts(State_GVB)
        do jb = 1, nBlock
          ib = iBlockList_I(jb)
          call cz_bind(ib)
          if (r_off + nb(1) <= k_join) cycle
          call ddr_block(State_GVB(:, :, :, psi_, ib), Dr_C)
          associate (th_ => Metric_DB(2, ib)%X_I)
            do i = max(1, k_join - r_off + 1), nb(1)
              State_GVB(i, 1:nb(2), 1:nb(3), br_, ib) = -Dr_C(i, :, :)
            end do
            call arbor_deriv_cell(State_GVB, psi_, 2, ib, D_C)
            do i = max(1, k_join - r_off + 1), nb(1)
              State_GVB(i, 1:nb(2), 1:nb(3), bt_, ib) = -D_C(i, :, :)/r_c(i)
            end do
            call arbor_deriv_cell(State_GVB, psi_, 3, ib, D_C)
            do k = 1, nb(3)
              do j = 1, nb(2)
                do i = max(1, k_join - r_off + 1), nb(1)
                  State_GVB(i, j, k, bp_, ib) = -D_C(i, j, k)/(r_c(i)*sin(th_(j)))
                end do
              end do
            end do
          end associate
        end do
        call cz_unbind()
      end if
      call arbor_fill_ghosts(State_GVB)
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        if (r_off + nb(1) <= k_join) cycle
        B_GV = State_GVB(:, :, :, br_:bp_, ib)
        call curl_block(B_GV, ib, J_CV)
        call div_block(State_GVB(:, :, :, br_:bp_, ib), ib, divB_C)
        dth = Metric_DB(2, ib)%Dx_I(1)
        dph = Metric_DB(3, ib)%Dx_I(1)
        do k = 1, nb(3)
          do j = 1, nb(2)
            th = Metric_DB(2, ib)%X_I(j)
            ph = Metric_DB(3, ib)%X_I(k)
            ireg = 1
            if (th < 0.25_dp*pi + 2.5_dp*dth .or. th > 0.75_dp*pi - 2.5_dp*dth .or. &
                abs(ph) > 0.75_dp*pi - 2.5_dp*dph) ireg = 2
            do i = 1, nb(1)
              if (r_off + i < k_join + 3 .or. r_c(i) > 0.95_dp*rss) cycle
              z = r_c(i) - rj
              ih = 1
              if (z > 5.0e8_dp) ih = 2
              if (z > 5.0e9_dp) ih = 3
              b2 = sum(State_GVB(i, j, k, br_:bp_, ib)**2)
              if (b2 <= 0.0_dp) cycle
              c2 = sum(J_CV(i, j, k, :)**2)*(r_c(i)*dth)**2
              d2 = (divB_C(i, j, k)*r_c(i)*dth)**2
              q = sqrt(c2/b2)
              acc(:, ireg, ih, m) = [acc(1, ireg, ih, m) + b2, acc(2, ireg, ih, m) + c2, &
                  max(acc(3, ireg, ih, m), q), acc(4, ireg, ih, m) + d2]
            end do
          end do
        end do
      end do
      call cz_unbind()
    end do
    do m = 1, 2
      do ih = 1, 3
        do ireg = 1, 2
          acc(1, ireg, ih, m) = arbor_sum(acc(1, ireg, ih, m))
          acc(2, ireg, ih, m) = arbor_sum(acc(2, ireg, ih, m))
          acc(3, ireg, ih, m) = arbor_max(acc(3, ireg, ih, m))
          acc(4, ireg, ih, m) = arbor_sum(acc(4, ireg, ih, m))
          write(msg, '(A,A,1X,A,1X,A,A,ES10.3,A,ES10.3,A,ES10.3)') '  dpf_check ', fname(m), rname(ireg), hname(ih), &
              ': rms |curl B| r dth/|B| ', sqrt(acc(2, ireg, ih, m)/max(acc(1, ireg, ih, m), tiny(1.0_dp))), &
              '  max ', acc(3, ireg, ih, m), '  rms div B r dth/|B| ', sqrt(acc(4, ireg, ih, m)/max(acc(1, ireg, ih, m), tiny(1.0_dp)))
          call arbor_print(msg)
        end do
      end do
    end do
    ! the analytic field back (psi 0 above the join, as cz_seed_pfss sets it)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (r_off + nb(1) <= k_join) cycle
      do i = max(1, k_join - r_off + 1), nb(1)
        State_GVB(i, 1:nb(2), 1:nb(3), br_:bp_, ib) = Ban(i, :, :, :, ib)
        State_GVB(i, 1:nb(2), 1:nb(3), psi_, ib) = 0.0_dp
      end do
    end do
    call cz_unbind()
    call arbor_fill_ghosts(State_GVB)
  end subroutine dpf_check_run

  ! b_seed = 'dpf': the discrete potential field above the join. Unknown Phi
  ! in the cells k_join+1..k_ss (whole radial blocks; k_ss a block face, the
  ! source surface). Compact finite-volume Laplacian (V times it, symmetric
  ! within a patch): sum over the faces of c_f (Phi_nb - Phi), c_f = area/
  ! distance; at the join face the flux of the join cell's B_r (r^2 B_r kept)
  ! on the right-hand side, at the k_ss face Phi = 0 (odd ghost). The
  ! neighbours come from ghost fills (Yin-Yang: linear interpolation,
  ! LIMITER_NONE, during the solve). BiCGSTAB, right-preconditioned by a radial
  ! tridiagonal solve per column within each block, started from the PFSS
  ! potential phi_B. Then B = -grad Phi with the code's stencils in the cells
  ! k_join+1..k_ss (Phi below the join extended with the CZ's B_r, odd above
  ! k_ss), and B_r ~ r^-2, B_h = 0 above k_ss; psi = 0 above the join.
  subroutine dpf_solve(phi_B, k_ss)
    real(dp), intent(in) :: phi_B(:,:,:,:)
    integer, intent(in) :: k_ss
    real(dp), allocatable :: X(:,:,:,:,:), Cf(:,:,:,:,:), xv(:,:,:,:), rv(:,:,:,:), rh(:,:,:,:), pv(:,:,:,:), &
        vv(:,:,:,:), sv(:,:,:,:), tv(:,:,:,:), bv(:,:,:,:), ph_(:,:,:,:), sh(:,:,:,:), yc(:,:,:,:)
    logical, allocatable :: edge(:,:,:)
    real(dp) :: rho_, rho_old, alpha, omega, beta, bnorm, rnorm, t0, dO, drr, brf
    integer :: it, jb, ib, i, j, k, lim_save, layer
    character(len=200) :: msg
    t0 = arbor_wtime()
    lim_save = iLimiter
    iLimiter = LIMITER_NONE
    allocate(X(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng, 1, nBlockMax), Cf(nb(1), nb(2), nb(3), 6, nBlockMax))
    allocate(xv(nb(1), nb(2), nb(3), nBlockMax))
    allocate(rv, rh, pv, vv, sv, tv, bv, ph_, sh, yc, mold=xv)
    allocate(edge(nb(2), nb(3), nBlockMax))
    edge = .false.
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      do k = 1, nb(3)
        do j = 1, nb(2)
          associate (th => Metric_DB(2, ib)%X_I(j), ph => Metric_DB(3, ib)%X_I(k), &
              dth => Metric_DB(2, ib)%Dx_I(j), dph => Metric_DB(3, ib)%Dx_I(k))
            edge(j, k, ib) = dpf_operator == 'wide' .and. dpf_edge_cells > 0.0_dp .and. &
                (th < 0.25_dp*pi + dpf_edge_cells*dth .or. th > 0.75_dp*pi - dpf_edge_cells*dth &
                .or. abs(ph) > 0.75_dp*pi - dpf_edge_cells*dph)
          end associate
        end do
      end do
    end do
    X = 0.0_dp;  Cf = 0.0_dp
    xv = 0.0_dp;  rv = 0.0_dp;  rh = 0.0_dp;  pv = 0.0_dp;  vv = 0.0_dp;  sv = 0.0_dp;  tv = 0.0_dp;  bv = 0.0_dp
    ph_ = 0.0_dp;  sh = 0.0_dp
    ! the face coefficients (1: r below, 2: r above, 3: theta below, 4: theta above, 5/6: phi) and b
    call arbor_fill_ghosts(State_GVB)                    ! (the join cell's B_r in the ghosts of the bottom blocks)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (.not. in_solve(ib)) cycle
      associate (th => Metric_DB(2, ib)%X_I, thf => Metric_DB(2, ib)%X_F, dth => Metric_DB(2, ib)%Dx_I, &
          ph => Metric_DB(3, ib)%X_I, dph => Metric_DB(3, ib)%Dx_I, vt => Metric_DB(2, ib)%Vol_I, &
          vp => Metric_DB(3, ib)%Vol_I)
        do k = 1, nb(3)
          do j = 1, nb(2)
            dO = vt(j)*vp(k)
            do i = 1, nb(1)
              drr = r_f(i) - r_f(i-1)
              if (r_off + i > k_join + 1) Cf(i, j, k, 1, ib) = r_f(i-1)**2*dO/(r_c(i) - r_c(i-1))
              Cf(i, j, k, 2, ib) = r_f(i)**2*dO/(r_c(i+1) - r_c(i))
              if (r_off + i == k_ss) Cf(i, j, k, 2, ib) = 2.0_dp*Cf(i, j, k, 2, ib)     ! (Phi = 0 at the face: odd ghost)
              Cf(i, j, k, 3, ib) = drr*sin(thf(j-1))*dph(k)/(th(j) - th(j-1))
              Cf(i, j, k, 4, ib) = drr*sin(thf(j))*dph(k)/(th(j+1) - th(j))
              Cf(i, j, k, 5, ib) = drr*dth(j)/(sin(th(j))*(ph(k) - ph(k-1)))
              Cf(i, j, k, 6, ib) = drr*dth(j)/(sin(th(j))*(ph(k+1) - ph(k)))
              if (r_off + i == k_join + 1) bv(i, j, k, ib) = -dO*r_c(i-1)**2*State_GVB(i-1, j, k, br_, ib)
            end do
          end do
        end do
      end associate
      xv(:, :, :, ib) = phi_B(:, :, :, ib)
    end do
    call cz_unbind()
    if (dpf_operator == 'wide') then                 ! b = V div B(Phi = 0): the CZ's field below the join
      sv = 0.0_dp
      call wide_op(sv, tv, .true.)
      do jb = 1, nBlock                                ! (bv held the compact operator's b)
        ib = iBlockList_I(jb)
        do k = 1, nb(3)
          do j = 1, nb(2)
            if (.not. edge(j, k, ib)) bv(:, j, k, ib) = -tv(:, j, k, ib) + max(dpf_eps, 0.0_dp)*bv(:, j, k, ib)
          end do
        end do
      end do
      tv = 0.0_dp
    end if
    ! BiCGSTAB (right preconditioned): r = b - A x
    call apply_A(xv, vv)
    rv = bv - vv
    rh = rv
    bnorm = sqrt(dot(bv, bv))
    rnorm = sqrt(dot(rv, rv))
    write(msg, '(A,I0,A,ES10.3,A,ES10.3)') '  dpf: source surface at face ', k_ss, ', |b| ', bnorm, &
        ', relative residual of the PFSS start ', rnorm/bnorm
    call arbor_print(msg)
    rho_old = 1.0_dp;  alpha = 1.0_dp;  omega = 1.0_dp
    pv = 0.0_dp;  vv = 0.0_dp
    do it = 1, dpf_maxit
      rho_ = dot(rh, rv)
      beta = (rho_/rho_old)*(alpha/omega)
      pv = rv + beta*(pv - omega*vv)
      call precond(pv, ph_)
      call apply_A(ph_, vv)
      alpha = rho_/dot(rh, vv)
      sv = rv - alpha*vv
      call precond(sv, sh)
      call apply_A(sh, tv)
      omega = dot(tv, sv)/dot(tv, tv)
      xv = xv + alpha*ph_ + omega*sh
      rv = sv - omega*tv
      rho_old = rho_
      rnorm = sqrt(dot(rv, rv))
      if (mod(it, 50) == 0 .or. rnorm <= dpf_tol*bnorm) then
        write(msg, '(A,I6,A,ES10.3)') '  dpf: iteration ', it, ', relative residual ', rnorm/bnorm
        call arbor_print(msg)
      end if
      if (rnorm <= dpf_tol*bnorm) exit
    end do
    ! B = -grad Phi (Phi in psi; its ghosts below the join and above k_ss set here)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (in_solve(ib)) State_GVB(1:nb(1), 1:nb(2), 1:nb(3), psi_, ib) = xv(:, :, :, ib)
    end do
    call arbor_fill_ghosts(State_GVB)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (.not. in_solve(ib)) cycle
      if (r_off == k_join) then                        ! Phi continued with the CZ's B_r (ghosts 0, -1: cells k_join, k_join - 1)
        State_GVB(0, :, :, psi_, ib) = State_GVB(1, :, :, psi_, ib) &
            + State_GVB(0, :, :, br_, ib)*(r_c(0)/r_f(0))**2*(r_c(1) - r_c(0))
        State_GVB(-1, :, :, psi_, ib) = State_GVB(0, :, :, psi_, ib) + State_GVB(0, :, :, br_, ib)*(r_c(0) - r_c(-1))
      end if
      if (r_off + nb(1) == k_ss) then                  ! odd about the k_ss face
        State_GVB(nb(1)+1, :, :, psi_, ib) = -State_GVB(nb(1), :, :, psi_, ib)
        State_GVB(nb(1)+2, :, :, psi_, ib) = -State_GVB(nb(1)-1, :, :, psi_, ib)
      end if
      call ddr_block(State_GVB(:, :, :, psi_, ib), Dr_C)
      State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_, ib) = -Dr_C
      call arbor_deriv_cell(State_GVB, psi_, 2, ib, D_C)
      do i = 1, nb(1)
        State_GVB(i, 1:nb(2), 1:nb(3), bt_, ib) = -D_C(i, :, :)/r_c(i)
      end do
      call arbor_deriv_cell(State_GVB, psi_, 3, ib, D_C)
      associate (th => Metric_DB(2, ib)%X_I)
        do j = 1, nb(2)
          do i = 1, nb(1)
            State_GVB(i, j, 1:nb(3), bp_, ib) = -D_C(i, j, :)/(r_c(i)*sin(th(j)))
          end do
        end do
      end associate
    end do
    call cz_unbind()
    ! above k_ss: radial, block layer by layer (each from the ghosts of the one below)
    do layer = 0, (n_r - k_ss)/nb(1) - 1
      call arbor_fill_ghosts(State_GVB)
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        if (r_off /= k_ss + layer*nb(1)) cycle
        do k = 1, nb(3)
          do j = 1, nb(2)
            if (layer == 0) then                       ! B_r at the k_ss face from Phi (0 there)
              brf = 2.0_dp*State_GVB(0, j, k, psi_, ib)/(r_c(1) - r_c(0))*(r_f(0)/r_c(1))**0
              brf = brf*r_f(0)**2
            else
              brf = State_GVB(0, j, k, br_, ib)*r_c(0)**2
            end if
            do i = 1, nb(1)
              State_GVB(i, j, k, br_, ib) = brf/r_c(i)**2
            end do
            State_GVB(1:nb(1), j, k, bt_:bp_, ib) = 0.0_dp
          end do
        end do
      end do
      call cz_unbind()
    end do
    ! psi 0 above the join
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (r_off + nb(1) <= k_join) cycle
      State_GVB(1:nb(1), 1:nb(2), 1:nb(3), psi_, ib) = 0.0_dp
    end do
    call cz_unbind()
    call arbor_fill_ghosts(State_GVB)
    iLimiter = lim_save
    write(msg, '(A,I0,A,ES10.3,A,F8.2,A)') '  dpf: ', min(it, dpf_maxit), ' iterations, relative residual ', &
        rnorm/bnorm, ', ', arbor_wtime() - t0, ' s'
    call arbor_print(msg)
    call b_field_stats('dpf', r_f(k_ss))
    call arbor_fill_ghosts(State_GVB)
  contains
    logical function in_solve(ib)
      integer, intent(in) :: ib
      in_solve = Key_IB(KI+1, ib)*nb(1) >= k_join .and. Key_IB(KI+1, ib)*nb(1) + nb(1) <= k_ss
    end function in_solve
    real(dp) function dot(a, b)
      real(dp), intent(in) :: a(:,:,:,:), b(:,:,:,:)
      integer :: jb_, ib_
      dot = 0.0_dp
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        if (in_solve(ib_)) dot = dot + sum(a(:, :, :, ib_)*b(:, :, :, ib_))
      end do
      dot = arbor_sum(dot)
    end function dot
    ! y = A x (ghosts of x filled first)
    subroutine apply_A(xin, y)
      real(dp), intent(in) :: xin(:,:,:,:)
      real(dp), intent(out) :: y(:,:,:,:)
      integer :: jb_, ib_, i_, j_, k_
      if (dpf_operator == 'wide') call wide_op(xin, y, .false.)
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        X(1:nb(1), 1:nb(2), 1:nb(3), 1, ib_) = xin(:, :, :, ib_)
      end do
      call arbor_fill_ghosts(X)
      yc = 0.0_dp
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        if (.not. in_solve(ib_)) cycle
        do k_ = 1, nb(3)
          do j_ = 1, nb(2)
            do i_ = 1, nb(1)
              associate (c => Cf(i_, j_, k_, :, ib_), x0 => X(i_, j_, k_, 1, ib_))
                yc(i_, j_, k_, ib_) = c(1)*(X(i_-1, j_, k_, 1, ib_) - x0) &
                    + merge(-c(2)*x0, c(2)*(X(i_+1, j_, k_, 1, ib_) - x0), Key_IB(KI+1, ib_)*nb(1) + i_ == k_ss) &
                    + c(3)*(X(i_, j_-1, k_, 1, ib_) - x0) + c(4)*(X(i_, j_+1, k_, 1, ib_) - x0) &
                    + c(5)*(X(i_, j_, k_-1, 1, ib_) - x0) + c(6)*(X(i_, j_, k_+1, 1, ib_) - x0)
              end associate
            end do
          end do
        end do
      end do
      if (dpf_operator == 'compact') then
        y = yc
      else                                           ! wide (+ eps compact); compact rows in the edge band
        do jb_ = 1, nBlock
          ib_ = iBlockList_I(jb_)
          do k_ = 1, nb(3)
            do j_ = 1, nb(2)
              if (edge(j_, k_, ib_)) then
                y(:, j_, k_, ib_) = yc(:, j_, k_, ib_)
              else
                y(:, j_, k_, ib_) = y(:, j_, k_, ib_) + dpf_eps*yc(:, j_, k_, ib_)
              end if
            end do
          end do
        end do
      end if
    end subroutine apply_A
    ! The wide operator: B = -grad x with the code's stencils in the solve
    ! cells (x in psi; below the join x continued, with the CZ's B_r if
    ! full, flat if not; odd above k_ss), its ghosts filled (the CZ's B below
    ! the join if full, 0 if not; above k_ss B_r even, B_h odd), then
    ! y = -V div B (V the cell volume): y = V div(grad x) for the homogeneous
    ! part. full: x = 0 gives -V div B of the CZ's field alone (-b)
    subroutine wide_op(xin, y, full)
      real(dp), intent(in) :: xin(:,:,:,:)
      real(dp), intent(out) :: y(:,:,:,:)
      logical, intent(in) :: full
      integer :: jb_, ib_, i_, j_, k_, m_
      real(dp) :: vol
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        if (in_solve(ib_)) then
          State_GVB(1:nb(1), 1:nb(2), 1:nb(3), psi_, ib_) = xin(:, :, :, ib_)
        else if (Key_IB(KI+1, ib_)*nb(1) >= k_join) then
          State_GVB(1:nb(1), 1:nb(2), 1:nb(3), psi_, ib_) = 0.0_dp
        end if
      end do
      call arbor_fill_ghosts(State_GVB, iVarMin=psi_, iVarMax=psi_)
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        call cz_bind(ib_)
        if (.not. in_solve(ib_)) cycle
        if (r_off == k_join) then
          if (full) then
            State_GVB(0, :, :, psi_, ib_) = State_GVB(1, :, :, psi_, ib_) &
                + State_GVB(0, :, :, br_, ib_)*(r_c(0)/r_f(0))**2*(r_c(1) - r_c(0))
            State_GVB(-1, :, :, psi_, ib_) = State_GVB(0, :, :, psi_, ib_) + State_GVB(0, :, :, br_, ib_)*(r_c(0) - r_c(-1))
          else
            State_GVB(0, :, :, psi_, ib_) = State_GVB(1, :, :, psi_, ib_)
            State_GVB(-1, :, :, psi_, ib_) = State_GVB(1, :, :, psi_, ib_)
          end if
        end if
        if (r_off + nb(1) == k_ss) then
          State_GVB(nb(1)+1, :, :, psi_, ib_) = -State_GVB(nb(1), :, :, psi_, ib_)
          State_GVB(nb(1)+2, :, :, psi_, ib_) = -State_GVB(nb(1)-1, :, :, psi_, ib_)
        end if
        call ddr_block(State_GVB(:, :, :, psi_, ib_), Dr_C)
        State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_, ib_) = -Dr_C
        call arbor_deriv_cell(State_GVB, psi_, 2, ib_, D_C)
        do i_ = 1, nb(1)
          State_GVB(i_, 1:nb(2), 1:nb(3), bt_, ib_) = -D_C(i_, :, :)/r_c(i_)
        end do
        call arbor_deriv_cell(State_GVB, psi_, 3, ib_, D_C)
        associate (th => Metric_DB(2, ib_)%X_I)
          do j_ = 1, nb(2)
            do i_ = 1, nb(1)
              State_GVB(i_, j_, 1:nb(3), bp_, ib_) = -D_C(i_, j_, :)/(r_c(i_)*sin(th(j_)))
            end do
          end do
        end associate
      end do
      call cz_unbind()
      call arbor_fill_ghosts(State_GVB, iVarMin=br_, iVarMax=bp_)
      y = 0.0_dp
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        call cz_bind(ib_)
        if (.not. in_solve(ib_)) cycle
        if (r_off == k_join .and. .not. full) State_GVB(-1:0, :, :, br_:bp_, ib_) = 0.0_dp
        if (r_off + nb(1) == k_ss) then
          do m_ = 1, 2
            State_GVB(nb(1)+m_, :, :, br_, ib_) = State_GVB(nb(1)+1-m_, :, :, br_, ib_)
            State_GVB(nb(1)+m_, :, :, bt_:bp_, ib_) = -State_GVB(nb(1)+1-m_, :, :, bt_:bp_, ib_)
          end do
        end if
        call div_block(State_GVB(:, :, :, br_:bp_, ib_), ib_, divB_C)
        associate (vt => Metric_DB(2, ib_)%Vol_I, vp => Metric_DB(3, ib_)%Vol_I)
          do k_ = 1, nb(3)
            do j_ = 1, nb(2)
              do i_ = 1, nb(1)
                vol = (r_f(i_)**3 - r_f(i_-1)**3)/3.0_dp*vt(j_)*vp(k_)
                y(i_, j_, k_, ib_) = -vol*divB_C(i_, j_, k_)
              end do
            end do
          end do
        end associate
      end do
      call cz_unbind()
    end subroutine wide_op
    ! z = M^-1 r: per column of each block, the radial tridiagonal part of A
    ! (the couplings out of the block dropped, the diagonal kept whole)
    subroutine precond(rin, z)
      real(dp), intent(in) :: rin(:,:,:,:)
      real(dp), intent(out) :: z(:,:,:,:)
      real(dp) :: dd(nb(1)), lo(nb(1)), up(nb(1)), cp(nb(1)), dp_(nb(1)), m_
      integer :: jb_, ib_, i_, j_, k_
      z = 0.0_dp
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        if (.not. in_solve(ib_)) cycle
        do k_ = 1, nb(3)
          do j_ = 1, nb(2)
            do i_ = 1, nb(1)
              associate (c => Cf(i_, j_, k_, :, ib_))
                dd(i_) = -(c(1) + c(2) + c(3) + c(4) + c(5) + c(6))
                lo(i_) = merge(c(1), 0.0_dp, i_ > 1)
                up(i_) = merge(c(2), 0.0_dp, i_ < nb(1))
              end associate
            end do
            ! Thomas
            cp(1) = up(1)/dd(1)
            dp_(1) = rin(1, j_, k_, ib_)/dd(1)
            do i_ = 2, nb(1)
              m_ = dd(i_) - lo(i_)*cp(i_-1)
              cp(i_) = up(i_)/m_
              dp_(i_) = (rin(i_, j_, k_, ib_) - lo(i_)*dp_(i_-1))/m_
            end do
            z(nb(1), j_, k_, ib_) = dp_(nb(1))
            do i_ = nb(1) - 1, 1, -1
              z(i_, j_, k_, ib_) = dp_(i_) - cp(i_)*z(i_+1, j_, k_, ib_)
            end do
          end do
        end do
      end do
    end subroutine precond
  end subroutine dpf_solve

  ! The field's curl and divergence (the code's stencils, ghosts filled) as
  ! q = |curl B| r dtheta/|B| and |div B| r dtheta/|B| (B^2-weighted rms;
  ! max of q) in the atmosphere cells whose stencils stay above the join and
  ! below 0.95 rss, patch interior / patch edge (within 2.5 cells of the
  ! Yin-Yang boundary), z < 5 Mm, 5-50 Mm, > 50 Mm
  subroutine b_field_stats(label, rss)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: rss
    real(dp) :: acc(4, 2, 3), th, ph, dth, dph, q, b2, z, c2, d2, rj
    integer :: jb, ib, i, j, k, ireg, ih
    character(len=200) :: msg
    character(len=8), parameter :: rname(2) = ['interior', 'edge    ']
    character(len=10), parameter :: hname(3) = ['z<5Mm     ', '5-50Mm    ', 'z>50Mm    ']
    call cz_unbind()
    rj = r_f(k_join)
    call arbor_fill_ghosts(State_GVB)
    acc = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (r_off + nb(1) <= k_join) cycle
      B_GV = State_GVB(:, :, :, br_:bp_, ib)
      call curl_block(B_GV, ib, J_CV)
      call div_block(State_GVB(:, :, :, br_:bp_, ib), ib, divB_C)
      dth = Metric_DB(2, ib)%Dx_I(1)
      dph = Metric_DB(3, ib)%Dx_I(1)
      do k = 1, nb(3)
        do j = 1, nb(2)
          th = Metric_DB(2, ib)%X_I(j)
          ph = Metric_DB(3, ib)%X_I(k)
          ireg = 1
          if (th < 0.25_dp*pi + 2.5_dp*dth .or. th > 0.75_dp*pi - 2.5_dp*dth .or. &
              abs(ph) > 0.75_dp*pi - 2.5_dp*dph) ireg = 2
          do i = 1, nb(1)
            if (r_off + i < k_join + 3 .or. r_c(i) > 0.95_dp*rss) cycle
            z = r_c(i) - rj
            ih = 1
            if (z > 5.0e8_dp) ih = 2
            if (z > 5.0e9_dp) ih = 3
            b2 = sum(State_GVB(i, j, k, br_:bp_, ib)**2)
            if (b2 <= 0.0_dp) cycle
            c2 = sum(J_CV(i, j, k, :)**2)*(r_c(i)*dth)**2
            d2 = (divB_C(i, j, k)*r_c(i)*dth)**2
            q = sqrt(c2/b2)
            acc(:, ireg, ih) = [acc(1, ireg, ih) + b2, acc(2, ireg, ih) + c2, max(acc(3, ireg, ih), q), acc(4, ireg, ih) + d2]
          end do
        end do
      end do
    end do
    call cz_unbind()
    do ih = 1, 3
      do ireg = 1, 2
        acc(1, ireg, ih) = arbor_sum(acc(1, ireg, ih))
        acc(2, ireg, ih) = arbor_sum(acc(2, ireg, ih))
        acc(3, ireg, ih) = arbor_max(acc(3, ireg, ih))
        acc(4, ireg, ih) = arbor_sum(acc(4, ireg, ih))
        write(msg, '(A,A,1X,A,1X,A,A,ES10.3,A,ES10.3,A,ES10.3)') '  field ', label, rname(ireg), hname(ih), &
            ': rms |curl B| r dth/|B| ', sqrt(acc(2, ireg, ih)/max(acc(1, ireg, ih), tiny(1.0_dp))), &
            '  max ', acc(3, ireg, ih), '  rms |div B| r dth/|B| ', sqrt(acc(4, ireg, ih)/max(acc(1, ireg, ih), tiny(1.0_dp)))
        call arbor_print(msg)
      end do
    end do
  end subroutine b_field_stats

  ! Magnetofriction (mf_iter > 0, after b_seed = 'dpf'; see cz3d_params_mod):
  ! the field in the atmosphere cells below rss moved towards zero J x B as
  ! the code computes it, the gas frozen. The displacement d is stored in the
  ! velocity slots (0 elsewhere) for the ghost fill; v and psi restored after.
  subroutine mf_relax(rss)
    real(dp), intent(in) :: rss
    real(dp), allocatable :: Vsave(:,:,:,:,:)
    real(dp) :: b(3), jj(3), f(3), h(3), b2, b2m, dmax, fsum(4), fe, t0
    integer :: it, jb, ib, i, j, k
    logical :: edg
    character(len=200) :: msg
    t0 = arbor_wtime()
    call cz_unbind()
    allocate(Vsave(nb(1), nb(2), nb(3), 3, nBlockMax))
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      Vsave(:, :, :, :, ib) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), vr_:vp_, ib)
    end do
    call b_field_stats('before mf', rss)
    do it = 1, mf_iter
      call arbor_fill_ghosts(State_GVB)
      dmax = 0.0_dp
      fsum = 0.0_dp                ! sum B^2 and sum (|F| h/B^2)^2 B^2, interior / edge
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        State_GVB(:, :, :, vr_:vp_, ib) = 0.0_dp
        if (r_off + nb(1) <= k_join) cycle
        B_GV = State_GVB(:, :, :, br_:bp_, ib)
        call curl_block(B_GV, ib, J_CV)
        associate (th => Metric_DB(2, ib)%X_I, dth => Metric_DB(2, ib)%Dx_I, ph => Metric_DB(3, ib)%X_I, &
            dph => Metric_DB(3, ib)%Dx_I)
          b2m = 0.0_dp
          do k = 1, nb(3)
            do j = 1, nb(2)
              edg = th(j) < 0.25_dp*pi + 2.5_dp*dth(j) .or. th(j) > 0.75_dp*pi - 2.5_dp*dth(j) &
                  .or. abs(ph(k)) > 0.75_dp*pi - 2.5_dp*dph(k)
              do i = 1, nb(1)
                if (r_off + i <= k_join .or. r_c(i) > rss) cycle
                b = State_GVB(i, j, k, br_:bp_, ib)
                jj = J_CV(i, j, k, :)
                b2 = sum(b**2) + 1.0e-6_dp
                f = [jj(2)*b(3) - jj(3)*b(2), jj(3)*b(1) - jj(1)*b(3), jj(1)*b(2) - jj(2)*b(1)]
                h = [r_f(i) - r_f(i-1), r_c(i)*dth(j), r_c(i)*sin(th(j))*dph(k)]
                State_GVB(i, j, k, vr_:vp_, ib) = mf_alpha*minval(h)**2*f/b2   ! (isotropic: the step limited by the smallest scale)
                dmax = max(dmax, maxval(abs(State_GVB(i, j, k, vr_:vp_, ib)))/minval(h))
                fe = sum(f**2)*minval(h)**2/b2         ! (|F| h_min/B^2)^2 B^2
                if (edg) then
                  fsum(3:4) = fsum(3:4) + [b2, fe]
                else
                  fsum(1:2) = fsum(1:2) + [b2, fe]
                end if
              end do
            end do
          end do
        end associate
      end do
      call cz_unbind()
      call arbor_fill_ghosts(State_GVB, iVarMin=vr_, iVarMax=vp_)
      ! B += curl(d x B) in the relaxed cells
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        if (r_off + nb(1) <= k_join) cycle
        do k = 1 - ng, nb(3) + ng
          do j = 1 - ng, nb(2) + ng
            do i = 1 - ng, nb(1) + ng
              jj = State_GVB(i, j, k, vr_:vp_, ib)
              b = State_GVB(i, j, k, br_:bp_, ib)
              B_GV(i, j, k, :) = [jj(2)*b(3) - jj(3)*b(2), jj(3)*b(1) - jj(1)*b(3), jj(1)*b(2) - jj(2)*b(1)]
            end do
          end do
        end do
        call curl_block(B_GV, ib, C_CV)
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, nb(1)
              if (r_off + i <= k_join .or. r_c(i) > rss) cycle
              State_GVB(i, j, k, br_:bp_, ib) = State_GVB(i, j, k, br_:bp_, ib) + C_CV(i, j, k, :)
            end do
          end do
        end do
      end do
      call cz_unbind()
      if (it == 1 .or. mod(it, 25) == 0 .or. it == mf_iter) then
        fsum(1) = arbor_sum(fsum(1));  fsum(2) = arbor_sum(fsum(2))
        fsum(3) = arbor_sum(fsum(3));  fsum(4) = arbor_sum(fsum(4))
        dmax = arbor_max(dmax)
        write(msg, '(A,I5,A,ES10.3,A,ES10.3,A,ES10.3)') '  mf: iteration ', it, ': rms |J x B| h/B^2 interior ', &
            sqrt(fsum(2)/max(fsum(1), tiny(1.0_dp))), ', edge ', sqrt(fsum(4)/max(fsum(3), tiny(1.0_dp))), &
            ', max step/h ', dmax
        call arbor_print(msg)
      end if
    end do
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      State_GVB(1:nb(1), 1:nb(2), 1:nb(3), vr_:vp_, ib) = Vsave(:, :, :, :, ib)
    end do
    call arbor_fill_ghosts(State_GVB)
    write(msg, '(A,I0,A,F8.2,A)') '  mf: ', mf_iter, ' iterations, ', arbor_wtime() - t0, ' s'
    call arbor_print(msg)
    call b_field_stats('after mf', rss)
  end subroutine mf_relax

  ! ---------------------------------------------------------------- state transfer

  ! Every local block's interior state (all nv variables), one file per rank:
  ! out_dir/state_export_rNNNN.bin, stream access: n_block, nb(1:3), nv; then
  ! per block its key (KEY_LEN integers) and State(1:nb(1), 1:nb(2), 1:nb(3), 1:nv).
  ! The root also writes out_dir/state_export_meta.bin: t, the step, and the
  ! restart numbers (profile counter and schedule, the shell clocks: as
  ! cz_restart_save), so an import continues the run.
  subroutine export_state_write(t, iStep)
    real(dp), intent(in) :: t
    integer, intent(in) :: iStep
    integer :: u, jb, ib
    character(len=300) :: name
    real(dp), allocatable :: Data_I(:)
    if (.not. export_state) return
    write(name, '(A,A,I4.4,A)') trim(out_dir), '/state_export_r', my_rank, '.bin'
    open(newunit=u, file=trim(name), access='stream', form='unformatted', status='replace', action='write')
    write(u) nBlock, nb, nv
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      write(u) Key_IB(:, ib), State_GVB(1:nb(1), 1:nb(2), 1:nb(3), 1:nv, ib)
    end do
    close(u)
    if (.not. is_root) return
    Data_I = [real(i_prof, dp), t_next_prof, t_prof0, tloc]
    open(newunit=u, file=trim(out_dir)//'/state_export_meta.bin', access='stream', form='unformatted', &
        status='replace', action='write')
    write(u) t, iStep, size(Data_I), Data_I
    close(u)
  end subroutine export_state_write

  ! Index of the export files in import_dir (every rank: the keys only, read
  ! by seeking): file, position and key of every exported block
  subroutine import_open()
    integer :: u, n, nbf(3), nvf, f, k, m
    integer(i8) :: rec, pos
    character(len=300) :: name
    logical :: exists
    integer, allocatable :: key_K(:,:), file_K(:)
    integer(i8), allocatable :: pos_K(:)
    if (allocated(imp_key)) return
    m = 0
    allocate(key_K(KEY_LEN, 0), file_K(0), pos_K(0))
    do f = 0, 99999
      write(name, '(A,A,I4.4,A)') trim(import_dir), '/state_export_r', f, '.bin'
      inquire(file=trim(name), exist=exists)
      if (.not. exists) exit
      open(newunit=u, file=trim(name), access='stream', form='unformatted', status='old', action='read')
      read(u) n, nbf, nvf
      if (f == 0) then
        imp_nb = nbf
        imp_nv = nvf
      else if (any(nbf /= imp_nb) .or. nvf /= imp_nv) then
        call arbor_abort('cz3d: import: the export files differ in block size or variables')
      end if
      rec = 4_i8*KEY_LEN + 8_i8*product(int(nbf, i8))*nvf
      key_K = reshape([key_K, [(0, k = 1, KEY_LEN*n)]], [KEY_LEN, m + n])
      file_K = [file_K, [(f, k = 1, n)]]
      pos_K = [pos_K, [(0_i8, k = 1, n)]]
      do k = 1, n
        pos = 21_i8 + int(k - 1, i8)*rec                ! (stream positions count from 1; the header is 5 integers)
        read(u, pos=pos) key_K(:, m + k)
        pos_K(m + k) = pos
      end do
      m = m + n
      close(u)
    end do
    if (m == 0) call arbor_abort('cz3d: import: no state_export files in '//trim(import_dir))
    if (any(imp_nb(2:3) /= nb(2:3))) call arbor_abort('cz3d: import: the horizontal block size differs')
    if (mod(imp_nb(1), nb(1)) /= 0) call arbor_abort('cz3d: import: the radial block size must divide the exported one')
    call move_alloc(key_K, imp_key)
    call move_alloc(file_K, imp_file)
    call move_alloc(pos_K, imp_pos)
  end subroutine import_open

  ! The finest horizontal level of the exported blocks
  integer function import_max_level() result(l)
    l = maxval(imp_key(KL+2, :))
  end function import_max_level

  ! Refinement criterion for the import tree (arbor_tree_refine): refine a
  ! block until an exported block has its horizontal position and level
  subroutine import_criterion(key, x_lo, x_hi, flags)
    integer, intent(in) :: key(KEY_LEN)
    real(dp), intent(in) :: x_lo(3), x_hi(3)
    integer, intent(out) :: flags(3)
    flags = 0
    if (import_record(key, .false.) == 0) flags(2:3) = 1
  end subroutine import_criterion

  ! The exported block holding block key (horizontal position and level
  ! equal; with radial: also the radial block that contains it), 0 if none
  integer function import_record(key, radial) result(r)
    integer, intent(in) :: key(KEY_LEN)
    logical, intent(in) :: radial
    integer :: m
    r = 0
    do m = 1, size(imp_key, 2)
      if (imp_key(KP, m) /= key(KP)) cycle
      if (any(imp_key(KL+2:KL+3, m) /= key(KL+2:KL+3))) cycle
      if (any(imp_key(KI+2:KI+3, m) /= key(KI+2:KI+3))) cycle
      if (radial) then
        if (imp_key(KI+1, m) /= (key(KI+1)*nb(1))/imp_nb(1)) cycle
      end if
      r = m
      return
    end do
  end function import_record

  ! The block's state from the export files in import_dir: its radial part of
  ! the exported block at the same place (the export may hold whole columns
  ! or larger radial blocks); the variables both runs have are copied, the
  ! others stay as set (B, psi: 0). With the atmosphere, a block above the
  ! join without an exported block keeps its initial state (a CZ-only export,
  ! e.g. CZ_3D's, under the atmosphere of the init file)
  subroutine import_block(iBlock)
    integer, intent(in) :: iBlock
    real(dp), allocatable :: buf(:,:,:,:)
    integer :: m, u, i0
    character(len=300) :: name
    call import_open()
    m = import_record(Key_IB(:, iBlock), .true.)
    if (m == 0 .and. use_atmosphere .and. Key_IB(KI+1, iBlock)*nb(1) >= k_join) return
    if (m == 0) call arbor_abort('cz3d: import: no exported block at this place (a different grid?)')
    allocate(buf(imp_nb(1), imp_nb(2), imp_nb(3), imp_nv))
    write(name, '(A,A,I4.4,A)') trim(import_dir), '/state_export_r', imp_file(m), '.bin'
    open(newunit=u, file=trim(name), access='stream', form='unformatted', status='old', action='read')
    read(u, pos=imp_pos(m) + 4_i8*KEY_LEN) buf
    close(u)
    i0 = Key_IB(KI+1, iBlock)*nb(1) - imp_key(KI+1, m)*imp_nb(1)
    State_GVB(1:nb(1), 1:nb(2), 1:nb(3), 1:min(nv, imp_nv), iBlock) = buf(i0+1:i0+nb(1), :, :, 1:min(nv, imp_nv))
  end subroutine import_block

  ! The time, step and restart numbers of the export (state_export_meta.bin),
  ! if there: t and iStep set, the profile counter and the shell clocks
  ! restored; found false (nothing changed) without the file
  subroutine import_meta(t, iStep, found)
    real(dp), intent(inout) :: t
    integer, intent(inout) :: iStep
    logical, intent(out) :: found
    integer :: u, nd
    real(dp), allocatable :: Data_I(:)
    inquire(file=trim(import_dir)//'/state_export_meta.bin', exist=found)
    if (.not. found) return
    open(newunit=u, file=trim(import_dir)//'/state_export_meta.bin', access='stream', form='unformatted', &
        status='old', action='read')
    read(u) t, iStep, nd
    allocate(Data_I(nd))
    read(u) Data_I
    close(u)
    i_prof = nint(Data_I(1))
    t_next_prof = Data_I(2)
    t_prof0 = Data_I(3)
    if (nd == 3 + n_r) then
      tloc = Data_I(4:)
    else
      tloc = t                                            ! (a different radial grid: all clocks at the run's time)
    end if
  end subroutine import_meta

  ! Uniform in [-1, 1] from (patch, i, j, k): Park-Miller steps (no overflow)
  real(dp) function urand(p, i, j, k)
    integer, intent(in) :: p, i, j, k
    integer(i8), parameter :: a = 48271_i8, m = 2147483647_i8
    integer(i8) :: x
    integer :: n
    x = 12345_i8 + 999331_i8*int(p, i8)
    do n = 1, 3
      x = mod(a*mod(x + int(i, i8), m), m)
      x = mod(a*mod(x + 7919_i8*int(j, i8), m), m)
      x = mod(a*mod(x + 104729_i8*int(k, i8), m), m)
    end do
    urand = 2.0_dp*real(x, dp)/real(m, dp) - 1.0_dp
  end function urand

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
    if (at_top .and. wall_top /= 'outflow') then
      if (wall_top == 'face') then
        W(m+1, :, :) = -W(m, :, :);  W(m+2, :, :) = -W(m-1, :, :)
      else
        W(m, :, :) = 0.0_dp;  W(m+1, :, :) = -W(m-1, :, :);  W(m+2, :, :) = -W(m-2, :, :)
      end if
    end if
  end subroutine odd_wall_flux

  ! ================================================================ equation of state

  ! rho and p1 = p0 expm1(ln p - ln p0) on the whole block (ghosts included);
  ! with_cs: also the sound speed and T. The equation of state is evaluated a
  ! radial column at a time (eos_hhe's column evaluators: vectorized loops,
  ! the same values as the scalar calls).
  subroutine derived(ib, with_cs)
    integer, intent(in) :: ib
    logical, intent(in) :: with_cs
    real(dp) :: rp(1-ng:nb(1)+ng), s(1-ng:nb(1)+ng), q(1-ng:nb(1)+ng, 3), d(1-ng:nb(1)+ng), &
        e1(1-ng:nb(1)+ng), e2(1-ng:nb(1)+ng), e3(1-ng:nb(1)+ng)
    integer :: i, j, k, n
    n = nb(1) + 2*ng
    do k = 1 - ng, nb(3) + ng
      do j = 1 - ng, nb(2) + ng
        do i = 1 - ng, nb(1) + ng
          rho_G(i, j, k) = rho0(i) + State_GVB(i, j, k, rho1_, ib)
          rp(i) = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)                     ! physical density
          s(i) = s0(i) + State_GVB(i, j, k, s1_, ib)
        end do
        if (with_cs) then
          call eos_fields_n(tab, [F_LNP, F_LNT, F_GAMMA1], n, rp, s, q)
          d = q(:, 1) - lnp0(1-ng:nb(1)+ng)
          call vexpm1(n, d, e1)
          call vexp(n, q(:, 1), e2)
          call vexp(n, q(:, 2), e3)
          do i = 1 - ng, nb(1) + ng
            p1_G(i, j, k) = p0(i)*e1(i)
            cs_G(i, j, k) = sqrt(q(i, 3)*e2(i)/rp(i))/xi(i)                        ! reduced (RSST)
            T_G(i, j, k) = e3(i)
          end do
          ! with B: the fast speed's bound c_red + v_A (the diffusion's speed; v_A capped, see va_eff)
          if (use_B) then
            do i = 1 - ng, nb(1) + ng
              cs_G(i, j, k) = cs_G(i, j, k) + va_eff(sum(State_GVB(i, j, k, br_:bp_, ib)**2), rp(i), &
                  cell_cap(i, j, k, ib, cs_G(i, j, k)))
            end do
          end if
        else
          call eos_field_n(tab, F_LNP, n, rp, s, q(:, 1))
          d = q(:, 1) - lnp0(1-ng:nb(1)+ng)
          call vexpm1(n, d, e1)
          do i = 1 - ng, nb(1) + ng
            p1_G(i, j, k) = p0(i)*e1(i)
          end do
        end if
      end do
    end do
  end subroutine derived

  ! ================================================================ time step

  ! dt < cfl / sum_d (|v_d| + cs/xi)/l_d, l_r = CZTR's dr/dxi
  real(dp) function cz_stable_dt(self, iBlock) result(dt)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iBlock
    real(dp) :: cs, rate, rr, sth, rc, rp(nb(1)), s(nb(1)), q(nb(1), 2), e(nb(1))
    integer :: i, j, k
    call cz_bind(iBlock)
    rate = 0.0_dp
    do k = 1, nb(3)
      do j = 1, nb(2)
        sth = sin(Metric_DB(2, iBlock)%X_I(j))
        do i = 1, nb(1)
          rp(i) = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, iBlock)
          s(i) = s0(i) + State_GVB(i, j, k, s1_, iBlock)
        end do
        call eos_fields_n(tab, [F_LNP, F_GAMMA1], nb(1), rp, s, q)
        call vexp(nb(1), q(:, 1), e)
        do i = 1, nb(1)
          cs = sqrt(q(i, 2)*e(i)/rp(i))/xi(i)
          if (use_B) cs = cs + va_eff(sum(State_GVB(i, j, k, br_:bp_, iBlock)**2), rp(i), cell_cap(i, j, k, iBlock, cs))
          rr = r_c(i)
          rc = (abs(State_GVB(i, j, k, vr_, iBlock)) + cs)/drdxi_c(i) &
              + (abs(State_GVB(i, j, k, vt_, iBlock)) + cs)/(rr*Metric_DB(2, iBlock)%Dx_I(j)) &
              + (abs(State_GVB(i, j, k, vp_, iBlock)) + cs)/(rr*sth*Metric_DB(3, iBlock)%Dx_I(k))
          rate = max(rate, rc)
          rate_r(i) = max(rate_r(i), rc)
        end do
      end do
    end do
    dt = min(cfl/rate, dt_max)
    if (use_atmosphere) dt = min(dt, atm_dt(iBlock, rate_r(1:nb(1))))
    call cz_unbind()
  end function cz_stable_dt

  subroutine cz_advance(self, iBlocks, t, dt)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iBlocks(:)
    real(dp), intent(in) :: t, dt
    integer :: stage, jb, ib, i, j, k, v
    real(dp) :: cdt
    real(dp), parameter :: c_stage(4) = [0.0_dp, 0.25_dp, 1.0_dp/3.0_dp, 0.5_dp]
    call local_steps(dt)
    dt_now = dt
    if (freeze_cz) then
      do jb = 1, size(iBlocks)
        if (frz_have(iBlocks(jb))) cycle
        frz_B(:, :, :, :, iBlocks(jb)) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), 1:nv, iBlocks(jb))
        frz_have(iBlocks(jb)) = .true.
      end do
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
    do stage = 1, 4
      call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.true., time=t + c_stage(stage)*dt)
      do jb = 1, size(iBlocks)
        ib = iBlocks(jb)
        call cz_bind(ib)
        want_work = budget_on .and. stage == 1
        call rhs(ib)
        want_work = .false.
        ! the RK update, with each shell's step factor tsc (local time stepping):
        ! one pass in memory order (the same arithmetic as scaling R_CV first)
        cdt = dt/real(5 - stage, dp)
        if (dt1_now > 0.0_dp) then
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
      end do
      call cz_unbind()
    end do
    if (budget_on) then
      do jb = 1, size(iBlocks)
        call budget_hydro(iBlocks(jb), dt)
      end do
    end if
    if (use_atmosphere) then
      do jb = 1, size(iBlocks)
        call cz_bind(iBlocks(jb))
        call atm_floors(iBlocks(jb))
      end do
      call cz_unbind()
    end if
    call cz_transport_all(iBlocks, dt, .false., t + dt)
    if (use_atmosphere) then
      ! conduction (the ghosts in every direction), losses, heating, then the floors
      call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.true., time=t + dt)
      do jb = 1, size(iBlocks)
        call cz_bind(iBlocks(jb))
        call atm_block(iBlocks(jb), dt, at_bottom, at_top, shell_steps())
        call atm_floors(iBlocks(jb))
      end do
      call cz_unbind()
    end if
    if (use_diffusion) then
      call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.true., time=t + dt)
      do jb = 1, size(iBlocks)
        call cz_bind(iBlocks(jb))
        call diffusion(iBlocks(jb), dt)
        if (use_atmosphere) call atm_floors(iBlocks(jb))
      end do
      call cz_unbind()
    end if
    ! freeze_cz: the CZ cells back to their stored state
    if (freeze_cz) then
      do jb = 1, size(iBlocks)
        call cz_bind(iBlocks(jb))
        do i = 1, min(nb(1), k_join - r_off)
          State_GVB(i, 1:nb(2), 1:nb(3), 1:nv, iBlocks(jb)) = frz_B(i, :, :, :, iBlocks(jb))
        end do
      end do
      call cz_unbind()
    end if
    ! friction on v above the join (v_damp_n > 0), for the run's first steps (n_adv: local_steps)
    if (v_damp_n > 0.0_dp .and. use_atmosphere) then
      if (n_adv <= v_damp_steps) then
        cdt = 1.0_dp
      else if (n_adv <= v_damp_steps + v_damp_ramp) then
        cdt = 1.0_dp - real(n_adv - v_damp_steps, dp)/real(v_damp_ramp + 1, dp)
      else
        cdt = 0.0_dp
      end if
      if (cdt > 0.0_dp) then
        cdt = 1.0_dp - cdt/(v_damp_n + 1.0_dp)
        do jb = 1, size(iBlocks)
          call cz_bind(iBlocks(jb))
          do i = max(1, k_join - r_off + 1), nb(1)
            ! (v_r towards the init file's profile: the wind keeps flowing)
            State_GVB(i, 1:nb(2), 1:nb(3), vr_, iBlocks(jb)) = v_init(i) &
                + cdt*(State_GVB(i, 1:nb(2), 1:nb(3), vr_, iBlocks(jb)) - v_init(i))
            State_GVB(i, 1:nb(2), 1:nb(3), vt_:vp_, iBlocks(jb)) = cdt*State_GVB(i, 1:nb(2), 1:nb(3), vt_:vp_, iBlocks(jb))
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
    dt_last = dt
    if (UseMultirate) then
      tloc = max(tloc, t + dt)                          ! (every level ends at the level-1 step's end)
    else
      tloc = tloc + tsc_g*dt
    end if
  contains
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
  ! > 0: dt1 scaled by max(0, 1 - n/dt1_ramp_steps), n this run's steps.
  subroutine local_steps(dt)
    real(dp), intent(in) :: dt
    real(dp) :: a(1:n_r)
    integer :: i
    dt1_now = dt1
    if (dt1_ramp_steps > 0) dt1_now = dt1*max(0.0_dp, 1.0_dp - real(n_adv, dp)/real(dt1_ramp_steps, dp))
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
    call derived(ib, .false.)
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

      ! ---- momentum and entropy: advection -(v.grad) q for q = v_r, v_theta, v_phi, s1
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

      ! ---- pressure gradient, gravity on rho1, curvature terms, ds0/dr
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
    real(dp) :: b(3), v(3), jj(3), ch, rho
    ! J = curl B
    B_GV = State_GVB(:, :, :, br_:bp_, ib)
    call curl_block(B_GV, ib, J_CV)
    ! E = v x B (with the ghosts), its curl
    do k = 1 - ng, nb(3) + ng
      do j = 1 - ng, nb(2) + ng
        do i = 1 - ng, nb(1) + ng
          v = State_GVB(i, j, k, vr_:vp_, ib)
          b = State_GVB(i, j, k, br_:bp_, ib)
          B_GV(i, j, k, :) = [v(2)*b(3) - v(3)*b(2), v(3)*b(1) - v(1)*b(3), v(1)*b(2) - v(2)*b(1)]
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
    associate (th => Metric_DB(2, ib)%X_I)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            R_CV(i, j, k, bp_) = C_CV(i, j, k, 3) - D_C(i, j, k)/(r_c(i)*sin(th(j)))
            b = State_GVB(i, j, k, br_:bp_, ib)
            jj = J_CV(i, j, k, :)
            rho = rho_G(i, j, k)
            R_CV(i, j, k, vr_:vp_) = R_CV(i, j, k, vr_:vp_) &
                + lorentz_f(sum(b**2), rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib), &
                cell_cap(i, j, k, ib, cell_cs(i, j, k, rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)))) &
                *[jj(2)*b(3) - jj(3)*b(2), jj(3)*b(1) - jj(1)*b(3), jj(1)*b(2) - jj(2)*b(1)]/(four_pi*rho)
            ch = ch_cell(i, j, k, ib)
            R_CV(i, j, k, psi_) = -ch**2*divB_C(i, j, k)
          end do
        end do
      end do
    end associate
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
    associate (th => Metric_DB(2, ib)%X_I)
      ! C_r
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          W_G(:, j, k) = sin(th(j))*A_GV(:, j, k, 3)
        end do
      end do
      call arbor_deriv_cell(W_G, 2, ib, D_C)
      W_G = A_GV(:, :, :, 2)
      call arbor_deriv_cell(W_G, 3, ib, Dr_C)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            C(i, j, k, 1) = (D_C(i, j, k) - Dr_C(i, j, k))/(r_c(i)*sin(th(j)))
          end do
        end do
      end do
      ! C_theta
      W_G = A_GV(:, :, :, 1)
      call arbor_deriv_cell(W_G, 3, ib, D_C)
      do i = 1 - ng, nb(1) + ng
        W_G(i, :, :) = r_c(i)*A_GV(i, :, :, 3)
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
      do i = 1 - ng, nb(1) + ng
        W_G(i, :, :) = r_c(i)*A_GV(i, :, :, 2)
      end do
      call ddr_block(W_G, Dr_C)
      W_G = A_GV(:, :, :, 1)
      call arbor_deriv_cell(W_G, 2, ib, D_C)
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
      do i = 1 - ng, nb(1) + ng
        W_G(i, :, :) = r_c(i)**2*A_GV(i, :, :, 1)
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
      W_G = A_GV(:, :, :, 3)
      call arbor_deriv_cell(W_G, 3, ib, Dr_C)
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
    ch = cs0(i)/xi(i) + va_eff(sum(State_GVB(i, j, k, br_:bp_, ib)**2), rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib), &
        cell_cap(i, j, k, ib, cs0(i)/xi(i)))
    if (glm_ch /= 'cfl' .or. dt_now <= 0.0_dp) return
    s_inv = 1.0_dp/drdxi_c(i) + 1.0_dp/(r_c(i)*Metric_DB(2, ib)%Dx_I(j)) &
        + 1.0_dp/(r_c(i)*sin(Metric_DB(2, ib)%X_I(j))*Metric_DB(3, ib)%Dx_I(k))
    dti = tsc(i)*dt_now
    ch = max(ch, glm_cfl_frac*cfl/(s_inv*dti))
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
    f = 1.0_dp
    if (c > 0.0_dp) f = 1.0_dp/sqrt(1.0_dp + (b2/(four_pi*rho*c**2))**2)
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
    c = va_max
    if (va_dt_target <= 0.0_dp) return
    dxm = min(drdxi_c(i), r_c(i)*Metric_DB(2, ib)%Dx_I(j), &
        r_c(i)*sin(Metric_DB(2, ib)%X_I(j))*Metric_DB(3, ib)%Dx_I(k))
    fl = 1.0e5_dp
    if (use_atmosphere) fl = max(fl, va_floor_cs*cs*atm_switch(z_c(i)))
    c = max(cfl*dxm/va_dt_target - (cs + norm2(State_GVB(i, j, k, vr_:vp_, ib))), fl)
  end function cell_cap

  ! The cell's signal speed from the pressure p0 + p1 (Gamma1 ~ 5/3; for the
  ! cap only), reduced by xi
  real(dp) function cell_cs(i, j, k, rho) result(cs)
    integer, intent(in) :: i, j, k
    real(dp), intent(in) :: rho
    cs = sqrt(5.0_dp/3.0_dp*max(p0(i) + p1_G(i, j, k), 1.0e-30_dp)/rho)/xi(i)
  end function cell_cs

  ! Longitude means in the global frame (rotation axis = global z), per radial
  ! cell and latitude bin (n_lat_bin bins of equal latitude width): the zonal
  ! velocity <v_phi_g>, Omega - Omega_0 = <v_phi_g>/(r cos lat) [nHz], and
  ! the meridional velocity <v_theta_g> (> 0 southward); solid-angle weights
  ! with the Yin-Yang overlap counted once. out_dir/rotation_NNNN.dat (the
  ! profile's number).
  subroutine rotation_profile(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    integer, parameter :: n_lat_bin = 36
    real(dp) :: A_IV(1:n_r, 3*n_lat_bin), e(3,3), x(3), vc(3), w, th, ph, th2, ph2, rs, lat, eph(3), eth(3), lc
    integer :: jb, ib, i, j, k, m, u
    character(len=300) :: name
    A_IV = 0.0_dp
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
          x = e(:, 1)
          if (Key_IB(KP, ib) == 2) x = yy_cart(x)
          rs = max(hypot(x(1), x(2)), 1.0e-12_dp)
          lat = asin(max(-1.0_dp, min(1.0_dp, x(3))))
          m = min(n_lat_bin, 1 + int((lat + 0.5_dp*pi)/pi*n_lat_bin))
          eph = [-x(2), x(1), 0.0_dp]/rs
          eth = [x(3)*x(1)/rs, x(3)*x(2)/rs, -rs]
          do i = 1, nb(1)
            vc = matmul(e, State_GVB(i, j, k, vr_:vp_, ib))
            if (Key_IB(KP, ib) == 2) vc = yy_cart(vc)
            A_IV(r_off + i, m) = A_IV(r_off + i, m) + w
            A_IV(r_off + i, n_lat_bin + m) = A_IV(r_off + i, n_lat_bin + m) + w*dot_product(vc, eph)
            A_IV(r_off + i, 2*n_lat_bin + m) = A_IV(r_off + i, 2*n_lat_bin + m) + w*dot_product(vc, eth)
          end do
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(A_IV)
    if (.not. is_root) return
    write(name, '(A,A,I4.4,A)') trim(out_dir), '/rotation_', i_prof, '.dat'
    open(newunit=u, file=trim(name), status='replace', action='write')
    write(u, '(A,ES14.6,A,I0,A,ES12.4,A)') '# t = ', t, ' s; step ', iStep, '; omega_rot = ', omega_rot, &
        ' rad/s; longitude means in the global frame'
    write(u, '(A)') '# i  r[cm]  z[km]  lat[deg]  <v_phi>[cm/s]  Omega-Omega0[nHz]  <v_theta>[cm/s] (>0 southward)'
    do i = 1, n_r
      do m = 1, n_lat_bin
        if (A_IV(i, m) <= 0.0_dp) cycle
        lc = -0.5_dp*pi + (m - 0.5_dp)*pi/n_lat_bin
        write(u, '(I5,ES16.8,F12.1,F8.2,3ES14.5)') i, r_c(i), z_c(i)/1.0e5_dp, lc*180.0_dp/pi, &
            A_IV(i, n_lat_bin + m)/A_IV(i, m), A_IV(i, n_lat_bin + m)/A_IV(i, m)/(r_c(i)*cos(lc))/(2.0_dp*pi)*1.0e9_dp, &
            A_IV(i, 2*n_lat_bin + m)/A_IV(i, m)
      end do
    end do
    close(u)
  contains
    logical function core(tt, pp)
      real(dp), intent(in) :: tt, pp
      core = tt >= 0.25_dp*pi .and. tt <= 0.75_dp*pi .and. abs(pp) <= 0.75_dp*pi
    end function core
  end subroutine rotation_profile

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

  ! Remove the non-Coriolis part of L_z's change since the reference: delta_L
  ! = L_z - (lz_ref + lz_imp), delta_Omega = delta_L / I_z (I_z = int rho s^2
  ! dV), v -= delta_Omega z x x (global frame, back to each patch's components)
  subroutine lz_correct(t)
    real(dp), intent(in) :: t
    integer :: jb, ib, i, j, k
    real(dp) :: e(3,3), x(3), dv(3), w, th, ph, th2, ph2, lz, iz, dom
    character(len=200) :: msg
    if (.not. lz_ref_set) return
    call angular_momentum(lz)
    iz = 0.0_dp
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
          x = e(:, 1)
          if (Key_IB(KP, ib) == 2) x = yy_cart(x)
          do i = 1, nb(1)
            iz = iz + w*Metric_DB(1, ib)%Vol_I(i)*(rho0(i) + State_GVB(i, j, k, rho1_, ib))*r_c(i)**2*(x(1)**2 + x(2)**2)
          end do
        end do
      end do
    end do
    call cz_unbind()
    iz = arbor_sum(iz)
    lz_ref = lz_ref + lz_imp
    lz_imp = 0.0_dp
    dom = (lz - lz_ref)/iz
    ! v -= dom (z x x) everywhere (all cells, the overlap included)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
          do i = 1, nb(1)
            x = r_c(i)*e(:, 1)
            if (Key_IB(KP, ib) == 2) x = yy_cart(x)
            dv = -dom*[-x(2), x(1), 0.0_dp]
            if (Key_IB(KP, ib) == 2) dv = yy_cart(dv)
            State_GVB(i, j, k, vr_:vp_, ib) = State_GVB(i, j, k, vr_:vp_, ib) + matmul(transpose(e), dv)
          end do
        end do
      end do
    end do
    call cz_unbind()
    if (lz_correct_dt > 0.0_dp) then
      do while (t_next_lzc <= t)
        t_next_lzc = t_next_lzc + lz_correct_dt
      end do
    end if
    if (abs(dom) > 1.0e-3_dp*max(abs(omega_rot), 1.0e-6_dp)) then
      write(msg, '(A,ES11.3,A)') 'cz3d: large L_z correction, delta_Omega = ', dom, ' rad/s'
      call arbor_print(msg)
    end if
  end subroutine lz_correct

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

  ! Energy exchange between the flow and the field over the sphere (ghosts
  ! filled): the Lorentz force's work wl = int v.(J x B)/(4 pi) dV and the
  ! induction term's change of the magnetic energy wi = int B.curl(v x B)/(4 pi)
  ! dV [erg/s]. With walls that pass no Poynting flux (conductor: B_r = v_r =
  ! 0; radial: B_h = v_r = 0), wl + wi = 0 up to the discretisation: a check of
  ! the signs and factors of both terms. Yin-Yang overlap weights as the budget's.
  subroutine mhd_exchange(wl, wi)
    real(dp), intent(out) :: wl, wi
    integer :: jb, ib, i, j, k
    real(dp) :: w, th, ph, th2, ph2, v(3), b(3), jj(3), vol
    wl = 0.0_dp;  wi = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      B_GV = State_GVB(:, :, :, br_:bp_, ib)
      call curl_block(B_GV, ib, J_CV)
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          do i = 1 - ng, nb(1) + ng
            v = State_GVB(i, j, k, vr_:vp_, ib)
            b = State_GVB(i, j, k, br_:bp_, ib)
            B_GV(i, j, k, :) = [v(2)*b(3) - v(3)*b(2), v(3)*b(1) - v(1)*b(3), v(1)*b(2) - v(2)*b(1)]
          end do
        end do
      end do
      call curl_block(B_GV, ib, C_CV)
      do k = 1, nb(3)
        do j = 1, nb(2)
          th = Metric_DB(2, ib)%X_I(j);  ph = Metric_DB(3, ib)%X_I(k)
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          if (.not. core(th, ph)) cycle
          call yy_angles(th, ph, th2, ph2)
          if (core(th2, ph2)) w = 0.5_dp*w
          do i = 1, nb(1)
            vol = w*Metric_DB(1, ib)%Vol_I(i)
            v = State_GVB(i, j, k, vr_:vp_, ib)
            b = State_GVB(i, j, k, br_:bp_, ib)
            jj = J_CV(i, j, k, :)
            wl = wl + vol*lorentz_f(sum(b**2), rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib), &
                cell_cap(i, j, k, ib, cs0(i)/xi(i))) &
                *dot_product(v, [jj(2)*b(3) - jj(3)*b(2), jj(3)*b(1) - jj(1)*b(3), jj(1)*b(2) - jj(2)*b(1)])/four_pi
            wi = wi + vol*dot_product(b, C_CV(i, j, k, :))/four_pi
          end do
        end do
      end do
    end do
    call cz_unbind()
    wl = arbor_sum(wl);  wi = arbor_sum(wi)
  contains
    logical function core(t, p)
      real(dp), intent(in) :: t, p
      core = t >= 0.25_dp*pi .and. t <= 0.75_dp*pi .and. abs(p) <= 0.75_dp*pi
    end function core
  end subroutine mhd_exchange

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
      if (present(time)) then
        call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.true., time=time)
      else
        call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.true.)
      end if
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
  end subroutine cz_transport_all

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
    ! the top-layer blend weights at the radial faces (wr) and the cells (wh: horizontal faces)
    wr = 0.0_dp;  wh = 0.0_dp
    if (diff_top_km > 0.0_dp) then
      do i = 0, nb(1)
        wr(i) = top_weight(r_top_cz - r_f(i))
        if (i >= 1) wh(i) = top_weight(r_top_cz - r_c(i))
      end do
    end if
    call derived(ib, .true.)
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
      ! all fluxes with the donor values
      call faces_b(1, 0, 0, F1_FV, C1_F)
      call faces_b(0, 1, 0, F2_FV, C2_F)
      call faces_b(0, 0, 1, F3_FV, C3_F)
    end associate
    call wall_fluxes(1, nVarCz)
    call radial_energy_flux(ib)
    if (diffusion_heating) then
      ! the energy each face's fluxes carry (variable 6), with h + Phi per cell
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          do i = 1 - ng, nb(1) + ng
            hP_G(i, j, k) = h0_r(i) + T_G(i, j, k)*State_GVB(i, j, k, s1_, ib) + p1_G(i, j, k)/rho_G(i, j, k) + phi_r(i)
          end do
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
    call arbor_flux_divergence(ib, F1_FV, F2_FV, F3_FV, Div_CV)
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, nb(1)
          dti = tsc(i)*dt                                   ! (the shell's own step)
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
          ho = h0_r(i) + T_G(i, j, k)*s1o + p1_G(i, j, k)/rho_G(i, j, k) + phi_r(i)   ! h + Phi
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
            F(i, j, k, rho1_) = -C(i, j, k)*jump4(State_GVB(i-di, j-dj, k-dk, rho1_, ib), State_GVB(i, j, k, rho1_, ib), &
                State_GVB(i+di, j+dj, k+dk, rho1_, ib), State_GVB(i+2*di, j+2*dj, k+2*dk, rho1_, ib))
          end do
        end do
      end do
      ! the top layer: blended toward the full jump (diff_top_km)
      if (diff_top_km > 0.0_dp) then
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
      if (diff_top_km > 0.0_dp) then
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
      ! the magnetic field and psi: -c J(q) (no density weight)
      do v = nVarCz + 1, nv
        do k = 1 - dk, nb(3)
          do j = 1 - dj, nb(2)
            do i = 1 - di, nb(1)
              F(i, j, k, v) = -b_diff_factor*C(i, j, k)*jump4(State_GVB(i-di, j-dj, k-dk, v, ib), State_GVB(i, j, k, v, ib), &
                  State_GVB(i+di, j+dj, k+dk, v, ib), State_GVB(i+2*di, j+2*dj, k+2*dk, v, ib))
            end do
          end do
        end do
      end do
    end subroutine faces_b

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

  ! ================================================================ energy budget

  ! h0(r) and Phi(r) of the background, the overlap weights, the shells'
  ! energy at the start, the output files. Call after cz_setup at the run's
  ! start (also after a restart: the budget counts from there).
  subroutine budget_start(t)
    real(dp), intent(in) :: t
    real(dp) :: th2, ph2
    integer :: j, k, jb, ib, p
    if (dt_budget <= 0.0_dp) return
    if (rsst_h2015) call arbor_abort('cz3d: the energy budget is built for rsst_form = hotta2012')
    budget_on = .true.
    allocate(wyy_B(nb(2), nb(3), nBlockMax), acc_IV(1:n_r, nAcc), E0_I(1:n_r), bgf_r(1:n_r))
    allocate(work_CB(nb(1), nb(2), nb(3), nBlockMax), workh_CB(nb(1), nb(2), nb(3), nBlockMax))
    work_CB = 0.0_dp
    workh_CB = 0.0_dp
    bgf_r = ddr_c(p0) + rho0(1:n_r)*g_c(1:n_r)
    ! weight: solid angle; 1/2 where the cell centre is in both patches'
    ! cores (theta in [pi/4, 3pi/4], phi in [-3pi/4, 3pi/4]), 0 in a margin
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      p = Key_IB(KP, ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          wyy_B(j, k, ib) = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          if (.not. in_core(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))) then
            wyy_B(j, k, ib) = 0.0_dp
          else
            call yy_angles(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k), th2, ph2)
            if (in_core(th2, ph2)) wyy_B(j, k, ib) = 0.5_dp*wyy_B(j, k, ib)
          end if
        end do
      end do
    end do
    acc_IV = 0.0_dp
    call shell_energy(E0_I)
    t_bud0 = t
    t_next_bud = t + dt_budget
    if (is_root .and. restarted) then
      ! a continued run appends (its records count from the restart time t0)
      open(newunit=u_bud, file=trim(out_dir)//'/budget_shells.dat', status='unknown', position='append', action='write')
      write(u_bud, '(A,ES16.8)') '# continued from a restart file; changes below count from t0 = ', t
      open(newunit=u_bsum, file=trim(out_dir)//'/budget_log.txt', status='unknown', position='append', action='write')
      write(u_bsum, '(A,ES16.8,A)') '# continued from a restart file: mean rates since t0 = ', t, ' s'
      open(newunit=u_bkh, file=trim(out_dir)//'/budget_kh.dat', status='unknown', position='append', action='write')
      open(newunit=u_bsp, file=trim(out_dir)//'/budget_split.dat', status='unknown', position='append', action='write')
      write(u_bsp, '(A,ES16.8)') '# continued from a restart file; rates below count from t0 = ', t
    else if (is_root) then
      open(newunit=u_bud, file=trim(out_dir)//'/budget_shells.dat', status='replace', action='write')
      write(u_bud, '(A)') '# energy budget per radial shell [erg], changes since the start (t0 below): records of'
      write(u_bud, '(A)') '# i  z[km]  dE_total  dE_hydro  dE_radiation  dE_cooling  dE_diffusion  E_rsst_source  ' // &
          'dE_total-(sum of operators)  W_background_force  E_rsst_source_fluct   (the RSST source and the ' // &
          'background-force work are part of dE_hydro)'
      open(newunit=u_bsum, file=trim(out_dir)//'/budget_log.txt', status='replace', action='write')
      write(u_bsum, '(A,ES16.8,A)') '# energy budget of the whole domain since t0 = ', t, ' s, as mean rates in Lsun:'
      write(u_bsum, '(A)') '# step  t  dE_total/dt  hydro  radiation  cooling  diffusion  rsst_source  residual  ' // &
          'bg_force_work  rsst_fluct  hydro-rsst-bg_work (scheme error)'
      open(newunit=u_bkh, file=trim(out_dir)//'/budget_kh.dat', status='replace', action='write')
      write(u_bkh, '(A)') '# horizontal kinetic energy K_h = rho v_h^2/2 per radial shell, mean rates in Lsun since t0 ' // &
          '(records): i  z[km]  dK_h(hydro)  W_h = -v_h.grad_h p1  F_Kh(c) = rho v_r v_h^2/2 over the sphere  ' // &
          'dK_h(diffusion)  dK(diffusion)  dK(hydro)  W  F_K(c)'
      open(newunit=u_bsp, file=trim(out_dir)//'/budget_split.dat', status='replace', action='write')
      write(u_bsp, '(A)') '# hydro operator split into kinetic and thermal energy per radial shell, mean rates in ' // &
          'Lsun since t0 (records): i  z[km]  dE_kin  dE_th  W  rsst_kin  rsst_total  bg_force_work  L_res_total(c)  ' // &
          'L_res_kin(c)   (dE_th = rho T ds + (h + Phi) drho; W = -v.grad p1 - g rho1 v_r; L_res at the cell centre: ' // &
          'rho v_r (h + v^2/2 + Phi) and rho v_r v^2/2 over the sphere; see scripts/analyze_split.py)'
    end if
  contains
    logical function in_core(th, ph)
      real(dp), intent(in) :: th, ph
      in_core = th >= 0.25_dp*pi .and. th <= 0.75_dp*pi .and. abs(ph) <= 0.75_dp*pi
    end function in_core
  end subroutine budget_start

  ! Background enthalpy h0(r) (EOS core) and potential Phi(r) (Phi(1) = 0)
  subroutine background_h_phi()
    type(eos_state) :: st
    real(dp) :: out(n_field)
    integer :: i
    allocate(h0_r_g(-1:n_r+2), phi_r_g(-1:n_r+2))
    h0_r => h0_r_g
    phi_r => phi_r_g
    do i = -1, n_r + 2
      call eos_eval(tab, rho0(i), s0(i), out)
      call eos_get_from_rho_T(rho0(i), exp(out(F_LNT)), st, tab%comp, tab%coulomb)
      h0_r(i) = st%e + st%p/st%rho
    end do
    phi_r(1) = 0.0_dp
    do i = 1, n_r + 1
      phi_r(i+1) = phi_r(i) + 0.5_dp*(g_c(i) + g_c(i+1))*(r_c(i+1) - r_c(i))
    end do
    do i = 0, -1, -1
      phi_r(i) = phi_r(i+1) - 0.5_dp*(g_c(i) + g_c(i+1))*(r_c(i+1) - r_c(i))
    end do
  end subroutine background_h_phi

  ! E per shell (e from the EOS core at the cell's rho, T)
  subroutine shell_energy(E_I)
    real(dp), intent(out) :: E_I(1:n_r)
    type(eos_state) :: st
    real(dp) :: out(n_field), rho, w
    integer :: jb, ib, i, j, k
    E_I = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          w = wyy_B(j, k, ib)
          if (w == 0.0_dp) cycle
          do i = 1, nb(1)
            rho = rho0(i) + State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            call eos_get_from_rho_T(rho, exp(out(F_LNT)), st, tab%comp, tab%coulomb)
            E_I(r_off + i) = E_I(r_off + i) + w*Metric_DB(1, ib)%Vol_I(i)*(rho*(st%e + phi_r(i)) &
                + 0.5_dp*rho*sum(State_GVB(i, j, k, vr_:vp_, ib)**2))
          end do
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(E_I)
  end subroutine shell_energy

  ! Energy change of the RK step of block ib (StateOld -> State) and the
  ! RSST source from the step's drho
  subroutine budget_hydro(ib, dt)
    integer, intent(in) :: ib
    real(dp), intent(in) :: dt
    real(dp) :: out(n_field), rho_o, rho_n, T_o, h_o, w, vol, ke_o, ke_n
    integer :: i, j, k
    do k = 1, nb(3)
      do j = 1, nb(2)
        w = wyy_B(j, k, ib)
        if (w == 0.0_dp) cycle
        do i = 1, nb(1)
          rho_o = rho0(i) + StateOld_GVB(i, j, k, rho1_, ib)
          rho_n = rho0(i) + State_GVB(i, j, k, rho1_, ib)
          call eos_eval(tab, rho_o, s0(i) + StateOld_GVB(i, j, k, s1_, ib), out)
          T_o = exp(out(F_LNT))
          h_o = h0_r(i) + T_o*StateOld_GVB(i, j, k, s1_, ib) + p0(i)*expm1(out(F_LNP) - lnp0(i))/rho_o
          ke_o = 0.5_dp*sum(StateOld_GVB(i, j, k, vr_:vp_, ib)**2)
          ke_n = 0.5_dp*rho_n*sum(State_GVB(i, j, k, vr_:vp_, ib)**2)
          vol = w*Metric_DB(1, ib)%Vol_I(i)
          acc_IV(i, 1) = acc_IV(i, 1) + vol*(rho_o*T_o*(State_GVB(i, j, k, s1_, ib) - StateOld_GVB(i, j, k, s1_, ib)) &
              + (h_o + phi_r(i))*(rho_n - rho_o) + ke_n - rho_o*ke_o)
          ! RSST source over the step: drho = -dt div(rho v)/xi^2 (the RK
          ! step's actual change), so S dt = -(xi^2 - 1)(h + v^2/2 + Phi) drho
          acc_IV(i, 5) = acc_IV(i, 5) - vol*(xi(i)**2 - 1.0_dp)*(h_o + ke_o + phi_r(i))*(rho_n - rho_o)
          acc_IV(i, 7) = acc_IV(i, 7) - vol*(xi(i)**2 - 1.0_dp)*(h_o - h0_r(i) + ke_o)*(rho_n - rho_o)
          acc_IV(i, 6) = acc_IV(i, 6) + vol*dt*StateOld_GVB(i, j, k, vr_, ib)*bgf_r(i)
          acc_IV(i, 8) = acc_IV(i, 8) + vol*(ke_n - rho_o*ke_o)
          acc_IV(i, 9) = acc_IV(i, 9) + vol*(rho_o*T_o*(State_GVB(i, j, k, s1_, ib) - StateOld_GVB(i, j, k, s1_, ib)) &
              + (h_o + phi_r(i))*(rho_n - rho_o))
          acc_IV(i, 10) = acc_IV(i, 10) + vol*dt*work_CB(i, j, k, ib)
          acc_IV(i, 11) = acc_IV(i, 11) - vol*(xi(i)**2 - 1.0_dp)*ke_o*(rho_n - rho_o)
          acc_IV(i, 12) = acc_IV(i, 12) + w*dt*r_c(i)**2*rho_o*StateOld_GVB(i, j, k, vr_, ib)*(h_o + ke_o + phi_r(i))
          acc_IV(i, 13) = acc_IV(i, 13) + w*dt*r_c(i)**2*rho_o*StateOld_GVB(i, j, k, vr_, ib)*ke_o
          acc_IV(i, 14) = acc_IV(i, 14) + vol*0.5_dp*(rho_n*sum(State_GVB(i, j, k, vt_:vp_, ib)**2) &
              - rho_o*sum(StateOld_GVB(i, j, k, vt_:vp_, ib)**2))
          acc_IV(i, 15) = acc_IV(i, 15) + vol*dt*workh_CB(i, j, k, ib)
          acc_IV(i, 16) = acc_IV(i, 16) + w*dt*r_c(i)**2*rho_o*StateOld_GVB(i, j, k, vr_, ib) &
              *0.5_dp*sum(StateOld_GVB(i, j, k, vt_:vp_, ib)**2)
        end do
      end do
    end do
  end subroutine budget_hydro

  ! One budget record: per shell, and the whole domain as mean rates
  subroutine budget_write(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    real(dp) :: E_I(1:n_r), A_IV(1:n_r, nAcc), dE(1:n_r), f
    integer :: i
    call shell_energy(E_I)
    A_IV = acc_IV
    call arbor_sum_inplace(A_IV)
    t_next_bud = t_next_bud + dt_budget
    if (.not. is_root) return
    dE = E_I - E0_I
    write(u_bud, '(A,I0,A,ES16.8,A,ES16.8)') '# step= ', iStep, ' t= ', t, ' t0= ', t_bud0
    do i = 1, n_r
      write(u_bud, '(I5,F12.1,9ES16.7E3)') i, z_c(i)/1.0e5_dp, dE(i), A_IV(i, 1:5), dE(i) - sum(A_IV(i, 1:4)), &
          A_IV(i, 6:7)
    end do
    write(u_bud, '(A)') ''
    flush(u_bud)
    f = 1.0_dp/(max(t - t_bud0, 1.0e-30_dp)*Lsun)
    write(u_bsum, '(I9,ES16.8,10F12.5)') iStep, t, sum(dE)*f, (sum(A_IV(:, i))*f, i = 1, 5), &
        (sum(dE) - sum(A_IV(:, 1:4)))*f, sum(A_IV(:, 6))*f, sum(A_IV(:, 7))*f, &
        (sum(A_IV(:, 1)) - sum(A_IV(:, 5)) - sum(A_IV(:, 6)))*f
    flush(u_bsum)
    write(u_bsp, '(A,I0,A,ES16.8,A,ES16.8)') '# step= ', iStep, ' t= ', t, ' t0= ', t_bud0
    do i = 1, n_r
      write(u_bsp, '(I5,F12.1,8ES16.7E3)') i, z_c(i)/1.0e5_dp, A_IV(i, 8)*f, A_IV(i, 9)*f, A_IV(i, 10)*f, &
          A_IV(i, 11)*f, A_IV(i, 5)*f, A_IV(i, 6)*f, A_IV(i, 12)*f, A_IV(i, 13)*f
    end do
    write(u_bsp, '(A)') ''
    flush(u_bsp)
    write(u_bkh, '(A,I0,A,ES16.8,A,ES16.8)') '# step= ', iStep, ' t= ', t, ' t0= ', t_bud0
    do i = 1, n_r
      write(u_bkh, '(I5,F12.1,8ES16.7E3)') i, z_c(i)/1.0e5_dp, A_IV(i, 14)*f, A_IV(i, 15)*f, A_IV(i, 16)*f, &
          A_IV(i, 17)*f, A_IV(i, 18)*f, A_IV(i, 8)*f, A_IV(i, 10)*f, A_IV(i, 13)*f
    end do
    write(u_bkh, '(A)') ''
    flush(u_bkh)
  end subroutine budget_write

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

  ! ================================================================ diagnostics

  subroutine cz_after_step(self, iStep, t)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    if (probe_radius_deg > 0.0_dp) call probe_write(iStep, t)
    call nan_check(iStep, t)
    if (lz_correct_steps > 0 .or. lz_correct_dt > 0.0_dp) then
      if (lz_correct_steps > 0) then
        if (mod(iStep, lz_correct_steps) == 0) call lz_correct(t)
      end if
      if (lz_correct_dt > 0.0_dp .and. t >= t_next_lzc) call lz_correct(t)
    end if
    if (mod(iStep, n_log) == 0 .or. t >= t_end*(1.0_dp - 1.0e-12_dp)) call cz_log(iStep, t)
    if (prof_every > 0) then
      if (mod(iStep, prof_every) == 0) call cz_profile(iStep, t)
    else if (dt_prof > 0.0_dp .and. t >= t_next_prof - 1.0e-9_dp*dt_prof) then
      call cz_profile(iStep, t)
    end if
    if (budget_on .and. (t >= t_next_bud - 1.0e-9_dp*dt_budget .or. t >= t_end*(1.0_dp - 1.0e-12_dp))) &
        call budget_write(iStep, t)
  end subroutine cz_after_step

  ! Debugging probe (probe_radius_deg > 0): the top probe_cells cells of the
  ! columns near (probe_lat_deg, probe_lon_deg), every step, one file per rank:
  ! step t i patch level j k lat lon rho1/rho0 s1 v_r v_theta v_phi T c_s Mach
  subroutine probe_write(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    integer, save :: u = 0
    real(dp) :: x0(3), x(3), e(3,3), la, lo, out(n_field), rho, cs, Tc
    integer :: jb, ib, i, j, k
    character(len=300) :: name
    la = probe_lat_deg*pi/180.0_dp;  lo = probe_lon_deg*pi/180.0_dp
    x0 = [cos(la)*cos(lo), cos(la)*sin(lo), sin(la)]
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (.not. arbor_physical_side(ib, 1, 2)) cycle          ! top blocks only
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
          x = e(:, 1)
          if (Key_IB(KP, ib) == 2) x = yy_cart(x)
          if (dot_product(x, x0) < cos(probe_radius_deg*pi/180.0_dp)) cycle
          if (u == 0) then
            write(name, '(A,A,I4.4,A)') trim(out_dir), '/probe_r', my_rank, '.dat'
            open(newunit=u, file=trim(name), status='replace', action='write')
            write(u, '(A)') '# step t i patch level j k lat lon rho1/rho0 s1 v_r v_theta v_phi T c_s Mach'
          end if
          do i = max(1, nb(1) - probe_cells + 1), nb(1)
            rho = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            Tc = exp(out(F_LNT));  cs = sqrt(out(F_GAMMA1)*exp(out(F_LNP))/rho)
            write(u, '(I9,ES15.7,I5,4I3,2F9.3,9ES14.5)') iStep, t, r_off + i, Key_IB(KP, ib), Key_IB(KL+2, ib), j, k, &
                asin(max(-1.0_dp, min(1.0_dp, x(3))))*180.0_dp/pi, atan2(x(2), x(1))*180.0_dp/pi, &
                State_GVB(i, j, k, rho1_, ib)/rho0(i), State_GVB(i, j, k, s1_, ib), State_GVB(i, j, k, vr_:vp_, ib), &
                Tc, cs, norm2(State_GVB(i, j, k, vr_:vp_, ib))/cs
          end do
        end do
      end do
      call cz_unbind()
    end do
    if (u /= 0) flush(u)
  end subroutine probe_write

  ! Stop with a profile if the state holds a NaN
  subroutine nan_check(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    logical :: bad
    integer :: jb, ib, i, j, k, v
    real(dp) :: e(3,3), x(3)
    character(len=300) :: msg
    bad = .false.
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (.not. any(State_GVB(1:nb(1), 1:nb(2), 1:nb(3), :, ib) /= State_GVB(1:nb(1), 1:nb(2), 1:nb(3), :, ib))) cycle
      bad = .true.
      ! where: the first NaN cell of the block (radial index, depth, global lat/lon, level)
      call cz_bind(ib)
      cells: do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            do v = 1, nv
              if (State_GVB(i, j, k, v, ib) /= State_GVB(i, j, k, v, ib)) exit cells
            end do
          end do
        end do
      end do cells
      e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
      x = e(:, 1)
      if (Key_IB(KP, ib) == 2) x = yy_cart(x)
      write(msg, '(A,I0,A,I0,A,F9.1,A,F7.2,A,F8.2,A,I0,A,I0,A,2I4)') 'cz3d: NaN on rank ', my_rank, ': cell ', r_off + i, &
          ' (z ', z_c(i)/1.0e5_dp, ' km), lat ', asin(max(-1.0_dp, min(1.0_dp, x(3))))*180.0_dp/pi, ', lon ', &
          atan2(x(2), x(1))*180.0_dp/pi, ', patch ', Key_IB(KP, ib), ', level ', Key_IB(KL+2, ib), ', j k in block', j, k
      write(0, '(A)') trim(msg)
      call cz_unbind()
    end do
    if (.not. arbor_any(bad)) return
    call cz_log(iStep, t)
    call arbor_abort('cz3d: NaN in the state')
  end subroutine nan_check

  ! Horizontal sums per radial cell (weighted by solid angle; the Yin-Yang
  ! overlap counts twice, which does not matter for means of a uniform
  ! state): returns Sum_IV(i, v) for v = 1 weight, 2 rho1, 3 (rho1 - <rho1>)^2,
  ! 4 s1, 5 (s1 - <s1>)^2, 6 v_r, 7 v_r^2, 8 v_h^2, 9 T, 10 L_rad and
  ! 11 L_conv at face i, 12 rho v_r h1, 13 rho v_r, 14 h1 (h1 = T s1 + p1/rho,
  ! the enthalpy perturbation to first order), 15 rho v^2 v_r / 2
  ! 16 the diffusion's energy flux at face i, 17 rho |v_r| (the up/down mass exchange),
  ! 18 B_r^2, 19 B_theta^2 + B_phi^2 (use_B)
  ! (i = 0..n_r; the face values in row i, the cell values in rows 1..n_r).
  ! The deviations in a second pass, about the means of the first (a
  ! one-pass variance loses the digits). The resolved enthalpy flux is
  ! <rho v_r h1> - <rho v_r><h1> (columns 12-14).
  ! with_T is kept for the interface; T and p1 are always evaluated.
  subroutine horizontal_sums(Sum_IV, with_T, only_level)
    real(dp), intent(out) :: Sum_IV(0:n_r, 19)
    logical, intent(in) :: with_T
    integer, intent(in), optional :: only_level     ! only the blocks of this horizontal level (a refined patch)
    real(dp) :: w, out(n_field), rho, T, p1, h1
    real(dp), allocatable :: Dev_IV(:,:)
    integer :: jb, ib, i, j, k
    Sum_IV = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (present(only_level)) then
        if (Key_IB(KL+2, ib) /= only_level) cycle
      end if
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          ! (the faces between radial blocks from the block above)
          do i = merge(0, 1, at_bottom), nb(1)
            Sum_IV(r_off + i, 10) = Sum_IV(r_off + i, 10) + w*qrad_FB(i, j, k, ib)
            if (qconv_face) then
              Sum_IV(r_off + i, 11) = Sum_IV(r_off + i, 11) + w*FaceState_FVDB(i, j, k, 1, 1, ib)
            else
              Sum_IV(r_off + i, 11) = Sum_IV(r_off + i, 11) + w*qconv_FB(i, j, k, ib)
            end if
            Sum_IV(r_off + i, 16) = Sum_IV(r_off + i, 16) + w*qdif_FB(i, j, k, ib)
          end do
          do i = 1, nb(1)
            associate (s => State_GVB(i, j, k, :, ib), Sm => Sum_IV(r_off + i, :))
              Sm(1) = Sm(1) + w
              Sm(2) = Sm(2) + w*s(rho1_)
              Sm(4) = Sm(4) + w*s(s1_)
              Sm(6) = Sm(6) + w*s(vr_)
              Sm(7) = Sm(7) + w*s(vr_)**2
              Sm(8) = Sm(8) + w*(s(vt_)**2 + s(vp_)**2)
              rho = rho0(i) + s(rho1_)                                    ! inertia (mass flux) density
              call eos_eval(tab, rho0(i) + rq(i)*s(rho1_), s0(i) + s(s1_), out)
              T = exp(out(F_LNT))
              p1 = p0(i)*expm1(out(F_LNP) - lnp0(i))
              h1 = T*s(s1_) + p1/(rho0(i) + rq(i)*s(rho1_))
              Sm(9) = Sm(9) + w*T
              Sm(12) = Sm(12) + w*rho*s(vr_)*h1
              Sm(13) = Sm(13) + w*rho*s(vr_)
              Sm(14) = Sm(14) + w*h1
              Sm(15) = Sm(15) + w*0.5_dp*rho*sum(s(vr_:vp_)**2)*s(vr_)
              Sm(17) = Sm(17) + w*rho*abs(s(vr_))
              if (use_B) then
                Sm(18) = Sm(18) + w*s(br_)**2
                Sm(19) = Sm(19) + w*(s(bt_)**2 + s(bp_)**2)
              end if
            end associate
          end do
          if (at_bottom) Sum_IV(0, 1) = Sum_IV(0, 1) + w
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(Sum_IV)
    allocate(Dev_IV(0:n_r, 2))
    Dev_IV = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (present(only_level)) then
        if (Key_IB(KL+2, ib) /= only_level) cycle
      end if
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          do i = 1, nb(1)
            Dev_IV(r_off + i, 1) = Dev_IV(r_off + i, 1) + w*(State_GVB(i, j, k, rho1_, ib) - Sum_IV(r_off + i, 2)/Sum_IV(r_off + i, 1))**2
            Dev_IV(r_off + i, 2) = Dev_IV(r_off + i, 2) + w*(State_GVB(i, j, k, s1_, ib) - Sum_IV(r_off + i, 4)/Sum_IV(r_off + i, 1))**2
          end do
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(Dev_IV)
    Sum_IV(:, 3) = Dev_IV(:, 1)
    Sum_IV(:, 5) = Dev_IV(:, 2)
  end subroutine horizontal_sums

  ! One log line: max Mach, max |v_r|, max |v_h|, the horizontally averaged
  ! luminosity (as CZTR: at the top face and its range below), and the
  ! largest horizontal rms deviations of s1 and rho1/rho0 (0 for a uniform state)
  subroutine cz_log(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    real(dp) :: S(0:n_r, 19), L(0:n_r), Lres(1:n_r), mach, mach_r, dmin, vrmax, vhmax, cs, ds_h, dr_h, out(n_field), rho
    real(dp) :: bmax, emag, ekin, eratio, divb_rel, sd, sb, dxm, wl, wi, lz, tcor, vamax, n_cap, n_all, trmin, trmax
    integer :: nt                                   ! the top of the CZ (L_top, the L range): the join with the atmosphere
    integer :: jb, ib, i, j, k
    character(len=300) :: line
    mach = 0.0_dp;  mach_r = 0.0_dp;  vrmax = 0.0_dp;  vhmax = 0.0_dp;  dmin = huge(1.0_dp);  bmax = 0.0_dp
    vamax = 0.0_dp;  n_cap = 0.0_dp;  n_all = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            rho = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            cs = sqrt(out(F_GAMMA1)*exp(out(F_LNP))/rho)
            mach = max(mach, norm2(State_GVB(i, j, k, vr_:vp_, ib))/cs)
            mach_r = max(mach_r, xi(i)*norm2(State_GVB(i, j, k, vr_:vp_, ib))/cs)
            dmin = min(dmin, State_GVB(i, j, k, rho1_, ib)/rho0(i))
            vrmax = max(vrmax, abs(State_GVB(i, j, k, vr_, ib)))
            vhmax = max(vhmax, hypot(State_GVB(i, j, k, vt_, ib), State_GVB(i, j, k, vp_, ib)))
            if (use_B) then
              bmax = max(bmax, norm2(State_GVB(i, j, k, br_:bp_, ib)))
              vamax = max(vamax, norm2(State_GVB(i, j, k, br_:bp_, ib))/sqrt(four_pi*rho))
              if (lorentz_f(sum(State_GVB(i, j, k, br_:bp_, ib)**2), rho, cell_cap(i, j, k, ib, cs)) < 0.9_dp) &
                  n_cap = n_cap + 1.0_dp
              n_all = n_all + 1.0_dp
            end if
          end do
        end do
      end do
    end do
    call cz_unbind()
    mach = arbor_max(mach);  mach_r = arbor_max(mach_r);  dmin = arbor_min(dmin);  vrmax = arbor_max(vrmax);  vhmax = arbor_max(vhmax)
    bmax = arbor_max(bmax);  vamax = arbor_max(vamax);  n_cap = arbor_sum(n_cap);  n_all = arbor_sum(n_all)
    nt = merge(k_join, n_r, use_atmosphere)
    call angular_momentum(lz, tcor)
    ! the transition region's height range (T = 1e5 K crossings below 10 Mm; -1 without)
    trmin = huge(1.0_dp);  trmax = -huge(1.0_dp)
    if (use_atmosphere) then
      do jb = 1, nBlock
        call cz_bind(iBlockList_I(jb))
        call atm_tr_heights(iBlockList_I(jb), 1.0e4_dp, trmin, trmax)
      end do
      call cz_unbind()
    end if
    trmin = arbor_min(trmin);  trmax = arbor_max(trmax)
    if (trmin > 1.0e30_dp) trmin = -1.0_dp
    if (trmax < -1.0e30_dp) trmax = -1.0_dp
    divb_rel = 0.0_dp;  wl = 0.0_dp;  wi = 0.0_dp
    if (use_B) then
      ! rms(div B dx_min) / rms|B| over all cells (ghosts filled first)
      call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.true.)
      sd = 0.0_dp;  sb = 0.0_dp
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        call div_block(State_GVB(:, :, :, br_:bp_, ib), ib, divB_C)
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, nb(1)
              dxm = min(drdxi_c(i), r_c(i)*Metric_DB(2, ib)%Dx_I(j), r_c(i)*sin(Metric_DB(2, ib)%X_I(j))*Metric_DB(3, ib)%Dx_I(k))
              sd = sd + (divB_C(i, j, k)*dxm)**2
              sb = sb + sum(State_GVB(i, j, k, br_:bp_, ib)**2)
            end do
          end do
        end do
      end do
      call cz_unbind()
      sd = arbor_sum(sd);  sb = arbor_sum(sb)
      divb_rel = sqrt(sd/max(sb, tiny(1.0_dp)))
      call mhd_exchange(wl, wi)
    end if
    call horizontal_sums(S, .false.)
    ! magnetic over kinetic energy (the shells' horizontal means; the kinetic with rho0)
    emag = 0.0_dp;  ekin = 0.0_dp
    do i = 1, n_r
      emag = emag + vol_c(i)*(S(i, 18) + S(i, 19))/S(i, 1)/(2.0_dp*four_pi)
      ekin = ekin + vol_c(i)*0.5_dp*rho0(i)*(S(i, 7) + S(i, 8))/S(i, 1)
    end do
    eratio = emag/max(ekin, tiny(1.0_dp))
    ! radiation, MLT, diffusion and the cooling layer at the faces
    L = 4.0_dp*pi*area_f*(S(:, 10) + S(:, 11) + S(:, 16) + S(0, 1)*cool_f)/S(0, 1)/Lsun
    call resolved_L(S, Lres)
    ds_h = 0.0_dp;  dr_h = 0.0_dp
    do i = 1, n_r
      ds_h = max(ds_h, sqrt(S(i, 5)/S(i, 1)))
      dr_h = max(dr_h, sqrt(S(i, 3)/S(i, 1))/rho0(i))
    end do
    if (.not. is_root) return
    if (u_log == 0 .and. restarted) then
      open(newunit=u_log, file=trim(out_dir)//'/log.txt', status='unknown', position='append', action='write')
      write(u_log, '(A,I0)') '# continued from a restart file at step ', iStep
      wall0 = arbor_wtime()
    else if (u_log == 0) then
      open(newunit=u_log, file=trim(out_dir)//'/log.txt', status='replace', action='write')
      write(u_log, '(A)') '# step  t[s]  dt[s]  max_Mach  max|v_r|  max|v_h| [cm/s]  L_top/Lsun  min_L  max_L  ' // &
          'max_rms_h(s1)[erg/g/K]  max_rms_h(rho1)/rho0  wall[s]  max_Mach(reduced c_s)  ' // &
          'max|L_enth+L_kin|  min_L_total  max_L_total   (L_total = L_rad + L_conv + L_diff + L_cool + L_enth + L_kin)' // &
          '  min(rho1/rho0)  max_tsc  t_local_bottom[s]  t_local_mid[s]  max|B|[G]  E_mag/E_kin  rms(divB dx)/rms|B|' // &
          '  W_Lorentz[Lsun]  dE_mag/dt(induction)[Lsun] (their sum is ~0 for walls without Poynting flux)' // &
          '  L_z[g cm2/s] (angular momentum of the flow about the global z axis, in the rotating frame)' // &
          '  T_Coriolis[g cm2/s2] (-2 Omega int rho s v_s dV: the Coriolis torque on L_z)' // &
          '  max_v_A[cm/s] (uncapped)  cap_frac (fraction of cells with the Lorentz factor f < 0.9, see va_max)' // &
          '  TR_min[km]  TR_max[km] (heights of the T = 1e5 K crossings, use_atmosphere)'
      wall0 = arbor_wtime()
    end if
    write(u_log, '(I9,ES14.6,ES11.3,3ES11.3,3F10.5,2ES11.3,F10.1,ES11.3,3F10.5,ES11.3,F9.3,2ES14.6,3ES11.3,2ES12.4,2ES14.6,2ES11.3,2F10.1)') iStep, t, &
        dt_last, mach, vrmax, vhmax, L(nt), minval(L(1:nt-1)), maxval(L(1:nt-1)), ds_h, dr_h, arbor_wtime() - wall0, &
        mach_r, maxval(abs(Lres(1:nt))), minval(L(1:nt-1) + Lres(1:nt-1)), maxval(L(1:nt-1) + Lres(1:nt-1)), dmin, maxval(tsc), &
        tloc(1), tloc(n_r/2), bmax, eratio, divb_rel, wl/Lsun, wi/Lsun, lz, tcor, vamax, n_cap/max(n_all, 1.0_dp), trmin, trmax
    flush(u_log)
    write(line, '(A,I8,A,ES11.4,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,F8.5,A,F8.5,A,F8.5,A,ES9.2)') 'step ', iStep, &
        ' t=', t, ' dt=', dt_last, ' Mach=', mach, ' (red. ', mach_r, ') |vr|=', vrmax, ' |vh|=', vhmax, ' L_top=', &
        L(nt), ' L ', minval(L(1:nt-1) + Lres(1:nt-1)), '..', maxval(L(1:nt-1) + Lres(1:nt-1)), &
        ' rms_h(s1)=', ds_h
    call arbor_print(line)
  end subroutine cz_log

  ! Horizontally averaged radial profile, out_dir/profile_NNNN.dat
  subroutine cz_profile(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    real(dp) :: S(0:n_r, 19)
    character(len=300) :: name
    call horizontal_sums(S, .true.)
    write(name, '(A,A,I4.4,A)') trim(out_dir), '/profile_', i_prof, '.dat'
    call write_profile(S, name, 'horizontal means (solid-angle weighted) and rms deviations')
    if (patch_level > refine_h) then
      ! the same over the refined patch alone (its finest blocks)
      call horizontal_sums(S, .true., patch_level)
      write(name, '(A,A,I4.4,A)') trim(out_dir), '/profile_patch_', i_prof, '.dat'
      call write_profile(S, name, 'means over the refined patch (blocks of level patch_level) and rms deviations')
    end if
    if (omega_rot /= 0.0_dp) call rotation_profile(iStep, t)
    i_prof = i_prof + 1
    t_next_prof = t_prof0 + i_prof*dt_prof
  contains
    subroutine write_profile(S, name, what)
      real(dp), intent(in) :: S(0:n_r, 19)
      character(*), intent(in) :: name, what
      real(dp) :: w, m, Lres(1:n_r)
      integer :: u, i
      call resolved_L(S, Lres)
      if (.not. is_root) return
      open(newunit=u, file=trim(name), status='replace', action='write')
      write(u, '(A,ES14.6,A,I0,A,A)') '# t = ', t, ' s; step ', iStep, '; ', what
      write(u, '(A)') '# i  r[cm]  z[km]  <rho>  <T>  <s1>  rms_h(s1)  <v_r>  rms(v_r)  rms(v_h)  ' // &
          'L_rad/Lsun  L_conv/Lsun  L_enth+L_kin/Lsun  L_kin/Lsun  L_diff/Lsun  L_cool/Lsun   (L_rad, L_conv, ' // &
          'L_diff (numerical diffusion), L_cool (cooling layer): through the face above the cell; L_enth, L_kin: ' // &
          'resolved, at the cell centre)  t_local[s]  tsc (the shell''s own clock and step factor; dt1 > 0)  ' // &
          '<rho v_r>  <rho |v_r|> [g/cm2/s] (the mean mass flux, the up/down exchange)  rms(B_r)  rms(B_h) [G]'
      do i = 1, n_r
        w = S(i, 1)
        m = S(i, 4)/w
        write(u, '(I5,ES20.12,20ES16.7E3)') i, r_c(i), z_c(i)/1.0e5_dp, rho0(i) + rq(i)*S(i, 2)/w, S(i, 9)/w, m, &
            sqrt(S(i, 5)/w), S(i, 6)/w, sqrt(S(i, 7)/w), sqrt(S(i, 8)/w), &
            4.0_dp*pi*area_f(i)*S(i, 10)/S(0, 1)/Lsun, 4.0_dp*pi*area_f(i)*S(i, 11)/S(0, 1)/Lsun, Lres(i), &
            4.0_dp*pi*r_c(i)**2*S(i, 15)/w/Lsun, 4.0_dp*pi*area_f(i)*S(i, 16)/S(0, 1)/Lsun, &
            4.0_dp*pi*area_f(i)*cool_f(i)/Lsun, tloc(i), tsc(i), S(i, 13)/w, S(i, 17)/w, sqrt(S(i, 18)/w), sqrt(S(i, 19)/w)
      end do
      close(u)
    end subroutine write_profile
  end subroutine cz_profile

  ! Resolved luminosity at the cell centres: 4 pi r^2 (<rho v_r h1> -
  ! <rho v_r><h1> + <rho v^2 v_r/2>), in Lsun
  subroutine resolved_L(S, Lres)
    real(dp), intent(in) :: S(0:n_r, 19)
    real(dp), intent(out) :: Lres(1:n_r)
    integer :: i
    real(dp) :: w
    do i = 1, n_r
      w = S(i, 1)
      Lres(i) = 4.0_dp*pi*r_c(i)**2*(S(i, 12)/w - (S(i, 13)/w)*(S(i, 14)/w) + S(i, 15)/w)/Lsun
    end do
  end subroutine resolved_L

  ! Every cell of the radial indices dump_shells, one file per rank
  ! (out_dir/shell_dump_rNNNN.dat): i, patch, the position as a global unit
  ! vector (x, y, z; Yin = global frame, Yang mapped by (x,y,z) -> (-x,z,y)),
  ! the overlap weight (solid angle, 1/2 where both patches' cores cover the
  ! point, 0 in a margin), v_r, the horizontal velocity as a global Cartesian
  ! vector, rho1/rho0, s1
  subroutine dump_shells_write()
    real(dp) :: th, ph, th2, ph2, e(3,3), x(3), vh(3), w
    integer :: u, jb, ib, j, k, n, i
    character(len=300) :: name
    if (all(dump_shells <= 0)) return
    write(name, '(A,A,I4.4,A)') trim(out_dir), '/shell_dump_r', my_rank, '.dat'
    open(newunit=u, file=trim(name), status='replace', action='write')
    write(u, '(A)') '# i  patch  x  y  z  weight  v_r  vh_x  vh_y  vh_z [cm/s]  rho1/rho0  s1'
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          th = Metric_DB(2, ib)%X_I(j)
          ph = Metric_DB(3, ib)%X_I(k)
          e = sph_basis(th, ph)
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          if (.not. core(th, ph)) then
            w = 0.0_dp
          else
            call yy_angles(th, ph, th2, ph2)
            if (core(th2, ph2)) w = 0.5_dp*w
          end if
          do n = 1, size(dump_shells)
            i = dump_shells(n) - r_off
            if (i < 1 .or. i > nb(1)) cycle
            x = e(:, 1)
            vh = State_GVB(i, j, k, vt_, ib)*e(:, 2) + State_GVB(i, j, k, vp_, ib)*e(:, 3)
            if (Key_IB(KP, ib) == 2) then
              x = yy_cart(x)
              vh = yy_cart(vh)
            end if
            write(u, '(2I5,3F12.8,ES14.6,4ES14.6,2ES14.6)') r_off + i, Key_IB(KP, ib), x, w, State_GVB(i, j, k, vr_, ib), vh, &
                State_GVB(i, j, k, rho1_, ib)/rho0(i), State_GVB(i, j, k, s1_, ib)
          end do
        end do
      end do
    end do
    call cz_unbind()
    close(u)
  contains
    logical function core(t, p)
      real(dp), intent(in) :: t, p
      core = t >= 0.25_dp*pi .and. t <= 0.75_dp*pi .and. abs(p) <= 0.75_dp*pi
    end function core
  end subroutine dump_shells_write

  ! The profile counter and schedule go into restart files
  subroutine cz_restart_save(self, Data_I)
    class(cz3d_app_t), intent(inout) :: self
    real(dp), allocatable, intent(out) :: Data_I(:)
    Data_I = [real(i_prof, dp), t_next_prof, t_prof0, tloc]
  end subroutine cz_restart_save

  subroutine cz_restart_load(self, Data_I)
    class(cz3d_app_t), intent(inout) :: self
    real(dp), intent(in) :: Data_I(:)
    if (size(Data_I) /= 2 .and. size(Data_I) /= 3 .and. size(Data_I) /= 3 + n_r) &
        call arbor_abort('cz3d: restart file from another application')
    i_prof = nint(Data_I(1))
    t_next_prof = Data_I(2)
    t_prof0 = 0.0_dp
    if (size(Data_I) >= 3) t_prof0 = Data_I(3)           ! (files before 2026-10-11 hold two numbers)
    if (.not. allocated(tloc)) allocate(tloc(1:n_r))
    tloc = -1.0_dp                                        ! (older files: set to the run's time by the caller)
    if (size(Data_I) == 3 + n_r) tloc = Data_I(4:)
  end subroutine cz_restart_load

  ! Order-independent checksum of the interior state of all blocks, to
  ! compare runs (e.g. continued vs uninterrupted) bitwise
  integer(i8) function state_checksum() result(c)
    integer(i8) :: x
    integer :: jb, ib, i, j, k, v
    c = 0
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      x = int(iLeaf_B(ib), i8)
      do v = 1, nv
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, nb(1)
              x = ieor(x, transfer(State_GVB(i, j, k, v, ib), 1_i8))
              x = ieor(x, ishft(x, 13));  x = ieor(x, ishft(x, -7));  x = ieor(x, ishft(x, 17))
            end do
          end do
        end do
      end do
      c = c + iand(x, 2_i8**40 - 1)
    end do
    c = arbor_sum(c)
  end function state_checksum

  ! The cell of the largest Mach number |v|/c_s: printed (radius, global
  ! latitude and longitude, z, T, rho, v_r, v_h, c_s, |B|) by the rank that holds it
  subroutine max_mach_report(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    real(dp) :: out(n_field), rhop, cs, m, mbest, gbest, e(3,3), x(3), loc(8)
    integer :: jb, ib, i, j, k
    character(len=400) :: msg
    mbest = -1.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            rhop = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rhop, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            cs = sqrt(out(F_GAMMA1)*exp(out(F_LNP))/rhop)
            m = norm2(State_GVB(i, j, k, vr_:vp_, ib))/cs
            if (m > mbest) then
              mbest = m
              e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
              x = e(:, 1)
              if (Key_IB(KP, ib) == 2) x = yy_cart(x)
              loc = [r_c(i), asin(max(-1.0_dp, min(1.0_dp, x(3))))*180.0_dp/pi, atan2(x(2), x(1))*180.0_dp/pi, &
                  exp(out(F_LNT)), rhop, State_GVB(i, j, k, vr_, ib), hypot(State_GVB(i, j, k, vt_, ib), &
                  State_GVB(i, j, k, vp_, ib)), cs]
              if (use_B) then
                write(msg, '(A,F9.4,A,F7.2,A,F8.2,A,F10.1,A,ES10.3,A,ES10.3,A,2F9.2,A,F8.2,A,F9.2,A)') 'max Mach: r/R ', &
                    loc(1)/6.957e10_dp, '  lat ', loc(2), '  lon ', loc(3), '  z[km] ', z_c(i)/1.0e5_dp, '  T ', loc(4), &
                    '  rho ', loc(5), '  v_r, v_h [km/s] ', loc(6)/1.0e5_dp, loc(7)/1.0e5_dp, '  c_s ', loc(8)/1.0e5_dp, &
                    '  |B| [G] ', norm2(State_GVB(i, j, k, br_:bp_, ib)), ''
              else
                write(msg, '(A,F9.4,A,F7.2,A,F8.2,A,F10.1,A,ES10.3,A,ES10.3,A,2F9.2,A,F8.2)') 'max Mach: r/R ', &
                    loc(1)/6.957e10_dp, '  lat ', loc(2), '  lon ', loc(3), '  z[km] ', z_c(i)/1.0e5_dp, '  T ', loc(4), &
                    '  rho ', loc(5), '  v_r, v_h [km/s] ', loc(6)/1.0e5_dp, loc(7)/1.0e5_dp, '  c_s ', loc(8)/1.0e5_dp
              end if
            end if
          end do
        end do
      end do
    end do
    call cz_unbind()
    gbest = arbor_max(mbest)
    if (mbest == gbest) then
      write(0, '(A,I0,A,ES12.5,A,F8.2,A,A)') 'step ', iStep, '  t ', t, '  Mach ', mbest, '  ', trim(msg)
      flush(0)
    end if
  end subroutine max_mach_report

  ! Tecplot file of the state in out_dir: rho1, v_x, v_y, v_z, s1, [B_x, B_y,
  ! B_z,] T - T0, rho, T and the Mach number |v|/c_s (cell centres; ARBOR
  ! writes node values). The vectors as Cartesian components in the global
  ! (Yin) frame, the frame of the X, Y, Z coordinates (Yang's rotated into it).
  ! The plot quantities sit in slots beyond 8 and only those are passed to the
  ! writer: its ghost fill then treats them as scalars (it rotates the slots
  ! declared as spherical vectors, v at vr_ and B at br_, which hold nothing here)
  subroutine cz_write_plot(self, dom, iStep, t)
    class(cz3d_app_t), intent(inout) :: self
    type(arbor_domain_t), intent(in) :: dom
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    integer, parameter :: o = 8
    character(len=300) :: name
    real(dp), allocatable :: P_GVB(:,:,:,:,:)
    real(dp) :: out(n_field), rhop, e(3,3), vc(3), bc(3)
    integer :: jb, ib, i, j, k, n, q
    n = merge(12, 9, use_B)
    allocate(P_GVB(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng, o + n, nBlockMax))
    P_GVB = 0.0_dp
    q = merge(3, 0, use_B)                                ! the scalars after the field
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
          do i = 1, nb(1)
            vc = matmul(e, State_GVB(i, j, k, vr_:vp_, ib))  ! Cartesian in the patch frame
            if (Key_IB(KP, ib) == 2) vc = yy_cart(vc)       ! Yang -> the global frame
            P_GVB(i, j, k, o+1, ib) = State_GVB(i, j, k, rho1_, ib)
            P_GVB(i, j, k, o+2:o+4, ib) = vc
            P_GVB(i, j, k, o+5, ib) = State_GVB(i, j, k, s1_, ib)
            if (use_B) then
              bc = matmul(e, State_GVB(i, j, k, br_:bp_, ib))
              if (Key_IB(KP, ib) == 2) bc = yy_cart(bc)
              P_GVB(i, j, k, o+6:o+8, ib) = bc
            end if
            rhop = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rhop, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            P_GVB(i, j, k, o+q+6, ib) = exp(out(F_LNT)) - T0(i)
            P_GVB(i, j, k, o+q+7, ib) = rhop
            P_GVB(i, j, k, o+q+8, ib) = exp(out(F_LNT))
            P_GVB(i, j, k, o+q+9, ib) = norm2(State_GVB(i, j, k, vr_:vp_, ib))/sqrt(out(F_GAMMA1)*exp(out(F_LNP))/rhop)
          end do
        end do
      end do
    end do
    call cz_unbind()
    write(name, '(A,A,I8.8,A)') trim(out_dir), '/plot_', iStep, '.plt'
    if (use_B) then
      call arbor_plot_tecplot(dom, name, P_GVB, iVarMin=o+1, iVarMax=o+n, names=['rho1', 'vx  ', 'vy  ', 'vz  ', &
          's1  ', 'Bx  ', 'By  ', 'Bz  ', 'T1  ', 'rho ', 'T   ', 'Mach'], t=t, title='cztr3d')
    else
      call arbor_plot_tecplot(dom, name, P_GVB, iVarMin=o+1, iVarMax=o+n, names=['rho1', 'vx  ', 'vy  ', 'vz  ', &
          's1  ', 'T1  ', 'rho ', 'T   ', 'Mach'], t=t, title='cztr3d')
    end if
  end subroutine cz_write_plot

end module cz3d_app_mod
