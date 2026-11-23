module cztr3d_atm_mod
  ! The atmosphere physics of CZTR (cztr_physics_mod's atm_step and
  ! cztr_advance_mod's floors), ported from CZTR_2D's cztr2d_atm_mod to the
  ! 3D Yin-Yang shell with radial blocks and multirate (use_atmosphere, see
  ! cz3d_params_mod):
  !   - hyperbolic Spitzer conduction (HYDRAD's H I term above, Rempel
  !     saturation 1/(1 + |F|/(1.5 rho c_s^3))) at the faces in r, theta and
  !     phi, times the switch S at the face; each face flux relaxes to its
  !     target in tau = max(D/budget^2, 4 dt), budget = f_cfl_cond dx/dt - |v|,
  !     as CZTR. The relaxing fluxes are kept per block (qr_B, qt_B, qp_B); a
  !     face shared by two blocks has a copy in each, computed from the same
  !     temperatures through the ghost cells (CZTR_2D's rule for theta faces,
  !     now in every direction; with multirate the two copies relax with their
  !     own steps). No flux through the bottom and top walls. The heat goes in
  !     as CZTR's add_heat, ds = cv ln(1 + dQ/(rho cv T)), dQ = -dt div q
  !     (ARBOR's flux divergence: the exact areas);
  !   - Carlsson & Leenaarts + coronal losses (blended at 2.4e4 K): the column
  !     mass above a cell from p/g (hydrostatic estimate, no column sums
  !     across blocks), the neutral-H column integrated from the block's top
  !     down plus x_n n_H H_p of the cell above the block (its ghost); times S;
  !   - the heating S H_peak x e^(1-x), x = max(z, 0)/z_heat;
  !   - floors: rho >= rho_floor everywhere, T >= T_floor where S > 0.5.
  ! Losses and heating use the state before the conduction (as CZTR). The
  ! routines work on the bound block (cz_bind: the radial arrays are views
  ! of its cells, indexed 1..nb(1) with ghosts).
  use arbor
  use cz3d_params_mod
  use cz3d_bg_mod, only: n_r, k_join, r_c, r_f, z_c, rho0, s0, p0, g_c, tab, ddr_c, r_top_cz
  use radiative_loss_mod
  use eos_table, only: eos_fields_n, eos_field_n, eos_eval, eos_s_from_rho_T, n_field, F_LNP, F_LNT, F_GAMMA1, F_LNCV, F_LNNE
  implicit none
  private

  real(dp), parameter :: kappa0_spitzer = 7.8e-7_dp, T_SWITCH = 2.4e4_dp, T_SWITCH_WIDTH = 3.0e3_dp
  real(dp), parameter :: Rsun_atm = 6.957e10_dp, kB_atm = 1.380649e-16_dp
  integer, parameter :: RHO1 = 1, VR = 2, VT = 3, VP = 4, S1 = 5, BR = 6, BT = 7, BP = 8

  type(species_table_t) :: h1_tab, ca2_tab, mg2_tab
  type(coronal_table_t) :: coronal_tab
  type(hconduct_table_t) :: hcond_tab

  ! the background-force switch (cells 1:n_r of the whole column: 0 up to
  ! k_join - 2, 1 above) and the background force dp0/dr + rho0 g
  real(dp), allocatable, public :: S_bgf(:), bgforce(:)
  ! relaxing conductive fluxes (faces in r, theta, phi) and the last losses per block slot
  real(dp), allocatable :: qr_B(:,:,:,:), qt_B(:,:,:,:), qp_B(:,:,:,:)
  real(dp), allocatable, public :: Qloss_B(:,:,:,:), Qheat_B(:,:,:,:)   ! the last step's losses and heating [erg/cm3/s] (box output)
  ! the last step's s1 change by the conduction (box_fine_from: allocated by the caller)
  real(dp), allocatable, public :: ds1c_B(:,:,:,:)
  ! cond_limiter = 'mc' in use (with cond_along_b and the field): the states and
  ! the fill include the edge ghost cells
  logical, public :: cond_mc_on = .false.
  ! cond_hot_to_cold: faces with a conductive flux, and those set to 0 (counts,
  ! and the sums of |q| [erg/cm2/s]); summed over the run (floors report)
  real(dp), public :: hc_stat(4) = 0.0_dp
  ! the floors' action (atm_floors), summed over calls, per global radial cell:
  ! (1, i) cells checked, (2, i) rho floored, (3, i) T floored (timers_report prints them)
  ! (4, i) the energy the T floor added [erg] (rho T ds dV; Yin-Yang overlap cells count twice)
  real(dp), allocatable, public :: floor_cnt(:,:)
  ! one-block scratch: cells 0:nb(1)+1, 0:nb(2)+1, 0:nb(3)+1 (the faces' ghost cross only)
  real(dp), allocatable :: rho_A(:,:,:), T_A(:,:,:), cv_A(:,:,:), cs_A(:,:,:), ne_A(:,:,:), p_A(:,:,:)
  real(dp), allocatable :: F1(:,:,:,:), F2(:,:,:,:), F3(:,:,:,:), Dv(:,:,:,:)

  integer :: tm_st = 0, tm_co = 0, tm_lo = 0     ! (ARBOR timers of atm_block's parts)
  real(dp), parameter :: exp_m025 = exp(-0.25_dp) ! (face_flux: the relaxation factor at tau = 4 dt)
  ! the interior columns' states (with radial ghosts) per block, from states_interior:
  ! (0:nb(1)+1, j, k, 1 rho 2 T 3 cv 4 ne 5 p 6 cs, block)
  real(dp), allocatable :: stI_B(:,:,:,:,:)
  public :: atm_init, atm_setup, atm_block, atm_losses, states_interior, atm_active, atm_floors, atm_dt, atm_tr_heights, atm_switch, trac_f

contains

  ! The switch S at height z [cm] above the CZ top
  elemental real(dp) function atm_switch(z) result(s)
    real(dp), intent(in) :: z
    s = 0.5_dp*(1.0_dp + tanh((z - z_atm_km*1.0e5_dp)/(w_atm_km*1.0e5_dp)))
  end function atm_switch

  ! Once, after bg_load (whole column): tables, background force
  subroutine atm_init()
    real(dp) :: aH, aCa, aMg
    call load_abundance_mixture(trim(atomic_dir)//'/abundances/asplund.ab', trim(atomic_dir)//'/masses/masses.amu', &
        aH, aCa, aMg)
    call load_species_table('h_1', trim(atomic_dir)//'/OpticallyThick', aH, h1_tab)
    call load_species_table('ca_2', trim(atomic_dir)//'/OpticallyThick', aCa, ca2_tab)
    call load_species_table('mg_2', trim(atomic_dir)//'/OpticallyThick', aMg, mg2_tab)
    call load_coronal_table(coronal_file, coronal_tab)
    call load_hconduct_table(trim(atomic_dir)//'/OpticallyThick/thermal_conductivities/h_1.tc', hcond_tab)
    allocate(S_bgf(1:n_r), bgforce(1:n_r))
    ! left out (well balanced) only where the radial stencil (+-2 cells) sees
    ! the CZ background alone (CZTR_2D: in the last two CZ cells the left-out
    ! force was +1.0% and +2.6% of rho0 g, the release kick of 2026-10-30)
    S_bgf = 1.0_dp
    S_bgf(1:min(k_join - 2, n_r)) = 0.0_dp
    bgforce = ddr_c(p0) + rho0(1:n_r)*g_c(1:n_r)
  end subroutine atm_init

  ! After the blocks are (re)built: the per-block arrays (fluxes from 0) and scratch
  subroutine atm_setup()
    if (allocated(qr_B)) deallocate(qr_B, qt_B, qp_B, Qloss_B, Qheat_B)
    allocate(qr_B(0:nb(1), nb(2), nb(3), nBlockMax), qt_B(nb(1), 0:nb(2), nb(3), nBlockMax), &
        qp_B(nb(1), nb(2), 0:nb(3), nBlockMax), Qloss_B(nb(1), nb(2), nb(3), nBlockMax), Qheat_B(nb(1), nb(2), nb(3), nBlockMax))
    qr_B = 0.0_dp;  qt_B = 0.0_dp;  qp_B = 0.0_dp;  Qloss_B = 0.0_dp;  Qheat_B = 0.0_dp
    cond_mc_on = cond_limiter == 'mc' .and. cond_along_b .and. use_B
    ! the ARBOR timers on every rank (arbor_timer_report is collective: a rank
    ! without atmosphere blocks must hold the same timers)
    if (tm_st == 0) then
      tm_st = arbor_timer_id('atm: states');  tm_co = arbor_timer_id('atm: conduction')
      tm_lo = arbor_timer_id('atm: losses, heating')
    end if
    if (allocated(rho_A)) return
    allocate(rho_A(0:nb(1)+1, 0:nb(2)+1, 0:nb(3)+1))
    allocate(T_A, cv_A, cs_A, ne_A, p_A, mold=rho_A)
    rho_A = 1.0_dp;  T_A = 1.0_dp;  cv_A = 1.0_dp;  cs_A = 1.0_dp;  ne_A = 0.0_dp;  p_A = 1.0_dp
    allocate(F1(0:nb(1), nb(2), nb(3), 1), F2(nb(1), 0:nb(2), nb(3), 1), F3(nb(1), nb(2), 0:nb(3), 1), &
        Dv(nb(1), nb(2), nb(3), 1))
  end subroutine atm_setup

  ! The state of one radial run of cells i0..i1 at (j, k) of block ib
  subroutine column_states(ib, i0, i1, j, k)
    integer, intent(in) :: ib, i0, i1, j, k
    real(dp) :: rp(i0:i1), s(i0:i1), q(i0:i1, 5)
    integer :: n
    n = i1 - i0 + 1
    rp = rho0(i0:i1) + State_GVB(i0:i1, j, k, RHO1, ib)
    s = s0(i0:i1) + State_GVB(i0:i1, j, k, S1, ib)
    call eos_fields_n(tab, [F_LNT, F_LNCV, F_LNNE, F_LNP, F_GAMMA1], n, rp, s, q)
    rho_A(i0:i1, j, k) = rp
    T_A(i0:i1, j, k) = exp(q(:, 1))
    cv_A(i0:i1, j, k) = exp(q(:, 2))
    ne_A(i0:i1, j, k) = exp(q(:, 3))
    p_A(i0:i1, j, k) = exp(q(:, 4))
    cs_A(i0:i1, j, k) = sqrt(q(:, 5)*p_A(i0:i1, j, k)/rp)
  end subroutine column_states

  ! The cell states of block ib: the interior with its radial ghosts
  ! (0, nb(1)+1), and the theta and phi ghost layers next to the interior
  ! (the corners are not used and may be stale)
  subroutine states(ib)
    integer, intent(in) :: ib
    integer :: j, k
    do k = 1, nb(3)
      do j = 1, nb(2)
        call column_states(ib, 0, nb(1) + 1, j, k)
      end do
    end do
    call ghost_states(ib)
  end subroutine states

  ! The horizontal ghost columns' states (cells 1..nb(1)); with cond_mc_on
  ! also their radial ghosts and the theta-phi edge columns (the limited
  ! gradients' diagonal neighbours)
  subroutine ghost_states(ib)
    integer, intent(in) :: ib
    integer :: j, k, i0, i1
    i0 = merge(0, 1, cond_mc_on);  i1 = merge(nb(1) + 1, nb(1), cond_mc_on)
    do k = 1, nb(3)
      call column_states(ib, i0, i1, 0, k)
      call column_states(ib, i0, i1, nb(2) + 1, k)
    end do
    do j = 1, nb(2)
      call column_states(ib, i0, i1, j, 0)
      call column_states(ib, i0, i1, j, nb(3) + 1)
    end do
    if (cond_mc_on) then
      do k = 0, nb(3) + 1, nb(3) + 1
        do j = 0, nb(2) + 1, nb(2) + 1
          call column_states(ib, 1, nb(1), j, k)
        end do
      end do
    end if
  end subroutine ghost_states

  ! Conduction, losses and heating of block ib over dt (ghosts filled;
  ! at_bot/at_top: the block holds the bottom/top wall of the domain);
  ! stf(0:nb(1)+1): the step factors of the cells (local time stepping; with
  ! the cells below and above the block): cell i advances stf(i) dt, a radial
  ! face the smaller step of its two cells
  subroutine atm_block(ib, dt_in, at_bot, at_top, stf, losses_ready)
    integer, intent(in) :: ib
    real(dp), intent(in) :: dt_in
    logical, intent(in) :: at_bot, at_top
    real(dp), intent(in) :: stf(0:)
    logical, intent(in), optional :: losses_ready   ! Qloss_B, Qheat_B already this state's (atm_losses), states stored
    logical :: ready
    ! per block: the switch at the radial faces (sfr) and the cells (scc), the
    ! unresolved field's dB_r^2 there (dbf, dbc), sin theta (sth_j), and each
    ! interior cell's T gradient (gc), once instead of per face
    real(dp) :: sfr(0:nb(1)), scc(1:nb(1)), dbf(0:nb(1)), dbc(1:nb(1)), sth_j(1:nb(2)), &
        gc(3, nb(1), nb(2), nb(3))
    ! cond_limiter = 'mc': the one-sided T gradients between neighbouring cells,
    ! once per block (OsR: from i to i+1, OsT: from j to j+1, OsP: from k to k+1)
    real(dp) :: OsR(0:nb(1), 0:nb(2)+1, 0:nb(3)+1), OsT(0:nb(1)+1, 0:nb(2), 0:nb(3)+1), OsP(0:nb(1)+1, 0:nb(2)+1, 0:nb(3))
    real(dp) :: Tf, rf, Tg, csf, cvf, dxf, vf, x, dQ, Sc, Sf, sth, dt, fS, qcl, gt(3), bf(3), aniso
    real(dp) :: colH(nb(1)), dr, qd, ionH, lT, w, nHt, QH, iH, QCa, iCa, QMg, iMg, chromo, coron, Qheat, Qlos, colm
    integer :: i, j, k, m
    m = nb(1)
    if (maxval(atm_switch(z_c(1:m))) < 1.0e-8_dp .and. atm_switch(r_f(m) - r_top_cz) < 1.0e-8_dp) return  ! (deep CZ block)
    call arbor_timer_start(tm_st)
    ready = .false.
    if (present(losses_ready)) ready = losses_ready
    if (ready) then
      call states_from_store(ib, at_bot, at_top)
    else
      call states(ib)
    end if
    call arbor_timer_stop(tm_st)
    call arbor_timer_start(tm_co)
    aniso = merge(1.0_dp, 0.0_dp, cond_along_b .and. use_B)
    associate (X2 => Metric_DB(2, ib)%X_I, X3 => Metric_DB(3, ib)%X_I)
    if (use_conduction) then
      F1 = 0.0_dp;  F2 = 0.0_dp;  F3 = 0.0_dp
      do i = 0, m
        sfr(i) = atm_switch(r_f(i) - r_top_cz)
        dbf(i) = db2_at(r_f(i) - r_top_cz)
      end do
      do i = 1, m
        scc(i) = atm_switch(z_c(i))
        dbc(i) = db2_at(z_c(i))
      end do
      do j = 1, nb(2)
        sth_j(j) = sin(X2(j))
      end do
      if (aniso > 0.0_dp .and. cond_mc_on) then
        do k = 0, nb(3) + 1
          do j = 0, nb(2) + 1
            do i = 0, m
              OsR(i, j, k) = (T_A(i+1, j, k) - T_A(i, j, k))/(r_c(i+1) - r_c(i))
            end do
            if (j <= nb(2)) then
              do i = 0, m + 1
                OsT(i, j, k) = (T_A(i, j+1, k) - T_A(i, j, k))/(r_c(i)*(X2(j+1) - X2(j)))
              end do
            end if
            if (k <= nb(3)) then
              do i = 0, m + 1
                OsP(i, j, k) = (T_A(i, j, k+1) - T_A(i, j, k))/(r_c(i)*sin(X2(j))*(X3(k+1) - X3(k)))
              end do
            end if
          end do
        end do
      end if
      if (aniso > 0.0_dp .and. .not. cond_mc_on) then
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, m
              gc(:, i, j, k) = grad_cell(i, j, k)
            end do
          end do
        end do
      end if
      do k = 1, nb(3)
        do j = 1, nb(2)
          ! radial faces 0..m (none through the domain's bottom and top walls)
          do i = 0, m
            Sf = sfr(i)
            if (Sf < 1.0e-8_dp .or. (i == 0 .and. at_bot) .or. (i == m .and. at_top)) then
              qr_B(i, j, k, ib) = 0.0_dp
              ! an outflow top: the top cell passes the conductive energy on
              ! (the same F A through its top face as through the face below)
              if (i == m .and. at_top .and. wall_top == 'outflow') then
                qr_B(m, j, k, ib) = qr_B(m-1, j, k, ib)*(r_f(m-1)/r_f(m))**2
                F1(m, j, k, 1) = qr_B(m, j, k, ib)
              end if
              cycle
            end if
            dt = dt_in*min(stf(i), stf(i+1))
            dxf = r_c(i+1) - r_c(i)
            Tf = 0.5_dp*(T_A(i, j, k) + T_A(i+1, j, k))
            rf = 0.5_dp*(rho_A(i, j, k) + rho_A(i+1, j, k))
            csf = 0.5_dp*(cs_A(i, j, k) + cs_A(i+1, j, k))
            cvf = 0.5_dp*(cv_A(i, j, k) + cv_A(i+1, j, k))
            Tg = (T_A(i+1, j, k) - T_A(i, j, k))/dxf
            if (aniso > 0.0_dp .and. cond_mc_on) then
              ! along the field, the components along the face limited (cond_limiter)
              gt(1) = Tg
              gt(2) = mc(mc(OsT(i, j, k), OsT(i, j-1, k)), mc(OsT(i+1, j, k), OsT(i+1, j-1, k)))
              gt(3) = mc(mc(OsP(i, j, k), OsP(i, j, k-1)), mc(OsP(i+1, j, k), OsP(i+1, j, k-1)))
              bf = 0.5_dp*(State_GVB(i, j, k, BR:BP, ib) + State_GVB(i+1, j, k, BR:BP, ib))
              Tg = along_b(gt, bf, 1, dbf(i))
            else if (aniso > 0.0_dp) then
              ! along the field: the face's (r, theta, phi) gradient (tangential
              ! parts from its interior cells' central differences)
              gt = 0.0_dp;  x = 0.0_dp
              if (i >= 1) then;  gt = gt + gc(:, i, j, k);  x = x + 1.0_dp;  end if
              if (i + 1 <= m) then;  gt = gt + gc(:, i+1, j, k);  x = x + 1.0_dp;  end if
              gt = gt/max(x, 1.0_dp);  gt(1) = Tg
              bf = 0.5_dp*(State_GVB(i, j, k, BR:BP, ib) + State_GVB(i+1, j, k, BR:BP, ib))
              Tg = along_b(gt, bf, 1, dbf(i))
            end if
            vf = 0.5_dp*abs(State_GVB(i, j, k, VR, ib) + State_GVB(i+1, j, k, VR, ib))
            fS = 1.0_dp;  qcl = 0.0_dp
            if (cond_r_rsun > 0.0_dp) then
              ! AWSoM: Spitzer -> the collisionless flux (3/2) alpha p_e v_r far out
              fS = 1.0_dp/(1.0_dp + (r_f(i)/(cond_r_rsun*Rsun_atm))**2)
              qcl = 1.5_dp*alpha_cl*kB_atm*Tf*0.5_dp*(ne_A(i, j, k) + ne_A(i+1, j, k)) &
                  *0.5_dp*(State_GVB(i, j, k, VR, ib) + State_GVB(i+1, j, k, VR, ib))
            end if
            call face_flux(Sf, Tf, rf, csf, cvf, Tg, dxf, vf, qr_B(i, j, k, ib), fS, qcl)
            if (cond_hot_to_cold) call hot_to_cold(qr_B(i, j, k, ib), T_A(i+1, j, k) - T_A(i, j, k))
            F1(i, j, k, 1) = qr_B(i, j, k, ib)
          end do
        end do
        ! theta faces 0..nb(2) (the switch of the cell's height)
        do j = 0, nb(2)
          dxf = X2(j+1) - X2(j)
          do i = 1, m
            Sc = scc(i)
            if (Sc < 1.0e-8_dp) then
              qt_B(i, j, k, ib) = 0.0_dp
              cycle
            end if
            dt = dt_in*stf(i)
            Tf = 0.5_dp*(T_A(i, j, k) + T_A(i, j+1, k))
            rf = 0.5_dp*(rho_A(i, j, k) + rho_A(i, j+1, k))
            csf = 0.5_dp*(cs_A(i, j, k) + cs_A(i, j+1, k))
            cvf = 0.5_dp*(cv_A(i, j, k) + cv_A(i, j+1, k))
            Tg = (T_A(i, j+1, k) - T_A(i, j, k))/(r_c(i)*dxf)
            if (aniso > 0.0_dp .and. cond_mc_on) then
              gt(2) = Tg
              gt(1) = mc(mc(OsR(i, j, k), OsR(i-1, j, k)), mc(OsR(i, j+1, k), OsR(i-1, j+1, k)))
              gt(3) = mc(mc(OsP(i, j, k), OsP(i, j, k-1)), mc(OsP(i, j+1, k), OsP(i, j+1, k-1)))
              bf = 0.5_dp*(State_GVB(i, j, k, BR:BP, ib) + State_GVB(i, j+1, k, BR:BP, ib))
              Tg = along_b(gt, bf, 2, dbc(i))
            else if (aniso > 0.0_dp) then
              gt = 0.0_dp;  x = 0.0_dp
              if (j >= 1) then;  gt = gt + gc(:, i, j, k);  x = x + 1.0_dp;  end if
              if (j + 1 <= nb(2)) then;  gt = gt + gc(:, i, j+1, k);  x = x + 1.0_dp;  end if
              gt = gt/max(x, 1.0_dp);  gt(2) = Tg
              bf = 0.5_dp*(State_GVB(i, j, k, BR:BP, ib) + State_GVB(i, j+1, k, BR:BP, ib))
              Tg = along_b(gt, bf, 2, dbc(i))
            end if
            vf = 0.5_dp*abs(State_GVB(i, j, k, VT, ib) + State_GVB(i, j+1, k, VT, ib))
            fS = 1.0_dp
            if (cond_r_rsun > 0.0_dp) fS = 1.0_dp/(1.0_dp + (r_c(i)/(cond_r_rsun*Rsun_atm))**2)
            call face_flux(Sc, Tf, rf, csf, cvf, Tg, r_c(i)*dxf, vf, qt_B(i, j, k, ib), fS, 0.0_dp)
            if (cond_hot_to_cold) call hot_to_cold(qt_B(i, j, k, ib), T_A(i, j+1, k) - T_A(i, j, k))
            F2(i, j, k, 1) = qt_B(i, j, k, ib)
          end do
        end do
      end do
      ! phi faces 0..nb(3)
      do k = 0, nb(3)
        dxf = X3(k+1) - X3(k)
        do j = 1, nb(2)
          sth = sth_j(j)
          do i = 1, m
            Sc = scc(i)
            if (Sc < 1.0e-8_dp) then
              qp_B(i, j, k, ib) = 0.0_dp
              cycle
            end if
            dt = dt_in*stf(i)
            Tf = 0.5_dp*(T_A(i, j, k) + T_A(i, j, k+1))
            rf = 0.5_dp*(rho_A(i, j, k) + rho_A(i, j, k+1))
            csf = 0.5_dp*(cs_A(i, j, k) + cs_A(i, j, k+1))
            cvf = 0.5_dp*(cv_A(i, j, k) + cv_A(i, j, k+1))
            Tg = (T_A(i, j, k+1) - T_A(i, j, k))/(r_c(i)*sth*dxf)
            if (aniso > 0.0_dp .and. cond_mc_on) then
              gt(3) = Tg
              gt(1) = mc(mc(OsR(i, j, k), OsR(i-1, j, k)), mc(OsR(i, j, k+1), OsR(i-1, j, k+1)))
              gt(2) = mc(mc(OsT(i, j, k), OsT(i, j-1, k)), mc(OsT(i, j, k+1), OsT(i, j-1, k+1)))
              bf = 0.5_dp*(State_GVB(i, j, k, BR:BP, ib) + State_GVB(i, j, k+1, BR:BP, ib))
              Tg = along_b(gt, bf, 3, dbc(i))
            else if (aniso > 0.0_dp) then
              gt = 0.0_dp;  x = 0.0_dp
              if (k >= 1) then;  gt = gt + gc(:, i, j, k);  x = x + 1.0_dp;  end if
              if (k + 1 <= nb(3)) then;  gt = gt + gc(:, i, j, k+1);  x = x + 1.0_dp;  end if
              gt = gt/max(x, 1.0_dp);  gt(3) = Tg
              bf = 0.5_dp*(State_GVB(i, j, k, BR:BP, ib) + State_GVB(i, j, k+1, BR:BP, ib))
              Tg = along_b(gt, bf, 3, dbc(i))
            end if
            vf = 0.5_dp*abs(State_GVB(i, j, k, VP, ib) + State_GVB(i, j, k+1, VP, ib))
            fS = 1.0_dp
            if (cond_r_rsun > 0.0_dp) fS = 1.0_dp/(1.0_dp + (r_c(i)/(cond_r_rsun*Rsun_atm))**2)
            call face_flux(Sc, Tf, rf, csf, cvf, Tg, r_c(i)*sth*dxf, vf, qp_B(i, j, k, ib), fS, 0.0_dp)
            if (cond_hot_to_cold) call hot_to_cold(qp_B(i, j, k, ib), T_A(i, j, k+1) - T_A(i, j, k))
            F3(i, j, k, 1) = qp_B(i, j, k, ib)
          end do
        end do
      end do
      call arbor_flux_divergence(ib, F1, F2, F3, Dv)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, m
            dQ = -dt_in*stf(i)*Dv(i, j, k, 1)
            if (dQ == 0.0_dp) cycle
            x = dQ/(rho_A(i, j, k)*cv_A(i, j, k)*T_A(i, j, k))
            State_GVB(i, j, k, S1, ib) = State_GVB(i, j, k, S1, ib) + cv_A(i, j, k)*log(1.0_dp + max(x, -0.9_dp))
            if (allocated(ds1c_B)) ds1c_B(i, j, k, ib) = cv_A(i, j, k)*log(1.0_dp + max(x, -0.9_dp))
          end do
        end do
      end do
    end if

    call arbor_timer_stop(tm_co)
    call arbor_timer_start(tm_lo)
    ! losses and heating (the state before the conduction: computed by
    ! atm_losses, during the ghost fill when losses_ready), applied now
    if (.not. ready) call atm_losses(ib, at_top)
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, m
          State_GVB(i, j, k, S1, ib) = State_GVB(i, j, k, S1, ib) &
              + dt_in*stf(i)*(Qheat_B(i, j, k, ib) - Qloss_B(i, j, k, ib))/(rho_A(i, j, k)*T_A(i, j, k))
        end do
      end do
    end do
    call arbor_timer_stop(tm_lo)
    end associate

  contains

    ! cond_hot_to_cold: the face flux q (positive towards the cell with the
    ! larger index) set to 0 if it carries heat from the colder cell to the
    ! hotter one (dT = T(next) - T(this)); counted in hc_stat
    subroutine hot_to_cold(q, dT)
      real(dp), intent(inout) :: q
      real(dp), intent(in) :: dT
      if (q == 0.0_dp) return
      hc_stat(1) = hc_stat(1) + 1.0_dp;  hc_stat(3) = hc_stat(3) + abs(q)
      if (q*dT > 0.0_dp) then
        hc_stat(2) = hc_stat(2) + 1.0_dp;  hc_stat(4) = hc_stat(4) + abs(q)
        q = 0.0_dp
      end if
    end subroutine hot_to_cold

    ! The (r, theta, phi) gradient of T at interior cell (i, j, k) of the block
    ! (central differences; the ghost layers next to the interior are filled)
    function grad_cell(i, j, k) result(g)
      integer, intent(in) :: i, j, k
      real(dp) :: g(3)
      associate (Y2 => Metric_DB(2, ib)%X_I, Y3 => Metric_DB(3, ib)%X_I)
        g(1) = (T_A(i+1, j, k) - T_A(i-1, j, k))/(r_c(i+1) - r_c(i-1))
        g(2) = (T_A(i, j+1, k) - T_A(i, j-1, k))/(r_c(i)*(Y2(j+1) - Y2(j-1)))
        g(3) = (T_A(i, j, k+1) - T_A(i, j, k-1))/(r_c(i)*sin(Y2(j))*(Y3(k+1) - Y3(k-1)))
      end associate
    end function grad_cell

    ! The gradient that drives the face flux along direction d at height zf
    ! [cm] (cond_along_b): (<b b> g)_d with the unresolved vertical field of
    ! rms dB_r(z) and the isotropic part b_iso_G,
    !   <b b> = (B B + dB_r^2 e_r e_r + b_iso^2 I)/(B^2 + dB_r^2 + b_iso^2)
    ! dB_r(z) = cond_dB0_G (1 - tanh((z - cond_dB_zj_km)/cond_dB_w_km))/2
    ! (0 when cond_dB0_G = 0: then w b_d (b.g) + (1 - w) g_d, w = B^2/(B^2 + b_iso^2))
    ! (db2: dB_r^2 at the face's height, db2_at)
    real(dp) function along_b(g, bv, d, db2) result(gd)
      real(dp), intent(in) :: g(3), bv(3), db2
      integer, intent(in) :: d
      real(dp) :: dn
      dn = sum(bv**2) + db2 + b_iso_G**2
      gd = (bv(d)*dot_product(bv, g) + b_iso_G**2*g(d))/dn
      if (d == 1) gd = gd + db2*g(1)/dn
    end function along_b

    ! dB_r(z)^2 of along_b at height zf [cm]
    real(dp) function db2_at(zf) result(db2)
      real(dp), intent(in) :: zf
      db2 = 0.0_dp
      if (cond_dB0_G > 0.0_dp) &
          db2 = (cond_dB0_G*0.5_dp*(1.0_dp - tanh((zf - cond_dB_zj_km*1.0e5_dp)/(cond_dB_w_km*1.0e5_dp))))**2
    end function db2_at

    ! CZTR's hyperbolic face flux: target -fsat k0 T^2.5 dT/dx, relaxed;
    ! fS: the Spitzer weight, qcl: the collisionless flux (weight 1 - fS)
    subroutine face_flux(Sf, Tf, rf, csf, cvf, Tg, dx, vf, q, fS, qcl)
      real(dp), intent(in) :: Sf, Tf, rf, csf, cvf, Tg, dx, vf, fS, qcl
      real(dp), intent(inout) :: q
      real(dp) :: k0, Fc, qt, D, budget, tau, t25, wb
      ! (the neutral-H term only where its blend weight is nonzero: 0 x it is 0;
      ! the weight 0.5 (1 - tanh(x)) is exactly 0 for x > 19.1 in double precision)
      if ((Tf - T_SWITCH)/T_SWITCH_WIDTH > 20.0_dp) then
        wb = 0.0_dp
      else
        wb = blend_weight(Tf, T_SWITCH, T_SWITCH_WIDTH)
      end if
      if (wb > 0.0_dp) then
        k0 = Sf*(kappa0_spitzer/trac_f(Tf) + wb*extra_h_kappa0(hcond_tab, log10(Tf)))
      else
        k0 = Sf*(kappa0_spitzer/trac_f(Tf) + 0.0_dp)
      end if
      t25 = Tf**2.5_dp
      Fc = k0*t25*Tg
      qt = -Fc/(1.0_dp + abs(Fc)/(1.5_dp*rf*csf**3))
      qt = fS*qt + (1.0_dp - fS)*Sf*qcl
      D = fS*k0*t25/(rf*cvf)
      ! (cond_tau_fsat: the saturation factor also in tau, Rempel 2017 eq. 16 with
      ! f_sat: lower tau where the flux saturates, the flux follows its target)
      if (cond_tau_fsat) D = D/(1.0_dp + abs(Fc)/(1.5_dp*rf*csf**3))
      budget = max(f_cfl_cond*dx/dt - vf, 1.0e-3_dp*dx/dt)
      tau = max(D/budget**2, 4.0_dp*dt)
      ! (at tau's floor 4 dt, dt/tau is exactly 0.25)
      if (tau == 4.0_dp*dt) then
        q = qt + (q - qt)*exp_m025
      else
        q = qt + (q - qt)*exp(-dt/tau)
      end if
    end subroutine face_flux
  end subroutine atm_block

  ! The losses Qloss_B and the heating Qheat_B [erg cm^-3 s^-1] of block ib's
  ! interior cells from the block's states (rho_A, T_A, p_A, ne_A of the
  ! interior and its radial ghosts: states_interior or states), the state
  ! before the conduction. at_top: the block holds the domain's top wall
  subroutine atm_losses(ib, at_top)
    integer, intent(in) :: ib
    logical, intent(in) :: at_top
    real(dp) :: colH(nb(1)), dr, qd, ionH, lT, w, nHt, QH, iH, QCa, iCa, QMg, iMg, chromo, coron, Qheat, Qlos, colm, &
        Sc, x
    integer :: i, j, k, m
    m = nb(1)
    ! losses and heating (the state before the conduction)
    do k = 1, nb(3)
      do j = 1, nb(2)
        if (use_losses) then
          ! the neutral-H column: from above the block (the ghost cell's
          ! x_n n_H H_p; 0 at the top wall), then down through the block
          colH = 0.0_dp
          if (.not. at_top) then
            call species_emiss_ionfrac(h1_tab, log10(T_A(m+1, j, k)), qd, ionH)
            colH(m) = ionH*N_H_per_gram_global*rho_A(m+1, j, k)*p_A(m+1, j, k)/(rho_A(m+1, j, k)*g_c(m+1))
          end if
          do i = m, 1, -1
            dr = r_f(i) - r_f(i-1)
            call species_emiss_ionfrac(h1_tab, log10(T_A(i, j, k)), qd, ionH)
            if (i < m) colH(i) = colH(i+1)
            colH(i) = colH(i) + ionH*N_H_per_gram_global*rho_A(i, j, k)*dr
          end do
        end if
        do i = 1, m
          Qlos = 0.0_dp
          Qheat = 0.0_dp
          Sc = atm_switch(z_c(i))
          if (Sc >= 1.0e-8_dp) then
            if (use_losses) then
              lT = log10(T_A(i, j, k))
              w = blend_weight(T_A(i, j, k), T_SWITCH, T_SWITCH_WIDTH)
              nHt = N_H_per_gram_global*rho_A(i, j, k)
              chromo = 0.0_dp
              if (w > 1.0e-6_dp) then
                colm = p_A(i, j, k)/g_c(i)                ! column mass above (hydrostatic)
                call species_emiss_ionfrac(h1_tab, lT, QH, iH)
                call species_emiss_ionfrac(ca2_tab, lT, QCa, iCa)
                call species_emiss_ionfrac(mg2_tab, lT, QMg, iMg)
                chromo = (QH*species_escape_prob(h1_tab, log10(max(4.0e-14_dp*colH(i), 1.0e-30_dp))) &
                    + QCa*species_escape_prob(ca2_tab, log10(max(colm, 1.0e-30_dp))) &
                    + QMg*species_escape_prob(mg2_tab, log10(max(colm, 1.0e-30_dp))))*nHt*ne_A(i, j, k)
              end if
              coron = 0.0_dp
              if (w < 1.0_dp - 1.0e-6_dp) coron = ne_A(i, j, k)*nHt*coronal_lambda(coronal_tab, lT)
              ! TR broadening: the whole loss times f (1 below trac_Tlow)
              Qlos = Sc*(w*chromo + (1.0_dp - w)*coron)*trac_f(T_A(i, j, k))
            end if
            if (use_heating) then
              x = max(z_c(i), 0.0_dp)/(z_heat_km*1.0e5_dp)
              Qheat = Sc*H_peak*x*exp(1.0_dp - x)
              if (H_c > 0.0_dp) Qheat = Qheat + Sc*H_c*exp(-max(z_c(i), 0.0_dp)/(lambda_c_km*1.0e5_dp))
              if (H_ext > 0.0_dp) Qheat = Qheat + Sc*H_ext*exp(-max(z_c(i), 0.0_dp)/(lambda_ext_rsun*Rsun_atm))
              if (H_B > 0.0_dp .and. use_B) Qheat = Qheat + Sc*H_B*bsat(norm2(State_GVB(i, j, k, BR:BP, ib))) &
                  *exp(-max(z_c(i), 0.0_dp)/(lambda_B_km*1.0e5_dp))/(1.0_dp + (ne_A(i, j, k)/n_gate)**4) &
                  *merge(min(1.0_dp, rho_A(i, j, k)/rho0(i))**H_B_rho_exp, 1.0_dp, H_B_rho_exp > 0.0_dp)
            end if
          end if
          Qloss_B(i, j, k, ib) = Qlos
          Qheat_B(i, j, k, ib) = Qheat
        end do
      end do
    end do
  end subroutine atm_losses

  ! The states of block ib's interior columns with their radial ghosts (what
  ! atm_losses reads; the radial neighbours are on this rank with si_radial)
  ! They are kept per block (stI_B) for atm_block (losses_ready): it then
  ! evaluates only the horizontal ghost columns
  subroutine states_interior(ib)
    integer, intent(in) :: ib
    integer :: j, k, m
    m = nb(1)
    if (.not. allocated(stI_B)) allocate(stI_B(0:m+1, nb(2), nb(3), 6, size(State_GVB, 5)))
    if (size(stI_B, 5) /= size(State_GVB, 5)) then
      deallocate(stI_B)
      allocate(stI_B(0:m+1, nb(2), nb(3), 6, size(State_GVB, 5)))
    end if
    do k = 1, nb(3)
      do j = 1, nb(2)
        call column_states(ib, 0, m + 1, j, k)
        stI_B(:, j, k, 1, ib) = rho_A(0:m+1, j, k);  stI_B(:, j, k, 2, ib) = T_A(0:m+1, j, k)
        stI_B(:, j, k, 3, ib) = cv_A(0:m+1, j, k);  stI_B(:, j, k, 4, ib) = ne_A(0:m+1, j, k)
        stI_B(:, j, k, 5, ib) = p_A(0:m+1, j, k);  stI_B(:, j, k, 6, ib) = cs_A(0:m+1, j, k)
      end do
    end do
  end subroutine states_interior

  ! states with the interior columns from states_interior's store (the same
  ! values), the horizontal ghost columns evaluated, and the radial ghost
  ! at a physical wall (at_bot/at_top) too: ARBOR applies the boundary of a
  ! block that waits for messages only at the fill's end
  subroutine states_from_store(ib, at_bot, at_top)
    integer, intent(in) :: ib
    logical, intent(in) :: at_bot, at_top
    integer :: j, k, m
    m = nb(1)
    do k = 1, nb(3)
      do j = 1, nb(2)
        rho_A(0:m+1, j, k) = stI_B(:, j, k, 1, ib);  T_A(0:m+1, j, k) = stI_B(:, j, k, 2, ib)
        cv_A(0:m+1, j, k) = stI_B(:, j, k, 3, ib);  ne_A(0:m+1, j, k) = stI_B(:, j, k, 4, ib)
        p_A(0:m+1, j, k) = stI_B(:, j, k, 5, ib);  cs_A(0:m+1, j, k) = stI_B(:, j, k, 6, ib)
        if (at_bot) call column_states(ib, 0, 0, j, k)
        if (at_top) call column_states(ib, m + 1, m + 1, j, k)
      end do
    end do
    call ghost_states(ib)
  end subroutine states_from_store

  ! Does block ib have atmosphere work (atm_block's own test)?
  logical function atm_active(ib)
    integer, intent(in) :: ib
    integer :: m
    m = nb(1)
    atm_active = .not. (maxval(atm_switch(z_c(1:m))) < 1.0e-8_dp .and. atm_switch(r_f(m) - r_top_cz) < 1.0e-8_dp)
  end function atm_active

  ! The monotonized-central limiter: minmod(2a, 2b, (a + b)/2)
  pure real(dp) function mc(a, b)
    real(dp), intent(in) :: a, b
    mc = merge(sign(min(2.0_dp*abs(a), 2.0_dp*abs(b), 0.5_dp*abs(a + b)), a), 0.0_dp, a*b > 0.0_dp)
  end function mc

  ! |B| for the field-strength heating, saturated at H_B_sat_G (> 0)
  real(dp) function bsat(b)
    real(dp), intent(in) :: b
    bsat = b
    if (H_B_sat_G > 0.0_dp) bsat = min(b, H_B_sat_G)
  end function bsat

  ! Floors of block ib (interior cells): rho >= rho_floor; T >= T_floor where S > 0.5.
  ! T by the EOS's column evaluator (ln T only; the same values as eos_eval)
  ! over the row's cells with S > 0.5; the inversion only where T < T_floor
  subroutine atm_floors(ib)
    integer, intent(in) :: ib
    real(dp) :: rh, rp(nb(1)), ss(nb(1)), lt(nb(1))
    integer :: i, j, k, n, idx(nb(1)), m, r0
    logical :: on(nb(1))
    if (.not. allocated(floor_cnt)) then
      allocate(floor_cnt(4, n_r))
      floor_cnt = 0.0_dp
    end if
    r0 = Key_IB(KI+1, ib)*nb(1)                      ! (the block's radial offset)
    m = 0
    do i = 1, nb(1)
      on(i) = atm_switch(z_c(i)) > 0.5_dp
      if (on(i)) then
        m = m + 1
        idx(m) = i
      end if
    end do
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, nb(1)
          if (State_GVB(i, j, k, RHO1, ib) < rho_floor - rho0(i)) floor_cnt(2, r0 + i) = floor_cnt(2, r0 + i) + 1.0_dp
          State_GVB(i, j, k, RHO1, ib) = max(State_GVB(i, j, k, RHO1, ib), rho_floor - rho0(i))
          floor_cnt(1, r0 + i) = floor_cnt(1, r0 + i) + 1.0_dp
        end do
        if (m == 0) cycle
        do n = 1, m
          i = idx(n)
          rp(n) = rho0(i) + State_GVB(i, j, k, RHO1, ib)
          ss(n) = s0(i) + State_GVB(i, j, k, S1, ib)
        end do
        call eos_field_n(tab, F_LNT, m, rp, ss, lt)
        do n = 1, m
          if (exp(lt(n)) < T_floor) then
            i = idx(n)
            floor_cnt(3, r0 + i) = floor_cnt(3, r0 + i) + 1.0_dp
            rh = rp(n)
            State_GVB(i, j, k, S1, ib) = eos_s_from_rho_T(tab, rh, T_floor, s_guess=ss(n)) - s0(i)
            floor_cnt(4, r0 + i) = floor_cnt(4, r0 + i) + rh*0.5_dp*(exp(lt(n)) + T_floor) &
                *(s0(i) + State_GVB(i, j, k, S1, ib) - ss(n)) &
                *Metric_DB(1, ib)%Vol_I(i)*Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          end if
        end do
      end do
    end do
  end subroutine atm_floors

  ! The TR broadening factor f(T) <= 1 (trac_Tc, trac_Tlow; 1 when off)
  pure real(dp) function trac_f(T) result(f)
    real(dp), intent(in) :: T
    real(dp) :: x, r
    f = 1.0_dp
    if (trac_Tc <= 0.0_dp .or. T >= trac_Tc .or. T <= trac_Tlow) return
    x = min(1.0_dp, log(T/trac_Tlow)/log(1.5_dp))
    r = x*x*(3.0_dp - 2.0_dp*x)
    f = exp(2.5_dp*r*log(T/trac_Tc))
  end function trac_f

  ! dt_rad_factor x the shortest loss time rho cv T/Q_loss of block ib (the
  ! last step's losses; huge if none); rate_r (optional): per radial cell,
  ! raised to cfl over the cell's own loss-time limit (local time stepping)
  real(dp) function atm_dt(ib, rate_r, i_min) result(dt)
    integer, intent(in) :: ib
    real(dp), intent(inout), optional :: rate_r(:)
    integer, intent(out), optional :: i_min        ! the block's cell (radial index) of the limit
    real(dp) :: rh, d, rpw(nb(1)), sw(nb(1)), q2(nb(1), 2)
    integer :: i, j, k, n, l, idx(nb(1))
    logical :: on(nb(1))
    if (present(i_min)) i_min = 0
    dt = huge(1.0_dp)
    if (.not. (use_losses .and. allocated(Qloss_B))) return
    ! ln cv, ln T by the column evaluator over each row's cells with S > 0.5
    ! and a loss (the same values as eos_eval), the height switch per index
    do i = 1, nb(1)
      on(i) = atm_switch(z_c(i)) > 0.5_dp
    end do
    do k = 1, nb(3)
      do j = 1, nb(2)
        n = 0
        do i = 1, nb(1)
          if (.not. on(i) .or. Qloss_B(i, j, k, ib) <= 0.0_dp) cycle
          n = n + 1
          idx(n) = i
          rpw(n) = rho0(i) + State_GVB(i, j, k, RHO1, ib)
          sw(n) = s0(i) + State_GVB(i, j, k, S1, ib)
        end do
        if (n == 0) cycle
        call eos_fields_n(tab, [F_LNCV, F_LNT], n, rpw(1:n), sw(1:n), q2(1:n, :))
        do l = 1, n
          i = idx(l)
          rh = rpw(l)
          d = rh*exp(q2(l, 1))*exp(q2(l, 2))/Qloss_B(i, j, k, ib)
          if (d < dt) then
            dt = d
            if (present(i_min)) i_min = i
          end if
          if (present(rate_r)) rate_r(i) = max(rate_r(i), cfl*Qloss_B(i, j, k, ib) &
              /(dt_rad_factor*rh*exp(q2(l, 1))*exp(q2(l, 2))))
        end do
      end do
    end do
    dt = dt_rad_factor*dt
  end function atm_dt

  ! The transition region of block ib's columns: the heights [km] of the
  ! upward crossings of T = 1e5 K with z between 0 and z_max_km (a cell above
  ! 1e5 K whose lower neighbour, possibly the ghost, is below); hmin/hmax are
  ! updated (start them at huge/-huge); the lowest crossing per column only
  subroutine atm_tr_heights(ib, z_max_km, hmin, hmax)
    integer, intent(in) :: ib
    real(dp), intent(in) :: z_max_km
    real(dp), intent(inout) :: hmin, hmax
    real(dp) :: out(n_field), Tlo, Thi
    integer :: i, j, k
    do k = 1, nb(3)
      do j = 1, nb(2)
        call eos_eval(tab, rho0(0) + State_GVB(0, j, k, RHO1, ib), s0(0) + State_GVB(0, j, k, S1, ib), out)
        Tlo = exp(out(F_LNT))
        do i = 1, nb(1)
          call eos_eval(tab, rho0(i) + State_GVB(i, j, k, RHO1, ib), s0(i) + State_GVB(i, j, k, S1, ib), out)
          Thi = exp(out(F_LNT))
          if (z_c(i) > 0.0_dp .and. z_c(i) < z_max_km*1.0e5_dp .and. Thi > 1.0e5_dp .and. Tlo <= 1.0e5_dp) then
            hmin = min(hmin, z_c(i)/1.0e5_dp)
            hmax = max(hmax, z_c(i)/1.0e5_dp)
            exit
          end if
          Tlo = Thi
        end do
      end do
    end do
  end subroutine atm_tr_heights

end module cztr3d_atm_mod
