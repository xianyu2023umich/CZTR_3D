module cz3d_transport_mod
  ! Radial energy transport of one column (cells 1..n_r, radial ghosts
  ! -1, 0, n_r+1, n_r+2), as CZTR's transport_step with step 1's settings:
  !   - mixing-length flux at the faces 0..n_r (CZ_1D's targets, relaxed
  !     hyperbolically), times the subgrid factor f_mlt of the face;
  !     divergence -> ds1 = dt (A F)_div/(V rho T);
  !   - grey Eddington radiation (grey_core, implicit coupling) over cells
  !     1..n_act (the last opaque cell, at most n_r), the Eddington top
  !     condition with optical depth grey_tau_above above cell n_act (CZTR's
  !     grey_top_at_join), heating Q per cell;
  !   - the bottom face carries L_bottom, split between radiation and
  !     convection as their targets.
  ! Each block works on its own part of the column (its radial ghosts filled
  ! first); a face between radial blocks is computed by both (the relaxing
  ! convective flux kept by each, not exchanged), the implicit grey solve
  ! only in the top block (grey_top_km), explicit diffusion below. Opacity: Rosseland table, zero outside its coverage.
  use kinds_mod, only: dp
  use const_mod, only: a_rad, c_light, pi
  use cz3d_bg_mod, only: n_r, n_c, r_c, r_f, area_f, vol_c, g_f, s0, ddr_f, face_value
  use opacity_table_mod, only: opac_table_t, opac_load, opac_get_logkappa
  use vmath_mod, only: vexp
  use mlt_mod, only: mlt_flux
  use grey_rt_mod, only: grey_core
  implicit none
  private

  type(opac_table_t) :: opac
  real(dp), parameter :: kr_floor = 1.0e-30_dp

  ! top_energy = 'cooling': the cooling rate per volume of each cell (> 0 removes heat)
  real(dp), pointer, public :: cool_c(:) => null()       ! cooling rate per volume of each cell (0 without cooling)
  real(dp), pointer, public :: cool_f(:) => null()  ! the cooling layer's flux at faces 0:n_r (0 without cooling)
  real(dp), pointer :: cool_c_g(:) => null(), cool_f_g(:) => null()
  ! the current column's walls: its block holds the bottom (top) of the domain
  logical, public :: col_bottom = .true., col_top = .true.
  ! first cell of the implicit grey solve (set from grey_top_km): below it the
  ! radiative flux is the explicit diffusion flux (the optically thick limit)
  integer, public :: grey_i0 = 1
  logical :: top_cooling = .false.

  public :: transport_init, set_top_cooling, transport_column, kappa, tr_view

contains

  subroutine transport_init(opacity_file)
    character(*), intent(in) :: opacity_file
    call opac_load(opac, opacity_file)
    allocate(cool_f(0:n_r), cool_c(1:n_r))
    cool_f = 0.0_dp
    cool_c = 0.0_dp
  end subroutine transport_init

  ! Cooling layer below the top face: flux F(r) = F_top exp(-(r_f(n_r) - r)/w),
  ! F_top = L/(4 pi r_top^2), so the layer removes L (to exp(-depth/w));
  ! per cell the conservative difference over its faces
  subroutine set_top_cooling(L, w)
    real(dp), intent(in) :: L, w
    real(dp) :: F(0:n_r)
    integer :: i
    top_cooling = .true.
    do i = 0, n_r
      F(i) = L/(4.0_dp*pi*area_f(n_r))*exp(-(r_f(n_r) - r_f(i))/w)
    end do
    cool_f = F
    do i = 1, n_r
      cool_c(i) = (area_f(i)*F(i) - area_f(i-1)*F(i-1))/vol_c(i)
    end do
  end subroutine set_top_cooling

  ! Rosseland opacity, zero outside the table's coverage
  real(dp) function kappa(rho, T)
    real(dp), intent(in) :: rho, T
    kappa = kappa_ln(log(rho), log(T))
  end function kappa

  ! The same from ln rho and ln T (the equation of state's own variables: no
  ! log of an exp): one table lookup, log10 by a product, kappa by exp
  real(dp) function kappa_ln(lnrho, lnT)
    real(dp), intent(in) :: lnrho, lnT
    real(dp), parameter :: ln10 = log(10.0_dp), rln10 = 1.0_dp/log(10.0_dp)
    real(dp) :: logT, logR, Tmax, Rmax, logk
    logical :: oot
    logT = lnT*rln10
    logR = lnrho*rln10 - 3.0_dp*(logT - 6.0_dp)
    Tmax = opac%Tmin + opac%dT*real(opac%nT - 1, dp)
    Rmax = opac%Rmin + opac%dR*real(opac%nR - 1, dp)
    kappa_ln = 0.0_dp
    if (logT < opac%Tmin .or. logT > Tmax .or. logR < opac%Rmin .or. logR > Rmax) return
    call opac_get_logkappa(opac, logT, logR, logk, oot)
    if (.not. oot) kappa_ln = exp(logk*ln10)
  end function kappa_ln

  ! One column. State arrays on cells -1:n_c+2 (ghosts filled); s1 updated
  ! in cells 1..n_c. q_conv(0:n_c) is the column's relaxing convective flux
  ! (kept between steps); first: set it to its targets (no relaxation).
  ! f_mlt(0:n_c): subgrid factor of the mixing-length flux. Out: the
  ! radiative flux q_rad(0:n_c) and heating Qrad(1:n_c).
  ! top_wall_face: the top face is a wall (no convective flux through it);
  ! else, as CZTR step 1, the top face lies inside the 1D grid.
  ! tsc(0:n_c+1) (optional, local time stepping): cell i advances tsc(i) dt,
  ! a face min of its two cells; tsc(0), tsc(n_c+1): the neighbouring cells
  ! below and above the column part (another radial block's, or the end
  ! cell's own at a wall); the grey coupling sees rho cv/tsc (the same
  ! implicit factor eta as with the cell's own step)
  ! lnrho, lnT, lnP, lncp: the logarithms of rho, T, P, cp (from the equation
  ! of state; the face values are interpolated in them)
  subroutine transport_column(dt, rho, T, P, cp, cv, Q, lnrho, lnT, lnP, lncp, v, s1, q_conv, f_mlt, alpha, cfl, &
      L_bottom, tau_above, use_rad, use_conv, use_hyp, first, q_rad, Qrad, top_wall_face, tsc)
    real(dp), intent(in) :: dt, alpha, cfl, L_bottom, tau_above
    real(dp), intent(in) :: rho(-1:), T(-1:), P(-1:), cp(-1:), cv(-1:), Q(-1:), v(-1:), f_mlt(0:)
    real(dp), intent(in) :: lnrho(-1:), lnT(-1:), lnP(-1:), lncp(-1:)
    real(dp), intent(inout) :: s1(-1:), q_conv(0:)
    logical, intent(in) :: use_rad, use_conv, use_hyp, first, top_wall_face
    real(dp), intent(out) :: q_rad(0:), Qrad(1:)
    real(dp), intent(in), optional :: tsc(0:)
    real(dp) :: sc(1:n_c), dtc(1:n_c), dtf(0:n_c), rho_f(0:n_c), T_f(0:n_c), P_f(0:n_c), cp_f(0:n_c), Q_f(0:n_c), dTdr(0:n_c), dsdr(0:n_c)
    real(dp) :: stot(-1:n_c+2), Fc_t(0:n_c), K_s(0:n_c), Fr_t0, kap, y, Hp, budget, one_over_tau, vf, F_bot
    real(dp) :: lnrho_f(0:n_c), lnT_f(0:n_c), lnP_f(0:n_c), lncp_f(0:n_c)
    integer :: i

    sc = 1.0_dp
    if (present(tsc)) sc = tsc(1:n_c)
    dtc = dt*sc
    dtf(0) = dtc(1)
    dtf(n_c) = dtc(n_c)
    if (present(tsc)) then
      dtf(0) = dt*min(tsc(0), tsc(1))
      dtf(n_c) = dt*min(tsc(n_c), tsc(n_c+1))
    end if
    do i = 1, n_c - 1
      dtf(i) = min(dtc(i), dtc(i+1))
    end do
    lnrho_f = face_value(lnrho)
    lnT_f = face_value(lnT)
    call vexp(n_c + 1, lnrho_f, rho_f)
    call vexp(n_c + 1, lnT_f, T_f)
    lnP_f = face_value(lnP)
    call vexp(n_c + 1, lnP_f, P_f)
    lncp_f = face_value(lncp)
    call vexp(n_c + 1, lncp_f, cp_f)
    Q_f = face_value(Q)
    dTdr = ddr_f(T)
    stot = s0 + s1
    dsdr = ddr_f(stot)

    ! convective flux targets
    Fc_t = 0.0_dp
    K_s = 0.0_dp
    if (use_conv) then
      do i = 0, n_c
        kap = kappa_ln(lnrho_f(i), lnT_f(i))
        if (kap <= 0.0_dp) cycle
        Hp = P_f(i)/(rho_f(i)*g_f(i))
        y = -(Hp/cp_f(i))*dsdr(i)
        call mlt_flux(y, rho_f(i), T_f(i), P_f(i), g_f(i), kap, cp_f(i), Q_f(i), alpha, Fc_t(i), K_s(i))
        Fc_t(i) = f_mlt(i)*Fc_t(i)
        K_s(i) = f_mlt(i)*K_s(i)
      end do
    end if
    if (first .or. .not. use_hyp) then
      q_conv = Fc_t
    else if (dt > 0.0_dp) then                  ! (dt = 0: the flux as it is, e.g. after a restart)
      do i = 0, n_c
        vf = max(abs(v(i)), abs(v(i+1)))
        budget = max(cfl*(r_c(i+1) - r_c(i))/dtf(i) - vf, 1.0e-3_dp*(r_c(i+1) - r_c(i))/dtf(i))
        one_over_tau = 1.0_dp/(4.0_dp*dtf(i))
        if (K_s(i) > 0.0_dp) one_over_tau = min(budget**2*rho_f(i)*T_f(i)/K_s(i), one_over_tau)
        q_conv(i) = Fc_t(i) + (q_conv(i) - Fc_t(i))*exp(-dtf(i)*one_over_tau)
      end do
    end if

    if (top_wall_face .and. col_top) q_conv(n_c) = 0.0_dp

    ! bottom face: L_bottom split as the radiative (diffusion) and convective targets
    F_bot = L_bottom/(4.0_dp*pi*area_f(0))
    kap = kappa_ln(lnrho_f(0), lnT_f(0))
    Fr_t0 = -4.0_dp*a_rad*c_light*T_f(0)**3/(3.0_dp*kap*rho_f(0))*dTdr(0)
    if (.not. use_rad) Fr_t0 = 0.0_dp
    if (col_bottom) then
      if (abs(Fr_t0 + Fc_t(0)) > 0.0_dp) then
        q_conv(0) = F_bot*Fc_t(0)/(Fr_t0 + Fc_t(0))
      else
        q_conv(0) = F_bot
      end if
    end if

    do i = 1, n_c
      s1(i) = s1(i) + dtc(i)*(area_f(i-1)*q_conv(i-1) - area_f(i)*q_conv(i))/(vol_c(i)*rho(i)*T(i))
    end do
    if (top_cooling) then
      do i = 1, n_c
        s1(i) = s1(i) - dtc(i)*cool_c(i)/(rho(i)*T(i))
      end do
    end if

    q_rad = 0.0_dp
    Qrad = 0.0_dp
    if (use_rad) then
      call grey_column(dt, sc, F_bot - q_conv(0), rho, T, cv, lnrho, lnT, tau_above, q_rad, Qrad)
      do i = 1, n_c
        s1(i) = s1(i) + dtc(i)*Qrad(i)/(rho(i)*T(i))
      end do
    end if
  end subroutine transport_column

  subroutine grey_column(dt, sc, F0, rho, T, cv, lnrho, lnT, tau_above, q_rad, Qrad)
    real(dp), intent(in) :: dt, sc(1:), F0, tau_above
    real(dp), intent(in) :: rho(-1:), T(-1:), cv(-1:), lnrho(-1:), lnT(-1:)
    real(dp), intent(inout) :: q_rad(0:), Qrad(1:)
    real(dp) :: kr(n_c), dtf(n_c), J(n_c), rcv(n_c), top_coef, lnr_f(0:n_c), lnT_f(0:n_c), dTdr_f(0:n_c)
    integer :: i, n_act, i0
    n_act = 1
    do i = 1, n_c
      kr(i) = kappa_ln(lnrho(i), lnT(i))*rho(i)
      if (kr(i) > 0.0_dp) n_act = i
      kr(i) = max(kr(i), kr_floor)
      rcv(i) = rho(i)*cv(i)/sc(i)                ! (sc: the cell's step is sc dt)
    end do
    dtf = 0.0_dp
    do i = 1, n_act - 1
      dtf(i) = (r_c(i+1) - r_c(i))/(0.5_dp*(1.0_dp/kr(i) + 1.0_dp/kr(i+1)))
    end do
    top_coef = 4.0_dp*pi/(3.0_dp*(tau_above + kr(n_act)*(r_f(n_act) - r_c(n_act))) + 2.0_dp)
    if (top_cooling) top_coef = 0.0_dp            ! no radiative flux through the top face
    if (col_bottom .and. col_top .and. (grey_i0 <= 1 .or. grey_i0 >= n_act)) then
      call grey_core(n_act, area_f(0:n_act), vol_c(1:n_act), kr(1:n_act), dtf(1:n_act), T(1:n_act), &
          rcv(1:n_act), dt, F0, top_coef, q_rad(0:n_act), J(1:n_act), Qrad(1:n_act))
    else
      ! cells 1..grey_i0-1: the explicit diffusion flux F = -4 a c T^3/(3 kappa rho) dT/dr
      ! at faces 1..grey_i0-1 (F0 through face 0); radiation there is slow diffusion
      ! i0: the first cell of the implicit solve (only in the top block; none below)
      i0 = n_act + 1
      if (col_top) i0 = max(1, grey_i0)
      lnr_f = face_value(lnrho)
      lnT_f = face_value(lnT)
      dTdr_f = ddr_f(T)
      do i = 0, i0 - 1
        q_rad(i) = -4.0_dp*a_rad*c_light*exp(3.0_dp*lnT_f(i))/(3.0_dp*kappa_ln(lnr_f(i), lnT_f(i))*exp(lnr_f(i))) &
            *dTdr_f(i)
      end do
      if (col_bottom) q_rad(0) = F0
      do i = 1, i0 - 1
        Qrad(i) = -(area_f(i)*q_rad(i) - area_f(i-1)*q_rad(i-1))/vol_c(i)
      end do
      ! cells i0..n_act: the implicit grey solve, its bottom flux the one at face i0-1
      if (i0 <= n_act) call grey_core(n_act - i0 + 1, area_f(i0-1:n_act), vol_c(i0:n_act), kr(i0:n_act), &
          dtf(i0:n_act), T(i0:n_act), rcv(i0:n_act), dt, q_rad(i0-1), top_coef, &
          q_rad(i0-1:n_act), J(i0:n_act), Qrad(i0:n_act))
    end if
  end subroutine grey_column

  ! The cooling arrays of the m cells o+1..o+m (as bg_view)
  subroutine tr_view(o, m)
    integer, intent(in) :: o, m
    if (.not. associated(cool_c_g)) then
      cool_c_g => cool_c; cool_f_g => cool_f
    end if
    cool_c(1:m) => cool_c_g(o+1:o+m)
    cool_f(0:m) => cool_f_g(o:o+m)
  end subroutine tr_view

end module cz3d_transport_mod
