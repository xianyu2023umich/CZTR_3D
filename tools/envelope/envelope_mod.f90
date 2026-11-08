module envelope_mod
  ! Surface starting point (grey Eddington T-tau relation + self-consistent
  ! hydrostatic pressure using the real EOS/opacity tables) and an RK4
  ! integrator in ln(P), stepping down from the photosphere into the deep
  ! convection zone. NOTE: the near-surface radiative transfer is still the
  ! simplified grey-Eddington starting condition, not yet the tunable-n_mu
  ! explicit transfer discussed earlier -- that refinement is a follow-up,
  ! confined to the outermost few points where it matters.
  use kinds_mod, only: dp
  use const_mod, only: G_grav
  use cz_eos_mod, only: cz_eos_t, cz_eos_state_t, cz_eos_PT, cz_eos_rhoT
  use opacity_table_mod, only: opac_table_t, opac_get_kappa_rhoT
  use structure_mod, only: structure_point_t, eval_rhs
  implicit none
  private

  public :: envelope_profile_t, surface_start, integrate_envelope

  type :: envelope_profile_t
    integer :: n = 0
    real(dp), allocatable :: r(:), T(:), P(:), rho(:), kappa(:), Mr(:)
    real(dp), allocatable :: grad(:), grad_ad(:), grad_rad(:), F_conv(:)
    logical,  allocatable :: is_convective(:)
    real(dp) :: logS_final
  end type envelope_profile_t

contains

  subroutine surface_start(eos, opac, M, Teff, R, tau_start, T0, P0, rho0)
    type(cz_eos_t),     intent(in) :: eos
    type(opac_table_t), intent(in) :: opac
    real(dp), intent(in)  :: M, Teff, R, tau_start
    real(dp), intent(out) :: T0, P0, rho0

    real(dp) :: grav, kappa, rho_guess
    logical :: oot
    integer :: it

    T0 = Teff*(0.75_dp*(tau_start + 2.0_dp/3.0_dp))**0.25_dp
    grav = G_grav*M/R**2
    rho_guess = 1.0e-8_dp

    do it = 1, 10
      kappa = opac_get_kappa_rhoT(opac, rho_guess, T0, oot)
      P0 = grav*tau_start/kappa
      rho_guess = rho_PT(P0, T0, rho_guess)
    end do
    rho0 = rho_guess

  contains

    function rho_PT(P, T, guess) result(rho)
      real(dp), intent(in) :: P, T, guess
      real(dp) :: rho
      type(cz_eos_state_t) :: st
      st = cz_eos_PT(eos, P, T, guess)
      rho = st%rho
    end function rho_PT

  end subroutine surface_start

  subroutine integrate_envelope(eos, opac, M_total, L, alpha, T0, P0, R0, lnP_end, n_steps, prof)
    type(cz_eos_t),     intent(in) :: eos
    type(opac_table_t), intent(in) :: opac
    real(dp), intent(in) :: M_total, L, alpha, T0, P0, R0, lnP_end
    integer,  intent(in) :: n_steps
    type(envelope_profile_t), intent(out) :: prof

    real(dp) :: lnP0, dlnP, lnP_i, r, T, Mr, P, rho_guess
    real(dp) :: k1r, k1T, k1M, k2r, k2T, k2M, k3r, k3T, k3M, k4r, k4T, k4M
    real(dp) :: r2, T2, Mr2, r3, T3, Mr3, r4, T4, Mr4, Pm, Pf
    type(structure_point_t) :: pt
    integer :: i

    allocate(prof%r(n_steps+1), prof%T(n_steps+1), prof%P(n_steps+1), prof%rho(n_steps+1))
    allocate(prof%kappa(n_steps+1), prof%grad(n_steps+1), prof%grad_ad(n_steps+1))
    allocate(prof%grad_rad(n_steps+1), prof%F_conv(n_steps+1), prof%is_convective(n_steps+1))
    allocate(prof%Mr(n_steps+1))

    lnP0 = log(P0)
    dlnP = (lnP_end - lnP0)/real(n_steps, dp)

    r = R0
    T = T0
    Mr = M_total   ! essentially exact: negligible mass lies above tau_start
    P = P0
    rho_guess = 1.0e-8_dp

    call eval_rhs(eos, opac, Mr, L, alpha, r, T, P, rho_guess, k1r, k1T, k1M, pt)
    rho_guess = pt%rho
    call store(1, r, T, P, Mr, pt)

    do i = 1, n_steps
      lnP_i = lnP0 + real(i-1, dp)*dlnP
      Pm = exp(lnP_i + 0.5_dp*dlnP)
      Pf = exp(lnP_i + dlnP)

      call eval_rhs(eos, opac, Mr, L, alpha, r, T, P, rho_guess, k1r, k1T, k1M, pt)
      rho_guess = pt%rho

      r2 = r + 0.5_dp*dlnP*k1r;  T2 = T + 0.5_dp*dlnP*k1T;  Mr2 = Mr + 0.5_dp*dlnP*k1M
      call eval_rhs(eos, opac, Mr2, L, alpha, r2, T2, Pm, rho_guess, k2r, k2T, k2M, pt)

      r3 = r + 0.5_dp*dlnP*k2r;  T3 = T + 0.5_dp*dlnP*k2T;  Mr3 = Mr + 0.5_dp*dlnP*k2M
      call eval_rhs(eos, opac, Mr3, L, alpha, r3, T3, Pm, rho_guess, k3r, k3T, k3M, pt)

      r4 = r + dlnP*k3r;  T4 = T + dlnP*k3T;  Mr4 = Mr + dlnP*k3M
      call eval_rhs(eos, opac, Mr4, L, alpha, r4, T4, Pf, rho_guess, k4r, k4T, k4M, pt)

      r  = r  + (dlnP/6.0_dp)*(k1r + 2.0_dp*k2r + 2.0_dp*k3r + k4r)
      T  = T  + (dlnP/6.0_dp)*(k1T + 2.0_dp*k2T + 2.0_dp*k3T + k4T)
      Mr = Mr + (dlnP/6.0_dp)*(k1M + 2.0_dp*k2M + 2.0_dp*k3M + k4M)
      P = Pf

      call eval_rhs(eos, opac, Mr, L, alpha, r, T, P, rho_guess, k1r, k1T, k1M, pt)
      rho_guess = pt%rho
      call store(i+1, r, T, P, Mr, pt)
    end do

    prof%n = n_steps + 1
    block
      type(cz_eos_state_t) :: st_final
      st_final = cz_eos_rhoT(eos, prof%rho(prof%n), prof%T(prof%n))
      prof%logS_final = st_final%logS
    end block

    print '(A,F8.5,A)', '  Mr/M_total at domain base = ', prof%Mr(prof%n)/M_total, &
         '  (1.0 would be the M_r=M approximation we replaced)'

  contains

    subroutine store(idx, rr, TT, PP, MMr, ptl)
      integer, intent(in) :: idx
      real(dp), intent(in) :: rr, TT, PP, MMr
      type(structure_point_t), intent(in) :: ptl
      prof%r(idx) = rr; prof%T(idx) = TT; prof%P(idx) = PP; prof%rho(idx) = ptl%rho
      prof%Mr(idx) = MMr
      prof%kappa(idx) = ptl%kappa; prof%grad(idx) = ptl%grad; prof%grad_ad(idx) = ptl%grad_ad
      prof%grad_rad(idx) = ptl%grad_rad; prof%F_conv(idx) = ptl%F_conv
      prof%is_convective(idx) = ptl%is_convective
    end subroutine store

  end subroutine integrate_envelope

end module envelope_mod
