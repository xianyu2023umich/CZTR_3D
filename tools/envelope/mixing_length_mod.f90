module mixing_length_mod
  ! Standard local mixing-length closure (Boehm-Vitense 1958; see textbook
  ! treatments e.g. Kippenhahn, Weigert & Weiss 2012 Ch.7, or Hansen, Kawaler
  ! & Trimble 2004 Ch.5). Derivation used here, for the record:
  !
  !   H_P   = P/(rho g)                       pressure scale height
  !   ell   = alpha * H_P                     mixing length
  !   Q     = (dlnrho/dlnT)_P (sign-flipped)  = chiT/chiRho from the EOS
  !
  ! A bubble displaced by x=ell picks up, from buoyancy work integrated over
  ! its path (standard mixing-length energetics):
  !   v_conv(x)^2 = (g Q ell^2)/(8 H_P) * x ,   x >= 0
  ! where x == grad - grad_element is the (unknown) local superadiabatic
  ! temperature-gradient excess driving convection.
  !
  ! The bubble loses heat radiatively while it travels; the ratio of the
  ! radiative-diffusion cooling time to the convective turnover time defines
  ! the convective efficiency:
  !   Gamma(x) = (3 kappa rho^2 Cp ell)/(16 sigma_SB T^3) * v_conv(x)
  ! with the standard closure (Gamma->inf: fully efficient/adiabatic bubble;
  ! Gamma->0: fully inefficient, bubble thermalises instantly):
  !   x = [Gamma/(1+Gamma)] * (grad - grad_ad)
  !   =>  grad = grad_ad + x + x/Gamma(x)
  !
  ! Flux conservation (grad_rad is, by construction, the gradient that would
  ! carry the whole local flux radiatively, so F_rad(grad) = F_tot*grad/grad_rad
  ! exactly in the diffusion approximation) then gives one equation in the one
  ! unknown x:
  !   F_conv(x) = rho Cp T v_conv(x) (ell/H_P) x
  !   G(x) = F_conv(x) - F_tot*(grad_rad - grad(x))/grad_rad = 0
  ! G(0+) < 0 and G(x) -> +inf as x -> infinity (F_conv ~ x^1.5 dominates the
  ! linear term), so G has a unique positive root; solved here by bracketed
  ! bisection rather than the classical cubic-in-xi reduction, to avoid
  ! hand-deriving that algebra. O(1) geometric prefactors in the energetics
  ! above are a known convention-dependent ambiguity in mixing-length theory;
  ! since alpha itself is calibrated (not taken from the literature), this
  ! ambiguity is absorbed into the calibrated alpha and does not need to match
  ! any particular textbook's exact constants.
  use kinds_mod, only: dp
  use const_mod, only: sigma_SB
  implicit none
  private

  public :: mlt_result_t, solve_mlt

  type :: mlt_result_t
    real(dp) :: grad        ! actual dlnT/dlnP
    real(dp) :: F_conv      ! convective flux, erg cm^-2 s^-1
    real(dp) :: v_conv      ! cm/s
    real(dp) :: gamma_eff   ! convective efficiency parameter
    logical  :: is_convective
  end type mlt_result_t

contains

  subroutine solve_mlt(rho, T, P, grav, kappa, Cp, Q, alpha, grad_ad, grad_rad, res)
    real(dp), intent(in) :: rho, T, P, grav, kappa, Cp, Q, alpha, grad_ad, grad_rad
    type(mlt_result_t), intent(out) :: res

    real(dp) :: Hp, ell, Ftot, x_lo, x_hi, G_lo, G_hi, x_mid, G_mid
    integer :: it
    integer, parameter :: MAXIT_BRACKET = 60, MAXIT_BISECT = 80
    real(dp), parameter :: REL_TOL = 1.0e-8_dp

    if (grad_rad <= grad_ad) then
      res%grad = grad_rad
      res%F_conv = 0.0_dp
      res%v_conv = 0.0_dp
      res%gamma_eff = 0.0_dp
      res%is_convective = .false.
      return
    end if

    Hp  = P/(rho*grav)
    ell = alpha*Hp
    ! Ftot = F_rad(grad_rad) = local total flux, via the diffusion-approx
    ! relation F_rad(grad) = 4 a_rad c g T^4 grad /(3 kappa P); see structure_mod
    ! for the grad_rad formula this is paired with.
    Ftot = flux_of_grad(grad_rad, T, P, grav, kappa)

    x_lo = 1.0e-10_dp*(grad_rad - grad_ad)
    x_hi = (grad_rad - grad_ad)
    G_lo = residual(x_lo)
    G_hi = residual(x_hi)

    it = 0
    do while (G_hi < 0.0_dp .and. it < MAXIT_BRACKET)
      x_hi = x_hi*2.0_dp
      G_hi = residual(x_hi)
      it = it + 1
    end do

    do it = 1, MAXIT_BISECT
      x_mid = 0.5_dp*(x_lo + x_hi)
      G_mid = residual(x_mid)
      if (G_mid < 0.0_dp) then
        x_lo = x_mid; G_lo = G_mid
      else
        x_hi = x_mid; G_hi = G_mid
      end if
      if (abs(x_hi - x_lo) < REL_TOL*x_hi) exit
    end do

    x_mid = 0.5_dp*(x_lo + x_hi)
    call fill_result(x_mid, res)

  contains

    function flux_of_grad(grad, TT, PP, gg, kap) result(F)
      use const_mod, only: a_rad, c_light
      real(dp), intent(in) :: grad, TT, PP, gg, kap
      real(dp) :: F
      F = (4.0_dp*a_rad*c_light*gg*TT**4*grad)/(3.0_dp*kap*PP)
    end function flux_of_grad

    function v_conv_of_x(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: v
      v = sqrt(max(0.0_dp, (grav*Q*ell**2)/(8.0_dp*Hp)*x))
    end function v_conv_of_x

    function gamma_of_x(x) result(gam)
      real(dp), intent(in) :: x
      real(dp) :: gam
      gam = (3.0_dp*kappa*rho**2*Cp*ell)/(16.0_dp*sigma_SB*T**3) * v_conv_of_x(x)
    end function gamma_of_x

    function grad_of_x(x) result(gr)
      real(dp), intent(in) :: x
      real(dp) :: gr, gam
      gam = gamma_of_x(x)
      gr = grad_ad + x + x/gam
    end function grad_of_x

    function residual(x) result(G)
      real(dp), intent(in) :: x
      real(dp) :: G, v, Fconv_x, grad_x
      v = v_conv_of_x(x)
      Fconv_x = rho*Cp*T*v*(ell/Hp)*x
      grad_x = grad_of_x(x)
      G = Fconv_x - Ftot*(grad_rad - grad_x)/grad_rad
    end function residual

    subroutine fill_result(x, r)
      real(dp), intent(in) :: x
      type(mlt_result_t), intent(out) :: r
      r%v_conv = v_conv_of_x(x)
      r%gamma_eff = gamma_of_x(x)
      r%grad = grad_of_x(x)
      r%F_conv = rho*Cp*T*r%v_conv*(ell/Hp)*x
      r%is_convective = .true.
    end subroutine fill_result

  end subroutine solve_mlt

end module mixing_length_mod
