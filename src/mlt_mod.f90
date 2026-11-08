module mlt_mod
  ! Mixing-length convective flux with the same relations as the static CZ
  ! code (CZ/srcFortran/mixing_length_mod.f90):
  !
  !   H_P = P/(rho g),  ell = alpha H_P,  Q = chi_T/chi_rho
  !   v_conv(x)^2 = g Q ell^2 x /(8 H_P)
  !   Gamma(x)   = 3 kappa rho^2 cp ell /(16 sigma T^3) * v_conv(x)
  !   grad - grad_ad = x + x/Gamma(x)
  !   F_conv     = rho cp T v_conv (ell/H_P) x
  !
  ! The static code knows grad_rad and solves the flux balance for x. Here the
  ! superadiabaticity y = grad - grad_ad is known from the current state (from
  ! the entropy gradient, y = -(H_P/cp) ds/dr), and x follows in closed form:
  ! with Gamma = A sqrt(x), y = u^2 + u/A for u = sqrt(x), so
  !   u = 2y / (1/A + sqrt(1/A^2 + 4y)).
  use kinds_mod, only: dp
  use const_mod, only: sigma_SB
  implicit none
  private

  public :: mlt_flux

contains

  ! F: convective flux [erg cm^-2 s^-1]. K_s: effective conductivity for
  ! entropy, dF/d(-ds/dr), used to set the hyperbolic relaxation time.
  pure subroutine mlt_flux(y, rho, T, P, g, kappa, cp, Q, alpha, F, K_s)
    real(dp), intent(in) :: y, rho, T, P, g, kappa, cp, Q, alpha
    real(dp), intent(out) :: F, K_s
    real(dp) :: Hp, ell, A, u, x, v, dF_dx, dx_dy

    F = 0.0_dp
    K_s = 0.0_dp
    if (y <= 0.0_dp .or. Q <= 0.0_dp) return

    Hp = P/(rho*g)
    ell = alpha*Hp
    A = (3.0_dp*kappa*rho**2*cp*ell/(16.0_dp*sigma_SB*T**3)) * sqrt(g*Q*ell**2/(8.0_dp*Hp))
    u = 2.0_dp*y/(1.0_dp/A + sqrt(1.0_dp/A**2 + 4.0_dp*y))
    x = u*u
    v = sqrt(g*Q*ell**2*x/(8.0_dp*Hp))
    F = rho*cp*T*v*(ell/Hp)*x

    ! F ~ x^(3/2); dy/dx = 1 + 1/(2 A u); dy/d(-ds/dr) = H_P/cp.
    dF_dx = 1.5_dp*F/x
    dx_dy = 1.0_dp/(1.0_dp + 1.0_dp/(2.0_dp*A*u))
    K_s = dF_dx*dx_dy*Hp/cp
  end subroutine mlt_flux

end module mlt_mod
