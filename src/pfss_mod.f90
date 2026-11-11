module pfss_mod
  ! Potential field source surface (PFSS) pieces: orthonormal associated
  ! Legendre functions and the radial factors of the PFSS solution.
  !
  ! B_r(R0, theta, phi) = sum_lm p_lm(theta) (c_lm cos(m phi) + s_lm sin(m phi))
  ! with Y_lm = p_lm e^(i m phi), int |Y_lm|^2 dOmega = 1 (no Condon-Shortley
  ! phase). The potential Phi_l ~ (r/r_ss)^l - (r_ss/r)^(l+1) (zero at the
  ! source surface r_ss), B = -grad Phi, scaled so that B_r(R0) is the given
  ! field:
  !   B_r     = sum_l fr_l(r) S_l,  S_l = sum_m p_lm (c cos + s sin)
  !   B_theta = sum_l gh_l(r) T_l,  T_l = sum_m dp_lm/dtheta (c cos + s sin)
  !   B_phi   = sum_l gh_l(r) U_l,  U_l = sum_m p_lm m (-c sin + s cos)/sin(theta)
  ! with q = R0/r_ss, x = R0/r, y = r/r_ss, D_l = (l+1) + l q^(2l+1):
  !   fr_l = [(l+1) x^(l+2) + l y^(l-1) q^(l+2)]/D_l      (1 at R0)
  !   gh_l = x [q^(l+1) y^l - x^(l+1)]/D_l                 (0 at r_ss)
  ! (every power <= 1: no overflow at high l). For l = 1 this is
  ! cz3d_app_mod's seed_dipole.
  use kinds_mod, only: dp
  implicit none
  private
  public :: pfss_legendre, pfss_radial

contains

  ! p(l, m) and, if present, dp/dtheta (l, m) for 0 <= m <= l <= lmx at the
  ! colatitude with cos ct and sin st (st > 0); the standard stable
  ! recurrences. Values below 1e-280 (high m near the poles) are left 0.
  subroutine pfss_legendre(lmx, ct, st, p, dpt)
    integer, intent(in) :: lmx
    real(dp), intent(in) :: ct, st
    real(dp), intent(out) :: p(0:lmx, 0:lmx)
    real(dp), intent(out), optional :: dpt(0:lmx, 0:lmx)
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp) :: pmm, sti
    integer :: l, m, mtop
    p = 0.0_dp
    pmm = 1.0_dp/sqrt(4.0_dp*pi)
    mtop = lmx
    do m = 0, lmx
      if (m > 0) pmm = sqrt(real(2*m + 1, dp)/real(2*m, dp))*st*pmm
      if (abs(pmm) < 1.0e-280_dp) then
        mtop = m - 1
        exit
      end if
      p(m, m) = pmm
      if (m < lmx) p(m+1, m) = sqrt(real(2*m + 3, dp))*ct*pmm
      do l = m + 2, lmx
        p(l, m) = sqrt(real(4*l*l - 1, dp)/real(l*l - m*m, dp)) &
            *(ct*p(l-1, m) - sqrt(real((l-1)**2 - m*m, dp)/real(4*(l-1)**2 - 1, dp))*p(l-2, m))
      end do
    end do
    if (.not. present(dpt)) return
    ! sin(theta) dp_lm/dtheta = l cos(theta) p_lm - sqrt((2l+1)/(2l-1) (l^2 - m^2)) p_(l-1)m
    dpt = 0.0_dp
    sti = 1.0_dp/max(st, 1.0e-12_dp)
    do m = 0, mtop
      do l = max(m, 1), lmx
        dpt(l, m) = sti*(l*ct*p(l, m) - sqrt(real(2*l + 1, dp)/real(2*l - 1, dp)*real(l*l - m*m, dp))*p(l-1, m))
      end do
    end do
  end subroutine pfss_legendre

  ! fr(0:lmx), gh(0:lmx) at radius r (R0 <= r <= r_ss)
  subroutine pfss_radial(lmx, r, R0, rss, fr, gh)
    integer, intent(in) :: lmx
    real(dp), intent(in) :: r, R0, rss
    real(dp), intent(out) :: fr(0:lmx), gh(0:lmx)
    real(dp) :: q, x, y, xl, yl, ql, d
    integer :: l
    q = R0/rss
    x = R0/r
    y = r/rss
    xl = x        ! x^(l+1)
    yl = 1.0_dp   ! y^l
    ql = q        ! q^(l+1)
    do l = 0, lmx
      d = real(l + 1, dp) + real(l, dp)*ql*ql/q            ! (l+1) + l q^(2l+1)
      if (l == 0) then
        fr(l) = x*x
      else
        fr(l) = (real(l + 1, dp)*xl*x + real(l, dp)*(yl/y)*ql*q)/d
      end if
      gh(l) = x*(ql*yl - xl)/d
      xl = xl*x
      yl = yl*y
      ql = ql*q
    end do
  end subroutine pfss_radial

end module pfss_mod
