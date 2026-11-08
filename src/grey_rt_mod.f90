module grey_rt_mod
  ! The grey Eddington moment solve, shared by CZ_1Dsimulation (transport_mod)
  ! and the TR code (tr_equilibrium.f90, --radiation=grey): one copy of the
  ! radiation scheme for the merged code. Depends only on the kind parameter.
  !
  !   F = -(4 pi/3) dJ/dtau        (faces; dtau between cell centres given)
  !   div F = 4 pi kappa rho eta (B - J),  B = sigma T^4/pi        (cells)
  !   eta = rho cv/(rho cv + 16 sigma T^3 kappa rho dt)   (implicit coupling)
  !
  ! with the flux F0 through face 0 and F = top_coef J at the top face.
  ! The heating Q (optional output) is -div F where a cell is optically thick
  ! and the equal local form 4 pi kappa rho eta (J - B) where it is thin
  ! (tau per cell < 1): in a thick cell J - B is below round-off, in a thin
  ! (nearly transparent) cell the face coefficient 4 pi/(3 dtau) is huge and
  ! -div F amplifies round-off in J. The two agree when the system is solved
  ! exactly.
  use kinds_mod, only: dp
  implicit none
  private
  real(dp), parameter :: sigma_SB = 5.670374e-5_dp            ! as CZ's const_mod
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  public :: grey_core

contains

  ! The grey moment solve on given arrays (also used by tests/test_grey_re):
  ! cells 1..nn with face areas af(0:nn) and volumes vc(1:nn) (per steradian;
  ! af = 1, vc = dr for plane-parallel), kappa rho kr, the optical depth
  ! dtau_f(i) between cell centres i and i+1, temperature T, rho cv rcv. F0 is
  ! the flux through face 0, top_coef = F(nn)/J(nn) the top condition (2 pi:
  ! vacuum at the top cell). dt = 0: instantaneous field; dt -> huge:
  ! radiative equilibrium div F = 0.
  subroutine grey_core(nn, af, vc, kr, dtau_f, T, rcv, dt, F0, top_coef, F, J, Q)
    integer, intent(in) :: nn
    real(dp), intent(in) :: af(0:nn), vc(1:nn), kr(1:nn), dtau_f(1:nn), T(1:nn), rcv(1:nn)
    real(dp), intent(in) :: dt, F0, top_coef
    real(dp), intent(inout) :: F(0:nn)
    real(dp), intent(out) :: J(1:nn)
    real(dp), intent(out), optional :: Q(1:nn)
    real(dp) :: D(0:nn), a(1:nn), b(1:nn), c(1:nn), rhs(1:nn), cpr(1:nn), dpr(1:nn), ex(1:nn), m, tau_cell
    integer :: i

    ! Face coefficient 4 pi/(3 dtau).
    D = 0.0_dp
    do i = 1, nn - 1
      D(i) = 4.0_dp*pi/(3.0_dp*dtau_f(i))
    end do
    do i = 1, nn
      ex(i) = 4.0_dp*pi*kr(i)*rcv(i)/(rcv(i) + 16.0_dp*sigma_SB*T(i)**3*kr(i)*dt)   ! 4 pi kappa rho eta
      a(i) = 0.0_dp
      c(i) = 0.0_dp
      b(i) = ex(i)
      rhs(i) = ex(i)*sigma_SB*T(i)**4/pi
      if (i > 1) then
        a(i) = -af(i-1)*D(i-1)/vc(i)
        b(i) = b(i) + af(i-1)*D(i-1)/vc(i)
      else
        rhs(i) = rhs(i) + af(0)*F0/vc(1)
      end if
      if (i < nn) then
        c(i) = -af(i)*D(i)/vc(i)
        b(i) = b(i) + af(i)*D(i)/vc(i)
      else
        b(i) = b(i) + af(nn)*top_coef/vc(nn)
      end if
    end do
    ! Thomas algorithm (diagonally dominant).
    cpr(1) = c(1)/b(1)
    dpr(1) = rhs(1)/b(1)
    do i = 2, nn
      m = b(i) - a(i)*cpr(i-1)
      cpr(i) = c(i)/m
      dpr(i) = (rhs(i) - a(i)*dpr(i-1))/m
    end do
    J(nn) = dpr(nn)
    do i = nn - 1, 1, -1
      J(i) = dpr(i) - cpr(i)*J(i+1)
    end do
    F(0) = F0
    do i = 1, nn - 1
      F(i) = -D(i)*(J(i+1) - J(i))
    end do
    F(nn) = top_coef*J(nn)
    if (present(Q)) then
      do i = 1, nn
        if (i == 1) then
          tau_cell = dtau_f(1)
        else if (i == nn) then
          tau_cell = dtau_f(nn-1)
        else
          tau_cell = 0.5_dp*(dtau_f(i-1) + dtau_f(i))
        end if
        if (tau_cell < 1.0_dp) then
          Q(i) = ex(i)*(J(i) - sigma_SB*T(i)**4/pi)
        else
          Q(i) = -(af(i)*F(i) - af(i-1)*F(i-1))/vc(i)
        end if
      end do
    end if
  end subroutine grey_core

end module grey_rt_mod
