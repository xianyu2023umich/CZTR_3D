module thermo_mod
  ! Invert the EOS table for density given (P, T): P_total is assumed equal
  ! to the table's gas pressure alone (radiation pressure is <0.2% of P_gas
  ! even at the base of the CZ, ~2e6 K, and far smaller near the surface --
  ! negligible at the precision this project targets, and the FreeEOS table
  ! only reports P_gas in any case).
  use kinds_mod, only: dp
  use eos_table_mod, only: eos_table_t, eos_state_t, eos_get_rhoT
  implicit none
  private
  public :: rho_of_PT

contains

  function rho_of_PT(eos, P, T, rho_guess, st_out) result(rho)
    type(eos_table_t), intent(in) :: eos
    real(dp), intent(in) :: P, T, rho_guess
    type(eos_state_t), intent(out), optional :: st_out
    real(dp) :: rho

    real(dp) :: logrho, logP_target, resid, dlogrho, step
    type(eos_state_t) :: st
    integer :: it
    integer, parameter :: MAXIT = 40
    real(dp), parameter :: TOL = 1.0e-10_dp
    real(dp), parameter :: MAX_STEP = 0.3_dp   ! clamp Newton step in log10(rho)

    logrho = log10(rho_guess)
    logP_target = log10(P)

    do it = 1, MAXIT
      st = eos_get_rhoT(eos, 10.0_dp**logrho, T)
      resid = st%logPgas - logP_target
      if (abs(resid) < TOL) exit
      ! Newton step using chiRho = dlnP/dlnrho|T
      step = -resid/max(st%chiRho, 1.0e-3_dp)
      step = max(-MAX_STEP, min(MAX_STEP, step))
      dlogrho = step
      logrho = logrho + dlogrho
    end do

    rho = 10.0_dp**logrho
    if (present(st_out)) st_out = st
  end function rho_of_PT

end module thermo_mod
