module structure_mod
  ! Assembles the envelope structure equations, independent variable ln(P):
  !   dr/dlnP  = -H_P                  (hydrostatic equilibrium)
  !   dT/dlnP  = T * grad              (grad = grad_rad if stable,
  !                                      grad = MLT solution if unstable)
  !   dMr/dlnP = -4*pi*r^2*rho*H_P     (mass continuity, dMr/dr=4*pi*r^2*rho,
  !                                      combined with dr/dlnP=-H_P)
  ! Mr is the mass interior to r, integrated alongside (r,T) rather than
  ! approximated as the fixed total mass M. Only L is still held fixed
  ! (L_r = L, no nuclear energy generation in this envelope -- an excellent
  ! approximation, unlike M_r=M which is off by ~2.5% at the CZ base).
  use kinds_mod, only: dp
  use const_mod, only: pi, a_rad, c_light, G_grav
  use cz_eos_mod, only: cz_eos_t, cz_eos_state_t, cz_eos_PT
  use opacity_table_mod, only: opac_table_t, opac_get_kappa_rhoT
  use mixing_length_mod, only: mlt_result_t, solve_mlt
  implicit none
  private

  public :: structure_point_t, eval_rhs

  type :: structure_point_t
    real(dp) :: rho, kappa, grad_rad, grad_ad, grad, F_conv, v_conv, Hp
    logical  :: is_convective, opacity_out_of_table
  end type structure_point_t

contains

  pure function grad_rad_of(kappa, L, Mr, P, T) result(gr)
    real(dp), intent(in) :: kappa, L, Mr, P, T
    real(dp) :: gr
    gr = 3.0_dp*kappa*L*P/(16.0_dp*pi*a_rad*c_light*G_grav*Mr*T**4)
  end function grad_rad_of

  subroutine eval_rhs(eos, opac, Mr, L, alpha, r, T, P, rho_guess, &
                       dr_dlnP, dT_dlnP, dMr_dlnP, pt)
    type(cz_eos_t),     intent(in)  :: eos
    type(opac_table_t), intent(in)  :: opac
    real(dp), intent(in)  :: Mr, L, alpha, r, T, P, rho_guess
    real(dp), intent(out) :: dr_dlnP, dT_dlnP, dMr_dlnP
    type(structure_point_t), intent(out) :: pt

    real(dp) :: grav, Q
    type(cz_eos_state_t) :: st
    type(mlt_result_t) :: mlt

    st = cz_eos_PT(eos, P, T, rho_guess)
    pt%rho = st%rho
    pt%kappa = opac_get_kappa_rhoT(opac, pt%rho, T, pt%opacity_out_of_table)

    grav = G_grav*Mr/r**2
    pt%Hp = P/(pt%rho*grav)
    pt%grad_ad = st%grad_ad
    pt%grad_rad = grad_rad_of(pt%kappa, L, Mr, P, T)
    Q = st%chiT/st%chiRho

    call solve_mlt(pt%rho, T, P, grav, pt%kappa, st%Cp, Q, alpha, pt%grad_ad, pt%grad_rad, mlt)
    pt%grad = mlt%grad
    pt%F_conv = mlt%F_conv
    pt%v_conv = mlt%v_conv
    pt%is_convective = mlt%is_convective

    dr_dlnP = -pt%Hp
    dT_dlnP = T*pt%grad
    dMr_dlnP = -4.0_dp*pi*r**2*pt%rho*pt%Hp
  end subroutine eval_rhs

end module structure_mod
