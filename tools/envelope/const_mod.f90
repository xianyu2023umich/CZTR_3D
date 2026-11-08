module const_mod
  ! Physical constants (cgs) and adopted solar parameters.
  use kinds_mod, only: dp
  implicit none
  private

  ! --- physical constants (cgs) ---
  real(dp), parameter, public :: G_grav   = 6.67430e-8_dp      ! cm^3 g^-1 s^-2
  real(dp), parameter, public :: sigma_SB = 5.670374e-5_dp     ! erg cm^-2 s^-1 K^-4
  real(dp), parameter, public :: c_light  = 2.99792458e10_dp   ! cm s^-1
  real(dp), parameter, public :: a_rad    = 4.0_dp*sigma_SB/c_light  ! radiation constant

  ! --- adopted solar parameters (IAU nominal values) ---
  real(dp), parameter, public :: Msun = 1.98892e33_dp   ! g
  real(dp), parameter, public :: Rsun = 6.95700e10_dp   ! cm
  real(dp), parameter, public :: Teff_sun = 5772.0_dp   ! K
  real(dp), parameter, public :: Lsun = 4.0_dp*3.14159265358979323846_dp*Rsun**2*sigma_SB*Teff_sun**4
  real(dp), parameter, public :: g_sun = G_grav*Msun/Rsun**2   ! cm s^-2, surface gravity

  ! --- target for the deep-CZ shooting condition (log10 specific entropy, erg g^-1 K^-1) ---
  ! logS_deep_target is the FreeEOS table's entropy at the deep-CZ point
  ! (logQ_target, logT_target), i.e. rho = 0.2085 g cm^-3, T = 2.326e6 K
  ! (read_eos.ipynb). Other EOS backends evaluate their own entropy at the
  ! same point, since entropy zero points differ between EOS codes.
  real(dp), parameter, public :: logS_deep_target = 9.252751889803044_dp
  real(dp), parameter, public :: logQ_target = -1.41394277_dp
  real(dp), parameter, public :: logT_target = 6.36653994_dp

  real(dp), parameter, public :: pi = 3.14159265358979323846_dp

end module const_mod
