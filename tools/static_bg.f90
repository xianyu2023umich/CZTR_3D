program static_bg
  ! The static CZ envelope (CZ/srcFortran: hydrostatic + MLT + radiative
  ! diffusion, integrated down from tau = 0.01 at R_sun) as a background for
  ! cz3d.x, deeper than CZ's own output: alpha is found exactly as CZ's
  ! main.f90 does (shooting to the deep-CZ entropy target at P = 6e13 with
  ! eos_hhe), then the envelope is integrated once more with that alpha down
  ! to p_deep, below the 3D grid's bottom ghost cells.
  !
  ! Usage (from CZ_3D): build/static_bg.x out.dat p_deep --eos=hhe     (built and run by `make data`)
  ! The envelope modules are copies of CZ's (tools/envelope/).
  use kinds_mod, only: dp
  use const_mod, only: Msun, Rsun, Teff_sun, Lsun
  use cz_eos_mod, only: cz_eos_t, cz_eos_from_args, cz_eos_logS_target, cz_eos_tag
  use opacity_table_mod, only: opac_table_t, opac_load
  use envelope_mod, only: envelope_profile_t, surface_start, integrate_envelope
  implicit none
  type(cz_eos_t) :: eos
  type(opac_table_t) :: opac
  type(envelope_profile_t) :: prof
  real(dp), parameter :: TAU_START = 0.01_dp, P_END = 6.0e13_dp
  integer, parameter :: N_STEPS = 4000, MAXIT = 40
  real(dp), parameter :: TOL = 1.0e-6_dp
  real(dp) :: a1, a2, a3, r1, r2, r3, alpha, logS_target, p_deep
  character(len=256) :: out, arg
  integer :: it, u, j, n_deep

  call get_command_argument(1, out)
  p_deep = 2.5e14_dp
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg)
    read(arg, *) p_deep
  end if
  ! eos_hhe, X = 0.75, Z = 0.02 (as CZ's solar_profile_hhe.dat and the 3D EOS
  ! table): the backend comes from --eos=hhe on the command line
  call cz_eos_from_args(eos, 'none')
  if (trim(cz_eos_tag(eos)) /= '_hhe') stop 'static_bg: run with --eos=hhe (the FreeEOS backend is not supported here)'
  logS_target = cz_eos_logS_target(eos)
  call opac_load(opac, 'data/opacity_merged_x70_z02.dat')

  a1 = 1.2_dp;  a2 = 2.4_dp
  r1 = residual(a1, P_END, N_STEPS);  r2 = residual(a2, P_END, N_STEPS)
  alpha = a2
  do it = 1, MAXIT
    if (abs(r2 - r1) < 1.0e-14_dp) exit
    a3 = max(0.3_dp, min(4.0_dp, a2 - r2*(a2 - a1)/(r2 - r1)))
    r3 = residual(a3, P_END, N_STEPS)
    alpha = a3
    if (abs(r3) < TOL) exit
    a1 = a2;  r1 = r2;  a2 = a3;  r2 = r3
  end do
  print '(A,F10.6)', 'alpha = ', alpha
  ! the same step in ln P as CZ's main down to P_END, continued to p_deep
  n_deep = nint(N_STEPS*log(p_deep/prof_p0())/log(P_END/prof_p0()))
  r3 = residual(alpha, p_deep, n_deep)
  print '(A,ES12.4,A,I0,A,F9.5)', 'deep integration to P = ', p_deep, ': ', prof%n, ' points, base r/Rsun = ', &
      prof%r(prof%n)/Rsun
  open(newunit=u, file=trim(out), status='replace', action='write')
  write(u, '(A,F10.6)') '# static CZ envelope (CZ/srcFortran envelope_mod, eos_hhe X=0.75 Z=0.02), alpha = ', alpha
  write(u, '(A)') '# r[cm]  T[K]  P[dyn/cm2]  rho[g/cm3]  Mr/Msun  grad  grad_ad  is_convective'
  do j = 1, prof%n
    write(u, '(7(ES24.16,1X),L1)') prof%r(j), prof%T(j), prof%P(j), prof%rho(j), prof%Mr(j)/Msun, prof%grad(j), &
        prof%grad_ad(j), prof%is_convective(j)
  end do
  close(u)
  print *, 'wrote ', trim(out)

contains

  real(dp) function prof_p0()
    real(dp) :: T0, P0, rho0
    call surface_start(eos, opac, Msun, Teff_sun, Rsun, TAU_START, T0, P0, rho0)
    prof_p0 = P0
  end function prof_p0

  real(dp) function residual(alpha, p_end, n)
    real(dp), intent(in) :: alpha, p_end
    integer, intent(in) :: n
    real(dp) :: T0, P0, rho0
    call surface_start(eos, opac, Msun, Teff_sun, Rsun, TAU_START, T0, P0, rho0)
    call integrate_envelope(eos, opac, Msun, Lsun, alpha, T0, P0, Rsun, log(p_end), n, prof)
    residual = prof%logS_final - logS_target
  end function residual

end program static_bg
