program grid_gen
  ! The radial grid of CZ_3D (and CZTR): the cell density n(r) of
  ! tools/grid/grid_mod.f90 (5 cells per H_p, the atmosphere law above the
  ! photosphere, a smooth 30 km floor, and optionally dr <= aspect * r dtheta)
  ! integrated into a smooth map r(xi), from the base of the profile (cut at
  ! the pressure p_base) to z_top_km above the photosphere.
  !
  ! Usage: build/grid_gen.x <static profile> <grid namelist> <grid file>
  !   static profile: tools/static_bg.f90's output (r T P rho Mr/Msun ...,
  !                   from the photosphere down; the last line deepest)
  !   grid namelist:  &grid (see grid_mod: cells_per_Hp, dr_min_km, horiz_deg,
  !                   aspect, p_base, ...); missing items keep their defaults
  !   grid file:      '#' comments, n, then 'i r_c' for cells -1..n+2, then
  !                   'i r_f' for faces 0..n [cm] (as CZTR's export_grids)
  use kinds_mod, only: dp
  use grid_mod
  implicit none
  real(dp), parameter :: G_grav = 6.67430e-8_dp, Msun = 1.98892e33_dp, Rsun = 6.957e10_dp
  real(dp), parameter :: r_show(4) = [0.72_dp, 0.85_dp, 0.96_dp, 0.99_dp]
  type(grid_params) :: par
  character(len=512) :: prof_file, nml_file, out_file, line
  real(dp), allocatable :: rp(:), Pp(:), rhop(:), Mp(:), Hp(:), r_c(:), r_f(:)
  real(dp) :: tmp(5), w
  integer :: u, ios, m, n, k, i, j

  if (command_argument_count() < 3) then
    print '(A)', 'usage: grid_gen.x <static profile> <grid namelist> <grid file>'
    stop 1
  end if
  call get_command_argument(1, prof_file)
  call get_command_argument(2, nml_file)
  call get_command_argument(3, out_file)
  call grid_read_nml(nml_file, par)

  ! the profile, top first in the file; keep P <= p_base
  allocate(rp(40000), Pp(40000), rhop(40000), Mp(40000))
  open(newunit=u, file=trim(prof_file), status='old', action='read')
  m = 0
  do
    read(u, '(A)', iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
    read(line, *, iostat=ios) tmp
    if (ios /= 0) cycle
    if (m > 0) then
      if (tmp(3) > par%p_base) then
        ! the base point at exactly p_base: linear in ln P
        w = log(par%p_base/Pp(m))/log(tmp(3)/Pp(m))
        m = m + 1
        rp(m) = rp(m-1) + w*(tmp(1) - rp(m-1))
        Pp(m) = par%p_base
        rhop(m) = exp(log(rhop(m-1)) + w*(log(tmp(4)) - log(rhop(m-1))))
        Mp(m) = Mp(m-1) + w*(tmp(5)*Msun - Mp(m-1))
        exit
      end if
    end if
    m = m + 1
    rp(m) = tmp(1);  Pp(m) = tmp(3);  rhop(m) = tmp(4);  Mp(m) = tmp(5)*Msun
    if (Pp(m) == par%p_base) exit
  end do
  close(u)
  if (abs(Pp(m)/par%p_base - 1.0_dp) > 1.0e-9_dp) stop 'grid_gen: the profile does not reach p_base'
  rp(1:m) = rp(m:1:-1);  Pp(1:m) = Pp(m:1:-1);  rhop(1:m) = rhop(m:1:-1);  Mp(1:m) = Mp(m:1:-1)
  allocate(Hp(m))
  Hp = Pp(1:m)/(rhop(1:m)*G_grav*Mp(1:m)/rp(1:m)**2)

  call grid_build(par, rp(1:m), Hp, n, r_c, r_f)
  k = count(r_c(1:n) <= rp(m))

  open(newunit=u, file=trim(out_file), status='replace', action='write')
  write(u, '(4A)') '# merged grid, all cells: tools/grid_gen.f90 ', trim(prof_file), ' ', trim(nml_file)
  write(u, '(A,F6.2,A,ES10.3,A,F7.4,A,F6.3)') '# cells_per_Hp ', par%cells_per_Hp, '  dr_min_km ', par%dr_min_km, &
      '  horiz_deg ', par%horiz_deg, '  aspect ', par%aspect
  write(u, '(A,I0,A)') '# merged cells 1..', n, '; then n; cells -1..n+2 (i r_c); faces 0..n (i r_f) [cm]'
  write(u, '(I0)') n
  do i = -1, n + 2
    write(u, '(I6,ES26.17)') i, r_c(i)
  end do
  do i = 0, n
    write(u, '(I6,ES26.17)') i, r_f(i)
  end do
  close(u)
  print '(A,A,A,I0,A,I0,A,F8.5,A)', 'grid_gen: ', trim(out_file), ': ', n, ' cells, ', k, &
      ' to the photosphere; base face r = ', r_f(0)/Rsun, ' R_sun'
  do j = 1, 4
    w = r_show(j)*Rsun
    i = minloc(abs(r_c(1:n) - w), 1)
    print '(A,F6.3,A,F9.1,A)', '  r = ', r_c(i)/Rsun, ' R_sun: dr = ', (r_f(i) - r_f(i-1))/1.0e5_dp, ' km'
  end do
end program grid_gen
