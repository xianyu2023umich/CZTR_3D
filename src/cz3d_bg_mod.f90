module cz3d_bg_mod
  ! The 1D background, initial state and radial metric of the 3D run, by
  ! radial cell index, taken from the merged 1D grid and a CZTR state file
  ! (init_full.dat format, e.g. step 1's final_init.dat) exactly as CZTR
  ! builds them (cztr_bg_mod, grid_mod):
  !   - cells -1..n_r+2 and faces 0..n_r of the 1D grid; the ghost cells are
  !     the 1D grid's own neighbours (below the CZ base, and the atmosphere
  !     cells above the join), so radial stencils are those of CZTR;
  !   - rho0, T0, g from the file; s0 = s(rho0, T0) from the EOS table, in the
  !     bottom ghosts extended linearly from cells 1, 2; p0 = p(rho0, s0);
  !   - ds0/dr and the 4th-order radial metric dr/dxi of CZTR's grid_mod;
  !   - the initial state rho, T (-> s1 = s(rho, T) - s0) and v_r.
  ! Every rank reads the files. Also holds the EOS table.
  use kinds_mod, only: dp
  use eos_table, only: eos_table_t, eos_table_read, eos_eval, eos_s_from_rho_T, n_field, F_LNP, F_LNT, F_GAMMA1
  use eos_constants, only: k_B, m_H
  implicit none
  private

  type(eos_table_t), public :: tab
  integer, public :: n_r = 0, n_full = 0, k_join = 0
  ! The radial arrays are pointers: after bg_load they hold the whole column
  ! (cells -1:n_r+2, faces 0:n_r); bg_view(o, m) points them at the m cells
  ! o+1..o+m of one radial block (with its ghosts), indexed 1..m, so that
  ! block code reads them by its local cell index. n_c: the cells of the
  ! current view (n_r for the whole column).
  integer, public :: n_c = 0
  real(dp), pointer, public :: r_c(:) => null(), drdxi_c(:) => null()                  ! cells -1:n_r+2
  real(dp), pointer, public :: r_f(:) => null(), drdxi_f(:) => null(), area_f(:) => null()   ! faces 0:n_r (area r_f^2)
  real(dp), pointer, public :: vol_c(:) => null()                                      ! cells 1:n_r, volume per steradian
  real(dp), pointer, public :: rho0(:) => null(), s0(:) => null(), T0(:) => null(), p0(:) => null(), &
      lnp0(:) => null(), g_c(:) => null(), z_c(:) => null()                            ! cells -1:n_r+2
  real(dp), pointer, public :: g_f(:) => null()                                        ! faces 0:n_r
  real(dp), pointer, public :: ds0dr(:) => null()                                      ! cells 1:n_r
  real(dp), pointer, public :: cs0(:) => null(), xi(:) => null()   ! background sound speed, RSST factor; cells -1:n_r+2
  real(dp), pointer, public :: rho1_init(:) => null(), s1_init(:) => null(), v_init(:) => null()   ! cells -1:n_r+2
  ! the whole column (the views' targets)
  real(dp), pointer :: r_c_g(:), drdxi_c_g(:), r_f_g(:), drdxi_f_g(:), area_f_g(:), vol_c_g(:), rho0_g(:), &
      s0_g(:), T0_g(:), p0_g(:), lnp0_g(:), g_c_g(:), z_c_g(:), g_f_g(:), ds0dr_g(:), cs0_g(:), xi_g(:), &
      rho1_init_g(:), s1_init_g(:), v_init_g(:)
  real(dp), public :: r_top_cz = 0.0_dp

  public :: bg_load, set_rsst, ddr_c, ddr_f, face_value, bg_view

contains

  ! n_r_in: radial cells of the 3D run (0: the join cell k of the state file)
  ! r_top_rsun > 0: n_r = the cells with their centre below r_top_rsun x R_sun
  ! g_balanced: g = -(dp0/dr)/rho0 in cells 1..n_r (CZTR's ddr_c), so the
  ! background is exactly hydrostatic in the discrete momentum equation
  subroutine bg_load(grid_file, init_file, table_file, n_r_in, r_top_rsun, g_balanced)
    use const_mod, only: Rsun
    character(*), intent(in) :: grid_file, init_file, table_file
    integer, intent(in) :: n_r_in
    real(dp), intent(in) :: r_top_rsun
    logical, intent(in), optional :: g_balanced
    real(dp), allocatable :: rc_all(:), rf_all(:), row(:,:)
    real(dp) :: w, out(n_field)
    integer :: ierr, j, i

    call eos_table_read(tab, table_file, ierr)
    if (ierr /= 0) then
      print '(A,A,A,I0,A)', 'cz3d: cannot read the EOS table ', trim(table_file), ' (error ', ierr, &
          '); generate it with `make data` (eos_hhe''s eos_table_gen and data/eos_table.nml)'
      stop 1
    end if

    call read_grid(grid_file, n_full, rc_all, rf_all)
    call read_state(init_file, n_full, k_join, row)
    n_r = n_r_in
    if (n_r <= 0) n_r = k_join
    if (r_top_rsun > 0.0_dp) n_r = count(rc_all(1:k_join) <= r_top_rsun*Rsun)
    if (n_r + 2 > n_full) stop 'cz3d: n_r too large for the grid'
    n_c = n_r
    do j = -1, n_full + 2
      if (abs(row(1, j) - rc_all(j)) > 1.0e-6_dp*(rf_all(1) - rf_all(0))) stop 'cz3d: state and grid cells differ'
    end do

    allocate(r_c(-1:n_r+2), drdxi_c(-1:n_r+2), r_f(0:n_r), drdxi_f(0:n_r), area_f(0:n_r), vol_c(1:n_r))
    allocate(rho0(-1:n_r+2), s0(-1:n_r+2), T0(-1:n_r+2), p0(-1:n_r+2), lnp0(-1:n_r+2), g_c(-1:n_r+2), &
        z_c(-1:n_r+2), g_f(0:n_r), ds0dr(1:n_r), rho1_init(-1:n_r+2), s1_init(-1:n_r+2), v_init(-1:n_r+2), &
        cs0(-1:n_r+2), xi(-1:n_r+2))
    r_c = rc_all(-1:n_r+2)
    r_f = rf_all(0:n_r)
    ! metric as CZTR's grid_mod (4th order, exact for f linear in r)
    drdxi_c = 0.0_dp
    do i = 1, n_r
      drdxi_c(i) = (r_c(i-2) - 8.0_dp*r_c(i-1) + 8.0_dp*r_c(i+1) - r_c(i+2))/12.0_dp
    end do
    do i = 0, n_r
      drdxi_f(i) = (r_c(i-1) - 27.0_dp*r_c(i) + 27.0_dp*r_c(i+1) - r_c(i+2))/24.0_dp
    end do
    area_f = r_f**2
    do i = 1, n_r
      vol_c(i) = (r_f(i)**3 - r_f(i-1)**3)/3.0_dp
    end do

    ! background (as cztr_bg_mod)
    z_c = row(2, -1:n_r+2)*1.0e5_dp
    rho0 = row(6, -1:n_r+2)
    T0 = row(7, -1:n_r+2)
    g_c = row(8, -1:n_r+2)
    r_top_cz = r_c(1) - z_c(1)
    do j = -1, n_r + 2
      if (j == -1) then
        s0(j) = eos_s_from_rho_T(tab, rho0(j), T0(j), s_guess=30.0_dp*k_B/m_H)
      else
        s0(j) = eos_s_from_rho_T(tab, rho0(j), T0(j), s_guess=s0(j-1))
      end if
    end do
    do j = -1, 0
      s0(j) = s0(1) + (r_c(j) - r_c(1))*(s0(2) - s0(1))/(r_c(2) - r_c(1))
    end do
    do j = -1, n_r + 2
      call eos_eval(tab, rho0(j), s0(j), out)
      lnp0(j) = out(F_LNP)
      p0(j) = exp(lnp0(j))
      T0(j) = exp(out(F_LNT))
      cs0(j) = sqrt(out(F_GAMMA1)*p0(j)/rho0(j))
    end do
    if (present(g_balanced)) then
      if (g_balanced) g_c(1:n_r) = -ddr_c(p0)/rho0(1:n_r)
    end if
    do j = 0, n_r
      w = (r_f(j) - r_c(j))/(r_c(j+1) - r_c(j))
      g_f(j) = exp((1.0_dp - w)*log(g_c(j)) + w*log(g_c(j+1)))
    end do
    ds0dr = ddr_c(s0)
    xi = 1.0_dp

    ! initial state (cells 1..n_r; ghosts are set by the boundary conditions)
    rho1_init = 0.0_dp
    s1_init = 0.0_dp
    v_init = 0.0_dp
    do j = 1, n_r
      rho1_init(j) = row(3, j) - rho0(j)
      s1_init(j) = eos_s_from_rho_T(tab, row(3, j), row(4, j), s_guess=s0(j)) - s0(j)
      v_init(j) = row(5, j)
    end do
  end subroutine bg_load

  ! RSST factor xi = c_s0/c_red, c_red = min(c_s0, c_target (c_s0/c_s0(n_r))^q)
  ! (c_target <= 0: xi = 1; q = 0: xi = max(1, c_s0/c_target), as CZTR)
  ! Reduced sound speed c_red = min(c_s0, c_ref (c_s0/c_scale)^q), xi = c_s0/c_red.
  ! depth_km <= 0: c_ref = c_target, c_scale = c_s0 of the top cell (Hotta
  ! 2014 for q = 0, uniform xi for q = 1). depth_km > 0: c_ref = c_scale =
  ! c_s0 at that depth below the photosphere (c_target unused), so xi = 1
  ! wherever c_s0 <= c_ref (above that depth) and (c_s0/c_ref)^(1-q) below.
  ! Never any reduction above the photosphere (z > 0).
  subroutine set_rsst(c_target, q, depth_km)
    real(dp), intent(in) :: c_target
    real(dp), intent(in), optional :: q, depth_km
    real(dp) :: qq, c_ref, c_scale, zd, w
    integer :: k
    xi = 1.0_dp
    qq = 0.0_dp
    if (present(q)) qq = q
    c_ref = c_target
    c_scale = cs0(n_r)
    if (present(depth_km)) then
      if (depth_km > 0.0_dp) then
        zd = -depth_km*1.0e5_dp
        if (zd < z_c(1) .or. zd > z_c(n_r)) stop 'set_rsst: rsst_depth_km outside the domain'
        do k = 1, n_r - 1
          if (z_c(k+1) >= zd) exit
        end do
        w = (zd - z_c(k))/(z_c(k+1) - z_c(k))
        c_ref = (1.0_dp - w)*cs0(k) + w*cs0(k+1)
        c_scale = c_ref
      end if
    end if
    if (c_ref <= 0.0_dp) return
    xi = cs0/min(cs0, c_ref*(cs0/c_scale)**qq)
    where (z_c > 0.0_dp) xi = 1.0_dp
    end subroutine set_rsst

  ! CZTR grid file: '#' lines, n, then 'i r_c' for i = -1..n+2, 'i r_f' for i = 0..n
  subroutine read_grid(path, n, rc, rf)
    character(*), intent(in) :: path
    integer, intent(out) :: n
    real(dp), allocatable, intent(out) :: rc(:), rf(:)
    character(len=256) :: line
    integer :: u, ios, i, j
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'cz3d: cannot open the grid file'
    do
      read(u, '(A)') line
      if (line(1:1) /= '#') exit
    end do
    read(line, *) n
    allocate(rc(-1:n+2), rf(0:n))
    do i = -1, n + 2
      read(u, *) j, rc(i)
    end do
    do i = 0, n
      read(u, *) j, rf(i)
    end do
    close(u)
  end subroutine read_grid

  ! CZTR state file (init_full.dat format): n, k, then rows
  ! 'i r z[km] rho T v rho0 T0 g' for i = -1..n+2; row(1:8, i)
  subroutine read_state(path, n, k, row)
    character(*), intent(in) :: path
    integer, intent(in) :: n
    integer, intent(out) :: k
    real(dp), allocatable, intent(out) :: row(:,:)
    character(len=512) :: line
    integer :: u, ios, i, j, nn
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'cz3d: cannot open the state file'
    do
      read(u, '(A)') line
      if (line(1:1) /= '#') exit
    end do
    read(line, *) nn, k
    if (nn /= n) stop 'cz3d: state file and grid have different cell counts'
    allocate(row(8, -1:n+2))
    do j = -1, n + 2
      read(u, *) i, row(:, j)
      if (i /= j) stop 'cz3d: state file: unexpected cell index'
    end do
    close(u)
  end subroutine read_state

  ! Radial derivative at cells 1..n_r of f(-1:n_r+2) (CZTR's ddr_c)
  pure function ddr_c(f) result(df)
    real(dp), intent(in) :: f(-1:)
    real(dp) :: df(1:n_c)
    integer :: i
    do i = 1, n_c
      df(i) = (f(i-2) - 8.0_dp*f(i-1) + 8.0_dp*f(i+1) - f(i+2))/(12.0_dp*drdxi_c(i))
    end do
  end function ddr_c

  ! Radial derivative at faces 0..n_r (CZTR's ddr_f)
  pure function ddr_f(f) result(df)
    real(dp), intent(in) :: f(-1:)
    real(dp) :: df(0:n_c)
    integer :: i
    do i = 0, n_c
      df(i) = (f(i-1) - 27.0_dp*f(i) + 27.0_dp*f(i+1) - f(i+2))/(24.0_dp*drdxi_f(i))
    end do
  end function ddr_f

  ! 4th-order face values at faces 0..n_r (CZTR's face_value)
  pure function face_value(f) result(ff)
    real(dp), intent(in) :: f(-1:)
    real(dp) :: ff(0:n_c)
    integer :: i
    do i = 0, n_c
      ff(i) = (-f(i-1) + 9.0_dp*f(i) + 9.0_dp*f(i+1) - f(i+2))/16.0_dp
    end do
  end function face_value

  ! Point the radial arrays at the m cells o+1..o+m (one radial block, its
  ! ghosts -1, 0, m+1, m+2 the neighbouring cells), indexed 1..m; bg_view(0,
  ! n_r): the whole column. The first call stores the whole column.
  subroutine bg_view(o, m)
    integer, intent(in) :: o, m
    if (.not. associated(r_c_g)) then
      r_c_g => r_c; drdxi_c_g => drdxi_c; r_f_g => r_f; drdxi_f_g => drdxi_f; area_f_g => area_f
      vol_c_g => vol_c; rho0_g => rho0; s0_g => s0; T0_g => T0; p0_g => p0; lnp0_g => lnp0; g_c_g => g_c
      z_c_g => z_c; g_f_g => g_f; ds0dr_g => ds0dr; cs0_g => cs0; xi_g => xi
      rho1_init_g => rho1_init; s1_init_g => s1_init; v_init_g => v_init
    end if
    r_c(-1:m+2) => r_c_g(o-1:o+m+2);  drdxi_c(-1:m+2) => drdxi_c_g(o-1:o+m+2)
    rho0(-1:m+2) => rho0_g(o-1:o+m+2);  s0(-1:m+2) => s0_g(o-1:o+m+2);  T0(-1:m+2) => T0_g(o-1:o+m+2)
    p0(-1:m+2) => p0_g(o-1:o+m+2);  lnp0(-1:m+2) => lnp0_g(o-1:o+m+2);  g_c(-1:m+2) => g_c_g(o-1:o+m+2)
    z_c(-1:m+2) => z_c_g(o-1:o+m+2);  cs0(-1:m+2) => cs0_g(o-1:o+m+2);  xi(-1:m+2) => xi_g(o-1:o+m+2)
    rho1_init(-1:m+2) => rho1_init_g(o-1:o+m+2);  s1_init(-1:m+2) => s1_init_g(o-1:o+m+2)
    v_init(-1:m+2) => v_init_g(o-1:o+m+2)
    r_f(0:m) => r_f_g(o:o+m);  drdxi_f(0:m) => drdxi_f_g(o:o+m);  area_f(0:m) => area_f_g(o:o+m)
    g_f(0:m) => g_f_g(o:o+m)
    vol_c(1:m) => vol_c_g(o+1:o+m);  ds0dr(1:m) => ds0dr_g(o+1:o+m)
    n_c = m
  end subroutine bg_view

end module cz3d_bg_mod
