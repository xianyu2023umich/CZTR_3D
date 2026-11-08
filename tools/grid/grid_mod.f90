module grid_mod
  ! CZ_3D's copy of CZTR's src/cztr_grid_mod.f90 (2026-10-22), with one more
  ! density term: n_h = 1/(aspect r dtheta) for a horizontal spacing dtheta
  ! (horiz_deg > 0), so that dr <= aspect * r dtheta everywhere; and the base
  ! pressure p_base used by tools/grid_gen.f90 to cut the profile.
  !
  ! One grid for the whole CZ + atmosphere domain, from the CZ base to z_top
  ! above the CZ top (tau = 0.01, z = 0). The grid is defined by a target cell
  ! spacing dr(r) built from the cell densities (cells per unit length) that
  ! each region asks for, added together:
  !
  !   n(r) = n_CZ + n_atm + n_TR
  !   n_CZ  = cells_per_Hp / H_p(r)            (CZ_1D's ln P grid: k cells per H_p)
  !   n_atm = S(z) / (atm_dr0 (1 + z+/atm_za))  (TR's grid: spacing growing
  !                                              linearly with height)
  !   n_TR  = (1/tr_band_dr) B(z)               (optional band for the TR;
  !                                              off by default: with it 30 km
  !                                              cells, 421 cells, dt 0.15 s)
  !
  ! CZ_3D: with horiz_deg > 0, n is then raised to at least
  !
  !   n_h = 1 / (aspect r dtheta),  dtheta = horiz_deg
  !
  ! by a smooth maximum, n -> (n^8 + n_h^8)^(1/8), so that radial cells are no
  ! taller than aspect x the horizontal spacing (and are unchanged where they
  ! already are);
  ! and the smallest spacing capped smoothly at dr_min:
  !
  !   dr = (dr_min^p + n^-p)^(1/p)
  !
  ! S(z) = (1 + tanh(z/w_s))/2 switches the atmosphere term on above the CZ
  ! top, z+ = w_s ln(1 + exp(z/w_s)) is a smooth max(z, 0), and B(z) is a
  ! smooth top-hat between tr_band_lo and tr_band_hi. The cell index is
  ! xi(r) = integral of dr'/dr(r'), rescaled so the domain holds an integer
  ! number of cells; faces sit at integer xi, cells at half-integers.
  ! Everything is smooth (dr is C-infinity in r apart from the tabulated H_p),
  ! so the map r(xi) is smooth: no seam between the CZ and atmosphere parts.
  !
  ! H_p: below the CZ top from the CZ profile; above it from a simple model
  ! atmosphere (T and mu stepping by tanh from chromospheric to coronal values
  ! at z_TRH), H = k_B T/(mu m_H g), scaled to equal the profile's H_p at the
  ! CZ top. ln H_p is then smoothed with a Gaussian of width smooth_km, which
  ! removes the derivative mismatch at the join. The grid therefore does not
  ! depend on the (still to be designed) atmosphere background.
  !
  ! Defaults reproduce the two separate grids in their own regions: the CZ
  ! with 5 cells per H_p (CZ_1D grid_type 'lnp'), the corona with TR's law
  ! (N = 256: 52 km at its base, spacing proportional to z + 0.01 Rsun), and
  ! no spacing below ~30 km: 394 cells, CFL step 0.34 s (2026-10-02).
  use kinds_mod, only: dp
  implicit none
  private

  real(dp), parameter :: km = 1.0e5_dp
  real(dp), parameter :: G_grav = 6.67430e-8_dp, Msun = 1.98892e33_dp
  real(dp), parameter :: k_B = 1.380649e-16_dp, m_H = 1.6735575e-24_dp

  type, public :: grid_params
    real(dp) :: cells_per_Hp  = 5.0_dp       ! CZ / chromosphere: cells per pressure scale height
    real(dp) :: dr_min_km     = 30.0_dp      ! smooth floor of the spacing
    integer  :: p_cap         = 4            ! sharpness of the floor
    real(dp) :: atm_dr0_km    = 52.3_dp      ! atmosphere law: spacing at z = 0 (TR, N = 256)
    real(dp) :: atm_za_km     = 6957.0_dp    ! ... grows as (1 + z/atm_za) (TR: 0.01 Rsun)
    real(dp) :: atm_switch_km = 200.0_dp     ! width of the switch-on of the atmosphere term
    logical  :: use_tr_band   = .false.      ! off: the TR gets ~40-75 km cells (as TR N=256), dt 0.34 s
    real(dp) :: tr_band_lo_km = 1400.0_dp    ! band where the TR may sit
    real(dp) :: tr_band_hi_km = 3000.0_dp
    real(dp) :: tr_band_w_km  = 200.0_dp     ! edge width of the band
    real(dp) :: tr_band_dr_km = 30.0_dp      ! spacing asked for inside the band
    real(dp) :: z_top_km      = 40000.0_dp   ! domain top above the CZ top
    ! optional refinement band below the photosphere (the convection-to-
    ! radiation handover layer): spacing cz_band_dr_km between cz_band_lo_km
    ! and cz_band_hi_km (heights z, CZ top = 0), the floor dr_min lowered there
    logical  :: use_cz_band   = .false.
    real(dp) :: cz_band_lo_km = -700.0_dp
    real(dp) :: cz_band_hi_km = 50.0_dp
    real(dp) :: cz_band_w_km  = 100.0_dp     ! edge width of the band
    real(dp) :: cz_band_dr_km = 12.0_dp
    ! model atmosphere for H_p above the CZ top
    real(dp) :: T_chrom = 4871.75_dp, T_cor = 5.0e5_dp
    real(dp) :: mu_chrom = 1.25_dp, mu_cor = 0.6_dp
    real(dp) :: z_TRH_km = 2000.0_dp, w_TRH_km = 300.0_dp
    real(dp) :: smooth_km     = 50.0_dp      ! Gaussian smoothing of ln H_p
    real(dp) :: table_dr_km   = 1.0_dp       ! resolution of the internal tables
    ! CZ_3D: radial cells at most aspect x the horizontal spacing r dtheta
    real(dp) :: horiz_deg     = 0.0_dp       ! horizontal spacing dtheta [deg]; 0: term off
    real(dp) :: aspect        = 2.0_dp
    integer  :: n_cz          = 0            ! CZ_3D: cells below the CZ top (0: as the density gives)
    real(dp) :: p_base        = 6.0e13_dp    ! domain base: the profile is cut at this pressure (grid_gen)
  end type grid_params

  ! Internal tables, kept after grid_build for grid_target.
  type(grid_params) :: par_s
  real(dp) :: r_top_s = 0.0_dp, hscale_s = 1.0_dp
  real(dp), allocatable :: tab_r(:), tab_lnH(:)

  public :: grid_build, grid_target, atm_T_mu, grid_read_nml

contains

  ! Build the grid. r_prof(:), Hp_prof(:): the CZ profile ordered by
  ! increasing r, its last point being the CZ top. Returns n cells, cell
  ! centres r_c(-1:n+2) (two ghosts each side) and faces r_f(0:n).
  subroutine grid_build(par, r_prof, Hp_prof, n, r_c, r_f)
    type(grid_params), intent(in) :: par
    real(dp), intent(in) :: r_prof(:), Hp_prof(:)
    integer, intent(out) :: n
    real(dp), allocatable, intent(out) :: r_c(:), r_f(:)
    real(dp), allocatable :: raw(:), dr(:), xi(:), slope(:)
    real(dp) :: r_bot, r_end, h, sig, w, ws, scale, s0
    integer :: m, i, j, kh

    par_s = par
    r_bot = r_prof(1)
    r_top_s = r_prof(size(r_prof))
    r_end = r_top_s + par%z_top_km*km
    hscale_s = 1.0_dp
    hscale_s = Hp_prof(size(Hp_prof))/H_model(0.0_dp)

    m = int((r_end - r_bot)/(par%table_dr_km*km)) + 1
    h = (r_end - r_bot)/real(m - 1, dp)
    if (allocated(tab_r)) deallocate(tab_r, tab_lnH)
    allocate(tab_r(m), tab_lnH(m), raw(m), dr(m), xi(m), slope(m))
    do i = 1, m
      tab_r(i) = r_bot + real(i - 1, dp)*h
      raw(i) = lnH_raw(tab_r(i), r_prof, Hp_prof)
    end do

    ! Gaussian smoothing of ln H_p (kernel renormalized at the domain ends).
    sig = par%smooth_km*km
    kh = max(1, nint(4.0_dp*sig/h))
    do i = 1, m
      s0 = 0.0_dp
      ws = 0.0_dp
      do j = max(1, i - kh), min(m, i + kh)
        w = exp(-0.5_dp*((real(j - i, dp)*h)/sig)**2)
        s0 = s0 + w*raw(j)
        ws = ws + w
      end do
      tab_lnH(i) = s0/ws
    end do

    ! Target spacing and the cell index xi(r) = int dr/dr_target.
    do i = 1, m
      dr(i) = target_spacing(tab_r(i), tab_lnH(i))
    end do
    xi(1) = 0.0_dp
    do i = 2, m
      xi(i) = xi(i-1) + 0.5_dp*h*(1.0_dp/dr(i-1) + 1.0_dp/dr(i))
    end do
    n = max(8, nint(xi(m)))
    scale = xi(m)/real(n, dp)               ! actual spacing = target * scale (~1)
    if (par%n_cz > 0) then
      ! CZ_3D: exactly n_cz cells below the CZ top (the CZ top on face n_cz); the
      ! atmosphere above takes the count the same scaled density gives
      scale = interp_xi(r_top_s)/real(par%n_cz, dp)
      n = int(xi(m)/scale)
    end if
    xi = xi/scale
    slope = dr*scale                        ! dr/dxi at the table nodes

    allocate(r_f(0:n), r_c(-1:n+2))
    do i = 0, n
      r_f(i) = r_of_xi(real(i, dp))
    end do
    do i = -1, n + 2
      r_c(i) = r_of_xi(real(i, dp) - 0.5_dp)
    end do

  contains

    ! xi(r) on the table (linear interpolation), unscaled
    real(dp) function interp_xi(rr) result(x)
      real(dp), intent(in) :: rr
      integer :: ii
      real(dp) :: w
      ii = min(m - 1, max(1, int((rr - tab_r(1))/h) + 1))
      w = (rr - tab_r(ii))/h
      x = xi(ii) + w*(xi(ii+1) - xi(ii))
    end function interp_xi

    ! Inverse of xi(r): cubic Hermite on the table (exact slopes), linear
    ! extension beyond the ends (ghost cells).
    real(dp) function r_of_xi(x) result(r)
      real(dp), intent(in) :: x
      integer :: lo, hi, mid
      real(dp) :: dx, s, h00, h10, h01, h11
      if (x <= 0.0_dp) then
        r = tab_r(1) + x*slope(1)
        return
      end if
      if (x >= xi(m)) then
        r = tab_r(m) + (x - xi(m))*slope(m)
        return
      end if
      lo = 1
      hi = m
      do while (hi - lo > 1)
        mid = (lo + hi)/2
        if (xi(mid) <= x) then
          lo = mid
        else
          hi = mid
        end if
      end do
      dx = xi(lo+1) - xi(lo)
      s = (x - xi(lo))/dx
      h00 = (1.0_dp + 2.0_dp*s)*(1.0_dp - s)**2
      h10 = s*(1.0_dp - s)**2
      h01 = s**2*(3.0_dp - 2.0_dp*s)
      h11 = s**2*(s - 1.0_dp)
      r = h00*tab_r(lo) + h10*dx*slope(lo) + h01*tab_r(lo+1) + h11*dx*slope(lo+1)
    end function r_of_xi

  end subroutine grid_build

  ! Target spacing at r and its ingredients (each as a spacing 1/n_i; huge if
  ! that term is off), plus the smoothed H_p used. Valid after grid_build.
  subroutine grid_target(r, dr, dr_cz, dr_atm, dr_tr, Hp)
    real(dp), intent(in) :: r
    real(dp), intent(out) :: dr, dr_cz, dr_atm, dr_tr, Hp
    real(dp) :: lnH, ncz, natm, ntr
    lnH = interp_lnH(r)
    Hp = exp(lnH)
    call densities(r, lnH, ncz, natm, ntr)
    dr = cap(with_horiz(ncz + natm + ntr, r), r - r_top_s)
    dr_cz = 1.0_dp/max(ncz, 1.0e-300_dp)
    dr_atm = 1.0_dp/max(natm, 1.0e-300_dp)
    dr_tr = 1.0_dp/max(ntr, 1.0e-300_dp)
  end subroutine grid_target

  ! Model atmosphere temperature and mean molecular weight at height z [cm].
  pure subroutine atm_T_mu(par, z, T, mu)
    type(grid_params), intent(in) :: par
    real(dp), intent(in) :: z
    real(dp), intent(out) :: T, mu
    real(dp) :: st
    st = 0.5_dp*(1.0_dp + tanh((z - par%z_TRH_km*km)/(par%w_TRH_km*km)))
    T = par%T_chrom + (par%T_cor - par%T_chrom)*st
    mu = par%mu_chrom + (par%mu_cor - par%mu_chrom)*st
  end subroutine atm_T_mu

  ! Read grid parameters from a namelist &grid in path (missing items keep
  ! their defaults). Used by tools/grid_gen.f90.
  subroutine grid_read_nml(path, par)
    character(*), intent(in) :: path
    type(grid_params), intent(inout) :: par
    real(dp) :: cells_per_Hp, dr_min_km, atm_dr0_km, atm_za_km, tr_band_lo_km, tr_band_hi_km, tr_band_dr_km, &
        z_top_km, cz_band_lo_km, cz_band_hi_km, cz_band_w_km, cz_band_dr_km, horiz_deg, aspect, p_base
    integer :: n_cz
    logical :: use_tr_band, use_cz_band
    integer :: u, ios
    namelist /grid/ cells_per_Hp, dr_min_km, atm_dr0_km, atm_za_km, use_tr_band, tr_band_lo_km, tr_band_hi_km, &
        tr_band_dr_km, z_top_km, use_cz_band, cz_band_lo_km, cz_band_hi_km, cz_band_w_km, cz_band_dr_km, &
        horiz_deg, aspect, p_base, n_cz
    cells_per_Hp = par%cells_per_Hp;  dr_min_km = par%dr_min_km;  atm_dr0_km = par%atm_dr0_km
    atm_za_km = par%atm_za_km;  use_tr_band = par%use_tr_band;  tr_band_lo_km = par%tr_band_lo_km
    tr_band_hi_km = par%tr_band_hi_km;  tr_band_dr_km = par%tr_band_dr_km;  z_top_km = par%z_top_km
    use_cz_band = par%use_cz_band;  cz_band_lo_km = par%cz_band_lo_km;  cz_band_hi_km = par%cz_band_hi_km
    cz_band_w_km = par%cz_band_w_km;  cz_band_dr_km = par%cz_band_dr_km
    horiz_deg = par%horiz_deg;  aspect = par%aspect;  p_base = par%p_base;  n_cz = par%n_cz
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print '(2A)', 'Cannot open grid namelist ', trim(path)
      stop 1
    end if
    read(u, nml=grid)
    close(u)
    par%cells_per_Hp = cells_per_Hp;  par%dr_min_km = dr_min_km;  par%atm_dr0_km = atm_dr0_km
    par%atm_za_km = atm_za_km;  par%use_tr_band = use_tr_band;  par%tr_band_lo_km = tr_band_lo_km
    par%tr_band_hi_km = tr_band_hi_km;  par%tr_band_dr_km = tr_band_dr_km;  par%z_top_km = z_top_km
    par%use_cz_band = use_cz_band;  par%cz_band_lo_km = cz_band_lo_km;  par%cz_band_hi_km = cz_band_hi_km
    par%cz_band_w_km = cz_band_w_km;  par%cz_band_dr_km = cz_band_dr_km
    par%horiz_deg = horiz_deg;  par%aspect = aspect;  par%p_base = p_base;  par%n_cz = n_cz
  end subroutine grid_read_nml

  ! ---- internals ----

  real(dp) function target_spacing(r, lnH) result(dr)
    real(dp), intent(in) :: r, lnH
    real(dp) :: ncz, natm, ntr
    call densities(r, lnH, ncz, natm, ntr)
    dr = cap(with_horiz(ncz + natm + ntr, r), r - r_top_s)
  end function target_spacing

  subroutine densities(r, lnH, ncz, natm, ntr)
    real(dp), intent(in) :: r, lnH
    real(dp), intent(out) :: ncz, natm, ntr
    real(dp) :: z, ws, zp, wb, bb
    z = r - r_top_s
    ncz = par_s%cells_per_Hp/exp(lnH)
    ws = par_s%atm_switch_km*km
    zp = ws*log1p_exp(z/ws)
    natm = 0.5_dp*(1.0_dp + tanh(z/ws))/(par_s%atm_dr0_km*km*(1.0_dp + zp/(par_s%atm_za_km*km)))
    ntr = 0.0_dp
    if (par_s%use_tr_band) then
      wb = par_s%tr_band_w_km*km
      ntr = 0.5_dp*(tanh((z - par_s%tr_band_lo_km*km)/wb) - tanh((z - par_s%tr_band_hi_km*km)/wb)) &
          /(par_s%tr_band_dr_km*km)
    end if
    if (par_s%use_cz_band) then
      bb = cz_band_weight(z)
      ntr = ntr + bb/(par_s%cz_band_dr_km*km)
    end if
  end subroutine densities

  ! CZ_3D: the total density with the horizontal-spacing constraint, a smooth
  ! maximum (power mean, p = 8) of the sum of the other terms and
  ! n_h = 1/(aspect r dtheta): n_h only raises the density where the other
  ! terms ask for cells taller than aspect x r dtheta, and adds nothing where
  ! they are already finer (a plain sum would refine both regions).
  real(dp) function with_horiz(n, r) result(nt)
    real(dp), intent(in) :: n, r
    real(dp) :: nh
    nt = n
    if (par_s%horiz_deg <= 0.0_dp) return
    nh = 1.0_dp/(par_s%aspect*r*par_s%horiz_deg*acos(-1.0_dp)/180.0_dp)
    nt = n*(1.0_dp + (nh/n)**8)**0.125_dp
  end function with_horiz

  ! Smooth top-hat of the CZ-top refinement band at height z [cm].
  real(dp) function cz_band_weight(z) result(bb)
    real(dp), intent(in) :: z
    real(dp) :: w
    w = par_s%cz_band_w_km*km
    bb = 0.5_dp*(tanh((z - par_s%cz_band_lo_km*km)/w) - tanh((z - par_s%cz_band_hi_km*km)/w))
  end function cz_band_weight

  ! Smooth floor on the spacing 1/n: dr_min, lowered to cz_band_dr inside the
  ! CZ-top band (z = height above the CZ top).
  real(dp) function cap(ntot, z) result(dr)
    real(dp), intent(in) :: ntot, z
    real(dp) :: d0
    integer :: p
    p = par_s%p_cap
    d0 = par_s%dr_min_km*km
    if (par_s%use_cz_band) d0 = d0 - (par_s%dr_min_km - par_s%cz_band_dr_km)*km*cz_band_weight(z)
    dr = (d0**p + (1.0_dp/ntot)**p)**(1.0_dp/real(p, dp))
  end function cap

  ! ln(1 + e^x) without overflow.
  pure real(dp) function log1p_exp(x)
    real(dp), intent(in) :: x
    if (x > 30.0_dp) then
      log1p_exp = x
    else
      log1p_exp = log(1.0_dp + exp(x))
    end if
  end function log1p_exp

  ! Model-atmosphere H_p at height z (unscaled).
  real(dp) function H_model(z)
    real(dp), intent(in) :: z
    real(dp) :: T, mu, g
    call atm_T_mu(par_s, z, T, mu)
    g = G_grav*Msun/(r_top_s + z)**2
    H_model = k_B*T/(mu*m_H*g)
  end function H_model

  ! Unsmoothed ln H_p: the CZ profile below the top, the scaled model above.
  real(dp) function lnH_raw(r, r_prof, Hp_prof)
    real(dp), intent(in) :: r, r_prof(:), Hp_prof(:)
    integer :: lo, hi, mid, m
    real(dp) :: w
    if (r > r_top_s) then
      lnH_raw = log(hscale_s*H_model(r - r_top_s))
      return
    end if
    m = size(r_prof)
    lo = 1
    hi = m
    do while (hi - lo > 1)
      mid = (lo + hi)/2
      if (r_prof(mid) <= r) then
        lo = mid
      else
        hi = mid
      end if
    end do
    w = (r - r_prof(lo))/(r_prof(lo+1) - r_prof(lo))
    lnH_raw = log(Hp_prof(lo)) + w*(log(Hp_prof(lo+1)) - log(Hp_prof(lo)))
  end function lnH_raw

  real(dp) function interp_lnH(r)
    real(dp), intent(in) :: r
    integer :: i, m
    real(dp) :: w
    m = size(tab_r)
    i = int((r - tab_r(1))/(tab_r(2) - tab_r(1))) + 1
    i = max(1, min(m - 1, i))
    w = (r - tab_r(i))/(tab_r(i+1) - tab_r(i))
    interp_lnH = tab_lnH(i) + w*(tab_lnH(i+1) - tab_lnH(i))
  end function interp_lnH

end module grid_mod
