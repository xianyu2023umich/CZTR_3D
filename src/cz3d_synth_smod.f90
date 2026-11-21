submodule (cz3d_app_mod) cz3d_synth_smod
  ! synthetic observables: EUV (optically thin) and continuum (1.5D LTE) images
  ! (part of cz3d_app_mod: the module's data and routines by host association)
  implicit none

contains

  ! Box output (box_radius_deg > 0; see cz3d_params_mod): the cells near a
  ! point (e.g. where a run failed), full state, one file per rank
  ! Synthetic EUV images (euv_dt, euv_at_end; see cz3d_params_mod). The
  ! response table (euv_table, read once): '#' comment lines, then
  ! "n_T n_channel", the channel wavelengths [A], then n_T rows
  ! "log10 T  K_1 .. K_n" (uniform in log10 T; K [DN cm^5 s^-1 pix^-1] per
  ! n_e n_H). The emission n_e n_H K(T) is sampled at sub-points of each
  ! atmosphere cell, uniform in (r, theta, phi): per direction
  ! ceil(3 extent/pixel) points (at most 12 across, 32 radially), radially
  ! also ceil(dlogT/0.02)
  ! where log10 T changes by dlogT to a radial neighbour (at most 32), so the
  ! transition region is integrated through. At each sub-point log10 T and
  ! ln(n_e n_H) are trilinear between the cell centres (one ghost layer, the
  ! state's ghosts filled first), each carries the cell volume / its count;
  ! with euv_trac the emission times the TR broadening factor f(T) (trac_f).
  ! A sub-point is weighted by its Yin-Yang core cover (1/2 where both
  ! patches' cores hold it) and dropped when the sphere at the join hides it
  ! from the observer; it is deposited as a square of side its horizontal
  ! spacing, split over the pixels it overlaps (so a cell's points tile the
  ! image without the moire a point deposit gives). The images are summed over the ranks and written by the root:
  ! out_dir/euv_NNNNNNNN.bin, stream: int32 n_view, n_channel, n_pix; real64
  ! t, fov [R_sun], lon(n_view), lat(n_view) [deg]; int32 channels(n_channel);
  ! real32 image(n_pix, n_pix, n_channel, n_view) (x fastest: increasing
  ! longitude at the sub-observer point; y: north; pixel (1,1) centred at
  ! (-fov + pix/2, -fov + pix/2)).
  subroutine euv_coronal(iStep, t)
    use const_mod, only: Rsun
    use radiative_loss_mod, only: N_H_per_gram_global
    use eos_table, only: F_LNNE
    use mpi_f08, only: MPI_SUM, MPI_Reduce
    use iso_fortran_env, only: int32, real32
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    real(dp), parameter :: dlt_sub = 0.02_dp
    real(dp), allocatable :: img(:,:,:,:), em(:), lt(:,:,:), ln(:,:,:)
    real(dp) :: d(3, euv_max_view), ex(3, euv_max_view), ey(3, euv_max_view), lo, la, pix, rocc, fov, &
        out(n_field), rho, lgT, f, w, dvs, rc, dth, dph, ra, rb, rs, ths, phs, p(3), q(3), u, v, &
        sdot, rho2, th2, ph2, ext, t0, mx, fx, fy, wr, wt, wp, ltp, lnp, dlt, lt_lo, lt_hi, hp, wx(2), wy(2)
    real(dp) :: rsub(32), thsub(12), phsub(12)
    integer :: nch, nT, jb, ib, i, j, k, m, iv, ir, it, ip, nsr, nst, nsp, ix, iy, jx, jy, u_, n_img, &
        i0, j0, k0
    character(len=300) :: name
    character(len=200) :: msg
    t0 = arbor_wtime()
    if (.not. allocated(euv_K)) call euv_load()
    nch = size(euv_ch)
    nT = size(euv_lgT)
    lt_lo = euv_lgT(1);  lt_hi = euv_lgT(nT)
    allocate(img(euv_npix, euv_npix, nch, euv_n_view), em(nch))
    allocate(lt(0:nb(1)+1, 0:nb(2)+1, 0:nb(3)+1), ln(0:nb(1)+1, 0:nb(2)+1, 0:nb(3)+1))
    img = 0.0_dp
    do iv = 1, euv_n_view
      lo = euv_lon_deg(iv)*pi/180.0_dp;  la = euv_lat_deg(iv)*pi/180.0_dp
      d(:, iv) = [cos(la)*cos(lo), cos(la)*sin(lo), sin(la)]                ! towards the observer
      ex(:, iv) = [-sin(lo), cos(lo), 0.0_dp]                                ! increasing longitude
      ey(:, iv) = [-sin(la)*cos(lo), -sin(la)*sin(lo), cos(la)]              ! north
    end do
    fov = euv_fov_rsun*Rsun
    pix = 2.0_dp*fov/euv_npix
    call arbor_fill_ghosts(State_GVB, n_ghost=1, faces_only=.false.)
    call cz_unbind()
    rocc = r_f(k_join)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (r_off + nb(1) <= k_join) cycle
      associate (th => Metric_DB(2, ib)%X_I, dth_ => Metric_DB(2, ib)%Dx_I, ph => Metric_DB(3, ib)%X_I, &
                 dph_ => Metric_DB(3, ib)%Dx_I, vth => Metric_DB(2, ib)%Vol_I, vph => Metric_DB(3, ib)%Vol_I)
      if (lbound(th, 1) > 0 .or. ubound(th, 1) < nb(2) + 1 .or. lbound(ph, 1) > 0 .or. ubound(ph, 1) < nb(3) + 1) &
          call arbor_abort('euv_write: the block metric has no ghost centres')
      ! log10 T and ln(n_e n_H) at the cell centres and one ghost layer
      do k = 0, nb(3) + 1
        do j = 0, nb(2) + 1
          do i = 0, nb(1) + 1
            rho = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, max(rho, 1.0e-30_dp), s0(i) + State_GVB(i, j, k, s1_, ib), out)
            lt(i, j, k) = out(F_LNT)/log(10.0_dp)
            ln(i, j, k) = out(F_LNNE) + log(N_H_per_gram_global*max(rho, 1.0e-30_dp))
          end do
        end do
      end do
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            if (r_off + i <= k_join) cycle
            ! nothing in the table's range in the cell or its neighbours: no emission
            if (maxval(lt(i-1:i+1, j-1:j+1, k-1:k+1)) <= lt_lo .or. minval(lt(i-1:i+1, j-1:j+1, k-1:k+1)) >= lt_hi) cycle
            ra = r_f(i-1);  rb = r_f(i)
            rc = r_c(i);  dth = dth_(j);  dph = dph_(k)
            ! the centre's projection: skip the cell when it lies outside every image by more than its size
            p = rc*[sin(th(j))*cos(ph(k)), sin(th(j))*sin(ph(k)), cos(th(j))]
            if (Key_IB(KP, ib) == 2) p = yy_cart(p)
            ext = (rb - ra) + rc*(dth + dph)
            if (all([(max(abs(dot_product(p, ex(:, iv))), abs(dot_product(p, ey(:, iv)))) > fov + ext, &
                iv = 1, euv_n_view)])) cycle
            dlt = max(abs(lt(i+1, j, k) - lt(i, j, k)), abs(lt(i, j, k) - lt(i-1, j, k)))
            nsr = min(32, max(1, ceiling(3.0_dp*(rb - ra)/pix), ceiling(dlt/dlt_sub)))
            nst = min(12, max(1, ceiling(3.0_dp*rc*dth/pix)))
            nsp = min(12, max(1, ceiling(3.0_dp*rc*sin(th(j))*dph/pix)))
            do ir = 1, nsr
              rsub(ir) = ra + (ir - 0.5_dp)*(rb - ra)/nsr
            end do
            do it = 1, nst
              thsub(it) = th(j) + ((it - 0.5_dp)/nst - 0.5_dp)*dth
            end do
            do ip = 1, nsp
              phsub(ip) = ph(k) + ((ip - 0.5_dp)/nsp - 0.5_dp)*dph
            end do
            dvs = (rb**3 - ra**3)/3.0_dp*vth(j)*vph(k)/real(nsr*nst*nsp, dp)/pix**2
            hp = min(1.0_dp, rc*sqrt(dth*sin(th(j))*dph/real(nst*nsp, dp))/pix)   ! the horizontal spacing [pixels]
            do ip = 1, nsp
              phs = phsub(ip)
              k0 = merge(k - 1, k, phs < ph(k))
              wp = (phs - ph(k0))/(ph(k0+1) - ph(k0))
              do it = 1, nst
                ths = thsub(it)
                w = 1.0_dp
                if (.not. yy_core(ths, phs)) cycle
                call yy_angles(ths, phs, th2, ph2)
                if (yy_core(th2, ph2)) w = 0.5_dp
                j0 = merge(j - 1, j, ths < th(j))
                wt = (ths - th(j0))/(th(j0+1) - th(j0))
                q = [sin(ths)*cos(phs), sin(ths)*sin(phs), cos(ths)]
                if (Key_IB(KP, ib) == 2) q = yy_cart(q)
                do ir = 1, nsr
                  rs = rsub(ir)
                  i0 = merge(i - 1, i, rs < r_c(i))
                  wr = (rs - r_c(i0))/(r_c(i0+1) - r_c(i0))
                  ltp = tri(lt)
                  if (ltp <= lt_lo .or. ltp >= lt_hi) cycle
                  lnp = tri(ln)
                  f = (ltp - lt_lo)/(euv_lgT(2) - lt_lo)
                  m = min(nT - 1, int(f) + 1)
                  f = f - (m - 1)
                  em = w*exp(lnp)*((1.0_dp - f)*euv_K(m, :) + f*euv_K(m+1, :))*dvs
                  if (euv_trac) em = em*trac_f(10.0_dp**ltp)        ! (the stretched TR's emission measure)
                  p = rs*q
                  do iv = 1, euv_n_view
                    u = dot_product(p, ex(:, iv));  v = dot_product(p, ey(:, iv))
                    if (abs(u) >= fov .or. abs(v) >= fov) cycle
                    sdot = dot_product(p, d(:, iv))
                    rho2 = rs*rs - sdot*sdot
                    if (rho2 < rocc*rocc) then
                      if (sdot < sqrt(rocc*rocc - rho2)) cycle      ! behind the sphere at the join
                    end if
                    ! the sub-point as a square of side hp (its spacing, in pixels), split
                    ! over the pixels it overlaps (at most 2 x 2; pixel n covers [n-1, n])
                    fx = (u + fov)/pix - 0.5_dp*hp
                    fy = (v + fov)/pix - 0.5_dp*hp
                    ix = floor(fx);  iy = floor(fy)
                    wx(1) = min(1.0_dp, (ix + 1 - fx)/hp);  wx(2) = 1.0_dp - wx(1)
                    wy(1) = min(1.0_dp, (iy + 1 - fy)/hp);  wy(2) = 1.0_dp - wy(1)
                    do jy = 1, 2
                      if (wy(jy) <= 0.0_dp .or. iy + jy < 1 .or. iy + jy > euv_npix) cycle
                      do jx = 1, 2
                        if (wx(jx) <= 0.0_dp .or. ix + jx < 1 .or. ix + jx > euv_npix) cycle
                        img(ix+jx, iy+jy, :, iv) = img(ix+jx, iy+jy, :, iv) + wx(jx)*wy(jy)*em
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      end associate
    end do
    call cz_unbind()
    n_img = size(img)
    if (is_root) then
      call MPI_Reduce(MPI_IN_PLACE, img, n_img, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD)
    else
      call MPI_Reduce(img, img, n_img, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD)
    end if
    if (.not. is_root) return
    write(name, '(A,A,I8.8,A)') trim(out_dir), '/euv_', iStep, '.bin'
    open(newunit=u_, file=trim(name), access='stream', form='unformatted', status='replace', action='write')
    write(u_) int(euv_n_view, int32), int(nch, int32), int(euv_npix, int32), t, euv_fov_rsun, &
        euv_lon_deg(1:euv_n_view), euv_lat_deg(1:euv_n_view), int(euv_ch, int32), real(img, real32)
    close(u_)
    mx = maxval(img(:, :, 1, 1))
    write(msg, '(A,I0,A,ES11.4,A,I0,A,ES10.3,A,F7.2,A)') '  euv: step ', iStep, ', t ', t, ' s: ', euv_n_view, &
        ' view(s); max (first channel, first view) ', mx, ' DN/s/pix; ', arbor_wtime() - t0, ' s'
    call arbor_print(msg)
  contains
    logical function yy_core(t_, p_)
      real(dp), intent(in) :: t_, p_
      yy_core = t_ >= 0.25_dp*pi .and. t_ <= 0.75_dp*pi .and. abs(p_) <= 0.75_dp*pi
    end function yy_core
    ! trilinear between the centres (i0..i0+1, j0..j0+1, k0..k0+1), weights wr, wt, wp
    real(dp) function tri(a)
      real(dp), intent(in) :: a(0:, 0:, 0:)
      tri = (1.0_dp - wp)*((1.0_dp - wt)*((1.0_dp - wr)*a(i0, j0, k0) + wr*a(i0+1, j0, k0)) &
                         + wt*((1.0_dp - wr)*a(i0, j0+1, k0) + wr*a(i0+1, j0+1, k0))) &
          + wp*((1.0_dp - wt)*((1.0_dp - wr)*a(i0, j0, k0+1) + wr*a(i0+1, j0, k0+1)) &
                + wt*((1.0_dp - wr)*a(i0, j0+1, k0+1) + wr*a(i0+1, j0+1, k0+1)))
    end function tri
  end subroutine euv_coronal

  ! The images of euv_dt/euv_at_end: the coronal EUV (use_atmosphere) and the
  ! continuum (phot_n_lambda > 0)
  module subroutine euv_write(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    if (use_atmosphere) call euv_coronal(iStep, t)
    if (phot_on) call phot_write(iStep, t)
  end subroutine euv_write

  ! Synthetic continuum images (phot_n_lambda > 0; see cz3d_params_mod), the
  ! views and pixels of the EUV images. Each column (whole on its rank: the
  ! si_radial distribution) is sampled at sub-columns, ceil(3 extent/pixel)
  ! per angle (at most 12); along a sub-column ln T and ln kappa (evaluated
  ! at the centres) are bilinear in angle between the cell centres (one
  ! ghost layer) at each cell centre with phot_z_km(1) < z < phot_z_km(2). For a view with
  ! mu = r_hat . d > 0 the emergent intensity of the plane-parallel column,
  ! I = int S exp(-tau/mu) dtau/mu, S = B_lambda(T) linear in tau between the
  ! centres, from the top down to tau/mu = 40 (or the range's bottom, where
  ! S exp(-tau/mu) is added: the diffusion limit). The sub-column is placed
  ! at the radius of tau/mu = 1 (the first wavelength) and deposited as the
  ! EUV sub-points, carrying I times its projected area mu r^2 dOmega / A_pixel
  ! (Yin-Yang cover weights as EUV); each pixel's sum divided by its deposited
  ! area times min(1, area): the pixel-mean intensity (the boxes' few-% tiling
  ! unevenness removed, the limb pixels partly covered).
  ! The opacity (cm^-1, absorption only, stimulated emission included):
  ! H- bound-free and free-free per neutral H atom and electron pressure
  ! (John 1988, A&A 193, 189; the fits' temperature clipped to 1400..10080 K),
  ! hydrogen bound-free from n = 1..8 and free-free (hydrogenic, Gaunt 1),
  ! H ionization by Saha with the EOS n_e (the metals' electrons included).
  subroutine phot_write(iStep, t)
    use const_mod, only: Rsun
    use radiative_loss_mod, only: N_H_per_gram_global
    use eos_table, only: F_LNNE
    use mpi_f08, only: MPI_SUM, MPI_Reduce
    use iso_fortran_env, only: int32, real32
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    real(dp), parameter :: tau_max = 40.0_dp, mu_min = 0.01_dp
    real(dp), allocatable :: img(:,:,:,:), area(:,:,:), a(:,:,:,:), ak(:,:,:,:), rg(:), zg(:), kap(:,:), src(:,:), tb(:), &
        Iout(:), cv(:), z1(:), zs(:,:)
    real(dp) :: d(3, euv_max_view), ex(3, euv_max_view), ey(3, euv_max_view), lo, la, pix, fov, out(n_field), &
        rho, ths, phs, wt, wp, w, q(3), p(3), mu, tau, dtau, e1, ww, r1, u, v, fx, fy, hp, &
        wx(2), wy(2), dA, t0, Tg, ne, nh, th2, ph2, mx, ic_max
    integer, allocatable :: cb(:,:)
    integer :: nl, g_lo, g_hi, kb_lo, kb_hi, nk, ncol, c, jb, ib, ib1, i, j, k, g, it, ip, nst, nsp, j0, k0, &
        iv, l, ix, iy, jx, jy, n_img, u_, kk, nc
    character(len=300) :: name
    character(len=200) :: msg
    t0 = arbor_wtime()
    if (.not. allocated(ph_lam)) call phot_setup()
    nl = size(ph_lam)                ! the wavelengths
    nc = size(ph_id)                 ! the channels (bands or wavelengths)
    call cz_unbind()
    ! the radial range: whole radial blocks covering phot_z_km
    g_lo = 1;  g_hi = n_r
    do g = 1, n_r
      if (z_c(g) > phot_z_km(1)*1.0e5_dp) then
        g_lo = g
        exit
      end if
    end do
    do g = n_r, 1, -1
      if (z_c(g) < phot_z_km(2)*1.0e5_dp) then
        g_hi = g
        exit
      end if
    end do
    kb_lo = (g_lo - 1)/nb(1);  kb_hi = (g_hi - 1)/nb(1)
    nk = kb_hi - kb_lo + 1
    ! the local columns: cb(radial block - kb_lo + 1, column)
    allocate(cb(nk, nBlock))
    cb = 0
    ncol = 0
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (Key_IB(KI+1, ib) < kb_lo .or. Key_IB(KI+1, ib) > kb_hi) cycle
      c = 0
      do kk = 1, ncol
        ib1 = maxval(cb(:, kk))
        if (all(Key_IB([KP, KI+2, KI+3], ib1) == Key_IB([KP, KI+2, KI+3], ib))) then
          c = kk
          exit
        end if
      end do
      if (c == 0) then
        ncol = ncol + 1
        c = ncol
      end if
      cb(Key_IB(KI+1, ib) - kb_lo + 1, c) = ib
    end do
    if (any(cb(:, 1:ncol) == 0)) call arbor_abort('phot_write: a column is not whole on its rank (needs si_radial)')
    do iv = 1, euv_n_view
      lo = euv_lon_deg(iv)*pi/180.0_dp;  la = euv_lat_deg(iv)*pi/180.0_dp
      d(:, iv) = [cos(la)*cos(lo), cos(la)*sin(lo), sin(la)]
      ex(:, iv) = [-sin(lo), cos(lo), 0.0_dp]
      ey(:, iv) = [-sin(la)*cos(lo), -sin(la)*sin(lo), cos(la)]
    end do
    fov = euv_fov_rsun*Rsun
    pix = 2.0_dp*fov/euv_npix
    allocate(img(euv_npix, euv_npix, nc, euv_n_view), area(euv_npix, euv_npix, euv_n_view), &
        a(3, g_lo:g_hi, 0:nb(2)+1, 0:nb(3)+1), ak(nl, g_lo:g_hi, 0:nb(2)+1, 0:nb(3)+1), rg(g_lo:g_hi), zg(g_lo:g_hi), kap(nl, g_lo:g_hi), src(nl, g_lo:g_hi), &
        tb(g_lo:g_hi), Iout(nl), cv(nc), z1(nl), zs(nl, 3))
    zs = 0.0_dp
    img = 0.0_dp
    area = 0.0_dp
    ic_max = 0.0_dp
    call arbor_fill_ghosts(State_GVB, n_ghost=1, faces_only=.false.)
    do c = 1, ncol
      ! ln T, ln n_e, ln n_H at the centres (one horizontal ghost layer)
      do kk = 1, nk
        ib = cb(kk, c)
        call cz_bind(ib)
        do i = 1, nb(1)
          g = r_off + i
          if (g < g_lo .or. g > g_hi) cycle
          rg(g) = r_c(i);  zg(g) = z_c(i)
          do k = 0, nb(3) + 1
            do j = 0, nb(2) + 1
              rho = max(rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib), 1.0e-30_dp)
              call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
              a(:, g, j, k) = [out(F_LNT), out(F_LNNE), log(N_H_per_gram_global*rho)]
              do l = 1, nl
                ak(l, g, j, k) = log(max(kappa_cont(ph_lam(l), exp(out(F_LNT)), exp(out(F_LNNE)), &
                    N_H_per_gram_global*rho), 1.0e-300_dp))
              end do
            end do
          end do
        end do
      end do
      ib = cb(nk, c)
      associate (th => Metric_DB(2, ib)%X_I, dth_ => Metric_DB(2, ib)%Dx_I, ph => Metric_DB(3, ib)%X_I, &
                 dph_ => Metric_DB(3, ib)%Dx_I, vth => Metric_DB(2, ib)%Vol_I, vph => Metric_DB(3, ib)%Vol_I)
      if (lbound(th, 1) > 0 .or. ubound(th, 1) < nb(2) + 1 .or. lbound(ph, 1) > 0 .or. ubound(ph, 1) < nb(3) + 1) &
          call arbor_abort('phot_write: the block metric has no ghost centres')
      do k = 1, nb(3)
        do j = 1, nb(2)
          nst = min(12, max(1, ceiling(3.0_dp*Rsun*dth_(j)/pix)))
          nsp = min(12, max(1, ceiling(3.0_dp*Rsun*sin(th(j))*dph_(k)/pix)))
          dA = vth(j)*vph(k)/real(nst*nsp, dp)                         ! solid angle per sub-column
          do ip = 1, nsp
            phs = ph(k) + ((ip - 0.5_dp)/nsp - 0.5_dp)*dph_(k)
            k0 = merge(k - 1, k, phs < ph(k))
            wp = (phs - ph(k0))/(ph(k0+1) - ph(k0))
            do it = 1, nst
              ths = th(j) + ((it - 0.5_dp)/nst - 0.5_dp)*dth_(j)
              if (.not. yy_core(ths, phs)) cycle
              w = 1.0_dp
              call yy_angles(ths, phs, th2, ph2)
              if (yy_core(th2, ph2)) w = 0.5_dp
              j0 = merge(j - 1, j, ths < th(j))
              wt = (ths - th(j0))/(th(j0+1) - th(j0))
              q = [sin(ths)*cos(phs), sin(ths)*sin(phs), cos(ths)]
              if (Key_IB(KP, ib) == 2) q = yy_cart(q)
              if (all([(dot_product(q, d(:, iv)) <= mu_min, iv = 1, euv_n_view)])) cycle
              ! the sub-column's opacity and source function
              do g = g_lo, g_hi
                Tg = exp(bil(1, g))
                tb(g) = Tg
                do l = 1, nl
                  kap(l, g) = exp((1.0_dp - wp)*((1.0_dp - wt)*ak(l, g, j0, k0) + wt*ak(l, g, j0+1, k0)) &
                                  + wp*((1.0_dp - wt)*ak(l, g, j0, k0+1) + wt*ak(l, g, j0+1, k0+1)))
                  src(l, g) = planck_A(ph_lam(l), Tg)
                end do
              end do
              do iv = 1, euv_n_view
                mu = dot_product(q, d(:, iv))
                if (mu <= mu_min) cycle
                z1 = huge(1.0_dp)
                do l = 1, nl
                  tau = 0.0_dp
                  Iout(l) = 0.0_dp
                  do g = g_hi, g_lo + 1, -1
                    dtau = 0.5_dp*(kap(l, g) + kap(l, g-1))*(rg(g) - rg(g-1))/mu
                    e1 = -expm1(-dtau)
                    if (dtau > 1.0e-4_dp) then
                      ww = (e1 - dtau*exp(-dtau))/dtau
                    else
                      ww = 0.5_dp*dtau
                    end if
                    Iout(l) = Iout(l) + exp(-tau)*(src(l, g)*e1 + (src(l, g-1) - src(l, g))*ww)
                    if (z1(l) == huge(1.0_dp) .and. tau + dtau >= 1.0_dp) &
                        z1(l) = rg(g) - (1.0_dp - tau)/dtau*(rg(g) - rg(g-1))       ! (the radius of tau/mu = 1)
                    tau = tau + dtau
                    if (tau > tau_max) exit
                  end do
                  if (tau <= tau_max) Iout(l) = Iout(l) + exp(-tau)*src(l, g_lo)
                end do
                where (z1 == huge(1.0_dp)) z1 = rg(g_lo)
                r1 = z1(1)
                cv = 0.0_dp
                do l = 1, nl
                  cv(ph_chn(l)) = cv(ph_chn(l)) + ph_wgt(l)*Iout(l)
                end do
                if (mu > 0.999_dp) then          ! disk centre: the tau = 1 heights, I and its spread
                  ic_max = max(ic_max, cv(1))
                  zs(:, 1) = zs(:, 1) + 1.0_dp
                  zs(:, 2) = zs(:, 2) + (z1 - rg(g_lo)) + zg(g_lo)
                  zs(:, 3) = zs(:, 3) + Iout
                end if
                p = r1*q
                u = dot_product(p, ex(:, iv));  v = dot_product(p, ey(:, iv))
                if (abs(u) >= fov .or. abs(v) >= fov) cycle
                hp = min(1.0_dp, r1*sqrt(dA)/pix)
                fx = (u + fov)/pix - 0.5_dp*hp
                fy = (v + fov)/pix - 0.5_dp*hp
                ix = floor(fx);  iy = floor(fy)
                wx(1) = min(1.0_dp, (ix + 1 - fx)/hp);  wx(2) = 1.0_dp - wx(1)
                wy(1) = min(1.0_dp, (iy + 1 - fy)/hp);  wy(2) = 1.0_dp - wy(1)
                do jy = 1, 2
                  if (wy(jy) <= 0.0_dp .or. iy + jy < 1 .or. iy + jy > euv_npix) cycle
                  do jx = 1, 2
                    if (wx(jx) <= 0.0_dp .or. ix + jx < 1 .or. ix + jx > euv_npix) cycle
                    ww = wx(jx)*wy(jy)*w*mu*r1*r1*dA/pix**2
                    img(ix+jx, iy+jy, :, iv) = img(ix+jx, iy+jy, :, iv) + ww*cv
                    area(ix+jx, iy+jy, iv) = area(ix+jx, iy+jy, iv) + ww
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      end associate
    end do
    call cz_unbind()
    n_img = size(img)
    if (is_root) then
      call MPI_Reduce(MPI_IN_PLACE, img, n_img, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD)
      call MPI_Reduce(MPI_IN_PLACE, area, size(area), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD)
    else
      call MPI_Reduce(img, img, n_img, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD)
      call MPI_Reduce(area, area, size(area), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD)
    end if
    ! the mean intensity of what each pixel covers, times the covered fraction
    ! (at most 1): the sub-columns' boxes do not tile the pixels exactly (a few %)
    do iv = 1, euv_n_view
      do l = 1, nc
        where (area(:, :, iv) > 0.0_dp) img(:, :, l, iv) = img(:, :, l, iv)/area(:, :, iv)*min(1.0_dp, area(:, :, iv))
      end do
    end do
    ic_max = arbor_max(ic_max)
    call arbor_sum_inplace(zs)
    if (.not. is_root) return
    if (len_trim(phot_band_file) > 0) then
      write(name, '(A,A,I8.8,A)') trim(out_dir), '/photb_', iStep, '.bin'
    else
      write(name, '(A,A,I8.8,A)') trim(out_dir), '/phot_', iStep, '.bin'
    end if
    open(newunit=u_, file=trim(name), access='stream', form='unformatted', status='replace', action='write')
    write(u_) int(euv_n_view, int32), int(nc, int32), int(euv_npix, int32), t, euv_fov_rsun, &
        euv_lon_deg(1:euv_n_view), euv_lat_deg(1:euv_n_view), int(ph_id, int32), real(img, real32)
    close(u_)
    mx = maxval(img(:, :, 1, 1))
    write(msg, '(A,I0,A,ES11.4,A,I0,A,ES10.3,A,ES10.3,A,A,A,F7.2,A)') '  phot: step ', iStep, ', t ', t, ' s: ', &
        ph_id(1), ', max at mu = 1 ', ic_max, ', image max ', mx, ' ', &
        trim(merge('DN/s/pix       ', 'erg/s/cm2/sr/A ', len_trim(phot_band_file) > 0)), '; ', arbor_wtime() - t0, ' s'
    call arbor_print(msg)
    ! the disk-centre tau = 1 heights and brightness temperatures (the channel's
    ! weighted mean over its wavelengths)
    if (zs(1, 1) > 0.0_dp) then
      do l = 1, nc
        write(msg, '(A,I6,A,F8.1,A,F8.0,A)') '    phot ', ph_id(l), ': disk centre tau = 1 at z ', &
            sum(ph_wgt*zs(:, 2), mask=ph_chn == l)/sum(ph_wgt, mask=ph_chn == l)/zs(1, 1)/1.0e5_dp, ' km, T_b ', &
            sum(ph_wgt*tbright(ph_lam, zs(:, 3)/zs(1, 1)), mask=ph_chn == l)/sum(ph_wgt, mask=ph_chn == l), ' K'
        call arbor_print(msg)
      end do
    end if
  contains
    logical function yy_core(t_, p_)
      real(dp), intent(in) :: t_, p_
      yy_core = t_ >= 0.25_dp*pi .and. t_ <= 0.75_dp*pi .and. abs(p_) <= 0.75_dp*pi
    end function yy_core
    real(dp) function bil(m, g_)
      integer, intent(in) :: m, g_
      bil = (1.0_dp - wp)*((1.0_dp - wt)*a(m, g_, j0, k0) + wt*a(m, g_, j0+1, k0)) &
          + wp*((1.0_dp - wt)*a(m, g_, j0, k0+1) + wt*a(m, g_, j0+1, k0+1))
    end function bil
  end subroutine phot_write

  ! The continuum channels: phot_lambda_A(1:phot_n_lambda) (weight 1, the
  ! intensity), or the bands of phot_band_file (tools/aia_uv_response.py):
  ! '#' comment lines, n_band, then per band "id n" and n lines "lambda w"
  ! (DN/s/pix = sum w I(lambda))
  subroutine phot_setup()
    integer :: u, ios, nb_, m, n, l, id
    character(len=512) :: line
    if (len_trim(phot_band_file) == 0) then
      allocate(ph_lam(phot_n_lambda), ph_wgt(phot_n_lambda), ph_chn(phot_n_lambda), ph_id(phot_n_lambda))
      ph_lam = phot_lambda_A(1:phot_n_lambda)
      ph_wgt = 1.0_dp
      ph_chn = [(l, l = 1, phot_n_lambda)]
      ph_id = nint(ph_lam)
      return
    end if
    open(newunit=u, file=trim(phot_band_file), status='old', action='read', iostat=ios)
    if (ios /= 0) call arbor_abort('phot: cannot open phot_band_file '//trim(phot_band_file))
    do
      read(u, '(A)') line
      if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
    end do
    read(line, *) nb_
    allocate(ph_id(nb_), ph_lam(0), ph_wgt(0), ph_chn(0))
    do m = 1, nb_
      read(u, *) id, n
      ph_id(m) = id
      do l = 1, n
        read(u, '(A)') line
        block
          real(dp) :: x(2)
          read(line, *) x
          if (x(1) < 1300.0_dp) call arbor_abort('phot: band wavelengths must be >= 1300 A')
          ph_lam = [ph_lam, x(1)];  ph_wgt = [ph_wgt, x(2)];  ph_chn = [ph_chn, m]
        end block
      end do
    end do
    close(u)
  end subroutine phot_setup

  ! the brightness temperature [K] of the intensity I [erg s^-1 cm^-2 sr^-1 A^-1] at lambda [A]
  elemental real(dp) function tbright(lam, I)
    real(dp), intent(in) :: lam, I
    real(dp) :: lc
    lc = lam*1.0e-8_dp
    tbright = 1.4387770_dp/lc/log(1.0_dp + 1.1910430e-5_dp/lc**5*1.0e-8_dp/max(I, tiny(1.0_dp)))
  end function tbright

  ! B_lambda(T) [erg s^-1 cm^-2 sr^-1 A^-1], lambda [A]
  pure real(dp) function planck_A(lam, T)
    real(dp), intent(in) :: lam, T
    real(dp), parameter :: c1 = 1.1910430e-5_dp, c2 = 1.4387770_dp      ! 2 h c^2 [erg cm^2 s^-1], h c / k [cm K]
    real(dp) :: lc
    lc = lam*1.0e-8_dp
    planck_A = c1/lc**5/expm1(c2/(lc*T))*1.0e-8_dp
  end function planck_A

  ! LTE continuum absorption coefficient [cm^-1] at lambda [A] (>= 1300) for
  ! T [K], n_e and n_H (all hydrogen) [cm^-3]: H- bound-free and free-free
  ! (John 1988; the free-free fit only above 3645 A, below it is negligible
  ! against the metals), H bound-free (n = 1..8) and free-free (hydrogenic,
  ! Gaunt 1), metal bound-free below 2600 A (kappa_metal); n_HI from Saha
  ! with n_e (see phot_write)
  pure real(dp) function kappa_cont(lam, T, ne, nh)
    real(dp), intent(in) :: lam, T, ne, nh
    real(dp), parameter :: kB = 1.380649e-16_dp, lam0 = 1.6419_dp, alpha = 1.439e4_dp, chi_eV = 13.598_dp, &
        kB_eV = 8.617333e-5_dp, cl = 2.99792458e10_dp
    real(dp), parameter :: cb(6) = [152.519_dp, 49.534_dp, -118.858_dp, 92.536_dp, -34.194_dp, 4.982_dp]
    real(dp), parameter :: ff(5, 6) = reshape([ &
        2483.346_dp, -3449.889_dp, 2200.040_dp, -696.271_dp, 88.283_dp, &
        285.827_dp, -1158.382_dp, 2427.719_dp, -1841.400_dp, 444.517_dp, &
        -2054.291_dp, 8746.523_dp, -13651.105_dp, 8624.970_dp, -1863.864_dp, &
        2827.776_dp, -11485.632_dp, 16755.524_dp, -10051.530_dp, 2095.288_dp, &
        -1341.537_dp, 5303.609_dp, -7510.494_dp, 4400.067_dp, -901.788_dp, &
        208.952_dp, -812.939_dp, 1132.738_dp, -655.020_dp, 132.985_dp], [5, 6])
    real(dp) :: lu, x, sig, f, Tf, th, kbf, kff, saha, nhi, np, pe, nu, stim, hb, chin
    integer :: n
    lu = lam*1.0e-4_dp                                  ! [micron]
    Tf = min(max(T, 1400.0_dp), 10080.0_dp)
    th = 5040.0_dp/Tf
    saha = 2.4148e15_dp*T**1.5_dp*exp(-157807.0_dp/T)   ! n_p n_e / n_HI
    nhi = nh*ne/(ne + saha)
    np = nh - nhi
    pe = ne*kB*T
    ! H- bound-free [cm^4 dyn^-1 per H atom]
    kbf = 0.0_dp
    if (lu < lam0) then
      x = 1.0_dp/lu - 1.0_dp/lam0
      f = 0.0_dp
      do n = 1, 6
        f = f + cb(n)*x**(0.5_dp*(n - 1))
      end do
      sig = 1.0e-18_dp*lu**3*x**1.5_dp*f
      kbf = 0.750_dp*Tf**(-2.5_dp)*exp(alpha/(lam0*Tf))*(1.0_dp - exp(-alpha/(lu*Tf)))*sig
    end if
    ! H- free-free (lambda > 0.3645 micron)
    kff = 0.0_dp
    if (lu > 0.3645_dp) then
    do n = 1, 6
      kff = kff + th**(0.5_dp*(n + 1))*(ff(1, n) + ff(2, n)/lu + ff(3, n)/lu**2 + ff(4, n)/lu**3 + ff(5, n)/lu**4)
    end do
    end if
    kff = 1.0e-29_dp*kff
    kappa_cont = (kbf + max(kff, 0.0_dp))*pe*nhi
    ! hydrogen bound-free (levels whose edge lies beyond lambda) and free-free
    nu = cl/(lam*1.0e-8_dp)
    stim = -expm1(-1.4387770e8_dp/(lam*T))
    hb = 0.0_dp
    do n = 1, 8
      if (911.75_dp*n*n < lam) cycle
      chin = chi_eV*(1.0_dp - 1.0_dp/(n*n))
      hb = hb + 2.815e29_dp/(real(n, dp)**5*nu**3)*real(n*n, dp)*exp(-chin/(kB_eV*T))
    end do
    kappa_cont = kappa_cont + hb*nhi*stim + 3.69e8_dp/sqrt(T)/nu**3*ne*np*stim
    ! metal bound-free (lambda < 2600 A)
    if (lam < 2600.0_dp) kappa_cont = kappa_cont + kappa_metal(lam, T, ne, nh)*stim
  end function kappa_cont

  ! Metal bound-free absorption [cm^-1] (without stimulated emission) for the
  ! UV continua (AIA 1600/1700): the low levels of Si I, Mg I, Al I, Fe I,
  ! Ca I, C I with Boltzmann populations over the listed levels (their sum
  ! the partition function U_I), Saha with n_e and a constant U_II for the
  ! neutral fraction, Asplund 2009 abundances, and the cross section
  ! sigma = sigma_0 (lambda/lambda_edge)^p for lambda below the level's edge.
  ! sigma_0 and p are rough threshold values (Si I 3P/1D/1S ~ 37/35/40 Mb in
  ! the range of R-matrix results; Fe I a flat 4 Mb in place of its
  ! resonances): the opacity is good to a factor ~2, i.e. the tau = 1
  ! height to ~0.7 scale heights. Si I 1D (edge 1682 A) and 3P (1521 A)
  ! dominate the 1600/1700 bands.
  pure real(dp) function kappa_metal(lam, T, ne, nh)
    real(dp), intent(in) :: lam, T, ne, nh
    integer, parameter :: n_el = 6, n_lev = 13
    real(dp), parameter :: kB_eV = 8.617333e-5_dp
    !                                     Si       Mg       Al       Fe       Ca       C
    real(dp), parameter :: abund(n_el) = [7.51_dp, 7.60_dp, 6.45_dp, 7.50_dp, 6.34_dp, 8.43_dp]
    real(dp), parameter :: ip_eV(n_el) = [8.1517_dp, 7.6462_dp, 5.9858_dp, 7.9024_dp, 6.1132_dp, 11.2603_dp]
    real(dp), parameter :: u2(n_el) = [5.7_dp, 2.0_dp, 1.0_dp, 42.0_dp, 2.2_dp, 5.9_dp]
    ! the levels: element, excitation [eV], g, sigma_0 [Mb], p
    integer, parameter :: lev_el(n_lev) = [1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 6, 6, 6]
    real(dp), parameter :: lev_e(n_lev) = [0.0_dp, 0.781_dp, 1.909_dp, 0.0_dp, 2.709_dp, 0.0_dp, 0.0_dp, 0.86_dp, &
        1.49_dp, 0.0_dp, 0.0_dp, 1.264_dp, 2.684_dp]
    real(dp), parameter :: lev_g(n_lev) = [9.0_dp, 5.0_dp, 1.0_dp, 1.0_dp, 9.0_dp, 6.0_dp, 25.0_dp, 35.0_dp, 21.0_dp, &
        1.0_dp, 9.0_dp, 5.0_dp, 1.0_dp]
    real(dp), parameter :: lev_s(n_lev) = [37.0_dp, 35.0_dp, 40.0_dp, 1.5_dp, 20.0_dp, 65.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, &
        2.0_dp, 12.0_dp, 10.0_dp, 9.0_dp]
    real(dp), parameter :: lev_p(n_lev) = [1.5_dp, 1.5_dp, 1.5_dp, 0.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.5_dp, 1.5_dp, 1.5_dp]
    real(dp) :: u1(n_el), fn(n_el), kT, edge
    integer :: e, m
    kT = kB_eV*T
    u1 = 0.0_dp
    do m = 1, n_lev
      u1(lev_el(m)) = u1(lev_el(m)) + lev_g(m)*exp(-lev_e(m)/kT)
    end do
    ! n_I / n_H (the neutral fraction times the abundance) / U_I
    do e = 1, n_el
      fn(e) = 10.0_dp**(abund(e) - 12.0_dp)*ne/(ne + 2.4148e15_dp*T**1.5_dp*2.0_dp*u2(e)/u1(e)*exp(-ip_eV(e)/kT))/u1(e)
    end do
    kappa_metal = 0.0_dp
    do m = 1, n_lev
      e = lev_el(m)
      edge = 12398.42_dp/(ip_eV(e) - lev_e(m))
      if (lam > edge) cycle
      kappa_metal = kappa_metal + fn(e)*lev_g(m)*exp(-lev_e(m)/kT)*lev_s(m)*1.0e-18_dp*(lam/edge)**lev_p(m)
    end do
    kappa_metal = kappa_metal*nh
  end function kappa_metal

  ! euv_table: the AIA temperature responses (see euv_coronal)
  subroutine euv_load()
    integer :: u, nT, nch, m, ios
    character(len=512) :: line
    open(newunit=u, file=trim(euv_table), status='old', action='read', iostat=ios)
    if (ios /= 0) call arbor_abort('euv: cannot open euv_table '//trim(euv_table))
    do
      read(u, '(A)') line
      if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
    end do
    read(line, *) nT, nch
    allocate(euv_lgT(nT), euv_K(nT, nch), euv_ch(nch))
    read(u, *) euv_ch
    do m = 1, nT
      read(u, *) euv_lgT(m), euv_K(m, :)
    end do
    close(u)
    if (nT < 2 .or. maxval(abs((euv_lgT(2:) - euv_lgT(:nT-1)) - (euv_lgT(2) - euv_lgT(1)))) > 1.0e-6_dp) &
        call arbor_abort('euv: euv_table must be uniform in log10 T')
  end subroutine euv_load

end submodule cz3d_synth_smod
