submodule (cz3d_app_mod) cz3d_setup_smod
  ! setup: initialisation, field seeds (PFSS, DPF, magnetofriction, smoothing), state import/export and matching, the column distribution, restart files
  ! (part of cz3d_app_mod: the module's data and routines by host association)
  implicit none

contains

  ! The number of state variables (5 hydro; 9 with the magnetic field)
  module integer function cz_n_var()
    cz_n_var = merge(nVarB, nVarCz, use_B)
  end function cz_n_var

  ! After arbor_block_sync: per-block arrays, the subgrid factor, scratch
  ! Refinement to refine_h levels in theta and phi everywhere and to
  ! patch_level near the patch centre (arbor_regrid)
  module subroutine cz_flag_block(self, iBlock, flags)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iBlock
    integer, intent(out) :: flags(3)
    real(dp) :: x0(3), x(3), e(3,3), la, lo
    integer :: j, k
    flags = 0
    if (Key_IB(KL+2, iBlock) < refine_h) flags(2:3) = 1
    if (Key_IB(KL+2, iBlock) >= patch_level) return
    la = patch_lat_deg*pi/180.0_dp;  lo = patch_lon_deg*pi/180.0_dp
    x0 = [cos(la)*cos(lo), cos(la)*sin(lo), sin(la)]
    do k = 1, nb(3)
      do j = 1, nb(2)
        e = sph_basis(Metric_DB(2, iBlock)%X_I(j), Metric_DB(3, iBlock)%X_I(k))
        x = e(:, 1)
        if (Key_IB(KP, iBlock) == 2) x = yy_cart(x)
        if (dot_product(x, x0) >= cos(patch_radius_deg*pi/180.0_dp)) then
          flags(2:3) = 1
          return
        end if
      end do
    end do
  end subroutine cz_flag_block

  ! After arbor_block_sync or a regrid (blocks may have changed): the
  ! per-block arrays (zeroed: the convective flux restarts at its targets),
  ! the subgrid factor, the scratch arrays (once)
  module subroutine cz_setup()
    integer :: i, lev
    real(dp) :: Hp, R
    if (mod(n_r, nb(1)) /= 0) call arbor_abort('cz3d: the radial block size must divide the radial cells')
    if (nb(1) < n_r .and. dt_budget > 0.0_dp) call arbor_abort('cz3d: the energy budget needs whole columns per block')
    if (mr_levels > 1 .and. (dt1 > 0.0_dp .or. dt_budget > 0.0_dp .or. lz_correct_steps > 0 .or. lz_correct_dt > 0.0_dp)) &
        call arbor_abort('cz3d: multirate needs dt1 = 0 and no energy budget or L_z correction')
    call cz_unbind()
    rsst_h2012 = rsst_form == 'hotta2012'
    glm_cfl = glm_ch == 'cfl'
    grey_i0_g = 1                                 ! first cell of the implicit grey solve (whole column)
    if (grey_top_km > 0.0_dp) then
      ! (a deep top, the cascade: no cell within grey_top_km of the photosphere;
      ! explicit radiative diffusion in the whole column)
      grey_i0_g = n_r + 1
      do i = 1, n_r
        if (-z_c(i) < grey_top_km*1.0e5_dp) then
          grey_i0_g = i
          exit
        end if
      end do
    end if
    if (use_atmosphere .and. mod(k_join, nb(1)) /= 0) call arbor_abort('cztr3d: the join cell must end a radial block '// &
        '(k_join a multiple of the radial block size)')
    if (nb(1) < merge(k_join, n_r, use_atmosphere) .and. grey_i0_g <= merge(k_join, n_r, use_atmosphere) - nb(1)) &
        call arbor_abort('cz3d: with radial blocks the implicit grey layer '// &
        '(grey_top_km) must lie within the top block')
    grey_i0 = grey_i0_g
    rsst_h2015 = rsst_form == 'hotta2015'
    if (use_atmosphere) call atm_setup()
    if (allocated(frz_B)) deallocate(frz_B, frz_have)
    if (freeze_side /= 'none') then
      allocate(frz_B(nb(1), nb(2), nb(3), cz_n_var(), nBlockMax), frz_have(nBlockMax))
      frz_have = .false.
    end if
    if (allocated(qrad_FB)) deallocate(qrad_FB, qdif_FB)
    if (allocated(qconv_FB)) deallocate(qconv_FB)
    if (associated(f_mlt_g)) deallocate(f_mlt_g)
    if (.not. qconv_face) then
      allocate(qconv_FB(0:nb(1), nb(2), nb(3), nBlockMax))
      qconv_FB = 0.0_dp
    end if
    allocate(qrad_FB(0:nb(1), nb(2), nb(3), nBlockMax))
    allocate(qdif_FB(0:nb(1), nb(2), nb(3), nBlockMax))
    qrad_FB = 0.0_dp
    qdif_FB = 0.0_dp
    ! the subgrid factor per horizontal level (the cell size halves per level)
    allocate(f_mlt_g(0:n_r, 0:max(refine_h, patch_level)))
    f_mlt_g = 1.0_dp
    do lev = 0, ubound(f_mlt_g, 2)
      dtheta = 0.5_dp*pi/real(n_theta - 2*merge(yy_margin, ng + 1, yy_margin >= 0), dp)/2.0_dp**lev
      if (.not. mlt_switch) cycle
      do i = 0, n_r
        ! background pressure scale height at the face per horizontal cell size
        Hp = 0.5_dp*(p0(i)/(rho0(i)*g_c(i)) + p0(i+1)/(rho0(i+1)*g_c(i+1)))
        R = Hp/(r_f(i)*dtheta)
        f_mlt_g(i, lev) = 1.0_dp/(1.0_dp + (R/mlt_R0)**mlt_p)
        ! under the top wall (mlt_wall_hp): the larger of the two
        if (mlt_wall_hp > 0.0_dp) f_mlt_g(i, lev) = max(f_mlt_g(i, lev), 1.0_dp/(1.0_dp + ((r_f(n_r) - r_f(i)) &
            /(mlt_wall_hp*p0(n_r)/(rho0(n_r)*g_c(n_r))))**mlt_p))
      end do
    end do
    lev_view = refine_h
    f_mlt(0:n_r) => f_mlt_g(0:n_r, refine_h)
    if (.not. allocated(tloc)) then
      allocate(tloc(1:n_r))
      tloc = 0.0_dp
    end if
    if (si_radial .and. si_cols) call si_setup()     ! (after cz_columns_to_ranks)
    if (allocated(rho_G)) return
    allocate(tsc_g(1:n_r), rate_r_g(1:n_r))
    tsc => tsc_g
    rate_r => rate_r_g
    tsc = 1.0_dp
    rate_r = 0.0_dp
    allocate(rq_g(-1:n_r+2))
    rq => rq_g
    rq = 1.0_dp
    if (rsst_h2015) rq = 1.0_dp/xi**2
    call background_h_phi()
    allocate(rho0_w(-1:n_r+2), rq_w(-1:n_r+2), s0_w(-1:n_r+2), r_w(-1:n_r+2))
    rho0_w = rho0(-1:n_r+2)
    s0_w = s0(-1:n_r+2)
    r_w = r_c(-1:n_r+2)
    rq_w = rq(-1:n_r+2)
    allocate(rho_G(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng))
    allocate(p1_G, cs_G, W_G, lim_G, T_G, hP_G, mold=rho_G)
    nv = cz_n_var()
    ie_ = nv + 1
    allocate(D_C(nb(1), nb(2), nb(3)), R_CV(nb(1), nb(2), nb(3), nv))
    allocate(Dr_C, mold=D_C)
    ! (variable ie_ = nv + 1: the diffusion's energy flux)
    allocate(F1_FV(0:nb(1), nb(2), nb(3), ie_), F2_FV(nb(1), 0:nb(2), nb(3), ie_), &
        F3_FV(nb(1), nb(2), 0:nb(3), ie_), Div_CV(nb(1), nb(2), nb(3), ie_))
    allocate(C1_F(0:nb(1), nb(2), nb(3)), C2_F(nb(1), 0:nb(2), nb(3)), C3_F(nb(1), nb(2), 0:nb(3)))
    if (use_B) allocate(B_GV(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng, 3), C_CV(nb(1), nb(2), nb(3), 3), &
        J_CV(nb(1), nb(2), nb(3), 3), divB_C(nb(1), nb(2), nb(3)), vA_G(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng))
  end subroutine cz_setup

  ! Initial state: the 1D state in every column, plus an optional random s1
  ! (from each cell's global index, so independent of the number of ranks)
  module subroutine cz_init_block(self, iBlock)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iBlock
    integer :: i, j, k
    call cz_bind(iBlock)
    State_GVB(:, :, :, :, iBlock) = 0.0_dp
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, nb(1)
          State_GVB(i, j, k, rho1_, iBlock) = rho1_init(i)/rq(i)
          State_GVB(i, j, k, vr_, iBlock) = v_init(i)
          State_GVB(i, j, k, s1_, iBlock) = s1_init(i)
          if (ds1_random > 0.0_dp .and. (.not. use_atmosphere .or. r_off + i <= k_join)) &
              State_GVB(i, j, k, s1_, iBlock) = s1_init(i) + ds1_random &
              *urand(Key_IB(KP, iBlock), r_off + i, Key_IB(KI+2, iBlock)*nb(2) + j, Key_IB(KI+3, iBlock)*nb(3) + k)
        end do
      end do
    end do
    if (len_trim(import_dir) > 0) call import_block(iBlock)
    if (use_B .and. b_seed == 'toroidal') call seed_toroidal(iBlock)
    if (use_B .and. b_seed == 'dipole') call seed_dipole(iBlock)
    call cz_unbind()
  end subroutine cz_init_block

  ! Toroidal seed field about the global z axis, added to the block's field:
  ! B = b_seed_G w(r) sin(theta_g) e_phi_g = b_seed_G w(r) (-y, x, 0) for the
  ! unit position (x, y, z) in the global (Yin) frame; divergence free. w = 1
  ! below b_seed_r1_rsun, 0 above b_seed_r2_rsun, sin^2 between.
  subroutine seed_toroidal(iBlock)
    use const_mod, only: Rsun
    integer, intent(in) :: iBlock
    real(dp) :: e(3,3), x(3), bg(3), bl(3), w, r1, r2
    integer :: i, j, k
    r1 = b_seed_r1_rsun*Rsun
    r2 = b_seed_r2_rsun*Rsun
    do k = 1, nb(3)
      do j = 1, nb(2)
        e = sph_basis(Metric_DB(2, iBlock)%X_I(j), Metric_DB(3, iBlock)%X_I(k))
        x = e(:, 1)
        if (Key_IB(KP, iBlock) == 2) x = yy_cart(x)          ! the position in the global frame
        bg = [-x(2), x(1), 0.0_dp]
        if (Key_IB(KP, iBlock) == 2) bg = yy_cart(bg)        ! back to this patch's frame
        bl = matmul(transpose(e), bg)                         ! (B_r, B_theta, B_phi) per unit amplitude
        do i = 1, nb(1)
          if (r_c(i) <= r1) then
            w = 1.0_dp
          else if (r_c(i) >= r2) then
            w = 0.0_dp
          else
            w = cos(0.5_dp*pi*(r_c(i) - r1)/(r2 - r1))**2
          end if
          State_GVB(i, j, k, br_:bp_, iBlock) = State_GVB(i, j, k, br_:bp_, iBlock) + b_seed_G*w*bl
        end do
      end do
    end do
  end subroutine seed_toroidal

  ! Dipole seed (b_seed = 'dipole'): the l = 1 potential field about the
  ! global z axis with a source surface at r_ss = b_seed_rss_rsun, Phi = A
  ! cos(theta) [(r_ss/r)^2 - r/r_ss]: B_r = A cos(theta) (2 r_ss^2/r^3 +
  ! 1/r_ss), B_theta = A sin(theta) (r_ss^2/r^3 - 1/r_ss) (0 at r_ss); above r_ss
  ! radial, B_r(r_ss, theta) (r_ss/r)^2; A from b_seed_G = B_r at the pole at
  ! r_top_cz. Added to the block's field (theta: the global colatitude)
  subroutine seed_dipole(iBlock)
    use const_mod, only: Rsun
    integer, intent(in) :: iBlock
    real(dp) :: e(3,3), x(3), bg(3), bl(3), rss, A, ct, st, r, br, bth, rr, thg(3)
    integer :: i, j, k
    rss = b_seed_rss_rsun*Rsun
    A = b_seed_G/(2.0_dp*rss**2/r_top_cz**3 + 1.0_dp/rss)
    do k = 1, nb(3)
      do j = 1, nb(2)
        e = sph_basis(Metric_DB(2, iBlock)%X_I(j), Metric_DB(3, iBlock)%X_I(k))
        x = e(:, 1)
        if (Key_IB(KP, iBlock) == 2) x = yy_cart(x)          ! the position in the global frame
        ct = x(3)
        st = sqrt(max(1.0_dp - ct*ct, 0.0_dp))
        thg = [0.0_dp, 0.0_dp, -1.0_dp]                       ! e_theta of the global frame (at the poles: any)
        if (st > 1.0e-12_dp) thg = [ct*x(1)/st, ct*x(2)/st, -st]
        do i = 1, nb(1)
          r = r_c(i)
          if (r <= rss) then
            br = A*ct*(2.0_dp*rss**2/r**3 + 1.0_dp/rss)
            bth = A*st*(rss**2/r**3 - 1.0_dp/rss)
          else
            rr = rss/r
            br = A*ct*3.0_dp/rss*rr**2
            bth = 0.0_dp
          end if
          bg = br*x + bth*thg                                 ! global Cartesian
          if (Key_IB(KP, iBlock) == 2) bg = yy_cart(bg)      ! to this patch's frame
          bl = matmul(transpose(e), bg)
          State_GVB(i, j, k, br_:bp_, iBlock) = State_GVB(i, j, k, br_:bp_, iBlock) + bl
        end do
      end do
    end do
  end subroutine seed_dipole

  ! A hydrodynamic restart read into an MHD run (arbor_restart_read with
  ! n_var_file < nv): B and psi came back zero; add the seed field to every
  ! block (init_block is not called on a restart)
  module subroutine cz_field_from_hydro()
    integer :: jb
    character(len=200) :: msg
    if (.not. use_B) return
    do jb = 1, nBlock
      call cz_bind(iBlockList_I(jb))
      State_GVB(:, :, :, br_:psi_, iBlockList_I(jb)) = 0.0_dp
      if (b_seed == 'toroidal') call seed_toroidal(iBlockList_I(jb))
      if (b_seed == 'dipole') call seed_dipole(iBlockList_I(jb))
    end do
    call cz_unbind()
    write(msg, '(A,A,A,F8.2,A)') '  magnetic field added to a hydrodynamic restart: seed ', trim(b_seed), ' (', b_seed_G, ' G)'
    call arbor_print(msg)
  end subroutine cz_field_from_hydro

  ! b_seed = 'pfss': the field above the join (cells > k_join) replaced by the
  ! PFSS field of the join cell's B_r (pfss_mod; R0 = that cell's centre, the
  ! source surface at b_seed_rss_rsun, radial above it). The coefficients
  ! come from a midpoint quadrature over the join shell (solid-angle weights,
  ! 1/2 where both Yin-Yang cores cover a point, as dump_shells_write), refined
  ! pfss_iter times on the residual, summed over the ranks; the monopole is
  ! dropped. Each block column then sums the harmonics once (O(L^2)) and each
  ! cell its l series (O(L)), with l cut where (R0/r)^l < 1e-10. psi = 0
  ! above the join. pfss_check: the CZ field set to the analytic dipole
  ! first, and the result compared with it.
  module subroutine cz_seed_pfss()
    use const_mod, only: Rsun
    real(dp), allocatable :: c(:,:), s(:,:), dc(:,:), ds(:,:), p(:,:), dpt(:,:), Sl(:), Tl(:), Ul(:), &
        fr(:,:), gh(:,:), frs(:), ghs(:), cmp(:), smp(:), Bsave(:,:,:,:), phi_B(:,:,:,:)
    real(dp) :: R0, rss, dth, x(3), e(3,3), eth(3), eph(3), bg(3), ct, st, ph, w, b, syn, a, &
        st_(5), br, bt, bp, x0, err(3), wt, t0
    integer :: lmx, lcut, it, jb, ib, i, j, k, l, m, k_ss
    logical :: want_phi
    character(len=200) :: msg
    if (pfss_check) then
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        State_GVB(:, :, :, br_:psi_, ib) = 0.0_dp
        call seed_dipole(ib)
      end do
      call cz_unbind()
    end if
    call cz_unbind()
    t0 = arbor_wtime()
    R0 = r_c(k_join)                       ! (unbound view: the global grid)
    rss = b_seed_rss_rsun*Rsun
    k_ss = 0
    if (b_seed == 'dpf') then                          ! the source surface at the nearest radial-block face
      k_ss = k_join + nb(1)
      do i = k_join + nb(1), n_r, nb(1)
        if (abs(r_f(i) - rss) < abs(r_f(k_ss) - rss)) k_ss = i
      end do
      rss = r_f(k_ss)
    end if
    dth = huge(1.0_dp)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      dth = min(dth, minval(Metric_DB(2, ib)%Dx_I(1:nb(2))))
    end do
    dth = arbor_min(dth)
    lmx = pfss_lmax
    if (lmx <= 0) lmx = nint(pi/dth)
    allocate(c(0:lmx, 0:lmx), s(0:lmx, 0:lmx), dc(0:lmx, 0:lmx), ds(0:lmx, 0:lmx), p(0:lmx, 0:lmx), dpt(0:lmx, 0:lmx), &
        Sl(0:lmx), Tl(0:lmx), Ul(0:lmx), fr(0:lmx, nb(1)), gh(0:lmx, nb(1)), frs(0:lmx), ghs(0:lmx), cmp(0:lmx), smp(0:lmx))
    c = 0.0_dp
    s = 0.0_dp
    ! --- the fit: pass it synthesises the current coefficients at the join
    ! shell; the residual's quadrature corrects them (not after the last pass)
    do it = 0, pfss_iter
      dc = 0.0_dp
      ds = 0.0_dp
      st_ = 0.0_dp                          ! sum w, sum w B, sum w B^2, sum w res^2, max |res|
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        if (Key_IB(KI+1, ib)*nb(1) + nb(1) /= k_join) cycle
        do k = 1, nb(3)
          do j = 1, nb(2)
            w = shell_weight(ib, j, k)
            if (w <= 0.0_dp) cycle
            call global_angles(ib, j, k, x, ct, st, ph)
            call pfss_legendre(lmx, ct, st, p)
            do m = 0, lmx
              cmp(m) = cos(m*ph)
              smp(m) = sin(m*ph)
            end do
            b = State_GVB(nb(1), j, k, br_, ib)
            syn = 0.0_dp
            do m = 0, lmx
              do l = m, lmx
                syn = syn + p(l, m)*(c(l, m)*cmp(m) + s(l, m)*smp(m))
              end do
            end do
            st_ = [st_(1) + w, st_(2) + w*b, st_(3) + w*b*b, st_(4) + w*(b - syn)**2, max(st_(5), abs(b - syn))]
            if (it == pfss_iter) cycle
            do m = 0, lmx
              a = merge(1.0_dp, 2.0_dp, m == 0)*w*(b - syn)
              dc(m:lmx, m) = dc(m:lmx, m) + a*cmp(m)*p(m:lmx, m)
              ds(m:lmx, m) = ds(m:lmx, m) + a*smp(m)*p(m:lmx, m)
            end do
          end do
        end do
      end do
      call arbor_sum_inplace(st_(1:4))
      st_(5) = arbor_max(st_(5))
      write(msg, '(A,I0,A,I0,A,F9.5,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A)') '  pfss: l_max ', lmx, ', pass ', it, &
          ': sum of weights/4pi ', st_(1)/(4.0_dp*pi), ', <B_r> ', st_(2)/st_(1), ' G, rms B_r ', sqrt(st_(3)/st_(1)), &
          ' G, rms residual ', sqrt(st_(4)/st_(1)), ', max ', st_(5), ' G'
      call arbor_print(msg)
      if (it == pfss_iter) exit
      call arbor_sum_inplace(dc)
      call arbor_sum_inplace(ds)
      c = c + dc
      s = s + ds
    end do
    c(0, 0) = 0.0_dp                          ! (no monopole)
    ! --- the field above the join
    want_phi = dpf_check .or. b_seed == 'dpf'
    if (want_phi) then
      allocate(phi_B(nb(1), nb(2), nb(3), nBlockMax))
      phi_B = 0.0_dp
    end if
    err = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (r_off + nb(1) < k_join .or. (r_off + nb(1) == k_join .and. .not. want_phi)) cycle
      x0 = R0/min(r_f(0), rss)
      lcut = lmx
      if (x0 < 1.0_dp - 1.0e-12_dp) lcut = min(lmx, max(1, ceiling(log(1.0e-10_dp)/log(x0))))
      do i = 1, nb(1)
        call pfss_radial(lcut, min(r_c(i), rss), R0, rss, fr(0:lcut, i), gh(0:lcut, i))
      end do
      call pfss_radial(lcut, rss, R0, rss, frs(0:lcut), ghs(0:lcut))
      if (pfss_check) Bsave = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:bp_, ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          call global_angles(ib, j, k, x, ct, st, ph)
          call pfss_legendre(lcut, ct, st, p(0:lcut, 0:lcut), dpt(0:lcut, 0:lcut))
          Sl(0:lcut) = 0.0_dp
          Tl(0:lcut) = 0.0_dp
          Ul(0:lcut) = 0.0_dp
          do m = 0, lcut
            cmp(m) = cos(m*ph)
            smp(m) = sin(m*ph)
            do l = m, lcut
              a = c(l, m)*cmp(m) + s(l, m)*smp(m)
              Sl(l) = Sl(l) + p(l, m)*a
              Tl(l) = Tl(l) + dpt(l, m)*a
              Ul(l) = Ul(l) + p(l, m)*m*(s(l, m)*cmp(m) - c(l, m)*smp(m))
            end do
          end do
          Ul(0:lcut) = Ul(0:lcut)/max(st, 1.0e-12_dp)
          eth = [ct*cos(ph), ct*sin(ph), -st]
          eph = [-sin(ph), cos(ph), 0.0_dp]
          e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
          do i = 1, nb(1)
            if (want_phi) then                         ! the potential: Phi = -r sum gh_l S_l
              phi_B(i, j, k, ib) = 0.0_dp
              if (r_c(i) <= rss) phi_B(i, j, k, ib) = -r_c(i)*sum(gh(0:lcut, i)*Sl(0:lcut))
              if (r_off + i <= k_join) cycle           ! (the CZ cells: Phi only)
            end if
            if (r_c(i) <= rss) then
              br = sum(fr(0:lcut, i)*Sl(0:lcut))
              bt = sum(gh(0:lcut, i)*Tl(0:lcut))
              bp = sum(gh(0:lcut, i)*Ul(0:lcut))
            else
              br = sum(frs(0:lcut)*Sl(0:lcut))*(rss/r_c(i))**2
              bt = 0.0_dp
              bp = 0.0_dp
            end if
            bg = br*x + bt*eth + bp*eph                 ! global Cartesian
            if (Key_IB(KP, ib) == 2) bg = yy_cart(bg)  ! to this patch's frame
            State_GVB(i, j, k, br_:bp_, ib) = matmul(transpose(e), bg)
            State_GVB(i, j, k, psi_, ib) = 0.0_dp
            if (pfss_check) then
              wt = sum((State_GVB(i, j, k, br_:bp_, ib) - Bsave(i, j, k, :))**2)
              err = [err(1) + wt, err(2) + sum(Bsave(i, j, k, :)**2), max(err(3), sqrt(wt))]
            end if
          end do
        end do
      end do
    end do
    call cz_unbind()
    if (pfss_check) then
      call arbor_sum_inplace(err(1:2))
      err(3) = arbor_max(err(3))
      write(msg, '(A,ES10.3,A,ES10.3,A)') '  pfss_check: |B - B_dipole| above the join: rms/rms|B| ', &
          sqrt(err(1)/max(err(2), tiny(1.0_dp))), ', max ', err(3)/b_seed_G, ' b_seed_G'
      call arbor_print(msg)
    end if
    write(msg, '(A,ES10.3,A,F8.2,A)') '  pfss: field above the join set (source surface ', rss/Rsun, ' R_sun), ', &
        arbor_wtime() - t0, ' s'
    call arbor_print(msg)
    if (dpf_check) call dpf_check_run(phi_B, rss)
    if (b_seed == 'dpf') then
      call dpf_solve(phi_B, k_ss)
      if (mf_iter > 0) call mf_relax(r_f(k_ss))
      if (pfss_check) then                           ! the solved field against the analytic dipole (below 0.95 r_ss)
        err = 0.0_dp
        do jb = 1, nBlock
          ib = iBlockList_I(jb)
          call cz_bind(ib)
          if (r_off + nb(1) <= k_join) cycle
          Bsave = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:bp_, ib)
          State_GVB(:, :, :, br_:bp_, ib) = 0.0_dp
          call seed_dipole(ib)
          do k = 1, nb(3)
            do j = 1, nb(2)
              do i = 1, nb(1)
                if (r_c(i) > 0.95_dp*rss) cycle
                wt = sum((Bsave(i, j, k, :) - State_GVB(i, j, k, br_:bp_, ib))**2)
                err = [err(1) + wt, err(2) + sum(State_GVB(i, j, k, br_:bp_, ib)**2), max(err(3), sqrt(wt))]
              end do
            end do
          end do
          State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:bp_, ib) = Bsave
        end do
        call cz_unbind()
        call arbor_sum_inplace(err(1:2))
        err(3) = arbor_max(err(3))
        write(msg, '(A,ES10.3,A,ES10.3,A)') '  pfss_check (dpf): |B - B_dipole| below 0.95 r_ss: rms/rms|B| ', &
            sqrt(err(1)/max(err(2), tiny(1.0_dp))), ', max ', err(3)/b_seed_G, ' b_seed_G'
        call arbor_print(msg)
      end if
    end if
  contains
    ! the cell's solid angle, halved where both cores cover the point, 0 outside the core
    real(dp) function shell_weight(ib, j, k) result(w)
      integer, intent(in) :: ib, j, k
      real(dp) :: th, ph_, th2, ph2
      th = Metric_DB(2, ib)%X_I(j)
      ph_ = Metric_DB(3, ib)%X_I(k)
      w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
      if (.not. core(th, ph_)) then
        w = 0.0_dp
      else
        call yy_angles(th, ph_, th2, ph2)
        if (core(th2, ph2)) w = 0.5_dp*w
      end if
    end function shell_weight
    logical function core(t, p_)
      real(dp), intent(in) :: t, p_
      core = t >= 0.25_dp*pi .and. t <= 0.75_dp*pi .and. abs(p_) <= 0.75_dp*pi
    end function core
    ! the cell centre's global unit vector, cos/sin of the global colatitude and the longitude
    subroutine global_angles(ib, j, k, x, ct, st, ph)
      integer, intent(in) :: ib, j, k
      real(dp), intent(out) :: x(3), ct, st, ph
      real(dp) :: e(3,3)
      e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
      x = e(:, 1)
      if (Key_IB(KP, ib) == 2) x = yy_cart(x)
      ct = max(-1.0_dp, min(1.0_dp, x(3)))
      st = sqrt(max(1.0_dp - ct*ct, 0.0_dp))
      ph = atan2(x(2), x(1))
    end subroutine global_angles
  end subroutine cz_seed_pfss

  ! dpf_check: the curl of the analytic PFSS field (now in the state) and of
  ! B = -grad Phi from the code's own stencils (Phi = the PFSS potential,
  ! phi_B), both with filled ghosts, measured by the code's curl as
  ! q = |curl B| r dtheta/|B| (B^2-weighted rms, and the max) in atmosphere
  ! cells whose stencils stay above the join and below 0.95 r_ss, split into
  ! patch-interior and patch-edge cells (within 2.5 cells of the Yin-Yang
  ! patch boundary: their stencils reach interpolated ghosts) and three
  ! heights (z < 5 Mm, 5-50 Mm, > 50 Mm). Also rms(div B dx)/rms|B| of both.
  ! The state's field is restored afterwards (the run keeps the analytic PFSS).
  subroutine dpf_check_run(phi_B, rss)
    real(dp), intent(in) :: phi_B(:,:,:,:), rss
    real(dp), allocatable :: Ban(:,:,:,:,:)
    real(dp) :: acc(4, 2, 3, 2), th, ph, dth, dph, q, b2, z, c2, d2, rj
    integer :: jb, ib, i, j, k, m, ireg, ih
    character(len=200) :: msg
    character(len=8), parameter :: fname(2) = ['PFSS    ', '-grad Ph']
    character(len=8), parameter :: rname(2) = ['interior', 'edge    ']
    character(len=10), parameter :: hname(3) = ['z<5Mm     ', '5-50Mm    ', 'z>50Mm    ']
    call cz_unbind()
    rj = r_f(k_join)                                     ! the join face (unbound view)
    allocate(Ban(nb(1), nb(2), nb(3), 3, nBlockMax))
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      Ban(:, :, :, :, ib) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:bp_, ib)
    end do
    acc = 0.0_dp            ! (sum B^2, sum B^2 q^2, max q, sum (div B dx)^2) x region x height x field
    do m = 1, 2
      if (m == 2) then
        ! Phi into psi (above the join and the top CZ block), then B = -grad Phi above the join
        do jb = 1, nBlock
          ib = iBlockList_I(jb)
          State_GVB(1:nb(1), 1:nb(2), 1:nb(3), psi_, ib) = phi_B(:, :, :, ib)
        end do
        call arbor_fill_ghosts(State_GVB)
        do jb = 1, nBlock
          ib = iBlockList_I(jb)
          call cz_bind(ib)
          if (r_off + nb(1) <= k_join) cycle
          call ddr_block(State_GVB(:, :, :, psi_, ib), Dr_C)
          associate (th_ => Metric_DB(2, ib)%X_I)
            do i = max(1, k_join - r_off + 1), nb(1)
              State_GVB(i, 1:nb(2), 1:nb(3), br_, ib) = -Dr_C(i, :, :)
            end do
            call arbor_deriv_cell(State_GVB, psi_, 2, ib, D_C)
            do i = max(1, k_join - r_off + 1), nb(1)
              State_GVB(i, 1:nb(2), 1:nb(3), bt_, ib) = -D_C(i, :, :)/r_c(i)
            end do
            call arbor_deriv_cell(State_GVB, psi_, 3, ib, D_C)
            do k = 1, nb(3)
              do j = 1, nb(2)
                do i = max(1, k_join - r_off + 1), nb(1)
                  State_GVB(i, j, k, bp_, ib) = -D_C(i, j, k)/(r_c(i)*sin(th_(j)))
                end do
              end do
            end do
          end associate
        end do
        call cz_unbind()
      end if
      call arbor_fill_ghosts(State_GVB)
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        if (r_off + nb(1) <= k_join) cycle
        B_GV = State_GVB(:, :, :, br_:bp_, ib)
        call curl_block(B_GV, ib, J_CV)
        call div_block(State_GVB(:, :, :, br_:bp_, ib), ib, divB_C)
        dth = Metric_DB(2, ib)%Dx_I(1)
        dph = Metric_DB(3, ib)%Dx_I(1)
        do k = 1, nb(3)
          do j = 1, nb(2)
            th = Metric_DB(2, ib)%X_I(j)
            ph = Metric_DB(3, ib)%X_I(k)
            ireg = 1
            if (th < 0.25_dp*pi + 2.5_dp*dth .or. th > 0.75_dp*pi - 2.5_dp*dth .or. &
                abs(ph) > 0.75_dp*pi - 2.5_dp*dph) ireg = 2
            do i = 1, nb(1)
              if (r_off + i < k_join + 3 .or. r_c(i) > 0.95_dp*rss) cycle
              z = r_c(i) - rj
              ih = 1
              if (z > 5.0e8_dp) ih = 2
              if (z > 5.0e9_dp) ih = 3
              b2 = sum(State_GVB(i, j, k, br_:bp_, ib)**2)
              if (b2 <= 0.0_dp) cycle
              c2 = sum(J_CV(i, j, k, :)**2)*(r_c(i)*dth)**2
              d2 = (divB_C(i, j, k)*r_c(i)*dth)**2
              q = sqrt(c2/b2)
              acc(:, ireg, ih, m) = [acc(1, ireg, ih, m) + b2, acc(2, ireg, ih, m) + c2, &
                  max(acc(3, ireg, ih, m), q), acc(4, ireg, ih, m) + d2]
            end do
          end do
        end do
      end do
      call cz_unbind()
    end do
    do m = 1, 2
      do ih = 1, 3
        do ireg = 1, 2
          acc(1, ireg, ih, m) = arbor_sum(acc(1, ireg, ih, m))
          acc(2, ireg, ih, m) = arbor_sum(acc(2, ireg, ih, m))
          acc(3, ireg, ih, m) = arbor_max(acc(3, ireg, ih, m))
          acc(4, ireg, ih, m) = arbor_sum(acc(4, ireg, ih, m))
          write(msg, '(A,A,1X,A,1X,A,A,ES10.3,A,ES10.3,A,ES10.3)') '  dpf_check ', fname(m), rname(ireg), hname(ih), &
              ': rms |curl B| r dth/|B| ', sqrt(acc(2, ireg, ih, m)/max(acc(1, ireg, ih, m), tiny(1.0_dp))), &
              '  max ', acc(3, ireg, ih, m), '  rms div B r dth/|B| ', sqrt(acc(4, ireg, ih, m)/max(acc(1, ireg, ih, m), tiny(1.0_dp)))
          call arbor_print(msg)
        end do
      end do
    end do
    ! the analytic field back (psi 0 above the join, as cz_seed_pfss sets it)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (r_off + nb(1) <= k_join) cycle
      do i = max(1, k_join - r_off + 1), nb(1)
        State_GVB(i, 1:nb(2), 1:nb(3), br_:bp_, ib) = Ban(i, :, :, :, ib)
        State_GVB(i, 1:nb(2), 1:nb(3), psi_, ib) = 0.0_dp
      end do
    end do
    call cz_unbind()
    call arbor_fill_ghosts(State_GVB)
  end subroutine dpf_check_run

  ! b_seed = 'dpf': the discrete potential field above the join. Unknown Phi
  ! in the cells k_join+1..k_ss (whole radial blocks; k_ss a block face, the
  ! source surface). Compact finite-volume Laplacian (V times it, symmetric
  ! within a patch): sum over the faces of c_f (Phi_nb - Phi), c_f = area/
  ! distance; at the join face the flux of the join cell's B_r (r^2 B_r kept)
  ! on the right-hand side, at the k_ss face Phi = 0 (odd ghost). The
  ! neighbours come from ghost fills (Yin-Yang: linear interpolation,
  ! LIMITER_NONE, during the solve). BiCGSTAB, right-preconditioned by a radial
  ! tridiagonal solve per column within each block, started from the PFSS
  ! potential phi_B. Then B = -grad Phi with the code's stencils in the cells
  ! k_join+1..k_ss (Phi below the join extended with the CZ's B_r, odd above
  ! k_ss), and B_r ~ r^-2, B_h = 0 above k_ss; psi = 0 above the join.
  subroutine dpf_solve(phi_B, k_ss)
    real(dp), intent(in) :: phi_B(:,:,:,:)
    integer, intent(in) :: k_ss
    real(dp), allocatable :: X(:,:,:,:,:), Cf(:,:,:,:,:), xv(:,:,:,:), rv(:,:,:,:), rh(:,:,:,:), pv(:,:,:,:), &
        vv(:,:,:,:), sv(:,:,:,:), tv(:,:,:,:), bv(:,:,:,:), ph_(:,:,:,:), sh(:,:,:,:), yc(:,:,:,:)
    logical, allocatable :: edge(:,:,:)
    real(dp) :: rho_, rho_old, alpha, omega, beta, bnorm, rnorm, t0, dO, drr, brf
    integer :: it, jb, ib, i, j, k, lim_save, layer
    character(len=200) :: msg
    t0 = arbor_wtime()
    lim_save = iLimiter
    iLimiter = LIMITER_NONE
    allocate(X(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng, 1, nBlockMax), Cf(nb(1), nb(2), nb(3), 6, nBlockMax))
    allocate(xv(nb(1), nb(2), nb(3), nBlockMax))
    allocate(rv, rh, pv, vv, sv, tv, bv, ph_, sh, yc, mold=xv)
    allocate(edge(nb(2), nb(3), nBlockMax))
    edge = .false.
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      do k = 1, nb(3)
        do j = 1, nb(2)
          associate (th => Metric_DB(2, ib)%X_I(j), ph => Metric_DB(3, ib)%X_I(k), &
              dth => Metric_DB(2, ib)%Dx_I(j), dph => Metric_DB(3, ib)%Dx_I(k))
            edge(j, k, ib) = dpf_operator == 'wide' .and. dpf_edge_cells > 0.0_dp .and. &
                (th < 0.25_dp*pi + dpf_edge_cells*dth .or. th > 0.75_dp*pi - dpf_edge_cells*dth &
                .or. abs(ph) > 0.75_dp*pi - dpf_edge_cells*dph)
          end associate
        end do
      end do
    end do
    X = 0.0_dp;  Cf = 0.0_dp
    xv = 0.0_dp;  rv = 0.0_dp;  rh = 0.0_dp;  pv = 0.0_dp;  vv = 0.0_dp;  sv = 0.0_dp;  tv = 0.0_dp;  bv = 0.0_dp
    ph_ = 0.0_dp;  sh = 0.0_dp
    ! the face coefficients (1: r below, 2: r above, 3: theta below, 4: theta above, 5/6: phi) and b
    call arbor_fill_ghosts(State_GVB)                    ! (the join cell's B_r in the ghosts of the bottom blocks)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (.not. in_solve(ib)) cycle
      associate (th => Metric_DB(2, ib)%X_I, thf => Metric_DB(2, ib)%X_F, dth => Metric_DB(2, ib)%Dx_I, &
          ph => Metric_DB(3, ib)%X_I, dph => Metric_DB(3, ib)%Dx_I, vt => Metric_DB(2, ib)%Vol_I, &
          vp => Metric_DB(3, ib)%Vol_I)
        do k = 1, nb(3)
          do j = 1, nb(2)
            dO = vt(j)*vp(k)
            do i = 1, nb(1)
              drr = r_f(i) - r_f(i-1)
              if (r_off + i > k_join + 1) Cf(i, j, k, 1, ib) = r_f(i-1)**2*dO/(r_c(i) - r_c(i-1))
              Cf(i, j, k, 2, ib) = r_f(i)**2*dO/(r_c(i+1) - r_c(i))
              if (r_off + i == k_ss) Cf(i, j, k, 2, ib) = 2.0_dp*Cf(i, j, k, 2, ib)     ! (Phi = 0 at the face: odd ghost)
              Cf(i, j, k, 3, ib) = drr*sin(thf(j-1))*dph(k)/(th(j) - th(j-1))
              Cf(i, j, k, 4, ib) = drr*sin(thf(j))*dph(k)/(th(j+1) - th(j))
              Cf(i, j, k, 5, ib) = drr*dth(j)/(sin(th(j))*(ph(k) - ph(k-1)))
              Cf(i, j, k, 6, ib) = drr*dth(j)/(sin(th(j))*(ph(k+1) - ph(k)))
              if (r_off + i == k_join + 1) bv(i, j, k, ib) = -dO*r_c(i-1)**2*State_GVB(i-1, j, k, br_, ib)
            end do
          end do
        end do
      end associate
      xv(:, :, :, ib) = phi_B(:, :, :, ib)
    end do
    call cz_unbind()
    if (dpf_operator == 'wide') then                 ! b = V div B(Phi = 0): the CZ's field below the join
      sv = 0.0_dp
      call wide_op(sv, tv, .true.)
      do jb = 1, nBlock                                ! (bv held the compact operator's b)
        ib = iBlockList_I(jb)
        do k = 1, nb(3)
          do j = 1, nb(2)
            if (.not. edge(j, k, ib)) bv(:, j, k, ib) = -tv(:, j, k, ib) + max(dpf_eps, 0.0_dp)*bv(:, j, k, ib)
          end do
        end do
      end do
      tv = 0.0_dp
    end if
    ! BiCGSTAB (right preconditioned): r = b - A x
    call apply_A(xv, vv)
    rv = bv - vv
    rh = rv
    bnorm = sqrt(dot(bv, bv))
    rnorm = sqrt(dot(rv, rv))
    write(msg, '(A,I0,A,ES10.3,A,ES10.3)') '  dpf: source surface at face ', k_ss, ', |b| ', bnorm, &
        ', relative residual of the PFSS start ', rnorm/bnorm
    call arbor_print(msg)
    rho_old = 1.0_dp;  alpha = 1.0_dp;  omega = 1.0_dp
    pv = 0.0_dp;  vv = 0.0_dp
    do it = 1, dpf_maxit
      rho_ = dot(rh, rv)
      beta = (rho_/rho_old)*(alpha/omega)
      pv = rv + beta*(pv - omega*vv)
      call precond(pv, ph_)
      call apply_A(ph_, vv)
      alpha = rho_/dot(rh, vv)
      sv = rv - alpha*vv
      call precond(sv, sh)
      call apply_A(sh, tv)
      omega = dot(tv, sv)/dot(tv, tv)
      xv = xv + alpha*ph_ + omega*sh
      rv = sv - omega*tv
      rho_old = rho_
      rnorm = sqrt(dot(rv, rv))
      if (mod(it, 50) == 0 .or. rnorm <= dpf_tol*bnorm) then
        write(msg, '(A,I6,A,ES10.3)') '  dpf: iteration ', it, ', relative residual ', rnorm/bnorm
        call arbor_print(msg)
      end if
      if (rnorm <= dpf_tol*bnorm) exit
    end do
    ! B = -grad Phi (Phi in psi; its ghosts below the join and above k_ss set here)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (in_solve(ib)) State_GVB(1:nb(1), 1:nb(2), 1:nb(3), psi_, ib) = xv(:, :, :, ib)
    end do
    call arbor_fill_ghosts(State_GVB)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (.not. in_solve(ib)) cycle
      if (r_off == k_join) then                        ! Phi continued with the CZ's B_r (ghosts 0, -1: cells k_join, k_join - 1)
        State_GVB(0, :, :, psi_, ib) = State_GVB(1, :, :, psi_, ib) &
            + State_GVB(0, :, :, br_, ib)*(r_c(0)/r_f(0))**2*(r_c(1) - r_c(0))
        State_GVB(-1, :, :, psi_, ib) = State_GVB(0, :, :, psi_, ib) + State_GVB(0, :, :, br_, ib)*(r_c(0) - r_c(-1))
      end if
      if (r_off + nb(1) == k_ss) then                  ! odd about the k_ss face
        State_GVB(nb(1)+1, :, :, psi_, ib) = -State_GVB(nb(1), :, :, psi_, ib)
        State_GVB(nb(1)+2, :, :, psi_, ib) = -State_GVB(nb(1)-1, :, :, psi_, ib)
      end if
      call ddr_block(State_GVB(:, :, :, psi_, ib), Dr_C)
      State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_, ib) = -Dr_C
      call arbor_deriv_cell(State_GVB, psi_, 2, ib, D_C)
      do i = 1, nb(1)
        State_GVB(i, 1:nb(2), 1:nb(3), bt_, ib) = -D_C(i, :, :)/r_c(i)
      end do
      call arbor_deriv_cell(State_GVB, psi_, 3, ib, D_C)
      associate (th => Metric_DB(2, ib)%X_I)
        do j = 1, nb(2)
          do i = 1, nb(1)
            State_GVB(i, j, 1:nb(3), bp_, ib) = -D_C(i, j, :)/(r_c(i)*sin(th(j)))
          end do
        end do
      end associate
    end do
    call cz_unbind()
    ! above k_ss: radial, block layer by layer (each from the ghosts of the one below)
    do layer = 0, (n_r - k_ss)/nb(1) - 1
      call arbor_fill_ghosts(State_GVB)
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        if (r_off /= k_ss + layer*nb(1)) cycle
        do k = 1, nb(3)
          do j = 1, nb(2)
            if (layer == 0) then                       ! B_r at the k_ss face from Phi (0 there)
              brf = 2.0_dp*State_GVB(0, j, k, psi_, ib)/(r_c(1) - r_c(0))*(r_f(0)/r_c(1))**0
              brf = brf*r_f(0)**2
            else
              brf = State_GVB(0, j, k, br_, ib)*r_c(0)**2
            end if
            do i = 1, nb(1)
              State_GVB(i, j, k, br_, ib) = brf/r_c(i)**2
            end do
            State_GVB(1:nb(1), j, k, bt_:bp_, ib) = 0.0_dp
          end do
        end do
      end do
      call cz_unbind()
    end do
    ! psi 0 above the join
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (r_off + nb(1) <= k_join) cycle
      State_GVB(1:nb(1), 1:nb(2), 1:nb(3), psi_, ib) = 0.0_dp
    end do
    call cz_unbind()
    call arbor_fill_ghosts(State_GVB)
    iLimiter = lim_save
    write(msg, '(A,I0,A,ES10.3,A,F8.2,A)') '  dpf: ', min(it, dpf_maxit), ' iterations, relative residual ', &
        rnorm/bnorm, ', ', arbor_wtime() - t0, ' s'
    call arbor_print(msg)
    call b_field_stats('dpf', r_f(k_ss))
    call arbor_fill_ghosts(State_GVB)
  contains
    logical function in_solve(ib)
      integer, intent(in) :: ib
      in_solve = Key_IB(KI+1, ib)*nb(1) >= k_join .and. Key_IB(KI+1, ib)*nb(1) + nb(1) <= k_ss
    end function in_solve
    real(dp) function dot(a, b)
      real(dp), intent(in) :: a(:,:,:,:), b(:,:,:,:)
      integer :: jb_, ib_
      dot = 0.0_dp
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        if (in_solve(ib_)) dot = dot + sum(a(:, :, :, ib_)*b(:, :, :, ib_))
      end do
      dot = arbor_sum(dot)
    end function dot
    ! y = A x (ghosts of x filled first)
    subroutine apply_A(xin, y)
      real(dp), intent(in) :: xin(:,:,:,:)
      real(dp), intent(out) :: y(:,:,:,:)
      integer :: jb_, ib_, i_, j_, k_
      if (dpf_operator == 'wide') call wide_op(xin, y, .false.)
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        X(1:nb(1), 1:nb(2), 1:nb(3), 1, ib_) = xin(:, :, :, ib_)
      end do
      call arbor_fill_ghosts(X)
      yc = 0.0_dp
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        if (.not. in_solve(ib_)) cycle
        do k_ = 1, nb(3)
          do j_ = 1, nb(2)
            do i_ = 1, nb(1)
              associate (c => Cf(i_, j_, k_, :, ib_), x0 => X(i_, j_, k_, 1, ib_))
                yc(i_, j_, k_, ib_) = c(1)*(X(i_-1, j_, k_, 1, ib_) - x0) &
                    + merge(-c(2)*x0, c(2)*(X(i_+1, j_, k_, 1, ib_) - x0), Key_IB(KI+1, ib_)*nb(1) + i_ == k_ss) &
                    + c(3)*(X(i_, j_-1, k_, 1, ib_) - x0) + c(4)*(X(i_, j_+1, k_, 1, ib_) - x0) &
                    + c(5)*(X(i_, j_, k_-1, 1, ib_) - x0) + c(6)*(X(i_, j_, k_+1, 1, ib_) - x0)
              end associate
            end do
          end do
        end do
      end do
      if (dpf_operator == 'compact') then
        y = yc
      else                                           ! wide (+ eps compact); compact rows in the edge band
        do jb_ = 1, nBlock
          ib_ = iBlockList_I(jb_)
          do k_ = 1, nb(3)
            do j_ = 1, nb(2)
              if (edge(j_, k_, ib_)) then
                y(:, j_, k_, ib_) = yc(:, j_, k_, ib_)
              else
                y(:, j_, k_, ib_) = y(:, j_, k_, ib_) + dpf_eps*yc(:, j_, k_, ib_)
              end if
            end do
          end do
        end do
      end if
    end subroutine apply_A
    ! The wide operator: B = -grad x with the code's stencils in the solve
    ! cells (x in psi; below the join x continued, with the CZ's B_r if
    ! full, flat if not; odd above k_ss), its ghosts filled (the CZ's B below
    ! the join if full, 0 if not; above k_ss B_r even, B_h odd), then
    ! y = -V div B (V the cell volume): y = V div(grad x) for the homogeneous
    ! part. full: x = 0 gives -V div B of the CZ's field alone (-b)
    subroutine wide_op(xin, y, full)
      real(dp), intent(in) :: xin(:,:,:,:)
      real(dp), intent(out) :: y(:,:,:,:)
      logical, intent(in) :: full
      integer :: jb_, ib_, i_, j_, k_, m_
      real(dp) :: vol
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        if (in_solve(ib_)) then
          State_GVB(1:nb(1), 1:nb(2), 1:nb(3), psi_, ib_) = xin(:, :, :, ib_)
        else if (Key_IB(KI+1, ib_)*nb(1) >= k_join) then
          State_GVB(1:nb(1), 1:nb(2), 1:nb(3), psi_, ib_) = 0.0_dp
        end if
      end do
      call arbor_fill_ghosts(State_GVB, iVarMin=psi_, iVarMax=psi_)
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        call cz_bind(ib_)
        if (.not. in_solve(ib_)) cycle
        if (r_off == k_join) then
          if (full) then
            State_GVB(0, :, :, psi_, ib_) = State_GVB(1, :, :, psi_, ib_) &
                + State_GVB(0, :, :, br_, ib_)*(r_c(0)/r_f(0))**2*(r_c(1) - r_c(0))
            State_GVB(-1, :, :, psi_, ib_) = State_GVB(0, :, :, psi_, ib_) + State_GVB(0, :, :, br_, ib_)*(r_c(0) - r_c(-1))
          else
            State_GVB(0, :, :, psi_, ib_) = State_GVB(1, :, :, psi_, ib_)
            State_GVB(-1, :, :, psi_, ib_) = State_GVB(1, :, :, psi_, ib_)
          end if
        end if
        if (r_off + nb(1) == k_ss) then
          State_GVB(nb(1)+1, :, :, psi_, ib_) = -State_GVB(nb(1), :, :, psi_, ib_)
          State_GVB(nb(1)+2, :, :, psi_, ib_) = -State_GVB(nb(1)-1, :, :, psi_, ib_)
        end if
        call ddr_block(State_GVB(:, :, :, psi_, ib_), Dr_C)
        State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_, ib_) = -Dr_C
        call arbor_deriv_cell(State_GVB, psi_, 2, ib_, D_C)
        do i_ = 1, nb(1)
          State_GVB(i_, 1:nb(2), 1:nb(3), bt_, ib_) = -D_C(i_, :, :)/r_c(i_)
        end do
        call arbor_deriv_cell(State_GVB, psi_, 3, ib_, D_C)
        associate (th => Metric_DB(2, ib_)%X_I)
          do j_ = 1, nb(2)
            do i_ = 1, nb(1)
              State_GVB(i_, j_, 1:nb(3), bp_, ib_) = -D_C(i_, j_, :)/(r_c(i_)*sin(th(j_)))
            end do
          end do
        end associate
      end do
      call cz_unbind()
      call arbor_fill_ghosts(State_GVB, iVarMin=br_, iVarMax=bp_)
      y = 0.0_dp
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        call cz_bind(ib_)
        if (.not. in_solve(ib_)) cycle
        if (r_off == k_join .and. .not. full) State_GVB(-1:0, :, :, br_:bp_, ib_) = 0.0_dp
        if (r_off + nb(1) == k_ss) then
          do m_ = 1, 2
            State_GVB(nb(1)+m_, :, :, br_, ib_) = State_GVB(nb(1)+1-m_, :, :, br_, ib_)
            State_GVB(nb(1)+m_, :, :, bt_:bp_, ib_) = -State_GVB(nb(1)+1-m_, :, :, bt_:bp_, ib_)
          end do
        end if
        call div_block(State_GVB(:, :, :, br_:bp_, ib_), ib_, divB_C)
        associate (vt => Metric_DB(2, ib_)%Vol_I, vp => Metric_DB(3, ib_)%Vol_I)
          do k_ = 1, nb(3)
            do j_ = 1, nb(2)
              do i_ = 1, nb(1)
                vol = (r_f(i_)**3 - r_f(i_-1)**3)/3.0_dp*vt(j_)*vp(k_)
                y(i_, j_, k_, ib_) = -vol*divB_C(i_, j_, k_)
              end do
            end do
          end do
        end associate
      end do
      call cz_unbind()
    end subroutine wide_op
    ! z = M^-1 r: per column of each block, the radial tridiagonal part of A
    ! (the couplings out of the block dropped, the diagonal kept whole)
    subroutine precond(rin, z)
      real(dp), intent(in) :: rin(:,:,:,:)
      real(dp), intent(out) :: z(:,:,:,:)
      real(dp) :: dd(nb(1)), lo(nb(1)), up(nb(1)), cp(nb(1)), dp_(nb(1)), m_
      integer :: jb_, ib_, i_, j_, k_
      z = 0.0_dp
      do jb_ = 1, nBlock
        ib_ = iBlockList_I(jb_)
        if (.not. in_solve(ib_)) cycle
        do k_ = 1, nb(3)
          do j_ = 1, nb(2)
            do i_ = 1, nb(1)
              associate (c => Cf(i_, j_, k_, :, ib_))
                dd(i_) = -(c(1) + c(2) + c(3) + c(4) + c(5) + c(6))
                lo(i_) = merge(c(1), 0.0_dp, i_ > 1)
                up(i_) = merge(c(2), 0.0_dp, i_ < nb(1))
              end associate
            end do
            ! Thomas
            cp(1) = up(1)/dd(1)
            dp_(1) = rin(1, j_, k_, ib_)/dd(1)
            do i_ = 2, nb(1)
              m_ = dd(i_) - lo(i_)*cp(i_-1)
              cp(i_) = up(i_)/m_
              dp_(i_) = (rin(i_, j_, k_, ib_) - lo(i_)*dp_(i_-1))/m_
            end do
            z(nb(1), j_, k_, ib_) = dp_(nb(1))
            do i_ = nb(1) - 1, 1, -1
              z(i_, j_, k_, ib_) = dp_(i_) - cp(i_)*z(i_+1, j_, k_, ib_)
            end do
          end do
        end do
      end do
    end subroutine precond
  end subroutine dpf_solve

  ! The field's curl and divergence (the code's stencils, ghosts filled) as
  ! q = |curl B| r dtheta/|B| and |div B| r dtheta/|B| (B^2-weighted rms;
  ! max of q) in the atmosphere cells whose stencils stay above the join and
  ! below 0.95 rss, patch interior / patch edge (within 2.5 cells of the
  ! Yin-Yang boundary), z < 5 Mm, 5-50 Mm, > 50 Mm
  subroutine b_field_stats(label, rss)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: rss
    real(dp) :: acc(4, 2, 3), th, ph, dth, dph, q, b2, z, c2, d2, rj
    integer :: jb, ib, i, j, k, ireg, ih
    character(len=200) :: msg
    character(len=8), parameter :: rname(2) = ['interior', 'edge    ']
    character(len=10), parameter :: hname(3) = ['z<5Mm     ', '5-50Mm    ', 'z>50Mm    ']
    call cz_unbind()
    rj = r_f(k_join)
    call arbor_fill_ghosts(State_GVB)
    acc = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (r_off + nb(1) <= k_join) cycle
      B_GV = State_GVB(:, :, :, br_:bp_, ib)
      call curl_block(B_GV, ib, J_CV)
      call div_block(State_GVB(:, :, :, br_:bp_, ib), ib, divB_C)
      dth = Metric_DB(2, ib)%Dx_I(1)
      dph = Metric_DB(3, ib)%Dx_I(1)
      do k = 1, nb(3)
        do j = 1, nb(2)
          th = Metric_DB(2, ib)%X_I(j)
          ph = Metric_DB(3, ib)%X_I(k)
          ireg = 1
          if (th < 0.25_dp*pi + 2.5_dp*dth .or. th > 0.75_dp*pi - 2.5_dp*dth .or. &
              abs(ph) > 0.75_dp*pi - 2.5_dp*dph) ireg = 2
          do i = 1, nb(1)
            if (r_off + i < k_join + 3 .or. r_c(i) > 0.95_dp*rss) cycle
            z = r_c(i) - rj
            ih = 1
            if (z > 5.0e8_dp) ih = 2
            if (z > 5.0e9_dp) ih = 3
            b2 = sum(State_GVB(i, j, k, br_:bp_, ib)**2)
            if (b2 <= 0.0_dp) cycle
            c2 = sum(J_CV(i, j, k, :)**2)*(r_c(i)*dth)**2
            d2 = (divB_C(i, j, k)*r_c(i)*dth)**2
            q = sqrt(c2/b2)
            acc(:, ireg, ih) = [acc(1, ireg, ih) + b2, acc(2, ireg, ih) + c2, max(acc(3, ireg, ih), q), acc(4, ireg, ih) + d2]
          end do
        end do
      end do
    end do
    call cz_unbind()
    do ih = 1, 3
      do ireg = 1, 2
        acc(1, ireg, ih) = arbor_sum(acc(1, ireg, ih))
        acc(2, ireg, ih) = arbor_sum(acc(2, ireg, ih))
        acc(3, ireg, ih) = arbor_max(acc(3, ireg, ih))
        acc(4, ireg, ih) = arbor_sum(acc(4, ireg, ih))
        write(msg, '(A,A,1X,A,1X,A,A,ES10.3,A,ES10.3,A,ES10.3)') '  field ', label, rname(ireg), hname(ih), &
            ': rms |curl B| r dth/|B| ', sqrt(acc(2, ireg, ih)/max(acc(1, ireg, ih), tiny(1.0_dp))), &
            '  max ', acc(3, ireg, ih), '  rms |div B| r dth/|B| ', sqrt(acc(4, ireg, ih)/max(acc(1, ireg, ih), tiny(1.0_dp)))
        call arbor_print(msg)
      end do
    end do
  end subroutine b_field_stats

  ! Magnetofriction (mf_iter > 0, after b_seed = 'dpf'; see cz3d_params_mod):
  ! the field in the atmosphere cells below rss moved towards zero J x B as
  ! the code computes it, the gas frozen. The displacement d is stored in the
  ! velocity slots (0 elsewhere) for the ghost fill; v and psi restored after.
  subroutine mf_relax(rss)
    real(dp), intent(in) :: rss
    real(dp), allocatable :: Vsave(:,:,:,:,:)
    real(dp) :: b(3), jj(3), f(3), h(3), b2, b2m, dmax, fsum(4), fe, t0
    integer :: it, jb, ib, i, j, k
    logical :: edg
    character(len=200) :: msg
    t0 = arbor_wtime()
    call cz_unbind()
    allocate(Vsave(nb(1), nb(2), nb(3), 3, nBlockMax))
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      Vsave(:, :, :, :, ib) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), vr_:vp_, ib)
    end do
    call b_field_stats('before mf', rss)
    do it = 1, mf_iter
      call arbor_fill_ghosts(State_GVB)
      dmax = 0.0_dp
      fsum = 0.0_dp                ! sum B^2 and sum (|F| h/B^2)^2 B^2, interior / edge
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        State_GVB(:, :, :, vr_:vp_, ib) = 0.0_dp
        if (r_off + nb(1) <= k_join) cycle
        B_GV = State_GVB(:, :, :, br_:bp_, ib)
        call curl_block(B_GV, ib, J_CV)
        associate (th => Metric_DB(2, ib)%X_I, dth => Metric_DB(2, ib)%Dx_I, ph => Metric_DB(3, ib)%X_I, &
            dph => Metric_DB(3, ib)%Dx_I)
          b2m = 0.0_dp
          do k = 1, nb(3)
            do j = 1, nb(2)
              edg = th(j) < 0.25_dp*pi + 2.5_dp*dth(j) .or. th(j) > 0.75_dp*pi - 2.5_dp*dth(j) &
                  .or. abs(ph(k)) > 0.75_dp*pi - 2.5_dp*dph(k)
              do i = 1, nb(1)
                if (r_off + i <= k_join .or. r_c(i) > rss) cycle
                b = State_GVB(i, j, k, br_:bp_, ib)
                jj = J_CV(i, j, k, :)
                b2 = sum(b**2) + 1.0e-6_dp
                f = [jj(2)*b(3) - jj(3)*b(2), jj(3)*b(1) - jj(1)*b(3), jj(1)*b(2) - jj(2)*b(1)]
                h = [r_f(i) - r_f(i-1), r_c(i)*dth(j), r_c(i)*sin(th(j))*dph(k)]
                State_GVB(i, j, k, vr_:vp_, ib) = mf_alpha*minval(h)**2*f/b2   ! (isotropic: the step limited by the smallest scale)
                dmax = max(dmax, maxval(abs(State_GVB(i, j, k, vr_:vp_, ib)))/minval(h))
                fe = sum(f**2)*minval(h)**2/b2         ! (|F| h_min/B^2)^2 B^2
                if (edg) then
                  fsum(3:4) = fsum(3:4) + [b2, fe]
                else
                  fsum(1:2) = fsum(1:2) + [b2, fe]
                end if
              end do
            end do
          end do
        end associate
      end do
      call cz_unbind()
      call arbor_fill_ghosts(State_GVB, iVarMin=vr_, iVarMax=vp_)
      ! B += curl(d x B) in the relaxed cells
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        if (r_off + nb(1) <= k_join) cycle
        do k = 1 - ng, nb(3) + ng
          do j = 1 - ng, nb(2) + ng
            do i = 1 - ng, nb(1) + ng
              jj = State_GVB(i, j, k, vr_:vp_, ib)
              b = State_GVB(i, j, k, br_:bp_, ib)
              B_GV(i, j, k, :) = [jj(2)*b(3) - jj(3)*b(2), jj(3)*b(1) - jj(1)*b(3), jj(1)*b(2) - jj(2)*b(1)]
            end do
          end do
        end do
        call curl_block(B_GV, ib, C_CV)
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, nb(1)
              if (r_off + i <= k_join .or. r_c(i) > rss) cycle
              State_GVB(i, j, k, br_:bp_, ib) = State_GVB(i, j, k, br_:bp_, ib) + C_CV(i, j, k, :)
            end do
          end do
        end do
      end do
      call cz_unbind()
      if (it == 1 .or. mod(it, 25) == 0 .or. it == mf_iter) then
        fsum(1) = arbor_sum(fsum(1));  fsum(2) = arbor_sum(fsum(2))
        fsum(3) = arbor_sum(fsum(3));  fsum(4) = arbor_sum(fsum(4))
        dmax = arbor_max(dmax)
        write(msg, '(A,I5,A,ES10.3,A,ES10.3,A,ES10.3)') '  mf: iteration ', it, ': rms |J x B| h/B^2 interior ', &
            sqrt(fsum(2)/max(fsum(1), tiny(1.0_dp))), ', edge ', sqrt(fsum(4)/max(fsum(3), tiny(1.0_dp))), &
            ', max step/h ', dmax
        call arbor_print(msg)
      end if
    end do
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      State_GVB(1:nb(1), 1:nb(2), 1:nb(3), vr_:vp_, ib) = Vsave(:, :, :, :, ib)
    end do
    call arbor_fill_ghosts(State_GVB)
    write(msg, '(A,I0,A,F8.2,A)') '  mf: ', mf_iter, ' iterations, ', arbor_wtime() - t0, ' s'
    call arbor_print(msg)
    call b_field_stats('after mf', rss)
  end subroutine mf_relax

  ! Every local block's interior state (all nv variables), one file per rank:
  ! out_dir/state_export_rNNNN.bin, stream access: n_block, nb(1:3), nv; then
  ! per block its key (KEY_LEN integers) and State(1:nb(1), 1:nb(2), 1:nb(3), 1:nv).
  ! The root also writes out_dir/state_export_meta.bin: t, the step, and the
  ! restart numbers (profile counter and schedule, the shell clocks: as
  ! cz_restart_save), so an import continues the run.
  module subroutine export_state_write(t, iStep)
    real(dp), intent(in) :: t
    integer, intent(in) :: iStep
    integer :: u, jb, ib
    character(len=300) :: name
    real(dp), allocatable :: Data_I(:)
    if (.not. export_state) return
    write(name, '(A,A,I4.4,A)') trim(out_dir), '/state_export_r', my_rank, '.bin'
    open(newunit=u, file=trim(name), access='stream', form='unformatted', status='replace', action='write')
    write(u) nBlock, nb, nv
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      write(u) Key_IB(:, ib), State_GVB(1:nb(1), 1:nb(2), 1:nb(3), 1:nv, ib)
    end do
    close(u)
    if (.not. is_root) return
    Data_I = [real(i_prof, dp), t_next_prof, t_prof0, tloc]
    open(newunit=u, file=trim(out_dir)//'/state_export_meta.bin', access='stream', form='unformatted', &
        status='replace', action='write')
    write(u) t, iStep, size(Data_I), Data_I
    close(u)
  end subroutine export_state_write

  ! Index of the export files in import_dir (every rank: the keys only, read
  ! by seeking): file, position and key of every exported block
  module subroutine import_open()
    integer :: u, n, nbf(3), nvf, f, k, m
    integer(i8) :: rec, pos
    character(len=300) :: name
    logical :: exists
    integer, allocatable :: key_K(:,:), file_K(:)
    integer(i8), allocatable :: pos_K(:)
    if (allocated(imp_key)) return
    m = 0
    allocate(key_K(KEY_LEN, 0), file_K(0), pos_K(0))
    do f = 0, 99999
      write(name, '(A,A,I4.4,A)') trim(import_dir), '/state_export_r', f, '.bin'
      inquire(file=trim(name), exist=exists)
      if (.not. exists) exit
      open(newunit=u, file=trim(name), access='stream', form='unformatted', status='old', action='read')
      read(u) n, nbf, nvf
      if (f == 0) then
        imp_nb = nbf
        imp_nv = nvf
      else if (any(nbf /= imp_nb) .or. nvf /= imp_nv) then
        call arbor_abort('cz3d: import: the export files differ in block size or variables')
      end if
      rec = 4_i8*KEY_LEN + 8_i8*product(int(nbf, i8))*nvf
      key_K = reshape([key_K, [(0, k = 1, KEY_LEN*n)]], [KEY_LEN, m + n])
      file_K = [file_K, [(f, k = 1, n)]]
      pos_K = [pos_K, [(0_i8, k = 1, n)]]
      do k = 1, n
        pos = 21_i8 + int(k - 1, i8)*rec                ! (stream positions count from 1; the header is 5 integers)
        read(u, pos=pos) key_K(:, m + k)
        pos_K(m + k) = pos
      end do
      m = m + n
      close(u)
    end do
    if (m == 0) call arbor_abort('cz3d: import: no state_export files in '//trim(import_dir))
    if (any(imp_nb(2:3) /= nb(2:3))) call arbor_abort('cz3d: import: the horizontal block size differs')
    ! the export's top (cells): a taller run (the cascade) starts its blocks
    ! above it from the init file (import_block)
    imp_top = maxval(key_K(KI+1, 1:m) + 1)*imp_nb(1)
    if (mod(imp_nb(1), nb(1)) /= 0) call arbor_abort('cz3d: import: the radial block size must divide the exported one')
    call move_alloc(key_K, imp_key)
    call move_alloc(file_K, imp_file)
    call move_alloc(pos_K, imp_pos)
  end subroutine import_open

  ! The finest horizontal level of the exported blocks
  module integer function import_max_level() result(l)
    l = maxval(imp_key(KL+2, :))
  end function import_max_level

  ! Refinement criterion for the import tree (arbor_tree_refine): refine a
  ! block until an exported block has its horizontal position and level
  module subroutine import_criterion(key, x_lo, x_hi, flags)
    integer, intent(in) :: key(KEY_LEN)
    real(dp), intent(in) :: x_lo(3), x_hi(3)
    integer, intent(out) :: flags(3)
    flags = 0
    if (import_record(key, .false.) == 0) flags(2:3) = 1
  end subroutine import_criterion

  ! The exported block holding block key (horizontal position and level
  ! equal; with radial: also the radial block that contains it), 0 if none
  integer function import_record(key, radial) result(r)
    integer, intent(in) :: key(KEY_LEN)
    logical, intent(in) :: radial
    integer :: m
    r = 0
    do m = 1, size(imp_key, 2)
      if (imp_key(KP, m) /= key(KP)) cycle
      if (any(imp_key(KL+2:KL+3, m) /= key(KL+2:KL+3))) cycle
      if (any(imp_key(KI+2:KI+3, m) /= key(KI+2:KI+3))) cycle
      if (radial) then
        if (imp_key(KI+1, m) /= (key(KI+1)*nb(1))/imp_nb(1)) cycle
      end if
      r = m
      return
    end do
  end function import_record

  ! The block's state from the export files in import_dir: its radial part of
  ! the exported block at the same place (the export may hold whole columns
  ! or larger radial blocks); the variables both runs have are copied, the
  ! others stay as set (B, psi: 0). With the atmosphere, a block above the
  ! join without an exported block keeps its initial state (a CZ-only export,
  ! e.g. CZ_3D's, under the atmosphere of the init file)
  subroutine import_block(iBlock)
    integer, intent(in) :: iBlock
    real(dp), allocatable :: buf(:,:,:,:)
    integer :: m, u, i0
    character(len=300) :: name
    call import_open()
    m = import_record(Key_IB(:, iBlock), .true.)
    if (m == 0 .and. use_atmosphere .and. Key_IB(KI+1, iBlock)*nb(1) >= k_join) return
    ! the cascade: a block above the export's top keeps the init file's state
    ! (B = 0 with no seed: with a conductor top below, B_r ~ 0 at the old top)
    if (m == 0 .and. Key_IB(KI+1, iBlock)*nb(1) >= imp_top) return
    if (m == 0) call arbor_abort('cz3d: import: no exported block at this place (a different grid?)')
    allocate(buf(imp_nb(1), imp_nb(2), imp_nb(3), imp_nv))
    write(name, '(A,A,I4.4,A)') trim(import_dir), '/state_export_r', imp_file(m), '.bin'
    open(newunit=u, file=trim(name), access='stream', form='unformatted', status='old', action='read')
    read(u, pos=imp_pos(m) + 4_i8*KEY_LEN) buf
    close(u)
    i0 = Key_IB(KI+1, iBlock)*nb(1) - imp_key(KI+1, m)*imp_nb(1)
    State_GVB(1:nb(1), 1:nb(2), 1:nb(3), 1:min(nv, imp_nv), iBlock) = buf(i0+1:i0+nb(1), :, :, 1:min(nv, imp_nv))
  end subroutine import_block

  ! import_cz_on_restart: the CZ blocks (below the join) from the export
  module subroutine cz_import_cz()
    integer :: jb, ib, n
    character(len=200) :: msg
    n = 0
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (Key_IB(KI+1, ib)*nb(1) + nb(1) > k_join) cycle
      call import_block(ib)
      n = n + 1
    end do
    n = arbor_sum(n)
    write(msg, '(A,I0,A,A)') '  import_cz_on_restart: ', n, ' CZ blocks from ', trim(import_dir)
    call arbor_print(msg)
  end subroutine cz_import_cz

  ! atm_match_base: each column's lower atmosphere pressure-matched to its
  ! own top CZ cell (see cz3d_params_mod; after an import, before the field
  ! seed). f per column from the blocks holding cell k_join, gathered in a
  ! (patch, theta, phi) map at the run's horizontal level; the atmosphere
  ! cells below atm_match_z_km(2) then get init_file's density times
  ! f_eff(z) at init_file's T (s from the EOS), v unchanged
  module subroutine cz_match_base()
    real(dp), allocatable :: fm(:,:,:)
    real(dp) :: out(n_field), lnp_ref, rho, rp, Tp, sp, rn, sn, zk, w, fe, fmin, fmax, qmin
    integer :: nth, nph, jb, ib, i, j, k, jg, kg, ig, n
    character(len=200) :: msg
    nth = n_theta*2**refine_h
    nph = n_phi*2**refine_h
    allocate(fm(2, nth, nph))
    fm = 0.0_dp
    call cz_unbind()
    call eos_eval(tab, rho0(k_join) + rho1_init(k_join), s0(k_join) + s1_init(k_join), out)
    lnp_ref = out(F_LNP)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (Key_IB(KL+2, ib) /= refine_h .or. Key_IB(KL+3, ib) /= refine_h) &
          call arbor_abort('cz3d: atm_match_base needs every block at the horizontal level refine_h')
      i = k_join - Key_IB(KI+1, ib)*nb(1)
      if (i < 1 .or. i > nb(1)) cycle
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          rho = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
          call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
          fm(Key_IB(KP, ib), Key_IB(KI+2, ib)*nb(2) + j, Key_IB(KI+3, ib)*nb(3) + k) = exp(out(F_LNP) - lnp_ref)
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(fm)
    fmin = minval(fm);  fmax = maxval(fm)
    if (fmin <= 0.0_dp) call arbor_abort('cz3d: atm_match_base: a column without its top CZ cell')
    qmin = huge(1.0_dp);  n = 0
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do i = 1, nb(1)
        ig = Key_IB(KI+1, ib)*nb(1) + i
        zk = z_c(i)/1.0e5_dp
        if (ig <= k_join .or. zk >= atm_match_z_km(2)) cycle
        w = 1.0_dp
        if (zk > atm_match_z_km(1)) w = cos(0.5_dp*pi*(zk - atm_match_z_km(1))/(atm_match_z_km(2) - atm_match_z_km(1)))**2
        rp = rho0(i) + rho1_init(i)
        sp = s0(i) + s1_init(i)
        call eos_eval(tab, rp, sp, out)
        Tp = exp(out(F_LNT))
        do k = 1, nb(3)
          do j = 1, nb(2)
            jg = Key_IB(KI+2, ib)*nb(2) + j
            kg = Key_IB(KI+3, ib)*nb(3) + k
            fe = 1.0_dp + (fm(Key_IB(KP, ib), jg, kg) - 1.0_dp)*w
            rn = fe*rp
            sn = eos_s_from_rho_T(tab, rn, Tp, s_guess=sp)
            State_GVB(i, j, k, rho1_, ib) = (rn - rho0(i))/rq(i)
            State_GVB(i, j, k, s1_, ib) = sn - s0(i)
            qmin = min(qmin, rn/rho0(i))
            n = n + 1
          end do
        end do
      end do
    end do
    call cz_unbind()
    qmin = arbor_min(qmin);  n = arbor_sum(n)
    write(msg, '(A,F7.4,A,F7.4,A,I0,A,F6.0,A,F6.0,A,F7.4)') '  atm_match_base: f = p(top CZ cell)/p(init_file) in ', fmin, &
        ' .. ', fmax, '; ', n, ' atmosphere cells rescaled (full below ', atm_match_z_km(1), ' km, none above ', &
        atm_match_z_km(2), ' km); min rho/rho0 there ', qmin
    call arbor_print(msg)
    if (qmin < 0.9_dp) call arbor_print('  atm_match_base: WARNING rho < 0.9 rho0 in the lower atmosphere (the background rule rho0 <= rho)')
  end subroutine cz_match_base

  ! The time, step and restart numbers of the export (state_export_meta.bin),
  ! if there: t and iStep set, the profile counter and the shell clocks
  ! restored; found false (nothing changed) without the file
  module subroutine import_meta(t, iStep, found)
    real(dp), intent(inout) :: t
    integer, intent(inout) :: iStep
    logical, intent(out) :: found
    integer :: u, nd
    real(dp), allocatable :: Data_I(:)
    inquire(file=trim(import_dir)//'/state_export_meta.bin', exist=found)
    if (.not. found) return
    open(newunit=u, file=trim(import_dir)//'/state_export_meta.bin', access='stream', form='unformatted', &
        status='old', action='read')
    read(u) t, iStep, nd
    allocate(Data_I(nd))
    read(u) Data_I
    close(u)
    i_prof = nint(Data_I(1))
    t_next_prof = Data_I(2)
    t_prof0 = Data_I(3)
    if (nd == 3 + n_r) then
      tloc = Data_I(4:)
    else
      tloc = t                                            ! (a different radial grid: all clocks at the run's time)
    end if
  end subroutine import_meta

  ! Uniform in [-1, 1] from (patch, i, j, k): Park-Miller steps (no overflow)
  real(dp) function urand(p, i, j, k)
    integer, intent(in) :: p, i, j, k
    integer(i8), parameter :: a = 48271_i8, m = 2147483647_i8
    integer(i8) :: x
    integer :: n
    x = 12345_i8 + 999331_i8*int(p, i8)
    do n = 1, 3
      x = mod(a*mod(x + int(i, i8), m), m)
      x = mod(a*mod(x + 7919_i8*int(j, i8), m), m)
      x = mod(a*mod(x + 104729_i8*int(k, i8), m), m)
    end do
    urand = 2.0_dp*real(x, dp)/real(m, dp) - 1.0_dp
  end function urand

  ! si_radial: the columns (patch, theta block, phi block) to ranks, each
  ! column's radial blocks on one rank (the column solves need no messages).
  ! The columns in the order their first blocks appear on the tree's Hilbert
  ! curve, cut into n_ranks contiguous stretches; the blocks move with their
  ! data (arbor_redistribute). All blocks at one horizontal level.
  module subroutine cz_columns_to_ranks(tree)
    type(arbor_tree_t), intent(in) :: tree
    integer, allocatable :: cid(:), pos(:), own(:), lev(:), per(:)
    integer :: j, nth_b, nph_b, ncol, n
    character(len=200) :: msg
    associate (key => tree%key(:, 1:tree%n_leaf))
    if (any(key(KL+2, :) /= key(KL+2, 1)) .or. any(key(KL+3, :) /= key(KL+2, 1))) &
        call arbor_abort('cz3d: si_radial needs every block at one horizontal level')
    nth_b = maxval(key(KI+2, :)) + 1
    nph_b = maxval(key(KI+3, :)) + 1
    ncol = maxval(key(KP, :))*nth_b*nph_b
    end associate
    allocate(cid(tree%n_leaf), pos(0:ncol-1), own(tree%n_leaf), lev(tree%n_leaf), per(0:n_ranks-1))
    pos = -1
    n = 0
    do j = 1, tree%n_leaf
      cid(j) = ((tree%key(KP, j) - 1)*nth_b + tree%key(KI+2, j))*nph_b + tree%key(KI+3, j)
      if (pos(cid(j)) < 0) then
        pos(cid(j)) = n
        n = n + 1
      end if
    end do
    per = 0
    do j = 0, ncol - 1
      if (pos(j) >= 0) per(int((int(pos(j), 8)*n_ranks)/n)) = per(int((int(pos(j), 8)*n_ranks)/n)) + 1
    end do
    do j = 1, tree%n_leaf
      own(j) = int((int(pos(cid(j)), 8)*n_ranks)/n)
    end do
    lev = 1
    call arbor_redistribute(tree, own, lev)
    call arbor_ghost_rebuild(tree)
    si_cols = .true.
    write(msg, '(A,I0,A,I0,A,I0,A,I0,A)') '  si_radial: ', n, ' columns to ', n_ranks, ' ranks (', minval(per), &
        ' .. ', maxval(per), ' per rank)'
    call arbor_print(msg)
  end subroutine cz_columns_to_ranks

  ! The implicit range (whole radial blocks covering si_z_km), the local
  ! columns and the per-block arrays (after every (re)distribution)
  subroutine si_setup()
    integer :: g0, g1, i, jb, ib, c, q, nk
    character(len=200) :: msg
    call cz_unbind()
    g0 = 1;  g1 = n_r
    do i = 1, n_r
      if (z_c(i) >= si_z_km(1)*1.0e5_dp) then
        g0 = i
        exit
      end if
    end do
    do i = n_r, 1, -1
      if (z_c(i) <= si_z_km(2)*1.0e5_dp) then
        g1 = i
        exit
      end if
    end do
    if (g1 < g0) call arbor_abort('cz3d: si_z_km holds no cell')
    si_kb0 = (g0 - 1)/nb(1)
    si_kb1 = (g1 - 1)/nb(1)
    si_i0 = si_kb0*nb(1) + 1
    si_i1 = (si_kb1 + 1)*nb(1)
    nk = si_kb1 - si_kb0 + 1
    if (allocated(si_cb)) deallocate(si_cb)
    allocate(si_cb(nk, nBlock))
    si_cb = 0
    si_ncol = 0
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (.not. si_in(ib)) cycle
      c = 0
      do q = 1, si_ncol
        ib_c: block
          integer :: ib1
          ib1 = maxval(si_cb(:, q))
          if (all(Key_IB([KP, KI+2, KI+3], ib1) == Key_IB([KP, KI+2, KI+3], ib))) c = q
        end block ib_c
        if (c > 0) exit
      end do
      if (c == 0) then
        si_ncol = si_ncol + 1
        c = si_ncol
      end if
      si_cb(Key_IB(KI+1, ib) - si_kb0 + 1, c) = ib
    end do
    if (any(si_cb(:, 1:si_ncol) == 0)) &
        call arbor_abort('cz3d: si_radial: a column is not whole on its rank (cz_columns_to_ranks not called?)')
    if (allocated(siS_B)) deallocate(siS_B, siR_B)
    allocate(siS_B(nb(1), nb(2), nb(3), nBlockMax), siR_B(nb(1), nb(2), nb(3), nv, nBlockMax))
    siS_B = 0.0_dp
    if (allocated(siSa_B)) deallocate(siSa_B)
    if (si_alfven) then
      allocate(siSa_B(nb(1), nb(2), nb(3), nBlockMax))
      siSa_B = 0.0_dp
    end if
    if (allocated(siV_B)) deallocate(siV_B, siX_B)
    if (si_diff_precond) allocate(siV_B(nb(1), nb(2), nb(3), nBlockMax), siX_B(nb(1), nb(2), nb(3), nBlockMax))
    if (allocated(EtaB_GB)) deallocate(EtaB_GB)
    if (use_B .and. use_diffusion .and. b_diff_form == 'resistive') then
      if (si_i0 /= 1 .or. si_i1 /= n_r) &
          call arbor_abort("cz3d: b_diff_form = 'resistive' needs the whole column implicit (si_z_km default)")
      allocate(EtaB_GB(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng, 1, size(State_GVB, 5)))
      EtaB_GB = 0.0_dp
    end if
    write(msg, '(A,I0,A,I0,A,F9.1,A,F10.1,A)') '  si_radial: implicit cells ', si_i0, ' .. ', si_i1, ' (z ', &
        z_c(si_i0)/1.0e5_dp, ' .. ', z_c(si_i1)/1.0e5_dp, ' km; the two cells at each end explicit)'
    call arbor_print(msg)
  end subroutine si_setup

  ! b_smooth_cells > 0 (see cz3d_params_mod): the CZ's field smoothed before
  ! the PFSS/DPF seed, each spherical component of B on each shell by the same
  ! explicit horizontal diffusion in angle (independent of r, so it commutes
  ! with the radial part of div B; the curvature terms of the components are
  ! left out, an O(dtheta^2) error): per pass
  !   B += D [(1/s) d_th(s d_th B) + d_ph^2 B/s^2],  D = 0.1 min(dth, s dph)^2,  s = sin(theta)
  ! n = b_smooth_cells^2/0.2 passes (diffusion length b_smooth_cells cells).
  ! The CZ cells only (with the atmosphere: up to the join); psi set to 0.
  module subroutine cz_smooth_b()
    integer :: n, it, jb, ib, i, j, k, c
    real(dp) :: acc(4), t0, w, dth, dph, st, sp, sm, dd
    real(dp), allocatable :: B0(:,:,:,:,:), X(:,:,:)
    character(len=200) :: msg
    if (b_smooth_cells <= 0.0_dp) return
    t0 = arbor_wtime()
    n = max(1, ceiling(b_smooth_cells**2/0.2_dp))
    allocate(B0(nb(1), nb(2), nb(3), br_:bp_, size(State_GVB, 5)), X(nb(1), nb(2), nb(3)))
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      B0(:, :, :, :, ib) = State_GVB(1:nb(1), 1:nb(2), 1:nb(3), br_:bp_, ib)
    end do
    do it = 1, n
      call arbor_fill_ghosts(State_GVB, iVarMin=br_, iVarMax=bp_, n_ghost=1, faces_only=.true.)
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        if (use_atmosphere .and. Key_IB(KI+1, ib)*nb(1) >= k_join) cycle
        associate (th => Metric_DB(2, ib)%X_I, dth_ => Metric_DB(2, ib)%Dx_I, dph_ => Metric_DB(3, ib)%Dx_I)
        do c = br_, bp_
          do k = 1, nb(3)
            do j = 1, nb(2)
              dth = dth_(j)
              dph = dph_(k)
              st = sin(th(j))
              sp = sin(th(j) + 0.5_dp*dth)
              sm = sin(th(j) - 0.5_dp*dth)
              dd = 0.1_dp*min(dth, st*dph)**2
              X(:, j, k) = State_GVB(1:nb(1), j, k, c, ib) + dd*( &
                  (sp*(State_GVB(1:nb(1), j+1, k, c, ib) - State_GVB(1:nb(1), j, k, c, ib)) &
                  - sm*(State_GVB(1:nb(1), j, k, c, ib) - State_GVB(1:nb(1), j-1, k, c, ib)))/(st*dth**2) &
                  + (State_GVB(1:nb(1), j, k+1, c, ib) - 2.0_dp*State_GVB(1:nb(1), j, k, c, ib) &
                  + State_GVB(1:nb(1), j, k-1, c, ib))/(st*dph)**2)
            end do
          end do
          do i = 1, nb(1)
            if (use_atmosphere .and. Key_IB(KI+1, ib)*nb(1) + i > k_join) cycle
            State_GVB(i, 1:nb(2), 1:nb(3), c, ib) = X(i, :, :)
          end do
        end do
        end associate
      end do
    end do
    ! report: rms |B - B0| / rms |B| over the CZ's top 3 Mm and over the whole CZ (volume weighted)
    acc = 0.0_dp
    call cz_unbind()
    w = r_f(merge(k_join, n_r, use_atmosphere))
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      State_GVB(:, :, :, psi_, ib) = 0.0_dp
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            if (use_atmosphere .and. r_off + i > k_join) cycle
            dd = r_c(i)**2*(r_f(i) - r_f(i-1))*Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
            st = sum((State_GVB(i, j, k, br_:bp_, ib) - B0(i, j, k, :, ib))**2)
            sp = sum(B0(i, j, k, :, ib)**2)
            acc(3:4) = acc(3:4) + dd*[st, sp]
            if (w - r_c(i) < 3.0e8_dp) acc(1:2) = acc(1:2) + dd*[st, sp]
          end do
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(acc)
    write(msg, '(A,F6.2,A,I0,A,ES10.3,A,ES10.3,A,F8.2,A)') '  b_smooth: ', b_smooth_cells, ' cells, ', n, &
        ' passes: rms|dB|/rms|B| ', sqrt(acc(1)/max(acc(2), tiny(1.0_dp))), ' (top 3 Mm of the CZ), ', &
        sqrt(acc(3)/max(acc(4), tiny(1.0_dp))), ' (CZ); ', arbor_wtime() - t0, ' s'
    call arbor_print(msg)
  end subroutine cz_smooth_b

  ! Background enthalpy h0(r) (EOS core) and potential Phi(r) (Phi(1) = 0)
  subroutine background_h_phi()
    type(eos_state) :: st
    real(dp) :: out(n_field)
    integer :: i
    allocate(h0_r_g(-1:n_r+2), phi_r_g(-1:n_r+2))
    h0_r => h0_r_g
    phi_r => phi_r_g
    do i = -1, n_r + 2
      call eos_eval(tab, rho0(i), s0(i), out)
      call eos_get_from_rho_T(rho0(i), exp(out(F_LNT)), st, tab%comp, tab%coulomb)
      h0_r(i) = st%e + st%p/st%rho
    end do
    phi_r(1) = 0.0_dp
    do i = 1, n_r + 1
      phi_r(i+1) = phi_r(i) + 0.5_dp*(g_c(i) + g_c(i+1))*(r_c(i+1) - r_c(i))
    end do
    do i = 0, -1, -1
      phi_r(i) = phi_r(i+1) - 0.5_dp*(g_c(i) + g_c(i+1))*(r_c(i+1) - r_c(i))
    end do
  end subroutine background_h_phi

  ! The profile counter and schedule go into restart files
  module subroutine cz_restart_save(self, Data_I)
    class(cz3d_app_t), intent(inout) :: self
    real(dp), allocatable, intent(out) :: Data_I(:)
    Data_I = [real(i_prof, dp), t_next_prof, t_prof0, tloc]
  end subroutine cz_restart_save

  module subroutine cz_restart_load(self, Data_I)
    class(cz3d_app_t), intent(inout) :: self
    real(dp), intent(in) :: Data_I(:)
    if (size(Data_I) /= 2 .and. size(Data_I) /= 3 .and. size(Data_I) /= 3 + n_r) &
        call arbor_abort('cz3d: restart file from another application')
    i_prof = nint(Data_I(1))
    t_next_prof = Data_I(2)
    t_prof0 = 0.0_dp
    if (size(Data_I) >= 3) t_prof0 = Data_I(3)           ! (files before 2026-10-11 hold two numbers)
    if (.not. allocated(tloc)) allocate(tloc(1:n_r))
    tloc = -1.0_dp                                        ! (older files: set to the run's time by the caller)
    if (size(Data_I) == 3 + n_r) tloc = Data_I(4:)
  end subroutine cz_restart_load

end submodule cz3d_setup_smod
