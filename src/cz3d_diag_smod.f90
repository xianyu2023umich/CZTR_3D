submodule (cz3d_app_mod) cz3d_diag_smod
  ! diagnostics: after-step outputs, log, profiles, energy budget, probes and boxes, shell dumps, plots, checksum, the angular-momentum correction
  ! (part of cz3d_app_mod: the module's data and routines by host association)
  implicit none

contains

  ! Longitude means in the global frame (rotation axis = global z), per radial
  ! cell and latitude bin (n_lat_bin bins of equal latitude width): the zonal
  ! velocity <v_phi_g>, Omega - Omega_0 = <v_phi_g>/(r cos lat) [nHz], and
  ! the meridional velocity <v_theta_g> (> 0 southward); solid-angle weights
  ! with the Yin-Yang overlap counted once. out_dir/rotation_NNNN.dat (the
  ! profile's number).
  subroutine rotation_profile(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    integer, parameter :: n_lat_bin = 36
    real(dp) :: A_IV(1:n_r, 3*n_lat_bin), e(3,3), x(3), vc(3), w, th, ph, th2, ph2, rs, lat, eph(3), eth(3), lc
    integer :: jb, ib, i, j, k, m, u
    character(len=300) :: name
    A_IV = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          th = Metric_DB(2, ib)%X_I(j);  ph = Metric_DB(3, ib)%X_I(k)
          if (.not. core(th, ph)) cycle
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          call yy_angles(th, ph, th2, ph2)
          if (core(th2, ph2)) w = 0.5_dp*w
          e = sph_basis(th, ph)
          x = e(:, 1)
          if (Key_IB(KP, ib) == 2) x = yy_cart(x)
          rs = max(hypot(x(1), x(2)), 1.0e-12_dp)
          lat = asin(max(-1.0_dp, min(1.0_dp, x(3))))
          m = min(n_lat_bin, 1 + int((lat + 0.5_dp*pi)/pi*n_lat_bin))
          eph = [-x(2), x(1), 0.0_dp]/rs
          eth = [x(3)*x(1)/rs, x(3)*x(2)/rs, -rs]
          do i = 1, nb(1)
            vc = matmul(e, State_GVB(i, j, k, vr_:vp_, ib))
            if (Key_IB(KP, ib) == 2) vc = yy_cart(vc)
            A_IV(r_off + i, m) = A_IV(r_off + i, m) + w
            A_IV(r_off + i, n_lat_bin + m) = A_IV(r_off + i, n_lat_bin + m) + w*dot_product(vc, eph)
            A_IV(r_off + i, 2*n_lat_bin + m) = A_IV(r_off + i, 2*n_lat_bin + m) + w*dot_product(vc, eth)
          end do
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(A_IV)
    if (.not. is_root) return
    write(name, '(A,A,I4.4,A)') trim(out_dir), '/rotation_', i_prof, '.dat'
    open(newunit=u, file=trim(name), status='replace', action='write')
    write(u, '(A,ES14.6,A,I0,A,ES12.4,A)') '# t = ', t, ' s; step ', iStep, '; omega_rot = ', omega_rot, &
        ' rad/s; longitude means in the global frame'
    write(u, '(A)') '# i  r[cm]  z[km]  lat[deg]  <v_phi>[cm/s]  Omega-Omega0[nHz]  <v_theta>[cm/s] (>0 southward)'
    do i = 1, n_r
      do m = 1, n_lat_bin
        if (A_IV(i, m) <= 0.0_dp) cycle
        lc = -0.5_dp*pi + (m - 0.5_dp)*pi/n_lat_bin
        write(u, '(I5,ES16.8,F12.1,F8.2,3ES14.5)') i, r_c(i), z_c(i)/1.0e5_dp, lc*180.0_dp/pi, &
            A_IV(i, n_lat_bin + m)/A_IV(i, m), A_IV(i, n_lat_bin + m)/A_IV(i, m)/(r_c(i)*cos(lc))/(2.0_dp*pi)*1.0e9_dp, &
            A_IV(i, 2*n_lat_bin + m)/A_IV(i, m)
      end do
    end do
    close(u)
  contains
    logical function core(tt, pp)
      real(dp), intent(in) :: tt, pp
      core = tt >= 0.25_dp*pi .and. tt <= 0.75_dp*pi .and. abs(pp) <= 0.75_dp*pi
    end function core
  end subroutine rotation_profile

  ! Remove the non-Coriolis part of L_z's change since the reference: delta_L
  ! = L_z - (lz_ref + lz_imp), delta_Omega = delta_L / I_z (I_z = int rho s^2
  ! dV), v -= delta_Omega z x x (global frame, back to each patch's components)
  subroutine lz_correct(t)
    real(dp), intent(in) :: t
    integer :: jb, ib, i, j, k
    real(dp) :: e(3,3), x(3), dv(3), w, th, ph, th2, ph2, lz, iz, dom
    character(len=200) :: msg
    if (.not. lz_ref_set) return
    call angular_momentum(lz)
    iz = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          th = Metric_DB(2, ib)%X_I(j);  ph = Metric_DB(3, ib)%X_I(k)
          if (.not. yy_core(th, ph)) cycle
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          call yy_angles(th, ph, th2, ph2)
          if (yy_core(th2, ph2)) w = 0.5_dp*w
          e = sph_basis(th, ph)
          x = e(:, 1)
          if (Key_IB(KP, ib) == 2) x = yy_cart(x)
          do i = 1, nb(1)
            iz = iz + w*Metric_DB(1, ib)%Vol_I(i)*(rho0(i) + State_GVB(i, j, k, rho1_, ib))*r_c(i)**2*(x(1)**2 + x(2)**2)
          end do
        end do
      end do
    end do
    call cz_unbind()
    iz = arbor_sum(iz)
    lz_ref = lz_ref + lz_imp
    lz_imp = 0.0_dp
    dom = (lz - lz_ref)/iz
    ! v -= dom (z x x) everywhere (all cells, the overlap included)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
          do i = 1, nb(1)
            x = r_c(i)*e(:, 1)
            if (Key_IB(KP, ib) == 2) x = yy_cart(x)
            dv = -dom*[-x(2), x(1), 0.0_dp]
            if (Key_IB(KP, ib) == 2) dv = yy_cart(dv)
            State_GVB(i, j, k, vr_:vp_, ib) = State_GVB(i, j, k, vr_:vp_, ib) + matmul(transpose(e), dv)
          end do
        end do
      end do
    end do
    call cz_unbind()
    if (lz_correct_dt > 0.0_dp) then
      do while (t_next_lzc <= t)
        t_next_lzc = t_next_lzc + lz_correct_dt
      end do
    end if
    if (abs(dom) > 1.0e-3_dp*max(abs(omega_rot), 1.0e-6_dp)) then
      write(msg, '(A,ES11.3,A)') 'cz3d: large L_z correction, delta_Omega = ', dom, ' rad/s'
      call arbor_print(msg)
    end if
  end subroutine lz_correct

  ! Energy exchange between the flow and the field over the sphere (ghosts
  ! filled): the Lorentz force's work wl = int v.(J x B)/(4 pi) dV and the
  ! induction term's change of the magnetic energy wi = int B.curl(v x B)/(4 pi)
  ! dV [erg/s]. With walls that pass no Poynting flux (conductor: B_r = v_r =
  ! 0; radial: B_h = v_r = 0), wl + wi = 0 up to the discretisation: a check of
  ! the signs and factors of both terms. Yin-Yang overlap weights as the budget's.
  subroutine mhd_exchange(wl, wi)
    real(dp), intent(out) :: wl, wi
    integer :: jb, ib, i, j, k
    real(dp) :: w, th, ph, th2, ph2, v(3), b(3), jj(3), vol
    wl = 0.0_dp;  wi = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      B_GV = State_GVB(:, :, :, br_:bp_, ib)
      call curl_block(B_GV, ib, J_CV)
      do k = 1 - ng, nb(3) + ng
        do j = 1 - ng, nb(2) + ng
          do i = 1 - ng, nb(1) + ng
            v = State_GVB(i, j, k, vr_:vp_, ib)
            b = State_GVB(i, j, k, br_:bp_, ib)
            B_GV(i, j, k, :) = [v(2)*b(3) - v(3)*b(2), v(3)*b(1) - v(1)*b(3), v(1)*b(2) - v(2)*b(1)]
          end do
        end do
      end do
      call curl_block(B_GV, ib, C_CV)
      do k = 1, nb(3)
        do j = 1, nb(2)
          th = Metric_DB(2, ib)%X_I(j);  ph = Metric_DB(3, ib)%X_I(k)
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          if (.not. core(th, ph)) cycle
          call yy_angles(th, ph, th2, ph2)
          if (core(th2, ph2)) w = 0.5_dp*w
          do i = 1, nb(1)
            vol = w*Metric_DB(1, ib)%Vol_I(i)
            v = State_GVB(i, j, k, vr_:vp_, ib)
            b = State_GVB(i, j, k, br_:bp_, ib)
            jj = J_CV(i, j, k, :)
            wl = wl + vol*lorentz_f(sum(b**2), rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib), &
                cell_cap(i, j, k, ib, cs0(i)/xi(i))) &
                *dot_product(v, [jj(2)*b(3) - jj(3)*b(2), jj(3)*b(1) - jj(1)*b(3), jj(1)*b(2) - jj(2)*b(1)])/four_pi
            wi = wi + vol*dot_product(b, C_CV(i, j, k, :))/four_pi
          end do
        end do
      end do
    end do
    call cz_unbind()
    wl = arbor_sum(wl);  wi = arbor_sum(wi)
  contains
    logical function core(t, p)
      real(dp), intent(in) :: t, p
      core = t >= 0.25_dp*pi .and. t <= 0.75_dp*pi .and. abs(p) <= 0.75_dp*pi
    end function core
  end subroutine mhd_exchange

  ! h0(r) and Phi(r) of the background, the overlap weights, the shells'
  ! energy at the start, the output files. Call after cz_setup at the run's
  ! start (also after a restart: the budget counts from there).
  module subroutine budget_start(t)
    real(dp), intent(in) :: t
    real(dp) :: th2, ph2
    integer :: j, k, jb, ib, p
    if (dt_budget <= 0.0_dp) return
    if (rsst_h2015) call arbor_abort('cz3d: the energy budget is built for rsst_form = hotta2012')
    budget_on = .true.
    allocate(wyy_B(nb(2), nb(3), nBlockMax), acc_IV(1:n_r, nAcc), E0_I(1:n_r), bgf_r(1:n_r))
    allocate(work_CB(nb(1), nb(2), nb(3), nBlockMax), workh_CB(nb(1), nb(2), nb(3), nBlockMax))
    work_CB = 0.0_dp
    workh_CB = 0.0_dp
    bgf_r = ddr_c(p0) + rho0(1:n_r)*g_c(1:n_r)
    ! weight: solid angle; 1/2 where the cell centre is in both patches'
    ! cores (theta in [pi/4, 3pi/4], phi in [-3pi/4, 3pi/4]), 0 in a margin
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      p = Key_IB(KP, ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          wyy_B(j, k, ib) = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          if (.not. in_core(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))) then
            wyy_B(j, k, ib) = 0.0_dp
          else
            call yy_angles(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k), th2, ph2)
            if (in_core(th2, ph2)) wyy_B(j, k, ib) = 0.5_dp*wyy_B(j, k, ib)
          end if
        end do
      end do
    end do
    acc_IV = 0.0_dp
    call shell_energy(E0_I)
    t_bud0 = t
    t_next_bud = t + dt_budget
    if (is_root .and. restarted) then
      ! a continued run appends (its records count from the restart time t0)
      open(newunit=u_bud, file=trim(out_dir)//'/budget_shells.dat', status='unknown', position='append', action='write')
      write(u_bud, '(A,ES16.8)') '# continued from a restart file; changes below count from t0 = ', t
      open(newunit=u_bsum, file=trim(out_dir)//'/budget_log.txt', status='unknown', position='append', action='write')
      write(u_bsum, '(A,ES16.8,A)') '# continued from a restart file: mean rates since t0 = ', t, ' s'
      open(newunit=u_bkh, file=trim(out_dir)//'/budget_kh.dat', status='unknown', position='append', action='write')
      open(newunit=u_bsp, file=trim(out_dir)//'/budget_split.dat', status='unknown', position='append', action='write')
      write(u_bsp, '(A,ES16.8)') '# continued from a restart file; rates below count from t0 = ', t
    else if (is_root) then
      open(newunit=u_bud, file=trim(out_dir)//'/budget_shells.dat', status='replace', action='write')
      write(u_bud, '(A)') '# energy budget per radial shell [erg], changes since the start (t0 below): records of'
      write(u_bud, '(A)') '# i  z[km]  dE_total  dE_hydro  dE_radiation  dE_cooling  dE_diffusion  E_rsst_source  ' // &
          'dE_total-(sum of operators)  W_background_force  E_rsst_source_fluct   (the RSST source and the ' // &
          'background-force work are part of dE_hydro)'
      open(newunit=u_bsum, file=trim(out_dir)//'/budget_log.txt', status='replace', action='write')
      write(u_bsum, '(A,ES16.8,A)') '# energy budget of the whole domain since t0 = ', t, ' s, as mean rates in Lsun:'
      write(u_bsum, '(A)') '# step  t  dE_total/dt  hydro  radiation  cooling  diffusion  rsst_source  residual  ' // &
          'bg_force_work  rsst_fluct  hydro-rsst-bg_work (scheme error)'
      open(newunit=u_bkh, file=trim(out_dir)//'/budget_kh.dat', status='replace', action='write')
      write(u_bkh, '(A)') '# horizontal kinetic energy K_h = rho v_h^2/2 per radial shell, mean rates in Lsun since t0 ' // &
          '(records): i  z[km]  dK_h(hydro)  W_h = -v_h.grad_h p1  F_Kh(c) = rho v_r v_h^2/2 over the sphere  ' // &
          'dK_h(diffusion)  dK(diffusion)  dK(hydro)  W  F_K(c)'
      open(newunit=u_bsp, file=trim(out_dir)//'/budget_split.dat', status='replace', action='write')
      write(u_bsp, '(A)') '# hydro operator split into kinetic and thermal energy per radial shell, mean rates in ' // &
          'Lsun since t0 (records): i  z[km]  dE_kin  dE_th  W  rsst_kin  rsst_total  bg_force_work  L_res_total(c)  ' // &
          'L_res_kin(c)   (dE_th = rho T ds + (h + Phi) drho; W = -v.grad p1 - g rho1 v_r; L_res at the cell centre: ' // &
          'rho v_r (h + v^2/2 + Phi) and rho v_r v^2/2 over the sphere; see scripts/analyze_split.py)'
    end if
  contains
    logical function in_core(th, ph)
      real(dp), intent(in) :: th, ph
      in_core = th >= 0.25_dp*pi .and. th <= 0.75_dp*pi .and. abs(ph) <= 0.75_dp*pi
    end function in_core
  end subroutine budget_start

  ! E per shell (e from the EOS core at the cell's rho, T)
  subroutine shell_energy(E_I)
    real(dp), intent(out) :: E_I(1:n_r)
    type(eos_state) :: st
    real(dp) :: out(n_field), rho, w
    integer :: jb, ib, i, j, k
    E_I = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          w = wyy_B(j, k, ib)
          if (w == 0.0_dp) cycle
          do i = 1, nb(1)
            rho = rho0(i) + State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            call eos_get_from_rho_T(rho, exp(out(F_LNT)), st, tab%comp, tab%coulomb)
            E_I(r_off + i) = E_I(r_off + i) + w*Metric_DB(1, ib)%Vol_I(i)*(rho*(st%e + phi_r(i)) &
                + 0.5_dp*rho*sum(State_GVB(i, j, k, vr_:vp_, ib)**2))
          end do
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(E_I)
  end subroutine shell_energy

  ! Energy change of the RK step of block ib (StateOld -> State) and the
  ! RSST source from the step's drho
  module subroutine budget_hydro(ib, dt)
    integer, intent(in) :: ib
    real(dp), intent(in) :: dt
    real(dp) :: out(n_field), rho_o, rho_n, T_o, h_o, w, vol, ke_o, ke_n
    integer :: i, j, k
    do k = 1, nb(3)
      do j = 1, nb(2)
        w = wyy_B(j, k, ib)
        if (w == 0.0_dp) cycle
        do i = 1, nb(1)
          rho_o = rho0(i) + StateOld_GVB(i, j, k, rho1_, ib)
          rho_n = rho0(i) + State_GVB(i, j, k, rho1_, ib)
          call eos_eval(tab, rho_o, s0(i) + StateOld_GVB(i, j, k, s1_, ib), out)
          T_o = exp(out(F_LNT))
          h_o = h0_r(i) + T_o*StateOld_GVB(i, j, k, s1_, ib) + p0(i)*expm1(out(F_LNP) - lnp0(i))/rho_o
          ke_o = 0.5_dp*sum(StateOld_GVB(i, j, k, vr_:vp_, ib)**2)
          ke_n = 0.5_dp*rho_n*sum(State_GVB(i, j, k, vr_:vp_, ib)**2)
          vol = w*Metric_DB(1, ib)%Vol_I(i)
          acc_IV(i, 1) = acc_IV(i, 1) + vol*(rho_o*T_o*(State_GVB(i, j, k, s1_, ib) - StateOld_GVB(i, j, k, s1_, ib)) &
              + (h_o + phi_r(i))*(rho_n - rho_o) + ke_n - rho_o*ke_o)
          ! RSST source over the step: drho = -dt div(rho v)/xi^2 (the RK
          ! step's actual change), so S dt = -(xi^2 - 1)(h + v^2/2 + Phi) drho
          acc_IV(i, 5) = acc_IV(i, 5) - vol*(xi(i)**2 - 1.0_dp)*(h_o + ke_o + phi_r(i))*(rho_n - rho_o)
          acc_IV(i, 7) = acc_IV(i, 7) - vol*(xi(i)**2 - 1.0_dp)*(h_o - h0_r(i) + ke_o)*(rho_n - rho_o)
          acc_IV(i, 6) = acc_IV(i, 6) + vol*dt*StateOld_GVB(i, j, k, vr_, ib)*bgf_r(i)
          acc_IV(i, 8) = acc_IV(i, 8) + vol*(ke_n - rho_o*ke_o)
          acc_IV(i, 9) = acc_IV(i, 9) + vol*(rho_o*T_o*(State_GVB(i, j, k, s1_, ib) - StateOld_GVB(i, j, k, s1_, ib)) &
              + (h_o + phi_r(i))*(rho_n - rho_o))
          acc_IV(i, 10) = acc_IV(i, 10) + vol*dt*work_CB(i, j, k, ib)
          acc_IV(i, 11) = acc_IV(i, 11) - vol*(xi(i)**2 - 1.0_dp)*ke_o*(rho_n - rho_o)
          acc_IV(i, 12) = acc_IV(i, 12) + w*dt*r_c(i)**2*rho_o*StateOld_GVB(i, j, k, vr_, ib)*(h_o + ke_o + phi_r(i))
          acc_IV(i, 13) = acc_IV(i, 13) + w*dt*r_c(i)**2*rho_o*StateOld_GVB(i, j, k, vr_, ib)*ke_o
          acc_IV(i, 14) = acc_IV(i, 14) + vol*0.5_dp*(rho_n*sum(State_GVB(i, j, k, vt_:vp_, ib)**2) &
              - rho_o*sum(StateOld_GVB(i, j, k, vt_:vp_, ib)**2))
          acc_IV(i, 15) = acc_IV(i, 15) + vol*dt*workh_CB(i, j, k, ib)
          acc_IV(i, 16) = acc_IV(i, 16) + w*dt*r_c(i)**2*rho_o*StateOld_GVB(i, j, k, vr_, ib) &
              *0.5_dp*sum(StateOld_GVB(i, j, k, vt_:vp_, ib)**2)
        end do
      end do
    end do
  end subroutine budget_hydro

  ! One budget record: per shell, and the whole domain as mean rates
  subroutine budget_write(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    real(dp) :: E_I(1:n_r), A_IV(1:n_r, nAcc), dE(1:n_r), f
    integer :: i
    call shell_energy(E_I)
    A_IV = acc_IV
    call arbor_sum_inplace(A_IV)
    t_next_bud = t_next_bud + dt_budget
    if (.not. is_root) return
    dE = E_I - E0_I
    write(u_bud, '(A,I0,A,ES16.8,A,ES16.8)') '# step= ', iStep, ' t= ', t, ' t0= ', t_bud0
    do i = 1, n_r
      write(u_bud, '(I5,F12.1,9ES16.7E3)') i, z_c(i)/1.0e5_dp, dE(i), A_IV(i, 1:5), dE(i) - sum(A_IV(i, 1:4)), &
          A_IV(i, 6:7)
    end do
    write(u_bud, '(A)') ''
    flush(u_bud)
    f = 1.0_dp/(max(t - t_bud0, 1.0e-30_dp)*Lsun)
    write(u_bsum, '(I9,ES16.8,10F12.5)') iStep, t, sum(dE)*f, (sum(A_IV(:, i))*f, i = 1, 5), &
        (sum(dE) - sum(A_IV(:, 1:4)))*f, sum(A_IV(:, 6))*f, sum(A_IV(:, 7))*f, &
        (sum(A_IV(:, 1)) - sum(A_IV(:, 5)) - sum(A_IV(:, 6)))*f
    flush(u_bsum)
    write(u_bsp, '(A,I0,A,ES16.8,A,ES16.8)') '# step= ', iStep, ' t= ', t, ' t0= ', t_bud0
    do i = 1, n_r
      write(u_bsp, '(I5,F12.1,8ES16.7E3)') i, z_c(i)/1.0e5_dp, A_IV(i, 8)*f, A_IV(i, 9)*f, A_IV(i, 10)*f, &
          A_IV(i, 11)*f, A_IV(i, 5)*f, A_IV(i, 6)*f, A_IV(i, 12)*f, A_IV(i, 13)*f
    end do
    write(u_bsp, '(A)') ''
    flush(u_bsp)
    write(u_bkh, '(A,I0,A,ES16.8,A,ES16.8)') '# step= ', iStep, ' t= ', t, ' t0= ', t_bud0
    do i = 1, n_r
      write(u_bkh, '(I5,F12.1,8ES16.7E3)') i, z_c(i)/1.0e5_dp, A_IV(i, 14)*f, A_IV(i, 15)*f, A_IV(i, 16)*f, &
          A_IV(i, 17)*f, A_IV(i, 18)*f, A_IV(i, 8)*f, A_IV(i, 10)*f, A_IV(i, 13)*f
    end do
    write(u_bkh, '(A)') ''
    flush(u_bkh)
  end subroutine budget_write

  module subroutine cz_after_step(self, iStep, t)
    class(cz3d_app_t), intent(inout) :: self
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    call tic(T_AFTER)
    step_now = iStep
    if (probe_radius_deg > 0.0_dp) call probe_write(iStep, t)
    if (box_radius_deg > 0.0_dp .and. box_every > 0) then
      if (mod(iStep, box_every) == 0 .or. (box_fine_from > 0 .and. iStep >= box_fine_from)) call box_write(iStep, t)
    end if
    if (euv_dt > 0.0_dp) then
      if (t_next_euv == -huge(1.0_dp)) t_next_euv = t
      if (t >= t_next_euv) then
        call euv_write(iStep, t)
        do while (t_next_euv <= t)
          t_next_euv = t_next_euv + euv_dt
        end do
      end if
    end if
    call nan_check(iStep, t)
    if (lz_correct_steps > 0 .or. lz_correct_dt > 0.0_dp) then
      if (lz_correct_steps > 0) then
        if (mod(iStep, lz_correct_steps) == 0) call lz_correct(t)
      end if
      if (lz_correct_dt > 0.0_dp .and. t >= t_next_lzc) call lz_correct(t)
    end if
    if (mod(iStep, n_log) == 0 .or. t >= t_end*(1.0_dp - 1.0e-12_dp)) then
      call cz_log(iStep, t)
      call dt_report(iStep)
      if (use_atmosphere) call footpoint_report(iStep, t)
    end if
    if (prof_every > 0) then
      if (mod(iStep, prof_every) == 0) call cz_profile(iStep, t)
    else if (dt_prof > 0.0_dp .and. t >= t_next_prof - 1.0e-9_dp*dt_prof) then
      call cz_profile(iStep, t)
    end if
    if (budget_on .and. (t >= t_next_bud - 1.0e-9_dp*dt_budget .or. t >= t_end*(1.0_dp - 1.0e-12_dp))) &
        call budget_write(iStep, t)
    dtl = huge(1.0_dp)                    ! (the dt report: the next step's minima)
    call toc(T_AFTER)
  end subroutine cz_after_step

  ! The dt report (every log line, stdout): per constraint (CFL sum, wave cap,
  ! losses) the smallest dt over all ranks in this step and its cell
  subroutine dt_report(iStep)
    integer, intent(in) :: iStep
    real(dp) :: v(8), g
    integer :: c, r, rmin
    character(len=240) :: msg
    character(len=*), parameter :: nm(3) = ['CFL sum  ', 'wave cap ', 'losses   ']
    character(len=*), parameter :: dn(3) = ['radial', 'theta ', 'phi   ']
    write(msg, '(A,I0,A)') 'dt limits (step ', iStep, '):'
    call arbor_print(msg)
    do c = 1, 3
      g = arbor_min(dtl(c, 1))
      if (g >= huge(1.0_dp)) cycle
      r = merge(my_rank, n_ranks, dtl(c, 1) == g)
      rmin = arbor_min(r)
      v = 0.0_dp
      if (my_rank == rmin) v = dtl(c, :)
      v = arbor_sum(v)
      if (c == 1) then
        write(msg, '(4X,A,F8.4,A,I4,A,F10.1,A,F7.1,A,F7.1,A,F7.1,A,F7.1,A,A)') nm(c), g, ' s  cell ', nint(v(2)), &
            ' (z ', v(3), ' km, dr ', v(8), ' km): |v_r| ', v(4)/1.0e5_dp, ', radial wave ', v(5)/1.0e5_dp, &
            ', horizontal fast ', v(6)/1.0e5_dp, ' km/s; dominant ', trim(dn(nint(v(7))))
      else if (c == 2) then
        write(msg, '(4X,A,F8.4,A,I4,A,F10.1,A,F7.1,A,F7.1,A)') nm(c), g, ' s  cell ', nint(v(2)), ' (z ', v(3), &
            ' km, dr ', v(8), ' km): c ', v(5)/1.0e5_dp, ' km/s'
      else
        write(msg, '(4X,A,F8.4,A,I4,A,F10.1,A)') nm(c), g, ' s  cell ', nint(v(2)), ' (z ', v(3), ' km)'
      end if
      call arbor_print(msg)
    end do
  end subroutine dt_report

  ! The stage timers at the end of the run (stdout): per stage the mean, min
  ! and max over the ranks [s], the mean's share of the run's wall time, and
  ! max/mean (load imbalance; time waiting for other ranks shows up in the
  ! ghost fill and in 'other': ARBOR's loop, the dt reduction, restart and
  ! plot files)
  module subroutine timers_report(wall, n_step)
    real(dp), intent(in) :: wall
    integer, intent(in) :: n_step
    real(dp) :: x(n_tim + 2), mn, mean, mx
    integer :: k, kk
    ! the print order: the diffusion's parts after it
    integer, parameter :: order(n_tim + 2) = [T_ADV, T_GHOST, T_G_RK, T_G_ATM, T_G_DIF, T_G_JOIN, T_R_OVL, T_RHS, T_R_EOS, T_R_MASS, T_R_ADV, T_R_PRS, T_R_MCURL, T_R_MLOOP, T_R_SIC, T_R_UPD, T_SI, T_S_SOL, T_TRANS, T_G_TR, T_ATM, T_A_OVL, T_DIFF, T_D_SIP, T_D_ETA, T_G_ETA, T_D_PRE, &
        T_D_EOS, T_D_FA, T_D_FB, T_D_FR, T_D_FE, T_D_DIV, T_D_UPD, T_D_IMP, n_tim + 1, T_DT, T_AFTER, n_tim + 2]
    character(len=200) :: msg
    character(len=24) :: lab
    x(1:n_tim) = tim_acc
    x(n_tim + 1) = tim_acc(T_ADV) - sum(tim_acc(T_GHOST:T_DIFF))       ! the rest of advance
    x(n_tim + 2) = wall - tim_acc(T_ADV) - tim_acc(T_DT) - tim_acc(T_AFTER)
    write(msg, '(A,I0,A,F9.1,A,I0,A)') 'timers: ', n_step, ' steps, wall ', wall, ' s, ', n_ranks, &
        ' ranks; per stage [s]: mean, min, max over ranks, mean/wall, max/mean'
    call arbor_print(msg)
    do kk = 1, n_tim + 2
      k = order(kk)
      mean = arbor_sum(x(k))/n_ranks
      mn = arbor_min(x(k))
      mx = arbor_max(x(k))
      if (k <= n_tim) then
        lab = tim_name(k)
      else if (k == n_tim + 1) then
        lab = '  other (in advance)'
      else
        lab = 'other (ARBOR, I/O, MPI)'
      end if
      write(msg, '(2X,A24)') lab
      write(msg(27:), '(3F10.2,F8.1,A,F7.2)') mean, mn, mx, 100.0_dp*mean/max(wall, tiny(1.0_dp)), ' %', &
          mx/max(mean, tiny(1.0_dp))
      call arbor_print(msg)
    end do
    if (b_eta_explicit) then
      x(1:2) = eta_cap_n
      call arbor_sum_inplace(x(1:2))
      write(msg, '(A,ES10.3,A,ES10.3,A)') "b_eta_radial = 'explicit': eta capped in ", x(1), ' of ', x(2), ' cell-steps'
      call arbor_print(msg)
    end if
    call floors_report()
    call arbor_timer_report(6)                 ! (ARBOR's own: the ghost fill's halves, ...)
  end subroutine timers_report

  ! The floors' action over the run (atm_floors' counts): the fraction of the
  ! checked cell-calls that a floor changed, per height band of shells with any
  ! (one line per band of consecutive shells; atm_floors runs once per step)
  subroutine floors_report()
    use cztr3d_atm_mod, only: floor_cnt
    real(dp), allocatable :: c(:,:)
    integer :: i, i0, f
    character(len=200) :: msg
    if (.not. use_atmosphere) return
    allocate(c(4, n_r))
    c = 0.0_dp
    if (allocated(floor_cnt)) c = floor_cnt
    call arbor_sum_inplace(c)
    write(msg, '(A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A)') 'floors: rho_floor ', rho_floor, ' g/cm3 acted in ', &
        sum(c(2, :)), ' of ', sum(c(1, :)), ' cell-calls; T_floor ', T_floor, ' K:'
    call arbor_print(msg)
    ! the T floor's mean power over the run, as a flux through the photosphere (4 pi R^2)
    if (dt_run > 0.0_dp) then
      write(msg, '(A,ES10.3,A,ES10.3,A,F8.1,A)') '    T floor: energy added ', sum(c(4, :)), ' erg, mean power ', &
          sum(c(4, :))/dt_run, ' erg/s = ', sum(c(4, :))/dt_run/(4.0_dp*pi*r_f(k_join)**2)/1.0e3_dp, &
          ' x 1e3 erg/cm2/s over the sphere (Yin-Yang overlap counted twice, ~+6%)'
      call arbor_print(msg)
    end if
    if (cond_hot_to_cold) then
      block
        use cztr3d_atm_mod, only: hc_stat
        real(dp) :: hs(4)
        hs = hc_stat
        call arbor_sum_inplace(hs)
        write(msg, '(A,ES10.3,A,F8.5,A,F8.5,A)') '    cond_hot_to_cold: ', hs(1), ' face fluxes, set to 0: ', &
            hs(2)/max(hs(1), 1.0_dp), ' of them, ', hs(4)/max(hs(3), tiny(1.0_dp)), ' of their sum of |q|'
        call arbor_print(msg)
      end block
    end if
    ! for comparison: the radiative losses and the heating (the last step's
    ! Qloss_B, Qheat_B) integrated over the floor's heights, z 0 .. 1200 km
    block
      use cztr3d_atm_mod, only: Qloss_B, Qheat_B
      real(dp) :: ql(2), qs(n_r)
      integer :: jb, ib, j, k
      ql = 0.0_dp
      qs = 0.0_dp
      if (allocated(Qloss_B)) then
        do jb = 1, nBlock
          ib = iBlockList_I(jb)
          call cz_bind(ib)
          do k = 1, nb(3)
            do j = 1, nb(2)
              do i = 1, nb(1)
                if (r_off + i <= k_join .or. z_c(i) > 1.2e8_dp) cycle
                ql = ql + [Qloss_B(i, j, k, ib), Qheat_B(i, j, k, ib)]*Metric_DB(1, ib)%Vol_I(i) &
                    *Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
                qs(r_off + i) = qs(r_off + i) + Qloss_B(i, j, k, ib)*Metric_DB(1, ib)%Vol_I(i) &
                    *Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
              end do
            end do
          end do
        end do
        call cz_unbind()
      end if
      call arbor_sum_inplace(ql)
      write(msg, '(A,F8.1,A,F8.1,A)') '    z 0-1200 km, last step: losses ', ql(1)/(4.0_dp*pi*r_f(k_join)**2)/1.0e3_dp, &
          ', heating ', ql(2)/(4.0_dp*pi*r_f(k_join)**2)/1.0e3_dp, ' x 1e3 erg/cm2/s over the sphere'
      call arbor_print(msg)
      call arbor_sum_inplace(qs)
      call arbor_print('    per shell (z km: losses, T-floor power, x 1e3 erg/cm2/s):')
      msg = '     '
      do i = k_join + 1, n_r
        if (z_c(i) > 1.2e8_dp) exit
        if (len_trim(msg) > 150) then
          call arbor_print(msg)
          msg = '     '
        end if
        write(msg(len_trim(msg)+1:), '(1X,I0,A,F0.0,A,F0.0)') nint(z_c(i)/1.0e5_dp), ': ', &
            qs(i)/(4.0_dp*pi*r_f(k_join)**2)/1.0e3_dp, '/', c(4, i)/dt_run/(4.0_dp*pi*r_f(k_join)**2)/1.0e3_dp
      end do
      call arbor_print(msg)
    end block
    do f = 2, 3
      i = 1
      do while (i <= n_r)
        if (c(f, i) <= 0.0_dp) then
          i = i + 1
          cycle
        end if
        i0 = i
        do while (i < n_r)
          if (c(f, i + 1) <= 0.0_dp) exit
          i = i + 1
        end do
        write(msg, '(4X,A,A,I4,A,I4,A,F9.1,A,F9.1,A,F8.4,A,F8.4)') merge('rho', 'T  ', f == 2), ' floored in cells ', i0, &
            ' ..', i, ' (z ', z_c(i0)/1.0e5_dp, ' ..', z_c(i)/1.0e5_dp, ' km): fraction ', &
            minval(c(f, i0:i)/max(c(1, i0:i), 1.0_dp)), ' ..', maxval(c(f, i0:i)/max(c(1, i0:i), 1.0_dp))
        call arbor_print(msg)
        i = i + 1
      end do
    end do
    ! the T floor per shell where it acts in more than 1% of the cell-calls: z [km] and the fraction
    msg = '    T floor > 1% (z km: fraction):'
    do i = 1, n_r
      if (c(3, i) <= 0.01_dp*max(c(1, i), 1.0_dp)) cycle
      if (len_trim(msg) > 170) then
        call arbor_print(msg)
        msg = '     '
      end if
      write(msg(len_trim(msg)+1:), '(1X,I0,A,F5.2)') nint(z_c(i)/1.0e5_dp), ':', c(3, i)/c(1, i)
    end do
    call arbor_print(msg)
  end subroutine floors_report

  ! Footpoint monitor (use_atmosphere, every log line): out_dir/footpoints.txt
  ! with, over all columns, the lowest gas pressure in the shell nearest
  ! z = 1 Mm relative to that shell's mean (a drained chromosphere), the
  ! strongest downflow at z = 0.5-3 Mm (a siphon or drainage into a
  ! footpoint) and the largest Mach number at z = 0-3 Mm, each with its place
  ! (global latitude, longitude; z and T for the last two). The means count
  ! the cells of both Yin-Yang patches (the overlap twice).
  subroutine footpoint_report(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    integer, save :: u = 0, i1 = 0
    real(dp) :: out(n_field), rho, cs, p, zk, m, x(3), e(3,3), lat, lon, h
    real(dp) :: psum, pn, pmin, vmin, mmax, gp, gv, gm
    real(dp) :: lp(2), lv(4), lm(4)
    integer :: jb, ib, i, j, k, ig
    if (i1 == 0) i1 = minloc(abs(z_c(1:n_r) - 1.0e8_dp), 1)    ! (unbound views: the global shells)
    psum = 0.0_dp;  pn = 0.0_dp;  pmin = huge(1.0_dp);  vmin = huge(1.0_dp);  mmax = -1.0_dp
    lp = 0.0_dp;  lv = 0.0_dp;  lm = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do i = 1, nb(1)
        zk = z_c(i)/1.0e5_dp
        ig = Key_IB(KI+1, ib)*nb(1) + i
        if (ig /= i1 .and. (zk < 0.0_dp .or. zk > 3000.0_dp)) cycle
        do k = 1, nb(3)
          do j = 1, nb(2)
            rho = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            p = exp(out(F_LNP))
            cs = sqrt(out(F_GAMMA1)*p/rho)
            m = norm2(State_GVB(i, j, k, vr_:vp_, ib))/cs
            e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
            x = e(:, 1)
            if (Key_IB(KP, ib) == 2) x = yy_cart(x)
            lat = asin(max(-1.0_dp, min(1.0_dp, x(3))))*180.0_dp/pi
            lon = atan2(x(2), x(1))*180.0_dp/pi
            if (ig == i1) then
              psum = psum + p;  pn = pn + 1.0_dp
              if (p < pmin) then
                pmin = p;  lp = [lat, lon]
              end if
            end if
            if (zk >= 500.0_dp .and. zk <= 3000.0_dp .and. State_GVB(i, j, k, vr_, ib) < vmin) then
              vmin = State_GVB(i, j, k, vr_, ib);  lv = [lat, lon, zk, exp(out(F_LNT))]
            end if
            if (zk >= 0.0_dp .and. m > mmax) then
              mmax = m;  lm = [lat, lon, zk, exp(out(F_LNT))]
            end if
          end do
        end do
      end do
    end do
    call cz_unbind()
    psum = arbor_sum(psum);  pn = arbor_sum(pn)
    gp = arbor_min(pmin);  gv = arbor_min(vmin);  gm = arbor_max(mmax)
    ! the places, from the rank(s) that hold each extreme
    if (pmin /= gp) lp = 0.0_dp
    if (vmin /= gv) lv = 0.0_dp
    if (mmax /= gm) lm = 0.0_dp
    h = arbor_sum(merge(1.0_dp, 0.0_dp, pmin == gp))
    do i = 1, 2
      lp(i) = arbor_sum(lp(i))/max(h, 1.0_dp)
    end do
    h = arbor_sum(merge(1.0_dp, 0.0_dp, vmin == gv))
    do i = 1, 4
      lv(i) = arbor_sum(lv(i))/max(h, 1.0_dp)
    end do
    h = arbor_sum(merge(1.0_dp, 0.0_dp, mmax == gm))
    do i = 1, 4
      lm(i) = arbor_sum(lm(i))/max(h, 1.0_dp)
    end do
    if (.not. is_root) return
    if (u == 0) then
      open(newunit=u, file=trim(out_dir)//'/footpoints.txt', status='unknown', position='append', action='write')
      write(u, '(A,F8.1,A)') '# step  t[s]  p_min/<p> at z = ', z_c(i1)/1.0e5_dp, ' km, lat, lon | min v_r [km/s] at' // &
          ' z = 0.5-3 Mm, lat, lon, z[km], T | max Mach at z = 0-3 Mm, lat, lon, z[km], T'
    end if
    write(u, '(I9,ES14.6,ES11.3,2F8.2,ES11.3,2F8.2,F9.1,ES10.2,ES11.3,2F8.2,F9.1,ES10.2)') iStep, t, &
        gp/max(psum/max(pn, 1.0_dp), tiny(1.0_dp)), lp, gv/1.0e5_dp, lv, gm, lm
    flush(u)
  end subroutine footpoint_report

  subroutine box_write(iStep, t)
    use cztr3d_atm_mod, only: Qloss_B, Qheat_B
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    integer, save :: u = 0
    real(dp) :: x0(3), x(3), e(3,3), la, lo, out(n_field), rho, cs, Tc, p, b(3), v(3), b2, qh, ql, zz, ops(9)
    logical :: fine
    integer :: jb, ib, i, j, k
    character(len=300) :: name
    la = box_lat_deg*pi/180.0_dp;  lo = box_lon_deg*pi/180.0_dp
    x0 = [cos(la)*cos(lo), cos(la)*sin(lo), sin(la)]
    fine = ops_on .and. allocated(ops_B)
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      if (z_c(nb(1)) < box_z_km(1)*1.0e5_dp .or. z_c(1) > box_z_km(2)*1.0e5_dp) cycle
      do k = 1, nb(3)
        do j = 1, nb(2)
          e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
          x = e(:, 1)
          if (Key_IB(KP, ib) == 2) x = yy_cart(x)
          if (dot_product(x, x0) < cos(box_radius_deg*pi/180.0_dp)) cycle
          if (u == 0) then
            write(name, '(A,A,I4.4,A)') trim(out_dir), '/box_r', my_rank, '.dat'
            open(newunit=u, file=trim(name), status='replace', action='write')
            write(u, '(A)') '# step t[s] cell z[km] lat lon patch rho[g/cm3] T[K] vx vy vz[cm/s] Bx By Bz[G] (global Cartesian) ' // &
                'beta Mach Q_heat Q_loss[erg/cm3/s]  (box_fine_from, from that step on: + dt[s], the step''s change of s1 ' // &
                '[erg/g/K] by the RK hydro step, the transport, the conduction, the losses and heating, the diffusion, ' // &
                'the rest (floors, walls, GLM, friction), and of v_r [cm/s] by the RK step, the diffusion, the rest; then the ' // &
                'cell''s s1 and v_r (full precision))'
          end if
          do i = 1, nb(1)
            zz = z_c(i)/1.0e5_dp
            if (zz < box_z_km(1) .or. zz > box_z_km(2)) cycle
            rho = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            Tc = exp(out(F_LNT));  p = exp(out(F_LNP));  cs = sqrt(out(F_GAMMA1)*p/rho)
            v = matmul(e, State_GVB(i, j, k, vr_:vp_, ib))
            b = 0.0_dp
            if (use_B) b = matmul(e, State_GVB(i, j, k, br_:bp_, ib))
            if (Key_IB(KP, ib) == 2) then
              v = yy_cart(v);  b = yy_cart(b)
            end if
            b2 = sum(b**2)
            qh = 0.0_dp;  ql = 0.0_dp
            if (allocated(Qheat_B) .and. r_off + i > k_join) then
              qh = Qheat_B(i, j, k, ib);  ql = Qloss_B(i, j, k, ib)
            end if
            if (fine) then
              associate (o => ops_B(i, j, k, :, ib))
                ops = [o(2) - o(1), o(3) - o(2), ds1c_B(i, j, k, ib), o(4) - o(3) - ds1c_B(i, j, k, ib), o(5) - o(4), &
                    State_GVB(i, j, k, s1_, ib) - o(5), o(7) - o(6), o(8) - o(7), State_GVB(i, j, k, vr_, ib) - o(8)]
              end associate
              write(u, '(I9,ES15.7,I5,F11.1,2F9.3,I2,12ES12.4,ES11.3,9ES12.4,2ES24.15E3)') iStep, t, r_off + i, zz, &
                  asin(max(-1.0_dp, min(1.0_dp, x(3))))*180.0_dp/pi, atan2(x(2), x(1))*180.0_dp/pi, Key_IB(KP, ib), &
                  rho, Tc, v, b, p/max(b2/(8.0_dp*pi), tiny(1.0_dp)), norm2(v)/cs, qh, ql, dt_last, ops, &
                  State_GVB(i, j, k, s1_, ib), State_GVB(i, j, k, vr_, ib)
            else
              write(u, '(I9,ES15.7,I5,F11.1,2F9.3,I2,14ES12.4)') iStep, t, r_off + i, zz, &
                  asin(max(-1.0_dp, min(1.0_dp, x(3))))*180.0_dp/pi, atan2(x(2), x(1))*180.0_dp/pi, Key_IB(KP, ib), &
                  rho, Tc, v, b, p/max(b2/(8.0_dp*pi), tiny(1.0_dp)), norm2(v)/cs, qh, ql
            end if
          end do
        end do
      end do
    end do
    call cz_unbind()
    if (u /= 0) flush(u)
  end subroutine box_write

  ! Debugging probe (probe_radius_deg > 0): the top probe_cells cells of the
  ! columns near (probe_lat_deg, probe_lon_deg), every step, one file per rank:
  ! step t i patch level j k lat lon rho1/rho0 s1 v_r v_theta v_phi T c_s Mach
  subroutine probe_write(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    integer, save :: u = 0
    real(dp) :: x0(3), x(3), e(3,3), la, lo, out(n_field), rho, cs, Tc
    integer :: jb, ib, i, j, k
    character(len=300) :: name
    la = probe_lat_deg*pi/180.0_dp;  lo = probe_lon_deg*pi/180.0_dp
    x0 = [cos(la)*cos(lo), cos(la)*sin(lo), sin(la)]
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (.not. arbor_physical_side(ib, 1, 2)) cycle          ! top blocks only
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
          x = e(:, 1)
          if (Key_IB(KP, ib) == 2) x = yy_cart(x)
          if (dot_product(x, x0) < cos(probe_radius_deg*pi/180.0_dp)) cycle
          if (u == 0) then
            write(name, '(A,A,I4.4,A)') trim(out_dir), '/probe_r', my_rank, '.dat'
            open(newunit=u, file=trim(name), status='replace', action='write')
            write(u, '(A)') '# step t i patch level j k lat lon rho1/rho0 s1 v_r v_theta v_phi T c_s Mach'
          end if
          do i = max(1, nb(1) - probe_cells + 1), nb(1)
            rho = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            Tc = exp(out(F_LNT));  cs = sqrt(out(F_GAMMA1)*exp(out(F_LNP))/rho)
            write(u, '(I9,ES15.7,I5,4I3,2F9.3,9ES14.5)') iStep, t, r_off + i, Key_IB(KP, ib), Key_IB(KL+2, ib), j, k, &
                asin(max(-1.0_dp, min(1.0_dp, x(3))))*180.0_dp/pi, atan2(x(2), x(1))*180.0_dp/pi, &
                State_GVB(i, j, k, rho1_, ib)/rho0(i), State_GVB(i, j, k, s1_, ib), State_GVB(i, j, k, vr_:vp_, ib), &
                Tc, cs, norm2(State_GVB(i, j, k, vr_:vp_, ib))/cs
          end do
        end do
      end do
      call cz_unbind()
    end do
    if (u /= 0) flush(u)
  end subroutine probe_write

  ! Stop with a profile if the state holds a NaN
  subroutine nan_check(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    logical :: bad
    integer :: jb, ib, i, j, k, v
    real(dp) :: e(3,3), x(3)
    character(len=300) :: msg
    bad = .false.
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (.not. any(State_GVB(1:nb(1), 1:nb(2), 1:nb(3), :, ib) /= State_GVB(1:nb(1), 1:nb(2), 1:nb(3), :, ib))) cycle
      bad = .true.
      ! where: the first NaN cell of the block (radial index, depth, global lat/lon, level)
      call cz_bind(ib)
      cells: do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            do v = 1, nv
              if (State_GVB(i, j, k, v, ib) /= State_GVB(i, j, k, v, ib)) exit cells
            end do
          end do
        end do
      end do cells
      e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
      x = e(:, 1)
      if (Key_IB(KP, ib) == 2) x = yy_cart(x)
      write(msg, '(A,I0,A,I0,A,F9.1,A,F7.2,A,F8.2,A,I0,A,I0,A,2I4)') 'cz3d: NaN on rank ', my_rank, ': cell ', r_off + i, &
          ' (z ', z_c(i)/1.0e5_dp, ' km), lat ', asin(max(-1.0_dp, min(1.0_dp, x(3))))*180.0_dp/pi, ', lon ', &
          atan2(x(2), x(1))*180.0_dp/pi, ', patch ', Key_IB(KP, ib), ', level ', Key_IB(KL+2, ib), ', j k in block', j, k
      write(0, '(A)') trim(msg)
      call cz_unbind()
    end do
    if (.not. arbor_any(bad)) return
    call cz_log(iStep, t)
    call arbor_abort('cz3d: NaN in the state')
  end subroutine nan_check

  ! Horizontal sums per radial cell (weighted by solid angle; the Yin-Yang
  ! overlap counts twice, which does not matter for means of a uniform
  ! state): returns Sum_IV(i, v) for v = 1 weight, 2 rho1, 3 (rho1 - <rho1>)^2,
  ! 4 s1, 5 (s1 - <s1>)^2, 6 v_r, 7 v_r^2, 8 v_h^2, 9 T, 10 L_rad and
  ! 11 L_conv at face i, 12 rho v_r h1, 13 rho v_r, 14 h1 (h1 = T s1 + p1/rho,
  ! the enthalpy perturbation to first order), 15 rho v^2 v_r / 2
  ! 16 the diffusion's energy flux at face i, 17 rho |v_r| (the up/down mass exchange),
  ! 18 B_r^2, 19 B_theta^2 + B_phi^2 (use_B)
  ! (i = 0..n_r; the face values in row i, the cell values in rows 1..n_r).
  ! The deviations in a second pass, about the means of the first (a
  ! one-pass variance loses the digits). The resolved enthalpy flux is
  ! <rho v_r h1> - <rho v_r><h1> (columns 12-14).
  ! with_T is kept for the interface; T and p1 are always evaluated.
  subroutine horizontal_sums(Sum_IV, with_T, only_level)
    real(dp), intent(out) :: Sum_IV(0:n_r, 19)
    logical, intent(in) :: with_T
    integer, intent(in), optional :: only_level     ! only the blocks of this horizontal level (a refined patch)
    real(dp) :: w, out(n_field), rho, T, p1, h1
    real(dp), allocatable :: Dev_IV(:,:)
    integer :: jb, ib, i, j, k
    Sum_IV = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (present(only_level)) then
        if (Key_IB(KL+2, ib) /= only_level) cycle
      end if
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          ! (the faces between radial blocks from the block above)
          do i = merge(0, 1, at_bottom), nb(1)
            Sum_IV(r_off + i, 10) = Sum_IV(r_off + i, 10) + w*qrad_FB(i, j, k, ib)
            if (qconv_face) then
              Sum_IV(r_off + i, 11) = Sum_IV(r_off + i, 11) + w*FaceState_FVDB(i, j, k, 1, 1, ib)
            else
              Sum_IV(r_off + i, 11) = Sum_IV(r_off + i, 11) + w*qconv_FB(i, j, k, ib)
            end if
            Sum_IV(r_off + i, 16) = Sum_IV(r_off + i, 16) + w*qdif_FB(i, j, k, ib)
          end do
          do i = 1, nb(1)
            associate (s => State_GVB(i, j, k, :, ib), Sm => Sum_IV(r_off + i, :))
              Sm(1) = Sm(1) + w
              Sm(2) = Sm(2) + w*s(rho1_)
              Sm(4) = Sm(4) + w*s(s1_)
              Sm(6) = Sm(6) + w*s(vr_)
              Sm(7) = Sm(7) + w*s(vr_)**2
              Sm(8) = Sm(8) + w*(s(vt_)**2 + s(vp_)**2)
              rho = rho0(i) + s(rho1_)                                    ! inertia (mass flux) density
              call eos_eval(tab, rho0(i) + rq(i)*s(rho1_), s0(i) + s(s1_), out)
              T = exp(out(F_LNT))
              p1 = p0(i)*expm1(out(F_LNP) - lnp0(i))
              h1 = T*s(s1_) + p1/(rho0(i) + rq(i)*s(rho1_))
              Sm(9) = Sm(9) + w*T
              Sm(12) = Sm(12) + w*rho*s(vr_)*h1
              Sm(13) = Sm(13) + w*rho*s(vr_)
              Sm(14) = Sm(14) + w*h1
              Sm(15) = Sm(15) + w*0.5_dp*rho*sum(s(vr_:vp_)**2)*s(vr_)
              Sm(17) = Sm(17) + w*rho*abs(s(vr_))
              if (use_B) then
                Sm(18) = Sm(18) + w*s(br_)**2
                Sm(19) = Sm(19) + w*(s(bt_)**2 + s(bp_)**2)
              end if
            end associate
          end do
          if (at_bottom) Sum_IV(0, 1) = Sum_IV(0, 1) + w
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(Sum_IV)
    allocate(Dev_IV(0:n_r, 2))
    Dev_IV = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      if (present(only_level)) then
        if (Key_IB(KL+2, ib) /= only_level) cycle
      end if
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          do i = 1, nb(1)
            Dev_IV(r_off + i, 1) = Dev_IV(r_off + i, 1) + w*(State_GVB(i, j, k, rho1_, ib) - Sum_IV(r_off + i, 2)/Sum_IV(r_off + i, 1))**2
            Dev_IV(r_off + i, 2) = Dev_IV(r_off + i, 2) + w*(State_GVB(i, j, k, s1_, ib) - Sum_IV(r_off + i, 4)/Sum_IV(r_off + i, 1))**2
          end do
        end do
      end do
    end do
    call cz_unbind()
    call arbor_sum_inplace(Dev_IV)
    Sum_IV(:, 3) = Dev_IV(:, 1)
    Sum_IV(:, 5) = Dev_IV(:, 2)
  end subroutine horizontal_sums

  ! One log line: max Mach, max |v_r|, max |v_h|, the horizontally averaged
  ! luminosity (as CZTR: at the top face and its range below), and the
  ! largest horizontal rms deviations of s1 and rho1/rho0 (0 for a uniform state)
  module subroutine cz_log(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    real(dp) :: S(0:n_r, 19), L(0:n_r), Lres(1:n_r), mach, mach_r, dmin, vrmax, vhmax, cs, ds_h, dr_h, out(n_field), rho
    real(dp) :: bmax, emag, ekin, eratio, divb_rel, sd, sb, dxm, wl, wi, lz, tcor, vamax, n_cap, n_all, trmin, trmax
    integer :: nt                                   ! the top of the CZ (L_top, the L range): the join with the atmosphere
    integer :: jb, ib, i, j, k
    character(len=300) :: line
    mach = 0.0_dp;  mach_r = 0.0_dp;  vrmax = 0.0_dp;  vhmax = 0.0_dp;  dmin = huge(1.0_dp);  bmax = 0.0_dp
    vamax = 0.0_dp;  n_cap = 0.0_dp;  n_all = 0.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            rho = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rho, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            cs = sqrt(out(F_GAMMA1)*exp(out(F_LNP))/rho)
            mach = max(mach, norm2(State_GVB(i, j, k, vr_:vp_, ib))/cs)
            mach_r = max(mach_r, xi(i)*norm2(State_GVB(i, j, k, vr_:vp_, ib))/cs)
            dmin = min(dmin, State_GVB(i, j, k, rho1_, ib)/rho0(i))
            vrmax = max(vrmax, abs(State_GVB(i, j, k, vr_, ib)))
            vhmax = max(vhmax, hypot(State_GVB(i, j, k, vt_, ib), State_GVB(i, j, k, vp_, ib)))
            if (use_B) then
              bmax = max(bmax, norm2(State_GVB(i, j, k, br_:bp_, ib)))
              vamax = max(vamax, norm2(State_GVB(i, j, k, br_:bp_, ib))/sqrt(four_pi*rho))
              if (lorentz_f(sum(State_GVB(i, j, k, br_:bp_, ib)**2), rho, cell_cap(i, j, k, ib, cs)) < 0.9_dp) &
                  n_cap = n_cap + 1.0_dp
              n_all = n_all + 1.0_dp
            end if
          end do
        end do
      end do
    end do
    call cz_unbind()
    mach = arbor_max(mach);  mach_r = arbor_max(mach_r);  dmin = arbor_min(dmin);  vrmax = arbor_max(vrmax);  vhmax = arbor_max(vhmax)
    bmax = arbor_max(bmax);  vamax = arbor_max(vamax);  n_cap = arbor_sum(n_cap);  n_all = arbor_sum(n_all)
    nt = merge(k_join, n_r, use_atmosphere)
    call angular_momentum(lz, tcor)
    ! the transition region's height range (T = 1e5 K crossings below 10 Mm; -1 without)
    trmin = huge(1.0_dp);  trmax = -huge(1.0_dp)
    if (use_atmosphere) then
      do jb = 1, nBlock
        call cz_bind(iBlockList_I(jb))
        call atm_tr_heights(iBlockList_I(jb), 1.0e4_dp, trmin, trmax)
      end do
      call cz_unbind()
    end if
    trmin = arbor_min(trmin);  trmax = arbor_max(trmax)
    if (trmin > 1.0e30_dp) trmin = -1.0_dp
    if (trmax < -1.0e30_dp) trmax = -1.0_dp
    divb_rel = 0.0_dp;  wl = 0.0_dp;  wi = 0.0_dp
    if (use_B) then
      ! rms(div B dx_min) / rms|B| over all cells (ghosts filled first)
      call arbor_fill_ghosts(State_GVB, n_ghost=2, faces_only=.true.)
      sd = 0.0_dp;  sb = 0.0_dp
      do jb = 1, nBlock
        ib = iBlockList_I(jb)
        call cz_bind(ib)
        call div_block(State_GVB(:, :, :, br_:bp_, ib), ib, divB_C)
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, nb(1)
              dxm = min(drdxi_c(i), r_c(i)*Metric_DB(2, ib)%Dx_I(j), r_c(i)*sin(Metric_DB(2, ib)%X_I(j))*Metric_DB(3, ib)%Dx_I(k))
              sd = sd + (divB_C(i, j, k)*dxm)**2
              sb = sb + sum(State_GVB(i, j, k, br_:bp_, ib)**2)
            end do
          end do
        end do
      end do
      call cz_unbind()
      sd = arbor_sum(sd);  sb = arbor_sum(sb)
      divb_rel = sqrt(sd/max(sb, tiny(1.0_dp)))
      call mhd_exchange(wl, wi)
    end if
    call horizontal_sums(S, .false.)
    ! magnetic over kinetic energy (the shells' horizontal means; the kinetic with rho0)
    emag = 0.0_dp;  ekin = 0.0_dp
    do i = 1, n_r
      emag = emag + vol_c(i)*(S(i, 18) + S(i, 19))/S(i, 1)/(2.0_dp*four_pi)
      ekin = ekin + vol_c(i)*0.5_dp*rho0(i)*(S(i, 7) + S(i, 8))/S(i, 1)
    end do
    eratio = emag/max(ekin, tiny(1.0_dp))
    ! radiation, MLT, diffusion and the cooling layer at the faces
    L = 4.0_dp*pi*area_f*(S(:, 10) + S(:, 11) + S(:, 16) + S(0, 1)*cool_f)/S(0, 1)/Lsun
    call resolved_L(S, Lres)
    ds_h = 0.0_dp;  dr_h = 0.0_dp
    do i = 1, n_r
      ds_h = max(ds_h, sqrt(S(i, 5)/S(i, 1)))
      dr_h = max(dr_h, sqrt(S(i, 3)/S(i, 1))/rho0(i))
    end do
    if (.not. is_root) return
    if (u_log == 0 .and. restarted) then
      open(newunit=u_log, file=trim(out_dir)//'/log.txt', status='unknown', position='append', action='write')
      write(u_log, '(A,I0)') '# continued from a restart file at step ', iStep
      wall0 = arbor_wtime()
    else if (u_log == 0) then
      open(newunit=u_log, file=trim(out_dir)//'/log.txt', status='replace', action='write')
      write(u_log, '(A)') '# step  t[s]  dt[s]  max_Mach  max|v_r|  max|v_h| [cm/s]  L_top/Lsun  min_L  max_L  ' // &
          'max_rms_h(s1)[erg/g/K]  max_rms_h(rho1)/rho0  wall[s]  max_Mach(reduced c_s)  ' // &
          'max|L_enth+L_kin|  min_L_total  max_L_total   (L_total = L_rad + L_conv + L_diff + L_cool + L_enth + L_kin)' // &
          '  min(rho1/rho0)  max_tsc  t_local_bottom[s]  t_local_mid[s]  max|B|[G]  E_mag/E_kin  rms(divB dx)/rms|B|' // &
          '  W_Lorentz[Lsun]  dE_mag/dt(induction)[Lsun] (their sum is ~0 for walls without Poynting flux)' // &
          '  L_z[g cm2/s] (angular momentum of the flow about the global z axis, in the rotating frame)' // &
          '  T_Coriolis[g cm2/s2] (-2 Omega int rho s v_s dV: the Coriolis torque on L_z)' // &
          '  max_v_A[cm/s] (uncapped)  cap_frac (fraction of cells with the Lorentz factor f < 0.9, see va_max)' // &
          '  TR_min[km]  TR_max[km] (heights of the T = 1e5 K crossings, use_atmosphere)'
      wall0 = arbor_wtime()
    end if
    write(u_log, '(I9,ES14.6,ES11.3,3ES11.3,3F10.5,2ES11.3,F10.1,ES11.3,3F10.5,ES11.3,F9.3,2ES14.6,3ES11.3,2ES12.4,2ES14.6,2ES11.3,2F10.1)') iStep, t, &
        dt_last, mach, vrmax, vhmax, L(nt), minval(L(1:nt-1)), maxval(L(1:nt-1)), ds_h, dr_h, arbor_wtime() - wall0, &
        mach_r, maxval(abs(Lres(1:nt))), minval(L(1:nt-1) + Lres(1:nt-1)), maxval(L(1:nt-1) + Lres(1:nt-1)), dmin, maxval(tsc), &
        tloc(1), tloc(n_r/2), bmax, eratio, divb_rel, wl/Lsun, wi/Lsun, lz, tcor, vamax, n_cap/max(n_all, 1.0_dp), trmin, trmax
    flush(u_log)
    write(line, '(A,I8,A,ES11.4,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,F8.5,A,F8.5,A,F8.5,A,ES9.2)') 'step ', iStep, &
        ' t=', t, ' dt=', dt_last, ' Mach=', mach, ' (red. ', mach_r, ') |vr|=', vrmax, ' |vh|=', vhmax, ' L_top=', &
        L(nt), ' L ', minval(L(1:nt-1) + Lres(1:nt-1)), '..', maxval(L(1:nt-1) + Lres(1:nt-1)), &
        ' rms_h(s1)=', ds_h
    call arbor_print(line)
  end subroutine cz_log

  ! Horizontally averaged radial profile, out_dir/profile_NNNN.dat
  module subroutine cz_profile(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    real(dp) :: S(0:n_r, 19)
    character(len=300) :: name
    call horizontal_sums(S, .true.)
    write(name, '(A,A,I4.4,A)') trim(out_dir), '/profile_', i_prof, '.dat'
    call write_profile(S, name, 'horizontal means (solid-angle weighted) and rms deviations')
    if (patch_level > refine_h) then
      ! the same over the refined patch alone (its finest blocks)
      call horizontal_sums(S, .true., patch_level)
      write(name, '(A,A,I4.4,A)') trim(out_dir), '/profile_patch_', i_prof, '.dat'
      call write_profile(S, name, 'means over the refined patch (blocks of level patch_level) and rms deviations')
    end if
    if (omega_rot /= 0.0_dp) call rotation_profile(iStep, t)
    i_prof = i_prof + 1
    t_next_prof = t_prof0 + i_prof*dt_prof
  contains
    subroutine write_profile(S, name, what)
      real(dp), intent(in) :: S(0:n_r, 19)
      character(*), intent(in) :: name, what
      real(dp) :: w, m, Lres(1:n_r)
      integer :: u, i
      call resolved_L(S, Lres)
      if (.not. is_root) return
      open(newunit=u, file=trim(name), status='replace', action='write')
      write(u, '(A,ES14.6,A,I0,A,A)') '# t = ', t, ' s; step ', iStep, '; ', what
      write(u, '(A)') '# i  r[cm]  z[km]  <rho>  <T>  <s1>  rms_h(s1)  <v_r>  rms(v_r)  rms(v_h)  ' // &
          'L_rad/Lsun  L_conv/Lsun  L_enth+L_kin/Lsun  L_kin/Lsun  L_diff/Lsun  L_cool/Lsun   (L_rad, L_conv, ' // &
          'L_diff (numerical diffusion), L_cool (cooling layer): through the face above the cell; L_enth, L_kin: ' // &
          'resolved, at the cell centre)  t_local[s]  tsc (the shell''s own clock and step factor; dt1 > 0)  ' // &
          '<rho v_r>  <rho |v_r|> [g/cm2/s] (the mean mass flux, the up/down exchange)  rms(B_r)  rms(B_h) [G]'
      do i = 1, n_r
        w = S(i, 1)
        m = S(i, 4)/w
        write(u, '(I5,ES20.12,20ES16.7E3)') i, r_c(i), z_c(i)/1.0e5_dp, rho0(i) + rq(i)*S(i, 2)/w, S(i, 9)/w, m, &
            sqrt(S(i, 5)/w), S(i, 6)/w, sqrt(S(i, 7)/w), sqrt(S(i, 8)/w), &
            4.0_dp*pi*area_f(i)*S(i, 10)/S(0, 1)/Lsun, 4.0_dp*pi*area_f(i)*S(i, 11)/S(0, 1)/Lsun, Lres(i), &
            4.0_dp*pi*r_c(i)**2*S(i, 15)/w/Lsun, 4.0_dp*pi*area_f(i)*S(i, 16)/S(0, 1)/Lsun, &
            4.0_dp*pi*area_f(i)*cool_f(i)/Lsun, tloc(i), tsc(i), S(i, 13)/w, S(i, 17)/w, sqrt(S(i, 18)/w), sqrt(S(i, 19)/w)
      end do
      close(u)
    end subroutine write_profile
  end subroutine cz_profile

  ! Resolved luminosity at the cell centres: 4 pi r^2 (<rho v_r h1> -
  ! <rho v_r><h1> + <rho v^2 v_r/2>), in Lsun
  subroutine resolved_L(S, Lres)
    real(dp), intent(in) :: S(0:n_r, 19)
    real(dp), intent(out) :: Lres(1:n_r)
    integer :: i
    real(dp) :: w
    do i = 1, n_r
      w = S(i, 1)
      Lres(i) = 4.0_dp*pi*r_c(i)**2*(S(i, 12)/w - (S(i, 13)/w)*(S(i, 14)/w) + S(i, 15)/w)/Lsun
    end do
  end subroutine resolved_L

  ! Every cell of the radial indices dump_shells, one file per rank
  ! (out_dir/shell_dump_rNNNN.dat): i, patch, the position as a global unit
  ! vector (x, y, z; Yin = global frame, Yang mapped by (x,y,z) -> (-x,z,y)),
  ! the overlap weight (solid angle, 1/2 where both patches' cores cover the
  ! point, 0 in a margin), v_r, the horizontal velocity as a global Cartesian
  ! vector, rho1/rho0, s1
  module subroutine dump_shells_write()
    real(dp) :: th, ph, th2, ph2, e(3,3), x(3), vh(3), w, bvec(3)
    integer :: u, jb, ib, j, k, n, i
    character(len=300) :: name
    if (all(dump_shells <= 0)) return
    write(name, '(A,A,I4.4,A)') trim(out_dir), '/shell_dump_r', my_rank, '.dat'
    open(newunit=u, file=trim(name), status='replace', action='write')
    write(u, '(A)') '# i  patch  x  y  z  weight  v_r  vh_x  vh_y  vh_z [cm/s]  rho1/rho0  s1  (use_B: B_x  B_y  B_z [G])'
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          th = Metric_DB(2, ib)%X_I(j)
          ph = Metric_DB(3, ib)%X_I(k)
          e = sph_basis(th, ph)
          w = Metric_DB(2, ib)%Vol_I(j)*Metric_DB(3, ib)%Vol_I(k)
          if (.not. core(th, ph)) then
            w = 0.0_dp
          else
            call yy_angles(th, ph, th2, ph2)
            if (core(th2, ph2)) w = 0.5_dp*w
          end if
          do n = 1, size(dump_shells)
            i = dump_shells(n) - r_off
            if (i < 1 .or. i > nb(1)) cycle
            x = e(:, 1)
            vh = State_GVB(i, j, k, vt_, ib)*e(:, 2) + State_GVB(i, j, k, vp_, ib)*e(:, 3)
            if (Key_IB(KP, ib) == 2) then
              x = yy_cart(x)
              vh = yy_cart(vh)
            end if
            if (use_B) then
              bvec = State_GVB(i, j, k, br_, ib)*e(:, 1) + State_GVB(i, j, k, bt_, ib)*e(:, 2) &
                  + State_GVB(i, j, k, bp_, ib)*e(:, 3)
              if (Key_IB(KP, ib) == 2) bvec = yy_cart(bvec)
              write(u, '(2I5,3F12.8,ES14.6,4ES14.6,2ES14.6,3ES14.6)') r_off + i, Key_IB(KP, ib), x, w, &
                  State_GVB(i, j, k, vr_, ib), vh, State_GVB(i, j, k, rho1_, ib)/rho0(i), State_GVB(i, j, k, s1_, ib), bvec
            else
              write(u, '(2I5,3F12.8,ES14.6,4ES14.6,2ES14.6)') r_off + i, Key_IB(KP, ib), x, w, State_GVB(i, j, k, vr_, ib), vh, &
                  State_GVB(i, j, k, rho1_, ib)/rho0(i), State_GVB(i, j, k, s1_, ib)
            end if
          end do
        end do
      end do
    end do
    call cz_unbind()
    close(u)
  contains
    logical function core(t, p)
      real(dp), intent(in) :: t, p
      core = t >= 0.25_dp*pi .and. t <= 0.75_dp*pi .and. abs(p) <= 0.75_dp*pi
    end function core
  end subroutine dump_shells_write

  ! Order-independent checksum of the interior state of all blocks, to
  ! compare runs (e.g. continued vs uninterrupted) bitwise
  ! (the result typed in the body: ifx does not resolve the module's kind i8
  ! in a submodule function's prefix)
  module function state_checksum() result(c)
    integer(i8) :: c
    integer(i8) :: x
    integer :: jb, ib, i, j, k, v
    c = 0
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      x = int(iLeaf_B(ib), i8)
      do v = 1, nv
        do k = 1, nb(3)
          do j = 1, nb(2)
            do i = 1, nb(1)
              x = ieor(x, transfer(State_GVB(i, j, k, v, ib), 1_i8))
              x = ieor(x, ishft(x, 13));  x = ieor(x, ishft(x, -7));  x = ieor(x, ishft(x, 17))
            end do
          end do
        end do
      end do
      c = c + iand(x, 2_i8**40 - 1)
    end do
    c = arbor_sum(c)
  end function state_checksum

  ! The cell of the largest Mach number |v|/c_s: printed (radius, global
  ! latitude and longitude, z, T, rho, v_r, v_h, c_s, |B|) by the rank that holds it
  module subroutine max_mach_report(iStep, t)
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    real(dp) :: out(n_field), rhop, cs, m, mbest, gbest, e(3,3), x(3), loc(8)
    integer :: jb, ib, i, j, k
    character(len=400) :: msg
    mbest = -1.0_dp
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, nb(1)
            rhop = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rhop, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            cs = sqrt(out(F_GAMMA1)*exp(out(F_LNP))/rhop)
            m = norm2(State_GVB(i, j, k, vr_:vp_, ib))/cs
            if (m > mbest) then
              mbest = m
              e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
              x = e(:, 1)
              if (Key_IB(KP, ib) == 2) x = yy_cart(x)
              loc = [r_c(i), asin(max(-1.0_dp, min(1.0_dp, x(3))))*180.0_dp/pi, atan2(x(2), x(1))*180.0_dp/pi, &
                  exp(out(F_LNT)), rhop, State_GVB(i, j, k, vr_, ib), hypot(State_GVB(i, j, k, vt_, ib), &
                  State_GVB(i, j, k, vp_, ib)), cs]
              if (use_B) then
                write(msg, '(A,F9.4,A,F7.2,A,F8.2,A,F10.1,A,ES10.3,A,ES10.3,A,2F9.2,A,F8.2,A,F9.2,A)') 'max Mach: r/R ', &
                    loc(1)/6.957e10_dp, '  lat ', loc(2), '  lon ', loc(3), '  z[km] ', z_c(i)/1.0e5_dp, '  T ', loc(4), &
                    '  rho ', loc(5), '  v_r, v_h [km/s] ', loc(6)/1.0e5_dp, loc(7)/1.0e5_dp, '  c_s ', loc(8)/1.0e5_dp, &
                    '  |B| [G] ', norm2(State_GVB(i, j, k, br_:bp_, ib)), ''
              else
                write(msg, '(A,F9.4,A,F7.2,A,F8.2,A,F10.1,A,ES10.3,A,ES10.3,A,2F9.2,A,F8.2)') 'max Mach: r/R ', &
                    loc(1)/6.957e10_dp, '  lat ', loc(2), '  lon ', loc(3), '  z[km] ', z_c(i)/1.0e5_dp, '  T ', loc(4), &
                    '  rho ', loc(5), '  v_r, v_h [km/s] ', loc(6)/1.0e5_dp, loc(7)/1.0e5_dp, '  c_s ', loc(8)/1.0e5_dp
              end if
            end if
          end do
        end do
      end do
    end do
    call cz_unbind()
    gbest = arbor_max(mbest)
    if (mbest == gbest) then
      write(0, '(A,I0,A,ES12.5,A,F8.2,A,A)') 'step ', iStep, '  t ', t, '  Mach ', mbest, '  ', trim(msg)
      flush(0)
    end if
  end subroutine max_mach_report

  ! Tecplot file of the state in out_dir: rho1, v_x, v_y, v_z, s1, [B_x, B_y,
  ! B_z,] T - T0, rho, T and the Mach number |v|/c_s (cell centres; ARBOR
  ! writes node values). The vectors as Cartesian components in the global
  ! (Yin) frame, the frame of the X, Y, Z coordinates (Yang's rotated into it).
  ! The plot quantities sit in slots beyond 8 and only those are passed to the
  ! writer: its ghost fill then treats them as scalars (it rotates the slots
  ! declared as spherical vectors, v at vr_ and B at br_, which hold nothing here)
  module subroutine cz_write_plot(self, dom, iStep, t)
    class(cz3d_app_t), intent(inout) :: self
    type(arbor_domain_t), intent(in) :: dom
    integer, intent(in) :: iStep
    real(dp), intent(in) :: t
    integer, parameter :: o = 8
    character(len=300) :: name
    real(dp), allocatable :: P_GVB(:,:,:,:,:)
    real(dp) :: out(n_field), rhop, e(3,3), vc(3), bc(3)
    integer :: jb, ib, i, j, k, n, q
    n = merge(12, 9, use_B)
    allocate(P_GVB(1-ng:nb(1)+ng, 1-ng:nb(2)+ng, 1-ng:nb(3)+ng, o + n, nBlockMax))
    P_GVB = 0.0_dp
    q = merge(3, 0, use_B)                                ! the scalars after the field
    do jb = 1, nBlock
      ib = iBlockList_I(jb)
      call cz_bind(ib)
      do k = 1, nb(3)
        do j = 1, nb(2)
          e = sph_basis(Metric_DB(2, ib)%X_I(j), Metric_DB(3, ib)%X_I(k))
          do i = 1, nb(1)
            vc = matmul(e, State_GVB(i, j, k, vr_:vp_, ib))  ! Cartesian in the patch frame
            if (Key_IB(KP, ib) == 2) vc = yy_cart(vc)       ! Yang -> the global frame
            P_GVB(i, j, k, o+1, ib) = State_GVB(i, j, k, rho1_, ib)
            P_GVB(i, j, k, o+2:o+4, ib) = vc
            P_GVB(i, j, k, o+5, ib) = State_GVB(i, j, k, s1_, ib)
            if (use_B) then
              bc = matmul(e, State_GVB(i, j, k, br_:bp_, ib))
              if (Key_IB(KP, ib) == 2) bc = yy_cart(bc)
              P_GVB(i, j, k, o+6:o+8, ib) = bc
            end if
            rhop = rho0(i) + rq(i)*State_GVB(i, j, k, rho1_, ib)
            call eos_eval(tab, rhop, s0(i) + State_GVB(i, j, k, s1_, ib), out)
            P_GVB(i, j, k, o+q+6, ib) = exp(out(F_LNT)) - T0(i)
            P_GVB(i, j, k, o+q+7, ib) = rhop
            P_GVB(i, j, k, o+q+8, ib) = exp(out(F_LNT))
            P_GVB(i, j, k, o+q+9, ib) = norm2(State_GVB(i, j, k, vr_:vp_, ib))/sqrt(out(F_GAMMA1)*exp(out(F_LNP))/rhop)
          end do
        end do
      end do
    end do
    call cz_unbind()
    write(name, '(A,A,I8.8,A)') trim(out_dir), '/plot_', iStep, '.plt'
    if (use_B) then
      call arbor_plot_tecplot(dom, name, P_GVB, iVarMin=o+1, iVarMax=o+n, names=['rho1', 'vx  ', 'vy  ', 'vz  ', &
          's1  ', 'Bx  ', 'By  ', 'Bz  ', 'T1  ', 'rho ', 'T   ', 'Mach'], t=t, title='cztr3d')
    else
      call arbor_plot_tecplot(dom, name, P_GVB, iVarMin=o+1, iVarMax=o+n, names=['rho1', 'vx  ', 'vy  ', 'vz  ', &
          's1  ', 'T1  ', 'rho ', 'T   ', 'Mach'], t=t, title='cztr3d')
    end if
  end subroutine cz_write_plot

end submodule cz3d_diag_smod
