module cztr3d_atm_mod
  ! The atmosphere physics of CZTR (cztr_physics_mod's atm_step and
  ! cztr_advance_mod's floors), ported from CZTR_2D's cztr2d_atm_mod to the
  ! 3D Yin-Yang shell with radial blocks and multirate (use_atmosphere, see
  ! cz3d_params_mod):
  !   - hyperbolic Spitzer conduction (HYDRAD's H I term above, Rempel
  !     saturation 1/(1 + |F|/(1.5 rho c_s^3))) at the faces in r, theta and
  !     phi, times the switch S at the face; each face flux relaxes to its
  !     target in tau = max(D/budget^2, 4 dt), budget = f_cfl_cond dx/dt - |v|,
  !     as CZTR. The relaxing fluxes are kept per block (qr_B, qt_B, qp_B); a
  !     face shared by two blocks has a copy in each, computed from the same
  !     temperatures through the ghost cells (CZTR_2D's rule for theta faces,
  !     now in every direction; with multirate the two copies relax with their
  !     own steps). No flux through the bottom and top walls. The heat goes in
  !     as CZTR's add_heat, ds = cv ln(1 + dQ/(rho cv T)), dQ = -dt div q
  !     (ARBOR's flux divergence: the exact areas);
  !   - Carlsson & Leenaarts + coronal losses (blended at 2.4e4 K): the column
  !     mass above a cell from p/g (hydrostatic estimate, no column sums
  !     across blocks), the neutral-H column integrated from the block's top
  !     down plus x_n n_H H_p of the cell above the block (its ghost); times S;
  !   - the heating S H_peak x e^(1-x), x = max(z, 0)/z_heat;
  !   - floors: rho >= rho_floor everywhere, T >= T_floor where S > 0.5.
  ! Losses and heating use the state before the conduction (as CZTR). The
  ! routines work on the bound block (cz_bind: the radial arrays are views
  ! of its cells, indexed 1..nb(1) with ghosts).
  use arbor
  use cz3d_params_mod
  use cz3d_bg_mod, only: n_r, k_join, r_c, r_f, z_c, rho0, s0, p0, g_c, tab, ddr_c, r_top_cz
  use radiative_loss_mod
  use eos_table, only: eos_fields_n, eos_eval, eos_s_from_rho_T, n_field, F_LNP, F_LNT, F_GAMMA1, F_LNCV, F_LNNE
  implicit none
  private

  real(dp), parameter :: kappa0_spitzer = 7.8e-7_dp, T_SWITCH = 2.4e4_dp, T_SWITCH_WIDTH = 3.0e3_dp
  real(dp), parameter :: Rsun_atm = 6.957e10_dp, kB_atm = 1.380649e-16_dp
  integer, parameter :: RHO1 = 1, VR = 2, VT = 3, VP = 4, S1 = 5, BR = 6, BT = 7, BP = 8

  type(species_table_t) :: h1_tab, ca2_tab, mg2_tab
  type(coronal_table_t) :: coronal_tab
  type(hconduct_table_t) :: hcond_tab

  ! the background-force switch (cells 1:n_r of the whole column: 0 up to
  ! k_join - 2, 1 above) and the background force dp0/dr + rho0 g
  real(dp), allocatable, public :: S_bgf(:), bgforce(:)
  ! relaxing conductive fluxes (faces in r, theta, phi) and the last losses per block slot
  real(dp), allocatable :: qr_B(:,:,:,:), qt_B(:,:,:,:), qp_B(:,:,:,:), Qloss_B(:,:,:,:)
  ! one-block scratch: cells 0:nb(1)+1, 0:nb(2)+1, 0:nb(3)+1 (the faces' ghost cross only)
  real(dp), allocatable :: rho_A(:,:,:), T_A(:,:,:), cv_A(:,:,:), cs_A(:,:,:), ne_A(:,:,:), p_A(:,:,:)
  real(dp), allocatable :: F1(:,:,:,:), F2(:,:,:,:), F3(:,:,:,:), Dv(:,:,:,:)

  public :: atm_init, atm_setup, atm_block, atm_floors, atm_dt, atm_tr_heights, atm_switch

contains

  ! The switch S at height z [cm] above the CZ top
  elemental real(dp) function atm_switch(z) result(s)
    real(dp), intent(in) :: z
    s = 0.5_dp*(1.0_dp + tanh((z - z_atm_km*1.0e5_dp)/(w_atm_km*1.0e5_dp)))
  end function atm_switch

  ! Once, after bg_load (whole column): tables, background force
  subroutine atm_init()
    real(dp) :: aH, aCa, aMg
    call load_abundance_mixture(trim(atomic_dir)//'/abundances/asplund.ab', trim(atomic_dir)//'/masses/masses.amu', &
        aH, aCa, aMg)
    call load_species_table('h_1', trim(atomic_dir)//'/OpticallyThick', aH, h1_tab)
    call load_species_table('ca_2', trim(atomic_dir)//'/OpticallyThick', aCa, ca2_tab)
    call load_species_table('mg_2', trim(atomic_dir)//'/OpticallyThick', aMg, mg2_tab)
    call load_coronal_table(coronal_file, coronal_tab)
    call load_hconduct_table(trim(atomic_dir)//'/OpticallyThick/thermal_conductivities/h_1.tc', hcond_tab)
    allocate(S_bgf(1:n_r), bgforce(1:n_r))
    ! left out (well balanced) only where the radial stencil (+-2 cells) sees
    ! the CZ background alone (CZTR_2D: in the last two CZ cells the left-out
    ! force was +1.0% and +2.6% of rho0 g, the release kick of 2026-10-30)
    S_bgf = 1.0_dp
    S_bgf(1:min(k_join - 2, n_r)) = 0.0_dp
    bgforce = ddr_c(p0) + rho0(1:n_r)*g_c(1:n_r)
  end subroutine atm_init

  ! After the blocks are (re)built: the per-block arrays (fluxes from 0) and scratch
  subroutine atm_setup()
    if (allocated(qr_B)) deallocate(qr_B, qt_B, qp_B, Qloss_B)
    allocate(qr_B(0:nb(1), nb(2), nb(3), nBlockMax), qt_B(nb(1), 0:nb(2), nb(3), nBlockMax), &
        qp_B(nb(1), nb(2), 0:nb(3), nBlockMax), Qloss_B(nb(1), nb(2), nb(3), nBlockMax))
    qr_B = 0.0_dp;  qt_B = 0.0_dp;  qp_B = 0.0_dp;  Qloss_B = 0.0_dp
    if (allocated(rho_A)) return
    allocate(rho_A(0:nb(1)+1, 0:nb(2)+1, 0:nb(3)+1))
    allocate(T_A, cv_A, cs_A, ne_A, p_A, mold=rho_A)
    rho_A = 1.0_dp;  T_A = 1.0_dp;  cv_A = 1.0_dp;  cs_A = 1.0_dp;  ne_A = 0.0_dp;  p_A = 1.0_dp
    allocate(F1(0:nb(1), nb(2), nb(3), 1), F2(nb(1), 0:nb(2), nb(3), 1), F3(nb(1), nb(2), 0:nb(3), 1), &
        Dv(nb(1), nb(2), nb(3), 1))
  end subroutine atm_setup

  ! The state of one radial run of cells i0..i1 at (j, k) of block ib
  subroutine column_states(ib, i0, i1, j, k)
    integer, intent(in) :: ib, i0, i1, j, k
    real(dp) :: rp(i0:i1), s(i0:i1), q(i0:i1, 5)
    integer :: n
    n = i1 - i0 + 1
    rp = rho0(i0:i1) + State_GVB(i0:i1, j, k, RHO1, ib)
    s = s0(i0:i1) + State_GVB(i0:i1, j, k, S1, ib)
    call eos_fields_n(tab, [F_LNT, F_LNCV, F_LNNE, F_LNP, F_GAMMA1], n, rp, s, q)
    rho_A(i0:i1, j, k) = rp
    T_A(i0:i1, j, k) = exp(q(:, 1))
    cv_A(i0:i1, j, k) = exp(q(:, 2))
    ne_A(i0:i1, j, k) = exp(q(:, 3))
    p_A(i0:i1, j, k) = exp(q(:, 4))
    cs_A(i0:i1, j, k) = sqrt(q(:, 5)*p_A(i0:i1, j, k)/rp)
  end subroutine column_states

  ! The cell states of block ib: the interior with its radial ghosts
  ! (0, nb(1)+1), and the theta and phi ghost layers next to the interior
  ! (the corners are not used and may be stale)
  subroutine states(ib)
    integer, intent(in) :: ib
    integer :: j, k
    do k = 1, nb(3)
      do j = 1, nb(2)
        call column_states(ib, 0, nb(1) + 1, j, k)
      end do
      call column_states(ib, 1, nb(1), 0, k)
      call column_states(ib, 1, nb(1), nb(2) + 1, k)
    end do
    do j = 1, nb(2)
      call column_states(ib, 1, nb(1), j, 0)
      call column_states(ib, 1, nb(1), j, nb(3) + 1)
    end do
  end subroutine states

  ! Conduction, losses and heating of block ib over dt (ghosts filled;
  ! at_bot/at_top: the block holds the bottom/top wall of the domain);
  ! stf(0:nb(1)+1): the step factors of the cells (local time stepping; with
  ! the cells below and above the block): cell i advances stf(i) dt, a radial
  ! face the smaller step of its two cells
  subroutine atm_block(ib, dt_in, at_bot, at_top, stf)
    integer, intent(in) :: ib
    real(dp), intent(in) :: dt_in
    logical, intent(in) :: at_bot, at_top
    real(dp), intent(in) :: stf(0:)
    real(dp) :: Tf, rf, Tg, csf, cvf, dxf, vf, x, dQ, Sc, Sf, sth, dt, fS, qcl, gt(3), bf(3), aniso
    real(dp) :: colH(nb(1)), dr, qd, ionH, lT, w, nHt, QH, iH, QCa, iCa, QMg, iMg, chromo, coron, Qheat, Qlos, colm
    integer :: i, j, k, m
    m = nb(1)
    if (maxval(atm_switch(z_c(1:m))) < 1.0e-8_dp .and. atm_switch(r_f(m) - r_top_cz) < 1.0e-8_dp) return  ! (deep CZ block)
    call states(ib)
    aniso = merge(1.0_dp, 0.0_dp, cond_along_b .and. use_B)
    associate (X2 => Metric_DB(2, ib)%X_I, X3 => Metric_DB(3, ib)%X_I)
    if (use_conduction) then
      F1 = 0.0_dp;  F2 = 0.0_dp;  F3 = 0.0_dp
      do k = 1, nb(3)
        do j = 1, nb(2)
          ! radial faces 0..m (none through the domain's bottom and top walls)
          do i = 0, m
            Sf = atm_switch(r_f(i) - r_top_cz)
            if (Sf < 1.0e-8_dp .or. (i == 0 .and. at_bot) .or. (i == m .and. at_top)) then
              qr_B(i, j, k, ib) = 0.0_dp
              ! an outflow top: the top cell passes the conductive energy on
              ! (the same F A through its top face as through the face below)
              if (i == m .and. at_top .and. wall_top == 'outflow') then
                qr_B(m, j, k, ib) = qr_B(m-1, j, k, ib)*(r_f(m-1)/r_f(m))**2
                F1(m, j, k, 1) = qr_B(m, j, k, ib)
              end if
              cycle
            end if
            dt = dt_in*min(stf(i), stf(i+1))
            dxf = r_c(i+1) - r_c(i)
            Tf = 0.5_dp*(T_A(i, j, k) + T_A(i+1, j, k))
            rf = 0.5_dp*(rho_A(i, j, k) + rho_A(i+1, j, k))
            csf = 0.5_dp*(cs_A(i, j, k) + cs_A(i+1, j, k))
            cvf = 0.5_dp*(cv_A(i, j, k) + cv_A(i+1, j, k))
            Tg = (T_A(i+1, j, k) - T_A(i, j, k))/dxf
            if (aniso > 0.0_dp) then
              ! along the field: the face's (r, theta, phi) gradient (tangential
              ! parts from its interior cells' central differences)
              gt = 0.0_dp;  x = 0.0_dp
              if (i >= 1) then;  gt = gt + grad_cell(i, j, k);  x = x + 1.0_dp;  end if
              if (i + 1 <= m) then;  gt = gt + grad_cell(i+1, j, k);  x = x + 1.0_dp;  end if
              gt = gt/max(x, 1.0_dp);  gt(1) = Tg
              bf = 0.5_dp*(State_GVB(i, j, k, BR:BP, ib) + State_GVB(i+1, j, k, BR:BP, ib))
              Tg = along_b(gt, bf, 1)
            end if
            vf = 0.5_dp*abs(State_GVB(i, j, k, VR, ib) + State_GVB(i+1, j, k, VR, ib))
            fS = 1.0_dp;  qcl = 0.0_dp
            if (cond_r_rsun > 0.0_dp) then
              ! AWSoM: Spitzer -> the collisionless flux (3/2) alpha p_e v_r far out
              fS = 1.0_dp/(1.0_dp + (r_f(i)/(cond_r_rsun*Rsun_atm))**2)
              qcl = 1.5_dp*alpha_cl*kB_atm*Tf*0.5_dp*(ne_A(i, j, k) + ne_A(i+1, j, k)) &
                  *0.5_dp*(State_GVB(i, j, k, VR, ib) + State_GVB(i+1, j, k, VR, ib))
            end if
            call face_flux(Sf, Tf, rf, csf, cvf, Tg, dxf, vf, qr_B(i, j, k, ib), fS, qcl)
            F1(i, j, k, 1) = qr_B(i, j, k, ib)
          end do
        end do
        ! theta faces 0..nb(2) (the switch of the cell's height)
        do j = 0, nb(2)
          dxf = X2(j+1) - X2(j)
          do i = 1, m
            Sc = atm_switch(z_c(i))
            if (Sc < 1.0e-8_dp) then
              qt_B(i, j, k, ib) = 0.0_dp
              cycle
            end if
            dt = dt_in*stf(i)
            Tf = 0.5_dp*(T_A(i, j, k) + T_A(i, j+1, k))
            rf = 0.5_dp*(rho_A(i, j, k) + rho_A(i, j+1, k))
            csf = 0.5_dp*(cs_A(i, j, k) + cs_A(i, j+1, k))
            cvf = 0.5_dp*(cv_A(i, j, k) + cv_A(i, j+1, k))
            Tg = (T_A(i, j+1, k) - T_A(i, j, k))/(r_c(i)*dxf)
            if (aniso > 0.0_dp) then
              gt = 0.0_dp;  x = 0.0_dp
              if (j >= 1) then;  gt = gt + grad_cell(i, j, k);  x = x + 1.0_dp;  end if
              if (j + 1 <= nb(2)) then;  gt = gt + grad_cell(i, j+1, k);  x = x + 1.0_dp;  end if
              gt = gt/max(x, 1.0_dp);  gt(2) = Tg
              bf = 0.5_dp*(State_GVB(i, j, k, BR:BP, ib) + State_GVB(i, j+1, k, BR:BP, ib))
              Tg = along_b(gt, bf, 2)
            end if
            vf = 0.5_dp*abs(State_GVB(i, j, k, VT, ib) + State_GVB(i, j+1, k, VT, ib))
            fS = 1.0_dp
            if (cond_r_rsun > 0.0_dp) fS = 1.0_dp/(1.0_dp + (r_c(i)/(cond_r_rsun*Rsun_atm))**2)
            call face_flux(Sc, Tf, rf, csf, cvf, Tg, r_c(i)*dxf, vf, qt_B(i, j, k, ib), fS, 0.0_dp)
            F2(i, j, k, 1) = qt_B(i, j, k, ib)
          end do
        end do
      end do
      ! phi faces 0..nb(3)
      do k = 0, nb(3)
        dxf = X3(k+1) - X3(k)
        do j = 1, nb(2)
          sth = sin(X2(j))
          do i = 1, m
            Sc = atm_switch(z_c(i))
            if (Sc < 1.0e-8_dp) then
              qp_B(i, j, k, ib) = 0.0_dp
              cycle
            end if
            dt = dt_in*stf(i)
            Tf = 0.5_dp*(T_A(i, j, k) + T_A(i, j, k+1))
            rf = 0.5_dp*(rho_A(i, j, k) + rho_A(i, j, k+1))
            csf = 0.5_dp*(cs_A(i, j, k) + cs_A(i, j, k+1))
            cvf = 0.5_dp*(cv_A(i, j, k) + cv_A(i, j, k+1))
            Tg = (T_A(i, j, k+1) - T_A(i, j, k))/(r_c(i)*sth*dxf)
            if (aniso > 0.0_dp) then
              gt = 0.0_dp;  x = 0.0_dp
              if (k >= 1) then;  gt = gt + grad_cell(i, j, k);  x = x + 1.0_dp;  end if
              if (k + 1 <= nb(3)) then;  gt = gt + grad_cell(i, j, k+1);  x = x + 1.0_dp;  end if
              gt = gt/max(x, 1.0_dp);  gt(3) = Tg
              bf = 0.5_dp*(State_GVB(i, j, k, BR:BP, ib) + State_GVB(i, j, k+1, BR:BP, ib))
              Tg = along_b(gt, bf, 3)
            end if
            vf = 0.5_dp*abs(State_GVB(i, j, k, VP, ib) + State_GVB(i, j, k+1, VP, ib))
            fS = 1.0_dp
            if (cond_r_rsun > 0.0_dp) fS = 1.0_dp/(1.0_dp + (r_c(i)/(cond_r_rsun*Rsun_atm))**2)
            call face_flux(Sc, Tf, rf, csf, cvf, Tg, r_c(i)*sth*dxf, vf, qp_B(i, j, k, ib), fS, 0.0_dp)
            F3(i, j, k, 1) = qp_B(i, j, k, ib)
          end do
        end do
      end do
      call arbor_flux_divergence(ib, F1, F2, F3, Dv)
      do k = 1, nb(3)
        do j = 1, nb(2)
          do i = 1, m
            dQ = -dt_in*stf(i)*Dv(i, j, k, 1)
            if (dQ == 0.0_dp) cycle
            x = dQ/(rho_A(i, j, k)*cv_A(i, j, k)*T_A(i, j, k))
            State_GVB(i, j, k, S1, ib) = State_GVB(i, j, k, S1, ib) + cv_A(i, j, k)*log(1.0_dp + max(x, -0.9_dp))
          end do
        end do
      end do
    end if

    ! losses and heating (the state before the conduction)
    do k = 1, nb(3)
      do j = 1, nb(2)
        if (use_losses) then
          ! the neutral-H column: from above the block (the ghost cell's
          ! x_n n_H H_p; 0 at the top wall), then down through the block
          colH = 0.0_dp
          if (.not. at_top) then
            call species_emiss_ionfrac(h1_tab, log10(T_A(m+1, j, k)), qd, ionH)
            colH(m) = ionH*N_H_per_gram_global*rho_A(m+1, j, k)*p_A(m+1, j, k)/(rho_A(m+1, j, k)*g_c(m+1))
          end if
          do i = m, 1, -1
            dr = r_f(i) - r_f(i-1)
            call species_emiss_ionfrac(h1_tab, log10(T_A(i, j, k)), qd, ionH)
            if (i < m) colH(i) = colH(i+1)
            colH(i) = colH(i) + ionH*N_H_per_gram_global*rho_A(i, j, k)*dr
          end do
        end if
        do i = 1, m
          Qlos = 0.0_dp
          Qheat = 0.0_dp
          Sc = atm_switch(z_c(i))
          if (Sc >= 1.0e-8_dp) then
            if (use_losses) then
              lT = log10(T_A(i, j, k))
              w = blend_weight(T_A(i, j, k), T_SWITCH, T_SWITCH_WIDTH)
              nHt = N_H_per_gram_global*rho_A(i, j, k)
              chromo = 0.0_dp
              if (w > 1.0e-6_dp) then
                colm = p_A(i, j, k)/g_c(i)                ! column mass above (hydrostatic)
                call species_emiss_ionfrac(h1_tab, lT, QH, iH)
                call species_emiss_ionfrac(ca2_tab, lT, QCa, iCa)
                call species_emiss_ionfrac(mg2_tab, lT, QMg, iMg)
                chromo = (QH*species_escape_prob(h1_tab, log10(max(4.0e-14_dp*colH(i), 1.0e-30_dp))) &
                    + QCa*species_escape_prob(ca2_tab, log10(max(colm, 1.0e-30_dp))) &
                    + QMg*species_escape_prob(mg2_tab, log10(max(colm, 1.0e-30_dp))))*nHt*ne_A(i, j, k)
              end if
              coron = 0.0_dp
              if (w < 1.0_dp - 1.0e-6_dp) coron = ne_A(i, j, k)*nHt*coronal_lambda(coronal_tab, lT)
              Qlos = Sc*(w*chromo + (1.0_dp - w)*coron)
            end if
            if (use_heating) then
              x = max(z_c(i), 0.0_dp)/(z_heat_km*1.0e5_dp)
              Qheat = Sc*H_peak*x*exp(1.0_dp - x)
              if (H_c > 0.0_dp) Qheat = Qheat + Sc*H_c*exp(-max(z_c(i), 0.0_dp)/(lambda_c_km*1.0e5_dp))
              if (H_ext > 0.0_dp) Qheat = Qheat + Sc*H_ext*exp(-max(z_c(i), 0.0_dp)/(lambda_ext_rsun*Rsun_atm))
              if (H_B > 0.0_dp .and. use_B) Qheat = Qheat + Sc*H_B*norm2(State_GVB(i, j, k, BR:BP, ib)) &
                  *exp(-max(z_c(i), 0.0_dp)/(lambda_B_km*1.0e5_dp))/(1.0_dp + (ne_A(i, j, k)/n_gate)**4)
            end if
          end if
          Qloss_B(i, j, k, ib) = Qlos
          State_GVB(i, j, k, S1, ib) = State_GVB(i, j, k, S1, ib) + dt_in*stf(i)*(Qheat - Qlos)/(rho_A(i, j, k)*T_A(i, j, k))
        end do
      end do
    end do
    end associate

  contains

    ! The (r, theta, phi) gradient of T at interior cell (i, j, k) of the block
    ! (central differences; the ghost layers next to the interior are filled)
    function grad_cell(i, j, k) result(g)
      integer, intent(in) :: i, j, k
      real(dp) :: g(3)
      associate (Y2 => Metric_DB(2, ib)%X_I, Y3 => Metric_DB(3, ib)%X_I)
        g(1) = (T_A(i+1, j, k) - T_A(i-1, j, k))/(r_c(i+1) - r_c(i-1))
        g(2) = (T_A(i, j+1, k) - T_A(i, j-1, k))/(r_c(i)*(Y2(j+1) - Y2(j-1)))
        g(3) = (T_A(i, j, k+1) - T_A(i, j, k-1))/(r_c(i)*sin(Y2(j))*(Y3(k+1) - Y3(k-1)))
      end associate
    end function grad_cell

    ! The gradient that drives the face flux along direction d: w b_d (b.g) +
    ! (1 - w) g_d, b = B/|B|, w = B^2/(B^2 + b_iso_G^2) (cond_along_b)
    real(dp) function along_b(g, bv, d) result(gd)
      real(dp), intent(in) :: g(3), bv(3)
      integer, intent(in) :: d
      real(dp) :: b2, w
      b2 = sum(bv**2)
      w = b2/(b2 + b_iso_G**2)
      gd = (1.0_dp - w)*g(d)
      if (b2 > 0.0_dp) gd = gd + w*bv(d)*dot_product(bv, g)/b2
    end function along_b

    ! CZTR's hyperbolic face flux: target -fsat k0 T^2.5 dT/dx, relaxed;
    ! fS: the Spitzer weight, qcl: the collisionless flux (weight 1 - fS)
    subroutine face_flux(Sf, Tf, rf, csf, cvf, Tg, dx, vf, q, fS, qcl)
      real(dp), intent(in) :: Sf, Tf, rf, csf, cvf, Tg, dx, vf, fS, qcl
      real(dp), intent(inout) :: q
      real(dp) :: k0, Fc, qt, D, budget, tau
      k0 = Sf*(kappa0_spitzer + blend_weight(Tf, T_SWITCH, T_SWITCH_WIDTH)*extra_h_kappa0(hcond_tab, log10(Tf)))
      Fc = k0*Tf**2.5_dp*Tg
      qt = -Fc/(1.0_dp + abs(Fc)/(1.5_dp*rf*csf**3))
      qt = fS*qt + (1.0_dp - fS)*Sf*qcl
      D = fS*k0*Tf**2.5_dp/(rf*cvf)
      budget = max(f_cfl_cond*dx/dt - vf, 1.0e-3_dp*dx/dt)
      tau = max(D/budget**2, 4.0_dp*dt)
      q = qt + (q - qt)*exp(-dt/tau)
    end subroutine face_flux
  end subroutine atm_block

  ! Floors of block ib (interior cells): rho >= rho_floor; T >= T_floor where S > 0.5
  subroutine atm_floors(ib)
    integer, intent(in) :: ib
    real(dp) :: rh, out(n_field)
    integer :: i, j, k
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, nb(1)
          State_GVB(i, j, k, RHO1, ib) = max(State_GVB(i, j, k, RHO1, ib), rho_floor - rho0(i))
          if (atm_switch(z_c(i)) > 0.5_dp) then
            rh = rho0(i) + State_GVB(i, j, k, RHO1, ib)
            call eos_eval(tab, rh, s0(i) + State_GVB(i, j, k, S1, ib), out)
            if (exp(out(F_LNT)) < T_floor) State_GVB(i, j, k, S1, ib) = &
                eos_s_from_rho_T(tab, rh, T_floor, s_guess=s0(i) + State_GVB(i, j, k, S1, ib)) - s0(i)
          end if
        end do
      end do
    end do
  end subroutine atm_floors

  ! dt_rad_factor x the shortest loss time rho cv T/Q_loss of block ib (the
  ! last step's losses; huge if none); rate_r (optional): per radial cell,
  ! raised to cfl over the cell's own loss-time limit (local time stepping)
  real(dp) function atm_dt(ib, rate_r) result(dt)
    integer, intent(in) :: ib
    real(dp), intent(inout), optional :: rate_r(:)
    real(dp) :: out(n_field), rh
    integer :: i, j, k
    dt = huge(1.0_dp)
    if (.not. (use_losses .and. allocated(Qloss_B))) return
    do k = 1, nb(3)
      do j = 1, nb(2)
        do i = 1, nb(1)
          if (atm_switch(z_c(i)) <= 0.5_dp .or. Qloss_B(i, j, k, ib) <= 0.0_dp) cycle
          rh = rho0(i) + State_GVB(i, j, k, RHO1, ib)
          call eos_eval(tab, rh, s0(i) + State_GVB(i, j, k, S1, ib), out)
          dt = min(dt, rh*exp(out(F_LNCV))*exp(out(F_LNT))/Qloss_B(i, j, k, ib))
          if (present(rate_r)) rate_r(i) = max(rate_r(i), cfl*Qloss_B(i, j, k, ib) &
              /(dt_rad_factor*rh*exp(out(F_LNCV))*exp(out(F_LNT))))
        end do
      end do
    end do
    dt = dt_rad_factor*dt
  end function atm_dt

  ! The transition region of block ib's columns: the heights [km] of the
  ! upward crossings of T = 1e5 K with z between 0 and z_max_km (a cell above
  ! 1e5 K whose lower neighbour, possibly the ghost, is below); hmin/hmax are
  ! updated (start them at huge/-huge); the lowest crossing per column only
  subroutine atm_tr_heights(ib, z_max_km, hmin, hmax)
    integer, intent(in) :: ib
    real(dp), intent(in) :: z_max_km
    real(dp), intent(inout) :: hmin, hmax
    real(dp) :: out(n_field), Tlo, Thi
    integer :: i, j, k
    do k = 1, nb(3)
      do j = 1, nb(2)
        call eos_eval(tab, rho0(0) + State_GVB(0, j, k, RHO1, ib), s0(0) + State_GVB(0, j, k, S1, ib), out)
        Tlo = exp(out(F_LNT))
        do i = 1, nb(1)
          call eos_eval(tab, rho0(i) + State_GVB(i, j, k, RHO1, ib), s0(i) + State_GVB(i, j, k, S1, ib), out)
          Thi = exp(out(F_LNT))
          if (z_c(i) > 0.0_dp .and. z_c(i) < z_max_km*1.0e5_dp .and. Thi > 1.0e5_dp .and. Tlo <= 1.0e5_dp) then
            hmin = min(hmin, z_c(i)/1.0e5_dp)
            hmax = max(hmax, z_c(i)/1.0e5_dp)
            exit
          end if
          Tlo = Thi
        end do
      end do
    end do
  end subroutine atm_tr_heights

end module cztr3d_atm_mod
