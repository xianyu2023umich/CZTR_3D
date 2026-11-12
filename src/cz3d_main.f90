program cz3d_run
  ! CZTR step 1 in 3D: the convection zone on a Yin-Yang shell (ARBOR), the
  ! 1D grid in r (whole columns or nb_r cells per block), the 1D state as the
  ! initial condition in every column. See cz3d_app_mod.
  !
  ! Usage: mpirun -np N ./cztr3d.x [run.nml]
  use arbor
  use cz3d_params_mod
  use cz3d_bg_mod, only: bg_load, set_rsst, n_r, r_f, r_c, z_c, k_join, r_top_cz, xi, p0, rho0, T0, g_c, cs0, s0
  use cz3d_transport_mod, only: transport_init, set_top_cooling
  use const_mod, only: Lsun
  use cz3d_app_mod
  use cztr3d_atm_mod, only: atm_init
  implicit none
  type(arbor_domain_t) :: dom
  type(arbor_tree_t) :: tree
  type(cz3d_app_t) :: app
  type(arbor_map1d_t) :: map_r
  character(len=256) :: file
  character(len=300) :: msg
  real(dp) :: t, wall
  integer :: n_step, i_step, n_var_file, n_chunk, n_done
  real(dp) :: wall_start, elapsed
  logical :: stopped_wall
  character(len=300) :: last_rst
  logical :: have_rst, regridded
  character(len=32) :: env
  integer :: n_env, ios
  real(dp) :: wall_job, t_relevel

  call arbor_init()
  wall_start = arbor_wtime()
  file = 'cz3d.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, file)
  call read_params(file)
  ! batch jobs: the job's own remaining wall time [s] (CZ3D_WALL_TIME_MAX, set
  ! by tools/job.* from squeue at the launch) caps wall_time_max, so a run
  ! submitted with a shorter -t than the namelist assumed still stops with a
  ! restart (the smallest over the ranks, in case a rank lacks the variable)
  call get_environment_variable('CZ3D_WALL_TIME_MAX', env, length=n_env)
  wall_job = huge(1.0_dp)
  if (n_env > 0) then
    read(env, *, iostat=ios) wall_job
    if (ios /= 0 .or. wall_job <= 0.0_dp) wall_job = huge(1.0_dp)
  end if
  wall_job = arbor_min(wall_job)
  if (wall_job < huge(1.0_dp) .and. (wall_time_max <= 0.0_dp .or. wall_job < wall_time_max)) then
    wall_time_max = wall_job
    write(msg, '(A,F10.1,A,F8.1,A)') 'cz3d: the job''s wall time ', wall_time_max, ' s sets wall_time_max (margin ', &
        wall_time_margin, ' s)'
    call arbor_print(msg)
    if (wall_time_margin >= wall_time_max) call arbor_abort('wall_time_margin must be below the job''s wall time')
  end if
  call bg_load(grid_file, init_file, table_file, n_r_cells, r_top_rsun, g_from_background)
  if (use_atmosphere) call atm_init()
  call transport_init(opacity_file)
  if (cool_width_hp > 0.0_dp) cool_width_km = cool_width_hp*p0(n_r)/(rho0(n_r)*g_c(n_r))/1.0e5_dp
  if (top_energy == 'cooling') call set_top_cooling(L_bottom_Lsun*Lsun, cool_width_km*1.0e5_dp)
  call set_rsst(c_target, rsst_q, rsst_depth_km)

  map_r = arbor_map_table(r_f(0:n_r))
  call arbor_map_set_extension(map_r, MAP_EXT_LINEAR, MAP_EXT_LINEAR)
  if (nb_r <= 0) nb_r = n_r
  if (mod(n_r, nb_r) /= 0) call arbor_abort('cz3d: nb_r must divide the radial cells')
  if (yy_margin >= 0) then
    call arbor_domain_yinyang(dom, [nb_r, nb_theta, nb_phi], 2, map_r, n_theta, n_phi, margin_cells=yy_margin)
  else
    call arbor_domain_yinyang(dom, [nb_r, nb_theta, nb_phi], 2, map_r, n_theta, n_phi)
  end if
  if (use_B) then
    call arbor_set_vectors([vr_, br_])            ! velocity and field: rotated between the Yin-Yang frames
  else
    call arbor_set_vectors([vr_])
  end if
  call arbor_block_init(dom, n_var=cz_n_var(), n_var_face=merge(1, 0, qconv_face))
  if (len_trim(restart_from) > 0) then
    ! tree, blocks, time, step and the profile counter from the file
    ! with use_B the file may be hydrodynamic (fewer variables): B and psi
    ! then start at zero and get the seed below
    if (use_B) then
      call arbor_restart_read(app, restart_from, dom, tree, t, i_step, n_var_file=n_var_file, &
          domain_tol=restart_domain_tol)
    else
      call arbor_restart_read(app, restart_from, dom, tree, t, i_step, domain_tol=restart_domain_tol)
      n_var_file = cz_n_var()
    end if
    call arbor_ghost_build(tree, dom, bc=wall_bc)
    call cz_setup()
    if (n_var_file < cz_n_var()) call cz_field_from_hydro()
    if (tloc(1) < 0.0_dp) tloc = t                        ! (no shell clocks in the file: all at the run's time)
    restarted = .not. restart_new_output
    if (restart_new_output) then
      i_prof = 0
      t_next_prof = t
      t_prof0 = t
    end if
  else if (len_trim(import_dir) > 0) then
    ! a run continued from export files (any radial block size dividing the
    ! exported one): the tree refined to the exported blocks' levels, each
    ! block's state, time, step and shell clocks from the files
    call import_open()
    if (import_max_level() > max(refine_h, patch_level)) &
        call arbor_abort('cz3d: import: the export is refined beyond max(refine_h, patch_level)')
    call arbor_tree_init(tree, dom, [0, max(refine_h, patch_level), max(refine_h, patch_level)])
    call arbor_tree_refine(tree, dom, import_criterion)
    call arbor_block_sync(tree)
    call arbor_ghost_build(tree, dom, bc=wall_bc)
    call cz_setup()
    call arbor_init_blocks(app)
    t = 0.0_dp
    i_step = 0
    call import_meta(t, i_step, have_rst)
    write(msg, '(A,A,A,I0,A,I0,A,L1)') '  imported from ', trim(import_dir), ': ', tree%n_leaf, ' blocks, step ', i_step, &
        ', time and clocks restored: ', have_rst
    call arbor_print(msg)
  else
    call arbor_tree_init(tree, dom, [0, 0, 0])
    call arbor_block_sync(tree)
    call arbor_ghost_build(tree, dom, bc=wall_bc)
    call cz_setup()
    call arbor_init_blocks(app)
    t = 0.0_dp
    i_step = 0
  end if
  ! the PFSS field above the join from the CZ's B_r: on a fresh start, an
  ! import, or a hydrodynamic restart (not on a restart that has B)
  if (use_B .and. (b_seed == 'pfss' .or. b_seed == 'dpf') .and. (len_trim(restart_from) == 0 .or. n_var_file < cz_n_var())) &
      call cz_seed_pfss()
  if (is_root) call execute_command_line('mkdir -p '//trim(out_dir))
  call arbor_barrier()
  if (is_root) call write_background()

  write(msg, '(A,I0,A,I0,A,I0,A,F6.3,A,I0,A,I0,A)') 'cz3d: ', n_r, ' radial cells (1D join cell ', k_join, '), Yin-Yang ', &
      n_theta, ' x ', 90.0_dp/real(n_theta - 2*merge(yy_margin, 3, yy_margin >= 0), dp), ' deg, ', tree%n_leaf, ' blocks on ', n_ranks, ' ranks'
  call arbor_print(msg)
  if (use_atmosphere) then
    write(msg, '(A,L1,A,L1,A,L1,A,F8.1,A,F6.1,A,ES9.2,A,F8.1,A,L1)') '  atmosphere: conduction ', use_conduction, &
        ', losses ', use_losses, ', heating ', use_heating, '; switch at z = ', z_atm_km, ' km, width ', w_atm_km, &
        ' km; H_peak ', H_peak, ', z_heat ', z_heat_km, ' km; CZ frozen ', freeze_cz
    call arbor_print(msg)
  end if
  write(msg, '(A,F9.1,A,F9.1,A,A,A,A,A,L1,A,F5.2)') '  r from ', (r_f(0) - r_top_cz)/1.0e5_dp, ' to ', &
      (r_f(n_r) - r_top_cz)/1.0e5_dp, ' km; walls: bottom ', trim(wall_bottom), ', top ', trim(wall_top), &
      '; MLT switch ', mlt_switch, '; diffusion h = ', diffusion_h
  call arbor_print(msg)
  write(msg, '(A,L1,A,L1)') '  diffusion heating ', diffusion_heating, '; well-balanced gravity ', g_from_background
  call arbor_print(msg)
  write(msg, '(A,L1,A,L1,A,A,A,F8.1,A)') '  convection (MLT) ', use_convection, ', radiation ', use_radiation, &
      ', top energy condition: ', trim(top_energy), ' (cooling width ', cool_width_km, ' km)'
  call arbor_print(msg)
  if (rsst_depth_km > 0.0_dp) then
    write(msg, '(A,F8.1,A,F6.2,A,F8.2,A,F8.2)') '  reduced sound speed anchored at ', rsst_depth_km, ' km depth, q = ', &
        rsst_q, '; xi at the bottom ', xi(1), ', at the top ', xi(n_r)
    call arbor_print(msg)
  else if (c_target > 0.0_dp) then
    write(msg, '(A,F7.3,A,F8.2,A,F8.2)') '  reduced sound speed: c_target = ', c_target/1.0e5_dp, ' km/s; xi at the bottom ', &
        xi(1), ', at the top ', xi(n_r)
    call arbor_print(msg)
  end if

  regridded = .false.
  if (max(refine_h, patch_level) > 0) call refine_horizontally()
  if (mr_levels > 1 .and. (.not. UseMultirate .or. mr_relevel)) then
    ! rate levels from the blocks' stable dt (a multirate restart file brings
    ! its levels, unless mr_relevel); blocks may move between ranks: the
    ! per-block arrays again (the convective fluxes restart at their targets below)
    call arbor_multirate_setup(app, tree, t, ratio=2, n_level_max=mr_levels, extrapolate_faster=mr_extrapolate, &
        uniform_dir=merge(1, 0, mr_by_shell))
    call cz_setup()
  end if

  if (plot_blocks) then
    ! (each block drawn as 1 x 4 x 4 bricks so the curved blocks bend)
    call arbor_plot_blocks(tree, dom, trim(out_dir)//'/blocks.dat', sub=[1, 4, 4])
    call arbor_print('  block outlines: '//trim(out_dir)//'/blocks.dat')
  end if
  ! as CZTR: the convective fluxes start at their targets; after a restart
  ! they continue from the file's face state (qconv_face; not if the blocks changed)
  call cz_transport_all(iBlockList_I(1:nBlock), 0.0_dp, &
      .not. (qconv_face .and. len_trim(restart_from) > 0 .and. .not. regridded))
  call budget_start(t)
  call cz_log(i_step, t)
  if ((dt_prof > 0.0_dp .or. prof_every > 0) .and. .not. restarted .and. t >= t_next_prof) call cz_profile(i_step, t)

  wall = arbor_wtime()
  stopped_wall = .false.
  t_relevel = huge(1.0_dp)
  if (mr_relevel_dt > 0.0_dp .and. UseMultirate) t_relevel = t + mr_relevel_dt
  if (wall_time_max <= 0.0_dp .and. t_relevel == huge(1.0_dp)) then
    call arbor_run(app, t, t_end, dt_max=dt_max, n_step_max=n_step_max, n_step=n_step, tree=tree, dom=dom, &
        step0=i_step, plot_dt=plot_dt, restart_dt=restart_dt, restart_file=trim(out_dir)//'/restart')
  else
    ! in chunks of wall_check_steps steps, the wall time checked between them
    ! (the restart and plot schedules depend only on t and the global step),
    ! and the levels set again once t passes t_relevel (mr_relevel_dt)
    n_step = 0
    do
      n_chunk = min(wall_check_steps, n_step_max - n_step)
      if (n_chunk <= 0 .or. t >= t_end*(1.0_dp - 1.0e-12_dp)) exit
      call arbor_run(app, t, t_end, dt_max=dt_max, n_step_max=n_chunk, n_step=n_done, tree=tree, dom=dom, &
          step0=i_step + n_step, plot_dt=plot_dt, restart_dt=restart_dt, restart_file=trim(out_dir)//'/restart')
      n_step = n_step + n_done
      if (n_done == 0) exit
      if (t >= t_relevel .and. t < t_end*(1.0_dp - 1.0e-12_dp)) then
        call arbor_multirate_setup(app, tree, t, ratio=2, n_level_max=mr_levels, extrapolate_faster=mr_extrapolate, &
        uniform_dir=merge(1, 0, mr_by_shell))
        call cz_setup()
        t_relevel = t + mr_relevel_dt
      end if
      if (wall_time_max <= 0.0_dp) cycle
      elapsed = arbor_max(arbor_wtime() - wall_start)
      if (elapsed + wall_time_margin >= wall_time_max .and. t < t_end*(1.0_dp - 1.0e-12_dp) &
          .and. n_step < n_step_max) then
        stopped_wall = .true.
        write(msg, '(A,F10.1,A,F10.1,A,F8.1,A)') 'cz3d: wall time ', elapsed, ' s + margin ', wall_time_margin, &
            ' s reaches the limit ', wall_time_max, ' s: stopping with a restart file'
        call arbor_print(msg)
        exit
      end if
    end do
  end if
  wall = arbor_wtime() - wall
  ! the largest Mach number's cell (stderr), and a plot at the end if asked
  call max_mach_report(i_step + n_step, t)
  if (plot_at_end) call app%write_plot(dom, i_step + n_step, t)
  ! a restart file at the end too (unless the schedule just wrote one), so a
  ! run can always be continued from where it stopped; always when stopped
  ! for the wall time
  last_rst = arbor_restart_name(trim(out_dir)//'/restart', i_step + n_step)
  if (restart_dt > 0.0_dp .or. stopped_wall) then
    inquire(file=trim(last_rst), exist=have_rst)
    if (.not. have_rst) call arbor_restart_write(app, tree, dom, t, i_step + n_step, trim(last_rst))
  end if
  ! status for job chains: 'stopped_wall <restart>' or 'finished <restart or none>'
  if (is_root) then
    open(newunit=n_chunk, file=trim(out_dir)//'/run_status.txt', status='replace', action='write')
    if (stopped_wall) then
      write(n_chunk, '(A,1X,A)') 'stopped_wall', trim(last_rst)
    else if (restart_dt > 0.0_dp) then
      write(n_chunk, '(A,1X,A)') 'finished', trim(last_rst)
    else
      write(n_chunk, '(A)') 'finished none'
    end if
    write(n_chunk, '(A,I0,A,ES16.8)') 'step ', i_step + n_step, '  t ', t
    close(n_chunk)
  end if
  write(msg, '(A,I0,A,ES11.4,A,F10.1,A,F8.4,A)') 'done: ', n_step, ' steps, t = ', t, ' s, wall ', wall, ' s (', &
      1.0e6_dp*wall*n_ranks/(real(max(n_step, 1), dp)*real(tree%n_leaf, dp)*product(nb)), ' us per cell-step per rank)'
  call arbor_print(msg)
  call dump_shells_write()
  call export_state_write(t, i_step + n_step)
  write(msg, '(A,I0)') 'state checksum: ', state_checksum()
  call arbor_print(msg)
  call arbor_finalize()

contains

  ! The background and the RSST factor per radial cell, out_dir/background.dat
  ! (for scripts/estimate_steps.py)
  subroutine write_background()
    integer :: u, i
    open(newunit=u, file=trim(out_dir)//'/background.dat', status='replace', action='write')
    write(u, '(A,F10.5,A,ES12.4,A,I0,A,I0)') '# r_top_rsun ', r_top_rsun, '  c_target ', c_target, '  n_theta ', n_theta, &
        '  yy_margin ', yy_margin
    write(u, '(A)') '# i  r[cm]  z[km]  dr[cm]  rho0  p0  T0  g  cs0[cm/s]  xi  s0[erg/g/K]'
    do i = 1, n_r
      write(u, '(I5,ES20.12,9ES16.7E3)') i, r_c(i), z_c(i)/1.0e5_dp, r_f(i) - r_f(i-1), rho0(i), p0(i), T0(i), g_c(i), &
          cs0(i), xi(i), s0(i)
    end do
    close(u)
  end subroutine write_background

  ! The same blocks in a tree that allows max(refine_h, patch_level) levels
  ! in theta and phi, then regrids until every block has its level
  ! (cz_flag_block: refine_h everywhere, patch_level in the patch)
  subroutine refine_horizontally()
    type(arbor_tree_t) :: t2
    logical :: ok
    character(len=200) :: m
    integer :: lev, lmax, n_changed, n0
    real(dp) :: d0
    lmax = max(refine_h, patch_level)
    n0 = tree%n_leaf
    call arbor_tree_from_keys(t2, dom, [0, lmax, lmax], tree%key(:, 1:tree%n_leaf), ok, m)
    if (.not. ok) call arbor_abort('cz3d: refine: '//trim(m))
    if (t2%n_leaf /= tree%n_leaf) call arbor_abort('cz3d: refine: the tree changed')
    if (any(t2%key(:, 1:t2%n_leaf) /= tree%key(:, 1:tree%n_leaf))) &
        call arbor_abort('cz3d: refine: the leaf order changed with l_max')
    tree = t2
    do lev = 1, lmax
      call arbor_regrid(app, tree, dom, t, n_changed=n_changed)
      if (n_changed > 0) regridded = .true.
    end do
    call cz_setup()
    d0 = 90.0_dp/real(n_theta - 2*merge(yy_margin, 3, yy_margin >= 0), dp)
    write(msg, '(A,F7.4,A,I0,A,I0,A)') '  refined horizontally to ', d0/2.0_dp**refine_h, ' deg: ', &
        tree%n_leaf, ' blocks (', n0, ' before)'
    call arbor_print(msg)
    if (patch_level > refine_h) then
      write(msg, '(A,F7.4,A,F7.2,A,F8.2,A,F6.2,A,I0,A)') '  patch at ', d0/2.0_dp**patch_level, ' deg around lat ', &
          patch_lat_deg, ', lon ', patch_lon_deg, ' (radius ', patch_radius_deg, ' deg): ', &
          count_level(patch_level), ' blocks'
      call arbor_print(msg)
    end if
  end subroutine refine_horizontally

  integer function count_level(lev) result(n)
    integer, intent(in) :: lev
    n = count(tree%key(KL+2, 1:tree%n_leaf) == lev)
  end function count_level

end program cz3d_run
