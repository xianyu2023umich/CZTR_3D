module cz3d_params_mod
  ! Run-time parameters of the 3D convection-zone run (namelist &cz3d in the
  ! parameter file, read by the root and shared through ARBOR's config).
  use arbor
  implicit none
  public

  ! --- input / output ---
  ! inputs in data/ (`make data` generates the EOS table and the static
  ! background there; the grid and the opacity table are kept in the repository)
  character(len=256) :: grid_file    = 'data/grid_full.dat'   ! merged 1D grid (faces of cells 1..n_r)
  ! 1D state + background. Default: the static CZ envelope (hydrostatic, from
  ! tools/static_bg.x + tools/build_static_init.py; state = background, so
  ! the run starts from ds1_random). Until 2026-10-14 the runs used CZTR's
  ! step-1 state '../CZTR/output/step1_cz_alone_dtcz/final_init.dat' (the same
  ! background below ~1.8 Mm depth; above it, CZ_1D's diffusion-balanced
  ! profile, up to 4% off hydrostatic at ~300 km)
  character(len=256) :: init_file    = 'data/static_init.dat'
  character(len=256) :: table_file   = 'data/eos_table.bin'                  ! eos_hhe runtime table (make data)
  character(len=256) :: opacity_file = 'data/opacity_merged_x70_z02.dat'
  character(len=256) :: out_dir      = 'output/run'

  ! --- grid: n_r radial cells (cells 1..n_r of the 1D grid; 0 = the join
  !     cell k of init_file), Yin-Yang n_theta x n_phi per patch, blocks of
  !     nb_r x nb_theta x nb_phi; nb_r = 0: the whole radial column in one
  !     block (nb_r = n_r), else nb_r must divide n_r ---
  integer :: n_r_cells = 0
  real(dp) :: r_top_rsun = 0.0_dp       ! > 0: the top cell is the last with its centre below r_top_rsun x R_sun (overrides n_r_cells)
  integer :: n_theta = 54, n_phi = 150
  integer :: nb_r = 0, nb_theta = 9, nb_phi = 25
  ! multirate (ARBOR's rate levels, ratio 2, from each block's stable dt;
  ! time accurate, global steps): mr_levels > 1 levels at most (needs
  ! radial blocks to gain); mr_extrapolate: ghosts from faster neighbours
  ! extrapolated in time (2nd order) instead of frozen at the substep start
  integer :: mr_levels = 0
  logical :: mr_extrapolate = .false.
  ! re-levelling: the levels from the blocks' current stable dt again
  ! (otherwise fixed at the setup, and a multirate restart keeps its own):
  ! mr_relevel at a restart, and every mr_relevel_dt [s] of model time during
  ! the run (at the end of a level-1 step; 0 = never). Blocks may move between
  ! ranks (each level is balanced on its own)
  logical :: mr_relevel = .false.
  real(dp) :: mr_relevel_dt = 0.0_dp
  ! one level per radial shell of blocks (the shell's fastest block): the
  ! near-surface layer on a single step (per-block levels put neighbouring top
  ! columns on different steps; the near-surface state depends on dt: mhd4)
  logical :: mr_by_shell = .true.
  ! horizontal refinement: after the start or restart, every block is
  ! refined refine_h levels in theta and phi (ARBOR's conservative
  ! interpolation), 2^refine_h times finer angles
  integer :: refine_h = 0
  ! a refined patch: blocks with a cell centre within patch_radius_deg of
  ! (patch_lat_deg, patch_lon_deg) (the global frame, the Yin patch's centre at
  ! 0, 0) refined to patch_level (> refine_h); the neighbours follow (2:1)
  integer :: patch_level = 0
  real(dp) :: patch_lat_deg = 0.0_dp, patch_lon_deg = 0.0_dp, patch_radius_deg = 5.0_dp
  ! debugging probe: every step, the top probe_cells radial cells of every
  ! column within probe_radius_deg (> 0) of (probe_lat_deg, probe_lon_deg) to
  ! out_dir/probe_rNNNN.dat
  real(dp) :: probe_lat_deg = 0.0_dp, probe_lon_deg = 0.0_dp, probe_radius_deg = 0.0_dp
  integer :: probe_cells = 4
  integer :: yy_margin = -1               ! Yin-Yang margin cells per side (< 0: ARBOR's default, ng + 1)

  ! --- time stepping ---
  real(dp) :: cfl = 0.5_dp
  real(dp) :: t_end = 3600.0_dp
  real(dp) :: dt_max = 20.0_dp
  integer  :: n_step_max = huge(1)
  integer  :: n_log = 100                 ! log line every n_log steps
  real(dp) :: dt_prof = 600.0_dp          ! horizontally averaged profile every dt_prof s (0 = none)
  real(dp) :: plot_dt = 0.0_dp            ! Tecplot file of the state every plot_dt s (0 = none)
  ! restart files out_dir/restart_NNNNNNNN.rst every restart_dt s (0 = none);
  ! restart_from: continue from that file (log.txt is appended to). The
  ! mixing-length flux memory is not stored: it restarts at its targets.
  real(dp) :: restart_dt = 0.0_dp
  ! energy budget per radial shell (budget_shells.dat, budget_log.txt in
  ! out_dir) every dt_budget s (0 = off); see cz3d_app_mod
  real(dp) :: dt_budget = 0.0_dp
  character(len=256) :: restart_from = ''
  ! local time stepping for relaxation (not time accurate; 0 = off): each
  ! radial shell advances min(its own stable dt, dt1) per step (at least the
  ! global dt), so the deep layers run ahead of the top. t is then the top's
  ! clock only; the log and the profiles carry each shell's own clock.
  ! Read at every start, so it can be changed at a restart.
  real(dp) :: dt1 = 0.0_dp
  ! > 0: dt1 falls linearly to 0 over this many steps from the run's start
  ! (or restart), turning the local stepping off smoothly (alpha -> 1 from the
  ! top down as dt1 drops below each shell's own step)
  integer  :: dt1_ramp_steps = 0
  integer  :: prof_every = 0              ! > 0: a profile every prof_every steps (instead of dt_prof)
  ! > 0: accept a restart file written on another machine (e.g. 1e-10), whose
  ! grid (map ends, face positions) differs in the last bits; this run uses its
  ! own grid (bit differences: fine for plots and analysis; a continued run is
  ! not bitwise the same as one on the original machine). 0: exact match
  real(dp) :: restart_domain_tol = 0.0_dp
  logical :: restart_new_output = .false.  ! after a restart: a new profile series (from 0, the first now)
  ! batch jobs: stop cleanly, with a restart file, once the elapsed wall time
  ! (from the program's start) plus wall_time_margin reaches wall_time_max
  ! [s] (0 = no limit); checked every wall_check_steps steps. The margin must
  ! cover one check interval plus writing the restart. out_dir/run_status.txt
  ! then says 'stopped_wall' (continue from the restart it names) or 'finished'.
  real(dp) :: wall_time_max = 0.0_dp
  real(dp) :: wall_time_margin = 1800.0_dp
  integer  :: wall_check_steps = 20

  ! --- physics ---
  logical  :: use_radiation  = .true.     ! grey Eddington radiation, radial (per column)
  logical  :: use_convection = .true.     ! mixing-length flux, radial (per column)
  real(dp) :: grey_top_km = 0.0_dp        ! implicit grey radiation only in the top grey_top_km (0: the whole column); explicit diffusion below
  logical  :: use_hyperbolic = .true.     ! hyperbolic relaxation of the convective flux
  logical  :: qconv_face = .true.         ! the relaxing flux in ARBOR's face state (in restart files, moves with
                                          ! blocks); .false.: a per-rank array, reset at restarts (older restart files)
  real(dp) :: alpha_mlt = 1.314263_dp
  real(dp) :: L_bottom_Lsun = 1.0_dp      ! luminosity through the bottom face
  real(dp) :: grey_tau_above = 0.01_dp    ! optical depth above the top cell (Eddington top condition)
  ! top energy condition: 'eddington' (radiation escapes through the top
  ! face, CZTR step 1) or 'cooling' (no radiative flux through the top face;
  ! a cooling layer removes L_bottom, its flux F(r) = F_top exp(-(r_top - r)/w)
  ! with w = cool_width_km, as convection.F90; for a top deep in the CZ)
  character(len=16) :: top_energy = 'eddington'
  real(dp) :: cool_width_km = 1000.0_dp
  real(dp) :: cool_width_hp = 0.0_dp      ! > 0: width = cool_width_hp x H_p of the top cell (overrides cool_width_km)
  ! MLT as a subgrid model: f = 1/(1 + (R/mlt_R0)^mlt_p), R = H_p/(r dtheta)
  ! the number of pressure scale heights per horizontal cell; off: f = 1
  logical  :: mlt_switch = .false.
  real(dp) :: mlt_R0 = 4.0_dp, mlt_p = 4.0_dp

  ! --- well balanced gravity: g = -(dp0/dr)/rho0 with the momentum
  !     equation's radial stencil in cells 1..n_r, so the background force
  !     (dp0/dr + rho0 g), which the momentum equation leaves out, is exactly
  !     0 (off: g of the state file, as CZTR) ---
  logical :: g_from_background = .false.

  ! --- reduced sound speed (RSST, after Hotta et al.): at most c_target
  !     [cm/s] everywhere, xi = max(1, c_s0/c_target), the continuity
  !     equation divided by xi^2 (0 = off; as CZTR's c_target) ---
  real(dp) :: c_target = 0.0_dp
  ! depth profile of the reduced sound speed: c_red = min(c_s0, c_target
  ! (c_s0/c_s0,top)^rsst_q), xi = c_s0/c_red; 0: c_red = c_target everywhere
  ! (Hotta et al. 2014), 1: xi uniform (Hotta et al. 2012)
  real(dp) :: rsst_q = 0.0_dp
  ! the form of the RSST: 'hotta2012' (continuity divided by xi^2; energy not
  ! conserved, Hotta et al. 2014 App. B) or 'hotta2015' (Hotta et al. 2015
  ! App. A: rho1 is the scaled perturbation xi^2 rho1 with the ordinary
  ! continuity equation, the EOS and the buoyancy see rho1/xi^2; conserves
  ! rho T s1 + rho v^2/2 to linear order)
  character(len=16) :: rsst_form = 'hotta2012'
  ! > 0: anchor the reduced sound speed at this depth below the photosphere
  ! instead of c_target: c_ref = c_s0 there, c_red = min(c_s0, c_ref
  ! (c_s0/c_ref)^rsst_q), so no reduction above that depth (xi = 1) and xi =
  ! (c_s0/c_ref)^(1 - rsst_q) below (c_target is then not used)
  real(dp) :: rsst_depth_km = 0.0_dp
  ! diagnostic: radial cell indices whose every cell is written at the end of
  ! the run (out_dir/shell_dump_rNNNN.dat per rank; 0 = unused)
  integer :: dump_shells(8) = 0

  ! --- magnetic field (ideal MHD, after solar_model): B = (B_r, B_theta,
  !     B_phi) [G] and the GLM cleaning scalar psi as state variables 6-9;
  !     dB/dt = curl(v x B) - grad psi, the Lorentz force (curl B) x B/(4 pi
  !     rho) in the momentum equation, dpsi/dt = -c_h^2 div B with c_h = c_red
  !     + v_A, then psi damped by exp(-dt alpha_glm c_h/dx_min); the Alfven
  !     speed enters the time step and the diffusion's speed. B and psi are
  !     diffused like the velocity (without the density weight); with
  !     diffusion_heating the dissipated magnetic energy returns as heat ---
  ! --- rotation: the frame rotates with angular velocity omega_rot [rad/s]
  !     about the global z axis (Yin's z); the Coriolis acceleration -2
  !     Omega x v in the momentum equation (no centrifugal term: Omega^2 r/g
  !     ~ 2e-5). Solar: ~2.6e-6 (414 nHz); 0 = off ---
  real(dp) :: omega_rot = 0.0_dp
  ! how: 'split' (exact rotation of v about the axis by -2 Omega dt/2 before
  ! and after each step, Strang; no stability limit, |v| kept exactly) or
  ! 'rk' (the term in the RK stages; RK4-like: unstable for 2 Omega dt > ~2.8,
  ! damping/phase errors ~(2 Omega dt)^5-6 per step)
  character(len=8) :: coriolis = 'split'
  ! global angular-momentum correction: the scheme does not conserve the
  ! flow's L_z about the global axis; every lz_correct_steps steps and/or
  ! lz_correct_dt seconds (0 = off) the part of L_z's change not due to the
  ! Coriolis torque is removed by a rigid rotation delta_Omega x r (L_z then
  ! follows its value at the run's start plus the Coriolis impulse)
  integer  :: lz_correct_steps = 0
  real(dp) :: lz_correct_dt = 0.0_dp

  logical  :: use_B = .false.
  real(dp) :: alpha_glm = 0.5_dp
  ! the GLM cleaning speed c_h: 'fast' (c_s0/xi + v_A, as solar_model) or
  ! 'cfl' (the largest the cell's step allows: glm_cfl_frac cfl / [(1/dr +
  ! 1/(r dtheta) + 1/(r sin theta dphi)) dt_cell], at least c_s0/xi + v_A;
  ! much faster than 'fast' wherever the cell's step is below its own limit,
  ! e.g. deep cells under global stepping)
  character(len=8) :: glm_ch = 'fast'
  ! factor on the diffusion of B and psi only (1: as the velocity; 0: off,
  ! for tests)
  real(dp) :: b_diff_factor = 1.0_dp
  real(dp) :: glm_cfl_frac = 0.8_dp
  ! walls: 'conductor' (B_r = 0, B_h mirrored) or 'radial' (B_h = 0, B_r
  ! mirrored: the field vertical at the wall); Hotta et al. 2014: conductor at
  ! the bottom, vertical at the top
  character(len=12) :: b_bc_bottom = 'conductor', b_bc_top = 'radial'
  ! Alfven-speed cap (after Rempel 2009/2017) [cm/s]: the Lorentz force times
  ! f = 1/sqrt(1 + (v_A/va_max)^4), so the waves' speed sqrt(f) v_A stays
  ! below ~va_max; that speed (not v_A) enters the time step, the diffusion's
  ! speed and c_h. Not energy conserving where it acts (W_Lorentz in the log
  ! is the scaled work: its sum with dE_mag/dt is the cap's energy sink). 0 = off
  real(dp) :: va_max = 0.0_dp
  ! initial field added at the start (also to an imported state):
  ! 'none' or 'toroidal': B_phi = b_seed_G sin(theta) w(r) about the global z
  ! axis (divergence free), w = 1 below b_seed_r1_rsun, 0 above
  ! b_seed_r2_rsun, a sin^2 step between (keeps it off a vertical-field top)
  character(len=12) :: b_seed = 'none'
  real(dp) :: b_seed_G = 100.0_dp, b_seed_r1_rsun = 0.95_dp, b_seed_r2_rsun = 0.98_dp

  ! --- state transfer between runs with different variable sets (ARBOR's
  !     restart files hold a fixed number of variables): export_state writes
  !     every block's interior state at the end of the run
  !     (out_dir/state_export_rNNNN.bin, and state_export_meta.bin: time,
  !     step, shell clocks); import_dir continues from such files: the same
  !     radial grid and horizontal blocks, any radial block size (nb_r) that
  !     divides the exported one (whole columns -> radial blocks), the tree
  !     refined to the exported levels (refine_h or patch_level at least as
  !     high), the common variables copied, the rest (B, psi) zero before the
  !     seed; the convective flux restarts at its targets ---
  logical :: export_state = .false.
  ! out_dir/blocks.dat at the start (after any refinement): the block outlines
  ! as a Tecplot ASCII file (ARBOR's arbor_plot_blocks: leaf, patch, levels, rank)
  logical :: plot_blocks = .false.
  character(len=256) :: import_dir = ''

  ! --- numerical diffusion (Rempel form, as CZTR) ---
  logical  :: use_diffusion = .true.
  logical  :: use_minmod    = .true.
  real(dp) :: diffusion_h   = 1.0_dp
  ! top layer of full (unlimited, first-order) diffusion: depth < diff_top_km
  ! full, a cos^2 ramp to the normal limited scheme at 2 diff_top_km (0: off).
  ! Keeps the grid from overturning at its own scale where granulation is
  ! unresolved (the MLT carries the flux there).
  real(dp) :: diff_top_km   = 0.0_dp
  ! energy-conserving diffusion: the energy the diffusion's mass, momentum
  ! and entropy fluxes carry gets its own face flux, and each cell's entropy
  ! is corrected so its energy changes by exactly -dt div F_E (the kinetic
  ! energy it dissipates returns as heat; see cz3d_app_mod). Off: as CZTR.
  logical  :: diffusion_heating = .false.

  ! --- walls: 'face' (wall on the boundary face, ghosts mirror cells 1, 2)
  !     or 'cell' (wall at the boundary cell's centre, v_r = 0 there, ghosts
  !     mirror the two cells inside); CZTR step 1: bottom face, top cell ---
  character(len=8) :: wall_bottom = 'face', wall_top = 'cell'
  ! wall_top = 'outflow' (solar wind, the top far out, ~20 R_sun, where the
  ! flow is super-Alfvenic): the ghosts above the top cell extrapolated: the
  ! total density ~ r^-2, v_r copied but not inward (>= 0), v_theta, v_phi
  ! and the total entropy copied; B_r ~ r^-2, B_theta, B_phi, psi copied;
  ! the top face's fluxes from these ghosts (not from the wall symmetry)

  ! --- initial perturbation: random s1 of this amplitude [erg/g/K] (0 = none) ---
  real(dp) :: ds1_random = 0.0_dp

  ! --- atmosphere (CZTR's TR physics, cztr3d_atm_mod; ported from CZTR_2D),
  !     use_atmosphere: on above the CZ top through the switch
  !     S(z) = (1 + tanh((z - z_atm)/w_atm))/2: hyperbolic Spitzer conduction
  !     at the faces in r, theta and phi (HYDRAD's H I term, Rempel
  !     saturation, relaxed in tau >= 4 dt, signal speed <= f_cfl_cond dx/dt;
  !     each block keeps its own copy of a face shared with a neighbour),
  !     Carlsson & Leenaarts + coronal losses (the column mass from p/g, the
  !     neutral-H column integrated down within the block plus an estimate
  !     above it), the heating H_peak x e^(1-x), x = z/z_heat, a temperature
  !     floor where S > 0.5 and a density floor everywhere; dt <= dt_rad_factor
  !     x the shortest loss time. The background force (dp0/dr + rho0 g) acts
  !     above the join (the atmosphere's background is not hydrostatic); below
  !     it is left out (well balanced). Grey radiation and MLT stay in the CZ
  !     (cells <= k_join, the Eddington top at the photosphere as without the
  !     atmosphere). Works with radial blocks and multirate. ---
  logical  :: use_atmosphere = .false.
  ! more heating, S [H_c exp(-z/lambda_c) + H_ext exp(-z/lambda_ext)]
  ! [erg cm^-3 s^-1] (Downs et al. 2010: the quiet-Sun term H_c 7.28e-5 with
  ! lambda_c 40 Mm, the extended term H_ext 5e-7 with lambda_ext 0.7 R_sun,
  ! which keeps the corona hot far out and drives the solar wind); 0 = off
  real(dp) :: H_c = 0.0_dp, lambda_c_km = 4.0e4_dp, H_ext = 0.0_dp, lambda_ext_rsun = 0.7_dp
  ! conduction blended to a collisionless flux far out (AWSoM, van der Holst
  ! et al. 2014): q = f_S q_Spitzer + (1 - f_S) (3/2) alpha_cl p_e v_r (radial
  ! faces), f_S = 1/(1 + (r/cond_r_rsun)^2); cond_r_rsun = 0: Spitzer only
  real(dp) :: cond_r_rsun = 0.0_dp, alpha_cl = 1.05_dp
  logical  :: use_conduction = .true., use_losses = .true., use_heating = .true.
  real(dp) :: z_atm_km = 0.0_dp, w_atm_km = 3.0_dp
  real(dp) :: H_peak = 1.0e-3_dp, z_heat_km = 1000.0_dp
  real(dp) :: T_floor = 4400.0_dp, rho_floor = 1.0e-20_dp
  real(dp) :: f_cfl_cond = 0.8_dp, dt_rad_factor = 0.1_dp
  character(len=256) :: atomic_dir   = 'data/atomic'                        ! (copied from FALCON/TR/atomic_data)
  character(len=256) :: coronal_file = 'data/atomic/RadCoolCorona_8.0.dat'   ! (from FALCON/TR)
  ! the CZ cells 1..k_join held at their starting state after every step (CZTR step 2)
  logical  :: freeze_cz = .false.
  ! diff_top_km above the photosphere: the full-diffusion weight falls from 1
  ! at z = 0 to 0 at z = diff_top_above_km (cos^2), so it stays a band at the
  ! CZ top when the atmosphere is in the domain
  real(dp) :: diff_top_above_km = 50.0_dp

  namelist /cz3d/ H_c, lambda_c_km, H_ext, lambda_ext_rsun, cond_r_rsun, alpha_cl, use_atmosphere, use_conduction, use_losses, use_heating, z_atm_km, w_atm_km, H_peak, z_heat_km, &
      T_floor, rho_floor, f_cfl_cond, dt_rad_factor, atomic_dir, coronal_file, freeze_cz, diff_top_above_km, grid_file, init_file, table_file, opacity_file, out_dir, n_r_cells, n_theta, n_phi, nb_r, mr_levels, mr_extrapolate, mr_relevel, mr_relevel_dt, mr_by_shell, qconv_face, nb_theta, &
      nb_phi, cfl, t_end, dt_max, n_step_max, n_log, dt_prof, plot_dt, use_radiation, use_convection, &
      use_hyperbolic, alpha_mlt, L_bottom_Lsun, grey_tau_above, mlt_switch, mlt_R0, mlt_p, use_diffusion, &
      use_minmod, diffusion_h, diff_top_km, diffusion_heating, wall_bottom, wall_top, ds1_random, c_target, rsst_q, rsst_form, yy_margin, refine_h, patch_level, patch_lat_deg, patch_lon_deg, patch_radius_deg, probe_lat_deg, probe_lon_deg, probe_radius_deg, probe_cells, g_from_background, restart_dt, dt_budget, restart_from, restart_domain_tol, restart_new_output, wall_time_max, wall_time_margin, wall_check_steps, r_top_rsun, top_energy, cool_width_km, cool_width_hp, dt1, prof_every, dt1_ramp_steps, grey_top_km, rsst_depth_km, dump_shells, omega_rot, coriolis, lz_correct_steps, lz_correct_dt, use_B, alpha_glm, glm_ch, glm_cfl_frac, b_diff_factor, va_max, b_bc_bottom, b_bc_top, b_seed, b_seed_G, b_seed_r1_rsun, b_seed_r2_rsun, export_state, import_dir, plot_blocks

contains

  subroutine read_params(file)
    character(*), intent(in) :: file
    type(arbor_config_t) :: cfg
    integer :: ios
    character(len=256) :: msg
    call arbor_config_load(cfg, file)
    call arbor_read_params(cfg)
    read(cfg%lines, nml=cz3d, iostat=ios, iomsg=msg)
    call arbor_config_check(cfg, 'cz3d', ios, msg)
    if (wall_bottom /= 'face' .and. wall_bottom /= 'cell') call arbor_abort('wall_bottom must be face or cell')
    if (wall_top /= 'face' .and. wall_top /= 'cell' .and. wall_top /= 'outflow') &
        call arbor_abort('wall_top must be face, cell or outflow')
    if (rsst_form /= 'hotta2012' .and. rsst_form /= 'hotta2015') call arbor_abort('rsst_form must be hotta2012 or hotta2015')
    if (dt1 > 0.0_dp .and. dt_budget > 0.0_dp) call arbor_abort('dt1 (local time stepping) does not conserve energy: set dt_budget = 0')
    if (use_B .and. dt_budget > 0.0_dp) call arbor_abort('the energy budget does not include the magnetic field: set dt_budget = 0')
    if (use_B .and. rsst_form /= 'hotta2012') call arbor_abort('use_B needs rsst_form = hotta2012')
    if (b_bc_bottom /= 'conductor' .and. b_bc_bottom /= 'radial') call arbor_abort('b_bc_bottom must be conductor or radial')
    if (b_bc_top /= 'conductor' .and. b_bc_top /= 'radial') call arbor_abort('b_bc_top must be conductor or radial')
    if (wall_time_max > 0.0_dp .and. wall_time_margin >= wall_time_max) call arbor_abort('wall_time_margin must be below wall_time_max')
    if (wall_check_steps < 1) call arbor_abort('wall_check_steps must be >= 1')
    if (coriolis /= 'split' .and. coriolis /= 'rk') call arbor_abort('coriolis must be split or rk')
    if (glm_ch /= 'fast' .and. glm_ch /= 'cfl') call arbor_abort('glm_ch must be fast or cfl')
    if (b_seed /= 'none' .and. b_seed /= 'toroidal') call arbor_abort('b_seed must be none or toroidal')
    if (top_energy /= 'eddington' .and. top_energy /= 'cooling') call arbor_abort('top_energy must be eddington or cooling')
    ! (local time stepping, dt1 > 0, works with the atmosphere: relaxation, e.g. of the wind)
    if (use_atmosphere .and. dt_budget > 0.0_dp) call arbor_abort('use_atmosphere needs dt_budget = 0')
    if (use_atmosphere .and. top_energy /= 'eddington') call arbor_abort('use_atmosphere needs top_energy = eddington')
  end subroutine read_params

end module cz3d_params_mod
