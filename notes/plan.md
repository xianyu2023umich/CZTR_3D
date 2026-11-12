# CZTR_3D plan

Started 2026-11-07. Nothing beyond the copy is implemented yet.

## Goal (December 2026)

A working 3D model from the base of the convection zone to the corona at
0.47 deg, with a magnetic field and coronal heating, that runs stably for
days of model time and gives synthetic AIA EUV movies. It is a demonstration
to show people (support, resources), not a publication; robustness and
visuals first, convergence studies later. 0.234 deg only after it works.

## Code base

- A copy of CZ_3D at e0703a4 (Alfven cap va_max, multirate with per-shell
  levels and re-levelling, the job's own wall-time limit); CZ_3D itself is
  not modified. Smoke checksum identical to CZ_3D's (38440707489373).
- To port from CZTR_2D (2.5D, r-theta): the atmosphere (cztr2d_atm_mod: the
  background above the join, background-force and atmosphere switches,
  floors), radiative_loss_mod (optically thin losses, column integrals), the
  prescribed heating, hyperbolic Spitzer conduction, the three-step start
  (scripts/chain.py), its lessons (below). The phi terms come back (CZTR_2D
  dropped them); Yin-Yang stays.
- No separate corona-only code (TR_3D): CZTR's step 2 (the CZ frozen, the
  atmosphere evolving) is the corona-only phase; frozen CZ blocks should do
  no work at all then.

## Agreed decisions (2026-11-07)

1. Folder FALCON/CZTR_3D (this one).
2. Domain to ~2.5 R_sun. The top is source-surface-like: the field radial
   there (b_bc_top = 'radial'). Mass: a closed wall at first (a hydrostatic
   corona on open field; no wind), an outflow boundary later. Coarse radial
   cells far out (dr growing ~ r).
3. The atmosphere's initial field: a PFSS extrapolation of the surface B_r
   (source surface at the top), computed once (Python, spherical harmonics)
   and imported when the atmosphere is added.
4. Coronal heating matters most for the result. Start simple (CZTR's H(z)),
   then a heating function derivable from the run (literature survey,
   notes/heating.md); try AWSoM-like w+- heating this month
   (CZTR_2D/notes/alfven_wave_heating.md) once CZTR_3D runs.
5. Thermal conduction along the field (CZTR_2D's hyperbolic Spitzer
   conduction, projected on b).
6. The Alfven cap varies with height. Idea (the user): set it from a target
   dt: per cell va_max = cfl dx_min/dt_target - (c_s + |v|), with a floor
   (e.g. >= 2 c_s) so the coronal sound waves are not distorted; dt_target
   possibly per multirate level. The target dt is then the cost knob.
7. Cost: set by the transition region/corona cells and the cap; estimate
   once the grid exists.

## Lessons from CZTR / CZTR_2D (keep)

- Steps 1-3 must use one fixed dt: the near-surface equilibrium depends on
  dt through the operator splitting; a dt change at the release kicks it.
- diffusion_heating only in the CZ (and keep the energy flux at the join face).
- Each column's atmosphere pressure-matched to its CZ top.
- The full-diffusion top band (diff_top_km) stays out of the steps with an
  atmosphere (it kicks the photosphere).
- Background: smooth (C1), rho0 <= rho, don't follow the TR (kinks kill runs).
- Per-block multirate levels near the surface are harmful (CZ_3D mhd4): keep
  mr_by_shell.

## Open questions

- The radial grid: the 196 CZ cells unchanged + a TR band (30 km floor) +
  a stretched corona to 2.5 R; total a multiple of the radial block size (28).
- Grey radiation with radial blocks: the implicit grey layer must sit in one
  block (grey_top_km); where does the atmosphere's radiation live?
- PFSS on Yin-Yang: how to import (interpolate the extrapolated field to the
  cell centres, then clean div B with GLM?).
- The heating function: survey in notes/heating.md (first choice: (H_B|B| + H0) exp(-h/lambda), gated on density).
- Synthetic EUV: in-situ optically thin projection with AIA response
  functions (R(T) tables from aiapy/CHIANTI), one MPI reduction per image;
  304 A excluded (optically thick).

## Stages

1. Copy (done 2026-11-07).
2. The extended radial grid to 2.5 R; background above the join (CZTR_2D).
3. Hydro atmosphere at 1.875 deg on the laptop: losses, H(z), conduction
   (isotropic first), the three steps; frozen CZ blocks skipped in step 2.
4. Field in the atmosphere: PFSS import, height-dependent cap (dt-based),
   conduction along b, heating from the survey.
5. 0.47 deg from the latest CZ_3D mhd5 state on Stampede3.
6. Synthetic EUV images and movies.
7. AWSoM-like w+- heating.

## Status

- 2026-11-07: copy of CZ_3D e0703a4; builds; smoke checksum identical.
- 2026-11-07: radial grid data/grid_atm25.dat (tools/grid/extend_grid.py):
  448 cells = 16 x 28. Cells 1..252 (faces 0..252) identical to
  grid_015a1.dat: the 196 CZ cells and the chromosphere to z = 2 Mm (~31 km
  cells, 65 km at 2 Mm). Above: geometric growth (ramping to 5% per cell,
  neighbour ratio <= 1.051) until dr = 0.918 r dtheta (0.469 deg), then dr
  ~ r dtheta to 2.5 R (no extra refinement there, the user). dr: 202 km at
  5 Mm, 950 km at 20 Mm, 4.5 Mm at 100 Mm, 9 Mm at 0.5 R up, 12.7 Mm at
  1 R up. 252 cells above the photosphere (162 above 1.01 R, 87 above 1.3 R).
  (--smax 0.08 gives 420 cells.)
- 2026-11-07: init/background data/init_atm25.dat (tools/build_atm_init.py,
  figure notes/fig/init_atm25.png). CZ cells (<= 196): state and background
  exactly static_init_015a1.dat (CZ_3D restarts keep their meaning).
  Initial atmosphere: CZTR's 1D step-3 state to z = 3 Mm (rho x 1.0071 to
  match our photosphere), then T -> 1.5 MK (C1, w 10 Mm) and hydrostatic
  (mu 0.593) to 2.5 R: rho 9e-18 at 0.29 R up, 2.9e-19 at the top.
  Background above the join: rho0 = rho_A + rho_B; rho_A CZTR's exponential
  (H 119.7 km at the join from the CZ background's slope -> H_inf 103 km by
  bisection, blend 300 km), rho_B hydrostatic at 1.2 MK, 0.1 x the state at
  3 Mm (falls faster than the state: margin if the corona cools). T0
  continues the CZ slope, ~4780 K. rho0 <= rho in every cell (rho/rho0 1.3-2.2
  in chromosphere/TR, 3-10 in the corona); largest second difference of
  ln rho0 per cell 0.100 at 2.2 Mm (CZTR's own background: 0.108), 2e-3 at
  the join.
- 2026-11-08: atmosphere ported (src/cztr3d_atm_mod.f90, commit 8acec12),
  with radial blocks and multirate (CZTR_2D needed whole columns and single
  rate). Choices: conduction faces in r, theta, phi, each block its own copy
  of shared faces (ghosts); chromospheric column mass from p/g, neutral-H
  column within the block + ghost estimate above; grey radiation and MLT
  only in the CZ (the join block is the transport top, Eddington top at the
  photosphere; no grey in the atmosphere, only thin losses: also avoids
  CZTR_2D's grey + C&L double counting); background force above the join;
  diffusion heating only in the CZ; diff_top band above the photosphere
  (diff_top_above_km); freeze_cz. The join must end a radial block (196 =
  7 x 28). No change without use_atmosphere (smoke checksum identical).
  First test (1.875 deg, laptop, frozen CZ from the static start, CZTR's
  heating, 30 coarse steps = 1312 s model, 837 s wall, ~1.6 model s/s):
  8 levels (finest 0.37 s), TR 1.73-1.91 Mm (CZTR 1D 1.72-1.84), coronal
  flows up to ~50 km/s (Mach 0.3) from the non-equilibrium 1.5 MK start,
  rms_h(s1) growing to ~8e5 then flat and horizontal flows ~200 m/s: to
  locate (profiles). Run on the Studio / plugged in: tests take >10 min.
- 2026-11-10: 3D wind relaxation (runs/wind3d_relax.nml, Studio, 1.875 deg,
  20000 LTS steps, 7.8 h on 8 ranks): steady from ~10000 steps (rho v r^2
  flat to +-0.3% above 1.1 R), still spherically symmetric (rms v_r =
  |<v_r>|, horizontal s1 variation <= 2e-4 of the mean, v_h <= 65 m/s).
  It settles to a different steady state than CZTR_2D's 1D wind: 457 vs 600
  km/s at 20 R, sonic point 4.0 vs 3.0 R, n v at 1 AU 3.4e8 vs 2.7e8, T at
  3-5 R 2.3 vs 3.3 MK, base density x2.5, low-corona pressure x2.
  KNOWN ISSUE (kept for now, the user: "not that serious"): the cause is the
  port's choice to stop the grey radiation at the photosphere. In CZTR_2D the
  grey solve continues through the (grey-opaque) chromosphere and heats it;
  in 3D only the thin losses act there, so the low chromosphere (0.2-0.5 Mm)
  cools to T_floor (4400 K, vs 4500-4700 K in 1D), the chromosphere
  restructures, the TR pressure rises and more mass evaporates (Hansteen &
  Leer: the TR pressure sets the mass flux): a denser, cooler, slower wind.
  The 3D wind is arguably closer to the slow wind, and CZTR_2D double counts
  the chromospheric radiation (grey + C&L). Fix later if needed: optically
  thin absorption Q = 4 pi kappa rho (J - B(T)) with J the diluted
  photospheric radiation (no implicit solve across blocks).

- 2026-11-11: b_seed = 'pfss' (src/pfss_mod.f90, cz_seed_pfss in
  cz3d_app_mod): the field above the join from the join cell's B_r (the CZ's
  own field, no dipole added: the user), computed in the run after the blocks
  are set up (fresh start, import, hydro restart). Quadrature fit over the
  Yin-Yang join shell (MPI sums), refined pfss_iter (6) times; l_max default
  pi/(2 dtheta) (48 at 1.875 deg, 192 at 0.47 deg); source surface
  b_seed_rss_rsun, radial above. Checks: runs/pfss_check.nml (CZ set to the
  analytic dipole, atmosphere rebuilt from cell 196 alone): error 7e-8 rms,
  9e-7 max; at l_max = pi/dtheta the midpoint sums alias (1% rms after 2
  passes). Standalone test with random coefficients to l = 40: div B and
  curl B ~1e-8 of |grad B| (FD noise). Import: atmosphere blocks without an
  exported block keep the init file's state (a CZ_3D export under the
  atmosphere). Smoke checksum unchanged.
- 2026-11-12: real field at 1.875 deg (mhd5 coarsened, runs/real1875_relax.nml):
  two NaNs at the TR base. (1) l_max = 96: near-grid-scale PFSS components
  are not discretely curl-free; spurious Lorentz force, Alfvenic flows at ~5
  Mm (143 km/s; 37 at l = 48, 4 at l = 24; not the heating). (2) l_max = 24
  still crashed: the background rho0 hugged the initial TR (rho/rho0 1.15-1.9
  at 1.5-1.7 Mm, from the bisection for the largest H_inf); the field makes
  the TR height vary by column (1330-1975 km) and rho fell to ~0.1 rho0 where
  it sat lower. Fix: build_atm_init.py --tr-margin-km (rho0 below the state
  as if the TR were that much lower, from margin below the TR up):
  data/init_wind20_dip_m600.dat (H_inf 99.5 -> 52.2 km; rho/rho0 >= 26 at
  1.5-1.7 Mm, ~1e3-1e4 in the upper chromosphere; ln rho0 2nd difference
  0.090). Lesson: in 3D with field the background needs a TR margin.
