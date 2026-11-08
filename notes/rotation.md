# Rotation (omega_rot), added 2026-10-18

Coriolis acceleration -2 Omega x v in the momentum equation, Omega =
omega_rot along the global z axis (Yin's z; Yang's local y, since yy_cart maps
(0,0,1) to (0,1,0)); spherical components per cell from sph_basis, the cross
product in the orthonormal (r, theta, phi) basis (src/cz3d_app_mod.f90:
rhs_coriolis).

Two ways (`coriolis`): 'split' (default since 2026-10-18; the user's idea):
the exact solution of dv/dt = -2 Omega x v, a rotation of v about the axis by
-2 Omega dt/2, before and after each step (Strang; StateOld rotated too, as
the RK stages build on it; each shell by its own step under local stepping):
no stability limit, |v| kept exactly. 'rk': the term in the RK stages (RK4-like:
unstable for 2 Omega dt > ~2.8). No centrifugal term (Omega^2 r / g ~ 2e-5). Solar: 2.6e-6 rad/s
(414 nHz; Carrington sidereal 2.87e-6). Does no work; local stepping scales
it with the rest of the RHS. log.txt: last column L_z = int rho (x cross v)_z
dV of the flow (rotating frame, global axis).

## Tests (0.96 R, mlt096 day 150, 5 days, time accurate)

- Rotation off: deep_c2 reproduced exactly.
- L_z: without rotation it wanders 6.9e42 -> 5.0e42 g cm2/s in 4 days, with
  rotation 6.9e42 -> -2.0e43. Scale: the CZ's moment of inertia ~1e53 g cm2,
  so these are rigid-rotation equivalents of 1e-10–2e-10 rad/s, ~1e-4 of
  Omega: no significant spurious torque (a wrong Omega direction in one patch
  would give torques ~Omega x the convective momentum). The faster change with
  rotation is consistent with the Coriolis torque from mass redistribution,
  -Omega d(I_z)/dt.
- Effect after 5 days: deep rms v_r ~30% lower (120 Mm: 67 -> 45 m/s), v_h
  slightly higher; the deep convective flux drops (L_total 0.70–0.79 at
  90–155 Mm vs ~1.0): rotation suppresses the deep convection's heat
  transport (deep Rossby number ~0.1); the deep layers will heat up and
  re-adjust. Near the top (Ro > 1) little change.

## split vs rk (0.96 R, mlt096 day 150, dt ~330 s)

- Solar (Omega 2.6e-6, 2 Omega dt 1.7e-3), 1300 steps: max Mach and max
  abs v_h agree to 3–4 digits throughout; L_z to ~0.2% after 5 days (the
  remaining difference is the convection's chaotic divergence).
- Omega 1e-3 (2 Omega dt 0.66), 300 steps: both stable, differ by a few %
  (at this rate the flow is strongly disrupted in both; min L < 0).
- Omega 5e-3 (2 Omega dt 3.3): 'split' runs at dt 332 s with normal
  velocities; 'rk' is unstable: spurious growth (max v_h 4–7x larger) until
  the CFL step shrinks to 2 Omega dt ~2.3–3.1, the RK4 stability limit: the
  run limps along with wrong flows. Main-sequence fast rotators (P ~0.1–1 d,
  Omega up to ~7e-4) with deep local-stepping steps (~1000 s) reach 2 Omega
  dt ~0.15–1.5, so 'split' is the robust default.

## Discussion points (2026-10-18)

- Rossby numbers in our runs (v / (2 Omega L)): deep ~0.1 (v ~70 m/s, L ~100
  Mm), mid ~1 (190 m/s, 30 Mm), near top >~3: the deep CZ is rotationally
  constrained, the upper part not.
- Expected: columnar (banana) cells outside the tangent cylinder; Reynolds
  stresses and meridional circulation build differential rotation (low Ro
  favours a solar-like fast equator, high Ro anti-solar; the "convective
  conundrum"); suppressed convective efficiency at low Ro (seen); Taylor-
  Proudman (cylindrical) isorotation unless a latitudinal entropy gradient
  (thermal wind) breaks it, which needs e.g. a subadiabatic layer below the
  CZ (our bottom is a wall at the CZ base).
- Time: differential rotation and meridional flow take many rotation periods
  (28 d) to years; local stepping reaches the same angular-momentum steady
  state in principle (time-scale invariant), but distorts the transient.
- Numerics: the advective momentum form does not conserve angular momentum
  exactly (watch L_z); the numerical diffusion acts as a viscosity that may
  damp differential rotation strongly at 1.875 deg (a spin-down test would
  measure it); MLT is non-rotating (fine where Ro > 1).
- With B: rotation enables a large-scale (alpha-Omega) dynamo; Hotta et al.
  2022 (Nature Astronomy) found the small-scale field needed for a solar-like
  differential rotation at high resolution.

## rot096: 100 days of rotation from mlt096 day 150 (0.96 R, Omega 2.6e-6)

runs/rot096.nml, output/rot096 (rotation_NNNN.dat daily; scripts/plot_rotation.py).

- Thermal re-adjustment: the deep convective flux collapsed first (min L 0.09
  at day 10, L 0.14–0.17 at 120–155 Mm), recovered by day ~40–50; days 60–100
  L flat 0.97–1.04. Deep rms v_r roughly halved (136 Mm: 61 -> 26 m/s; 92 Mm:
  72 -> 39), v_h up (max 480 -> 500–680 m/s).
- Differential rotation (days 220–250 mean, N-S symmetric): ANTI-SOLAR, in the
  top ~0.06 R: equator -12 nHz, latitude 42–58 deg +13 nHz, ~72 deg up to
  +38 nHz (polar bins few cells, noisy); below ~0.9 R nearly rigid (+-3 nHz).
  Days 10–30: equator -20, lat 60 +7 nHz (still developing).
- Meridional flow: one cell per hemisphere, poleward near the top (20–45
  m/s), equatorward return below; solar-like direction, larger amplitude than
  the Sun's ~15–20 m/s.
- Interpretation: anti-solar is the expected regime when the convection that
  carries the angular momentum has Ro >~ 1 (here the top layers, Ro ~1–3):
  the "convective conundrum" in our setup; the deep rotationally constrained
  layer (Ro ~0.1) stays nearly rigid (maybe also the strong numerical
  viscosity).
- L_z drifts linearly, 6.9e42 -> -7.1e44 g cm2/s in 100 days (~-7e42/day):
  a mean rotation change of ~-7e-9 rad/s (~-1.1 nHz per 100 days, ~-4 nHz per
  year) for I ~1e53 g cm2. Small against the +-10–40 nHz differential
  rotation but systematic: either -Omega dI/dt from mass redistribution or
  the scheme's angular-momentum non-conservation (advective momentum form,
  the 2012 RSST mass budget). To check if runs get long.

## L_z budget (2026-10-18; runs/rot096_lz.nml, runs/norot096_lz.nml, 20 days each)

log.txt now also has T_Coriolis = -2 Omega int rho s v_s dV (= -Omega dI_z/dt).
Delta L_z = int T_cor dt + residual (the scheme's torque):

| run | Delta L_z [g cm2/s] | int T_cor dt | residual | residual as a mean rotation change |
|---|---|---|---|---|
| rot096 day 250–270 (Omega 2.6e-6) | -1.45e44 | -4.7e43 (33%) | -9.8e43 (67%) | -2.8 nHz/yr |
| mlt096 day 150–170 (no rotation) | -3.1e43 | 0 | -3.1e43 | -0.9 nHz/yr |

(I ~ 1e53 g cm2.) The residual grows steadily (not noise) and is ~3x larger
with rotation (zonal flows, shear). Sources: the Rempel diffusion acts on rho
v_theta, rho v_phi in their own (patch) frames, not on the angular momentum
density; the advective momentum form; the 2012 RSST mass budget. Compared
with the +-10–40 nHz differential rotation it matters after a few years of
simulated time. Possible fix: remove the non-Coriolis part of each step's
Delta L_z by a rigid-rotation velocity correction delta_Omega x r (global
axis; exact conservation, common in global codes), or angular-momentum-
conserving diffusion fluxes.

## Global L_z correction (lz_correct_steps / lz_correct_dt), added 2026-10-18

The Coriolis impulse is accumulated every step (each shell with its own
step); at a correction the measured L_z minus (reference + impulse) is
removed by a rigid rotation v -= delta_Omega z x x, delta_Omega = delta_L /
I_z; the reference is set at the run's start/restart. Angular-momentum-
conserving diffusion was postponed by the user (it makes more sense once the
advection is conservative too).

Test, rot096 days 250–270 (runs/rot096_lzc.nml, every 20 steps): residual
-7.98e40 g cm2/s (vs -9.75e43 uncorrected), i.e. 1e-4 nHz; L_z now follows
the Coriolis impulse (Delta L_z -4.89e43 vs int T_cor dt -4.89e43). Flows
unchanged (mean max Mach 8.1127e-3 vs 8.1122e-3, max v_h 5.6009e4 vs
5.6006e4, min L 0.981 both); differential rotation over the 20 days the same
(equator -11.1 vs -11.2 nHz, lat 60 +15.4 vs +15.3 nHz). No large-correction
warnings. Recommended for long rotating runs.

## Flux components with rotation, with and without MLT (2026-10-18)

rot096 (MLT) and rotvan096 (no MLT; runs/rotvan096.nml, from vanilla096 day
100, L_z correction every 20 steps), last 40 of 100 days, against their
non-rotating parents (mlt096 days 130–150, vanilla096 days 80–100):

| depth | vanilla096: res / diff | rotvan096: res / diff | mlt096: MLT / res / diff | rot096: MLT / res / diff |
|---|---|---|---|---|
| 196 Mm | 0.10 / 0.07 | 0.04 / 0.17 | 0.00 / 0.12 / 0.03 | 0.00 / 0.03 / 0.19 |
| 145 Mm | 0.86 / -0.11 | 0.34 / **0.38** | 0.01 / 0.92 / -0.17 | 0.10 / 0.24 / **0.42** |
| 105 Mm | 0.92 / -0.02 | 0.77 / 0.10 | 0.06 / 0.92 / -0.05 | **0.35** / 0.53 / 0.04 |
| 75 Mm | 0.96 / -0.03 | 0.95 / -0.02 | 0.30 / 0.73 / -0.01 | 0.49 / 0.52 / 0.00 |
| 53 Mm | 0.96 / 0.01 | 0.95 / 0.02 | 0.69 / 0.34 / -0.01 | 0.73 / 0.30 / -0.01 |
| 38 Mm | 0.59 / 0.36 | 0.63 / 0.31 | 0.84 / 0.09 / -0.02 | 0.83 / 0.10 / -0.02 |

(radiation as in all runs: 0.77 at the base, 0.21 at 145 Mm, ~0 above 100 Mm)

- Bottom ~50 Mm: with rotation, the numerical diffusion carries ~0.4 Lsun at
  145 Mm with or without MLT: the rotationally suppressed convection near the
  base (Ro ~0.1) cannot carry the flux; the steepened mean entropy gradient
  is carried by the Rempel diffusion. Not an MLT effect.
- 60–150 Mm with MLT: rotation moved flux from the resolved convection to MLT
  (105 Mm: MLT 0.06 -> 0.35, resolved 0.92 -> 0.53). MLT is non-rotating, so
  it carries flux there that rotation suppresses for the resolved flow;
  without MLT the resolved flow keeps 0.77 at 105 Mm. A rotation-aware MLT
  (reduced efficiency at low Ro) or a stronger switch-off deep would remove
  this.
- Near the top without MLT: the known numerical-diffusion flux (~0.3 at 38 Mm),
  rotating or not.
- Differential rotation (days 70–100 means, N-S symmetric): both anti-solar.
  No MLT: equator -24 nHz at 30 Mm, -16 at 60 Mm; lat 60 +6 nHz. With MLT:
  equator -12 at 30 Mm, -4 at 60 Mm; lat 60 +13. Without MLT the resolved
  near-top convection is stronger (v_h 293 vs 224 m/s at 38 Mm) and the
  anti-solar equator deficit about twice as large and deeper.

## Rotation rate scan without MLT (2026-10-18; rotvan096_half: Omega 1.3e-6)

Last 40 of 100 days from vanilla096 day 100 (Omega 0: vanilla096 days 80–100):

| depth | resolved (0 / 1.3e-6 / 2.6e-6) | numerical diffusion | rms v_r [m/s] | <s1> change per cell (outward) |
|---|---|---|---|---|
| 164 Mm | 0.61 / 0.37 / 0.19 | -0.04 / +0.26 / +0.41 | 41 / 24 / 16 | +2.8 / -7.3 / -12.9 |
| 145 Mm | 0.86 / 0.60 / 0.34 | -0.11 / +0.16 / +0.38 | 56 / 36 / 25 | +10.9 / -10.5 / -30.4 |
| 128 Mm | 0.92 / 0.79 / 0.52 | -0.08 / +0.06 / +0.27 | 65 / 44 / 32 | +19.3 / -6.8 / -46.5 |
| 105 Mm | 0.92 / 0.91 / 0.77 | -0.02 / +0.01 / +0.10 | 71 / 54 / 42 | +29.1 / +3.5 / -46.4 |
| 61 Mm | 1.01 / 0.99 / 0.97 | -0.05 / -0.03 / -0.02 | 77 / 69 / 60 | -1.4 / -18.4 / -57.4 |

- Monotonic in Omega: slower deep flows, a superadiabatic mean entropy
  gradient that steepens and reaches higher with Omega, and a numerical
  diffusion flux that grows with it (the mechanism: the limiter turns to
  first order wherever fluctuations make local extrema, so the mean gradient
  is diffused with c ~ c_red, independent of the flow). Even at half the solar
  rate the bottom ~50 Mm rely on it (0.16–0.26 Lsun).
- Differential rotation (days 70–100, N-S symmetric), at 30 Mm: half rate eq
  -13, lat 32 +4, lat 58 +9 nHz; full rate eq -24, lat 32 -4, lat 58 +6 nHz.
  Both anti-solar; contrast 22 vs 30 nHz, i.e. 0.11 vs 0.07 of Omega_0: the
  slower rotation (higher Ro) is relatively more anti-solar, as expected.
- The user considered diffusing only s1 - <s1> (no numerical conduction of the
  mean profile) "awkward"; not done.
