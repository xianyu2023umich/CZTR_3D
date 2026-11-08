# Growing near-surface horizontal flow in the 1.0 R run (open issue)

Recorded 2026-10-17. Status: open, postponed by the user. Run: output/full100
(runs/full100*.nml; 1.0 R, static background, grey top, subgrid MLT R0 = 1,
p = 4, RSST anchored at 500 km with q = 0, Rempel diffusion with heating,
stress-free top wall at the top cell's centre, Yin-Yang margin 0, 1.875 deg).

## Symptom

After the deep convection saturated (step ~50k), the horizontal flow in the
top ~2 Mm (MLT-dominated, convectively stable, low density) keeps growing,
almost purely horizontal (v_r ~1 m/s there):

| step | rms v_h at 0.1 / 1 / 5 / 10 Mm [m/s] | max abs v_h | max true Mach | min rho1/rho0 |
|---|---|---|---|---|
| 60k | 88 / 113 / 229 / 297 | 0.66 km/s | 0.064 | -2.1% |
| 70k | 231 / 259 / 337 / 309 | 1.33 | 0.183 | -4.1% |
| 80k | 418 / 408 / 376 / 314 | 2.23 | 0.305 | -5.8% |

Exponential at first, roughly linear (+18 m/s per 1000 steps) by 70–80k; it
spread from the top ~1 Mm down to ~5 Mm. 10 Mm and deeper are steady.

## Measured

- Not local time stepping: 2000 steps from 70k with dt1 = 0 and dt1 = inf
  grow identically (runs/full100_test_lts0, _ltsinf).
- Not a large-scale flow: rigid rotation and l = 1 poloidal each <= 0.2% of
  the shell's v_h energy; correlation length a few cells (5–10 deg)
  (scripts/analyze_shells.py; dump_shells option).
- Edge-enhanced: rms v_h within ~4 deg of the own Yin-Yang patch edge is
  ~2x the patch interior above ~2 Mm (12 km depth at 80k: 572 vs 307 m/s);
  no edge preference at 5 Mm and deeper. Interior grows too.
- Horizontal kinetic energy budget (time accurate, 70k–73k, budget_kh.dat,
  scripts/analyze_kh.py): in every shell the change equals the horizontal
  pressure work W_h = -v_h.grad_h p1 (within a few %); the radial flux of
  K_h and the exchange with vertical motion are ~2% each; the numerical
  diffusion removes ~80% of W_h. Top 6 Mm: W_h 4.5e-3 Lsun, diffusion
  -3.7e-3, net +1.0e-3.

## Interpretation (not proven)

The deep convection's pressure fluctuations (cells >= 20 Mm) reach the thin,
low-density top layer and push it sideways (acceleration grad p1 / rho). The
layer is stable (MLT removes the superadiabaticity), so the flow cannot turn
over; MLT carries heat but no momentum (no turbulent friction); only the
numerical diffusion damps it. The Yin-Yang edge enhancement suggests part of
the driving may be spurious (interpolation of steeply stratified p near the
patch edges, margin 0).

## Options discussed

- A. Drag in the top ~2–5 Mm (non-physical, simple). ~25–30% more
  dissipation than the numerical diffusion would balance W_h.
- B. MLT eddy viscosity nu_t ~ v_mlt H_p f_mlt (~1e12–1e13 cm2/s at the top):
  horizontally negligible against the numerical diffusion (~1e14–1e15) at
  23 Mm cells; vertically ~10 h over 2 Mm, may limit but likely not stop it.
- C. Check whether the driving is partly spurious: Yin-Yang margin > 0
  (new grid, restart from scratch) or the ghost interpolation of p1 (ARBOR).
- No-slip top for v_h: ill-posed for the inviscid equations (only v_n = 0
  is allowed at a wall; v_t = 0 overdetermines and makes a one-cell numerical
  boundary layer, a resolution-dependent drag). Well posed only with an
  explicit viscosity; then stress-free at the top is still the natural choice.
- Resolving supergranulation (0.47 deg) would thin the MLT layer and make the
  forcing physical but not remove the need for subgrid momentum transport
  (granules ~1 Mm need ~0.03 deg).

## Related

Also open at 80k: L_top = 0.959 < 1 (the static model's surface under grey
radiation); ~0.04 Lsun accumulates in the interior (L at 100 Mm 1.09).
