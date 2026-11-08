# Near-top energy non-conservation without subgrid MLT (open issue)

Recorded 2026-10-13. Status: open, postponed. Working belief (the user's):
a resolution issue. The chain of cause and effect is not established.

## Symptom

At 0.99 R (1.875 deg Yin-Yang, 30 km–11 Mm radial cells, 2012 RSST,
cooling top 1000 km, MLT off), the relaxed horizontally averaged total
luminosity is not flat:

| run | base | 150 Mm | 100 Mm | 60 Mm | 30 Mm | min–max |
|---|---|---|---|---|---|---|
| deep_c2 (0.99 R, q = 0, no MLT) | 0.82 | 0.53 | 0.56 | 0.59 | 0.62 | 0.52–1.23 |
| lts_c2_inf_mlt (0.99 R, q = 0, subgrid MLT) | 0.93 | ~1.0 | 0.99 | ~1.0 | 1.0 | 0.88–1.00 |
| vanilla096 (0.96 R, q = 1, no MLT) | 0.93 | 0.94 | 0.97 | 0.97 | 1.04 | 0.90–1.10 |
| mlt096 (0.96 R, q = 1, subgrid MLT) | 0.92 | 0.95 | 1.00 | 1.04 | 1.01 | 0.88–1.04 |

With the bottom (1 Lsun in) and the cooling layer (1 Lsun out) fixed, a
steady state needs the scheme's spurious sources to sum to zero:
L(r) = 1 + (net spurious energy below r). A source near the top therefore
comes with an equal sink somewhere; the deep deficit mirrors the top excess.

## What is measured

- Per-shell budget, deep_c2 days 34–60 (q = 0): sink -0.41 Lsun in the
  bottom 3 cells, source ~+0.8 spread over 150–10 Mm, a dipole at the
  cooling layer.
- q = 1 (uniform xi, deep_c2_q1): base sink gone (+0.004); bulk source
  unchanged (~+0.7) although xi^2 there dropped ~44x. So the base sink is
  the RSST (mass sloshing between shells with depth-varying xi^2 weight),
  the bulk source is the advective scheme. deep_c2_q1 was not relaxed (15 d
  after the switch), so which sink balances the source at q = 1 is unknown.
- Kinetic/thermal split (deep_c2_split, days 60–61.5, budget_split.dat,
  scripts/analyze_split.py): the kinetic-energy equation closes to +0.17
  Lsun (mostly top 15 Mm; W = +2.0 Lsun). The thermal side is swamped by
  -(xi^2-1)(h0+Phi) drho from radial mass sloshing (±500 Lsun over 1.5 d);
  needs a perturbation-form thermal budget (h1, background terms cancelled
  analytically) to be readable.
- Mean mass flux (output/massflux, scripts/compare_massflux.py): net
  <rho v_r> ≤ 2–3% of <rho|v_r|> except the top 3 cells (+0.33, -0.77,
  -1.0 without MLT). Energy via -v_r ds0/dr: ~0.02 Lsun. Ruled out.
- Up/down exchange <rho|v_r|>: does not scale with the error (the MLT run
  exchanges up to 1.5x more mass deep and is flat there). Ruled out as a
  simple proportionality.

## Near-top features, no MLT vs MLT (relaxed, top 7.1–12.4 Mm)

- L_total: dips to 0.75–0.8 at 10–12 Mm, 1.17 at 7.9 Mm (inside the cooling
  layer); with MLT 1.00 in every cell.
- Carriers: resolved enthalpy 0.14–0.78, numerical diffusion up to 0.36 in
  the cooling layer; with MLT, MLT 0.48–0.94, resolved ≤ 0.06.
- rms(s1) 1.6–1.9e5 erg/g/K (20–25x MLT's ~8e3).
- rms v_h ~800 m/s (MLT ~335), rms v_r 30–90 m/s (MLT 2–37); max |v_h|
  2.4 km/s (0.74), max |v_r| 820 m/s (370).
- Max reduced Mach 1.20 (MLT 0.37): beyond the RSST validity limit 0.7.
- Mean v_r in the top cells -35, -19, +13.5 m/s (MLT ≤ 1 m/s).
- Mean T ~0.14% lower; min rho1/rho0 -3.5e-3 (MLT -2e-4).
- <rho|v_r|> at 8 Mm 5x larger.

Picture: with the top H_p (~1–2 Mm) far below the horizontal cell (23 Mm),
the cooling layer drives sharp, strong structures that the grid cannot
resolve; most of the outgoing energy passes through the numerical
diffusion. Subgrid MLT carries that flux instead, and the error disappears
(also deep, which is not understood).

## Suspects not yet tested

- Advective-form entropy transport (-v.grad s1) of the sharp near-top
  structures, together with the Rempel diffusion acting on them.
- RSST beyond its validity (reduced Mach > 1) near the top.

## Next diagnostic when resumed

Perturbation-form thermal budget per shell (h1 = T s1 + p1/rho in the
energy change and the flux; (h0+Phi) drho, rho v_r T0 ds0/dr and the
background-force work handled analytically), from deep_c2's day-60
restart, ~2000 steps; compare with the same from lts_c2_inf_mlt.
