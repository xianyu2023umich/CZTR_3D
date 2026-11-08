# Local time stepping and horizontal refinement: an instability (2026-10-25/26)

## Symptom

`c196_m0` (1.0 R, 1.875 deg, 196 radial cells, uncapped local time stepping
`dt1 = 1e30`) refined once horizontally (`c196_m0_r1`, 0.94 deg) blows up
(NaN) within ~230 steps. The motion grows in the top layers: an oscillation
with a period of roughly 5-8 min, peaking at ~170 km depth, about 3 cells
across, propagating upward at 0.4-1 Mm (v_r and rho1 in phase), growing ~3%
per step; no Yin-Yang-edge preference. The MLT turns the large amplitude into
NaN (with frozen MLT it grows but does not blow up).

## Tests (refined restarts from c196_m0, 200-500 steps)

| test | result |
|---|---|
| no hyperbolic MLT | same growth |
| refine from 40k / 70k / 150k | all grow (40k from ~200x smaller seeds): linear |
| maximum diffusion (`use_minmod = .false.`) | stable |
| half the time step | the same growth **per step**: numerical |
| global stepping (`dt1 = 0`) | stable: **local time stepping is the cause** |
| tsc bump near the top flattened; global steps in the top 3 Mm | no effect / delayed (starts at the jump) |
| tsc removed from the hydro / transport / diffusion update separately | **hydro without tsc: stable**; the others unstable |
| radial transport off; diffusion heating off; top wall on a face | same growth |
| tsc frozen in time | same growth |
| horizontally uniform state plus noise | local stepping grows, global slowly: lives on the mean stratification |
| per-sub-step variances | global: hydro x1.085 / diffusion x0.92 per step (grid-scale content); local: a radially smooth mode the limited diffusion barely damps |

Global stepping is stable only through the diffusion: without diffusion it
blows up too (grid-scale noise). A 1D linear model (the code's stencils, walls,
RK4, tsc; one horizontal wavenumber) is dominated by grid-scale zig-zag modes
without damping and has neutral smooth modes with and without tsc; the limited
diffusion is nonlinear (minmod), so the mechanism is not settled by a matrix
model. Working picture: the per-shell scaling of the hydro update makes the
deep layers' clock run faster than the surface's; the coupled top-layer
acoustic mode that results is radially smooth and escapes the limiter.

## The stable range shrinks with resolution

Maximum step factor (`dt1` cap or `tsc**beta`; the outcome follows the maximum,
not the profile's shape):

| max factor | 0.94 deg | 0.47 deg |
|---|---|---|
| ~10 (`dt1 = 20 s`) | stable | marginal (levels off ~40 m/s) |
| ~20 (`dt1 = 40 s`) | marginal (levels off 50-65 m/s) | NaN in ~190 steps |
| ~37 (uncapped) | NaN in ~230 steps | - |

The safe factor roughly halves per refinement (as the numerical diffusion):
~2-3 at 0.234 deg, ~1 at 0.117 deg (extrapolated). Local time stepping helps at
the coarse levels only; the fine levels need (nearly) global steps.
