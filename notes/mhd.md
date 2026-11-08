# Magnetic field option (use_B), added 2026-10-18

Ideal MHD in the perturbation formulation, after solar_model (src/ModAdvection,
ModSource, ModDivB): the field is not split into a background and a
perturbation; it is a state variable of its own.

## Equations (src/cz3d_app_mod.f90: rhs_mhd, curl_block, div_block, glm_damp)

State variables 6–9 when `use_B = .true.`: B_r, B_theta, B_phi [G] and the
GLM scalar psi. ARBOR rotates (B_r, B_theta, B_phi) between the Yin-Yang
frames like the velocity (`arbor_set_vectors([vr_, br_])`).

    dv/dt   += (curl B) x B / (4 pi rho)            (rho = rho0 + rho1, the inertia density)
    dB/dt    = curl(v x B) - grad psi
    dpsi/dt  = -c_h^2 div B,     c_h = c_s0/xi + v_A
    then psi <- psi exp(-dt alpha_glm c_h / dx_min)  (after each step; alpha_glm = 0.5)

c_h (`glm_ch`): 'fast' (default, as solar_model) or 'cfl': the largest speed
the cell's step allows, glm_cfl_frac cfl / [(1/dr + 1/(r dtheta) + 1/(r sin
theta dphi)) dt_cell] (glm_cfl_frac = 0.8), at least c_s0/xi + v_A. It gains
wherever a cell steps below its own limit (deep cells under global stepping:
up to dt_own/dt); little under uncapped local stepping.

Curl and divergence are evaluated in spherical components: radial derivatives
with CZTR's 4th-order stencil and metric (ddr_block), angular ones with ARBOR's
(arbor_deriv_cell), as the hydro terms.

- Time step: the Alfven speed is added to the (reduced) sound speed in
  stable_dt, so local time stepping sees it too (each shell's own step).
- Diffusion: B and psi get the Rempel flux -c J(q) (minmod, h switch), no
  density weight; c includes v_A. With `diffusion_heating`, the face energy
  flux includes B_f . F_B / (4 pi) and the cell's entropy is corrected so the
  dissipated magnetic energy returns as heat (energy conserving like the
  hydro part).
- Walls (`b_bc_bottom`, `b_bc_top`): 'conductor' (B_r odd, B_h even: B_r = 0
  at the wall) or 'radial' (B_h odd: vertical field). Defaults as Hotta et
  al. 2014: conductor at the bottom, radial at the top. psi even. A
  cell-centred wall also zeroes the odd components in its cell after each step.
- Not included: the energy budget (refused with use_B), the hotta2015 RSST
  form (refused), explicit resistivity (only the numerical diffusion).

## Seed field (`b_seed = 'toroidal'`)

B = b_seed_G w(r) sin(theta_g) e_phi_g about the global z axis (divergence
free; vanishes at the poles, unlike a strictly uniform B_phi), w = 1 below
`b_seed_r1_rsun`, 0 above `b_seed_r2_rsun`, sin^2 between (keeps it off a
vertical-field top). Added in init_block, also to an imported state. Hotta et
al. 2014 added a uniform B_phi = 100 G to the relaxed hydro state.

## Starting MHD from a hydro run

Directly from a hydro restart (ARBOR's `arbor_restart_read(..., n_var_file=)`,
added 2026-10-18 in ARBOR; a run may have more variables than its file):
set `use_B = .true.`, a seed, and `restart_from = '<hydro restart>'`. The
hydro variables are read, B and psi come back zero, cz_field_from_hydro adds
the seed to every block; t and the step continue (one continuous record). An
MHD restart (9 variables) gets no seed.

The older route, for runs that are not restart compatible otherwise: a hydro
run (or a 0-step run from a hydro restart) with `export_state = .true.`
writes out_dir/state_export_rNNNN.bin; an MHD run with `import_dir` copies the
common variables by block key (the same grid), adds the seed; t starts at 0.

## Diagnostics

- log.txt: last three columns max|B| [G], E_mag/E_kin (shell means; E_kin
  with rho0), rms(div B dx_min)/rms|B| (ghosts filled at the log step).
- profiles: columns 21, 22 rms(B_r), rms(B_h) [G].
- Tecplot: Br, Bth, Bph added.

## Tests (2026-10-18, 0.99 R deep_c2 setup)

- Hydro unchanged: deep_c2 reproduced exactly at step 100 with the new code.
- use_B with no seed: identical to hydro (B stays 0).
- Toroidal seed 100 G, nearly static: max|B| 99.98 G, rms B_h = 81.5 G
  (= 100 sqrt(2/3)), B_r ~1e-16; div B relative 8e-6 at t = 0, 6.7e-6 after
  100 steps.
- Export from deep_c2 day 60 + import with the seed, 1000 time-accurate steps:
  the imported state matches deep_c2 exactly; the flows stretch the field
  (max |B| 100 -> 182 G, rms B_r up to ~6 G); stable.
  div B relative error after 1000 steps: 8.6e-4 with glm_ch = 'fast' (growing,
  slowly decelerating), 2.1e-4 with 'cfl' (4x lower, clearly levelling off:
  +3.7e-5 per 100 steps at first, +1.3e-5 at the end); the field and flows
  are otherwise the same to the printed digits.
- Direct hydro -> MHD restart (deep_c2 day 60, 1000 steps, glm_ch = 'cfl'):
  bit-identical to the export/import run (same state checksum), t continuing
  from day 60. MHD restart continuity: 500 + restart + 500 steps gives the
  same checksum as 1000 straight. Hydro still reproduces deep_c2. For comparison, GLM in
  solar_model (Cartesian) grows to ~2e-2 and saturates at ~1e-2 (the user's
  test; its normalisation may differ).
- Likely sources of div B: the component-wise Rempel diffusion of B (minmod,
  h switch; not divergence preserving), the Yin-Yang interpolation, the walls.
  The induction term itself does not create it in the interior: both
  derivatives are the same 4th-order centred 5-point stencil in index space
  with one-index metric factors, so div(curl) cancels exactly there.

## Verification of the implementation (2026-10-18, 0.96 R)

1. Energy exchange (300 steps from mhd096 day 210): the Lorentz force's work
   W_L = int v.(J x B)/(4 pi) and the induction term's dE_mag/dt cancel to
   ~1e-4 relative (-7.78e-5 vs +7.78e-5 Lsun; sum <= 4e-9 Lsun), as they
   must with walls that pass no Poynting flux: signs and factors of both
   terms consistent. In the 10 kG test they cancel at 0.34 Lsun to 3e-4.
   (log.txt columns W_Lorentz, dE_mag/dt(induction).)
2. B diffusion off (b_diff_factor = 0, 1000 steps from mlt096 + 100 G seed):
   div B rel 1.4e-4 after 1000 steps vs 1.44e-3 with the diffusion on: the
   component-wise diffusion of B makes ~90% of the div B error; the rest is
   the Yin-Yang interpolation and the walls. The induction/curl
   discretisation itself is clean, as argued.
3. Strong field (10 kG seed, E_mag/E_kin 0.60 -> 0.69 in 4.3 days): stable,
   no NaN; the flow gives 0.34 Lsun to the field; rms velocities 3–6% lower
   than the weak-field run after 4 days (braking starting); flux balance
   unchanged (L 0.91–1.02); div B 1.4e-3, like the weak field.

## Results so far: mhd096 (0.96 R, from mlt096 day 150, 60 days, time accurate)

Toroidal seed 100 G (full below 0.92 R, zero above 0.95 R), glm_ch = 'cfl',
conductor bottom, radial top, output/mhd096 (runs/mhd096*.nml).

| day | max abs B | E_mag/E_kin | div B rel | min / max L |
|---|---|---|---|---|
| 0 | 100 G | 6.03e-5 | 9e-6 | 0.90 / 1.03 |
| 10 | 329 | 7.14e-5 | 1.76e-3 | 0.91 / 1.04 |
| 30 | 507 | 7.50e-5 (peak) | 1.78e-3 | 0.91 / 1.03 |
| 45 | 705 (peak) | 7.31e-5 | 1.93e-3 | 0.90 / 1.03 |
| 60 | 544 | 6.69e-5 | 1.95e-3 | 0.89 / 1.09 |

- Seed stretched (+18% in 10 days), then flat, then slow decay (-11% over
  days 30–60): no dynamo at this resolution (expected: the effective magnetic
  Reynolds number with 1.875 deg cells and the Rempel diffusion of B is far
  below a small-scale dynamo's threshold; Hotta et al. 2014 had one at ~0.09
  deg). Local concentrations up to ~700 G.
- B_h piles up at the conducting bottom (81 -> 124 G rms; downward pumping),
  B_r up to ~30 G rms at 110–165 Mm; field spread into the top layers.
- Convection and flux balance unaffected (E_mag/E_kin ~ 7e-5).
- div B saturated at ~1.7–1.95e-3 (vs ~2e-4 at 0.99 R after 1000 steps; the
  reduced Mach number is ~3x larger at 0.96 R).

## Ideas for later (postponed by the user, 2026-10-18)

- Field validation needs a dynamo-capable resolution (HPC); at 1.875 deg the
  field is a slowly decaying passive field.
- Dynamo threshold test: a factor on the diffusion of B only (e.g. 0.3, 0.1)
  to see whether E_mag then grows: how close the setup is to a dynamo, and
  what resolution would be needed.
- Local stepping with a field: branch one MHD state (a relaxed local-stepping
  hydro state plus the seed) into (1) dt1 ramped to 0 then time accurate, (2)
  dt1 kept at inf. Compare per shell in each shell's own clock (E_mag, rms
  B_r/B_h) and the normalised depth distribution (pumping between shells
  stepping at different rates is where they may differ). Most useful with a
  saturated dynamo: can local stepping reach dynamo saturation faster?
- div B: sources are most likely the component-wise Rempel diffusion of B,
  the Yin-Yang interpolation and the walls; a divergence-aware diffusion
  (Rempel et al. 2009 appendix, as Hotta 2014) or a larger alpha_glm if it
  matters.
- Rotation (large-scale dynamo, differential rotation) is not in the code.
- Coronal coupling later needs the field above the photosphere (xi = 1 there
  already) and the atmosphere physics in 3D.
