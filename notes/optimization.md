# Optimization notes

Measured on the T2 test (1.875 deg MHD + atmosphere restart, 20 steps,
8 ranks on the laptop; the end-of-run timer tables: `timers_report`, ARBOR's
`arbor_timer_report`). Single laptop runs drift by ~10% with the machine's
temperature: compare versions back to back, or against an untouched stage.

## Done (2026-11-21, commits daf18d1 .. this one)

- Submodules (setup, diagnostics, synthetic observables); hot path in the parent.
- Diffusion ~2.8x: top-layer blend only where its weights are nonzero, the
  resistive field current from per-block derivative arrays, no string
  compares in face loops, `in_atm` skips, EOS only on the cells used,
  vectorized psi limiter, lean heating EOS.
- RHS ~1.8x: `cz_bind` caches the height switch and sin theta, branch-free
  cap helpers, the Lorentz/psi loop written out (vectorized), `va_cap`,
  no B copies in the curls, `si_coeff` without T and during the RK fills.
- Ghost fills: the transport's fill replaced by a radial copy, no psi in the
  atmosphere's fill; `derived(ib, need, cells)`; the interior EOS (RK
  stages) and the losses (atmosphere) during the fills, with an MPI_Iprobe
  progress poke (`overlap_fills`, default on; neutral on the laptop, to be
  A/B tested on Stampede3).
- Atmosphere ~3x: conduction ingredients per block (switch, dB_r^2, sin
  theta, the cells' T gradients) instead of per face, T^2.5 once, the
  neutral-H term and its tanh skipped where its weight is exactly 0,
  exp(-0.25) at tau's floor; floors once per step with the column EOS.

- stable_dt -30%: the speeds in a vectorized pass (va_cap), atm_dt with the
  column EOS.
- Radial implicit solves: all (j, k) columns of a stack at once (the Thomas
  recursion vectorized across columns), one elimination for v_theta and
  v_phi (the same matrix); in the implicit blocks s1, B, psi updated right
  after the RHS (rho1 and the solved velocities after the solves: the
  operator uses the stage's starting density). ~-27%.
- ARBOR b71f61b installed (Yin-Yang stencils per radial column): ghost fill -38%.

Bit-identity: checked with the smoke test and T2 checksums at every step.
Loops that newly vectorize differ only by FMA contraction (identical when
built with `-ffp-contract=off`).

## Possible next steps (not done)

- **Conduction face loops (atm_block)**: still scalar (per-face calls of
  face_flux, along_b, trac_f; T^2.5, the saturation, exp where tau is above
  its floor). Vectorizing them as the Lorentz loop (the arithmetic written
  out in the loop, branch-free, ivdep) might save 30-40% of the ~1.0 s per
  20 steps (~2% of a step).
- Column solves: the remaining ~1 s per 20 steps is mostly the transposition
  (gathering columns from (i, j, k) into (jk, n) and back, strided); storing
  the implicit tendencies column-major from the start would avoid it.
- RHS leftovers: the p1 shell and the curls' scaled copies cover edges and
  corners (-5% of cells at 0.47 deg); fuse the advection's update passes.
- Overlap for the diffusion's fill (eta, the v_r copy, the interior EOS):
  only if Stampede3 shows the overlaps pay.
- The EOS evaluation itself (eos_hhe's table interpolation): later.
- Alfven cap from the horizontal cells (va_cap_horizontal, 2026-11-22, off):
  the cap's dx_min includes the 30-70 km radial cells, so it sits at its floor
  1.5 c_s everywhere above the join (the Lorentz force reduced by
  (cap/v_A)^2, 1e2-1e4) although the radial waves are semi-implicit (MAS uses
  no cap). With the horizontal cells only (the GLM speed keeps the old cap),
  1.875 tests: switching it on in a relaxed state (T2) releases huge forces
  (Mach 360 in 20 steps); a fresh DPF start also grows strong flows (Mach 52
  in 40 steps at lat -4.7, lon 40.3, z 1.4 Mm) while the old cap stays calm:
  the seeded field is far from force balance at the join, which the old cap
  hides. Needs a relaxed initial field (magnetofriction, or a ramp of the cap)
  before it can be used.
