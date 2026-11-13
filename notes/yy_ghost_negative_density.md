# real1875_relax_cont NaN: a Yin-Yang ghost extrapolated to negative density (2026-11-13)

## Summary

The NaN at the transition-region base comes from ARBOR's Yin-Yang ghost fill,
not from the physics. `run_stencil` extrapolates a ghost value from a one-sided
slope when its source cell is the last cell of the source block, so it never
looks across that block edge. In the steep transition region that undershoots
every real neighbour. Here it drove one ghost cell's density below zero over
40 steps. The diffusion operator then computed the sound speed from that ghost,
got NaN, and put it into the flux of the cell next to it. Turning the friction
off at the restart only triggered it, by letting the transition-region
structure steepen faster.

Proposed fix (ARBOR, a few lines): when a stencil offset points past the source
block's edge, do not extrapolate (zero slope, or clip to the cells used).

## Symptom

`runs/real1875_relax_cont.nml` restarts `real1875_relax` from step 800513 with
the friction (`v_damp_*`) dropped. It reaches NaN at step 800554 (41 steps):

    cz3d: NaN on rank 6: cell 246 (z 1620.5 km), lat -31.09, lon 54.30, patch 2, level 0, j k in block 1 1

Earlier tries in the same series also crashed on a patch boundary: the third try
at a Yin patch corner (lat -44, lon 134) and the l = 96 try on a Yang patch edge
(see the header of `runs/real1875_relax.nml`).

## Reproduction

The run is deterministic. A rerun reproduces the crash cell and step exactly,
and its log lines at steps 800520 and 800540 match the original's.

- `runs/real1875_relax_cont_plot30.nml` is the cont namelist with its own
  `out_dir`, a log line every step, and a Tecplot file at step 800543 (30 steps
  in). Its output is in `output/real1875_relax_cont_plot30/`.
- With `n_step_max = 60` it crashes at step 800554, as the original did.
- About 80 s on 8 ranks: `mpirun -np 8 ./cztr3d.x runs/real1875_relax_cont_plot30.nml`

## The crash cell

- It is cell (i, j, k) = (22, 1, 1) of a radial block that holds cells 225-252,
  so i = 22 is global cell 246.
- That block is the Yang patch's corner block. Its patch-local position is
  θ = 45.94° (first cell, edge at 45°) and φ = -134.06° (first cell, edge at
  -135°).
- The cell stays calm up to the step before the NaN. At step 800553:
  ρ = 3.15e-15, T = 2.1e5 K, |v| ≈ 23 km/s, Mach 0.3, β ≈ 0.28.
- At step 800554 it is the only NaN cell within 6° and cells 215-300. No probed
  cell is extreme one step earlier.

## Which operator makes the NaN

We instrumented a copy of `cz_advance` with a NaN scan after every operator. In
step 41, the state is clean after:

- the RK stages,
- `atm_floors`,
- `cz_transport_all` (MLT and radiation),
- `atm_block` (conduction, losses and heating).

The first NaN appears after **`diffusion`**, in variable 1 (rho1) of cell
(22, 1, 1).

`diffusion` (`src/cz3d_app_mod.f90:2023`) first calls `derived(ib, .true.)`,
which sets ρ, T and c_s over the whole ghost range. Its face stencils
(`faces_a`, `faces_b`) then use two ghost cells across each face, and no
diagonal ghosts. A NaN in rho1 therefore needs a bad face ghost.

## The bad ghost

The ghosts are filled with `faces_only=.true.` (lines 1033, 1082, 1091). In the
corner block:

- **The diagonal ghosts, (j, k) = (0, 0), (-1, 0), (0, -1) and (-1, -1), are
  NaN throughout the run.** They are never filled, and no operator reads them,
  so they are harmless.
- **The φ-face ghost (j, k) = (1, 0) is wrong**, at every fill, in the
  transition-region cells. Below are the ghost and its neighbours at the first
  step (RK stage 1 fill; ρ0(246) = 4.70e-17):

  | i (global) | ghost (1,0) | ghost (2,0) | cell (1,1) |
  |---|---|---|---|
  | 1-15 (230-239) | within a few % of the neighbours | | |
  | 20 (244) | 8.0e-15 | 1.18e-14 | 1.14e-14 |
  | 21 (245) | 2.3e-15 | 4.6e-15 | 4.5e-15 |
  | 22 (246) | **1.0e-15** | 2.9e-15 | 2.9e-15 |

- **It drifts down smoothly, by about 2.5e-17 per step.** Its rho1 goes from
  1.02e-15 at the restart to -6.5e-18 at step 41's RK stage 1.
- **At the post-RK fill of step 41 it reaches -5.37e-17 < -ρ0**, which is a
  negative density.
- **That ghost maps well inside Yin**, at Yin θ = 120.0° and φ = 53.4°, more
  than 8° from Yin's edge. The Yin cells around that point have
  ρ = 2.4-3.1e-15. No cell within 6° has a rho1 history that matches the
  ghost's.

## Why: one-sided extrapolation in `run_stencil`

This is the stencil entry for that ghost (dumped from
`arbor_ghost: iStSrc_DI, StOff_DI`):

    dst (22, 1, 0)  src leaf 291, rank 1, src cell (22, 8, 5), offset (0, +0.4906, -0.0166)

Source cell j = 8 is the source block's last cell in θ (nb(2) = 8), so
`ARBOR/src/arbor_ghost.f90:582` takes the one-sided branch:

```fortran
else if (c(d) + 1 > nb(d)) then
  um = src(c - e, v)
  s0 = limited(u0 - um, um - src(c - 2*e, v))   ! both differences behind the cell
...
val(v) = val(v) + s0*StOff_DI(d, m)             ! pushed +0.49 cells across the edge
```

The slope comes from j = 6, 7 and 8 only, and is extrapolated across the block
boundary. The next Yin block's j = 1 is never read. In the transition region
this row has a density minimum right at the block edge:

| step | u(j=6) | u(j=7) | u(j=8) | next block j=1 | ghost predicted by the formula | ghost ARBOR wrote |
|---|---|---|---|---|---|---|
| 800514 | 8.41e-15 | 5.33e-15 | 2.44e-15 | 2.86e-15 | 1.018e-15 | 1.018e-15 |
| 800530 | 9.57e-15 | 6.01e-15 | 2.42e-15 | 2.92e-15 | 6.756e-16 | 6.757e-16 |
| 800540 | 1.08e-14 | 6.61e-15 | 2.39e-15 | 2.96e-15 | 3.268e-16 | 3.270e-16 |
| 800550 | 1.21e-14 | 7.09e-15 | 2.37e-15 | 3.00e-15 | 5.855e-17 | 5.877e-17 |
| 800553 | 1.24e-14 | 7.21e-15 | 2.37e-15 | 3.01e-15 | -6.67e-18 | -6.45e-18 |

All values are rho1 at i = 246. The "ARBOR wrote" column is the ghost at the
next step's RK stage 1 fill. The formula's prediction (minmod; the small k
offset included) matches ARBOR's values to 3-4 digits. As mass piles up at j = 6
and 7, the backward slope steepens and the extrapolated ghost crosses zero. The
true neighbour, at 3.0e-15, stays positive throughout.

### Why patch boundaries, and why this place

- Only Yin-Yang ghosts are interpolated. Box ghosts inside a patch are copies,
  so this can only happen on a patch boundary.
- The patch corner is not special in itself. What matters is a ghost whose
  source cell is at its block's edge, with an offset near ±1/2 toward that edge,
  in steep structure. In the transition region, density changes by 3-5× from
  one column to the next.
- In the Yang-to-Yin mapping, roughly every source-block edge along the
  interface has such ghosts.
- Possibly related, not checked: CZ_3D's mhd4 NaN, near the Yin corner at
  lat -45, lon 135.

## Role of the friction

- With the main run's friction schedule kept (`v_damp_n = 5`,
  `v_damp_steps = 863`, `v_damp_ramp = 3000`, i.e. the same `n_adv` phase), the
  same restart gets past step 800554 with no NaN (90 steps tested).
- Without friction, max |v_h| rises from 36 to 150 km/s within 30 steps, which
  piles up mass in the j = 6-7 columns.
- So the friction slows the steepening, but does not remove the flaw. We did
  not check how far the ghost drifts with friction on.

## Fix options (none applied)

1. **ARBOR `run_stencil`: do not extrapolate past the source block's edge**
   (recommended). In the one-sided branches, use s0 = 0 when the offset points
   toward the missing side: `c(d) = nb(d)` with offset > 0, or `c(d) = 1` with
   offset < 0. Alternatively, clip `val` to the range of the cells used. Either
   way the result stays bounded by real values and keeps positivity. It is first
   order in those cells only.
2. **ARBOR: a central slope from the source block's ghost across that edge**,
   as round 2 does for prolongation. This is more accurate, but needs the
   source's own face ghosts filled first (a second round or an ordering change).
3. **Application safety net**, which hides the symptom rather than fixing it:
   positivity or a floor on density in Yin-Yang ghosts before `diffusion` and
   other operators read them, or a guard in `derived` for the ghost range.

After a fix, rerun `runs/real1875_relax_cont_plot30.nml` with
`n_step_max = 60`. The run should get past step 800554, and the ghost (1, 0)
should stay near 2.4-3.0e-15.

## Side finding: a CZTR_3D restart is not an exact continuation

Restarting from step 800513 with the main run's friction schedule does not
reproduce the main run:

| step | main run | restart (friction schedule matched) |
|---|---|---|
| 800520 | dt = 6.77e-2, Mach 12.6 | dt = 6.80e-2, Mach 11.4 |
| 800600 | Mach 15.3 | Mach 2.82 |

The dt differs already in the first logged step. Some state that the time step
(or the local step factors) depends on is not saved in the restart, or is
rebuilt differently. This is unrelated to the NaN. Not investigated further.

## Debug tools (scratchpad, not in the repo)

These are in
`/private/tmp/claude-503/-Users-xianyu-Documents-FALCON/c2625cc0-9bf3-4400-b123-78c26c8e15b9/scratchpad/`:

- `cztr_dbg/` is a copy of the CZTR_3D source with:
  - `probe_write` extended to every block and the radial range PROBE_ILO to
    PROBE_IHI (environment variables), writing the full state per cell;
  - a NaN scan after every operator in `cz_advance` (`dbg_nan`);
  - a dump of the corner block's stencil entries (`dbg_ghosts`).

  Run it as
  `PROBE_ILO=215 PROBE_IHI=300 mpirun -x PROBE_ILO -x PROBE_IHI -np 8 <that>/cztr3d.x <nml>`,
  with `probe_lat_deg = -31.09, probe_lon_deg = 54.30, probe_radius_deg = 6`.
- `cztr_crash*/` holds the probe output, and `probe_load.py` reads it.
