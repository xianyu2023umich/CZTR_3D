# CZTR_3D

**Work in progress (started 2026-11-07): the convection zone plus the
transition region and corona up to ~2.5 R_sun in 3D. The plan is in
`notes/plan.md`. CZTR_3D began as a copy of CZ_3D (commit e0703a4; smoke
checksum identical); the text below still describes the CZ_3D part.**

Global 3D simulation of the solar convection zone, from its base (0.715 R_sun)
up to a top radius of choice (0.96 R_sun ... the photosphere at 1.0 R_sun), on
the ARBOR block library's Yin-Yang grid. The radial direction uses the 1D
merged grid of CZTR (30 km cells near the photosphere, ~11 Mm at the base);
a radial column lives in one block or in several radial blocks (`nb_r`).

Physics: perturbations (rho1, v, s1) about a fixed 1D background (well
balanced: the background force is dropped), the eos_hhe equation of state
(runtime table), reduced speed of sound (Hotta et al. 2012 form), Rempel-type
slope-limited diffusion with the dissipated energy returned as heat, radial
grey Eddington radiation and a hyperbolically relaxed mixing-length flux per
column (optionally as a subgrid model that switches off where the grid resolves
convection), a cooling layer or a radiating top, local time stepping for
relaxation (`dt1`), rotation (Coriolis force, exact split rotation), and ideal
MHD with GLM divergence cleaning (`use_B`). The parameters are documented in
`src/cz3d_params_mod.f90`; physics notes, tests and open issues in `notes/`.

## Layout

| path | contents |
|---|---|
| `src/` | the code (`cz3d_main.f90`, `cz3d_app_mod.f90`, `cz3d_params_mod.f90`, background, transport, MLT, grey radiation, opacity) |
| `runs/` | namelists of the runs made so far (the record; see below) |
| `data/` | inputs: radial grid, opacity table, EOS table parameters; `make data` generates the rest |
| `tools/` | the static envelope (background) generator: `static_bg.f90`, `build_static_init.py`, copies of CZ's envelope code in `tools/envelope/` |
| `scripts/` | analysis and plotting (Python 3 with numpy and matplotlib) |
| `notes/` | the MHD and rotation options, open issues, test results |

## Building on a new machine

Dependencies: an MPI Fortran compiler (Fortran 2018; gfortran or Intel), and two
libraries from the same project:

- **eos_hhe** (https://github.com/xianyu2023umich/EOS_for_FALCON): the equation of state.
- **ARBOR** (https://github.com/xianyu2023umich/ARBOR): the block library.

Each code is configured once per machine (`./configure`), which chooses the
compiler and writes a `config.mk` that the Makefile reads. All three must use
the same compiler family (Fortran module files of two compilers do not mix);
CZ_3D's configure takes the compiler from ARBOR's install and refuses an
eos_hhe build of another family.

```bash
# 1. eos_hhe, the OpenMP build (library and table generator in build-omp/)
git clone https://github.com/xianyu2023umich/EOS_for_FALCON.git EOS/eos_hhe
cd EOS/eos_hhe && ./configure && make omp && cd ../..

# 2. ARBOR
git clone https://github.com/xianyu2023umich/ARBOR.git ARBOR
make -C ARBOR install MODE=opt

# 3. CZ_3D (defaults: ../ARBOR/install/opt and ../EOS/eos_hhe/build-omp)
cd CZ_3D
./configure            # options: ./configure --help
make

# 4. generated inputs (not kept in the repository):
#    data/eos_table.bin   from eos_hhe's eos_table_gen and data/eos_table.nml (89 MB, ~10 s)
#    data/static_init.dat the static envelope (tools/static_bg.f90 + tools/build_static_init.py)
make data
```

Radial grids: `data/grid_full.dat` (CZTR's, 111 cells to the photosphere) is
the default. Other grids come from CZTR's cell-density map
(`tools/grid/grid_mod.f90`: 5 cells per H_p, a 30 km floor, and optionally
radial cells no taller than `aspect` x the horizontal spacing for
`horiz_deg`): a parameter file `data/grid_<name>.nml` gives the grid and its
init file with `make data/static_init_<name>.dat`; a run names both
(`grid_file`, `init_file`). `data/grid_015a2.nml` (0.15 deg, aspect 2: 141 cells to the photosphere) and
`data/grid_015a1.nml` (aspect 1: 196 cells, the flagship cascade's grid). `n_cz`
sets the number of cells below the photosphere exactly (the photosphere on face
`n_cz`; e.g. `data/grid_015a1_n192.nml`: 192 = 8 radial blocks of 24).

Radial blocks and multirate: `nb_r` cells per radial block (0: the whole
column; it must divide the radial cells) and `mr_levels > 1`: ARBOR's
multirate (rate levels of ratio 2 from each block's stable dt, time
accurate, Berger-Oliger order; at most `mr_levels` levels). Split columns
need `grey_top_km` (the implicit grey layer inside the top block); multirate
needs global steps (`dt1 = 0`) and no energy budget or L_z correction.
192 cells as 8 x 24 at 1.875 deg: 6 levels, 3.4x faster than single rate;
time averages agree with single rate to 1e-4 in L. Split and whole columns
are bitwise equal only without FMA contraction (`-ffp-contract=off` in both
codes; else 1 ulp at the block tops), with global and with local time
stepping (`dt1`).

Changing the radial block size of a running case (restart files hold their
block layout): export, then import with the new `nb_r`.

```bash
# 1. a 0-step run from the old restart, writing the state and the time/step/clocks
#    (the old namelist + restart_from, n_step_max = 0, export_state = .true.)
#    -> out_dir/state_export_rNNNN.bin, out_dir/state_export_meta.bin
# 2. the new run: import_dir = that out_dir, nb_r = 28 (must divide the exported
#    radial block size), refine_h / patch_level as before, no restart_from
```

The import refines the tree to the exported blocks' levels, reads only each
rank's own blocks, and restores the time, step, profile counter and shell
clocks (bitwise the same state; the convective flux restarts at its targets).

Check the build with the smoke test (`runs/smoke.nml`: the 1.0 R setup for 300
steps from the static start; 8 ranks, about a minute):

```bash
make smoke             # ends with: state checksum: 38440707489373
```

On a cluster (no MPI on login nodes), submit it as a short job instead, e.g.
`sbatch -N 1 -n 8 -t 00:15:00 tools/job.stampede3 runs/smoke.nml`.

The checksum is bitwise only with the same compiler, flags and math library
(gfortran, ARBOR's opt flags, macOS with Accelerate's array `exp`/`log` in
`src/vmath_accelerate.f90`; other systems use plain loops,
`src/vmath_plain.f90`); elsewhere compare the log instead. Results do
not depend on the number of ranks.

### Compilers

ARBOR's install sets the MPI wrapper and the optimization flags (its
`arbor.mk`). `-march=native` and `-x...` target one CPU: build on (or for) the
compute nodes.

The Intel compilers are a second family; Intel builds go to their own
directories, so both families can be kept side by side. On Stampede3 (TACC;
Intel MPI, whose `mpi_f08` module gfortran cannot read):

```bash
cd EOS/eos_hhe && ./configure --fc=ifx && make omp       # build-omp-intel/
make -C ARBOR install MODE=opt COMPILER=intel            # ARBOR/install/opt-intel
cd CZ_3D && ./configure --arbor=../ARBOR/install/opt-intel && make   # build-intel/, cztr3d-intel.x
```

`make data` needs only one family (the generated inputs are the same files).
`config.mk` also records the batch launcher (`ibrun` where present, else
`srun`) and the modules loaded at configure time, which `tools/chain.slurm`
loads; configure in the environment the jobs should run in. Run with
`ulimit -s unlimited` (automatic arrays are on the stack).

## Running

```bash
mpirun -np N ./cztr3d.x runs/<case>.nml
```

Output goes to `out_dir`: `log.txt`, horizontally averaged profiles
(`profile_NNNN.dat`), restarts (`restart_NNNNNNNN.rst`, readable on any number
of ranks), optional Tecplot files, rotation maps, shell dumps. Continue a run
with `restart_from`. A hydrodynamic restart can be continued with `use_B`
(the field and the seed are added). The relaxing convective flux is kept in
the restart file (`qconv_face = .true.`, ARBOR's face state), so a continued
run is bitwise the uninterrupted one; restart files from before 2026-10-26
have no face state: read them with `qconv_face = .false.` (the flux then
restarts at its targets).

### Long runs on a workstation (macOS)

Use one rank per performance core (`sysctl -n hw.perflevel0.physicalcpu`;
efficiency cores slow every step down to their pace), with a number of
blocks divisible by the ranks (72 blocks at 1.875 deg: 8, 9, 12, 18, 24).
Start the run detached and kept awake:

```bash
nohup caffeinate -is mpirun -np 8 ./cztr3d.x runs/<case>.nml > output/<case>.out 2>&1 &
disown
```

`caffeinate -is` keeps the machine from sleeping while the run lasts; a
locked screen does not stop it, logging out or a restart does (continue from
the last restart then). `disown`: the terminal can be closed. Status:

```bash
tail -3 output/<case>/log.txt | cut -c1-150   # step, t, dt, Mach, |v_r|, |v_h|, L_top, L range
pgrep -fl cztr3d.x | head -2                    # still running?
cat output/<case>/run_status.txt              # at the end: finished ...
pmset -g log | grep -E " Sleep | Wake " | tail -3   # did the machine sleep?
```

### Batch jobs and wall-time limits

`wall_time_max` (seconds, from the program's start) makes a run stop cleanly
with a restart file once the elapsed wall time plus `wall_time_margin`
reaches it (checked every `wall_check_steps` steps; the margin must cover one
check interval plus writing a restart). `out_dir/run_status.txt` then reads
`stopped_wall <restart>` (or `finished ...`). `tools/chain.slurm` is a SLURM
template that runs a job, reads this file and submits the next job from the
named restart until the run is finished:

```bash
sbatch tools/chain.slurm runs/<case>.nml
```

`tools/job.stampede3` (TACC) and `tools/job.aces` (TAMU) are the same, ready
to use: queues, cores per node, build notes. `tools/bench/submit.sh` submits
the strong-scaling benchmark (`tools/bench/bench192.nml` on several node counts).

The older namelists in `runs/` (before 2026-10-14) take their initial state
from CZTR's output (`../CZTR/output/...`), which is not part of this
repository; they are kept as the record of those runs. New runs use the static
background in `data/` by default.

## Analysis and plots

The scripts in `scripts/` need Python 3 with numpy and matplotlib (with
Homebrew's Python: `python3 -m venv ~/venv-plot && ~/venv-plot/bin/pip install
numpy matplotlib`, then `~/venv-plot/bin/python`). Each prints its usage in
its header; they write PNG files (default: into the run directory).

| script | what |
|---|---|
| **one run** | |
| `plot_run.py <run>` | overview: log time series and mean profiles |
| `plot_L_evolution.py <run> [every n] [out]` | luminosity profiles by component (radiation, MLT, resolved, numerical diffusion), coloured by step |
| `plot_mach_evolution.py <run> [every n] [out]` | rms Mach against depth over time, the log's maxima against step |
| `plot_flux_components.py <run> [first profile] [out]` | luminosity by component, averaged over profiles |
| `plot_rotation.py <run> [first] [last] [out]` | differential rotation and meridional flow (omega_rot /= 0) |
| `estimate_steps.py <run.nml> [key=value ...]` | steps to relax a run, estimated before running it |
| **comparisons** | |
| `compare_evolution.py <run> <run> ... [max=] [out=]` | log time series and space-time maps of L_total of several runs against step |
| `compare_Lprofiles.py <A> <B> [interval] [max step] [out]` | L_total profiles of two runs at the same steps |
| `compare_relaxed.py <A> a0 a1 <B> b0 b1 [out]` | relaxed states (step windows) of two runs |
| `plot_flux_before_after.py out.png <A> a0 a1 <B> b0 b1` | flux components and velocities of two runs over day windows |
| `compare_lts.py <lts run> <reference>` | local time stepping against a time-accurate run |
| `compare_massflux.py out.png <dir>:<label> ...` | mass exchange between up- and downflows |
| `compare_1d.py <3D run> <1D run>` | horizontal means against a CZTR (1D) run |
| **diagnostics** | |
| `analyze_shells.py <dump dir>` | horizontal flow on spherical shells (`dump_shells` files) |
| `analyze_kh.py <run> [depth]` | horizontal kinetic energy budget per shell (`budget_kh.dat`) |
| `analyze_split.py <run>` | where the hydro operator creates or destroys energy |
