#!/bin/bash
# Job chain for long cztr3d.x runs on SLURM machines (a template: adapt the
# #SBATCH lines, modules and paths to the machine).
#
#   sbatch tools/chain.slurm runs/<case>.nml
#
# The executable, the MPI launcher (ibrun on TACC) and the modules to load come
# from config.mk (./configure), so the job runs in the environment of the build.
# The namelist should set out_dir, restart_dt > 0 and wall_time_max a bit
# below the job's time limit (in seconds), e.g. for --time=72:00:00:
#   wall_time_max = 259200.0, wall_time_margin = 3600.0
# Each job runs until the wall-time limit, writes a restart and
# out_dir/run_status.txt ('stopped_wall <restart>'), then submits the next job,
# which continues from that restart (restart_from set in a copy of the
# namelist). The chain ends when the status is 'finished' (t_end or
# n_step_max reached) or a job fails. Bound a chained run by t_end:
# n_step_max counts each job's own steps. restart_from must sit on a line of
# its own (or first on its line) in the namelist.
#SBATCH --job-name=cz3d
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --time=72:00:00
#SBATCH --output=cz3d_%j.log

set -e
nml=${1:?usage: sbatch tools/chain.slurm <namelist>}
cd "$SLURM_SUBMIT_DIR"
cfg() { sed -n "s/^$1 *= *//p" config.mk | head -1; }
[ -f config.mk ] || { echo "chain: no config.mk (run ./configure)"; exit 1; }
mods=$(cfg MODULES | tr ':' ' ')
if [ -n "$mods" ]; then module purge; module load $mods; fi
ulimit -s unlimited                      # large automatic arrays on the stack
exe=./$(cfg EXE)
launcher=$(cfg LAUNCHER)

out_dir=$(sed -n "s/^ *out_dir *= *'\([^']*\)'.*/\1/p" "$nml" | head -1)
$launcher "$exe" "$nml"

status=$(head -1 "$out_dir/run_status.txt" 2>/dev/null || echo "missing")
set -- $status
if [ "$1" = "stopped_wall" ]; then
  next="${nml%.nml}_cont.nml"
  # continue from the restart: replace (or add) restart_from, keep the output series
  sed -e "/^ *restart_from *=/d" -e "/^ *restart_new_output *=/d" "$nml" \
    | sed "s#^&cz3d#\&cz3d\n  restart_from = '$2'#" > "$next"
  sbatch tools/chain.slurm "$next"   # not $0: SLURM runs a spooled copy
  echo "chain: continuing from $2 with $next"
else
  echo "chain: done ($status)"
fi
