#!/bin/bash
# Submit short test runs to Stampede3's skx-dev queue when it has idle nodes
# (often it has many, with nobody waiting, while spr has a queue). Run on the
# login node from the CZTR_3D directory:
#
#   bash tools/submit_dev.sh --check                    # skx-dev state only
#   bash tools/submit_dev.sh runs/a.nml [runs/b.nml ...]
#   bash tools/submit_dev.sh -N 12 -t 01:30:00 runs/a.nml
#
# Options:
#   -N n       nodes (default: all idle skx-dev nodes, at most 16)
#   -t time    wall time hh:mm:ss (default and max 02:00:00)
#   -j script  job script (default tools/job.stampede3_seq: the namelists run
#              one after another; tools/job.stampede3 continues a run from its
#              run_status.txt)
#   -n         dry run: print the sbatch command, don't submit
#   -f         submit even if no node is idle (the job then waits)
#
# skx-dev (TACC docs, Nov 2026; check with qlimits): SKX nodes, 48 cores and
# 192 GB each, at most 16 nodes and 2 h per job, 2 running jobs and 4 submitted
# per user, 1 SU per node-hour. Per core SKX is ~1.35x slower than SPR, so 16
# skx nodes (768 ranks) run about as fast as 6 spr nodes x 96 ranks (576) at
# 16 instead of 12 SU per hour: worth it when spr has a queue.
#
# The job's time limit reaches cz3d (CZ3D_WALL_TIME_MAX), so a run that needs
# more than 2 h stops cleanly with a restart, wall_time_margin (namelist)
# before the end: a margin of 1800 s leaves 1.5 h of a 2 h job.
#
# The executable must run on Skylake: ARBOR built with INTEL_ARCH=-xCORE-AVX512
# (ARBOR's Makefile defaults to -xSAPPHIRERAPIDS, which does not run on skx).
# This script checks ARBOR_FFLAGS in the arbor.mk that config.mk names.
set -e
nodes=""; time_lim="02:00:00"; job="tools/job.stampede3_seq"; dry=0; force=0; check=0
while [ $# -gt 0 ]; do
  case "$1" in
    -N) nodes=$2; shift 2 ;;
    -t) time_lim=$2; shift 2 ;;
    -j) job=$2; shift 2 ;;
    -n) dry=1; shift ;;
    -f) force=1; shift ;;
    --check) check=1; shift ;;
    -h|--help) sed -n '2,/^set -e/p' "$0" | sed '$d'; exit 0 ;;
    -*) echo "submit_dev: unknown option $1"; exit 1 ;;
    *) break ;;
  esac
done

# skx-dev state
idle=$(sinfo -h -p skx-dev -t idle -o %D 2>/dev/null | awk '{s += $1} END {print s + 0}')
total=$(sinfo -h -p skx-dev -o %D 2>/dev/null | awk '{s += $1} END {print s + 0}')
pend=$(squeue -h -p skx-dev -t PD 2>/dev/null | wc -l | tr -d ' ')
mine_run=$(squeue -h -u "$USER" -p skx-dev -t R 2>/dev/null | wc -l | tr -d ' ')
mine_all=$(squeue -h -u "$USER" -p skx-dev 2>/dev/null | wc -l | tr -d ' ')
echo "skx-dev: $idle of $total nodes idle, $pend jobs waiting; yours: $mine_run running, $mine_all submitted (limits 2 / 4)"
if [ $check -eq 1 ]; then
  sinfo -p skx-dev -o '%P %a %l %D %t' 2>/dev/null
  squeue -p skx-dev -o '%.10i %.9u %.3t %.10M %.10l %.5D %R' 2>/dev/null | head -20
  exit 0
fi

[ $# -ge 1 ] || { echo "usage: bash tools/submit_dev.sh [-N n] [-t hh:mm:ss] [-j job] [-n] [-f] <namelist> [...]"; exit 1; }
for nml in "$@"; do [ -f "$nml" ] || { echo "submit_dev: no namelist $nml"; exit 1; }; done
[ -f "$job" ] || { echo "submit_dev: no job script $job"; exit 1; }
[ -f config.mk ] || { echo "submit_dev: no config.mk (run from CZTR_3D after ./configure)"; exit 1; }

# the build must run on Skylake
mk=$(sed -n 's/^ARBOR_MK *= *//p' config.mk | head -1)
flags=$(sed -n 's/^ARBOR_FFLAGS *:*= *//p' "$mk" 2>/dev/null | head -1)
case "$flags" in
  *SAPPHIRERAPIDS*|*xHost*|*xhost*|*march=native*)
    echo "submit_dev: ARBOR was built with '$flags' ($mk): it does not run on skx."
    echo "  Rebuild with -xCORE-AVX512 (runs on skx, icx and spr):"
    echo "    make -C ../ARBOR install MODE=opt COMPILER=intel INTEL_ARCH=-xCORE-AVX512 && make clean && make"
    exit 1 ;;
  "") echo "submit_dev: warning: no ARBOR_FFLAGS found in '$mk'; make sure the build has no SPR-only flags" ;;
esac

# wall time: at most 2 h
secs=$(echo "$time_lim" | awk -F: '{s = 0; for (i = 1; i <= NF; i++) s = 60*s + $i; print s}')
if [ "$secs" -gt 7200 ]; then echo "submit_dev: -t $time_lim is over skx-dev's 2 h"; exit 1; fi

# nodes: all idle ones (at most 16) unless -N was given
if [ -z "$nodes" ]; then
  nodes=$idle; [ "$nodes" -gt 16 ] && nodes=16
  if [ "$nodes" -lt 1 ]; then
    [ $force -eq 1 ] || { echo "submit_dev: no idle skx-dev node (use -N n -f to queue anyway)"; exit 1; }
    nodes=16
  fi
fi
if [ "$nodes" -gt 16 ]; then echo "submit_dev: -N $nodes is over skx-dev's 16"; exit 1; fi
if [ "$nodes" -gt "$idle" ] && [ $force -eq 0 ]; then
  echo "submit_dev: $nodes nodes asked, $idle idle: the job would wait (use -f to queue anyway)"; exit 1
fi
[ "$mine_run" -ge 2 ] && echo "submit_dev: warning: you already run 2 skx-dev jobs: this one waits"

cmd="sbatch -p skx-dev -N $nodes --ntasks-per-node 48 -t $time_lim $job $*"
echo "$cmd   # $((nodes*48)) ranks, at most $(awk -v n=$nodes -v s=$secs 'BEGIN {printf "%.1f", n*s/3600}') SU"
[ $dry -eq 1 ] && exit 0
$cmd
