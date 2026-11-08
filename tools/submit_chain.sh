#!/bin/bash
# Submit a chain of N jobs from the login node, each starting after the
# previous one ended (any state); the job script continues the run from its
# last restart (out_dir/run_status.txt) or exits once the run is finished.
#
#   bash tools/submit_chain.sh <job script> <namelist> <N> [sbatch options]
#   e.g. bash tools/submit_chain.sh tools/job.stampede3_spr4 runs/x.nml 3 -N 6 --ntasks-per-node 96
set -e
job=${1:?job script}; nml=${2:?namelist}; n=${3:?number of jobs}; shift 3
dep=""
for i in $(seq 1 "$n"); do
  id=$(sbatch --parsable $dep "$@" "$job" "$nml")
  echo "job $i of $n: $id"
  dep="--dependency=afterany:$id"
done
