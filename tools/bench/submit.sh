#!/bin/bash
# Strong-scaling benchmark: tools/bench/bench192.nml on several node counts,
# one job each (run from the CZ_3D directory):
#
#   bash tools/bench/submit.sh [queue] [cores per node] [job script] [node counts] [namelist]
#     Stampede3: bash tools/bench/submit.sh skx 48 tools/job.stampede3 "1 2 4 8"   (default)
#     ACES:      bash tools/bench/submit.sh cpu 96 tools/job.aces "1 2 4"
#     the 0.47 deg state in radial blocks: ... "1 2 4" tools/bench/s3_047_lts.nml
#
# Results: grep -h "blocks on\|done:" output/<name>_N*/stdout.txt (name: the
# namelist's out_dir, e.g. bench192)
q=${1:-skx}; c=${2:-48}; job=${3:-tools/job.stampede3}; nodes=${4:-1 2 4 8}; nml=${5:-tools/bench/bench192.nml}
base=$(sed -n "s/^ *out_dir *= *'\([^']*\)'.*/\1/p" "$nml" | head -1)
for n in $nodes; do
  d=${base}_N$n; mkdir -p $d
  sed "s#$base'#$d'#" "$nml" > $d/bench.nml
  sbatch -J bench$n -p $q -N $n --ntasks-per-node $c -t 00:30:00 -o $d/stdout.txt -e $d/stderr.txt \
    $job $d/bench.nml
done
