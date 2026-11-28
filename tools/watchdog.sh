# Hang watchdog for the job scripts (sourced; bash). A run that stops making
# progress but keeps its ranks alive (2026-11-28, casc047_l4, job 3767221:
# one node leader lost an MPI_Allreduce result and the job idled ~2.5 h)
# burns the whole allocation without output. The watchdog watches
# out_dir/log.txt (written every n_log steps): if it has not changed for
# CZ3D_WATCHDOG_S seconds (default 1200; 0 = off), it writes
# out_dir/run_status.txt ('stopped_wall <restart>' with the newest complete
# restart, so the next tools/job.stampede3 continues from it) and kills the
# run. Before the run's first log line (import, PFSS/DPF setup) it allows
# CZ3D_WATCHDOG_START_S more (default 1800). Choose CZ3D_WATCHDOG_S well
# above the time between log lines (n_log steps) and a restart write.
#
#   source tools/watchdog.sh
#   ibrun ./cztr3d.x run.nml & pid=$!
#   cz_watchdog "$out_dir" $pid & wd=$!
#   wait $pid; kill $wd 2>/dev/null
#
# A restart counts as complete if it has the size of the one before it (a
# hang during a restart write leaves a short file); with no older restart to
# compare, the newest is used.

cz_newest_restart() {
  local d=$1 a b
  a=$(ls -t "$d"/restart_*.rst 2>/dev/null | sed -n 1p)
  b=$(ls -t "$d"/restart_*.rst 2>/dev/null | sed -n 2p)
  [ -n "$a" ] || return 1
  if [ -n "$b" ] && [ "$(stat -c %s "$a")" != "$(stat -c %s "$b")" ]; then echo "$b"; else echo "$a"; fi
}

cz_watchdog() {
  local out_dir=$1 pid=$2 lim=${CZ3D_WATCHDOG_S:-1200} start=${CZ3D_WATCHDOG_START_S:-1800}
  local poll=${CZ3D_WATCHDOG_POLL:-60} log=$1/log.txt t0 now m ref rst
  [ "$lim" -gt 0 ] 2>/dev/null || return 0
  t0=$(date +%s)
  while sleep "$poll"; do
    kill -0 "$pid" 2>/dev/null || return 0
    now=$(date +%s)
    m=$(stat -c %Y "$log" 2>/dev/null || echo 0)
    ref=$m
    [ "$m" -lt "$t0" ] && ref=$((t0 + start))       # no log line from this run yet
    [ $((now - ref)) -gt "$lim" ] || continue
    rst=$(cz_newest_restart "$out_dir")
    echo "job $SLURM_JOB_ID: $(date '+%H:%M:%S') watchdog: $log unchanged for $((now - ref)) s (limit $lim s): killing the run"
    if [ -n "$rst" ]; then
      printf 'stopped_wall %s\nwatchdog: no log for %s s (job %s)\n' "$rst" $((now - ref)) "$SLURM_JOB_ID" > "$out_dir/run_status.txt"
      echo "job $SLURM_JOB_ID: watchdog: run_status.txt -> stopped_wall $rst"
    else
      echo "job $SLURM_JOB_ID: watchdog: no restart file to continue from"
    fi
    kill -TERM "$pid" 2>/dev/null; sleep 60
    if kill -0 "$pid" 2>/dev/null; then
      echo "job $SLURM_JOB_ID: watchdog: the launcher did not stop: scancel"
      scancel "$SLURM_JOB_ID"
    fi
    return 0
  done
}
