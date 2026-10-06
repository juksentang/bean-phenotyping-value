#!/bin/bash
# Run a shard of the Layer-1 sweep inside ONE Slurm allocation: many single-core cells side by side.
#
#   slurm/run_cells.sh <n_reps> <shard_index> <n_shards>
#
# Cells are 0..N-1 (N from config/sensitivity_grid.yaml); shard i takes ids with id % n_shards == i.
# Cells whose output already exists are skipped, so a resubmitted shard resumes where it stopped.
# Output: layer1_genetic_sim/outputs/$TOOL_SWEEP_SUBDIR (default from config: sweep_v2).
set -uo pipefail
NREPS=${1:?n_reps}; SHARD=${2:?shard_index}; NSHARD=${3:?n_shards}
JOBS=${SLURM_CPUS_PER_TASK:-$(nproc)}
SUB=${TOOL_SWEEP_SUBDIR:-$(Rscript -e 'cat(yaml::read_yaml("config/sensitivity_grid.yaml")$output_subdir)' 2>/dev/null | tail -1)}
OUT=layer1_genetic_sim/outputs/$SUB
NCELL=$(Rscript -e 'source("layer1_genetic_sim/R/07_param_sweep.R"); cat(nrow(expand_sweep_grid(load_sweep_grid())))' 2>/dev/null | tail -1)
mkdir -p "$OUT" slurm/logs/$SUB
echo "shard $SHARD/$NSHARD: $NCELL cells in grid, $NREPS reps, $JOBS concurrent, → $OUT"

run_one() {  # one cell, single core (mclapply reads SLURM_CPUS_PER_TASK)
  local id=$1 f
  f=$(printf '%s/summaries/cell_%05d.csv' "$OUT" "$id")  # written last, so present = cell complete
  [ -s "$f" ] && return 0
  # 40-min cap per cell (normal ~12.5 min, max seen 20): a hung AlphaSimR
  # process is killed instead of blocking its slot until the job times out
  SLURM_CPUS_PER_TASK=1 timeout 2400 Rscript layer1_genetic_sim/R/07_param_sweep.R "$id" "$NREPS" \
    > "slurm/logs/$SUB/cell_$id.log" 2>&1 || { echo "FAIL cell $id"; return 1; }
}
export -f run_one; export OUT NREPS SUB TOOL_SWEEP_SUBDIR=$SUB

# Two passes: the second retries cells that failed or timed out in the first
for pass in 1 2; do
  seq 0 $((NCELL - 1)) | awk -v s="$SHARD" -v n="$NSHARD" '$1 % n == s' \
    | xargs -P "$JOBS" -I{} bash -c 'run_one {}'
  RC=$?
  [ $RC = 0 ] && break
  echo "pass $pass: some cells failed, retrying"
done
echo "shard $SHARD done: $(ls "$OUT"/summaries/cell_*.csv 2>/dev/null | wc -l) cells complete; xargs rc=$RC"
exit $RC
