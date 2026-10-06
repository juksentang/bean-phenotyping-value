#!/bin/bash
# v2 sweep extension (meta-analysis-driven low err/cost ranges) + h2/GxE robustness.
# Array tasks: 0-5 extA (960 cells), 6-17 extB (1800), 18-21 robust (540).
#SBATCH --account=<your-allocation>
#SBATCH --job-name=tool_sweep_v2_ext
#SBATCH --array=0-21
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=ext_%A_%a.log
module purge >/dev/null 2>&1; module load StdEnv/2023 gcc/12.3 r/4.4.0 >/dev/null 2>&1
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:-$PWD}"
t=$SLURM_ARRAY_TASK_ID
if   [ $t -lt 6 ];  then B=extA;   S=$t;        N=6
elif [ $t -lt 18 ]; then B=extB;   S=$((t-6));  N=12
else                     B=robust; S=$((t-18)); N=4; fi
export TOOL_SWEEP_GRID=config/sweeps/sweep_v2_$B.yaml TOOL_SWEEP_SUBDIR=sweep_v2_$B
echo "task $t → block $B shard $S/$N"
bash slurm/run_cells.sh 100 $S $N
