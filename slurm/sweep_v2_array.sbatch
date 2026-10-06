#!/bin/bash
#SBATCH --account=<your-allocation>
#SBATCH --job-name=tool_sweep_v2
#SBATCH --array=0-15
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=1:30:00
#SBATCH --output=arr_%A_%a.log
module purge >/dev/null 2>&1; module load StdEnv/2023 gcc/12.3 r/4.4.0 >/dev/null 2>&1
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
export TOOL_SWEEP_SUBDIR=sweep_v2
cd "${SLURM_SUBMIT_DIR:-$PWD}"
bash slurm/run_cells.sh 100 $SLURM_ARRAY_TASK_ID 16
