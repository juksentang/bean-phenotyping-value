#!/bin/bash
# v3 sweep: all blocks re-centered on the public-data genetic parameters (h2 0.10, kappa 1.5)
# plus the h2 x kappa robustness block. Array tasks:
#   0-15 main (3000 cells), 16-21 extA (960), 22-33 extB (1800), 34-39 gxe (960), 40 stress (6).
#SBATCH --account=<your-allocation>
#SBATCH --job-name=tool_sweep_v3
#SBATCH --array=0-40
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=arr_%A_%a.log
module purge >/dev/null 2>&1; module load StdEnv/2023 gcc/12.3 r/4.4.0 >/dev/null 2>&1
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:-$PWD}"
t=$SLURM_ARRAY_TASK_ID
if   [ $t -lt 16 ]; then B=main;   S=$t;        N=16; SUB=sweep_v3
elif [ $t -lt 22 ]; then B=extA;   S=$((t-16)); N=6;  SUB=sweep_v3_extA
elif [ $t -lt 34 ]; then B=extB;   S=$((t-22)); N=12; SUB=sweep_v3_extB
elif [ $t -lt 40 ]; then B=gxe;    S=$((t-34)); N=6;  SUB=sweep_v3_gxe
else                     B=stress; S=0;         N=1;  SUB=sweep_v3_stress; fi
export TOOL_SWEEP_GRID=config/sweeps/sweep_v3_$B.yaml TOOL_SWEEP_SUBDIR=$SUB
echo "task $t → block $B shard $S/$N → $SUB"
bash slurm/run_cells.sh 100 $S $N
