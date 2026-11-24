#!/bin/bash
#SBATCH --job-name=sagehen-calib
#SBATCH --account=rrg-mclark
#SBATCH --nodes=1
#SBATCH --ntasks=121
#SBATCH --mem-per-cpu=3G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x-%j.out

# One Sagehen calibration from a run directory made by prepare_sagehen_calib.sh: rank 0 coordinates,
# every other rank runs one trial at a time on one core, so ranks = population + 1.
#
#   sbatch sagehen_calib.fir.slurm <run>                  # pop 120
#   sbatch --ntasks=13 sagehen_calib.fir.slurm <run>      # pop 12
set -euo pipefail

RUN=${1:?"usage: sbatch sagehen_calib.fir.slurm <run name>"}

source "${WORKFLOW:-$HOME/models/summa_modflow_workflow}/hpc/fir/env_sagehen.fir.bash"
ulimit -s unlimited

EXE="$SUMMA_ROOT/bin/summa_modflow6_opt_sundials_mizuroute.exe"
cd "$RUNS/$RUN"

echo "Job:      $SLURM_JOB_ID on $SLURM_JOB_NODELIST"
echo "Run:      $RUNS/$RUN"
echo "Ranks:    $SLURM_NTASKS"
echo "Start:    $(date)"

srun --mpi=pmix \
     -n "$SLURM_NTASKS" \
     "$EXE" -c config.toml > run.log 2>&1

echo "End:      $(date)"
grep -E "^NSGA-II generation" run.log || true
