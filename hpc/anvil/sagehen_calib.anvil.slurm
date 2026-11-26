#!/bin/bash
#SBATCH --job-name=sagehen-calib
#SBATCH --account=ees250064
#SBATCH --partition=wholenode
#SBATCH --nodes=1
#SBATCH --ntasks=121
#SBATCH --time=4-00:00:00
#SBATCH --output=%x-%j.out

# One Sagehen calibration from a run directory made by prepare_sagehen_calib.sh: rank 0 coordinates, every
# other rank runs one trial at a time on one core, so ranks = population + 1. Anvil allows 96 h.
# Submit with submit (env_sagehen.anvil.bash), which keeps us to 8 jobs queued or running.
#   sbatch sagehen_calib.anvil.slurm <run>                                          # pop 120, one node
#   sbatch -p shared --ntasks=41 --time=12:00:00 sagehen_calib.anvil.slurm <run>    # pop 40, charged by the core
set -eo pipefail

RUN=${1:?"usage: sbatch sagehen_calib.anvil.slurm <run name>"}
source "${WORKFLOW:-$HOME/models/summa_modflow_workflow}/hpc/anvil/env_sagehen.anvil.bash"
ulimit -s unlimited
cd "$RUNS/$RUN"
echo "Run $RUNS/$RUN on $SLURM_JOB_NODELIST, $SLURM_NTASKS ranks, $(date)"
mpirun -np "$SLURM_NTASKS" --bind-to none "$SUMMA_ROOT/bin/summa_modflow6_opt_sundials_mizuroute.exe" -c config.toml > run.log 2>&1
echo "End $(date)"
grep -E "^NSGA-II generation" run.log || true
