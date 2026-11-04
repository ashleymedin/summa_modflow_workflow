#!/bin/bash
#SBATCH --job-name=sagehen9-timing
#SBATCH --account=rrg-mclark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=%x-%j.out

# Times the 3396-HRU grid (test_mflow/run_sagehen9.sh: 72 h from 7 Feb 2017) at several GRU thread counts,
# for the time of one calibration trial and the speedup threads give. Runs in a copy under $RUNS/timing;
# the first pass also runs the case's three-week spin-up, which is not timed.
#   sbatch sagehen9_timing.fir.slurm [thread counts]      # default "1 2 4 5 8"
set -euo pipefail

THREADS=${1:-"1 2 4 5 8"}

source "${WORKFLOW:-$HOME/models/summa_modflow_workflow}/hpc/fir/env_sagehen.fir.bash"
ulimit -s unlimited

EXE="$SUMMA_ROOT/bin/summa_modflow6_serial_sundials_mizuroute.exe"
T="$RUNS/timing"
mkdir -p "$T"
rsync -a "$WORKFLOW/test_mflow/coupler_commands.sh" "$WORKFLOW/test_mflow/domain_sagehen9" \
         "$WORKFLOW/test_mflow/ex-gwf-sagehen" "$T/"
cd "$T"

run() {
  ./coupler_commands.sh -c domain_sagehen9/summa_modflow6.config \
                        -r domain_sagehen9/simulations/spinup/heads_latflow.bin \
                        -s domain_sagehen9/settings/SUMMA/fileManager_latflow_spinup.txt \
                        ex-gwf-sagehen domain_sagehen9/settings/SUMMA/fileManager_latflow.txt "$EXE"
}

if [ ! -f domain_sagehen9/simulations/spinup/heads_latflow.bin ]; then
  echo "spin-up on $SLURM_CPUS_PER_TASK threads..."
  OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK run > spinup.log 2>&1
fi

echo "threads  seconds for 72 h"
for n in $THREADS; do
  start=$(date +%s)
  OMP_NUM_THREADS=$n run > "event_${n}threads.log" 2>&1
  echo "$n  $(( $(date +%s) - start ))"
done
