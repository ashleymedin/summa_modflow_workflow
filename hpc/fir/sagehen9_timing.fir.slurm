#!/bin/bash
#SBATCH --job-name=sagehen9-timing
#SBATCH --account=rrg-mclark
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --mem-per-cpu=3G
#SBATCH --time=4:00:00
#SBATCH --output=%x-%j.out

# Times the 3396-HRU grid (test_mflow/run_sagehen9.sh: 72 h from 7 Feb 2017) through the MPI coupler,
# which splits the 9 GRUs; one rank is the time of one calibration trial. Runs in a copy under $RUNS/timing;
# the first pass also runs the case's three-week spin-up, which is not timed.
#   sbatch sagehen9_timing.fir.slurm [rank counts]      # default "1 3 5 9"
set -euo pipefail

RANKS=${1:-"1 3 5 9"}

source "${WORKFLOW:-$HOME/models/summa_modflow_workflow}/hpc/fir/env_sagehen.fir.bash"
ulimit -s unlimited

BIN="$SUMMA_ROOT/bin"
T="$RUNS/timing"
mkdir -p "$T"
rsync -a "$WORKFLOW/test_mflow/coupler_commands.sh" "$WORKFLOW/test_mflow/domain_sagehen9" "$WORKFLOW/test_mflow/ex-gwf-sagehen" "$T/"
cd "$T"

HEADS=domain_sagehen9/simulations/spinup/heads_latflow.bin
if [ ! -f "$HEADS" ]; then
  echo "spin-up..."
  ./coupler_commands.sh -c domain_sagehen9/summa_modflow6.config -w "$HEADS" ex-gwf-sagehen \
      domain_sagehen9/settings/SUMMA/fileManager_latflow_spinup.txt "$BIN/summa_modflow6_sundials_mizuroute.exe" > spinup.log 2>&1
fi

# the MPI coupler takes no restart flags, only the namelist's
sed "s|^ */ *$|  head_restart_read = '../$HEADS'\n/|" domain_sagehen9/summa_modflow6.config > mpi.config
cd ex-gwf-sagehen
echo "ranks  seconds for 72 h"
for n in $RANKS; do
  start=$(date +%s)
  srun --mpi=pmix -n "$n" "$BIN/summa_modflow6_mpi_sundials_mizuroute.exe" \
       ../domain_sagehen9/settings/SUMMA/fileManager_latflow.txt ../mpi.config > "../event_${n}ranks.log" 2>&1
  echo "$n  $(( $(date +%s) - start ))"
done
