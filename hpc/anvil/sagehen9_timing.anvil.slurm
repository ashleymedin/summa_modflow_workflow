#!/bin/bash
#SBATCH --job-name=sagehen9-timing
#SBATCH --account=ees250064
#SBATCH --partition=wholenode
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=1:30:00
#SBATCH --output=%x-%j.out

# Times the 3396-HRU GWE grid (run_sagehen9_gwe.sh: 72 h from 7 Feb 2017): the MPI coupler, which splits
# the 9 GRUs, on 1 to 9 ranks alone on the node, then COPIES serial runs at once, as a full calibration node.
# Expects its spin-up (heads_gwe.bin, spinup_gwe_restart_*) in $RUNS/timing/domain_sagehen9/simulations/spinup.
set -eo pipefail

source "${WORKFLOW:-$HOME/models/summa_modflow_workflow}/hpc/anvil/env_sagehen.anvil.bash"
ulimit -s unlimited
COPIES=${1:-121}
BIN="$SUMMA_ROOT/bin"
T="$RUNS/timing"
mkdir -p "$T/domain_sagehen9/simulations/run_1/SUMMA"
rsync -a "$WORKFLOW/test_mflow/domain_sagehen9" "$WORKFLOW/test_mflow/ex-gwf-sagehen-gwe" "$T/"

FM=../domain_sagehen9/settings/SUMMA/fileManager_gwe.txt
CF=../domain_sagehen9/summa_modflow6_gwe.config
HR=../domain_sagehen9/simulations/spinup/heads_gwe.bin

# the MPI coupler takes no restart flags, only the namelist's
sed "s|^ */ *$|  head_restart_read = '$HR'\n/|" "$T/domain_sagehen9/summa_modflow6_gwe.config" > "$T/mpi.config"

echo "MPI ranks  seconds for 72 h"
cd "$T/ex-gwf-sagehen-gwe"
for n in 1 3 5 9; do
  s=$(date +%s.%N)
  mpirun -np $n --bind-to none "$BIN/summa_modflow6_mpi_sundials_mizuroute.exe" "$FM" "$T/mpi.config" > "$T/mpi_$n.log" 2>&1 || echo "failed, see mpi_$n.log"
  echo "$n  $(echo "$(date +%s.%N) - $s" | bc)"
done

echo "$COPIES serial copies at once"
L=$(mktemp -d /tmp/sagehen9.XXXXXX)
for i in $(seq "$COPIES"); do
  mkdir -p "$L/$i/domain_sagehen9/simulations/run_1/SUMMA"
  ln -s "$T/domain_sagehen9/simulations/spinup" "$L/$i/domain_sagehen9/simulations/spinup"
  for f in "$T"/domain_sagehen9/*; do [ "$(basename "$f")" = simulations ] || ln -s "$f" "$L/$i/domain_sagehen9/"; done
  cp -r "$T/ex-gwf-sagehen-gwe" "$L/$i/"
done
export L FM CF HR BIN
s=$(date +%s.%N)
srun -n "$COPIES" --cpu-bind=cores bash -c 'cd $L/$((SLURM_PROCID+1))/ex-gwf-sagehen-gwe; s=$(date +%s); /usr/bin/time -f "%e s %M kB" "$BIN/summa_modflow6_sundials_mizuroute.exe" "$FM" "$CF" --head-restart-read "$HR" > run.log 2> time.log || echo "copy failed $SLURM_PROCID"'
echo "all: $(echo "$(date +%s.%N) - $s" | bc) s"
tail -qn1 "$L"/*/ex-gwf-sagehen-gwe/time.log | sort -n | awk '{a[NR]=$0} END{print "fastest", a[1]; print "median ", a[int((NR+1)/2)]; print "slowest", a[NR]}'
cp "$L/1/ex-gwf-sagehen-gwe/run.log" "$T/copy1.log"
rm -rf "$L"
