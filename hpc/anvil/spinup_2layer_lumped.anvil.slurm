#!/bin/bash
#SBATCH --job-name=sagehen-spin-lumped
#SBATCH --account=ees250064
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --output=%x-%j.out

# Spins up the lumped two-layer domain: heads (spinup_sagehen9_heads.py), then both models' states over
# cycled water years 2017-2018 (spinup_sagehen9_cycles.py), in $RUNS/spinup_2layer.
set -eo pipefail

source "${WORKFLOW:-$HOME/models/summa_modflow_workflow}/hpc/anvil/env_sagehen.anvil.bash"
ulimit -s unlimited
PY=${PY:-$HOME/models/venv_mf6build/bin/python}
TOOLS="$WORKFLOW/test_mflow/tools"
EXE="$SUMMA_ROOT/bin/summa_modflow6_sundials_mizuroute.exe"
S="$RUNS/spinup_2layer"
mkdir -p "$S"

"$PY" "$TOOLS/build_sagehen9_calibration.py" lumped "$WORKFLOW/test_calibration/sagehen/forcing" "$S/lumped" 2
"$PY" "$TOOLS/spinup_sagehen9_heads.py" "$S/lumped" "$EXE" "$MF6_ROOT/bin/mf6" "$S/strt_heads_lumped1.txt"
for k in 1 2; do cp "$S/strt_heads_lumped$k.txt" "$S/lumped/mf6/strt$k.txt"; done
"$PY" "$TOOLS/spinup_sagehen9_cycles.py" "$S/lumped" "$EXE" "$S/strt_spunup_lumped1.txt" "$S/state_spunup_lumped.nc" 8
