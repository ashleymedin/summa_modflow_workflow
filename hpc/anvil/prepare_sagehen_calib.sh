#!/bin/bash
# Builds a Sagehen calibration run directory, $RUNS/<run name>, on an Anvil login node.
#
# Usage:  ./prepare_sagehen_calib.sh <lumped|grid> <1|2 layers> <population> <generations> <run name>
# A two-layer run starts from $RUNS/spinup_2layer (spinup_2layer_lumped.anvil.slurm, spinup_2layer_grid.anvil.slurm).
set -eo pipefail

source "$(cd "$(dirname "$0")" && pwd)/env_sagehen.anvil.bash"
CALIB="$WORKFLOW/test_calibration"
PY=${PY:-$HOME/models/venv_mf6build/bin/python}

LAYOUT=${1:?"layout: lumped or grid"}; LAYERS=${2:?"1 or 2"}
POP=${3:?"population"}; GEN=${4:?"generations"}; NAME=${5:?"run name"}
R="$RUNS/$NAME"
[ ! -e "$R" ] || { echo "$0: $R exists"; exit 1; }

mkdir -p "$R/output" "$R/state"
"$PY" "$WORKFLOW/test_mflow/tools/build_sagehen9_calibration.py" "$LAYOUT" "$CALIB/sagehen/forcing" "$R/domain" "$LAYERS"
cp "$R/domain/settings/SUMMA/coldState.nc" "$R/state/"
if [ "$LAYERS" = 2 ]; then
  for k in 1 2; do cp "$RUNS/spinup_2layer/strt_spunup_$LAYOUT$k.txt" "$R/domain/mf6/strt$k.txt"; done
  cp "$RUNS/spinup_2layer/state_spunup_$LAYOUT.nc" "$R/state/coldState.nc"
elif [ "$LAYOUT" = lumped ]; then
  cp "$CALIB/sagehen/strt_spunup_lumped.txt" "$R/domain/mf6/strt1.txt"
  cp "$CALIB/sagehen/state_spunup_lumped.nc" "$R/state/coldState.nc"
fi
cp -r "$CALIB/sagehen/observations" "$R/observations"
cp -r "$WORKFLOW/test_mflow/domain_sagehen9/settings/mizuroute" "$R/mizuroute"
source "$CALIB/sagehen_config.bash"
write_sagehen_config "$R/config.toml" "$R" "$R/observations" "$R/mizuroute" "$((POP * GEN))" "$POP" "$LAYERS"
echo "Ready: $R; run on $((POP + 1)) ranks:"
echo "  sbatch --ntasks=$((POP + 1)) $WORKFLOW/hpc/anvil/sagehen_calib.anvil.slurm $NAME"
