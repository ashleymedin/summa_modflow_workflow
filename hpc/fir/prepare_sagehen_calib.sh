#!/bin/bash
# ---------------------------------------------------------------------------------------
# Builds a Sagehen calibration run directory on this machine, with every path as it will be on Fir,
# for rsync to $RUNS there. The domain, starting state, targets and searched parameters are those of
# test_calibration/test_calibration_sagehen.sh (WY2017 spin-up, WY2018 scored).
#
# Needs a python3 with netCDF4 and pyproj (PYTHON=~/miniconda3/envs/modflow6/bin/python).
#
# Usage:  ./prepare_sagehen_calib.sh <lumped|grid> <population> <generations> <run name> [Fir RUNS dir]
#   e.g.  ./prepare_sagehen_calib.sh grid 120 8 grid_p120 /scratch/amedin/sagehen
# ---------------------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW="$(cd "${HERE}/../.." && pwd)"
CALIB="${WORKFLOW}/test_calibration"
MFLOW="${WORKFLOW}/test_mflow"
PYTHON=${PYTHON:-python3}

LAYOUT=${1:?"usage: $0 <lumped|grid> <population> <generations> <run name> [Fir RUNS dir]"}
POPULATION=${2:?"population"}
GENERATIONS=${3:?"generations"}
NAME=${4:?"run name"}
REMOTE_RUNS=${5:-/scratch/${USER}/sagehen}
[ "${LAYOUT}" = lumped ] || [ "${LAYOUT}" = grid ] || { echo "$0: layout is lumped or grid"; exit 1; }

LOCAL="${HERE}/runs/${NAME}"
REMOTE="${REMOTE_RUNS}/${NAME}"

rm -rf "${LOCAL}"
mkdir -p "${LOCAL}/output" "${LOCAL}/state"
echo "Building the ${LAYOUT} domain in ${LOCAL}..."
"${PYTHON}" "${MFLOW}/tools/build_sagehen9_calibration.py" \
            "${LAYOUT}" "${CALIB}/sagehen/forcing" "${LOCAL}/domain" 2>&1 | grep -v -i deprecat
cp "${LOCAL}/domain/settings/SUMMA/coldState.nc" "${LOCAL}/state/"
if [ "${LAYOUT}" = lumped ]; then
  cp "${CALIB}/sagehen/strt_spunup_lumped.txt" "${LOCAL}/domain/mf6/strt1.txt"
  cp "${CALIB}/sagehen/state_spunup_lumped.nc" "${LOCAL}/state/coldState.nc"
fi

# everything the run reads goes inside the run directory
cp -r "${CALIB}/sagehen/observations" "${LOCAL}/observations"
cp -r "${MFLOW}/domain_sagehen9/settings/mizuroute" "${LOCAL}/mizuroute"

source "${CALIB}/sagehen_config.bash"
write_sagehen_config "${LOCAL}/config.toml" "${REMOTE}" "${REMOTE}/observations" "${REMOTE}/mizuroute" \
                     "$((POPULATION * GENERATIONS))" "${POPULATION}"

# the builder writes this machine's paths into the coupler config
grep -rIl --exclude='*.lst' "${LOCAL}" "${LOCAL}" | while read -r f; do
  sed -i.bak "s|${LOCAL}|${REMOTE}|g" "${f}" && rm "${f}.bak"
done
if grep -rIl --exclude='*.lst' "${WORKFLOW}" "${LOCAL}" >/dev/null; then
  echo "WARNING: these still name a path on this machine:"
  grep -rIl --exclude='*.lst' "${WORKFLOW}" "${LOCAL}"
fi

echo
echo "Ready: ${LOCAL}"
echo "  population ${POPULATION} x ${GENERATIONS} generations; ${POPULATION} workers + 1 = $((POPULATION + 1)) MPI ranks"
echo "Copy it to Fir:"
echo "  rsync -a ${LOCAL}/ fir:${REMOTE}/"
echo "then on Fir:"
echo "  sbatch --ntasks=$((POPULATION + 1)) \$WORKFLOW/hpc/fir/sagehen_calib.fir.slurm ${NAME}"
