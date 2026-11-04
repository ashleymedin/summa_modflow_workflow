#!/bin/bash
# ---------------------------------------------------------------------------------------
# Coupled SUMMA / MODFLOW 6 / mizuRoute calibration on Sagehen Creek, scoring discharge and
# stream temperature at USGS 10343500 and the water level in a synthetic well.
#
# The domain is the 9 D8 subcatchments of ex-gwf-sagehen, each with a stream HRU for its reach,
# built by test_mflow/tools/build_sagehen9_calibration.py from the bundled basin-mean
# forcing (Oct 2016 - Sep 2018) in sagehen/forcing. NSGA-II scores water year 2018 after a
# water-year-2017 spin-up: discharge KGE and T_reach RMSE (degC) at reach 9, the outlet, and
# head RMSE (m) as a departure from its mean at the synthetic well's cell, read from MODFLOW's OBS6
# output, whose record is a reference run's head (make_sagehen_synthetic_well.py).  MODFLOW's K and Sy are
# searched with the SUMMA parameters, as multipliers on the model's input files.  UZF's VKS and GWE's solid
# conductivity are not: these targets cannot see a VKS below the true one, nor the conductivity at all.
#
#   lumped  one land HRU per GRU (18 HRUs), mapped onto its GRU's cells; the default
#   grid    one land HRU per MODFLOW cell (3396 HRUs), hours per trial
#
# Needs bin/summa_modflow6_opt_sundials_mizuroute.exe, built with
# -DUSE_MPI=ON -DUSE_MIZUROUTE=ON -DUSE_MODFLOW6=ON, and a python3 with netCDF4 and pyproj.
#
# Usage:  ./test_calibration_sagehen.sh [lumped|grid] [population] [generations] [n_ranks]
# ---------------------------------------------------------------------------------------
set -euo pipefail

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
SUMMA_ROOT="${SUMMA_ROOT:-/Users/amedin/Research/MODFLOW/modflow6/srcextern/summa}"
MFLOW="$(cd "${TEST_DIR}/../test_mflow" && pwd)"
DATA="${TEST_DIR}/sagehen"

LAYOUT=${1:-lumped}
POPULATION=${2:-6}
GENERATIONS=${3:-3}
N_RANKS=${4:-7}
N_SAMPLES=$((POPULATION * GENERATIONS))
WORK="${TEST_DIR}/sagehen_${LAYOUT}_run"
PYTHON=${PYTHON:-python3}

SUMMA_EXE="${SUMMA_ROOT}/bin/summa_modflow6_opt_sundials_mizuroute.exe"
if [ ! -x "${SUMMA_EXE}" ]; then
  echo "ERROR: ${SUMMA_EXE} not found"
  echo "       build it with -DUSE_MPI=ON -DUSE_MIZUROUTE=ON -DUSE_MODFLOW6=ON"
  exit 1
fi

echo "======================================================================"
echo "SUMMA coupled calibration test (Sagehen, discharge + stream temperature + well level)"
echo "======================================================================"
echo "  executable:  ${SUMMA_EXE}"
echo "  layout:      ${LAYOUT}"
echo "  work dir:    ${WORK}"
echo "  budget:      ${N_SAMPLES} trials (population ${POPULATION} x ${GENERATIONS} generations)"
echo "  ranks:       ${N_RANKS}"
echo

rm -rf "${WORK}"
mkdir -p "${WORK}/output" "${WORK}/state"
echo "Building the ${LAYOUT} domain..."
"${PYTHON}" "${MFLOW}/tools/build_sagehen9_calibration.py" \
            "${LAYOUT}" "${DATA}/forcing" "${WORK}/domain" 2>&1 | grep -v -i deprecat
cp "${WORK}/domain/settings/SUMMA/coldState.nc" "${WORK}/state/"
# start SUMMA and MODFLOW from WY2017-2018 cycled until it ends where it starts (spinup_sagehen9_cycles.py)
if [ "${LAYOUT}" = "lumped" ]; then
  cp "${DATA}/strt_spunup_lumped.txt" "${WORK}/domain/mf6/strt1.txt"
  cp "${DATA}/state_spunup_lumped.nc" "${WORK}/state/coldState.nc"
fi
echo

source "${TEST_DIR}/sagehen_config.bash"
write_sagehen_config "${WORK}/config.toml" "${WORK}" "${DATA}/observations" \
                     "${MFLOW}/domain_sagehen9/settings/mizuroute" "${N_SAMPLES}" "${POPULATION}"

echo "Running nsga2..."
set +e
mpirun -np "${N_RANKS}" "${SUMMA_EXE}" -c "${WORK}/config.toml" > "${WORK}/run.log" 2>&1
status=$?
set -e
if [ ${status} -ne 0 ]; then
  echo "FAIL: calibration exited ${status}"
  echo "--- last 25 lines of ${WORK}/run.log ---"
  tail -25 "${WORK}/run.log"
  exit 1
fi
grep -E "^NSGA-II generation" "${WORK}/run.log" | sed 's/^/  /' || true
echo

echo "Checking the NSGA-II front..."
"${PYTHON}" - "${WORK}" "${POPULATION}" "${GENERATIONS}" <<'PYEOF'
import os, sys
import numpy as np
import netCDF4 as nc

work, pop, ngen = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])

d = nc.Dataset(os.path.join(work, "output", "sagehen_run0001_calibration.nc"))
sense = [s.strip() for s in d["objective"].sense.split(",")]
obj = np.ma.filled(d["objective"][:].astype(float), np.nan).T   # (target, sample)
failed = np.asarray(d["trial_failed"][:]) == 1

fail = False
def check(ok, what):
    global fail
    print(f"  {'ok  ' if ok else 'FAIL'} {what}")
    fail |= not ok

check(obj.shape == (3, pop * ngen), f"objective is 3 targets x {pop * ngen} trials")
check(sense == ["maximize", "minimize", "minimize"], "discharge KGE is maximized, temperature and well RMSE minimized")
check(not bool(failed.all()), f"{int((~failed).sum())} of {failed.size} trials ran")
ran = obj[:, ~failed]
check(bool(np.all(np.isfinite(ran))), "every objective of a trial that ran is finite")
# a daily stream temperature within a few degC of a mountain creek's, and the RMSE of a run that moves
check(bool(np.all((ran[1] > 0.0) & (ran[1] < 10.0))), "temperature RMSE is between 0 and 10 degC")
check(ran.shape[1] > 1 and bool(np.ptp(ran[1]) > 0.0), "temperature RMSE changes with the parameters")
# the well's record is the model's own head, so a trial misses it by no more than the water table moves
check(bool(np.all((ran[2] >= 0.0) & (ran[2] < 5.0))), "well level RMSE is between 0 and 5 m")

# the MODFLOW multipliers are searched within their bounds, and a worker's model holds its last trial's K
for name, lo, hi in [("K_mult", 0.2, 5.0), ("Sy_mult", 0.5, 1.5)]:
    v = np.asarray(d[name][:])
    check(bool(np.all((v >= lo) & (v <= hi)) and np.ptp(v) > 0.0), f"{name} varies within [{lo}, {hi}]")
kh = np.loadtxt(os.path.join(work, "domain", "mf6", "kh1.txt")).ravel()
rank = np.asarray(d["worker_rank"][:])
t = np.asarray(d["completion_time"][:]).reshape(rank.size, -1).astype(np.int64)   # y m d utc h min s ms
done = ((((((t[:, 0] * 13 + t[:, 1]) * 32 + t[:, 2]) * 24 + t[:, 4]) * 60 + t[:, 5]) * 60 + t[:, 6]) * 1000 + t[:, 7])
for r in np.unique(rank):
    last = np.where(rank == r)[0][np.argmax(done[rank == r])]
    run = np.array(open(os.path.join(work, "output", f"modflow_rank{r:04d}", "kh1.txt")).read().split(), float)
    ratio = run[kh > 0] / kh[kh > 0]
    check(bool(np.allclose(ratio, d["K_mult"][last], rtol=1e-12)), f"rank {r}'s kh1.txt is its last trial's K x {d['K_mult'][last]:.4f}")

front = np.asarray(d["pareto_front"][:]) == 1
print()
print(f"  NSGA-II front ({int(front.sum())} trials):")
for i in np.where(front)[0][np.argsort(-obj[0, front])]:
    print(f"    sample {i + 1:3d}: KGE {obj[0, i]:7.4f}   T_reach RMSE {obj[1, i]:6.3f} degC"
          f"   well RMSE {obj[2, i]:6.3f} m")
print(f"  failed trials (saved in output/failed_trials): {int(failed.sum())}")

sys.exit(1 if fail else 0)
PYEOF

echo
echo "PASS"
