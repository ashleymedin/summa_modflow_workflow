#!/usr/bin/env python3
"""Spin the MODFLOW heads of a calibration Sagehen domain up to a year that ends where it starts.

The native ex-gwf-sagehen heads put the water table near land surface, and one spin-up year does not
drain them, so a calibration scoring the next year sees the aquifer still falling.  Cycling a year does not get
there either: the low-conductivity ridges drain over centuries.  So this solves for the heads the mean
recharge holds up, then settles them with coupled years:

  1. a coupled run of water years 2017-2018, recording each HRU's soil drainage, the recharge the
     coupler hands MODFLOW
  2. MODFLOW 6 alone at steady state under that drainage averaged over the two years
  3. water year 2017 repeated the way the calibration spins up, SUMMA from its cold state and MODFLOW
     from the current heads, each year's final heads the next year's start, until a year changes the
     heads little

The result replaces strt1.txt in a domain built by build_sagehen9_calibration.py, and is the strt1.txt
every ../ex-gwf-sagehen* model carries.

Usage:
    spinup_sagehen9_heads.py <domain dir> <coupler exe> <mf6 exe> <output strt file> [max years]

The coupler is the serial one, e.g. bin/summa_modflow6_sundials_mizuroute.exe.  Everything runs in a
scratch directory beside the output file and the domain is left untouched.
"""

import glob
import os
import shutil
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset

import build_sagehen9 as b9
from make_sagehen_synthetic_well import read_heads

WY2017 = ("2016-10-01 01:00", "2017-10-01 00:00")   # 8760 hourly steps
WY2018 = ("2016-10-01 01:00", "2018-09-30 23:00")   # to the end of the forcing, 17519 steps
TOL_MEDIAN, TOL_P95 = 0.02, 0.10   # m of head change over a year, over the active cells


def write_array(path, a):
    with open(path, "w") as f:
        for row in a:
            f.write(" ".join(f"{v:.6g}" for v in row) + "\n")


def write_oc(mf6, save):
    with open(os.path.join(mf6, "sagehen.oc"), "w") as f:
        f.write(f"BEGIN OPTIONS\n  HEAD FILEOUT sagehen.hds\nEND OPTIONS\n\nBEGIN PERIOD 1\n  SAVE HEAD {save}\nEND PERIOD\n")


def write_file_manager(path, settings, forcing, output, period):
    with open(path, "w") as f:
        f.write(f"""controlVersion       'SUMMA_FILE_MANAGER_V3.0.0'
simStartTime         '{period[0]}'
simEndTime           '{period[1]}'
tmZoneInfo           'utcTime'
outFilePrefix        'spinup'
settingsPath         '{settings}/'
forcingPath          '{forcing}/'
outputPath           '{output}/'
initConditionFile    'coldState.nc'
attributeFile        'attributes.nc'
trialParamFile       'trialParams.nc'
forcingListFile      'forcingFileList.txt'
decisionsFile        'modelDecisions.txt'
outputControlFile    'outputControl.txt'
globalHruParamFile   'localParamInfo.txt'
globalGruParamFile   'basinParamInfo.txt'
vegTableFile         'TBL_VEGPARM.TBL'
soilTableFile        'TBL_SOILPARM.TBL'
generalTableFile     'TBL_GENPARM.TBL'
noahmpTableFile      'TBL_MPTABLE.TBL'
""")


def run_coupled(exe, work, period):
    """Run the coupler over a period from the heads in work/mf6/strt1.txt; the saved heads."""
    for old in glob.glob(os.path.join(work, "out", "*.nc")):
        os.remove(old)
    write_file_manager(os.path.join(work, "fileManager.txt"), os.path.join(work, "settings"),
                       os.path.join(work, "forcing"), os.path.join(work, "out"), period)
    with open(os.path.join(work, "run.log"), "w") as log:
        subprocess.run([exe, os.path.join(work, "fileManager.txt"), os.path.join(work, "summa_modflow6.config")],
                       cwd=os.path.join(work, "mf6"), stdout=log, stderr=subprocess.STDOUT, check=True)
    return read_heads(os.path.join(work, "mf6", "sagehen.hds"))[1]


def mean_recharge(work, active):
    """The mean soil drainage (m s-1) of each HRU, spread over the cells it maps to."""
    with Dataset(glob.glob(os.path.join(work, "out", "*_day.nc"))[0]) as d:
        drainage = np.asarray(d["scalarSoilDrainage_mean"][:]).mean(axis=0).reshape(len(d.dimensions["hru"]), -1)[:, 0]
    recharge = np.zeros(active.shape)
    for line in open(os.path.join(work, "hru2cell_map.txt")):
        if line.strip() and not line.startswith("#"):
            hru, cell = (int(v) for v in line.split()[:2])
            recharge.flat[cell - 1] = max(drainage[hru - 1], 0.0)
    return recharge


def steady_heads(mf6_exe, work, strt, recharge):
    """MODFLOW 6 alone at steady state under a recharge array (m s-1)."""
    ss = os.path.join(work, "mf6_ss")
    shutil.rmtree(ss, ignore_errors=True)
    shutil.copytree(os.path.join(work, "mf6"), ss)
    write_array(os.path.join(ss, "strt1.txt"), strt)
    write_array(os.path.join(ss, "recharge.txt"), recharge)
    with open(os.path.join(ss, "sagehen.tdis"), "w") as f:
        f.write("BEGIN OPTIONS\n  TIME_UNITS SECONDS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  NPER 1\nEND DIMENSIONS\n\n"
                "BEGIN PERIODDATA\n  1.0 1 1.0\nEND PERIODDATA\n")
    with open(os.path.join(ss, "sagehen.sto"), "w") as f:
        f.write("BEGIN GRIDDATA\n  ICONVERT\n    CONSTANT 1\n  SS\n    CONSTANT 1.0e-7\n  SY\n    OPEN/CLOSE sy1.txt\n"
                "END GRIDDATA\n\nBEGIN PERIOD 1\n  STEADY-STATE\nEND PERIOD\n")
    with open(os.path.join(ss, "sagehen.rcha"), "w") as f:
        f.write("BEGIN OPTIONS\n  READASARRAYS\nEND OPTIONS\n\nBEGIN PERIOD 1\n  RECHARGE\n    OPEN/CLOSE recharge.txt\n"
                "END PERIOD\n")
    write_oc(ss, "LAST")
    # the flow model alone: a GWE model's SSM reads a recharge temperature this RCH no longer carries
    nam = os.path.join(ss, "mfsim.nam")
    lines = [ln for ln in open(nam) if "GWE" not in ln.upper() or ln.lstrip().startswith("#")]
    open(nam, "w").writelines(lines)
    with open(os.path.join(work, "mf6_ss.log"), "w") as log:
        subprocess.run([mf6_exe], cwd=ss, stdout=log, stderr=subprocess.STDOUT, check=True)
    return read_heads(os.path.join(ss, "sagehen.hds"))[1][-1]


def main(domain, exe, mf6_exe, output, max_years=20):
    domain, exe, mf6_exe, output = map(os.path.abspath, (domain, exe, mf6_exe, output))
    work = os.path.splitext(output)[0] + "_spinup"
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(os.path.join(work, "out"))
    shutil.copytree(os.path.join(domain, "mf6"), os.path.join(work, "mf6"))
    shutil.copytree(os.path.join(domain, "settings", "SUMMA"), os.path.join(work, "settings"))
    os.symlink(os.path.join(domain, "forcing"), os.path.join(work, "forcing"))
    for name in ("summa_modflow6.config", "hru2cell_map.txt"):
        shutil.copy(os.path.join(domain, name), work)
    with open(os.path.join(work, "settings", "outputControl.txt"), "a") as f:
        f.write("scalarSoilDrainage | 24 | mean\n")

    top, active = b9.read_grid()
    strt_path = os.path.join(work, "mf6", "strt1.txt")
    strt = np.loadtxt(strt_path).reshape(b9.NROW, b9.NCOL)

    # 1-2: the heads the mean recharge of water years 2017-2018 holds up
    write_oc(os.path.join(work, "mf6"), "STEPS 17519")
    run_coupled(exe, work, WY2018)
    recharge = mean_recharge(work, active)
    head = np.where(active, steady_heads(mf6_exe, work, strt, recharge), strt)
    print(f"steady state under {recharge[active].mean() * 3.1536e10:.0f} mm a year of recharge: depth to water "
          f"median {np.median((top - head)[active]):.2f} m", flush=True)
    strt = head
    write_array(strt_path, strt)

    # 3: settle them under the coupled years
    write_oc(os.path.join(work, "mf6"), "STEPS 8760")
    for year in range(1, int(max_years) + 1):
        head = np.where(active, run_coupled(exe, work, WY2017)[-1], strt)
        change = np.abs(head - strt)[active]
        print(f"year {year:2d}: head change median {np.median(change):.3f} m, 95th percentile "
              f"{np.percentile(change, 95):.3f} m, largest {change.max():.3f} m", flush=True)
        strt = head
        write_array(strt_path, strt)
        if np.median(change) < TOL_MEDIAN and np.percentile(change, 95) < TOL_P95:
            break
    else:
        print(f"not settled after {max_years} years; writing the last heads anyway")

    shutil.copy(strt_path, output)
    print(f"{output}: steady state, then {year} years of water year 2017")


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        sys.exit(__doc__)
    main(*sys.argv[1:])
