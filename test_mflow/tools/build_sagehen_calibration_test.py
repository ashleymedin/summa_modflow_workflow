#!/usr/bin/env python3
"""Build the small coupled-calibration test case that SUMMA bundles in utils/test/test_mflow/sagehen_calibration.

The lumped domain of build_sagehen9_calibration.py (9 GRUs, each a land HRU and a stream HRU, on
ex-gwf-sagehen-gwe), cut to the forcing of October 2016 to November 2017, so a calibration spins up over
water year 2017 and scores October-November 2017.  On top of that domain, so that every package the
coupler drives is exercised:

  EVT       an EVTA package taking SUMMA's aquifer transpiration demand (role gw_et), with roots reaching
            below the soil column (rootingDepth 6 m, doubleExp, both scale factors 1)
  GHBB/DRNB the base exchange, named in the coupler config; on a lumped HRU it is a test of the path,
            not a better closure, since its flux answers the mean head of its cells

The heads start from the spun-up lumped heads (strt_spunup_lumped.txt, spinup_sagehen9_heads.py).  Run
output is left out, and the coupler config names its map file relatively; the test points it at its copy.

Usage:
    build_sagehen_calibration_test.py <basin-mean forcing dir> <spun-up lumped heads> <output dir>
"""

import glob
import os
import re
import shutil
import sys
import tempfile

import numpy as np

import build_sagehen9_calibration as b9c

MONTHS = [f"2016{m:02d}" for m in (10, 11, 12)] + [f"2017{m:02d}" for m in range(1, 12)]
ROOTS = {"rootingDepth": 6.0, "rootScaleFactor1": 1.0, "rootScaleFactor2": 1.0}


def main(forcing_src, strt, out_dir):
    with tempfile.TemporaryDirectory() as tmp:
        for f in sorted(os.listdir(forcing_src)):
            if f.endswith(".nc") and any(m in f for m in MONTHS):
                os.symlink(os.path.join(forcing_src, f), os.path.join(tmp, f))
        b9c.main("lumped", tmp, out_dir)

    mf6 = os.path.join(out_dir, "mf6")
    shutil.copy(strt, os.path.join(mf6, "strt1.txt"))
    for pattern in ("*.lst", "*.hds", "*.cbc", "*.ucn", "*.stdout", "*.grb"):
        for f in glob.glob(os.path.join(mf6, pattern)):
            os.remove(f)

    write_evt(mf6)
    write_config(out_dir)
    deepen_roots(os.path.join(out_dir, "settings", "SUMMA"))


def write_evt(mf6):
    """EVTA with SURFACE at the soil-column base (DIS/TOP) and a DEPTH that never limits."""
    shutil.copy(os.path.join(mf6, "top1.txt"), os.path.join(mf6, "evtsurf1.txt"))
    with open(os.path.join(mf6, "sagehen.evta"), "w") as f:
        f.write("""# Groundwater ET - the coupler writes SUMMA's aquifer transpiration demand into RATE every step.
# SUMMA has already limited the demand by the water table, so DEPTH is large: EVT must not limit it again.
BEGIN OPTIONS
  READASARRAYS
  SAVE_FLOWS
END OPTIONS

BEGIN PERIOD 1
  IEVT
    CONSTANT  1
  SURFACE
    OPEN/CLOSE  evtsurf1.txt
  RATE
    CONSTANT  0.0
  DEPTH
    CONSTANT  1000.0
END PERIOD
""")
    nam = os.path.join(mf6, "sagehen.nam")
    text = open(nam).read().replace("  OC6", "  EVT6   sagehen.evta  EVTA\n  OC6", 1)
    open(nam, "w").write(text)


def write_config(out_dir):
    with open(os.path.join(out_dir, "summa_modflow6.config"), "w") as f:
        f.write("""&coupler
  mf6_model_name     = 'SAGEHEN'
  rch_package_name   = 'RCHA'
  base_ghb_package_name = 'GHBB'  ! the base exchange, solved by MODFLOW against its new head
  base_drn_package_name = 'DRNB'
  uzf_hold_days      = 7.0     ! a week of drainage per UZF wave, so deep cells keep few alive
  bnd_package_names  = 'SFR', 'DRN', 'EVTA'
  bnd_package_roles  = 'baseflow', 'surface_discharge', 'gw_et'
  evt_package_name   = 'EVTA'    ! takes SUMMA's aquifer transpiration demand; the supply returns as gw_et
  gwe_model_name     = 'SAGEHEN_GWE'   ! its water-table temperature is SUMMA's scalarAquiferTemp
  esl_package_name   = 'ESL'
  uze_package_name   = 'UZE'
  map_file           = 'hru2cell_map.txt'
  mf6_epsg           = 0
  feedback           = .true.
/
""")


def deepen_roots(settings):
    path = os.path.join(settings, "localParamInfo.txt")
    out = []
    for line in open(path):
        key = line.split("|")[0].strip() if "|" in line else ""
        if key in ROOTS:
            cols = line.split("|")
            line = f"{cols[0]}| {ROOTS[key]:12.4f} |{cols[2]}|{cols[3]}"
        out.append(line)
    open(path, "w").writelines(out)

    path = os.path.join(settings, "modelDecisions.txt")
    text = open(path).read()
    line = "rootProfil                      doubleExp       ! rooting profile (Zeng 2001 double exponential)\n"
    if re.search(r"^rootProfil\s", text, flags=re.M):
        text = re.sub(r"^rootProfil\s.*\n", line, text, flags=re.M)
    else:
        text += line
    open(path, "w").write(text)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*(os.path.abspath(a) for a in sys.argv[1:]))
