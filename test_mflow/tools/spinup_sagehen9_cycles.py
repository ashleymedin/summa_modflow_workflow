#!/usr/bin/env python3
"""Cycle a calibration Sagehen domain over water years 2017-2018 until storage ends where it starts.

spinup_sagehen9_heads.py settles the heads with SUMMA cold-started every year. This carries both
models' states: each cycle starts SUMMA from the last cycle's restart and MODFLOW from its head
restart, and the first starts from coldState.nc and the domain's strt1.txt. The low-conductivity
ridges drain slowly, so after cycle EXTRAPOLATE_AT each head jumps to the limit of its decay. It
stops when soil water and the aquifer's volume each change by less than TOL_MM over a cycle, and
writes the final heads as a STRT array and SUMMA's final state as an initial-condition file, which
the calibration tests copy in.
Each cycle's water-year budgets go to cycles.json in the scratch directory.

Usage:
    spinup_sagehen9_cycles.py <domain dir> <coupler exe> <output strt file> <output state file> [max cycles]

The coupler is the serial one, e.g. bin/summa_modflow6_sundials_mizuroute.exe. Everything runs in a
scratch directory beside the output strt file and the domain is left untouched.
"""

import glob
import json
import os
import shutil
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset, num2date

PERIOD = ("2016-10-01 01:00", "2018-09-30 23:00")
TOL_MM = 5.0          # storage change over a cycle, basin mean
EXTRAPOLATE_AT = 3    # cycle after which the heads jump to the limit of their decay
DECAY = 0.8           # the ridges' head change shrinks by about this a cycle
OUTPUT = """scalarRainfall | 24 | mean
scalarSnowfall | 24 | mean
scalarTotalET | 24 | mean
scalarSoilDrainage | 24 | mean
scalarSurfaceRunoff | 24 | mean
scalarSoilBaseflow | 24 | mean
scalarAquiferBaseflow | 24 | mean
scalarAquiferSeepage | 24 | mean
scalarTotalSoilWat | 24 | instant
scalarSWE | 24 | instant
scalarAquiferStorage | 24 | instant
averageRoutedRunoff | 24 | mean
"""


def write_file_manager(work, init):
    with open(os.path.join(work, "fileManager.txt"), "w") as f:
        f.write(f"""controlVersion       'SUMMA_FILE_MANAGER_V3.0.0'
simStartTime         '{PERIOD[0]}'
simEndTime           '{PERIOD[1]}'
tmZoneInfo           'utcTime'
outFilePrefix        'cycle'
settingsPath         '{work}/settings/'
forcingPath          '{work}/forcing/'
outputPath           '{work}/out/'
initConditionFile    '{init}'
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


def budget(day_file, area):
    """Basin totals (mm) per water year, and soil, snow and aquifer storage (mm) at its start and end."""
    d = Dataset(day_file)
    t = num2date(d["time"][:], d["time"].units, only_use_cftime_datetimes=False)
    wy = np.array([x.year + (1 if x.month >= 10 else 0) for x in t])

    def hru(v):  # (time, hru), the upland and stream domains summed, fill values dropped
        x = np.asarray(d[v][:], float)
        return np.where((x < -9000) | (np.abs(x) > 1e3), 0.0, x).sum(2)

    def mean(x):
        return (x * area).sum(-1) / area.sum()

    gru_area = area.reshape(-1, 2).sum(1)
    out = {}
    for y in sorted(set(wy)):
        m = wy == y
        i0, i1 = np.where(m)[0][[0, -1]]
        r = {k: float(mean(hru(v + "_mean")[m]).sum() * 86400 * s) for k, (v, s) in {
            "p": ("scalarRainfall", 1), "et": ("scalarTotalET", -1)}.items()}
        r["p"] += float(mean(hru("scalarSnowfall_mean")[m]).sum() * 86400)
        for k, v in (("drain", "scalarSoilDrainage"), ("surface", "scalarSurfaceRunoff"),
                     ("latflow", "scalarSoilBaseflow"), ("baseflow", "scalarAquiferBaseflow"),
                     ("seepage", "scalarAquiferSeepage")):
            x = hru(v + "_mean")
            x = np.where(np.abs(x) > 1e-2, 0.0, x)  # a first-step value before the coupler has run
            r[k] = float(mean(x[m]).sum() * 86400 * 1000)
        q = np.asarray(d["averageRoutedRunoff_mean"][:], float)
        r["q"] = float(((q * gru_area).sum(1) / gru_area.sum())[m].sum() * 86400 * 1000)
        for k, v, s in (("soil", "scalarTotalSoilWat", 1), ("swe", "scalarSWE", 1), ("aquifer", "scalarAquiferStorage", 1000)):
            x = mean(hru(v)) * s
            r[k] = (float(x[i0]), float(x[i1]))
        out[f"WY{y}"] = r
    return out


def read_restart_heads(path):
    """The raw bytes of a coupler head restart, its node count, and its heads in reduced node order."""
    raw = open(path, "rb").read()
    assert raw[:8] == b"SUMMF6HD", f"{path} is not a coupler head restart"
    n = int(np.frombuffer(raw[20:24], "<i4")[0])
    return raw, n, np.frombuffer(raw[24:24 + 8 * n], "<f8").copy()


def aquifer_change(prev, last, sy, top):
    """Water gained by the aquifer below the column base between two head restarts, basin mean mm."""
    h0, h1 = read_restart_heads(prev)[2], read_restart_heads(last)[2]
    return 1000 * float(np.mean(sy * (np.minimum(h1, top) - np.minimum(h0, top))))


def extrapolate(prev, last, r):
    """Move each head to the limit of a decay by r a cycle, h + (h - h_prev) r/(1 - r), in place."""
    raw, n, h1 = read_restart_heads(last)
    h0 = read_restart_heads(prev)[2]
    h = h1 + (h1 - h0) * r / (1 - r)
    open(last, "wb").write(raw[:24] + h.astype("<f8").tobytes() + raw[24 + 8 * n:])


def head_strt(restart, idomain, strt):
    """A coupler head restart (reduced node order) as a full STRT array; inactive cells keep strt."""
    head = strt.copy()
    head[idomain > 0] = read_restart_heads(restart)[2]
    return head


def main(domain, exe, out_strt, out_state, max_cycles=8):
    domain, exe, out_strt, out_state = map(os.path.abspath, (domain, exe, out_strt, out_state))
    work = os.path.splitext(out_strt)[0] + "_cycles"
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(os.path.join(work, "out"))
    shutil.copytree(os.path.join(domain, "mf6"), os.path.join(work, "mf6"))
    shutil.copytree(os.path.join(domain, "settings", "SUMMA"), os.path.join(work, "settings"))
    os.symlink(os.path.join(domain, "forcing"), os.path.join(work, "forcing"))
    shutil.copy(os.path.join(domain, "hru2cell_map.txt"), work)
    config = [f"  map_file           = '{work}/hru2cell_map.txt'" if ln.strip().startswith("map_file") else ln
              for ln in open(os.path.join(domain, "summa_modflow6.config")).read().splitlines()]
    open(os.path.join(work, "summa_modflow6.config"), "w").write("\n".join(config) + "\n")
    open(os.path.join(work, "settings", "outputControl.txt"), "w").write(OUTPUT)
    area = np.asarray(Dataset(os.path.join(work, "settings", "attributes.nc"))["HRUarea"][:], float)
    mf6 = os.path.join(work, "mf6")
    idomain = np.loadtxt(os.path.join(mf6, "idomain1.txt"))
    active = idomain.reshape(-1) > 0
    sy = np.loadtxt(os.path.join(mf6, "sy1.txt")).reshape(-1)[active]
    top = np.loadtxt(os.path.join(mf6, "top1.txt")).reshape(-1)[active]

    results = []
    for c in range(int(max_cycles)):
        for old in glob.glob(os.path.join(work, "out", "*")):
            os.remove(old)
        write_file_manager(work, "coldState.nc" if c == 0 else f"restart_{c - 1}.nc")
        cmd = [exe, os.path.join(work, "fileManager.txt"), os.path.join(work, "summa_modflow6.config"),
               "--head-restart-write", os.path.join(work, f"heads_{c}.bin")]
        if c > 0:
            cmd += ["--head-restart-read", os.path.join(work, f"heads_{c - 1}.bin")]
        with open(os.path.join(work, f"run_{c}.log"), "w") as log:
            subprocess.run(cmd, cwd=os.path.join(work, "mf6"), stdout=log, stderr=subprocess.STDOUT, check=True)
        shutil.copy(glob.glob(os.path.join(work, "out", "*restart*.nc"))[0],
                    os.path.join(work, "settings", f"restart_{c}.nc"))
        b = budget(glob.glob(os.path.join(work, "out", "*_day.nc"))[0], area)
        results.append(b)
        json.dump(results, open(os.path.join(work, "cycles.json"), "w"), indent=1)
        change = {k: b["WY2018"][k][1] - b["WY2017"][k][0] for k in ("soil", "swe")}
        if c > 0:  # MODFLOW's own volume: SUMMA's aquifer storage is a lumped HRU's mean head
            change["aquifer"] = aquifer_change(os.path.join(work, f"heads_{c - 1}.bin"),
                                               os.path.join(work, f"heads_{c}.bin"), sy, top)
        w = b["WY2018"]
        print(f"cycle {c}: storage change, mm: " + ", ".join(f"{k} {v:+.1f}" for k, v in change.items())
              + f"; WY2018 P {w['p']:.0f}, ET {w['et']:.0f}, Q {w['q']:.0f}", flush=True)
        if c > 0 and abs(change["soil"]) < TOL_MM and abs(change["aquifer"]) < TOL_MM:
            break
        if c == EXTRAPOLATE_AT:
            extrapolate(os.path.join(work, f"heads_{c - 1}.bin"), os.path.join(work, f"heads_{c}.bin"), DECAY)
    else:
        print(f"not settled after {max_cycles} cycles; writing the last state anyway")

    strt = np.loadtxt(os.path.join(mf6, "strt1.txt")).reshape(idomain.shape)
    head = head_strt(os.path.join(work, f"heads_{c}.bin"), idomain, strt)
    with open(out_strt, "w") as f:
        for row in head:
            f.write(" ".join(f"{v:.6g}" for v in row) + "\n")
    shutil.copy(os.path.join(work, "settings", f"restart_{c}.nc"), out_state)
    print(f"{out_strt}, {out_state}: {c + 1} cycles of water years 2017-2018")


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        sys.exit(__doc__)
    main(*sys.argv[1:])
