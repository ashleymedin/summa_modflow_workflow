#!/usr/bin/env python3
"""Show that a well's drawdown reaches SUMMA's lower boundary, from run_sagehen9_wel.sh against run_sagehen9_noLatflow.sh.

The well is in ex-gwf-sagehen-wel/sagehen.wel; the identity map gives its cell one HRU. SUMMA's
lowerBoundHead there is scalarAquiferStorage / Sy, the head above the soil-column base, and should be
MODFLOW's head less DIS/TOP one step earlier, since the exchange is explicit. Both runs must start from
the same spin-up.

Usage:
    check_sagehen_wel.py          (run from anywhere)
"""

import os
import sys

import numpy as np
from netCDF4 import Dataset

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
MODEL = os.path.join(HERE, "ex-gwf-sagehen")
OUT = os.path.join(HERE, "domain_sagehen9/simulations/run_1/SUMMA/run1_{}_G1-9_timestep.nc")
NROW, NCOL = 73, 81
CELL_AREA = 90.0 * 90.0   # m2


def arr(name):
    return np.loadtxt(os.path.join(MODEL, name)).reshape(NROW, NCOL)


def well():
    with open(os.path.join(HERE, "ex-gwf-sagehen-wel/sagehen.wel")) as f:
        lines = [l.split() for l in f if l.strip() and not l.lstrip().startswith("#")]
    lay, row, col, q = lines[lines.index(["BEGIN", "PERIOD", "1"]) + 1]
    return int(row) - 1, int(col) - 1, float(q)


def heads(path):
    """Every saved step of a one-layer DIS head file, (nstep, nrow, ncol)."""
    rec = 52 + 8 * NROW * NCOL
    b = open(path, "rb").read()
    return np.array([np.frombuffer(b[k * rec + 52:(k + 1) * rec], float).reshape(NROW, NCOL)
                     for k in range(len(b) // rec)])


def hru_of(row, col):
    cell = row * NCOL + col + 1
    for line in open(os.path.join(HERE, "domain_sagehen9/hru2cell_map.txt")):
        p = line.split()
        if p and not p[0].startswith("#") and int(p[1]) == cell:
            return int(p[0]) - 1
    sys.exit(f"cell {cell} is not mapped")


def main():
    row, col, q = well()
    i = hru_of(row, col)
    top, sy = arr("top1.txt")[row, col], arr("sy1.txt")[row, col]
    h = heads(os.path.join(HERE, "ex-gwf-sagehen-wel/sagehen.hds"))[:, row, col] - top

    runs = {}
    for name in ("wel", "noLatflow"):
        path = OUT.format(name)
        if not os.path.exists(path):
            sys.exit(f"missing {path}: run run_sagehen9_{name}.sh first")
        with Dataset(path) as d:
            runs[name] = {v: np.array(d[v][:, i, 0], dtype=float)
                          for v in ("scalarSoilDrainage", "scalarAquiferStorage") if v in d.variables}
    lbh = runs["wel"]["scalarAquiferStorage"] / sy
    q_wel = runs["wel"]["scalarSoilDrainage"] * 3.6e6
    q_ref = runs["noLatflow"]["scalarSoilDrainage"] * 3.6e6
    n = min(len(h), len(lbh))

    print(f"well at row {row + 1}, column {col + 1} (HRU {i + 1}), {-q * 86400:.0f} m3 a day")
    print(f"{'step':>4} {'MODFLOW h - top, m':>19} {'SUMMA lowerBoundHead, m':>24} {'drainage wel / no well, mm/h':>30}")
    for k in range(0, n, 6):
        prev = f"{h[k - 1]:19.3f}" if k else f"{'(restart)':>19}"
        print(f"{k + 1:4d} {prev} {lbh[k]:24.3f} {q_wel[k]:15.4f} {q_ref[k]:14.4f}")
    lag = np.abs(lbh[1:n] - h[:n - 1]).max()
    print(f"max |lowerBoundHead(k) - (h - top)(k-1)| over steps 2-{n}: {lag:.2e} m")
    print(f"drawdown below the base at the end: {-h[n - 1]:.3f} m; extra drainage at the HRU: "
          f"{(q_wel - q_ref).sum() * 1e-3 * CELL_AREA:.1f} m3")


if __name__ == "__main__":
    main()
