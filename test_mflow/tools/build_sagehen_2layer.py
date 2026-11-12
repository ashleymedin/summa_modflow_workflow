#!/usr/bin/env python3
"""Build ex-gwf-sagehen-2layer: ex-gwf-sagehen with GSFLOW's two-layer aquifer.

The grid, IDOMAIN, TOP (the soil-column base) and every boundary package are ex-gwf-sagehen's, so the
map files and coupler configs carry over. From the GSFLOW Sagehen MODFLOW-NWT input it takes the two
layer bottoms (DIS) and K, K33, Sy and Ss (UPW), layer 1 convertible and layer 2 confined as there.

A cell active here with no positive GSFLOW value takes the nearest one. Each
layer keeps at least MIN_THICK under the one above. Both layers start from strt1.txt.

Usage:
    build_sagehen_2layer.py [<GSFLOW modflow input dir> [<output model dir>]]
        (defaults ../../../sagehen-orig/input/modflow and ../ex-gwf-sagehen-2layer)
"""

import os
import shutil
import sys

import numpy as np

from build_sagehen9 import HERE, MF6, NCOL, NROW

GSFLOW = os.path.join(HERE, "..", "..", "..", "sagehen-orig", "input", "modflow")
MIN_THICK = 1.0       # m
PER_DAY = 1.0 / 86400.0


def read_arrays(path, n):
    """The NROW x NCOL arrays of a MODFLOW-2005 free-format file, in order, multipliers applied."""
    lines = [l.split() for l in open(path).read().replace("\r", "").splitlines()]
    arrays, i = [], 0
    while i < len(lines):
        t = lines[i]
        i += 1
        if not t or t[0].startswith("#"):
            continue
        if t[0].upper() == "CONSTANT":
            arrays.append(np.full((NROW, NCOL), float(t[1])))
        elif t[0].upper() == "INTERNAL":
            vals = []
            while len(vals) < NROW * NCOL:
                vals += [float(v) for v in lines[i][:NCOL]]     # BAS numbers each row in a last column
                i += 1
            arrays.append(float(t[1]) * np.array(vals[:NROW * NCOL]).reshape(NROW, NCOL))
        if len(arrays) == n:
            break
    return arrays


def fill_nearest(a, known, need):
    """a where known, else the value of the nearest known cell, over the need mask."""
    out = a.copy()
    kr, kc = np.nonzero(known)
    for r, c in zip(*np.nonzero(need & ~known)):
        j = np.argmin((kr - r) ** 2 + (kc - c) ** 2)
        out[r, c] = a[kr[j], kc[j]]
    return out


def write_layers(path, arrays, fmt):
    with open(path, "w") as f:
        for a in arrays:
            np.savetxt(f, a, fmt=fmt)


def main(gsflow, out_dir):
    bot1, bot2 = read_arrays(os.path.join(gsflow, "sagehen.dis"), 5)[3:]     # DELR, DELC, TOP first
    hk1, vk1, ss1, sy1, hk2, vk2, ss2 = read_arrays(os.path.join(gsflow, "sagehen.upw"), 7)
    gs_active = read_arrays(os.path.join(gsflow, "sagehen.bas"), 1)[0] != 0

    top = np.loadtxt(os.path.join(MF6, "top1.txt"))
    idomain = np.loadtxt(os.path.join(MF6, "idomain1.txt")).astype(int)
    active = idomain > 0
    fill = lambda a: np.where(active, fill_nearest(a, gs_active & (a > 0), active), 0.0)
    hk1, vk1, sy1, hk2, vk2 = map(fill, (hk1, vk1, sy1, hk2, vk2))
    bot1 = np.where(active, np.minimum(fill(bot1), top - MIN_THICK), 0.0)
    bot2 = np.where(active, np.minimum(fill(bot2), bot1 - MIN_THICK), 0.0)
    print(f"{(active & ~gs_active).sum()} cells active here only, "
          f"{(active & (top - bot1 < 10.0)).sum()} with layer 1 under 10 m")

    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    shutil.copytree(MF6, out_dir, ignore=shutil.ignore_patterns(
        "*.lst", "*.hds", "*.cbc", "*.grb", "*.stdout", ".DS_Store", "bot1.txt", "kh1.txt", "kv1.txt",
        "sy1.txt"))

    write_layers(os.path.join(out_dir, "bot.txt"), [bot1, bot2], "%10.3f")
    write_layers(os.path.join(out_dir, "kh.txt"), [hk1, hk2], "%12.5e")
    write_layers(os.path.join(out_dir, "kv.txt"), [vk1, vk2], "%12.5e")
    write_layers(os.path.join(out_dir, "sy.txt"), [sy1, np.zeros_like(sy1)], "%8.4f")

    with open(os.path.join(out_dir, "sagehen.dis"), "w") as f:
        f.write(f"""# Structured discretisation: 2 layers, {NROW} rows, {NCOL} cols, 90 m cells.
# TOP and IDOMAIN are ex-gwf-sagehen's; BOTM is GSFLOW's two layer bottoms (bot.txt, build_sagehen_2layer.py).
BEGIN OPTIONS
  LENGTH_UNITS  METERS
END OPTIONS

BEGIN DIMENSIONS
  NLAY  2
  NROW  {NROW}
  NCOL  {NCOL}
END DIMENSIONS

BEGIN GRIDDATA
  DELR
    CONSTANT  90.0
  DELC
    CONSTANT  90.0
  TOP
    OPEN/CLOSE  top1.txt
  BOTM
    OPEN/CLOSE  bot.txt
  IDOMAIN LAYERED
    OPEN/CLOSE  idomain1.txt
    OPEN/CLOSE  idomain1.txt
END GRIDDATA
""")
    with open(os.path.join(out_dir, "sagehen.npf"), "w") as f:
        f.write(f"""# Node property flow: GSFLOW's UPW conductivities in m/day (kh.txt, kv.txt), converted to m/s.
# Layer 1 convertible, layer 2 confined, as GSFLOW's LAYTYP.
BEGIN OPTIONS
  SAVE_FLOWS
END OPTIONS

BEGIN GRIDDATA
  ICELLTYPE LAYERED
    CONSTANT  1
    CONSTANT  0
  K
    OPEN/CLOSE  kh.txt  FACTOR  {PER_DAY:.15g}
  K33
    OPEN/CLOSE  kv.txt  FACTOR  {PER_DAY:.15g}
END GRIDDATA
""")
    with open(os.path.join(out_dir, "sagehen.sto"), "w") as f:
        f.write(f"""# Storage: GSFLOW's UPW specific yield (sy.txt, layer 1 only) and specific storage.
BEGIN OPTIONS
  SAVE_FLOWS
END OPTIONS

BEGIN GRIDDATA
  ICONVERT LAYERED
    CONSTANT  1
    CONSTANT  0
  SS LAYERED
    CONSTANT  {ss1[0, 0]:.6g}
    CONSTANT  {ss2[0, 0]:.6g}
  SY
    OPEN/CLOSE  sy.txt
END GRIDDATA

BEGIN PERIOD 1
  TRANSIENT
END PERIOD
""")
    with open(os.path.join(out_dir, "sagehen.ic"), "w") as f:
        f.write("""# Initial conditions: both layers start from strt1.txt (m).
BEGIN OPTIONS
END OPTIONS

BEGIN GRIDDATA
  STRT LAYERED
    OPEN/CLOSE  strt1.txt
    OPEN/CLOSE  strt1.txt
END GRIDDATA
""")
    print(f"wrote {out_dir}")


if __name__ == "__main__":
    gsflow = sys.argv[1] if len(sys.argv) > 1 else GSFLOW
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "ex-gwf-sagehen-2layer")
    main(os.path.abspath(gsflow), os.path.abspath(out))
