#!/usr/bin/env python3
"""Add the base exchange's GHB and DRN to a Sagehen MODFLOW 6 model.

One row of each per UZF cell, at DIS/TOP with zero conductance: the coupler sets the GHB's stage and
the DRN's elevation from the last head, and both conductances from SUMMA's base conductance, each
step, so MODFLOW solves the exchange across the soil-column base against its own new head. The
coupler config names them base_ghb_package_name = 'GHBB' and base_drn_package_name = 'DRNB'.

Usage:
    build_sagehen_base.py [<model dir> ...]      (default every ../ex-gwf-sagehen* with its own UZF)
"""

import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGES = (("ghbb", "GHB6", "GHBB", "stage"), ("drnb", "DRN6", "DRNB", "elevation"))


def uzf_cells(model_dir):
    """(row, col) of every UZF cell, in packagedata order."""
    cells, inside = [], False
    with open(os.path.join(model_dir, "sagehen.uzf")) as f:
        for line in f:
            t = line.split()
            if t[:2] == ["BEGIN", "PACKAGEDATA"]:
                inside = True
            elif t[:2] == ["END", "PACKAGEDATA"]:
                break
            elif inside and t and not t[0].startswith("#"):
                cells.append((int(t[2]), int(t[3])))
    return cells


def add_base(model_dir):
    top = np.loadtxt(os.path.join(model_dir, "top1.txt"))
    cells = uzf_cells(model_dir)
    for ext, _, _, what in PACKAGES:
        with open(os.path.join(model_dir, f"sagehen.{ext}"), "w") as f:
            f.write(f"# Base exchange {ext[:3].upper()}, written by tools/build_sagehen_base.py: one row per UZF cell, at\n"
                    f"# DIS/TOP with zero conductance; the coupler sets the {what} and conductance each step.\n"
                    "BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n\n"
                    f"BEGIN DIMENSIONS\n  MAXBOUND  {len(cells)}\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n")
            for r, c in cells:
                f.write(f"  1 {r:3d} {c:3d} {top[r - 1, c - 1]:.4f} 0.0\n")
            f.write("END PERIOD\n")
    add_to_nam(model_dir)
    print(f"{model_dir}: {len(cells)} rows each")


def add_to_nam(model_dir):
    nam = os.path.join(model_dir, "sagehen.nam")
    text = open(nam).read()
    for ext, ftype, pname, _ in PACKAGES:
        if f"sagehen.{ext}" not in text:
            text = text.replace("  OC6 ", f"  {ftype:6s} sagehen.{ext:5s} {pname}\n  OC6 ", 1)
    open(nam, "w").write(text)


if __name__ == "__main__":
    dirs = sys.argv[1:] or sorted(d for d in glob.glob(os.path.join(HERE, "..", "ex-gwf-sagehen*"))
                                  if os.path.isfile(os.path.join(d, "sagehen.uzf"))
                                  and not os.path.islink(os.path.join(d, "sagehen.uzf")))
    for d in dirs:
        add_base(os.path.normpath(d))
