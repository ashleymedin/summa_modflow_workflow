#!/usr/bin/env python3
"""Add a UZF package below SUMMA's soil column to a Sagehen MODFLOW 6 model.

One landflag-1 UZF cell per active cell, except the CHD outlet cells where the model has CHD, as the upstream
ex-gwf-sagehen does. Parameters are the upstream example's, but for 3 trailing waves: VKS the cell's vertical conductivity,
THTS from thts1.txt, THTR 0.1002 and EPS 4. THTI is the water content whose kinematic-wave flux
VKS ((theta - THTR)/(THTS - THTR))^EPS carries the long-term mean recharge, so the unsaturated
zone starts draining rather than filling. No ET, groundwater ET, seepage or mover: the coupler
writes SINF each step and returns REJINF to the soil column, and EVT keeps groundwater ET.

Usage:
    build_sagehen_uzf.py [<model dir> ...]      (default ../ex-gwf-sagehen and ../ex-gwf-sagehen-ss)
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
THTR = 0.1002            # upstream uzf_ts.dat extwc
EPS = 4.0                # Brooks-Corey epsilon
SURFDEP = 0.1            # m
RECH_MEAN = 0.116 / 3.15576e7   # m s-1, the mean recharge spinup_sagehen9_heads.py settles the heads under
NTRAIL, NWAVESETS = 3, 1000   # 3 trailing waves give the upstream 15's flux to 1e-4 at a fifth of the cost; daily input to deep cells keeps hundreds of sets alive


def load(model_dir, name):
    with open(os.path.join(model_dir, name)) as f:
        return np.array(f.read().split(), dtype=float)


def add_uzf(model_dir):
    idomain = np.loadtxt(os.path.join(model_dir, "idomain1.txt"))
    nrow, ncol = idomain.shape
    active = idomain > 0
    vks = load(model_dir, "kv1.txt").reshape(nrow, ncol) / 86400.0
    thts = load(model_dir, "thts1.txt").reshape(nrow, ncol)
    chd = set()
    has_chd = "CHD6" in open(os.path.join(model_dir, "sagehen.nam")).read()
    if has_chd and os.path.exists(os.path.join(model_dir, "chd1.txt")):
        for line in open(os.path.join(model_dir, "chd1.txt")):
            w = line.split()
            if len(w) >= 3 and not line.lstrip().startswith("#"):
                chd.add((int(w[1]) - 1, int(w[2]) - 1))
    cells = [(i, j) for i, j in zip(*np.nonzero(active)) if (i, j) not in chd]

    with open(os.path.join(model_dir, "sagehen.uzf"), "w") as f:
        f.write("# Unsaturated zone below SUMMA's soil column, written by tools/build_sagehen_uzf.py.\n"
                "# The coupler overwrites SINF_PVAR (m s-1) each step with SUMMA's drainage where the water table is\n"
                "# below the cell top, and returns REJINF to the soil column; the period FINF here is only a placeholder.\n"
                "# THTI carries the long-term mean recharge, "
                f"{RECH_MEAN * 3.15576e10:.0f} mm a year, at each cell's VKS.\n")
        f.write(f"BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  NUZFCELLS  {len(cells)}\n"
                f"  NTRAILWAVES  {NTRAIL}\n  NWAVESETS  {NWAVESETS}\nEND DIMENSIONS\n\n"
                "BEGIN PACKAGEDATA\n# iuzno lay row col landflag ivertcon surfdep vks thtr thts thti eps\n")
        for n, (i, j) in enumerate(cells, 1):
            ts = max(thts[i, j], THTR + 0.02)
            ti = THTR + (ts - THTR) * min(RECH_MEAN / vks[i, j], 1.0) ** (1.0 / EPS)
            f.write(f"  {n:4d} 1 {i + 1:3d} {j + 1:3d} 1 0 {SURFDEP} {vks[i, j]:.6e} {THTR} {ts:.4f} {ti:.4f} {EPS}\n")
        f.write("END PACKAGEDATA\n\nBEGIN PERIOD 1\n# iuzno finf pet extdp extwc ha hroot rootact\n")
        for n in range(1, len(cells) + 1):
            f.write(f"  {n:4d} 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")
        f.write("END PERIOD\n")

    nam = os.path.join(model_dir, "sagehen.nam")
    text = open(nam).read()
    if "sagehen.uzf" not in text:
        text = text.replace("  RCH6 ", "  UZF6   sagehen.uzf   UZF\n  RCH6 ")
        open(nam, "w").write(text)
    print(f"{model_dir}: {len(cells)} UZF cells")


if __name__ == "__main__":
    dirs = sys.argv[1:] or [os.path.join(HERE, "..", d) for d in ("ex-gwf-sagehen", "ex-gwf-sagehen-ss")]
    for d in dirs:
        add_uzf(os.path.abspath(d))
