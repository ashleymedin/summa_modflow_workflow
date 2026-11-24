#!/usr/bin/env python3
"""Build ex-gwf-sagehen-gwe-2layer: ex-gwf-sagehen-gwe (SFR and GWE) on GSFLOW's two-layer aquifer.

Bottoms, K, K33, Sy and Ss are ex-gwf-sagehen-2layer's (build_sagehen_2layer.py), one file a layer so a
calibration can scale each: kh1, kh2, kv1, kv2 (m/day), sy1, bot1, bot2, strt1, strt2. The stream, drain, UZF,
base exchange and recharge packages stay on layer 1. UZF's VKS and THTS are GSFLOW's, with THTR = THTS - Sy as
MODFLOW-2005 sets it, and SFR's streambed K and thickness are GSFLOW's on ex-gwf-sagehen-sfr's network; the drain
conductance is layer 1's kh. Layer 2 is confined, so it has no Sy; its GWE porosity is layer 1's Sy.

The starting temperatures are equilibrated as build_sagehen_gwe.py's, the geothermal flux entering layer 2.

Usage:
    build_sagehen_gwe_2layer.py [<output model dir>]      (default ../ex-gwf-sagehen-gwe-2layer)
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

import build_sagehen_gwe as g
from build_sagehen9 import CELL_AREA, HERE, NCOL, NROW
from build_sagehen_2layer import GSFLOW, fill_nearest, read_arrays
from build_sagehen_uzf import add_uzf

GWE1 = os.path.join(HERE, "..", "ex-gwf-sagehen-gwe")
GWF2 = os.path.join(HERE, "..", "ex-gwf-sagehen-2layer")
NLAY = 2
PER_DAY = 1.0 / 86400.0
BED_K, BED_THK = 5.0, 1.0      # GSFLOW SFR streambed, m/day and m


g.BED_THK = BED_THK               # the bed SFE conducts through


def layered(name, suffix=""):
    return "".join(f"    OPEN/CLOSE  {name}{k}.txt{suffix}\n" for k in range(1, NLAY + 1))


def rewrite(path, rows):
    """Replace a package's PERIOD 1 rows."""
    text = open(path).read()
    head, rest = text.split("BEGIN PERIOD 1\n", 1)
    tail = rest[rest.index("END PERIOD"):]
    g.write(path, head + "BEGIN PERIOD 1\n" + "".join(rows) + tail)


def read_uzf(path):
    """VKS (m/day) and THTS from GSFLOW's UZF input; UZFBND is the first array, EPS the third."""
    return [a for i, a in enumerate(read_arrays(path, 4)) if i in (1, 3)]


def equilibrate(model_dir):
    """build_sagehen_gwe.equilibrate on two layers."""
    active = np.loadtxt(os.path.join(model_dir, "idomain1.txt")) > 0
    with tempfile.TemporaryDirectory() as tmp:
        run = os.path.join(tmp, "eq")
        shutil.copytree(model_dir, run)
        g.write(os.path.join(run, "sagehen.tdis"),
                f"BEGIN OPTIONS\n  TIME_UNITS  SECONDS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  NPER  1\nEND DIMENSIONS\n\n"
                f"BEGIN PERIODDATA\n  {g.EQ_YEARS * 365.25 * 86400.0:.6e}  60  1.1\nEND PERIODDATA\n")
        sto = os.path.join(run, "sagehen.sto")
        g.write(sto, open(sto).read().replace("  TRANSIENT\n", "  STEADY-STATE\n"))
        g.add_esl(run, active, g.Q_GEO * CELL_AREA, "Basal geothermal flux (W per 90 m cell).")
        esl = os.path.join(run, "sagehen_gwe.esl")
        rewrite(esl, [f"  {NLAY} {r + 1:3d} {c + 1:3d} {g.Q_GEO * CELL_AREA:.4f}\n" for r, c in zip(*np.nonzero(active))])
        g.write(os.path.join(run, "sagehen_gwe.oc"),
                "BEGIN OPTIONS\n  TEMPERATURE  FILEOUT  sagehen_gwe.ucn\nEND OPTIONS\n\n"
                "BEGIN PERIOD 1\n  SAVE  TEMPERATURE  STEPS 59 60\nEND PERIOD\n")
        p = subprocess.run([g.MF6], cwd=run, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if p.returncode:
            sys.exit(p.stdout[-3000:])
        t = np.array(g.read_arrays(os.path.join(run, "sagehen_gwe.ucn"), NROW, NCOL)).reshape(2, NLAY, NROW, NCOL)
    return np.where(active, t[1], 0.0), np.abs(t[1] - t[0])[:, active].max()


def main(out):
    if os.path.exists(out):
        shutil.rmtree(out)
    shutil.copytree(GWE1, out, ignore=shutil.ignore_patterns(
        "*.lst", "*.hds", "*.cbc", "*.ucn", "*.grb", "*.bin", "*.stdout", ".DS_Store", "kh1.txt", "kv1.txt", "sy1.txt",
        "bot1.txt", "strt_gwe1.txt"))

    arr = {n: np.loadtxt(os.path.join(GWF2, f"{n}.txt")).reshape(NLAY, NROW, NCOL) for n in ("bot", "kh", "kv", "sy")}
    for k in range(NLAY):
        np.savetxt(os.path.join(out, f"bot{k + 1}.txt"), arr["bot"][k], fmt="%10.3f")
        np.savetxt(os.path.join(out, f"kh{k + 1}.txt"), arr["kh"][k], fmt="%12.5e")
        np.savetxt(os.path.join(out, f"kv{k + 1}.txt"), arr["kv"][k], fmt="%12.5e")
    np.savetxt(os.path.join(out, "sy1.txt"), arr["sy"][0], fmt="%8.4f")
    shutil.copy(os.path.join(out, "strt1.txt"), os.path.join(out, "strt2.txt"))
    kh = arr["kh"][0] * PER_DAY
    active = np.loadtxt(os.path.join(out, "idomain1.txt")) > 0
    sy1 = arr["sy"][0]
    vks, thts = read_uzf(os.path.join(GSFLOW, "sagehen.uzf"))
    known = np.loadtxt(os.path.join(GWF2, "ibnd1.txt")).reshape(NROW, NCOL) != 0
    vks, thts = (np.where(active, fill_nearest(a, known & (a > 0), active), 0.0) for a in (vks, thts))
    np.savetxt(os.path.join(out, "vks1.txt"), vks, fmt="%12.5e")
    np.savetxt(os.path.join(out, "thts1.txt"), thts, fmt="%7.4f")
    np.savetxt(os.path.join(out, "thtr1.txt"), np.where(active, np.maximum(thts - sy1, 0.0), 0.0), fmt="%7.4f")

    dis = open(os.path.join(GWF2, "sagehen.dis")).read()
    dis = dis.replace("  BOTM\n    OPEN/CLOSE  bot.txt\n", "  BOTM LAYERED\n" + layered("bot")).replace("(bot.txt, ", "(")
    g.write(os.path.join(out, "sagehen.dis"), dis)
    g.write(os.path.join(out, "sagehen_gwe.dis"),
            dis.replace("  LENGTH_UNITS  METERS\n", "  LENGTH_UNITS  METERS\n  NOGRB\n"))
    f = f"  FACTOR  {PER_DAY:.15g}"
    npf = open(os.path.join(GWF2, "sagehen.npf")).read()
    npf = re.sub(r"  K\n.*\n", "  K LAYERED\n" + layered("kh", f), npf)
    npf = re.sub(r"  K33\n.*\n", "  K33 LAYERED\n" + layered("kv", f), npf)
    g.write(os.path.join(out, "sagehen.npf"), npf.replace("(kh.txt, kv.txt)", "(kh1, kh2, kv1, kv2.txt)"))
    sto = open(os.path.join(GWF2, "sagehen.sto")).read()
    g.write(os.path.join(out, "sagehen.sto"),
            re.sub(r"  SY\n.*\n", "  SY LAYERED\n    OPEN/CLOSE  sy1.txt\n    CONSTANT  0.0\n", sto).replace("(sy.txt,", "(sy1.txt,"))
    g.write(os.path.join(out, "sagehen.ic"),
            f"# Initial conditions: starting heads (m), one file a layer.\nBEGIN GRIDDATA\n  STRT LAYERED\n{layered('strt')}END GRIDDATA\n")

    # the streambed, and the drain conductance from layer 1's kh
    sfr = os.path.join(out, "sagehen.sfr")
    lines, inblock = [], False
    for line in open(sfr):
        if line.startswith("BEGIN PACKAGEDATA"):
            inblock = True
        elif line.startswith("END PACKAGEDATA"):
            inblock = False
        elif inblock and not line.startswith("#"):
            w = line.split()
            w[8], w[9] = f"{BED_THK}", f"{BED_K * PER_DAY:.4e}"
            line = "  " + " ".join(w) + "\n"
        lines.append(line)
    g.write(sfr, "".join(lines))
    drn = os.path.join(out, "sagehen.drn")
    period = open(drn).read().split("BEGIN PERIOD 1\n")[1].split("END PERIOD")[0]
    rows = []
    for line in period.splitlines():
        w = line.split()
        if len(w) >= 5 and not line.lstrip().startswith("#"):
            w[4] = f"{kh[int(w[1]) - 1, int(w[2]) - 1] * CELL_AREA:.6e}"
            rows.append("  " + " ".join(w) + "\n")
    rewrite(drn, rows)
    add_uzf(out)

    est = os.path.join(out, "sagehen_gwe.est")
    g.write(est, open(est).read().replace("  POROSITY\n    OPEN/CLOSE  sy1.txt\n",
                                          "  POROSITY LAYERED\n    OPEN/CLOSE  sy1.txt\n    OPEN/CLOSE  sy1.txt\n"))
    g.write(os.path.join(out, "sagehen_gwe.ic"), f"BEGIN GRIDDATA\n  STRT\n    CONSTANT  {g.STRT}\nEND GRIDDATA\n")
    nam = os.path.join(out, "sagehen_gwe.nam")
    g.write(nam, re.sub(r"  ESL6.*\n", "", open(nam).read()))

    temp, change = equilibrate(out)
    for k in range(NLAY):
        g.write_array(os.path.join(out, f"strt_gwe{k + 1}.txt"), temp[k])
    g.write(os.path.join(out, "sagehen_gwe.ic"),
            "# Thermal equilibrium from tools/build_sagehen_gwe_2layer.py, degC.\n"
            f"BEGIN GRIDDATA\n  STRT LAYERED\n{layered('strt_gwe')}END GRIDDATA\n")
    top = np.loadtxt(os.path.join(out, "land1.txt"))
    trch = np.where(active, g.T_REF + g.GW_OFFSET - g.LAPSE * (top - g.Z_REF), 0.0)
    cells = g.uzf_cells(out)
    g.write_sfe(out, [temp[0][c] for c in g.sfr_cells(out)])
    g.write_uze(out, [temp[0][c] for c in cells], [trch[c] for c in cells])
    g.add_esl(out, active, 0.0, "The coupler overwrites SENERRATE (__INPUT__/SAGEHEN_GWE/ESL/SENERRATE, W) with\n"
              "# SUMMA's conduction out the soil base, at each cell's water-table node.")
    for k in range(NLAY):
        t = temp[k][active]
        print(f"layer {k + 1}: equilibrium {t.min():.2f} to {t.max():.2f} degC, mean {t.mean():.2f}")
    print(f"last step changed it by at most {change:.1e} degC; wrote {out}")


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "ex-gwf-sagehen-gwe-2layer")))
