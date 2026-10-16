#!/usr/bin/env python3
"""Build ex-gwf-sagehen-sfr: ex-gwf-sagehen with SFR in place of CHD.

One SFR reach per channel cell of the D8 network build_sagehen9.py cuts, so the
stream boundary follows the nine mizuRoute reaches cell by cell. The coupler
returns SFR's aquifer exchange with role = baseflow, as it did CHD's, and SUMMA
keeps the delivery; SFR's own channel flow is never read back.

DIS/TOP is ex-gwf-sagehen's, the soil-column base; the drain stays at land surface.

Streambed top is RBED_DEPTH below land surface, under the drain, and the bed
conductivity is the cell's vertical conductivity. Width is mizuRoute's.

Usage:
    build_sagehen_sfr.py [<output model dir>]      (default ../ex-gwf-sagehen-sfr)
"""

import os
import shutil
import sys

import numpy as np

from build_sagehen9 import (MF6, NCOL, NROW, CELL_AREA, WSCALE, DXY, read_grid, flow_directions,
                            accumulate, build_network)

RBED_DEPTH = 1.0      # m, streambed top below land surface
RBED_THK = 0.5        # m
MANNING = 0.04
MIN_SLOPE = 1e-4
def main(out_dir):
    top, active = read_grid()
    fill, rec = flow_directions(top, active)
    acc, order = accumulate(rec, active)
    links, down, sub, chan = build_network(top, active, rec, acc, order)
    with open(os.path.join(MF6, "kv1.txt")) as f:
        kv = np.array(f.read().split(), dtype=float).reshape(NROW, NCOL) / 86400.0

    cells = [c for cells in links for c in cells]
    rno = {c: n + 1 for n, c in enumerate(cells)}

    # each cell drains to its D8 receiver if that is a channel cell; the outlet cell to nothing
    dn = {c: (rno[tuple(rec[c])] if rec[c][0] >= 0 and tuple(rec[c]) in rno else 0) for c in cells}
    ups = {n: [] for n in rno.values()}
    for c, d in dn.items():
        if d:
            ups[d].append(rno[c])

    def step(c):
        if not dn[c]:
            return DXY
        r = cells[dn[c] - 1]
        return DXY * (np.sqrt(2.0) if abs(r[0] - c[0]) + abs(r[1] - c[1]) == 2 else 1.0)

    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    shutil.copytree(MF6, out_dir, ignore=shutil.ignore_patterns(
        "*.chd", "chd1.txt", "*.lst", "*.hds", "*.cbc", ".DS_Store", "mfsim.*.*"))

    with open(os.path.join(out_dir, "sagehen.sfr"), "w") as f:
        f.write("# Streamflow routing in place of CHD, written by tools/build_sagehen_sfr.py.\n")
        f.write(f"# One reach per D8 channel cell (acc >= threshold of build_sagehen9.py), {len(cells)} in all.\n")
        f.write(f"# Streambed top {RBED_DEPTH} m below land surface, {RBED_THK} m thick, K = the cell's kv (m/s).\n")
        f.write("# The coupler returns the reach-aquifer exchange to SUMMA as baseflow.\n")
        f.write("BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n\n")
        f.write(f"BEGIN DIMENSIONS\n  NREACHES  {len(cells)}\nEND DIMENSIONS\n\n")
        f.write("BEGIN PACKAGEDATA\n# rno lay row col rlen rwid rgrd rtp rbth rhk man ncon ustrf ndv\n")
        for c in cells:
            n = rno[c]
            ncon = len(ups[n]) + (1 if dn[c] else 0)
            width = WSCALE * np.sqrt(acc[c] * CELL_AREA)
            r = cells[dn[c] - 1] if dn[c] else c
            grad = max((top[c] - top[r]) / step(c), MIN_SLOPE)
            f.write(f"  {n:4d} 1 {c[0] + 1:3d} {c[1] + 1:3d} {step(c):8.3f} {width:7.4f} {grad:.6f} "
                    f"{top[c] - RBED_DEPTH:9.3f} {RBED_THK} {max(kv[c], 1e-9):.4e} {MANNING} {ncon} 1.0 0\n")
        f.write("END PACKAGEDATA\n\nBEGIN CONNECTIONDATA\n# rno  upstream (+)  downstream (-)\n")
        for c in cells:
            n = rno[c]
            f.write(f"  {n:4d} " + " ".join(str(u) for u in ups[n]) + (f" -{dn[c]}" if dn[c] else "") + "\n")
        f.write("END CONNECTIONDATA\n")

    nam = os.path.join(out_dir, "sagehen.nam")
    with open(nam) as f:
        text = f.read()
    text = text.replace("  CHD6   sagehen.chd   chd\n", "  SFR6   sagehen.sfr   SFR\n")
    with open(nam, "w") as f:
        f.write(text)

    nout = sum(1 for c in cells if not dn[c])
    print(f"{len(cells)} SFR reaches over {len(links)} links, {nout} outlet reach(es) at "
          f"{[(c[0] + 1, c[1] + 1) for c in cells if not dn[c]]}")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    main(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "ex-gwf-sagehen-sfr")))
