#!/usr/bin/env python3
"""Write the coupler's map_file from HRU polygons and a MODFLOW 6 DIS grid.

Each row is "iHRU cell area": the area (m2) of the HRU polygon that falls in the active cell.
The coupler normalises the weights per HRU and scatters by the area each HRU holds of a cell,
so recharge volume is conserved on a domain whose HRUs are not cells.

The grid is placed by the DIS XORIGIN/YORIGIN (lower-left corner) and must have ANGROT 0 and
constant DELR/DELC. HRU polygons are reprojected to --grid-epsg. iHRU is SUMMA's HRU index,
GRU-major, the row order of attributes.nc.

Usage:
    build_map_file.py --hru HRUS.gpkg --index-field IHRU --dis MODEL/x.dis --grid-epsg 32617 --out map.txt
"""

import argparse
import os
import re
import sys

import geopandas as gpd
import numpy as np
from shapely.geometry import box


def read_dis(path):
    """NROW, NCOL, DELR, DELC, XORIGIN, YORIGIN, IDOMAIN from a DIS file with constant spacing."""
    text = open(path).read()
    body = "\n".join(line.split("#")[0] for line in text.splitlines())

    def num(key, default=None):
        m = re.search(rf"\b{key}\s+([-+0-9.eE]+)", body, re.I)
        if m is None and default is None:
            sys.exit(f"{path}: no {key}")
        return float(m.group(1)) if m else default

    nrow, ncol = int(num("NROW")), int(num("NCOL"))
    if num("ANGROT", 0.0) != 0.0:
        sys.exit("ANGROT must be 0")
    delr = re.search(r"\bDELR\s+CONSTANT\s+([-+0-9.eE]+)", body, re.I)
    delc = re.search(r"\bDELC\s+CONSTANT\s+([-+0-9.eE]+)", body, re.I)
    if not (delr and delc):
        sys.exit("DELR and DELC must be CONSTANT")
    idom = np.ones((nrow, ncol), int)
    m = re.search(r"\bIDOMAIN\s+OPEN/CLOSE\s+(\S+)", body, re.I)
    if m:
        idom = np.loadtxt(os.path.join(os.path.dirname(path), m.group(1))).astype(int).reshape(nrow, ncol)
    return nrow, ncol, float(delr.group(1)), float(delc.group(1)), num("XORIGIN", 0.0), num("YORIGIN", 0.0), idom


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hru", required=True)
    ap.add_argument("--index-field", required=True)
    ap.add_argument("--dis", required=True)
    ap.add_argument("--grid-epsg", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    nrow, ncol, dx, dy, x0, y0, idom = read_dis(a.dis)
    hru = gpd.read_file(a.hru).to_crs(a.grid_epsg)
    hru = hru.dissolve(by=a.index_field).reset_index()

    ytop = y0 + nrow * dy
    rows, cols = np.nonzero(idom > 0)
    cells = gpd.GeoDataFrame(
        {"cell": rows * ncol + cols + 1},
        geometry=[box(x0 + j * dx, ytop - (i + 1) * dy, x0 + (j + 1) * dx, ytop - i * dy) for i, j in zip(rows, cols)],
        crs=a.grid_epsg)

    inter = gpd.overlay(hru[[a.index_field, "geometry"]], cells, how="intersection", keep_geom_type=True)
    inter["area"] = inter.geometry.area
    inter = inter[inter["area"] > 1e-6].sort_values([a.index_field, "cell"])

    with open(a.out, "w") as f:
        f.write(f"# HRU -> MODFLOW cell map, written by build_map_file.py from {os.path.basename(a.hru)}\n")
        f.write("# columns: hru_index  cell_index((i-1)*ncol+j)  intersected area (m2)\n")
        for r in inter.itertuples(index=False):
            f.write(f"{int(getattr(r, a.index_field))} {int(r.cell)} {r.area:.3f}\n")

    poly = hru.set_index(a.index_field).geometry.area
    got = inter.groupby(a.index_field)["area"].sum()
    out = (1 - got.reindex(poly.index).fillna(0) / poly).abs()
    print(f"{len(inter)} rows for {len(hru)} HRUs over {len(cells)} active cells")
    print(f"HRU area outside active cells: worst {100 * out.max():.3f}%, total "
          f"{100 * (1 - got.sum() / poly.sum()):.3f}%")


if __name__ == "__main__":
    main()
