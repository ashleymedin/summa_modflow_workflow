#!/usr/bin/env python3
"""Write the Sagehen vegetation and soil class of every MODFLOW cell, landcover/veg1.txt and soil1.txt.

Vegetation is the MODIFIED_IGBP_MODIS_NOAH class. NLCD 2021 land cover (30 m) gives the type, the mode of
the cell's 3 x 3 pixels. A forest cell takes its class from NLCD 2017 tree canopy cover, the cell mean, by
IGBP's own thresholds: over 60% (or unmapped) its forest type, 30-60% woody savanna, below that savanna.

Soil is the ROSETTA class of the USDA texture of POLARIS sand and clay, the cell mean averaged over 0-2 m.

Inactive cells are 0. The rasters are not in the repository; the arrays are.

Usage:
    build_sagehen_landcover.py [<data dir>]      (default ../../../sagehen_data: nlcd2021.tif, tcc2017.tif,
                                                  polaris/{sand,clay}_<top>_<bottom>.tif)
"""

import os
import sys

import numpy as np

from build_sagehen9 import DXY, EPSG, HERE, MF6, NCOL, NROW, XORIGIN, YTOP

DATA = os.path.join(HERE, "..", "..", "..", "sagehen_data")
OUT = os.path.join(HERE, "..", "landcover")
DEPTHS = [(0, 5), (5, 15), (15, 30), (30, 60), (60, 100), (100, 200)]   # cm, POLARIS layers

# NLCD class -> IGBP class; forest is split by canopy cover below
NLCD_IGBP = {11: 11, 21: 13, 22: 13, 23: 13, 24: 13, 31: 16, 41: 4, 42: 1, 43: 5, 52: 7, 71: 10,
             81: 12, 82: 12, 90: 11, 95: 11}
FOREST = (41, 42, 43)
WOODY_SAVANNA, SAVANNA = 8, 9

ROSETTA = {"clay": 1, "clay loam": 2, "loam": 3, "loamy sand": 4, "sand": 5, "sandy clay": 6,
           "sandy clay loam": 7, "sandy loam": 8, "silt": 9, "silty clay": 10, "silty clay loam": 11,
           "silt loam": 12}


def usda_texture(sand, clay):
    """USDA texture class of a soil from its sand and clay percentages."""
    silt = 100.0 - sand - clay
    if silt + 1.5 * clay < 15:
        return "sand"
    if silt + 2 * clay < 30:
        return "loamy sand"
    if clay >= 40:
        return "silty clay" if silt >= 40 else ("sandy clay" if sand > 45 else "clay")
    if clay >= 35 and sand > 45:
        return "sandy clay"
    if clay >= 27:
        return "silty clay loam" if sand <= 20 else ("clay loam" if sand <= 45 else "sandy clay loam")
    if clay >= 20 and sand > 45 and silt < 28:
        return "sandy clay loam"
    if silt >= 80 and clay < 12:
        return "silt"
    if silt >= 50:
        return "silt loam"
    if clay >= 7 and silt >= 28 and sand <= 52:
        return "loam"
    return "sandy loam"


def sample(path, active, reduce):
    """reduce() of the raster's pixels at a 3 x 3 grid of 30 m points in each active cell."""
    import rasterio
    from pyproj import Transformer
    out = np.zeros((NROW, NCOL))
    with rasterio.open(path) as r:
        a = r.read(1)
        nodata = r.nodata
        to = Transformer.from_crs(f"EPSG:{EPSG}", r.crs, always_xy=True)
        off = np.array([-DXY / 3, 0.0, DXY / 3])
        for i, j in zip(*np.nonzero(active)):
            x, y = np.meshgrid(XORIGIN + (j + 0.5) * DXY + off, YTOP - (i + 0.5) * DXY + off)
            xs, ys = to.transform(x.ravel(), y.ravel())
            rows, cols = rasterio.transform.rowcol(r.transform, xs, ys)
            v = a[np.array(rows), np.array(cols)]
            v = v[v != nodata] if nodata is not None else v
            out[i, j] = reduce(v)
    return out


def mode(v):
    u, c = np.unique(v, return_counts=True)
    return u[np.argmax(c)]


def main(data):
    active = np.loadtxt(os.path.join(MF6, "idomain1.txt")) > 0

    nlcd = sample(os.path.join(data, "nlcd2021.tif"), active, mode).astype(int)
    canopy = sample(os.path.join(data, "tcc2017.tif"), active, lambda v: v[v <= 100].mean() if (v <= 100).any() else np.nan)
    veg = np.zeros((NROW, NCOL), dtype=int)
    for i, j in zip(*np.nonzero(active)):
        k = nlcd[i, j]
        veg[i, j] = NLCD_IGBP[k] if k not in FOREST or not canopy[i, j] <= 60 else \
            (WOODY_SAVANNA if canopy[i, j] > 30 else SAVANNA)

    thick = np.array([b - t for t, b in DEPTHS], dtype=float)
    sand = sum(w * sample(os.path.join(data, "polaris", f"sand_{t}_{b}.tif"), active, np.mean)
               for w, (t, b) in zip(thick, DEPTHS)) / thick.sum()
    clay = sum(w * sample(os.path.join(data, "polaris", f"clay_{t}_{b}.tif"), active, np.mean)
               for w, (t, b) in zip(thick, DEPTHS)) / thick.sum()
    soil = np.zeros((NROW, NCOL), dtype=int)
    for i, j in zip(*np.nonzero(active)):
        soil[i, j] = ROSETTA[usda_texture(sand[i, j], clay[i, j])]

    os.makedirs(OUT, exist_ok=True)
    np.savetxt(os.path.join(OUT, "veg1.txt"), veg, fmt="%3d",
               header="MODIFIED_IGBP_MODIS_NOAH class per cell, build_sagehen_landcover.py")
    np.savetxt(os.path.join(OUT, "soil1.txt"), soil, fmt="%3d",
               header="ROSETTA soil class per cell, build_sagehen_landcover.py")
    for name, a in (("veg", veg), ("soil", soil)):
        u, c = np.unique(a[active], return_counts=True)
        print(f"{name}: " + ", ".join(f"{k} x {n}" for k, n in zip(u, c)))
    print(f"sand {np.percentile(sand[active], [5, 50, 95]).round(0)} %, "
          f"clay {np.percentile(clay[active], [5, 50, 95]).round(0)} % (5th, 50th, 95th)")


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else DATA))
