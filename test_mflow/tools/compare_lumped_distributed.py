#!/usr/bin/env python3
"""Compare the lumped domain_sagehen1 with the distributed domain_sagehen9 over the same event.

Both domains run the February 2017 event on the same decisions and parameters, over the same
MODFLOW model, so what differs is model structure. The two carry different areas -- sagehen1
is 5.47% larger than the grid -- so every quantity is per unit area: HRU fluxes are summed
over the run and weighted by HRU area, GRU runoff by GRU area, all in mm.

Usage:
    compare_lumped_distributed.py          (run from anywhere; reads the four run_sagehen*_wet*
                                            and run_sagehen9* outputs under this directory)
"""

import os
import sys

import numpy as np
from netCDF4 import Dataset

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = "domain_sagehen{}/simulations/run_1/SUMMA/{}_timestep.nc"
ATTR = "domain_sagehen{}/settings/SUMMA/attributes.nc"

# (label, domain, output prefix)
PAIRS = {
    "wet": [("lumped", 1, "run1_wet_G1-1"),
            ("distributed", 9, "run1_noLatflow_G1-9"),
            ("distributed + latflow", 9, "run1_latflow_G1-9"),
            ("latflow + SFR", 9, "run1_sfr_G1-9")],
    "deeproot": [("lumped", 1, "run1_wet_deeproot_G1-1"),
                 ("distributed", 9, "run1_deeproot_G1-9")],
}

# HRU fluxes in m s-1 of the HRU, summed over the run
HRU_FLUX = ["scalarRainPlusMelt", "scalarTotalET", "scalarCanopyTranspiration", "scalarInfiltration",
            "scalarSurfaceRunoff", "scalarSoilDrainage", "scalarAquiferBaseflow",
            "scalarAquiferSeepage", "scalarAquiferTranspire"]
GRU_FLUX = ["basin__TotalRunoff", "averageRoutedRunoff"]
HRU_FACTOR = ["scalarAquiferRootFrac", "scalarTranspireLimAqfr"]


def load(domain, prefix):
    path = os.path.join(HERE, OUT.format(domain, prefix))
    if not os.path.exists(path):
        sys.exit(f"missing {path}: run run_sagehen*.sh first")
    with Dataset(os.path.join(HERE, ATTR.format(domain))) as a:
        area = np.array(a["HRUarea"][:], dtype=float)
        hru2gru = np.array(a["hru2gruId"][:])
        gru_id = np.array(a["gruId"][:])
        land = np.array(a["vegTypeIndex"][:]) != 16 if "streamSegId" in a.variables else np.ones(len(area), bool)
        # the cascade's exits: land HRUs draining to a stream HRU or out of the GRU
        hru_id, down = np.array(a["hruId"][:]), np.array(a["downHRUindex"][:])
        exits = land & ~np.isin(down, hru_id[land])
    d = Dataset(path)
    dt = np.diff(np.array(d["time"][:], dtype=float)).min()      # SUMMA writes seconds
    gru_area = np.array([area[hru2gru == g].sum() for g in gru_id])
    return d, dt, area, land, exits, gru_area


def field(d, name):
    """The variable with -9999 as NaN (a coupled return is missing on the first step), the upland domain."""
    v = np.ma.filled(d[name][:], np.nan).astype(float)
    v[v == -9999.0] = np.nan
    return v[..., 0] if "dom" in d[name].dimensions else v


def to_mm(d, name, dt):
    return dt * (1000.0 if d[name].units == "m s-1" else 1.0)   # kg m-2 s-1 is mm s-1


def summarise(domain, prefix):
    d, dt, area, land, exits, gru_area = load(domain, prefix)
    w = area[land] / area[land].sum()
    row = {}
    for name in HRU_FLUX:
        if name in d.variables:
            row[name] = to_mm(d, name, dt) * np.nansum(field(d, name)[:, land], axis=0) @ w
    # scalarSoilBaseflow nets inflow against outflow and carries exfiltration, so export is taken from the layers
    if "mLayerColumnOutflow" in d.variables:
        row["lateral export"] = 1000.0 * dt * np.nansum(field(d, "mLayerColumnOutflow")[:, :, exits]) / area[land].sum()
    for name in GRU_FLUX:
        if name in d.variables:
            row[name] = to_mm(d, name, dt) * np.nansum(field(d, name), axis=0) @ (gru_area / gru_area.sum())
    for name in HRU_FACTOR:
        if name in d.variables:
            x = np.nanmean(field(d, name)[:, land], axis=0)
            row[name] = x @ w
            if land.sum() > 1:
                row[name + " p10-p90"] = (np.nanpercentile(x, 10), np.nanpercentile(x, 90))
    row["area km2"] = area[land].sum() / 1e6
    d.close()
    return row


def main():
    for case, runs in PAIRS.items():
        rows = [(label, summarise(dom, pre)) for label, dom, pre in runs]
        keys = list(dict.fromkeys(k for _, r in rows for k in r))
        print(f"\n{case}: mm over the run unless marked, area-weighted")
        print(f"{'':32s}" + "".join(f"{label:>24s}" for label, _ in rows))
        for k in keys:
            cells = []
            for _, r in rows:
                v = r.get(k)
                cells.append("" if v is None else f"{v[0]:.3f} - {v[1]:.3f}" if isinstance(v, tuple) else f"{v:.3f}")
            print(f"{k:32s}" + "".join(f"{c:>24s}" for c in cells))


if __name__ == "__main__":
    main()
