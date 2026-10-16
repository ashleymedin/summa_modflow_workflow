#!/usr/bin/env python3
"""Give a single-GRU SUMMA domain one reach: a stream HRU and a one-reach mizuRoute topology.

The stream HRU is appended after the land HRUs.  Its first domain is an empty upland, its second the
reach water column (two lake layers over a copy of the donor HRU's soil column), and any further
domains are unused.  Its parameters and forcing are copied from the donor HRU, its planform area is
length x width.  The land HRUs are unchanged; mizuRoute routes the GRU's runoff into the one reach,
which carries the GRU's id.

Reads the attributes, cold state, trial parameters and forcingFileList.txt from the settings directory
and writes attributes.nc, coldState.nc, trialParams.nc and forcingFileList.txt to <output>/summa_inputs,
the monthly forcing files from --first-month to --last-month to <output>/forcing, and <output>/topology.nc.
The text settings are not copied.

Wolverine Glacier, from the SYMFLUENCE domain, reach from the terminus to USGS 15236900:

    ./add_stream_hru.py domain_Wolverine/settings/SUMMA domain_Wolverine/forcing/SUMMA_input wolverine \\
        --attributes attributes_glac.nc --cold-state coldstate_glac.nc \\
        --first-month 2016-01 --last-month 2019-12 \\
        --length 2000 --slope 0.03 --width 5 --elevation 390 --latitude 60.374 --longitude -148.912
"""

import argparse
import os

import numpy as np
from netCDF4 import Dataset

UPLAND, STREAM = 1, 6              # globalData.f90 domain types
VEG_WATER = 16                     # USGS land class for open water
LAKE_DEPTHS = [0.2, 0.8]           # m, top to bottom
LAKE_TEMP_K = 274.15


def copy_with_hru(src, dst, fill_new, extra_dims=None):
    """Copy a NetCDF file, adding one HRU; fill_new(name, old, new) fills each hru-dimensioned variable's new slot."""
    extra_dims = extra_dims or {}
    with Dataset(src) as s, Dataset(dst, "w", format="NETCDF4") as d:
        d.setncatts({k: s.getncattr(k) for k in s.ncattrs()})
        size = {name: len(dim) + 1 if name == "hru" else extra_dims.get(name, len(dim))
                for name, dim in s.dimensions.items()}
        for name, dim in s.dimensions.items():
            d.createDimension(name, None if dim.isunlimited() else size[name])
        for name, v in s.variables.items():
            fill = v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None
            o = d.createVariable(name, v.dtype, v.dimensions, zlib=True, complevel=4, fill_value=fill)
            o.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"})
            old = np.ma.filled(v[:], fill if fill is not None else 0)
            if "hru" not in v.dimensions and not any(n in extra_dims for n in v.dimensions):
                o[:] = old
                continue
            new = np.zeros(tuple(size[n] for n in v.dimensions), dtype=old.dtype)
            new[tuple(slice(0, n) for n in old.shape)] = old
            fill_new(name, v.dimensions, old, new)
            o[:] = new


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("settings", help="directory of the domain's SUMMA settings")
    p.add_argument("forcing", help="directory of its monthly forcing files, named in forcingFileList.txt")
    p.add_argument("output", help="directory to write the domain to")
    p.add_argument("--attributes", default="attributes.nc")
    p.add_argument("--cold-state", default="coldState.nc")
    p.add_argument("--trial-params", default="trialParams.nc")
    p.add_argument("--first-month", required=True, help="YYYY-MM, first forcing month written")
    p.add_argument("--last-month", required=True, help="YYYY-MM, last forcing month written")
    p.add_argument("--donor-hru", type=int, default=1, help="hruId whose parameters, forcing and soil the stream takes")
    p.add_argument("--length", type=float, required=True, help="reach length, m")
    p.add_argument("--slope", type=float, required=True, help="reach slope, m m-1")
    p.add_argument("--width", type=float, required=True, help="reach width, m")
    p.add_argument("--elevation", type=float, required=True, help="mean reach elevation, m")
    p.add_argument("--latitude", type=float, required=True)
    p.add_argument("--longitude", type=float, required=True)
    a = p.parse_args()

    os.makedirs(os.path.join(a.output, "forcing"), exist_ok=True)
    settings = os.path.join(a.output, "summa_inputs")
    os.makedirs(settings, exist_ok=True)
    area = a.length * a.width

    with Dataset(os.path.join(a.settings, a.attributes)) as s:
        hru_ids = np.array(s["hruId"][:])
        gru_id = int(np.array(s["gruId"][:]).ravel()[0])
        if len(s.dimensions["gru"]) != 1:
            raise SystemExit("only a single-GRU domain is handled")
        land_area = float(np.sum(s["HRUarea"][:]))
        has_seg = "streamSegId" in s.variables
    donor = int(np.flatnonzero(hru_ids == a.donor_hru)[0])
    stream_id = int(hru_ids.max()) + 1

    # attributes: the stream HRU takes the donor's classes, its own geometry, and names its reach
    stream_attr = {"hruId": stream_id, "hru2gruId": gru_id, "downHRUindex": 0, "longitude": a.longitude,
                   "latitude": a.latitude, "elevation": a.elevation, "HRUarea": area, "tan_slope": a.slope,
                   "contourLength": a.length, "vegTypeIndex": VEG_WATER, "streamSegId": gru_id}

    def fill_attr(name, dims, old, new):
        new[-1] = stream_attr.get(name, old[donor])

    copy_with_hru(os.path.join(a.settings, a.attributes), os.path.join(settings, "attributes.nc"), fill_attr)
    with Dataset(os.path.join(settings, "attributes.nc"), "a") as d:
        if not has_seg:
            v = d.createVariable("streamSegId", "i4", ("hru",))
            v.long_name = "mizuRoute reach id of a stream HRU (0 = reach mapped from the GRU id)"
            v[:] = [0] * (len(hru_ids)) + [gru_id]

    # cold state: an empty upland, then the reach column, lake layers over the donor's soil
    with Dataset(os.path.join(a.settings, a.cold_state)) as s:
        n_soil = int(s["nSoil"][0, donor, 0])
        n_toto = len(s.dimensions["midToto"])
    n_lake = len(LAKE_DEPTHS)
    extra = {"midToto": max(n_toto, n_lake + n_soil), "ifcToto": max(n_toto, n_lake + n_soil) + 1}
    zero_in_stream = {"scalarCanopyIce", "scalarCanopyLiq", "scalarSnowDepth", "scalarSWE", "scalarSfcMeltPond",
                      "glacMass4AreaChange", "scalarGlceWE", "scalarAblFrac", "nSnow", "nGlce"}

    def fill_state(name, dims, old, new):
        if "hru" not in dims:
            return
        if dims[1:] != ("hru", "dom"):
            if dims == ("hru",):
                new[-1] = stream_id
            return
        col = old[:, donor, 0]
        new[:, -1, :] = 0
        new[:len(col), -1, 0] = col                               # the empty upland keeps the donor's column
        if name == "DOMarea":
            new[0, -1, 0] = 0.0
        if name == "DOMelev":
            new[0, -1, 2:] = -9999.0
        if name in ("mLayerDepth", "mLayerTemp", "mLayerVolFracLiq", "mLayerVolFracIce"):
            lake = {"mLayerDepth": LAKE_DEPTHS, "mLayerTemp": [LAKE_TEMP_K] * n_lake,
                    "mLayerVolFracLiq": [1.0] * n_lake, "mLayerVolFracIce": [0.0] * n_lake}[name]
            new[:n_lake + n_soil, -1, 1] = np.concatenate([lake, col[:n_soil]])
        elif name == "iLayerHeight":
            h = np.concatenate([[-sum(LAKE_DEPTHS)], -sum(LAKE_DEPTHS) + np.cumsum(LAKE_DEPTHS)])
            new[:n_lake + n_soil + 1, -1, 1] = np.concatenate([h, col[1:n_soil + 1]])
        elif name == "mLayerMatricHead":
            new[:n_soil, -1, 1] = col[:n_soil]
        elif name == "domType":
            new[0, -1, :2] = [UPLAND, STREAM]
        elif name == "nLake":
            new[0, -1, 1] = n_lake
        elif name == "nSoil":
            new[0, -1, 1] = n_soil
        elif name == "DOMarea":
            new[0, -1, 1] = area
        elif name == "DOMelev":
            new[0, -1, 1] = a.elevation
        elif name == "DOMtan_slope":
            new[0, -1, 1] = a.slope
        elif name == "DOMcontourLength":
            new[0, -1, 1] = a.length
        elif name in ("scalarCanairTemp", "scalarCanopyTemp"):
            new[0, -1, 1] = LAKE_TEMP_K
        elif name not in zero_in_stream and name != "DOMaspect":
            new[0, -1, 1] = col[0]

    copy_with_hru(os.path.join(a.settings, a.cold_state), os.path.join(settings, "coldState.nc"),
                  fill_state, extra)

    def fill_donor(name, dims, old, new):
        new[-1] = stream_id if name == "hruId" else old[donor]

    copy_with_hru(os.path.join(a.settings, a.trial_params), os.path.join(settings, "trialParams.nc"), fill_donor)

    # forcing: the months asked for, the stream HRU forced as the donor
    y0, m0 = map(int, a.first_month.split("-"))
    y1, m1 = map(int, a.last_month.split("-"))
    months = [f"{y:04d}{m:02d}" for y in range(y0, y1 + 1) for m in range(1, 13)
              if (y, m) >= (y0, m0) and (y, m) <= (y1, m1)]
    names = [line.strip().strip("'") for line in open(os.path.join(a.settings, "forcingFileList.txt")) if line.strip()]
    chosen = [n for n in names if any(n.endswith(f"_{ym}.nc") for ym in months)]
    if len(chosen) != len(months):
        raise SystemExit(f"{len(chosen)} forcing files found for {len(months)} months")

    def fill_forcing(name, dims, old, new):
        if name == "hruId":
            new[-1] = stream_id
        elif name == "latitude":
            new[-1] = a.latitude
        elif name == "longitude":
            new[-1] = a.longitude
        else:
            new[..., -1] = old[..., donor]

    for n in chosen:
        copy_with_hru(os.path.join(a.forcing, n), os.path.join(a.output, "forcing", n), fill_forcing)
    with open(os.path.join(settings, "forcingFileList.txt"), "w") as f:
        f.writelines(f"{n}\n" for n in chosen)

    # the one reach, carrying the GRU's id, fed by the whole GRU
    with Dataset(os.path.join(a.output, "topology.nc"), "w", format="NETCDF4") as d:
        d.createDimension("hru", 1)
        d.createDimension("seg", 1)
        for name, dim, val, dtype, units in (("hruId", "hru", gru_id, "i8", "-"),
                                             ("area", "hru", land_area + area, "f8", "m2"),
                                             ("hruToSegId", "hru", gru_id, "i8", "-"),
                                             ("segId", "seg", gru_id, "i8", "-"),
                                             ("downSegId", "seg", 0, "i8", "-"),
                                             ("length", "seg", a.length, "f8", "m"),
                                             ("slope", "seg", a.slope, "f8", "-")):
            v = d.createVariable(name, dtype, (dim,))
            v.units = units
            v[:] = [val]
        d.description = "one reach draining the whole GRU, written by add_stream_hru.py"

    print(f"{len(hru_ids)} land HRUs + stream HRU {stream_id} ({area:.0f} m2), {len(chosen)} forcing files")


if __name__ == "__main__":
    main()
