#!/usr/bin/env python3
"""Split tributary GRUs off a one-reach SUMMA domain and give every GRU a reach of a mizuRoute network.

Starts from a domain written by add_stream_hru.py: one GRU of land HRUs and one stream HRU, routed to
one reach.  A JSON spec names the reaches (length m, slope, mean elevation m, latitude, longitude, the
reach downstream, 0 at the outlet) and the new GRUs (id, the reach it drains to, area m2, mean elevation
m, latitude, longitude).  The original GRU drains to the reach with its own id.

Each new GRU is one upland HRU, a copy of the donor land HRU's upland domain, plus a stream HRU, a copy
of the original stream HRU; the donor's upland loses the area the new GRUs take, so the land area is
unchanged.  Forcing air temperature and cold-state temperatures shift from the donor's elevation to
each HRU's own at the spec's lapse rate (K m-1), every stream HRU's included.  A reach is wscale x
sqrt(upstream area) wide, as mizuRoute makes it, and its stream HRU is length x width.

Wolverine Creek above USGS 15236900, from NHDPlus HR: the glacier GRU's reach from the east-lobe junction
to the gauge-map tributary, a 2.27 km2 tundra branch entering at its head, the 0.50 km2 tributary, and
the last 454 m to the gauge:

    ./add_stream_network.py tests/wolverine wolverine_network wolverine_network.json

writes <output>/summa_inputs/{attributes,coldState,trialParams}.nc and forcingFileList.txt, the forcing
files to <output>/forcing, and <output>/topology.nc.  The text settings are not copied.
"""

import argparse
import json
import os

import numpy as np
from netCDF4 import Dataset

UPLAND, STREAM = 1, 6              # globalData.f90 domain types
TEMPERATURES = {"mLayerTemp", "scalarCanairTemp", "scalarCanopyTemp"}
GRU_COUNTS = {"nGlac", "nWtld", "nGrid"}


def remap(src, dst, hru_from, gru_from, edit):
    """Copy a NetCDF file whose hru (gru) axis takes the source rows hru_from (gru_from; None = zeros),
    then edit(name, dims, array) adjusts each hru- or gru-dimensioned variable."""
    with Dataset(src) as s, Dataset(dst, "w", format="NETCDF4") as d:
        d.setncatts({k: s.getncattr(k) for k in s.ncattrs()})
        size = {"hru": len(hru_from), "gru": len(gru_from)}
        for name, dim in s.dimensions.items():
            d.createDimension(name, None if dim.isunlimited() else size.get(name, len(dim)))
        for name, v in s.variables.items():
            fill = v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None
            o = d.createVariable(name, v.dtype, v.dimensions, zlib=True, complevel=4, fill_value=fill)
            o.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"})
            a = np.ma.filled(v[:], fill if fill is not None else 0)
            for axis_name, rows in (("hru", hru_from), ("gru", gru_from)):
                if axis_name not in v.dimensions:
                    continue
                ax = v.dimensions.index(axis_name)
                take = np.take(a, [r if r is not None else 0 for r in rows], axis=ax)
                for i, r in enumerate(rows):
                    if r is None:
                        idx = [slice(None)] * take.ndim
                        idx[ax] = i
                        take[tuple(idx)] = 0
                a = take
            if "hru" in v.dimensions or "gru" in v.dimensions:
                edit(name, v.dimensions, a)
            o[:] = a


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("domain", help="one-reach domain from add_stream_hru.py")
    p.add_argument("output", help="directory to write the network domain to")
    p.add_argument("spec", help="JSON spec of the reaches and the new GRUs")
    p.add_argument("--wscale", type=float, default=0.001, help="mizuRoute wscale, m m-1")
    a = p.parse_args()
    spec = json.load(open(a.spec))
    lapse = spec["lapse_rate"]
    reaches = {int(k): v for k, v in spec["reaches"].items()}
    new_grus = spec["grus"]

    src = os.path.join(a.domain, "summa_inputs")
    settings = os.path.join(a.output, "summa_inputs")
    os.makedirs(settings, exist_ok=True)
    os.makedirs(os.path.join(a.output, "forcing"), exist_ok=True)

    with Dataset(os.path.join(src, "attributes.nc")) as s:
        hru_ids = np.array(s["hruId"][:])
        gru_id = int(np.array(s["gruId"][:]).ravel()[0])
        if len(s.dimensions["gru"]) != 1:
            raise SystemExit("only a single-GRU domain is handled")
        elev = np.array(s["elevation"][:])
        stream_hru = int(np.flatnonzero(np.array(s["streamSegId"][:]) > 0)[0])
    if gru_id not in reaches:
        raise SystemExit(f"the original GRU {gru_id} needs a reach with its id")
    with Dataset(os.path.join(src, "coldState.nc")) as s:
        dom_type = np.array(s["domType"][0])
        dom_area = np.array(s["DOMarea"][0])
        n_soil = int(s["nSoil"][0, stream_hru, 1])
        n_lake = int(s["nLake"][0, stream_hru, 1])
    donor = int(np.flatnonzero(hru_ids == spec["donor_hru"])[0])
    if dom_type[donor, 0] != UPLAND:
        raise SystemExit("the donor HRU's first domain is not an upland")
    taken = sum(g["area"] for g in new_grus)
    if taken >= dom_area[donor, 0]:
        raise SystemExit(f"the new GRUs take {taken:.0f} m2, more than the donor's {dom_area[donor, 0]:.0f} m2 of upland")

    # every HRU in order: the original GRU's, then a land and a stream HRU per new GRU
    gru_of = {gru_id: gru_id}
    gru_of.update({g["id"]: g["id"] for g in new_grus})
    land_area = float(np.sum(dom_area[:stream_hru].sum(axis=1)))
    gru_land = {gru_id: land_area - taken}
    gru_land.update({g["id"]: g["area"] for g in new_grus})
    drains = {gru_id: gru_id}
    drains.update({g["id"]: g["reach"] for g in new_grus})
    upstream = {}
    def total_area(r):
        if r not in upstream:
            upstream[r] = sum(gru_land[g] for g, rr in drains.items() if rr == r) + \
                          sum(total_area(u) for u, v in reaches.items() if v["down"] == r)
        return upstream[r]
    width = {r: a.wscale * np.sqrt(total_area(r)) for r in reaches}
    stream_area = {r: reaches[r]["length"] * width[r] for r in reaches}

    next_id = int(hru_ids.max()) + 1
    rows = list(range(len(hru_ids)))                        # source HRU of each new HRU
    new = [dict(id=int(i), gru=gru_id, kind="old") for i in hru_ids]
    new[stream_hru]["kind"] = "stream"
    new[stream_hru]["reach"] = gru_id
    for g in new_grus:
        rows += [donor, stream_hru]
        new.append(dict(id=next_id, gru=g["id"], kind="land", area=g["area"], elevation=g["elevation"],
                        latitude=g["latitude"], longitude=g["longitude"]))
        new.append(dict(id=next_id + 1, gru=g["id"], kind="stream", reach=g["reach"]))
        next_id += 2
    for h in new:
        if h["kind"] == "stream":
            r = reaches[h["reach"]]
            h.update(area=stream_area[h["reach"]], elevation=r["elevation"], latitude=r["latitude"],
                     longitude=r["longitude"], length=r["length"], slope=r["slope"])
        h["dT"] = lapse * (elev[donor] - h["elevation"]) if h["kind"] in ("land", "stream") else 0.0
    gru_ids = [gru_id] + [g["id"] for g in new_grus]
    gru_rows = [0] + [None] * len(new_grus)

    # attributes
    def fill_attr(name, dims, arr):
        if dims == ("gru",):
            if name == "gruId":
                arr[:] = gru_ids
            elif name in GRU_COUNTS:
                arr[1:] = 0
            return
        if dims != ("hru",):
            return
        for i, h in enumerate(new):
            if name == "hruId":
                arr[i] = h["id"]
            elif name == "hru2gruId":
                arr[i] = h["gru"]
            elif name == "HRUarea" and i == donor:
                arr[i] -= taken
            elif h["kind"] == "old":
                continue
            elif name in ("HRUarea", "elevation", "latitude", "longitude"):
                arr[i] = h[name.replace("HRUarea", "area")]
            elif name == "streamSegId":
                arr[i] = h.get("reach", 0) if h["kind"] == "stream" else 0
            elif h["kind"] == "stream" and name in ("tan_slope", "contourLength"):
                arr[i] = h["slope"] if name == "tan_slope" else h["length"]
            elif name == "downHRUindex":
                arr[i] = 0

    remap(os.path.join(src, "attributes.nc"), os.path.join(settings, "attributes.nc"), rows, gru_rows, fill_attr)

    # cold state: a new land HRU is the donor's upland alone; temperatures move with elevation
    def fill_state(name, dims, arr):
        if name == "gruId":
            arr[:] = gru_ids
        if "hru" not in dims:
            return
        if dims == ("hru",):
            arr[:] = [h["id"] for h in new]
            return
        for i, h in enumerate(new):
            if h["kind"] == "old":
                if i == donor and name == "DOMarea":
                    arr[0, i, 0] -= taken
                continue
            if h["kind"] == "stream":
                if name == "mLayerTemp":
                    arr[n_lake:n_lake + n_soil, i, 1] += h["dT"]
                elif name == "DOMarea":
                    arr[0, i, 1] = h["area"]
                elif name == "DOMelev":
                    arr[0, i, 1] = h["elevation"]
                elif name == "DOMtan_slope":
                    arr[0, i, 1] = h["slope"]
                elif name == "DOMcontourLength":
                    arr[0, i, 1] = h["length"]
                continue
            # a new land HRU: the upland domain, the others unused
            arr[..., i, 1:] = -9999.0 if name == "DOMelev" else 0
            if name in TEMPERATURES:
                arr[..., i, 0] = np.where(arr[..., i, 0] > 0, arr[..., i, 0] + h["dT"], arr[..., i, 0])
            elif name == "DOMarea":
                arr[0, i, 0] = h["area"]
            elif name == "DOMelev":
                arr[0, i, 0] = h["elevation"]

    remap(os.path.join(src, "coldState.nc"), os.path.join(settings, "coldState.nc"), rows, gru_rows, fill_state)

    def fill_params(name, dims, arr):
        if name == "hruId":
            arr[:] = [h["id"] for h in new]
        elif name == "gruId":
            arr[:] = gru_ids
        elif dims == ("gru",):
            arr[1:] = arr[0]

    remap(os.path.join(src, "trialParams.nc"), os.path.join(settings, "trialParams.nc"), rows, [0] * len(gru_ids),
          fill_params)

    # forcing: each new HRU forced as its source, its air temperature moved to its own elevation
    names = [line.strip().strip("'") for line in open(os.path.join(src, "forcingFileList.txt")) if line.strip()]

    def fill_forcing(name, dims, arr):
        for i, h in enumerate(new):
            if name == "hruId":
                arr[i] = h["id"]
            elif h["kind"] == "old":
                continue
            elif name in ("latitude", "longitude"):
                arr[i] = h[name]
            elif name == "airtemp":
                arr[..., i] += h["dT"]

    for n in names:
        remap(os.path.join(a.domain, "forcing", n), os.path.join(a.output, "forcing", n), rows, [], fill_forcing)
    with open(os.path.join(settings, "forcingFileList.txt"), "w") as f:
        f.writelines(f"{n}\n" for n in names)

    # the network: each GRU feeds the reach it drains to
    segs = sorted(reaches)
    with Dataset(os.path.join(a.output, "topology.nc"), "w", format="NETCDF4") as d:
        d.createDimension("hru", len(gru_ids))
        d.createDimension("seg", len(segs))
        reach_of_gru = {g: drains[g] for g in gru_ids}
        hru_area = [gru_land[g] + sum(h["area"] for h in new if h["gru"] == g and h["kind"] == "stream")
                    for g in gru_ids]
        for name, dim, val, dtype, units in (("hruId", "hru", gru_ids, "i8", "-"),
                                             ("area", "hru", hru_area, "f8", "m2"),
                                             ("hruToSegId", "hru", [reach_of_gru[g] for g in gru_ids], "i8", "-"),
                                             ("segId", "seg", segs, "i8", "-"),
                                             ("downSegId", "seg", [reaches[r]["down"] for r in segs], "i8", "-"),
                                             ("length", "seg", [reaches[r]["length"] for r in segs], "f8", "m"),
                                             ("slope", "seg", [reaches[r]["slope"] for r in segs], "f8", "-")):
            v = d.createVariable(name, dtype, (dim,))
            v.units = units
            v[:] = val
        d.description = "reach network written by add_stream_network.py"

    for r in segs:
        print(f"reach {r} -> {reaches[r]['down']}: {reaches[r]['length']:.0f} m x {width[r]:.2f} m, "
              f"upstream {total_area(r) / 1e6:.2f} km2")
    print(f"{len(new)} HRUs in {len(gru_ids)} GRUs; donor upland {dom_area[donor, 0] - taken:.0f} m2 left")


if __name__ == "__main__":
    main()
