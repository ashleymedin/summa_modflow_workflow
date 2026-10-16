#!/usr/bin/env python3
"""Build domain_sagehen9: one SUMMA HRU per active ex-gwf-sagehen cell, in 9 GRUs.

Every HRU is one 90 m MODFLOW cell, so HRU area equals cell area by construction
and the HRU->cell map is the identity. The GRUs are the 9 subcatchments of the
channel network that D8 over DIS/TOP produces, each carrying a stream HRU for
its reach; downHRUindex is the D8 receiver, and a cell leaving its GRU drains to
that GRU's stream HRU, which run_oneGRU treats as a GRU outlet.

Flow directions come from a priority-flooded TOP -- the raw field has one
interior pit -- but every reported elevation is the raw DIS/TOP.

The forcing is the single-HRU series of domain_sagehen1 given to every HRU: the
source has no sub-GRU information, so there is nothing finer to distribute.

Each case starts the event from its own spin-up, the three weeks before the event,
which its run script makes on first use.

Usage:
    build_sagehen9.py [<output domain dir>]      (default ../domain_sagehen9)
"""

import heapq
import os
import shutil
import sys

import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
MF6 = os.path.join(HERE, "..", "ex-gwf-sagehen")
SRC = os.path.join(HERE, "..", "domain_sagehen1")

NROW, NCOL, DXY = 73, 81, 90.0
CELL_AREA = DXY * DXY
OUTLET = [(42, 79), (43, 79), (44, 79)]          # the three CHD cells, 0-based
OUTLET_CELL = (43, 79)                           # the lowest of them, the domain outlet
CHANNEL_THRESHOLD = 250                          # cells; the value that yields 9 links
NB = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
EPS = 1e-3

STREAM, UPLAND, VEG_WATER = 6, 1, 16             # globalData.f90 domain types, USGS water class
STREAM_ID_OFFSET = 100000
WSCALE = 0.001                                   # mizuRoute width = WSCALE * sqrt(upstream area)
LAKE_DEPTHS = [0.2, 0.8]
LAKE_TEMP_K = 278.0

# the wet period; bundled August 2019 is the driest month in the record and drains nothing
FORCING_FILES = ["NWAM_SUMMA_forcing_201701.nc", "NWAM_SUMMA_forcing_201702.nc"]
SIM_START, SIM_END = "2017-02-07 00:00", "2017-02-09 23:00"
# the spin-up ends the step before the event, where SUMMA writes its restart
SPIN_START, SPIN_END = "2017-01-15 00:00", "2017-02-06 23:00"

# grid georeference (EPSG:32611), from fitting the TOP valley line to the Copernicus
# mainstem; the DIS carries none. East-west is pinned by the outlet, north-south to
# about two cells, which only ever reaches the per-HRU latitude used for solar geometry.
XORIGIN, YTOP, EPSG = 214161.61, 4373231.81, 32611


def read_grid():
    def arr(name, dtype=float):
        with open(os.path.join(MF6, name)) as f:
            return np.array(f.read().split(), dtype=dtype).reshape(NROW, NCOL)
    top = arr("land1.txt")  # land surface; DIS/TOP is the soil-column base
    active = arr("idomain1.txt", int) == 1
    return top, active


def flow_directions(top, active):
    """Priority-flood from the outlet, then steepest descent on the filled surface."""
    fill = np.where(active, np.inf, np.nan)
    done = np.zeros((NROW, NCOL), bool)
    heap = []
    for i, j in OUTLET:
        fill[i, j] = top[i, j]
        heapq.heappush(heap, (top[i, j], i, j))
    while heap:
        e, i, j = heapq.heappop(heap)
        if done[i, j]:
            continue
        done[i, j] = True
        for di, dj in NB:
            a, b = i + di, j + dj
            if not (0 <= a < NROW and 0 <= b < NCOL) or not active[a, b] or done[a, b]:
                continue
            ne = max(top[a, b], e + EPS)
            if ne < fill[a, b]:
                fill[a, b] = ne
                heapq.heappush(heap, (ne, a, b))
    if done.sum() != active.sum():
        sys.exit("priority flood did not reach every active cell")

    rec = np.full((NROW, NCOL, 2), -1, int)
    for i in range(NROW):
        for j in range(NCOL):
            if not active[i, j] or (i, j) in OUTLET:
                continue
            lo, best = fill[i, j], None
            for di, dj in NB:                      # scan order as determine_runoff_conns_4mvr
                a, b = i + di, j + dj
                if not (0 <= a < NROW and 0 <= b < NCOL) or not active[a, b]:
                    continue
                if fill[a, b] < lo:
                    lo, best = fill[a, b], (a, b)
            if best is None:
                sys.exit(f"unresolved pit at ({i+1},{j+1})")
            rec[i, j] = best
    for c in OUTLET:                               # the flanking CHD cells follow the channel out
        if c != OUTLET_CELL:
            rec[c] = OUTLET_CELL
    return fill, rec


def accumulate(rec, active):
    indeg = np.zeros((NROW, NCOL), int)
    for i in range(NROW):
        for j in range(NCOL):
            if active[i, j] and rec[i, j, 0] >= 0:
                indeg[tuple(rec[i, j])] += 1
    queue = [(i, j) for i in range(NROW) for j in range(NCOL) if active[i, j] and indeg[i, j] == 0]
    acc, order, k = active.astype(float).copy(), [], 0
    while k < len(queue):
        i, j = queue[k]; k += 1
        order.append((i, j))
        r = rec[i, j]
        if r[0] >= 0:
            acc[r[0], r[1]] += acc[i, j]
            indeg[r[0], r[1]] -= 1
            if indeg[r[0], r[1]] == 0:
                queue.append((r[0], r[1]))
    if len(order) != active.sum():
        sys.exit("the D8 receiver graph contains a loop")
    return acc, order


def build_network(top, active, rec, acc, order):
    """Cut the channel network at its junctions; each link is a reach and owns a GRU."""
    chan = active & (acc >= CHANNEL_THRESHOLD)
    ninf = np.zeros((NROW, NCOL), int)
    for i in range(NROW):
        for j in range(NCOL):
            if chan[i, j] and rec[i, j, 0] >= 0 and chan[tuple(rec[i, j])]:
                ninf[tuple(rec[i, j])] += 1

    links = []
    for s in [(i, j) for i in range(NROW) for j in range(NCOL) if chan[i, j] and ninf[i, j] != 1]:
        cells, c = [s], s
        while True:
            r = tuple(rec[c]) if rec[c][0] >= 0 else None
            if r is None or not chan[r] or ninf[r] >= 2:
                break                              # a junction cell starts the link below
            cells.append(r)
            c = r
        links.append(cells)
    link_of = {c: k for k, cells in enumerate(links) for c in cells}

    sub = np.full((NROW, NCOL), -1, int)
    for c in reversed(order):                      # downstream-first, so the receiver is known
        sub[c] = link_of[c] if c in link_of else (sub[tuple(rec[c])] if rec[c][0] >= 0 else -1)

    down = []
    for cells in links:
        r = tuple(rec[cells[-1]])
        down.append(link_of[r] if rec[cells[-1]][0] >= 0 and r in link_of else -1)

    # number the reaches upstream-first, so the outlet reach is the last
    depth = []
    for k in range(len(links)):
        d, c = 0, k
        while down[c] >= 0:
            c, d = down[c], d + 1
        depth.append(d)
    new = {k: n for n, k in enumerate(sorted(range(len(links)), key=lambda k: -depth[k]))}
    links = [links[k] for k in sorted(new, key=lambda k: new[k])]
    down = [(-1 if down[k] < 0 else new[down[k]]) for k in sorted(new, key=lambda k: new[k])]
    sub = np.where(sub >= 0, np.vectorize(lambda k: new.get(k, -1))(sub), -1)
    return links, down, sub, chan


def reach_length(cells):
    if len(cells) < 2:
        return DXY
    return sum(DXY * (np.sqrt(2.0) if abs(a[0] - b[0]) + abs(a[1] - b[1]) == 2 else 1.0)
               for a, b in zip(cells, cells[1:]))


def main(out_dir):
    top, active = read_grid()
    fill, rec = flow_directions(top, active)
    acc, order = accumulate(rec, active)
    links, down, sub, chan = build_network(top, active, rec, acc, order)
    nGRU = len(links)

    # SUMMA indexes HRUs GRU by GRU, and the coupler's map_file is keyed by that index,
    # so the file is written in the same order: each GRU's cells, then its stream HRU
    cells, gru_of, slot = [], [], []
    for k in range(nGRU):
        ck = [(i, j) for i in range(NROW) for j in range(NCOL) if active[i, j] and sub[i, j] == k]
        slot += [("L", len(cells) + t) for t in range(len(ck))] + [("S", k)]
        cells += ck
        gru_of += [k] * len(ck)
    gru_of = np.array(gru_of)
    hru_id = np.array([(i + 1) * 100 + (j + 1) for i, j in cells], dtype="i8")
    ix_of = {c: n for n, c in enumerate(cells)}
    nLand = len(cells)
    gru_id = np.arange(1, nGRU + 1, dtype="i8")
    stream_id = gru_id + STREAM_ID_OFFSET
    land_slot = np.array([p for p, (t, n) in enumerate(slot) if t == "L"])
    stream_slot = np.array([p for p, (t, n) in enumerate(slot) if t == "S"])

    # reach geometry: length from the link's cell path, width from the area above its foot
    length = np.array([reach_length(c) for c in links])
    up_area = np.array([acc[c[-1]] * CELL_AREA for c in links])
    width = WSCALE * np.sqrt(up_area)
    stream_area = length * width
    slope = np.array([max((top[c[0]] - top[c[-1]]) / reach_length(c), 1e-4) for c in links])

    # A land HRU keeps its whole cell, and the reach is laid over the channel cells rather than
    # cut out of them: map weights are normalised per HRU, so a one-cell HRU always weighs 1
    # against a reach's 1/n and the areal share cannot be expressed. The reach is left out of
    # the map instead -- its water column belongs to mizuRoute, and SFR is what should carry
    # stream-aquifer exchange -- so every land HRU matches its cell exactly.
    area = np.full(nLand, CELL_AREA)

    # downslope HRU: the D8 receiver, or this GRU's stream HRU where the flow path leaves it
    down_hru = np.zeros(nLand, dtype="i8")
    for n, c in enumerate(cells):
        r = tuple(rec[c])
        if rec[c][0] < 0:
            down_hru[n] = stream_id[gru_of[n]]                 # the domain outlet cell
        elif sub[r] == gru_of[n]:
            down_hru[n] = hru_id[ix_of[r]]
        else:
            down_hru[n] = stream_id[gru_of[n]]

    elev = np.array([top[c] for c in cells], dtype="f8")
    tan_slope = np.empty(nLand)
    aspect = np.empty(nLand)
    for n, c in enumerate(cells):
        r = tuple(rec[c]) if rec[c][0] >= 0 else c
        di, dj = r[0] - c[0], r[1] - c[1]
        dist = DXY * np.hypot(di, dj) if (di or dj) else DXY
        tan_slope[n] = max((top[c] - top[r]) / dist, 1e-3)
        aspect[n] = np.degrees(np.arctan2(dj, -di)) % 360.0 if (di or dj) else 0.0

    xc = np.array([XORIGIN + (j + 0.5) * DXY for i, j in cells])
    yc = np.array([YTOP - (i + 0.5) * DXY for i, j in cells])
    from pyproj import Transformer
    lon, lat = Transformer.from_crs(f"EPSG:{EPSG}", "EPSG:4326", always_xy=True).transform(xc, yc)

    s_elev = np.array([np.mean([top[c] for c in cs]) for cs in links])
    s_lon = np.array([np.mean([lon[ix_of[c]] for c in cs]) for cs in links])
    s_lat = np.array([np.mean([lat[ix_of[c]] for c in cs]) for cs in links])

    settings = os.path.join(out_dir, "settings", "SUMMA")
    forcing = os.path.join(out_dir, "forcing", "SUMMA_input")
    os.makedirs(settings, exist_ok=True)
    os.makedirs(forcing, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "simulations", "run_1", "SUMMA"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "simulations", "spinup"), exist_ok=True)
    src_set = os.path.join(SRC, "settings", "SUMMA")
    nHRU = nLand + nGRU

    def cat(land, stream):
        land, stream = np.asarray(land), np.asarray(stream)
        out = np.empty(nHRU, dtype=np.result_type(land, stream))
        out[land_slot] = land
        out[stream_slot] = stream
        return out

    with Dataset(os.path.join(src_set, "attributes.nc")) as a:
        s_type = int(np.array(a.variables["soilTypeIndex"][:]).ravel()[0])
        v_type = int(np.array(a.variables["vegTypeIndex"][:]).ravel()[0])
        p_type = int(np.array(a.variables["slopeTypeIndex"][:]).ravel()[0])
        m_hgt = float(np.array(a.variables["mHeight"][:]).ravel()[0])

    with Dataset(os.path.join(settings, "attributes.nc"), "w", format="NETCDF4") as d:
        d.createDimension("hru", nHRU)
        d.createDimension("gru", nGRU)

        def put(name, dim, data, dtype, units, long_name):
            v = d.createVariable(name, dtype, (dim,), zlib=True, complevel=4)
            v.units, v.long_name = units, long_name
            v[:] = data

        put("hruId", "hru", cat(hru_id, stream_id), "i8", "-", "Index of hydrological response unit (HRU)")
        put("gruId", "gru", gru_id, "i8", "-", "Index of grouped response unit (GRU)")
        put("hru2gruId", "hru", cat(gru_id[gru_of], gru_id), "i8", "-", "Index of GRU to which the HRU belongs")
        put("downHRUindex", "hru", cat(down_hru, np.zeros(nGRU)), "i4", "-", "Index of downslope HRU (0 = basin outlet)")
        put("streamSegId", "hru", cat(np.zeros(nLand), gru_id), "i4", "-",
            "mizuRoute reach id of a stream HRU (0 = reach mapped from the GRU id)")
        put("longitude", "hru", cat(lon, s_lon), "f8", "Decimal degree east", "Longitude of HRU's centroid")
        put("latitude", "hru", cat(lat, s_lat), "f8", "Decimal degree north", "Latitude of HRU's centroid")
        put("elevation", "hru", cat(elev, s_elev), "f8", "m", "Mean elevation of HRU")
        put("HRUarea", "hru", cat(area, stream_area), "f8", "m^2", "Area of HRU")
        put("tan_slope", "hru", cat(tan_slope, slope), "f8", "m m-1", "Average tangent slope of HRU")
        put("contourLength", "hru", cat(np.full(nLand, DXY), length), "f8", "m", "Contour length of HRU")
        put("aspect", "hru", cat(aspect, np.zeros(nGRU)), "f8", "degrees",
            "Mean azimuth of HRU, degrees East of North")
        put("slopeTypeIndex", "hru", np.full(nHRU, p_type), "i4", "-", "Index defining slope")
        put("soilTypeIndex", "hru", np.full(nHRU, s_type), "i4", "-", "Index defining soil type")
        put("vegTypeIndex", "hru", cat(np.full(nLand, v_type), np.full(nGRU, VEG_WATER)), "i4", "-",
            "Index defining vegetation type")
        put("mHeight", "hru", np.full(nHRU, m_hgt), "f8", "m", "Measurement height above bare ground")
        d.description = ("one HRU per active ex-gwf-sagehen cell in 9 D8 subcatchments, "
                         "generated by build_sagehen9.py")

    write_cold_state(settings, src_set, cat, land_slot, stream_slot, nHRU,
                     cat(hru_id, stream_id), cat(area, stream_area), cat(elev, s_elev),
                     cat(tan_slope, slope), cat(aspect, np.zeros(nGRU)),
                     cat(np.full(nLand, DXY), length))
    write_trial_params(settings, src_set, cat(hru_id, stream_id))
    write_forcing(settings, forcing, src_set, SRC, cat(hru_id, stream_id),
                  cat(lon, s_lon), cat(lat, s_lat))
    write_map(out_dir, cells, land_slot)
    write_topology(out_dir, gru_id, down, length, slope,
                   np.array([(sub == k).sum() for k in range(nGRU)]) * CELL_AREA)
    write_text_settings(settings, src_set, os.path.basename(out_dir), nGRU)
    write_run_files(out_dir, os.path.basename(out_dir))

    print(f"{nLand} land HRUs + {nGRU} stream HRUs in {nGRU} GRUs")
    print(f"land HRU area {area.sum():.6e} m2 vs active cell area {active.sum() * CELL_AREA:.6e} m2 "
          f"(mismatch {100 * (area.sum() - active.sum() * CELL_AREA) / (active.sum() * CELL_AREA):+.3e} %)")
    print(f"reach planform {stream_area.sum():.4e} m2 on top, "
          f"{100 * stream_area.sum() / (active.sum() * CELL_AREA):.3f}% of the grid area")
    print(f"channel threshold {CHANNEL_THRESHOLD} cells; {chan.sum()} channel cells")
    print("\n gru  cells  area_km2  reach_m  width_m  slope   downSeg  HRUs_draining_out")
    for k in range(nGRU):
        n = int((sub == k).sum())
        nout = int((down_hru[gru_of == k] == stream_id[k]).sum())
        print(f"{k+1:4d} {n:6d} {n * CELL_AREA / 1e6:9.3f} {length[k]:8.0f} {width[k]:8.2f} "
              f"{slope[k]:7.4f} {down[k]+1:8d} {nout:18d}")
    print(f"\nelevation {elev.min():.0f}-{elev.max():.0f} m, tan_slope {tan_slope.min():.4f}-{tan_slope.max():.4f}")


def write_cold_state(settings, src_set, cat, land_slot, stream_slot, nHRU,
                     hru_id, area, elev, tan_slope, aspect, contour):
    with Dataset(os.path.join(src_set, "coldState.nc")) as c:
        nSoil = int(np.array(c.variables["nSoil"][:]).ravel()[0])
        nLake = len(LAKE_DEPTHS)
        nToto = nSoil + nLake
        col = {k: np.array(c.variables[k][:], dtype="f8").ravel() for k in
               ("mLayerDepth", "mLayerTemp", "mLayerVolFracLiq", "mLayerVolFracIce",
                "mLayerMatricHead", "iLayerHeight")}
        scal = {k: float(np.array(c.variables[k][:]).ravel()[0]) for k in
                ("dt_init", "scalarCanopyIce", "scalarCanopyLiq", "scalarSnowDepth", "scalarSWE",
                 "scalarSfcMeltPond", "scalarAquiferStorage", "scalarSnowAlbedo",
                 "scalarCanairTemp", "scalarCanopyTemp")}

    with Dataset(os.path.join(settings, "coldState.nc"), "w", format="NETCDF4") as d:
        d.createDimension("hru", nHRU)
        d.createDimension("dom", 2)
        d.createDimension("scalarv", 1)
        d.createDimension("midSoil", nSoil)
        d.createDimension("midToto", nToto)
        d.createDimension("ifcToto", nToto + 1)

        d.createVariable("hruId", "i8", ("hru",), zlib=True, complevel=4)[:] = hru_id
        v = d.createVariable("domType", "i4", ("hru", "dom"), zlib=True, complevel=4)
        v.long_name = "type defining the domain response unit type"
        dt = np.zeros((nHRU, 2), dtype="i4")
        dt[:, 0] = UPLAND
        dt[stream_slot, 1] = STREAM
        v[:] = dt
        for name, val in (("nSnow", 0), ("nSoil", nSoil), ("nGlce", 0)):
            d.createVariable(name, "i4", ("hru", "dom"), zlib=True, complevel=4)[:] = \
                np.full((nHRU, 2), val, dtype="i4")
        nl = np.zeros((nHRU, 2), dtype="i4")
        nl[stream_slot, 1] = nLake
        d.createVariable("nLake", "i4", ("hru", "dom"), zlib=True, complevel=4)[:] = nl

        def put_scal(name, land_val, stream_val, units="-"):
            v = d.createVariable(name, "f8", ("scalarv", "hru", "dom"), zlib=True, complevel=4)
            v.units = units
            arr = np.zeros((1, nHRU, 2))
            arr[0, :, :] = land_val                      # the upland domain everywhere
            arr[0, stream_slot, 1] = stream_val
            v[:] = arr

        for name, val in scal.items():
            stream = val
            if name in ("scalarCanopyIce", "scalarCanopyLiq", "scalarSnowDepth",
                        "scalarSWE", "scalarSfcMeltPond"):
                stream = 0.0
            if name in ("scalarCanairTemp", "scalarCanopyTemp"):
                stream = LAKE_TEMP_K
            put_scal(name, val, stream)

        def put_geom(name, value, units):
            # the stream HRU's upland domain is empty and its stream domain carries the value
            v = d.createVariable(name, "f8", ("scalarv", "hru", "dom"), zlib=True, complevel=4)
            v.units = units
            arr = np.zeros((1, nHRU, 2))
            arr[0, :, 0] = value
            arr[0, :, 1] = value
            arr[0, stream_slot, 1] = value[stream_slot]
            if name == "DOMarea":
                arr[0, stream_slot, 0] = 0.0
            v[:] = arr

        put_geom("DOMarea", area, "m2")
        put_geom("DOMelev", elev, "m")
        put_geom("DOMtan_slope", tan_slope, "-")
        put_geom("DOMaspect", aspect, "degrees")
        put_geom("DOMcontourLength", contour, "m")

        s_depth = np.concatenate([LAKE_DEPTHS, col["mLayerDepth"][:nSoil]])
        s_temp = np.concatenate([np.full(nLake, LAKE_TEMP_K), col["mLayerTemp"][:nSoil]])
        s_liq = np.concatenate([np.ones(nLake), col["mLayerVolFracLiq"][:nSoil]])
        s_ice = np.concatenate([np.zeros(nLake), col["mLayerVolFracIce"][:nSoil]])
        s_h = np.empty(nToto + 1)
        s_h[0] = -sum(LAKE_DEPTHS)
        for k in range(nToto):
            s_h[k + 1] = s_h[k] + s_depth[k]

        def put_lay(name, dim, nlay, land, stream, units):
            v = d.createVariable(name, "f8", (dim, "hru", "dom"), zlib=True, complevel=4)
            v.units = units
            arr = np.zeros((nlay, nHRU, 2))
            arr[:len(land), :, 0] = land[:, None]
            arr[:len(land), :, 1] = land[:, None]
            arr[:, stream_slot, 1] = stream[:, None]
            v[:] = arr

        put_lay("mLayerDepth", "midToto", nToto, col["mLayerDepth"][:nSoil], s_depth, "m")
        put_lay("mLayerTemp", "midToto", nToto, col["mLayerTemp"][:nSoil], s_temp, "K")
        put_lay("mLayerVolFracLiq", "midToto", nToto, col["mLayerVolFracLiq"][:nSoil], s_liq, "-")
        put_lay("mLayerVolFracIce", "midToto", nToto, col["mLayerVolFracIce"][:nSoil], s_ice, "-")
        put_lay("mLayerMatricHead", "midSoil", nSoil, col["mLayerMatricHead"][:nSoil],
                col["mLayerMatricHead"][:nSoil], "m")
        put_lay("iLayerHeight", "ifcToto", nToto + 1, col["iLayerHeight"][:nSoil + 1], s_h, "m")


def write_trial_params(settings, src_set, ids):
    with Dataset(os.path.join(src_set, "trialParams.nc")) as p, \
         Dataset(os.path.join(settings, "trialParams.nc"), "w", format="NETCDF4") as d:
        d.createDimension("hru", len(ids))
        for name, v in p.variables.items():
            out = d.createVariable(name, v.dtype, ("hru",), zlib=True, complevel=4)
            for att in v.ncattrs():
                out.setncattr(att, v.getncattr(att))
            out[:] = ids if name == "hruId" else np.full(len(ids), np.array(v[:]).ravel()[0])


def write_forcing(settings, forcing, src_set, src_domain, ids, lon, lat):
    n = len(ids)
    for fname in FORCING_FILES:
        with Dataset(os.path.join(src_domain, "forcing", "SUMMA_input", fname)) as s, \
             Dataset(os.path.join(forcing, fname), "w", format="NETCDF4") as d:
            for dname, dim in s.dimensions.items():
                d.createDimension(dname, n if dname == "hru" else
                                  (None if dim.isunlimited() else len(dim)))
            for name, v in s.variables.items():
                fill = v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None
                out = d.createVariable(name, v.dtype, v.dimensions, fill_value=fill,
                                       zlib=True, complevel=4)
                for att in v.ncattrs():
                    if att != "_FillValue":
                        out.setncattr(att, v.getncattr(att))
                if "hru" not in v.dimensions:
                    out[:] = v[:]
                elif name == "hruId":
                    out[:] = ids
                elif name == "hru":
                    out[:] = np.arange(n)
                elif name == "longitude":
                    out[:] = lon
                elif name == "latitude":
                    out[:] = lat
                else:
                    src = np.array(v[:])
                    ax = v.dimensions.index("hru")
                    out[:] = np.repeat(src, n, axis=ax) if src.shape[ax] == 1 else src
    with open(os.path.join(settings, "forcingFileList.txt"), "w") as f:
        f.writelines(f"'{fname}'\n" for fname in FORCING_FILES)


def write_map(out_dir, cells, land_slot):
    """The identity map: HRU index -> its own cell. A reach spreads over its channel cells."""
    with open(os.path.join(out_dir, "hru2cell_map.txt"), "w") as f:
        f.write("# HRU -> MODFLOW cell map for domain_sagehen9, written by build_sagehen9.py\n")
        f.write("# columns: hru_index  cell_index((i-1)*ncol+j)  weight\n")
        f.write("#\n# One land HRU per active cell, so this is the identity map and HRU area\n")
        f.write("# equals cell area exactly. The stream HRUs are deliberately absent: a reach's\n")
        f.write("# water column exchanges with the network, not with the cell beneath it, and\n")
        f.write("# co-mapping it would take an unowned share of that cell's recharge.\n")
        f.write("# hru_index is SUMMA's, which runs GRU by GRU; the attributes file is written\n")
        f.write("# in that order, so it is the row order of attributes.nc.\n")
        for n, (i, j) in enumerate(cells):
            f.write(f"{land_slot[n] + 1} {i * NCOL + j + 1} 1.0\n")


def write_topology(out_dir, gru_id, down, length, slope, area):
    path = os.path.join(out_dir, "topology.nc")
    n = len(gru_id)
    with Dataset(path, "w", format="NETCDF4") as d:
        d.createDimension("hru", n)
        d.createDimension("seg", n)

        def put(name, dim, data, dtype, units, long_name):
            v = d.createVariable(name, dtype, (dim,), zlib=True, complevel=4)
            v.units, v.long_name = units, long_name
            v[:] = data

        put("hruId", "hru", gru_id, "i8", "-", "routing HRU id, matched to SUMMA GRU ids")
        put("area", "hru", area, "f8", "m2", "routing HRU area")
        put("hruToSegId", "hru", gru_id, "i8", "-", "id of the stream segment below each HRU")
        put("segId", "seg", gru_id, "i8", "-", "stream segment id")
        put("downSegId", "seg", np.array([0 if k < 0 else k + 1 for k in down], dtype="i8"),
            "i8", "-", "downstream segment id (<=0 at the outlet)")
        put("length", "seg", length, "f8", "m", "segment length")
        put("slope", "seg", slope, "f8", "-", "segment slope")
        d.description = ("Sagehen channel network from D8 over the ex-gwf-sagehen DIS/TOP, "
                         "generated by build_sagehen9.py")


# modLatFlow constrains the conductivity profile and the infiltration closure
LATFLOW_DECISIONS = {"groundwatr": "modLatflow", "hc_profile": "exp_prof", "infRateMax": "topmodel_GA"}


def write_text_settings(settings, src_set, domain, nGRU):
    for fname in os.listdir(src_set):
        # the decisions and file manager are written per variant below, not copied
        if fname.endswith((".txt", ".TBL")) and "deeproot" not in fname and "_wet" not in fname \
                and fname not in ("forcingFileList.txt", "fileManager.txt", "modelDecisions.txt"):
            shutil.copy(os.path.join(src_set, fname), settings)
    shutil.copy(os.path.join(src_set, "localParamInfo_deeproot.txt"), settings)

    def decisions(name, changes, src="modelDecisions.txt"):
        out = []
        for line in open(os.path.join(src_set, src)):
            key = line.split()[0] if line.split() else ""
            if key in changes:
                head, sep, tail = line.partition("!")
                out.append(f"{key:<32}{changes[key]:<16}{sep}{tail}" if sep else
                           f"{key:<32}{changes[key]}\n")
            else:
                out.append(line)
        open(os.path.join(settings, name), "w").writelines(out)

    # the pair differs on the groundwatr line alone, so running both isolates lateral flow
    decisions("modelDecisions_latflow.txt", LATFLOW_DECISIONS)
    decisions("modelDecisions_noLatflow.txt", dict(LATFLOW_DECISIONS, groundwatr="modflow"))
    # deep roots over the EVT model, without lateral flow so it pairs with the lumped run_sagehen1_wet_deeproot
    decisions("modelDecisions_deeproot.txt", dict(LATFLOW_DECISIONS, groundwatr="modflow"),
              src="modelDecisions_deeproot.txt")

    # SUMMA names its restart <prefix>_restart_<YYYYMMDDHH of the last step>_G<first>-<last GRU>.nc
    stamp = SPIN_END.replace("-", "").replace(" ", "")[:10]
    for tag in ("latflow", "noLatflow", "deeproot"):
        restart = f"spinup_{tag}_restart_{stamp}_G1-{nGRU}.nc"
        for spin in (False, True):
            fm = []
            for line in open(os.path.join(src_set, "fileManager.txt")):
                line = line.replace("domain_sagehen1", domain)
                if line.startswith("outFilePrefix"):
                    line = f"outFilePrefix        '{'spinup' if spin else 'run1'}_{tag}' !\n"
                if line.startswith("simStartTime"):
                    line = f"simStartTime         '{SPIN_START if spin else SIM_START}' !\n"
                if line.startswith("simEndTime"):
                    line = (f"simEndTime           '{SPIN_END}' ! 552 hourly steps, matches mf6/sagehen.tdis\n"
                            if spin else f"simEndTime           '{SIM_END}' ! 72 hourly steps\n")
                if line.startswith("outputPath") and spin:
                    line = f"outputPath           '../{domain}/simulations/spinup/' ! the restart lands here too\n"
                if line.startswith("initConditionFile") and not spin:
                    line = f"initConditionFile    '../../simulations/spinup/{restart}' ! Relative to settingsPath\n"
                if line.startswith("outputControlFile") and spin:
                    line = "outputControlFile    'outputControl_spinup.txt' ! Relative to settingsPath\n"
                if line.startswith("decisionsFile"):
                    line = f"decisionsFile        'modelDecisions_{tag}.txt' ! Relative to settingsPath\n"
                if line.startswith("globalHruParamFile") and tag == "deeproot":
                    line = "globalHruParamFile   'localParamInfo_deeproot.txt' ! Relative to settingsPath\n"
                fm.append(line)
            name = f"fileManager_{tag}_spinup.txt" if spin else f"fileManager_{tag}.txt"
            open(os.path.join(settings, name), "w").writelines(fm)

    # every flux in the water balance each step, so run totals can be summed rather than sampled
    oc = os.path.join(settings, "outputControl.txt")
    want = ["scalarRainPlusMelt", "scalarTotalET", "scalarSurfaceRunoff", "scalarInfiltration",
            "scalarSoilDrainage", "scalarSoilBaseflow", "scalarAquiferBaseflow", "scalarTotalRunoff",
            "scalarAquiferSeepage", "scalarAquiferTranspire", "scalarTranspireLimAqfr",
            "scalarAquiferRootFrac", "scalarCanopyTranspiration", "mLayerColumnInflow", "mLayerColumnOutflow",
            "scalarStreamTemp", "scalarStreamRunoff", "basin__TotalRunoff"]
    lines = [ln for ln in open(oc) if ln.split("|")[0].strip() not in want]
    open(oc, "w").writelines(lines + [f"{v} | 1\n" for v in want])
    # the spin-up keeps the same variables as daily means: 552 steps of every HRU is a third of a gigabyte
    lines = [ln.split("|")[0].rstrip() + " | 24\n" if "|" in ln else ln for ln in open(oc)]
    open(os.path.join(settings, "outputControl_spinup.txt"), "w").writelines(lines)


def write_run_files(out_dir, domain):
    with open(os.path.join(out_dir, "summa_modflow6.config"), "w") as f:
        f.write(f"""&coupler
  mf6_model_name     = 'SAGEHEN' ! GWF model name in mfsim.nam
  rch_package_name   = 'RCHA'    ! RCH package name, as in the GWF name file (upper case)
  bnd_package_names  = 'CHD', 'DRN'
  bnd_package_roles  = 'baseflow', 'surface_discharge' ! the land-surface drain returns as surface runoff
  map_file           = '../{domain}/hru2cell_map.txt' ! the identity map: one HRU per active cell
  mf6_epsg           = 0         ! no reprojection needed, an explicit map_file is supplied
  feedback           = .true.
/
""")
    with open(os.path.join(out_dir, "summa_modflow6_deeproot.config"), "w") as f:
        f.write(f"""&coupler
  mf6_model_name     = 'SAGEHEN' ! GWF model name in mfsim.nam
  rch_package_name   = 'RCHA'    ! RCH package name, as in the GWF name file (upper case)
  bflow_package_name = ''        ! superseded by the bnd_package_* table below
  bnd_package_names  = 'CHD', 'DRN', 'EVTA'
  bnd_package_roles  = 'baseflow', 'surface_discharge', 'gw_et'
  evt_package_name   = 'EVTA'    ! takes SUMMA's aquifer transpiration demand; the supply returns as gw_et
  map_file           = '../{domain}/hru2cell_map.txt' ! the identity map: one HRU per active cell
  mf6_epsg           = 0         ! no reprojection needed, an explicit map_file is supplied
  feedback           = .true.
/
""")
    base = domain.replace("domain_", "")
    run = os.path.join(os.path.dirname(out_dir), f"run_{base}_deeproot.sh")
    with open(run, "w") as f:
        f.write(f"""#!/bin/bash
# One HRU per MODFLOW cell, 9 GRUs, 6 m roots drawing on MODFLOW's EVT, without lateral flow;
# the distributed half of run_sagehen1_wet_deeproot.sh. See {domain}/README.md.
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c {domain}/summa_modflow6_deeproot.config \\
                      -r {domain}/simulations/spinup/heads_deeproot.bin \\
                      -s {domain}/settings/SUMMA/fileManager_deeproot_spinup.txt \\
                      ex-gwf-sagehen-ss \\
                      {domain}/settings/SUMMA/fileManager_deeproot.txt
""")
    os.chmod(run, 0o755)
    for tag, what in (("latflow", "lateral flow routed by SUMMA (groundwatr = modLatflow)"),
                      ("noLatflow", "without lateral flow (groundwatr = modflow)")):
        run = os.path.join(os.path.dirname(out_dir),
                           f"run_{base}.sh" if tag == "latflow" else f"run_{base}_{tag}.sh")
        with open(run, "w") as f:
            f.write(f"""#!/bin/bash
# One HRU per MODFLOW cell, 9 GRUs, {what}; see {domain}/README.md.
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c {domain}/summa_modflow6.config \\
                      -r {domain}/simulations/spinup/heads_{tag}.bin \\
                      -s {domain}/settings/SUMMA/fileManager_{tag}_spinup.txt \\
                      ex-gwf-sagehen \\
                      {domain}/settings/SUMMA/fileManager_{tag}.txt
""")
        os.chmod(run, 0o755)


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "domain_sagehen9")))
