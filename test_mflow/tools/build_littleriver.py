#!/usr/bin/env python3
"""Build ex-gwf-littleriver and two SUMMA layouts of the Little River Experimental Watershed, Georgia.

The basin above USGS 02317797 (Little River near Tifton, 334 km2), on a 200 m grid in UTM 17N, as
SWAT+ gwflow ran it (Bailey et al., 2020, Hydrology 7:75). Raw data come from DATA_DIR (default
~/Research/MODFLOW/littleriver_data): the NLDI basin, a 30 m 3DEP DEM and SoilGrids depth to bedrock.

MODFLOW: one convertible layer, TOP the base of a 2 m soil column, BOTM land surface less the depth
to bedrock (the Hawthorn), K 5 m/d and SY 0.12 everywhere (the middle of Bailey's calibrated zones),
SS 2e-6 as GSFLOW. UZF below the column, a land-surface DRN, SFR along the D8 channel network, and
an OBS head series at each surficial-aquifer USGS well. A steady first period under 184 mm/yr of
recharge, Bailey's, then hourly transient steps.

SUMMA: GRUs are the D8 subcatchments of the channel network, at Sagehen's cells per GRU. Two layouts:
domain_littleriver has one HRU per cell, draining to its D8 receiver within the GRU, for a wet month;
domain_littleriver_lumped has one HRU per GRU, mapped by build_map_file.py's area intersection, for
WY2011-2020. Forcing and observations are written by make_littleriver_forcing.py.

Usage:
    build_littleriver.py [DATA_DIR]
"""

import heapq
import json
import os
import shutil
import subprocess
import sys

import geopandas as gpd
import numpy as np
import rasterio
from netCDF4 import Dataset
from pyproj import Transformer
from rasterio import features
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject
from shapely.geometry import box, shape
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
TEST = os.path.join(HERE, "..")
SRC = os.path.join(TEST, "domain_sagehen1")
MF6_DIR = os.path.join(TEST, "ex-gwf-littleriver")
MODEL = "LITTLERIVER"

EPSG = 32617
X0, YTOP, DX = 236200.0, 3517600.0, 200.0
NROW, NCOL = 168, 127
CELL_AREA = DX * DX
GAUGE = (31.481694, -83.584111)                 # USGS 02317797
NB = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
EPS = 1e-3

CHANNEL_THRESHOLD = 190        # cells; 21 GRUs of about 400 cells, as Sagehen's 9 hold about 380
SOIL_DEPTHS = [0.025, 0.075, 0.15, 0.25, 0.5, 1.0]   # 2 m column
SOIL_DEPTH = sum(SOIL_DEPTHS)

K_MD, KV_MD, SY, SS = 5.0, 5.0, 0.12, 2.0e-6
MIN_THICK = 5.0                # m of aquifer below TOP at least
RECH_MEAN = 0.184 / 3.15576e7  # m s-1
THTS, THTR, EPS_BC, SURFDEP = 0.35, 0.05, 4.0, 0.1
NTRAIL, NWAVESETS = 3, 1000
RBED_DEPTH, RBED_THK, MANNING, WSCALE = 1.0, 0.5, 0.04, 0.001

SOIL_TYPE, VEG_TYPE = 4, 14         # ROSETTA loamy sand; IGBP cropland/natural vegetation mosaic
M_HEIGHT = 30.0                     # above the 16-20 m forest canopy

# (spin-up start, end; run start, end): the lumped layout runs WY2011-2020, the USGS gauge record from 2010,
# and the grid February 2013, the wettest month in it
DECADE = ("2009-10-01 00:00", "2010-09-30 23:00", "2010-10-01 00:00", "2020-09-30 23:00")
EVENT = ("2013-01-01 00:00", "2013-01-31 23:00", "2013-02-01 00:00", "2013-02-28 23:00")
NSTEP = 87696                       # hours in the longer of the two runs, with a day to spare


def grids(data):
    """Land surface, basin fraction and depth to bedrock on the 200 m grid; the 30 m DEM."""
    dst = from_origin(X0, YTOP, DX, DX)
    with rasterio.open(os.path.join(data, "dem30.tif")) as r:
        dem30, t30 = r.read(1).astype("f8"), r.transform
    land = np.zeros((NROW, NCOL))
    reproject(dem30, land, src_transform=t30, src_crs=EPSG, dst_transform=dst, dst_crs=EPSG,
              resampling=Resampling.average)
    basin = gpd.read_file(os.path.join(data, "basin.geojson")).to_crs(EPSG)
    in30 = features.rasterize(basin.geometry, out_shape=dem30.shape, transform=t30).astype("f8")
    frac = np.zeros((NROW, NCOL))
    reproject(in30, frac, src_transform=t30, src_crs=EPSG, dst_transform=dst, dst_crs=EPSG,
              resampling=Resampling.average)
    with rasterio.open(os.path.join(data, "bdticm.tif")) as r:
        bd = r.read(1).astype("f8")
        bd[bd < 0] = np.nan
        depth = np.full((NROW, NCOL), np.nan)
        reproject(bd / 100.0, depth, src_transform=r.transform, src_crs=r.crs, dst_transform=dst,
                  dst_crs=EPSG, resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan)
    depth[np.isnan(depth)] = np.nanmedian(depth)
    return land, frac >= 0.5, depth, dem30, t30


def outlet_cell():
    x, y = Transformer.from_crs(4326, EPSG, always_xy=True).transform(GAUGE[1], GAUGE[0])
    return int((YTOP - y) // DX), int((x - X0) // DX)


def flow_directions(top, active, outlet):
    """Priority-flood from the outlet over the active cells, then steepest descent; unreached cells drop out."""
    fill = np.full((NROW, NCOL), np.inf)
    done = np.zeros((NROW, NCOL), bool)
    fill[outlet] = top[outlet]
    heap = [(top[outlet], *outlet)]
    while heap:
        e, i, j = heapq.heappop(heap)
        if done[i, j]:
            continue
        done[i, j] = True
        for di, dj in NB:
            a, b = i + di, j + dj
            if 0 <= a < NROW and 0 <= b < NCOL and active[a, b] and not done[a, b]:
                ne = max(top[a, b], e + EPS)
                if ne < fill[a, b]:
                    fill[a, b] = ne
                    heapq.heappush(heap, (ne, a, b))
    active = active & done
    rec = np.full((NROW, NCOL, 2), -1, int)
    for i, j in zip(*np.nonzero(active)):
        if (i, j) == outlet:
            continue
        lo, best = fill[i, j], None
        for di, dj in NB:
            a, b = i + di, j + dj
            if 0 <= a < NROW and 0 <= b < NCOL and active[a, b] and fill[a, b] < lo:
                lo, best = fill[a, b], (a, b)
        rec[i, j] = best
    return active, rec


def accumulate(rec, active):
    indeg = np.zeros((NROW, NCOL), int)
    for i, j in zip(*np.nonzero(active)):
        if rec[i, j, 0] >= 0:
            indeg[tuple(rec[i, j])] += 1
    queue = [(i, j) for i, j in zip(*np.nonzero(active)) if indeg[i, j] == 0]
    acc, order, k = active.astype(float).copy(), [], 0
    while k < len(queue):
        c = queue[k]; k += 1
        order.append(c)
        r = tuple(rec[c])
        if r[0] >= 0:
            acc[r] += acc[c]
            indeg[r] -= 1
            if indeg[r] == 0:
                queue.append(r)
    if len(order) != active.sum():
        sys.exit("the D8 receiver graph contains a loop")
    return acc, order


def build_network(active, rec, acc, order):
    """Links between junctions of the channel network, numbered upstream-first; each owns a GRU."""
    chan = active & (acc >= CHANNEL_THRESHOLD)
    ninf = np.zeros((NROW, NCOL), int)
    for c in zip(*np.nonzero(chan)):
        if rec[c][0] >= 0 and chan[tuple(rec[c])]:
            ninf[tuple(rec[c])] += 1
    links = []
    for s in zip(*np.nonzero(chan & (ninf != 1))):
        cells, c = [s], s
        while rec[c][0] >= 0:
            r = tuple(rec[c])
            if not chan[r] or ninf[r] >= 2:
                break
            cells.append(r)
            c = r
        links.append(cells)
    link_of = {c: k for k, cells in enumerate(links) for c in cells}
    down = [link_of.get(tuple(rec[cells[-1]]), -1) if rec[cells[-1]][0] >= 0 else -1 for cells in links]
    depth = []
    for k in range(len(links)):
        d, c = 0, k
        while down[c] >= 0:
            c, d = down[c], d + 1
        depth.append(d)
    new = {k: n for n, k in enumerate(sorted(range(len(links)), key=lambda k: -depth[k]))}
    links = [links[k] for k in sorted(new, key=new.get)]
    down = [(-1 if down[k] < 0 else new[down[k]]) for k in sorted(new, key=new.get)]
    link_of = {c: k for k, cells in enumerate(links) for c in cells}
    sub = np.full((NROW, NCOL), -1, int)
    for c in reversed(order):
        sub[c] = link_of[c] if c in link_of else (sub[tuple(rec[c])] if rec[c][0] >= 0 else -1)
    return links, down, sub, chan


def main(data):
    land, active, depth, dem30, t30 = grids(data)
    outlet = outlet_cell()
    if not active[outlet]:
        sys.exit(f"gauge cell {outlet} is not active")
    active, rec = flow_directions(land, active, outlet)
    acc, order = accumulate(rec, active)
    links, down, sub, chan = build_network(active, rec, acc, order)
    nGRU = len(links)
    top = land - SOIL_DEPTH
    bot = land - np.maximum(depth, SOIL_DEPTH + MIN_THICK)
    wells = surficial_wells(data, active)
    write_modflow(land, top, bot, active, rec, acc, links, wells)

    to_ll = Transformer.from_crs(EPSG, 4326, always_xy=True)
    chan_len = np.array([DX * len(cs) for cs in links])

    # grid: one HRU per active cell, GRU by GRU, draining to its D8 receiver within the GRU
    cells = [c for k in range(nGRU) for c in zip(*np.nonzero(sub == k))]
    hid = {c: (c[0] + 1) * 1000 + c[1] + 1 for c in cells}
    xs = np.array([X0 + (j + 0.5) * DX for i, j in cells]); ys = np.array([YTOP - (i + 0.5) * DX for i, j in cells])
    lon, lat = to_ll.transform(xs, ys)
    dn, sl = [], []
    for c in cells:
        r = tuple(rec[c])
        same = rec[c][0] >= 0 and sub[r] == sub[c]
        dn.append(hid[r] if same else 0)
        dist = DX * (np.sqrt(2.0) if rec[c][0] >= 0 and abs(r[0] - c[0]) + abs(r[1] - c[1]) == 2 else 1.0)
        sl.append(max((land[c] - land[r]) / dist, 1e-3) if rec[c][0] >= 0 else 1e-3)
    grid = dict(hruId=[hid[c] for c in cells], gru=[sub[c] + 1 for c in cells], area=np.full(len(cells), CELL_AREA),
                elevation=[land[c] for c in cells], tan_slope=sl, lon=lon, lat=lat, down=dn,
                contour=np.full(len(cells), DX))
    write_summa("domain_littleriver", grid, nGRU, EVENT)
    with open(os.path.join(TEST, "domain_littleriver", "hru2cell_map.txt"), "w") as f:
        f.write("# HRU -> MODFLOW cell map for domain_littleriver: one HRU per active cell, the identity\n")
        for n, (i, j) in enumerate(cells, 1):
            f.write(f"{n} {i * NCOL + j + 1} 1.0\n")

    # lumped: one HRU per GRU, mapped by area intersection of the GRU polygons
    gpoly = [unary_union([box(X0 + j * DX, YTOP - (i + 1) * DX, X0 + (j + 1) * DX, YTOP - i * DX)
                          for i, j in zip(*np.nonzero(sub == k))]) for k in range(nGRU)]
    lumped_dir = os.path.join(TEST, "domain_littleriver_lumped")
    os.makedirs(lumped_dir, exist_ok=True)
    g = gpd.GeoDataFrame({"ihru": np.arange(1, nGRU + 1), "gru": np.arange(1, nGRU + 1)}, geometry=gpoly, crs=EPSG)
    g.to_file(os.path.join(lumped_dir, "gru_polygons.gpkg"), driver="GPKG")
    cen = g.geometry.representative_point()
    glon, glat = to_ll.transform(cen.x.to_numpy(), cen.y.to_numpy())
    lumped = dict(hruId=np.arange(1, nGRU + 1), gru=np.arange(1, nGRU + 1), area=g.geometry.area.to_numpy(),
                  elevation=[land[sub == k].mean() for k in range(nGRU)],
                  tan_slope=[np.mean([s for s, gg in zip(sl, grid["gru"]) if gg == k + 1]) for k in range(nGRU)],
                  lon=glon, lat=glat, down=np.zeros(nGRU), contour=2.0 * chan_len)
    write_summa("domain_littleriver_lumped", lumped, nGRU, DECADE)
    subprocess.run([sys.executable, os.path.join(HERE, "build_map_file.py"), "--hru",
                    os.path.join(lumped_dir, "gru_polygons.gpkg"), "--index-field", "ihru", "--dis",
                    os.path.join(MF6_DIR, "littleriver.dis"), "--grid-epsg", str(EPSG), "--out",
                    os.path.join(lumped_dir, "hru2cell_map.txt")], check=True)

    print(f"{active.sum()} active cells ({active.sum() * CELL_AREA / 1e6:.1f} km2), outlet cell {outlet}")
    print(f"{nGRU} GRUs, {active.sum() / nGRU:.0f} cells each on average")
    print(f"land {land[active].min():.1f}-{land[active].max():.1f} m, aquifer below TOP "
          f"{(top - bot)[active].min():.1f}-{(top - bot)[active].max():.1f} m (mean {(top - bot)[active].mean():.1f})")
    print(f"{sum(len(c) for c in links)} SFR reaches; {len(wells)} wells observed: " + ", ".join(w["id"] for w in wells))


def surficial_wells(data, active):
    """USGS surficial-aquifer wells (110QRNR, 110SFCL) in an active cell, with their cell."""
    meta = json.load(open(os.path.join(data, "wells_meta.json")))
    to_utm = Transformer.from_crs(4326, EPSG, always_xy=True)
    wells = []
    for wid, v in sorted(meta.items()):
        if (v["p"].get("aquifer_code") or "") not in ("110QRNR", "110SFCL"):
            continue
        lon, lat = v["geom"]["coordinates"][:2]
        x, y = to_utm.transform(lon, lat)
        i, j = int((YTOP - y) // DX), int((x - X0) // DX)
        if 0 <= i < NROW and 0 <= j < NCOL and active[i, j]:
            wells.append(dict(id=wid.replace("USGS-", ""), row=i, col=j))
    return wells


def write_modflow(land, top, bot, active, rec, acc, links, wells):
    if os.path.exists(MF6_DIR):
        shutil.rmtree(MF6_DIR)
    os.makedirs(MF6_DIR)
    p = lambda name: os.path.join(MF6_DIR, name)

    def arr(name, a, fmt="%.3f"):
        np.savetxt(p(name), np.where(active, a, 0.0) if fmt != "%d" else a, fmt=fmt)

    arr("land1.txt", land)
    arr("top1.txt", top)
    arr("bot1.txt", bot)
    arr("idomain1.txt", active.astype(int), "%d")
    arr("strt1.txt", np.maximum(land - 1.0, bot + 1.0))
    cells = list(zip(*np.nonzero(active)))
    kms = 1.0 / 86400.0
    x_ll, y_ll = X0, YTOP - NROW * DX

    open(p("mfsim.nam"), "w").write(f"""# Little River, Georgia, for the SUMMA coupler; written by tools/build_littleriver.py
BEGIN OPTIONS
END OPTIONS

BEGIN TIMING
  TDIS6  littleriver.tdis
END TIMING

BEGIN MODELS
  GWF6  littleriver.nam  {MODEL}
END MODELS

BEGIN EXCHANGES
END EXCHANGES

BEGIN SOLUTIONGROUP 1
  IMS6  littleriver.ims  {MODEL}
END SOLUTIONGROUP 1
""")
    open(p("littleriver.nam"), "w").write("""BEGIN OPTIONS
  NEWTON  UNDER_RELAXATION
  SAVE_FLOWS
END OPTIONS

BEGIN PACKAGES
  DIS6   littleriver.dis   dis
  IC6    littleriver.ic    ic
  NPF6   littleriver.npf   npf
  STO6   littleriver.sto   sto
  UZF6   littleriver.uzf   UZF
  RCH6   littleriver.rcha  RCHA
  SFR6   littleriver.sfr   SFR
  DRN6   littleriver.drn   DRN
  OBS6   littleriver.obs   obs
  OC6    littleriver.oc    oc
END PACKAGES
""")
    open(p("littleriver.tdis"), "w").write(f"""# A steady first period under the mean recharge, then hourly steps; one per SUMMA data step.
BEGIN OPTIONS
  TIME_UNITS  SECONDS
END OPTIONS

BEGIN DIMENSIONS
  NPER  2
END DIMENSIONS

BEGIN PERIODDATA
   1.0  1  1.0
   {NSTEP * 3600.0:.1f}  {NSTEP}  1.0
END PERIODDATA
""")
    open(p("littleriver.dis"), "w").write(f"""# 200 m cells in EPSG:{EPSG}; TOP is the base of SUMMA's {SOIL_DEPTH} m soil column, land surface in land1.txt.
# BOTM is land surface less the SoilGrids depth to bedrock, at least {MIN_THICK} m below TOP.
BEGIN OPTIONS
  LENGTH_UNITS  METERS
  XORIGIN  {x_ll:.1f}
  YORIGIN  {y_ll:.1f}
  ANGROT  0.0
END OPTIONS

BEGIN DIMENSIONS
  NLAY  1
  NROW  {NROW}
  NCOL  {NCOL}
END DIMENSIONS

BEGIN GRIDDATA
  DELR
    CONSTANT  {DX}
  DELC
    CONSTANT  {DX}
  TOP
    OPEN/CLOSE  top1.txt
  BOTM
    OPEN/CLOSE  bot1.txt
  IDOMAIN
    OPEN/CLOSE  idomain1.txt
END GRIDDATA
""")
    open(p("littleriver.ic"), "w").write("BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN GRIDDATA\n  STRT\n    OPEN/CLOSE  strt1.txt\nEND GRIDDATA\n")
    open(p("littleriver.npf"), "w").write(f"""# K {K_MD} m/d horizontal and vertical, in m/s.
BEGIN OPTIONS
  SAVE_FLOWS
END OPTIONS

BEGIN GRIDDATA
  ICELLTYPE
    CONSTANT  1
  K
    CONSTANT  {K_MD * kms:.6e}
  K33
    CONSTANT  {KV_MD * kms:.6e}
END GRIDDATA
""")
    open(p("littleriver.sto"), "w").write(f"""BEGIN OPTIONS
  SAVE_FLOWS
END OPTIONS

BEGIN GRIDDATA
  ICONVERT
    CONSTANT  1
  SS
    CONSTANT  {SS}
  SY
    CONSTANT  {SY}
END GRIDDATA

BEGIN PERIOD 1
  STEADY-STATE
END PERIOD

BEGIN PERIOD 2
  TRANSIENT
END PERIOD
""")
    open(p("littleriver.rcha"), "w").write(f"""# The coupler overwrites RECHARGE each transient step; the steady period takes {RECH_MEAN * 3.15576e10:.0f} mm/yr.
BEGIN OPTIONS
  READASARRAYS
  SAVE_FLOWS
END OPTIONS

BEGIN PERIOD 1
  RECHARGE
    CONSTANT  {RECH_MEAN:.6e}
END PERIOD

BEGIN PERIOD 2
  RECHARGE
    CONSTANT  0.0
END PERIOD
""")
    vks = KV_MD * kms
    thti = THTR + (THTS - THTR) * min(RECH_MEAN / vks, 1.0) ** (1.0 / EPS_BC)
    with open(p("littleriver.uzf"), "w") as f:
        f.write("# Below SUMMA's soil column; the coupler writes SINF_PVAR and returns REJINF. No ET, seepage or mover.\n")
        f.write(f"BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  NUZFCELLS  {len(cells)}\n"
                f"  NTRAILWAVES  {NTRAIL}\n  NWAVESETS  {NWAVESETS}\nEND DIMENSIONS\n\n"
                "BEGIN PACKAGEDATA\n# iuzno lay row col landflag ivertcon surfdep vks thtr thts thti eps\n")
        for n, (i, j) in enumerate(cells, 1):
            f.write(f"  {n:5d} 1 {i + 1:3d} {j + 1:3d} 1 0 {SURFDEP} {vks:.6e} {THTR} {THTS} {thti:.4f} {EPS_BC}\n")
        f.write("END PACKAGEDATA\n\nBEGIN PERIOD 1\n# iuzno finf pet extdp extwc ha hroot rootact\n")
        for n in range(1, len(cells) + 1):
            f.write(f"  {n:5d} 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")
        f.write("END PERIOD\n")
    with open(p("littleriver.drn"), "w") as f:
        f.write("# Land-surface drain, returned to SUMMA as seepage (role surface_discharge).\n")
        f.write(f"BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND  {len(cells)}\nEND DIMENSIONS\n\n"
                "BEGIN PERIOD 1\n")
        for i, j in cells:
            f.write(f"  1 {i + 1:3d} {j + 1:3d} {land[i, j]:9.3f} {K_MD * kms * CELL_AREA / 1.0:.6e}\n")
        f.write("END PERIOD\n")

    rcells = [c for cs in links for c in cs]
    rno = {c: n + 1 for n, c in enumerate(rcells)}
    dn = {c: (rno[tuple(rec[c])] if rec[c][0] >= 0 and tuple(rec[c]) in rno else 0) for c in rcells}
    ups = {n: [] for n in rno.values()}
    for c, d in dn.items():
        if d:
            ups[d].append(rno[c])

    def step(c):
        if not dn[c]:
            return DX
        r = rcells[dn[c] - 1]
        return DX * (np.sqrt(2.0) if abs(r[0] - c[0]) + abs(r[1] - c[1]) == 2 else 1.0)

    with open(p("littleriver.sfr"), "w") as f:
        f.write(f"# One reach per D8 channel cell, {len(rcells)}; bed top {RBED_DEPTH} m below land, {RBED_THK} m thick, K = the cell's K33.\n")
        f.write(f"BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  NREACHES  {len(rcells)}\nEND DIMENSIONS\n\n"
                "BEGIN PACKAGEDATA\n# rno lay row col rlen rwid rgrd rtp rbth rhk man ncon ustrf ndv\n")
        for c in rcells:
            n = rno[c]
            r = rcells[dn[c] - 1] if dn[c] else c
            grad = max((land[c] - land[r]) / step(c), 1e-4)
            width = WSCALE * np.sqrt(acc[c] * CELL_AREA)
            f.write(f"  {n:4d} 1 {c[0] + 1:3d} {c[1] + 1:3d} {step(c):8.3f} {width:7.4f} {grad:.6f} "
                    f"{land[c] - RBED_DEPTH:9.3f} {RBED_THK} {KV_MD * kms:.4e} {MANNING} "
                    f"{len(ups[n]) + (1 if dn[c] else 0)} 1.0 0\n")
        f.write("END PACKAGEDATA\n\nBEGIN CONNECTIONDATA\n")
        for c in rcells:
            n = rno[c]
            f.write(f"  {n:4d} " + " ".join(str(u) for u in ups[n]) + (f" -{dn[c]}" if dn[c] else "") + "\n")
        f.write("END CONNECTIONDATA\n")

    with open(p("littleriver.obs"), "w") as f:
        f.write("# Head at the cell of each USGS surficial-aquifer well, every step.\n")
        f.write("BEGIN OPTIONS\n  DIGITS  10\nEND OPTIONS\n\nBEGIN CONTINUOUS  FILEOUT  littleriver.head.csv\n")
        for w in wells:
            f.write(f"  W{w['id']}  HEAD  1 {w['row'] + 1} {w['col'] + 1}\n")
        f.write("END CONTINUOUS\n")
    open(p("littleriver.oc"), "w").write("""BEGIN OPTIONS
  BUDGET  FILEOUT  littleriver.cbc
  HEAD    FILEOUT  littleriver.hds
END OPTIONS

BEGIN PERIOD 1
  SAVE  HEAD  LAST
END PERIOD

BEGIN PERIOD 2
  SAVE  HEAD    FREQUENCY 24
  SAVE  BUDGET  FREQUENCY 24
  SAVE  HEAD    LAST
END PERIOD
""")
    shutil.copy(os.path.join(TEST, "ex-gwf-sagehen", "sagehen.ims"), p("littleriver.ims"))
    open(p(".gitignore"), "w").write("*.lst\n*.hds\n*.cbc\n*.grb\n*.csv\n")


def write_summa(name, h, nGRU, period):
    spin_start, spin_end, sim_start, sim_end = period
    out = os.path.join(TEST, name)
    settings = os.path.join(out, "settings", "SUMMA")
    src_set = os.path.join(SRC, "settings", "SUMMA")
    for sub_dir in (settings, os.path.join(out, "forcing", "SUMMA_input"),
                    os.path.join(out, "simulations", "run_1", "SUMMA"), os.path.join(out, "simulations", "spinup")):
        os.makedirs(sub_dir, exist_ok=True)
    n = len(h["hruId"])
    # forcing is rebuilt from AORC by make_littleriver_forcing.py
    open(os.path.join(out, ".gitignore"), "w").write("forcing/SUMMA_input/*.nc\nsimulations/\n")

    with Dataset(os.path.join(settings, "attributes.nc"), "w", format="NETCDF4") as d:
        d.createDimension("hru", n)
        d.createDimension("gru", nGRU)

        def put(var, dim, data, dtype, units, long_name):
            v = d.createVariable(var, dtype, (dim,))
            v.units, v.long_name = units, long_name
            v[:] = np.asarray(data)

        put("hruId", "hru", h["hruId"], "i8", "-", "Index of hydrological response unit (HRU)")
        put("gruId", "gru", np.arange(1, nGRU + 1), "i8", "-", "Index of grouped response unit (GRU)")
        put("hru2gruId", "hru", h["gru"], "i8", "-", "Index of GRU to which the HRU belongs")
        put("downHRUindex", "hru", h["down"], "i8", "-", "Index of downslope HRU (0 = basin outlet)")
        put("longitude", "hru", h["lon"], "f8", "Decimal degree east", "Longitude of HRU's centroid")
        put("latitude", "hru", h["lat"], "f8", "Decimal degree north", "Latitude of HRU's centroid")
        put("elevation", "hru", h["elevation"], "f8", "m", "Mean elevation of HRU")
        put("HRUarea", "hru", h["area"], "f8", "m^2", "Area of HRU")
        put("tan_slope", "hru", h["tan_slope"], "f8", "m m-1", "Average tangent slope of HRU")
        put("contourLength", "hru", h["contour"], "f8", "m", "Contour length of HRU")
        put("aspect", "hru", np.zeros(n), "f8", "degrees", "Mean azimuth of HRU, degrees East of North")
        put("slopeTypeIndex", "hru", np.ones(n), "i4", "-", "Index defining slope")
        put("soilTypeIndex", "hru", np.full(n, SOIL_TYPE), "i4", "-", "Index defining soil type")
        put("vegTypeIndex", "hru", np.full(n, VEG_TYPE), "i4", "-", "Index defining vegetation type")
        put("mHeight", "hru", np.full(n, M_HEIGHT), "f8", "m", "Measurement height above bare ground")
        d.description = f"Little River GA, {name}, written by build_littleriver.py"
    hru_id = np.asarray(h["hruId"])
    nSoil = len(SOIL_DEPTHS)
    with Dataset(os.path.join(settings, "coldState.nc"), "w", format="NETCDF4") as d:
        for dim, size in (("hru", n), ("midSoil", nSoil), ("midToto", nSoil), ("ifcToto", nSoil + 1), ("scalarv", 1)):
            d.createDimension(dim, size)
        d.createVariable("hruId", "i8", ("hru",))[:] = hru_id
        scal = dict(dt_init=3600.0, nSoil=nSoil, nSnow=0, scalarCanopyIce=0.0, scalarCanopyLiq=0.0,
                    scalarSnowDepth=0.0, scalarSWE=0.0, scalarSfcMeltPond=0.0, scalarAquiferStorage=1.0,
                    scalarSnowAlbedo=0.0, scalarCanairTemp=293.15, scalarCanopyTemp=293.15)
        for k, v in scal.items():
            d.createVariable(k, "i4" if k in ("nSoil", "nSnow") else "f8", ("scalarv", "hru"))[:] = np.full((1, n), v)
        z = np.concatenate([[0.0], np.cumsum(SOIL_DEPTHS)])
        lay = lambda name, dim, vals: d.createVariable(name, "f8", (dim, "hru")).__setitem__(
            slice(None), np.repeat(np.asarray(vals, "f8")[:, None], n, axis=1))
        lay("mLayerDepth", "midToto", SOIL_DEPTHS)
        lay("iLayerHeight", "ifcToto", z)
        lay("mLayerTemp", "midToto", np.full(nSoil, 293.15))
        lay("mLayerVolFracIce", "midToto", np.zeros(nSoil))
        lay("mLayerVolFracLiq", "midToto", np.full(nSoil, 0.25))
        lay("mLayerMatricHead", "midSoil", np.full(nSoil, -1.0))

    with Dataset(os.path.join(settings, "trialParams.nc"), "w", format="NETCDF4") as d:
        d.createDimension("hru", n)
        d.createVariable("hruId", "i8", ("hru",))[:] = hru_id

    for fname in ("TBL_GENPARM.TBL", "TBL_MPTABLE.TBL", "TBL_SOILPARM.TBL", "TBL_VEGPARM.TBL",
                  "localParamInfo.txt", "basinParamInfo.txt"):
        shutil.copy(os.path.join(src_set, fname), settings)
    lp = os.path.join(settings, "localParamInfo.txt")
    lines = []
    for line in open(lp):
        key = line.split("|")[0].strip() if "|" in line else ""
        if key == "rootingDepth":
            line = "rootingDepth              |       1.5000 |       0.0100 |      10.0000\n"
        lines.append(line)
    open(lp, "w").writelines(lines)

    def decisions(name, changes):
        out = []
        for line in open(os.path.join(src_set, "modelDecisions.txt")):
            key = line.split()[0] if line.split() else ""
            if key in changes:
                head, sep, tail = line.partition("!")
                out.append(f"{key:<32}{changes[key]:<16}{sep}{tail}")
            else:
                out.append(line)
        open(os.path.join(settings, name), "w").writelines(out)

    latflow = {"groundwatr": "modLatflow", "hc_profile": "exp_prof", "infRateMax": "topmodel_GA", "tmZoneInfo": "ncTime"}
    decisions("modelDecisions_latflow.txt", latflow)
    decisions("modelDecisions_noLatflow.txt", dict(latflow, groundwatr="modflow"))

    want = ["averageRoutedRunoff", "basin__TotalRunoff", "scalarTotalET", "scalarRainPlusMelt", "scalarSurfaceRunoff",
            "scalarInfiltration", "scalarSoilDrainage", "scalarSoilBaseflow", "scalarAquiferBaseflow",
            "scalarAquiferSeepage", "scalarTotalSoilWat"]
    open(os.path.join(settings, "outputControl.txt"), "w").writelines(f"{v:<22}| 24 | mean\n" for v in want)

    stamp = spin_end.replace("-", "").replace(" ", "")[:10]
    for tag in ("latflow", "noLatflow"):
        restart = f"spinup_{tag}_restart_{stamp}_G{1:0{len(str(nGRU))}d}-{nGRU}.nc"   # SUMMA pads the first GRU to the last one's width
        for spin in (False, True):
            fm = f"""controlVersion       'SUMMA_FILE_MANAGER_V3.0.0' ! file manager version
simStartTime         '{spin_start if spin else sim_start}' !
simEndTime           '{spin_end if spin else sim_end}' !
tmZoneInfo           'utcTime' !
outFilePrefix        '{"spinup" if spin else "run1"}_{tag}' !
settingsPath         '../{name}/settings/SUMMA/' ! relative to MODFLOW_CASE
forcingPath          '../{name}/forcing/SUMMA_input/' !
outputPath           '../{name}/simulations/{"spinup" if spin else "run_1/SUMMA"}/' !
initConditionFile    '{"coldState.nc" if spin else "../../simulations/spinup/" + restart}' ! Relative to settingsPath
attributeFile        'attributes.nc' ! Relative to settingsPath
trialParamFile       'trialParams.nc' ! Relative to settingsPath
forcingListFile      'forcingFileList.txt' ! Relative to settingsPath
decisionsFile        'modelDecisions_{tag}.txt' ! Relative to settingsPath
outputControlFile    'outputControl.txt' ! Relative to settingsPath
globalHruParamFile   'localParamInfo.txt' ! Relative to settingsPath
globalGruParamFile   'basinParamInfo.txt' ! Relative to settingsPath
vegTableFile         'TBL_VEGPARM.TBL' ! Relative to settingsPath
soilTableFile        'TBL_SOILPARM.TBL' ! Relative to settingsPath
generalTableFile     'TBL_GENPARM.TBL' ! Relative to settingsPath
noahmpTableFile      'TBL_MPTABLE.TBL' ! Relative to settingsPath
"""
            fmname = f"fileManager_{tag}_spinup.txt" if spin else f"fileManager_{tag}.txt"
            open(os.path.join(settings, fmname), "w").write(fm)

    open(os.path.join(out, "summa_modflow6.config"), "w").write(f"""&coupler
  mf6_model_name     = '{MODEL}'
  rch_package_name   = 'RCHA'
  bnd_package_names  = 'SFR', 'DRN'
  bnd_package_roles  = 'baseflow', 'surface_discharge'
  map_file           = '../{name}/hru2cell_map.txt'
  mf6_epsg           = 0         ! an explicit map_file is supplied
  feedback           = .true.
/
""")
    base = name.replace("domain_littleriver", "run_littleriver")
    span = "WY2011-2020" if period == DECADE else "February 2013"
    for tag, what in (("latflow", "with lateral soil flow (groundwatr = modLatflow)"),
                      ("noLatflow", "without lateral flow (groundwatr = modflow)")):
        run = os.path.join(TEST, f"{base}.sh" if tag == "latflow" else f"{base}_{tag}.sh")
        open(run, "w").write(f"""#!/bin/bash
# Little River GA, {name}, {span}, {what}.
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c {name}/summa_modflow6.config \\
                      -r {name}/simulations/spinup/heads_{tag}.bin \\
                      -s {name}/settings/SUMMA/fileManager_{tag}_spinup.txt \\
                      ex-gwf-littleriver \\
                      {name}/settings/SUMMA/fileManager_{tag}.txt
""")
        os.chmod(run, 0o755)

if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Research/MODFLOW/littleriver_data")))
