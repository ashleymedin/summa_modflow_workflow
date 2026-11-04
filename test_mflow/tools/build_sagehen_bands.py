#!/usr/bin/env python3
"""Build a Sagehen domain of downslope bands: the 9 D8 GRUs of build_sagehen9.py, each cut into bands
of height above its channel (HAND), about BAND_CELLS cells a band, each band draining to the one below
and the lowest to the GRU's stream HRU.

It tests whether an HRU over a few tens of cells bounces under the explicit base exchange.  A
uniform flux over an HRU's cells sees about the harmonic mean of their storage, so a band holding
confined valley cells may sit past the C dt / S limit that the lumped builds stay well under.

  bands    every land HRU a band; the explicit exchange (no base GHB or DRN)
  lumped   one band per GRU, the 9 land HRUs of the calibration build; the explicit exchange
  hybrid   one HRU per cell whose spun-up head is within HYBRID_DEPTH of DIS/TOP, bands elsewhere;
           the base GHB and DRN keep rows for the one-to-one cells only, so only they are implicit

Each layout runs from its own copy of ex-gwf-sagehen, saving the head every step, beside the domain.

Usage:
    build_sagehen_bands.py <bands|lumped|hybrid> [<output domain dir>]   (default ../domain_sagehen_<layout>)
"""

import os
import shutil
import sys

import numpy as np
from netCDF4 import Dataset

import build_sagehen9 as b9

BAND_CELLS = 26        # cells a band, about GSFLOW's 128 HRUs over the grid
HYBRID_DEPTH = 0.1     # m below DIS/TOP, the soil-base flux blend, within which a cell is its own HRU
SPUN_HEADS = [os.path.join(b9.HERE, "..", "domain_sagehen9", "simulations", "spinup", f"heads_{t}.bin")
              for t in ("noLatflow", "latflow")]


def read_heads(path, active):
    """Heads of a coupler head restart file on the (row, column) grid; nan where inactive."""
    raw = open(path, "rb").read()
    assert raw[:8] == b"SUMMF6HD"
    n = int(np.frombuffer(raw[8:24], dtype="<i4")[3])
    h = np.full((b9.NROW, b9.NCOL), np.nan)
    h[active] = np.frombuffer(raw[24:24 + 8 * n], dtype="<f8")
    return h


def hand(top, active, rec, chan):
    """Height of each active cell above the channel cell its D8 path first reaches."""
    out = np.full((b9.NROW, b9.NCOL), np.nan)
    for c in zip(*np.nonzero(active)):
        path, r = [], c
        while not chan[r] and np.isnan(out[r]) and rec[r][0] >= 0:
            path.append(r)
            r = tuple(rec[r])
        base = top[r] - (0.0 if chan[r] or rec[r][0] < 0 else out[r])
        for p in path + [r]:
            if np.isnan(out[p]):
                out[p] = top[p] - base
    return out


def main(layout, out_dir):
    top, active = b9.read_grid()
    fill, rec = b9.flow_directions(top, active)
    acc, order = b9.accumulate(rec, active)
    links, down, sub, chan = b9.build_network(top, active, rec, acc, order)
    nGRU = len(links)
    hnd = hand(top, active, rec, chan)
    dis_top = np.loadtxt(os.path.join(b9.MF6, "top1.txt"))

    single = np.zeros_like(active)
    if layout == "hybrid":
        for path in SPUN_HEADS:
            single |= active & (read_heads(path, active) >= dis_top - HYBRID_DEPTH)

    # land HRUs GRU by GRU: one-to-one cells first, then bands from the channel up
    hrus, gru_of = [], []          # each HRU its list of cells
    for k in range(nGRU):
        ones = [c for c in zip(*np.nonzero(active & (sub == k) & single))]
        rest = sorted(zip(*np.nonzero(active & (sub == k) & ~single)), key=lambda c: hnd[c])
        nb = (1 if layout == "lumped" else max(1, round(len(rest) / BAND_CELLS))) if rest else 0
        bands = [list(b) for b in np.array_split(np.array(rest, dtype=int), nb)] if nb else []
        for c in ones:
            hrus.append([tuple(c)])
            gru_of.append(k)
        for b in bands:
            hrus.append([tuple(c) for c in b])
            gru_of.append(k)
    gru_of = np.array(gru_of)
    nLand = len(hrus)
    hru_of = {c: n for n, cs in enumerate(hrus) for c in cs}
    gru_id = np.arange(1, nGRU + 1, dtype="i8")
    stream_id = gru_id + b9.STREAM_ID_OFFSET
    land_id = np.empty(nLand, dtype="i8")
    for k in range(nGRU):
        idx = np.nonzero(gru_of == k)[0]
        land_id[idx] = (k + 1) * 10000 + np.arange(1, len(idx) + 1)

    # downslope: a one-to-one cell to its D8 receiver's HRU, a band to the next band down, the
    # lowest band, or a path leaving the GRU, to the stream HRU
    band_rank = {}
    for k in range(nGRU):
        bands = [n for n in np.nonzero(gru_of == k)[0] if len(hrus[n]) > 1 or not single[hrus[n][0]]]
        for r, n in enumerate(bands):
            band_rank[n] = (k, r, bands)
    down_hru = np.empty(nLand, dtype="i8")
    contour = np.zeros(nLand)
    for n, cs in enumerate(hrus):
        k = gru_of[n]
        if n in band_rank:
            _, r, bands = band_rank[n]
            down_hru[n] = land_id[bands[r - 1]] if r > 0 else stream_id[k]
        else:
            r = tuple(rec[cs[0]])
            inside = rec[cs[0]][0] >= 0 and sub[r] == k
            down_hru[n] = land_id[hru_of[r]] if inside else stream_id[k]
        # contour: the D8 edges leaving the HRU
        contour[n] = b9.DXY * sum(1 for c in cs if rec[c][0] < 0 or hru_of.get(tuple(rec[c])) != n)
    contour = np.maximum(contour, b9.DXY)

    from pyproj import Transformer
    tr = Transformer.from_crs(f"EPSG:{b9.EPSG}", "EPSG:4326", always_xy=True)

    def cell_geom(c):
        r = tuple(rec[c]) if rec[c][0] >= 0 else c
        di, dj = r[0] - c[0], r[1] - c[1]
        dist = b9.DXY * np.hypot(di, dj) if (di or dj) else b9.DXY
        tan = max((top[c] - top[r]) / dist, 1e-3)
        asp = np.degrees(np.arctan2(dj, -di)) % 360.0 if (di or dj) else 0.0
        lon, lat = tr.transform(b9.XORIGIN + (c[1] + 0.5) * b9.DXY, b9.YTOP - (c[0] + 0.5) * b9.DXY)
        return tan, asp, lon, lat

    area = np.array([len(cs) * b9.CELL_AREA for cs in hrus])
    elev = np.array([np.mean([top[c] for c in cs]) for cs in hrus])
    tan_slope, aspect, lon, lat = np.array([np.mean([cell_geom(c) for c in cs], axis=0) for cs in hrus]).T

    length = np.array([b9.reach_length(c) for c in links])
    width = b9.WSCALE * np.sqrt(np.array([acc[c[-1]] * b9.CELL_AREA for c in links]))
    stream_area = length * width
    slope = np.array([max((top[c[0]] - top[c[-1]]) / b9.reach_length(c), 1e-4) for c in links])
    s_elev = np.array([np.mean([top[c] for c in cs]) for cs in links])
    s_ll = np.array([np.mean([cell_geom(c)[2:] for c in cs], axis=0) for cs in links])

    land_slot, stream_slot, p = [], [], 0
    for k in range(nGRU):
        nk = int((gru_of == k).sum())
        land_slot += list(range(p, p + nk))
        stream_slot.append(p + nk)
        p += nk + 1
    land_slot, stream_slot = np.array(land_slot), np.array(stream_slot)
    nHRU = nLand + nGRU

    def cat(land, stream):
        land, stream = np.asarray(land), np.asarray(stream)
        out = np.empty(nHRU, dtype=np.result_type(land, stream))
        out[land_slot] = land
        out[stream_slot] = stream
        return out

    domain = os.path.basename(out_dir)
    settings = os.path.join(out_dir, "settings", "SUMMA")
    forcing = os.path.join(out_dir, "forcing", "SUMMA_input")
    for d in (settings, forcing, os.path.join(out_dir, "simulations", "run_1", "SUMMA"),
              os.path.join(out_dir, "simulations", "spinup")):
        os.makedirs(d, exist_ok=True)
    src_set = os.path.join(b9.SRC, "settings", "SUMMA")

    with Dataset(os.path.join(src_set, "attributes.nc")) as a:
        s_type = int(np.array(a.variables["soilTypeIndex"][:]).ravel()[0])
        v_type = int(np.array(a.variables["vegTypeIndex"][:]).ravel()[0])
        p_type = int(np.array(a.variables["slopeTypeIndex"][:]).ravel()[0])
        m_hgt = float(np.array(a.variables["mHeight"][:]).ravel()[0])

    ids = cat(land_id, stream_id)
    with Dataset(os.path.join(settings, "attributes.nc"), "w", format="NETCDF4") as d:
        d.createDimension("hru", nHRU)
        d.createDimension("gru", nGRU)

        def put(name, dim, data, dtype, units, long_name):
            v = d.createVariable(name, dtype, (dim,), zlib=True, complevel=4)
            v.units, v.long_name = units, long_name
            v[:] = data

        put("hruId", "hru", ids, "i8", "-", "Index of hydrological response unit (HRU)")
        put("gruId", "gru", gru_id, "i8", "-", "Index of grouped response unit (GRU)")
        put("hru2gruId", "hru", cat(gru_id[gru_of], gru_id), "i8", "-", "Index of GRU to which the HRU belongs")
        put("downHRUindex", "hru", cat(down_hru, np.zeros(nGRU)), "i4", "-", "Index of downslope HRU (0 = basin outlet)")
        put("streamSegId", "hru", cat(np.zeros(nLand), gru_id), "i4", "-",
            "mizuRoute reach id of a stream HRU (0 = reach mapped from the GRU id)")
        put("longitude", "hru", cat(lon, s_ll[:, 0]), "f8", "Decimal degree east", "Longitude of HRU's centroid")
        put("latitude", "hru", cat(lat, s_ll[:, 1]), "f8", "Decimal degree north", "Latitude of HRU's centroid")
        put("elevation", "hru", cat(elev, s_elev), "f8", "m", "Mean elevation of HRU")
        put("HRUarea", "hru", cat(area, stream_area), "f8", "m^2", "Area of HRU")
        put("tan_slope", "hru", cat(tan_slope, slope), "f8", "m m-1", "Average tangent slope of HRU")
        put("contourLength", "hru", cat(contour, length), "f8", "m", "Contour length of HRU")
        put("aspect", "hru", cat(aspect, np.zeros(nGRU)), "f8", "degrees", "Mean azimuth of HRU, degrees East of North")
        put("slopeTypeIndex", "hru", np.full(nHRU, p_type), "i4", "-", "Index defining slope")
        put("soilTypeIndex", "hru", np.full(nHRU, s_type), "i4", "-", "Index defining soil type")
        put("vegTypeIndex", "hru", cat(np.full(nLand, v_type), np.full(nGRU, b9.VEG_WATER)), "i4", "-",
            "Index defining vegetation type")
        put("mHeight", "hru", np.full(nHRU, m_hgt), "f8", "m", "Measurement height above bare ground")
        d.description = f"Sagehen, 9 D8 GRUs, {layout} land HRUs, by build_sagehen_bands.py"

    b9.write_cold_state(settings, src_set, cat, land_slot, stream_slot, nHRU, ids,
                        cat(area, stream_area), cat(elev, s_elev), cat(tan_slope, slope),
                        cat(aspect, np.zeros(nGRU)), cat(contour, length))
    b9.write_trial_params(settings, src_set, ids)
    b9.write_forcing(settings, forcing, src_set, b9.SRC, ids, cat(lon, s_ll[:, 0]), cat(lat, s_ll[:, 1]))
    b9.write_topology(out_dir, gru_id, down, length, slope, np.array([(sub == k).sum() for k in range(nGRU)]) * b9.CELL_AREA)
    b9.write_text_settings(settings, src_set, domain, nGRU)

    with open(os.path.join(out_dir, "hru2cell_map.txt"), "w") as f:
        f.write(f"# HRU -> MODFLOW cell map, {layout} layout, written by build_sagehen_bands.py\n")
        f.write("# columns: hru_index  cell_index((i-1)*ncol+j)  weight; stream HRUs are not mapped\n")
        for n, cs in enumerate(hrus):
            for i, j in cs:
                f.write(f"{land_slot[n] + 1} {i * b9.NCOL + j + 1} 1.0\n")

    mf6 = write_modflow(out_dir, layout, single)
    write_run_files(out_dir, domain, mf6, layout)

    nb = sum(1 for cs in hrus if len(cs) > 1 or not single[cs[0]])
    print(f"{layout}: {nLand} land HRUs ({nLand - nb} one-to-one, {nb} bands) + {nGRU} stream HRUs")
    sizes = [len(cs) for cs in hrus if len(cs) > 1 or not single[cs[0]]]
    print(f"band cells {min(sizes)}-{max(sizes)}, median {int(np.median(sizes))}; "
          f"land area {area.sum():.6e} m2 against {active.sum() * b9.CELL_AREA:.6e} m2 of cells")


def write_modflow(out_dir, layout, single):
    """ex-gwf-sagehen beside the domain, heads every step; in hybrid the base rows of one-to-one cells only."""
    mf6 = os.path.join(os.path.dirname(out_dir), f"ex-gwf-sagehen-{layout}")
    shutil.copytree(b9.MF6, mf6, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("*.hds", "*.cbc", "*.lst", "*.stdout", "*.grb", "*.csv"))
    with open(os.path.join(mf6, "sagehen.oc"), "w") as f:
        f.write("""# Output control - the head every step, to count hourly swings; the budget at the last step.
BEGIN OPTIONS
  BUDGET  FILEOUT  sagehen.cbc
  HEAD    FILEOUT  sagehen.hds
END OPTIONS

BEGIN PERIOD 1
  SAVE  HEAD    ALL
  SAVE  BUDGET  LAST
END PERIOD
""")
    if layout == "hybrid":
        for pkg in ("sagehen.ghbb", "sagehen.drnb"):
            lines = open(os.path.join(mf6, pkg)).read().splitlines(keepends=True)
            rows = [ln for ln in lines if ln.split() and ln.split()[0].isdigit()]
            keep = [ln for ln in rows if single[int(ln.split()[1]) - 1, int(ln.split()[2]) - 1]]
            out = []
            for ln in lines:
                if ln in rows:
                    continue
                out.append(f"  MAXBOUND  {len(keep)}\n" if ln.strip().startswith("MAXBOUND") else ln)
                if ln.strip().startswith("BEGIN PERIOD"):
                    out += keep
            open(os.path.join(mf6, pkg), "w").writelines(out)
    return mf6


def write_run_files(out_dir, domain, mf6, layout):
    base = ("  base_ghb_package_name = 'GHBB'  ! rows for the one-to-one cells only\n"
            "  base_drn_package_name = 'DRNB'\n") if layout == "hybrid" else ""
    with open(os.path.join(out_dir, "summa_modflow6.config"), "w") as f:
        f.write(f"""&coupler
  mf6_model_name     = 'SAGEHEN'
  rch_package_name   = 'RCHA'
{base}  uzf_hold_days      = 7.0
  bnd_package_names  = 'CHD', 'DRN'
  bnd_package_roles  = 'baseflow', 'surface_discharge'
  map_file           = '../{domain}/hru2cell_map.txt'
  mf6_epsg           = 0
  feedback           = .true.
/
""")
    for tag in ("latflow", "noLatflow"):
        run = os.path.join(os.path.dirname(out_dir), f"run_{domain.replace('domain_', '')}_{tag}.sh")
        with open(run, "w") as f:
            f.write(f"""#!/bin/bash
# Sagehen in downslope bands ({layout}), {tag}; see tools/build_sagehen_bands.py.
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c {domain}/summa_modflow6.config \\
                      -r {domain}/simulations/spinup/heads_{tag}.bin \\
                      -s {domain}/settings/SUMMA/fileManager_{tag}_spinup.txt \\
                      {os.path.basename(mf6)} \\
                      {domain}/settings/SUMMA/fileManager_{tag}.txt "$@"
""")
        os.chmod(run, 0o755)


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3) or sys.argv[1] not in ("bands", "lumped", "hybrid"):
        sys.exit(__doc__)
    layout = sys.argv[1]
    main(layout, os.path.abspath(sys.argv[2] if len(sys.argv) == 3 else
                                 os.path.join(b9.HERE, "..", f"domain_sagehen_{layout}")))
