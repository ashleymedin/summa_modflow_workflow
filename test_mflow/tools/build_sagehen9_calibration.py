#!/usr/bin/env python3
"""Build a domain_sagehen9 for calibration: the 9 D8 GRUs and reaches, over a long forcing period.

Two layouts of the same GRUs, stream HRUs, reaches and MODFLOW model:

  grid     one land HRU per active cell, as build_sagehen9.py (3387 land HRUs)
  lumped   one land HRU per GRU, draining to its stream HRU, mapped onto every cell of
           the GRU with equal weight (9 land HRUs)

The forcing is a directory of single-HRU basin-mean monthly files, tiled to every HRU.
The MODFLOW model is ex-gwf-sagehen-gwe (SFR streams, GWE heat transport), copied with a TDIS long
enough for the forcing, saving heads, temperatures and budgets at its last step only, and an OBS6
file writing the head at each observation cell every step.  GWE's water-table temperature is the
temperature of the baseflow SUMMA hands the stream (deepTherml = aquiferTemp).  CND reads the solid
thermal conductivity from kts1.txt, so a calibration can scale it.

Usage:
    build_sagehen9_calibration.py <grid|lumped> <basin-mean forcing dir> <output domain dir>
"""

import os
import re
import shutil
import sys

import numpy as np
from netCDF4 import Dataset

import build_sagehen9 as b9
from make_sagehen_synthetic_well import ROW as WELL_ROW, COL as WELL_COL

MF6_GWE = os.path.join(os.path.dirname(b9.MF6), "ex-gwf-sagehen-gwe")
# observation cells (1-based row, column), one per decile of depth to water off the channel, in every GRU;
# well1 is the synthetic well
OBS_CELLS = {"well1": (WELL_ROW, WELL_COL), "well2": (9, 73), "well3": (20, 57), "well4": (51, 34), "well5": (18, 45),
             "well6": (33, 39), "well7": (62, 67), "well8": (36, 67), "well9": (42, 79), "well10": (49, 7)}
DECISIONS = dict(b9.LATFLOW_DECISIONS, deepTherml="aquiferTemp", bcLowrTdyn="presTemp")


def main(layout, forcing_src, out_dir):
    top, active = b9.read_grid()
    fill, rec = b9.flow_directions(top, active)
    acc, order = b9.accumulate(rec, active)
    links, down, sub, chan = b9.build_network(top, active, rec, acc, order)
    nGRU = len(links)
    gru_id = np.arange(1, nGRU + 1, dtype="i8")
    stream_id = gru_id + b9.STREAM_ID_OFFSET

    # reach geometry, as build_sagehen9.py
    length = np.array([b9.reach_length(c) for c in links])
    width = b9.WSCALE * np.sqrt(np.array([acc[c[-1]] * b9.CELL_AREA for c in links]))
    stream_area = length * width
    slope = np.array([max((top[c[0]] - top[c[-1]]) / b9.reach_length(c), 1e-4) for c in links])

    # every active cell, GRU by GRU, with its lon/lat, local slope and aspect
    gcells = [[(i, j) for i in range(b9.NROW) for j in range(b9.NCOL) if active[i, j] and sub[i, j] == k]
              for k in range(nGRU)]
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

    if layout == "grid":
        cells = [c for cs in gcells for c in cs]
        gru_of = np.concatenate([[k] * len(cs) for k, cs in enumerate(gcells)])
        land_id = np.array([(i + 1) * 100 + (j + 1) for i, j in cells], dtype="i8")
        ix_of = {c: n for n, c in enumerate(cells)}
        area = np.full(len(cells), b9.CELL_AREA)
        elev = np.array([top[c] for c in cells])
        geom = np.array([cell_geom(c) for c in cells])
        contour = np.full(len(cells), b9.DXY)
        down_hru = np.empty(len(cells), dtype="i8")
        for n, c in enumerate(cells):
            r = tuple(rec[c])
            inside = rec[c][0] >= 0 and sub[r] == gru_of[n]
            down_hru[n] = land_id[ix_of[r]] if inside else stream_id[gru_of[n]]
        members = [[n] for n in range(len(cells))]
    else:
        gru_of = np.arange(nGRU)
        land_id = gru_id.copy()
        area = np.array([len(cs) * b9.CELL_AREA for cs in gcells])
        elev = np.array([np.mean([top[c] for c in cs]) for cs in gcells])
        geom = np.array([np.mean([cell_geom(c) for c in cs], axis=0) for cs in gcells])
        # the land HRU drains into its reach along both banks
        contour = 2.0 * length
        down_hru = stream_id.copy()
        members = [[k] for k in range(nGRU)]
    tan_slope, aspect, lon, lat = geom.T
    nLand = len(land_id)

    # SUMMA indexes HRUs GRU by GRU, each GRU's stream HRU last; the map is keyed by that index
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

    s_elev = np.array([np.mean([top[c] for c in cs]) for cs in links])
    s_ll = np.array([np.mean([cell_geom(c)[2:] for c in cs], axis=0) for cs in links])

    settings = os.path.join(out_dir, "settings", "SUMMA")
    forcing = os.path.join(out_dir, "forcing")
    for d in (settings, forcing):
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
        d.description = f"Sagehen, 9 D8 GRUs, {layout} land HRUs, by build_sagehen9_calibration.py"

    b9.write_cold_state(settings, src_set, cat, land_slot, stream_slot, nHRU, ids,
                        cat(area, stream_area), cat(elev, s_elev), cat(tan_slope, slope),
                        cat(aspect, np.zeros(nGRU)), cat(contour, length))
    b9.write_trial_params(settings, src_set, ids)
    nSteps = write_forcing(settings, forcing, forcing_src, ids, cat(lon, s_ll[:, 0]), cat(lat, s_ll[:, 1]))
    b9.write_topology(out_dir, gru_id, down, length, slope, area=np.array([len(cs) for cs in gcells]) * b9.CELL_AREA)
    write_map(out_dir, layout, gcells, gru_of, land_slot)
    write_settings(settings, src_set)
    write_modflow(out_dir, nSteps, layout)

    print(f"{layout}: {nLand} land HRUs + {nGRU} stream HRUs in {nGRU} GRUs, {nSteps} forcing steps")


def write_forcing(settings, forcing, forcing_src, ids, lon, lat):
    """Tile each single-HRU month to every HRU; returns the number of steps."""
    n, nSteps = len(ids), 0
    files = sorted(f for f in os.listdir(forcing_src) if f.endswith(".nc"))
    for fname in files:
        with Dataset(os.path.join(forcing_src, fname)) as s, \
             Dataset(os.path.join(forcing, fname), "w", format="NETCDF4") as d:
            nSteps += len(s.dimensions["time"])
            for dname, dim in s.dimensions.items():
                d.createDimension(dname, n if dname == "hru" else len(dim))
            for name, v in s.variables.items():
                fill = v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None
                out = d.createVariable(name, v.dtype, v.dimensions, fill_value=fill, zlib=True, complevel=4)
                for att in v.ncattrs():
                    if att != "_FillValue":
                        out.setncattr(att, v.getncattr(att))
                if "hru" not in v.dimensions:
                    out[:] = v[:]
                elif name in ("hruId", "hru", "longitude", "latitude"):
                    out[:] = {"hruId": ids, "hru": np.arange(n), "longitude": lon, "latitude": lat}[name]
                else:
                    out[:] = np.repeat(np.array(v[:]), n, axis=v.dimensions.index("hru"))
    with open(os.path.join(settings, "forcingFileList.txt"), "w") as f:
        f.writelines(f"'{fname}'\n" for fname in files)
    return nSteps


def write_map(out_dir, layout, gcells, gru_of, land_slot):
    """HRU -> cell map: a grid HRU is its own cell, a lumped HRU every cell of its GRU."""
    with open(os.path.join(out_dir, "hru2cell_map.txt"), "w") as f:
        f.write(f"# HRU -> MODFLOW cell map, {layout} layout, written by build_sagehen9_calibration.py\n")
        f.write("# columns: hru_index  cell_index((i-1)*ncol+j)  weight; stream HRUs are not mapped\n")
        if layout == "grid":
            cells = [c for cs in gcells for c in cs]
            for n, (i, j) in enumerate(cells):
                f.write(f"{land_slot[n] + 1} {i * b9.NCOL + j + 1} 1.0\n")
        else:
            for k, cs in enumerate(gcells):
                for i, j in cs:
                    f.write(f"{land_slot[k] + 1} {i * b9.NCOL + j + 1} 1.0\n")


def write_settings(settings, src_set):
    for fname in os.listdir(src_set):
        if fname.endswith((".txt", ".TBL")) and "deeproot" not in fname and "_wet" not in fname \
                and fname not in ("forcingFileList.txt", "fileManager.txt", "modelDecisions.txt", "outputControl.txt"):
            shutil.copy(os.path.join(src_set, fname), settings)
    out = []
    for line in open(os.path.join(src_set, "modelDecisions.txt")):
        key = line.split()[0] if line.split() else ""
        if key in DECISIONS:
            head, sep, tail = line.partition("!")
            line = f"{key:<32}{DECISIONS[key]:<16}{sep}{tail}" if sep else f"{key:<32}{DECISIONS[key]}\n"
        out.append(line)
    open(os.path.join(settings, "modelDecisions.txt"), "w").writelines(out)
    # daily means of what a calibration run is judged on; each trial writes its own file
    with open(os.path.join(settings, "outputControl.txt"), "w") as f:
        for v in ("scalarSWE", "scalarTotalRunoff", "averageRoutedRunoff", "scalarStreamTemp", "basin__TotalRunoff"):
            f.write(f"{v} | 24\n")


def write_zone_file(array_file, zone_file):
    """Zones of a text array by its distinct values, 1 the smallest; 0 where the value is 0 (inactive)."""
    a = np.loadtxt(array_file)
    classes = np.unique(a[a > 0])
    z = np.where(a > 0, np.searchsorted(classes, a) + 1, 0)
    np.savetxt(zone_file, z, fmt="%d")


def write_obs(mf6):
    """An OBS6 file with the head at each observation cell every step, listed in the GWF name file."""
    rows = "".join(f"  {name:<6} HEAD  1 {r} {c}\n" for name, (r, c) in OBS_CELLS.items())
    with open(os.path.join(mf6, "sagehen.obs"), "w") as f:
        f.write(f"""# Observations - the head at each observation cell, every step; well1 is the synthetic well.
BEGIN CONTINUOUS FILEOUT sagehen.obs.csv
{rows}END CONTINUOUS
""")
    nam = os.path.join(mf6, "sagehen.nam")
    text = open(nam).read().replace("END PACKAGES", "  OBS6   sagehen.obs   obs\nEND PACKAGES")
    open(nam, "w").write(text)


def write_cnd(mf6):
    """CND with the solid thermal conductivity as an array file, kts1.txt, so a calibration can scale it."""
    cnd = open(os.path.join(mf6, "sagehen_gwe.cnd")).read()
    kts = float(re.search(r"KTS\s+CONSTANT\s+(\S+)", cnd).group(1))
    np.savetxt(os.path.join(mf6, "kts1.txt"), np.full((b9.NROW, b9.NCOL), kts), fmt="%g")
    cnd = re.sub(r"KTS\s+CONSTANT\s+\S+", "KTS\n    OPEN/CLOSE  kts1.txt", cnd)
    open(os.path.join(mf6, "sagehen_gwe.cnd"), "w").write(cnd)


def write_modflow(out_dir, nSteps, layout):
    """ex-gwf-sagehen-gwe with one stress period of nSteps hours, saving its last step only.

    Only the grid layout solves the base exchange in MODFLOW: a lumped HRU's flux answers its cells' mean
    head, not each cell's, so a per-cell correction there is not its linearisation.
    """
    base = ("  base_ghb_package_name = 'GHBB'  ! the base exchange, solved by MODFLOW against its new head\n"
            "  base_drn_package_name = 'DRNB'\n") if layout == "grid" else ""
    mf6 = os.path.join(out_dir, "mf6")
    shutil.copytree(MF6_GWE, mf6, dirs_exist_ok=True)
    with open(os.path.join(mf6, "sagehen.tdis"), "w") as f:
        f.write(f"""# Time discretisation - SECONDS, one stress period of {nSteps} hourly steps, the length of the forcing.
BEGIN OPTIONS
  TIME_UNITS  SECONDS
END OPTIONS

BEGIN DIMENSIONS
  NPER  1
END DIMENSIONS

BEGIN PERIODDATA
#  perlen      nstp   tsmult
   {nSteps * 3600.0:.1f}   {nSteps}    1.0
END PERIODDATA
""")
    write_obs(mf6)
    write_cnd(mf6)
    write_zone_file(os.path.join(mf6, "kh1.txt"), os.path.join(mf6, "kzone1.txt"))
    write_zone_file(os.path.join(mf6, "sy1.txt"), os.path.join(mf6, "syzone1.txt"))
    with open(os.path.join(mf6, "sagehen.oc"), "w") as f:
        f.write("""# Output control - the last step only; a calibration trial is judged on SUMMA's output.
BEGIN OPTIONS
  BUDGET  FILEOUT  sagehen.cbc
  HEAD    FILEOUT  sagehen.hds
END OPTIONS

BEGIN PERIOD 1
  SAVE  HEAD    LAST
  SAVE  BUDGET  LAST
END PERIOD
""")
    with open(os.path.join(mf6, "sagehen_gwe.oc"), "w") as f:
        f.write("""# Output control - the last step only.
BEGIN OPTIONS
  BUDGET       FILEOUT  sagehen_gwe.cbc
  TEMPERATURE  FILEOUT  sagehen_gwe.ucn
END OPTIONS

BEGIN PERIOD 1
  SAVE  TEMPERATURE  LAST
  SAVE  BUDGET       LAST
END PERIOD
""")
    with open(os.path.join(out_dir, "summa_modflow6.config"), "w") as f:
        f.write(f"""&coupler
  mf6_model_name     = 'SAGEHEN'
  rch_package_name   = 'RCHA'
{base}  uzf_hold_days      = 7.0     ! a week of drainage per UZF wave, so deep cells keep few alive
  bnd_package_names  = 'SFR', 'DRN'
  bnd_package_roles  = 'baseflow', 'surface_discharge'
  gwe_model_name     = 'SAGEHEN_GWE'   ! its water-table temperature is SUMMA's scalarAquiferTemp
  esl_package_name   = 'ESL'
  uze_package_name   = 'UZE'
  map_file           = '{os.path.join(out_dir, "hru2cell_map.txt")}'
  mf6_epsg           = 0
  feedback           = .true.
/
""")


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in ("grid", "lumped"):
        sys.exit(__doc__)
    main(sys.argv[1], os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]))
