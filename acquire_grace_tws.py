#!/usr/bin/env python3
"""Fetch the JPL GRACE/GRACE-FO mass anomaly for a basin and write it as a SUMMA calibration observation file.

The source is the JPL RL06.3 mascon grid with the Coastline Resolution Improvement (CRI) filter,
from PO.DAAC.  Downloading it needs a free NASA Earthdata login and the `earthaccess` package; a
grid already on disk can be named with --grid instead.

The value is the 0.5-degree grid cell holding the basin: its centroid, or a GRU's, taken as the
area-weighted mean of its HRUs from a SUMMA attributes file.  That is the right reduction for a
basin much smaller than a mascon (about 3 degrees), which is every basin this was written for;
a larger one would want the cells averaged over its outline instead.

The two basins of the bundled glacier test domain:

    ATTR=../test/test_calibration/gulkana_wolverine/summa_inputs/attributes.nc
    ./acquire_grace_tws.py --attributes $ATTR --gru 1 gulkana_grace_tws.nc
    ./acquire_grace_tws.py --attributes $ATTR --gru 2 wolverine_grace_tws.nc

The file holds `tws_obs`, equivalent water thickness in mm, as the anomaly JPL publishes: relative
to the product's own 2004.0-2009.999 mean.  A calibration target reading it references both series
to a common period with its `baseline_start` and `baseline_end`.  Each value covers the days its
monthly solution spans, which are not always a calendar month.
"""

import argparse
import datetime as dt
import os
import sys

import netCDF4 as nc
import numpy as np

from observation_netcdf import write_observations

SHORT_NAME = "TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.3_V4"
CM_TO_MM = 10.0


def download_grid(directory):
    """The latest JPL mascon grid from PO.DAAC, which holds every month to date."""
    try:
        import earthaccess
    except ImportError:
        sys.exit("downloading needs the earthaccess package; or name a grid already on disk with --grid")
    earthaccess.login()
    granules = earthaccess.search_data(short_name=SHORT_NAME)
    if not granules:
        sys.exit(f"PO.DAAC returned no granules of {SHORT_NAME}")
    return earthaccess.download(granules[-1:], directory)[0]


def gru_centroid(attributes, gru_id):
    """Area-weighted mean latitude and longitude of the HRUs of one GRU."""
    with nc.Dataset(attributes) as attr:
        in_gru = attr["hru2gruId"][:] == gru_id
        if not in_gru.any():
            sys.exit(f"GRU {gru_id} is not in {attributes}")
        area = attr["HRUarea"][:][in_gru]
        lat = np.average(attr["latitude"][:][in_gru], weights=area)
        lon = np.average(attr["longitude"][:][in_gru], weights=area)
    return float(lat), float(lon)


def cell_series(grid, lat, lon, apply_scale_factor):
    """The grid cell holding (lat, lon): the spans of each monthly solution and its anomaly in mm."""
    with nc.Dataset(grid) as grace:
        lon360 = lon % 360.0
        lat_bnds = grace["lat_bounds"][:]
        lon_bnds = grace["lon_bounds"][:]
        iLat = np.flatnonzero((lat_bnds[:, 0] <= lat) & (lat < lat_bnds[:, 1]))
        iLon = np.flatnonzero((lon_bnds[:, 0] <= lon360) & (lon360 < lon_bnds[:, 1]))
        if iLat.size != 1 or iLon.size != 1:
            sys.exit(f"no single grid cell holds latitude {lat}, longitude {lon}")
        iLat, iLon = int(iLat[0]), int(iLon[0])

        if grace["lwe_thickness"].units != "cm":
            sys.exit(f"unexpected lwe_thickness units {grace['lwe_thickness'].units!r}")
        values = np.ma.filled(grace["lwe_thickness"][:, iLat, iLon].astype("f8"), np.nan) * CM_TO_MM
        if apply_scale_factor:
            values *= float(grace["scale_factor"][iLat, iLon])

        # the first and last day each solution includes, so a span ends the day after its last
        bounds = nc.num2date(grace["time_bounds"][:], grace["time_bounds"].units,
                             only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        starts = [dt.datetime(b.year, b.month, b.day) for b in bounds[:, 0]]
        ends = [dt.datetime(b.year, b.month, b.day) + dt.timedelta(days=1) for b in bounds[:, 1]]
        cell = (float(grace["lat"][iLat]), float(grace["lon"][iLon]))
    return starts, ends, values, cell


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    where = parser.add_mutually_exclusive_group(required=True)
    where.add_argument("--gru", type=int, help="GRU id, located from --attributes")
    where.add_argument("--latlon", type=float, nargs=2, metavar=("LAT", "LON"), help="basin centroid")
    parser.add_argument("--attributes", help="SUMMA attributes file, for --gru")
    parser.add_argument("--grid", help="JPL mascon grid already on disk, rather than downloading it")
    parser.add_argument("--download-dir", default=".", help="where a downloaded grid goes")
    parser.add_argument("--scale-factor", action="store_true",
                        help="apply the grid's CRI land scale factor, which restores signal the filtering damps")
    parser.add_argument("output", help="NetCDF file to write")
    args = parser.parse_args()

    if args.gru is not None:
        if not args.attributes:
            sys.exit("--gru needs --attributes")
        lat, lon = gru_centroid(args.attributes, args.gru)
        where_text = f"GRU {args.gru} of {os.path.basename(args.attributes)}"
    else:
        lat, lon = args.latlon
        where_text = f"latitude {lat}, longitude {lon}"

    grid = args.grid or download_grid(args.download_dir)
    starts, ends, values, cell = cell_series(grid, lat, lon, args.scale_factor)

    write_observations(
        args.output, starts, ends, values,
        varname="tws_obs",
        units="mm",
        long_name="GRACE terrestrial water storage anomaly",
        attrs={
            "title": "GRACE terrestrial water storage anomaly, for SUMMA calibration",
            "source": f"{os.path.basename(grid)}, cell at latitude {cell[0]}, longitude {cell[1]}"
                      + (", CRI scale factor applied" if args.scale_factor else ""),
            "comment": f"{where_text}. Equivalent water thickness relative to the product's 2004.0-2009.999 mean.",
            "history": f"created {dt.date.today()} by acquire_grace_tws.py",
        },
    )
    print(f"{args.output}: {np.isfinite(values).sum()} months, cell at {cell[0]}N {cell[1]}E")


if __name__ == "__main__":
    main()
