#!/usr/bin/env python3
"""Write SUMMA forcing for the Little River layouts from AORC v1.1 hourly files.

Each GRU takes the mean of the 1 km AORC pixels whose centres fall in it (the nearest pixel if none
does), and every HRU in a GRU takes its GRU's series. AORC stamps each hour at its end, as SUMMA does.
The lumped layout gets one file per calendar year; the grid layout only the months its event spans.
The AORC files are clipped by littleriver_data/fetch_aorc.py.

Usage:
    make_littleriver_forcing.py [DATA_DIR]
"""

import glob
import os
import sys

import geopandas as gpd
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from shapely.geometry import Point

HERE = os.path.dirname(os.path.abspath(__file__))
TEST = os.path.join(HERE, "..")
GRID_MONTHS = ("2013-01", "2013-02")
TIME_UNITS = "hours since 1990-01-01 00:00:00"


def gru_weights(lat, lon, gru_poly):
    """Pixel-membership matrix (nGRU, nlat*nlon), rows summing to one."""
    g = gru_poly.to_crs(4326)
    la, lo = np.meshgrid(lat, lon, indexing="ij")
    pts = gpd.GeoDataFrame(geometry=[Point(x, y) for x, y in zip(lo.ravel(), la.ravel())], crs=4326)
    hit = gpd.sjoin(pts, g[["gru", "geometry"]], predicate="within")
    w = np.zeros((len(g), lo.size))
    for idx, k in zip(hit.index, hit.gru):
        w[k - 1, idx] = 1.0
    for k in range(len(g)):
        if w[k].sum() == 0:
            c = g.geometry.iloc[k].centroid
            w[k, np.argmin((lo.ravel() - c.x) ** 2 + (la.ravel() - c.y) ** 2)] = 1.0
    return w / w.sum(axis=1, keepdims=True)


def summa_vars(ds, w):
    """GRU means in SUMMA's names and units, (time, nGRU)."""
    flat = lambda v: ds[v].values.reshape(ds.sizes["time"], -1) @ w.T
    u, v = flat("UGRD_10maboveground"), flat("VGRD_10maboveground")
    return {"pptrate": flat("APCP_surface") / 3600.0, "SWRadAtm": flat("DSWRF_surface"),
            "LWRadAtm": flat("DLWRF_surface"), "airtemp": flat("TMP_2maboveground"),
            "spechum": flat("SPFH_2maboveground"), "airpres": flat("PRES_surface"), "windspd": np.hypot(u, v)}


def write(path, times, data, hru_id, hru_gru, lon, lat):
    with Dataset(path, "w", format="NETCDF4") as d:
        d.createDimension("time", None)
        d.createDimension("hru", len(hru_id))
        t = d.createVariable("time", "f8", ("time",))
        t.units, t.calendar = TIME_UNITS, "proleptic_gregorian"
        t[:] = (times - np.datetime64("1990-01-01T00:00")) / np.timedelta64(1, "h")
        d.createVariable("hruId", "i8", ("hru",))[:] = hru_id
        d.createVariable("longitude", "f8", ("hru",))[:] = lon
        d.createVariable("latitude", "f8", ("hru",))[:] = lat
        for name, units in (("pptrate", "kg m-2 s-1"), ("SWRadAtm", "W m-2"), ("LWRadAtm", "W m-2"),
                            ("airtemp", "K"), ("spechum", "g g-1"), ("airpres", "Pa"), ("windspd", "m s-1")):
            v = d.createVariable(name, "f4", ("time", "hru"), zlib=True, complevel=4)
            v.units = units
            v[:] = data[name][:, hru_gru - 1]
        s = d.createVariable("data_step", "f8")
        s.units, s.long_name = "s", "data step length in seconds"
        s.assignValue(3600.0)
        d.description = "AORC v1.1 GRU means, written by make_littleriver_forcing.py"


def main(data):
    years = sorted(glob.glob(os.path.join(data, "aorc", "aorc_[0-9][0-9][0-9][0-9].nc")))
    if not years:
        sys.exit("no AORC files yet")
    gru = gpd.read_file(os.path.join(TEST, "domain_littleriver_lumped", "gru_polygons.gpkg"))
    with xr.open_dataset(years[0]) as ds:
        w = gru_weights(ds.latitude.values, ds.longitude.values, gru)

    for layout in ("domain_littleriver_lumped", "domain_littleriver"):
        settings = os.path.join(TEST, layout, "settings", "SUMMA")
        forcing = os.path.join(TEST, layout, "forcing", "SUMMA_input")
        with Dataset(os.path.join(settings, "attributes.nc")) as a:
            hru_id, hru_gru = a["hruId"][:], a["hru2gruId"][:].astype(int)
            lon, lat = a["longitude"][:], a["latitude"][:]
        names = []
        for path in years:
            year = os.path.basename(path)[5:9]
            with xr.open_dataset(path) as ds:
                if layout == "domain_littleriver":
                    ds = ds.sel(time=slice(GRID_MONTHS[0], f"{GRID_MONTHS[-1]}-28T23:00"))
                    if ds.sizes["time"] == 0:
                        continue
                ds = ds.load()
                name = f"littleriver_forcing_{year}.nc"
                write(os.path.join(forcing, name), ds.time.values, summa_vars(ds, w), hru_id, hru_gru, lon, lat)
            names.append(name)
        with open(os.path.join(settings, "forcingFileList.txt"), "w") as f:
            f.writelines(f"'{n}'\n" for n in names)
        print(f"{layout}: {len(names)} forcing files, {len(hru_id)} HRUs")


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Research/MODFLOW/littleriver_data")))
