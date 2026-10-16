#!/usr/bin/env python3
"""Extract a basin-mean wet-period forcing file for the coupled Sagehen test cases.

The bundled August 2019 month is the driest in the record -- 3.8 mm over the month,
2.1 mm in its wettest 72 hours -- so nothing drains laterally and the cascade the
cases exist to guard carries no water. February 2017, the Sierra atmospheric river,
gives 128 mm over 72 hours at a mean 3.4 C, so almost all of it is rain.

Writes a one-HRU basin mean, which domain_sagehen4 and build_sagehen9.py tile out
to their own HRUs exactly as they did the August file. January 2017 is extracted
the same way for the spin-up that runs up to the event.

The bundled months come from SagehenCk/domain_sagehen/forcing/4_SUMMA_input, not
domain_sagehen_prcp_corrected, whose precipitation differs.

Usage:
    make_wet_forcing.py <source month .nc> <output .nc>          (source has one column per HRU)
"""

import sys

import numpy as np
from netCDF4 import Dataset

FORCING = ("airpres", "LWRadAtm", "SWRadAtm", "spechum", "windspd", "pptrate", "airtemp")


def main(src_path, out_path):
    with Dataset(src_path) as s, Dataset(out_path, "w", format="NETCDF4") as d:
        nTime = len(s.dimensions["time"])
        d.createDimension("time", nTime)
        d.createDimension("hru", 1)

        def copy(name, dims, data, src=None):
            ref = src if src is not None else s.variables.get(name)
            v = d.createVariable(name, "f8" if ref is None else ref.dtype, dims,
                                 zlib=True, complevel=4)
            if ref is not None:
                for att in ref.ncattrs():
                    if att != "_FillValue":
                        v.setncattr(att, ref.getncattr(att))
            v[:] = data

        copy("time", ("time",), np.array(s.variables["time"][:]))
        copy("hru", ("hru",), np.array([0]))
        copy("hruId", ("hru",), np.array([1.0]))
        # the basin mean sits at the area-weighted centre, which the source does not carry,
        # so take the plain mean of the HRU centres as the bundled August file does
        for name in ("latitude", "longitude"):
            copy(name, ("hru",), np.array([np.mean(np.array(s.variables[name][:]))]))
        for name in FORCING:
            col = np.ma.filled(s.variables[name][:], np.nan).astype(float)
            copy(name, ("time", "hru"), np.nanmean(col, axis=1)[:, None])
        copy("data_step", (), np.array(s.variables["data_step"][:]))

        d.description = ("basin-mean Sagehen forcing for the coupled wet-period tests, "
                         "written by make_wet_forcing.py from " + src_path.split("/")[-1])

    with Dataset(out_path) as d:
        p = np.array(d.variables["pptrate"][:]).ravel()
        t = np.array(d.variables["airtemp"][:]).ravel()
        # rank on rain, not precipitation: snow sits on the surface and drains nothing,
        # and these cases start from a snow-free cold state
        rain = np.where(t > 274.16, p, 0.0)
        c = np.convolve(rain, np.ones(72), "valid") * 3600.0
        k = int(c.argmax())
        print(f"wrote {out_path}: {len(p)} steps, {p.sum() * 3600:.1f} mm precipitation total")
        print(f"  rainiest 72 h starts at hour {k}: {c[k]:.1f} mm rain of "
              f"{p[k:k + 72].sum() * 3600:.1f} mm precipitation, mean air temperature "
              f"{t[k:k + 72].mean() - 273.15:.1f} C")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
