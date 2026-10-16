"""Write an observed series in the NetCDF form every SUMMA calibration target reads.

The form is the one the bundled streamflow observations already use: a `time` coordinate with CF
units, and one named variable carrying its own units, missing values as NaN.  Keeping every kind of
observation in that one shape is what lets the model read them all with one routine.

Each time stamps the END of the span the value covers, on the clock the SUMMA forcing is on, which is
how SUMMA stamps its own time steps and how the calibration aligns a series: an observation is compared
against the simulation averaged over the span since the observation before it.  The span itself is
written too, as CF `time_bnds`, so the file says what it means without the reader having to assume it.
"""

import datetime as dt

import netCDF4 as nc
import numpy as np

# The reference the time coordinate is written against; this one matches the bundled streamflow
# observations, so the files read alike.
TIME_UNITS = "minutes since 1950-01-01"
TIME_REFERENCE = dt.datetime(1950, 1, 1)


def _minutes(times):
    return np.asarray([(t - TIME_REFERENCE).total_seconds() / 60.0 for t in times], dtype="i8")


def write_observations(path, starts, ends, values, varname, units, long_name, attrs=None, time_zone="UTC",
                       cell_methods="time: mean"):
    """Write one observed series.

    starts, ends: datetimes bounding the span each value covers (end exclusive), on the model's clock
    values:       the observations, NaN where there is none
    attrs:        global attributes, typically title / source / history
    time_zone:    the clock the times are on, which has to be the one the SUMMA forcing is on
    cell_methods: how each value covers its span, a mean over it or a change across it ("time: sum")
    """
    starts, ends = list(starts), list(ends)
    values = np.asarray(values, dtype="f8")
    if not (len(starts) == len(ends) == len(values)):
        raise ValueError("starts, ends and values differ in length")
    if any(e <= s for s, e in zip(starts, ends)):
        raise ValueError("every span has to end after it starts")
    if any(b <= a for a, b in zip(ends, ends[1:])):
        raise ValueError("spans have to be in increasing order")

    with nc.Dataset(path, "w", format="NETCDF4") as out:
        out.createDimension("time", len(values))
        out.createDimension("nbnds", 2)

        time_var = out.createVariable("time", "i8", ("time",))
        time_var.standard_name = "time"
        time_var.units = TIME_UNITS
        time_var.calendar = "standard"
        time_var.bounds = "time_bnds"
        time_var.time_zone = time_zone
        time_var.comment = "end of the span each value covers"
        time_var[:] = _minutes(ends)

        bnds_var = out.createVariable("time_bnds", "i8", ("time", "nbnds"))
        bnds_var.units = TIME_UNITS
        bnds_var[:, 0] = _minutes(starts)
        bnds_var[:, 1] = _minutes(ends)

        obs_var = out.createVariable(varname, "f8", ("time",), fill_value=np.nan)
        obs_var.units = units
        obs_var.long_name = long_name
        obs_var.cell_methods = cell_methods
        obs_var[:] = values

        for key, value in (attrs or {}).items():
            out.setncattr(key, value)
