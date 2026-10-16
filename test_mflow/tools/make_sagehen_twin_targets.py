#!/usr/bin/env python3
"""Write the twin experiment's synthetic targets from a reference calibration trial.

The reference is one trial of the calibration driver at known parameters, writing its hourly
SUMMA/mizuRoute output (write_timeseries, write_Qreach) and MODFLOW's OBS6 CSV.  This takes the
discharge and water temperature leaving one reach, and the head at every OBS6 observation, as daily
means over UTC days, each stamped at the end of its day, the way the calibration's observation
files are.  The CSV's last rows are the trial's time steps, as the driver reads them.

    q_obs        discharge, m3 s-1
    t_obs        water temperature, degC
    h_<name>     head at OBS6 observation <name>, m above the model datum

Usage:
    make_sagehen_twin_targets.py <timestep output .nc> <sagehen.obs.csv> <reach id> <output .nc>
"""

import datetime as dt
import sys

import numpy as np
from netCDF4 import Dataset, num2date

EPOCH = dt.datetime(1950, 1, 1)


def daily_means(times, values):
    """Means over UTC days of the values stamped in (midnight, next midnight]; whole days only."""
    days = {}
    for t, v in zip(times, values):
        days.setdefault((t - dt.timedelta(microseconds=1)).date(), []).append(v)
    per_day = max(len(v) for v in days.values())
    full = sorted(d for d, v in days.items() if len(v) == per_day)
    return full, {d: float(np.mean(days[d])) for d in full}


def main(timestep_nc, obs_csv, reach, output):
    with Dataset(timestep_nc) as d:
        tv = d["time"]
        times = [dt.datetime(*t.timetuple()[:6]) for t in num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"))]
        iseg = int(np.where(np.asarray(d["seg"][:]) == int(reach))[0][0])
        q = np.asarray(d["Q_reach"][:])[:, iseg, 0]
        temp = np.asarray(d["T_reach"][:])[:, iseg] - 273.15

    names = open(obs_csv).readline().strip().split(",")[1:]
    rows = np.loadtxt(obs_csv, delimiter=",", skiprows=1, ndmin=2)[-len(times):, 1:]
    if len(rows) != len(times):
        sys.exit(f"{obs_csv} holds {len(rows)} rows for {len(times)} time steps")

    series = {"q_obs": (q, "m3 s-1", "synthetic streamflow"),
              "t_obs": (temp, "degC", "synthetic stream water temperature")}
    for k, name in enumerate(names):
        series[f"h_{name.lower()}"] = (rows[:, k], "m", f"synthetic head at MODFLOW observation {name.lower()}")

    days = None
    with Dataset(output, "w", format="NETCDF4") as o:
        for vname, (values, units, long_name) in series.items():
            full, mean = daily_means(times, values)
            if days is None:
                days = full
                o.createDimension("time", len(days))
                o.createDimension("nbnds", 2)
                end = [int((dt.datetime.combine(d, dt.time()) + dt.timedelta(days=1) - EPOCH).total_seconds() // 60)
                       for d in days]
                t = o.createVariable("time", "i8", ("time",))
                t.standard_name, t.units, t.calendar = "time", "minutes since 1950-01-01", "standard"
                t.bounds, t.time_zone, t.comment = "time_bnds", "UTC+0", "end of the span each value covers"
                t[:] = end
                b = o.createVariable("time_bnds", "i8", ("time", "nbnds"))
                b.units = "minutes since 1950-01-01"
                b[:] = np.column_stack([np.asarray(end) - 1440, end])
            v = o.createVariable(vname, "f8", ("time",), fill_value=np.nan)
            v.units, v.long_name, v.cell_methods = units, long_name, "time: mean"
            v[:] = [mean.get(d, np.nan) for d in days]
        o.title = "Synthetic Sagehen twin-experiment targets, for SUMMA calibration"
        o.source = f"reference trial {timestep_nc}, {obs_csv}, reach {reach}"
        o.history = f"created {dt.date.today()} by make_sagehen_twin_targets.py"
    print(f"{output}: {len(days)} days of {', '.join(series)}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
