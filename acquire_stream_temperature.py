#!/usr/bin/env python3
"""Fetch daily mean water temperature for a USGS site and write it as a SUMMA calibration observation file.

The source is the U.S. Geological Survey Water Data API daily values (api.waterdata.usgs.gov), parameter
00010, statistic 00003.  As for streamflow, each daily mean covers a day of the site's local standard time;
--model-utc-offset states the clock the SUMMA forcing is on, and each day is moved onto it and stamped at its
end.  The site's zone is looked up unless --gauge-utc-offset gives it.

Sagehen Creek near Truckee, the outlet of the bundled Sagehen domains, whose forcing is in UTC:

    ./acquire_stream_temperature.py 10343500 2016-10-01 2018-09-30 USGS_10343500_daily_temperature.nc \
        --model-utc-offset 0

The file holds `t_obs`, the daily mean in degC, which a target scores against T_reach at the reach the
site sits on.  A day with no value is written as NaN rather than dropped.
"""

import argparse
import datetime as dt
import sys

import numpy as np

from acquire_streamflow import fetch_features, usgs_utc_offset, SOURCES
from observation_netcdf import write_observations


def fetch_usgs_temperature(site, start, end):
    """Daily mean water temperature (degC) by date for a USGS site."""
    features = fetch_features(SOURCES["usgs"]["url"], {
        "monitoring_location_id": f"USGS-{site}",
        "parameter_code": "00010",   # water temperature
        "statistic_id": "00003",     # daily mean
        "time": f"{start}/{end}",
        "f": "json",
        "limit": 10000,
    })
    temps = {}
    for f in features:
        props = f["properties"]
        if props["value"] is None:
            continue
        if props["unit_of_measure"] != "degC":
            sys.exit(f"unexpected USGS unit {props['unit_of_measure']!r}")
        temps[dt.date.fromisoformat(props["time"][:10])] = float(props["value"])
    return temps


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("site", help="USGS site number")
    parser.add_argument("start", type=dt.date.fromisoformat, help="first day, YYYY-MM-DD")
    parser.add_argument("end", type=dt.date.fromisoformat, help="last day, YYYY-MM-DD")
    parser.add_argument("output", help="NetCDF file to write")
    parser.add_argument("--model-utc-offset", type=float, required=True,
                        help="UTC offset, hours, of the clock the SUMMA forcing is on (0 for UTC)")
    parser.add_argument("--gauge-utc-offset", type=float,
                        help="standard-time UTC offset, hours, of the site's days, when it cannot be looked up")
    args = parser.parse_args()

    gauge_offset = args.gauge_utc_offset
    if gauge_offset is None:
        gauge_offset = usgs_utc_offset(args.site)
    if gauge_offset is None:
        sys.exit(f"the time zone of {args.site} could not be looked up; give it with --gauge-utc-offset")

    temps = fetch_usgs_temperature(args.site, args.start, args.end)
    if not temps:
        sys.exit(f"no daily mean water temperature for {args.site} between {args.start} and {args.end}")

    # a site day starts at its local midnight, which on the model's clock is shifted by the zones' difference
    shift = dt.timedelta(hours=args.model_utc_offset - gauge_offset)
    days = [args.start + dt.timedelta(days=i) for i in range((args.end - args.start).days + 1)]
    starts = [dt.datetime.combine(day, dt.time()) + shift for day in days]
    ends = [start + dt.timedelta(days=1) for start in starts]
    values = np.asarray([temps.get(day, np.nan) for day in days], dtype="f8")

    write_observations(
        args.output, starts, ends, values,
        varname="t_obs",
        units="degC",
        long_name="observed stream water temperature",
        attrs={
            "title": f"Daily mean water temperature at USGS {args.site}, for SUMMA calibration",
            "source": SOURCES["usgs"]["name"],
            "comment": f"Daily means over the site's standard-time days (UTC{gauge_offset:+g}), "
                       f"on a clock at UTC{args.model_utc_offset:+g}, each stamped at the end of its day.",
            "history": f"created {dt.date.today()} by acquire_stream_temperature.py",
        },
        time_zone=f"UTC{args.model_utc_offset:+g}",
    )
    print(f"{args.output}: {np.isfinite(values).sum()} of {len(values)} days")


if __name__ == "__main__":
    main()
