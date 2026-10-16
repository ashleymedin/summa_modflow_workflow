#!/usr/bin/env python3
"""Fetch groundwater levels for a USGS well and write them as a SUMMA calibration observation file.

The source is the U.S. Geological Survey Water Data API (api.waterdata.usgs.gov), parameter 72019,
depth to water below land surface.  Two records exist for a well:

  daily   daily values from a continuous recorder, one per day of the site's local standard time
  field   field measurements, a tape or sounder reading on each visit, often a few a year

Each is put on a daily grid on the model's clock, as for streamflow: --model-utc-offset states the clock
the SUMMA forcing is on, and each site day is moved onto it and stamped at its end.  A field reading
stands for the site day it was taken on; readings taken while the well was pumping, dry or flowing are
left out.  A day with no value is written as NaN.  The site's zone is looked up unless
--gauge-utc-offset gives it.

The file holds `h_obs`, the water level relative to land surface in m, positive up, so it rises as the
water table rises.  Its datum is not the model's, so a target scores it as an anomaly against
lowerBoundHead at the HRU the well sits in, with baseline_start and baseline_end.

A recorder in Michigan that publishes only the daily maximum depth, and a well near Truckee measured by
hand, both for a model forced in UTC:

    ./acquire_groundwater_level.py 421641085350601 2016-10-01 2018-09-30 USGS_421641085350601_daily_level.nc \\
        --model-utc-offset 0 --statistic 00001
    ./acquire_groundwater_level.py 391900120100001 2016-10-01 2023-09-30 USGS_391900120100001_field_level.nc \\
        --model-utc-offset 0 --record field

--features reads a saved API response (a GeoJSON FeatureCollection of the chosen record) in place of
the network, which is how the synthetic well of the Sagehen test case is read.
"""

import argparse
import datetime as dt
import json
import sys

import numpy as np

from acquire_streamflow import fetch_features, usgs_utc_offset
from observation_netcdf import write_observations

FT_TO_M = 0.3048

RECORDS = {
    "daily": {
        "url": "https://api.waterdata.usgs.gov/ogcapi/v0/collections/daily/items",
        "name": "U.S. Geological Survey, Water Data API daily values",
    },
    "field": {
        "url": "https://api.waterdata.usgs.gov/ogcapi/v0/collections/field-measurements/items",
        "name": "U.S. Geological Survey, Water Data API field measurements",
    },
}

# qualifiers of a reading that is not the resting water level
DISTURBED = ("pump", "dry", "flowing", "obstruct")


def fetch_usgs_level(record, site, start, end, statistic):
    """Every feature of a well's depth-to-water record between two days."""
    params = {
        "monitoring_location_id": f"USGS-{site}",
        "parameter_code": "72019",   # depth to water below land surface
        "time": f"{start}/{end}",
        "f": "json",
        "limit": 10000,
    }
    if record == "daily":
        params["statistic_id"] = statistic
    return fetch_features(RECORDS[record]["url"], params)


def daily_levels(features, record, statistic, site_offset):
    """Water level relative to land surface (m, positive up) by site day, readings on a day averaged."""
    readings = {}
    for f in features:
        props = f["properties"]
        if props["value"] is None or props.get("parameter_code", "72019") != "72019":
            continue
        if record == "daily" and props.get("statistic_id", statistic) != statistic:
            continue
        if any(word in q.lower() for q in (props.get("qualifier") or []) for word in DISTURBED):
            continue
        if props["unit_of_measure"] != "ft":
            sys.exit(f"unexpected USGS unit {props['unit_of_measure']!r}")
        if record == "daily":
            day = dt.date.fromisoformat(props["time"][:10])
        else:
            taken = dt.datetime.fromisoformat(props["time"]).astimezone(dt.timezone.utc).replace(tzinfo=None)
            day = (taken + dt.timedelta(hours=site_offset)).date()
        readings.setdefault(day, []).append(-float(props["value"]) * FT_TO_M)
    return {day: float(np.mean(values)) for day, values in readings.items()}


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
    parser.add_argument("--record", choices=sorted(RECORDS), default="daily",
                        help="continuous-recorder daily values or field measurements (default daily)")
    parser.add_argument("--statistic", default="00003",
                        help="daily statistic: 00003 mean (default), 00001 maximum depth, 00002 minimum depth")
    parser.add_argument("--features", help="a saved API response to read in place of the network")
    args = parser.parse_args()

    site_offset = args.gauge_utc_offset
    if site_offset is None:
        if args.features:
            sys.exit("give the site's standard-time UTC offset with --gauge-utc-offset when reading --features")
        site_offset = usgs_utc_offset(args.site)
    if site_offset is None:
        sys.exit(f"the time zone of {args.site} could not be looked up; give it with --gauge-utc-offset")

    if args.features:
        with open(args.features) as f:
            features = json.load(f)["features"]
        source = f"{RECORDS[args.record]['name']}, as saved in {args.features}"
    else:
        features = fetch_usgs_level(args.record, args.site, args.start, args.end, args.statistic)
        source = RECORDS[args.record]["name"]

    levels = daily_levels(features, args.record, args.statistic, site_offset)
    if not levels:
        hint = f" with statistic {args.statistic}; many recorders publish only 00001 or 00002" \
            if args.record == "daily" else ""
        sys.exit(f"no {args.record} water levels for {args.site} between {args.start} and {args.end}{hint}")

    # a site day starts at its local midnight, which on the model's clock is shifted by the zones' difference
    shift = dt.timedelta(hours=args.model_utc_offset - site_offset)
    days = [args.start + dt.timedelta(days=i) for i in range((args.end - args.start).days + 1)]
    starts = [dt.datetime.combine(day, dt.time()) + shift for day in days]
    ends = [start + dt.timedelta(days=1) for start in starts]
    values = np.asarray([levels.get(day, np.nan) for day in days], dtype="f8")

    what = "Daily mean" if args.record == "daily" else "Field-measured"
    where = args.site if args.features else f"USGS {args.site}"
    write_observations(
        args.output, starts, ends, values,
        varname="h_obs",
        units="m",
        long_name="observed water level relative to land surface, positive up",
        attrs={
            "title": f"{what} groundwater level at {where}, for SUMMA calibration",
            "source": source,
            "comment": f"Negated depth to water below land surface (72019) over the site's standard-time days "
                       f"(UTC{site_offset:+g}), on a clock at UTC{args.model_utc_offset:+g}, each stamped at "
                       f"the end of its day.",
            "history": f"created {dt.date.today()} by acquire_groundwater_level.py",
        },
        time_zone=f"UTC{args.model_utc_offset:+g}",
    )
    print(f"{args.output}: {np.isfinite(values).sum()} of {len(values)} days")


if __name__ == "__main__":
    main()
