#!/usr/bin/env python3
"""Fetch daily streamflow for a gauge and write it as a SUMMA calibration observation file.

Two sources, both public and needing no account:

  wsc   Water Survey of Canada, through the MSC GeoMet API (api.weather.gc.ca)
  usgs  U.S. Geological Survey, through the Water Data API (api.waterdata.usgs.gov)

Both publish each daily mean over a day of the gauge's local standard time.  The calibration compares
it against the simulation over the same span on the model's clock, which is whatever clock the SUMMA
forcing is on, so --model-utc-offset states that clock and each day is moved onto it and stamped at
its end, as SUMMA stamps its own time steps.  The gauge's zone is looked up (WSC's real-time station
list, the USGS site record) unless --gauge-utc-offset gives it.

The Bow at Banff, the gauge of the bundled test domain, whose forcing is on Mountain Standard Time:

    ./acquire_streamflow.py wsc 05BB001 1950-01-01 2021-12-31 CAN_05BB001_daily_flow_observations.nc \
        --model-utc-offset -7

A USGS gauge, given by site number, for a model forced in UTC:

    ./acquire_streamflow.py usgs 09380000 2000-01-01 2010-12-31 USGS_09380000_daily_flow_observations.nc \
        --model-utc-offset 0

The file holds `q_obs`, the daily mean in m3 s-1; USGS values arrive in ft3 s-1 and are converted.
A day with no value is written as NaN rather than dropped, so the record stays evenly spaced.
"""

import argparse
import csv
import datetime as dt
import io
import json
import sys
import urllib.parse
import urllib.request

import numpy as np

from observation_netcdf import write_observations

CFS_TO_CMS = 0.028316846592  # ft3 s-1 to m3 s-1

WSC_STATION_LIST = "https://dd.weather.gc.ca/today/hydrometric/doc/hydrometric_StationList.csv"
USGS_SITE = "https://api.waterdata.usgs.gov/ogcapi/v0/collections/monitoring-locations/items/USGS-{site}?f=json"

# standard-time UTC offsets, hours, of the zones USGS names its sites by
USGS_ZONES = {"AST": -4, "EST": -5, "CST": -6, "MST": -7, "PST": -8, "AKST": -9, "HST": -10,
              "SST": -11, "ChST": 10, "UTC": 0, "GMT": 0}

SOURCES = {
    "wsc": {
        "url": "https://api.weather.gc.ca/collections/hydrometric-daily-mean/items",
        "name": "Water Survey of Canada, MSC GeoMet hydrometric-daily-mean",
    },
    "usgs": {
        "url": "https://api.waterdata.usgs.gov/ogcapi/v0/collections/daily/items",
        "name": "U.S. Geological Survey, Water Data API daily values",
    },
}


def fetch_features(url, params):
    """Every feature of an OGC API query, following its `next` links."""
    features = []
    url = url + "?" + urllib.parse.urlencode(params)
    while url:
        with urllib.request.urlopen(url, timeout=120) as response:
            page = json.load(response)
        features.extend(page.get("features", []))
        url = next((link["href"] for link in page.get("links", []) if link.get("rel") == "next"), None)
    return features


def wsc_utc_offset(station):
    """Standard-time UTC offset, hours, of a WSC station, from the real-time station list."""
    with urllib.request.urlopen(WSC_STATION_LIST, timeout=120) as response:
        rows = csv.reader(io.TextIOWrapper(response, encoding="utf-8-sig"))
        next(rows)
        for row in rows:
            if row and row[0] == station:
                zone = row[-1].strip()                     # e.g. UTC-07:00
                sign = -1 if zone[3] == "-" else 1
                hours, minutes = zone[4:].split(":")
                return sign * (int(hours) + int(minutes) / 60.0)
    return None


def usgs_utc_offset(site):
    """Standard-time UTC offset, hours, of a USGS site, from its site record."""
    with urllib.request.urlopen(USGS_SITE.format(site=site), timeout=120) as response:
        zone = json.load(response)["properties"].get("time_zone_abbreviation")
    return USGS_ZONES.get(zone)


def fetch_wsc(station, start, end):
    """Daily discharge (m3 s-1) by date for a WSC station."""
    features = fetch_features(SOURCES["wsc"]["url"], {
        "STATION_NUMBER": station,
        "datetime": f"{start}/{end}",
        "f": "json",
        "limit": 10000,
    })
    return {
        dt.date.fromisoformat(f["properties"]["DATE"]): f["properties"]["DISCHARGE"]
        for f in features
        if f["properties"]["DISCHARGE"] is not None
    }


def fetch_usgs(site, start, end):
    """Daily mean discharge (m3 s-1) by date for a USGS site."""
    features = fetch_features(SOURCES["usgs"]["url"], {
        "monitoring_location_id": f"USGS-{site}",
        "parameter_code": "00060",   # discharge
        "statistic_id": "00003",     # daily mean
        "time": f"{start}/{end}",
        "f": "json",
        "limit": 10000,
    })
    flows = {}
    for f in features:
        props = f["properties"]
        if props["value"] is None:
            continue
        if props["unit_of_measure"] != "ft^3/s":
            sys.exit(f"unexpected USGS unit {props['unit_of_measure']!r}")
        flows[dt.date.fromisoformat(props["time"][:10])] = float(props["value"]) * CFS_TO_CMS
    return flows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", choices=sorted(SOURCES))
    parser.add_argument("gauge", help="WSC station number or USGS site number")
    parser.add_argument("start", type=dt.date.fromisoformat, help="first day, YYYY-MM-DD")
    parser.add_argument("end", type=dt.date.fromisoformat, help="last day, YYYY-MM-DD")
    parser.add_argument("output", help="NetCDF file to write")
    parser.add_argument("--model-utc-offset", type=float, required=True,
                        help="UTC offset, hours, of the clock the SUMMA forcing is on (0 for UTC)")
    parser.add_argument("--gauge-utc-offset", type=float,
                        help="standard-time UTC offset, hours, of the gauge's days, when it cannot be looked up")
    args = parser.parse_args()

    gauge_offset = args.gauge_utc_offset
    if gauge_offset is None:
        gauge_offset = (wsc_utc_offset if args.source == "wsc" else usgs_utc_offset)(args.gauge)
    if gauge_offset is None:
        sys.exit(f"the time zone of {args.gauge} could not be looked up; give it with --gauge-utc-offset")

    fetch = fetch_wsc if args.source == "wsc" else fetch_usgs
    flows = fetch(args.gauge, args.start, args.end)
    if not flows:
        sys.exit(f"no discharge for {args.gauge} between {args.start} and {args.end}")

    # a gauge day starts at its local midnight, which on the model's clock is shifted by the zones' difference
    shift = dt.timedelta(hours=args.model_utc_offset - gauge_offset)
    days = [args.start + dt.timedelta(days=i) for i in range((args.end - args.start).days + 1)]
    starts = [dt.datetime.combine(day, dt.time()) + shift for day in days]
    ends = [start + dt.timedelta(days=1) for start in starts]
    values = np.asarray([flows.get(day, np.nan) for day in days], dtype="f8")

    write_observations(
        args.output, starts, ends, values,
        varname="q_obs",
        units="m3 s-1",
        long_name="observed streamflow",
        attrs={
            "title": f"Daily streamflow at {args.source.upper()} {args.gauge}, for SUMMA calibration",
            "source": SOURCES[args.source]["name"],
            "comment": f"Daily means over the gauge's standard-time days (UTC{gauge_offset:+g}), "
                       f"on a clock at UTC{args.model_utc_offset:+g}, each stamped at the end of its day.",
            "history": f"created {dt.date.today()} by acquire_streamflow.py",
        },
        time_zone=f"UTC{args.model_utc_offset:+g}",
    )
    print(f"{args.output}: {np.isfinite(values).sum()} of {len(values)} days")


if __name__ == "__main__":
    main()
