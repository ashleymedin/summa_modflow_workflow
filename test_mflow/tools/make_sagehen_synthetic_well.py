#!/usr/bin/env python3
"""Write a synthetic Sagehen observation well, from a reference run's MODFLOW heads, as a saved USGS response.

Sagehen has no observation well: the one USGS well in the basin has no water levels, and Niswonger &
Prudic calibrated ex-gwf-sagehen to discharge alone.  So the well is made up.  It sits on the lower valley side
two cells from the channel above the outlet gauge (row 41, column 75, 1-based, in GRU 7), where the
water table is about 15 m down and moves with the seasons; the valley floor below it is held at land
surface by the drains.  Its record is the head of a reference coupled run there, sampled the way a
USGS recorder reports it: depth to water below land surface (72019) in ft to 0.01 ft, the mean over each
Pacific-standard day.  The calibration that scores it therefore knows the true answer, the parameters
of the reference run.

The output is a GeoJSON FeatureCollection shaped like the Water Data API's daily values, read by
acquire_groundwater_level.py --features as a real well's would be.  The site id says it is synthetic.

Usage:
    make_sagehen_synthetic_well.py <sagehen.hds> <run start, UTC, "YYYY-MM-DD HH:MM"> <output.json>

The heads file has to hold the one layer at several times a day, e.g. SAVE HEAD FREQUENCY 6 in its OC.
"""

import datetime as dt
import json
import os
import sys

import numpy as np

import build_sagehen9 as b9

SITE = "SYNTHETIC-SAGEHEN-1"
ROW, COL = 41, 75          # 1-based MODFLOW row and column of the well
SITE_UTC_OFFSET = -8       # Pacific standard time
FT_TO_M = 0.3048


def read_heads(path):
    """(seconds since the run start, head field) for every record of a one-layer MODFLOW 6 heads file."""
    header = np.dtype([("kstp", "<i4"), ("kper", "<i4"), ("pertim", "<f8"), ("totim", "<f8"),
                       ("text", "S16"), ("ncol", "<i4"), ("nrow", "<i4"), ("ilay", "<i4")])
    raw = open(path, "rb").read()
    times, fields, pos = [], [], 0
    while pos < len(raw):
        h = np.frombuffer(raw, header, 1, pos)[0]
        pos += header.itemsize
        n = int(h["nrow"]) * int(h["ncol"])
        fields.append(np.frombuffer(raw, "<f8", n, pos).reshape(int(h["nrow"]), int(h["ncol"])))
        times.append(float(h["totim"]))
        pos += 8 * n
    return np.asarray(times), fields


def main(hds, start, output):
    top, active = b9.read_grid()
    if not active[ROW - 1, COL - 1]:
        sys.exit(f"row {ROW} column {COL} is not an active cell")
    land = top[ROW - 1, COL - 1]

    seconds, fields = read_heads(hds)
    run_start = dt.datetime.strptime(start, "%Y-%m-%d %H:%M")
    utc = [run_start + dt.timedelta(seconds=s) for s in seconds]
    head = np.asarray([f[ROW - 1, COL - 1] for f in fields])

    # a site day runs from its local midnight to the next; a head saved at the midnight closing it belongs to it
    depth_by_day = {}
    for t, h in zip(utc, head):
        day = (t + dt.timedelta(hours=SITE_UTC_OFFSET) - dt.timedelta(microseconds=1)).date()
        depth_by_day.setdefault(day, []).append((land - h) / FT_TO_M)
    per_day = max(len(v) for v in depth_by_day.values())
    full = {day: v for day, v in depth_by_day.items() if len(v) == per_day}

    features = [{
        "type": "Feature",
        "properties": {
            "monitoring_location_id": SITE,
            "parameter_code": "72019",
            "statistic_id": "00003",
            "time": day.isoformat(),
            "value": f"{np.mean(v):.2f}",
            "unit_of_measure": "ft",
            "approval_status": "Synthetic",
            "qualifier": None,
        },
        "geometry": None,
    } for day, v in sorted(full.items())]

    with open(output, "w") as f:
        json.dump({
            "type": "FeatureCollection",
            "comment": f"Synthetic well: head of a reference coupled SUMMA / MODFLOW 6 run ({os.path.basename(hds)}) "
                       f"at row {ROW}, column {COL} of ex-gwf-sagehen, land surface {land:g} m, written by "
                       f"make_sagehen_synthetic_well.py.  Not a USGS record.",
            "features": features,
        }, f, indent=1)
    print(f"{output}: {len(features)} days of {per_day} heads each, {len(depth_by_day) - len(full)} partial days left out")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
