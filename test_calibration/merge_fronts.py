#!/usr/bin/env python3
"""Merge the Pareto fronts of several NSGA-II runs of one search (seeds) and pick from the merged front.

Prints each run's front size and how many of its members survive, then the best member of each objective and the
balanced one, nearest the ideal corner with each objective scaled by the merged front's span, and, for each
parameter, the merged front's range as a share of all the runs' samples' range: near 1 means the targets do not pin it.

Usage:
    merge_fronts.py <run1 calibration.nc> <run2 calibration.nc> ...
"""
import sys

import numpy as np
from netCDF4 import Dataset

SKIP = {"sample", "worker_rank", "dispatch_time", "completion_time", "objective", "target_weight", "target_boxcox_lambda",
        "trial_failed", "target_name", "birth_generation", "population", "population_rank", "population_crowding", "pareto_front"}


def main(paths):
    cost, obj, src, values, sampled, targets = [], [], [], [], [], None
    for path in paths:
        with Dataset(path) as d:
            sense = [s.strip() for s in d["objective"].sense.split(",")]
            o = np.ma.filled(d["objective"][:].astype(float), np.nan)
            params = [v for v in d.variables if v not in SKIP]
            x = np.array([np.ma.filled(d[p][:].astype(float), np.nan) for p in params]).T
            front = np.flatnonzero((np.asarray(d["pareto_front"][:]) == 1) & np.isfinite(o).all(1))
        if targets is None:
            targets = (sense, params)
        elif targets != (sense, params):
            sys.exit(f"{path} has other targets or parameters than {paths[0]}")
        sign = np.array([-1.0 if s == "maximize" else 1.0 for s in sense])
        obj.append(o[front]); cost.append(o[front] * sign); values.append(x[front]); sampled.append(x)
        src += [(path, int(i) + 1) for i in front]
    obj, cost, values, sampled = map(np.concatenate, (obj, cost, values, sampled))
    keep = [i for i in range(len(cost)) if not ((cost <= cost[i]).all(1) & (cost < cost[i]).any(1)).any()]
    for path in paths:
        n = [i for i in range(len(src)) if src[i][0] == path]
        print(f"{path}: front of {len(n)}, {len(set(n) & set(keep))} survive the merge")
    obj, cost, values, src = obj[keep], cost[keep], values[keep], [src[i] for i in keep]
    span = np.ptp(cost, axis=0)
    span[span == 0] = 1.0
    balanced = int(np.argmin(np.linalg.norm((cost - cost.min(0)) / span, axis=1)))
    picks = [(f"best objective {j + 1}", int(np.argmin(cost[:, j]))) for j in range(cost.shape[1])] + [("balanced", balanced)]
    print(f"merged front of {len(src)}")
    for name, i in picks:
        print(f"  {name:17s} sample {src[i][1]:4d} of {src[i][0]}: {np.array2string(obj[i], precision=3)}")
    full = np.nanmax(sampled, 0) - np.nanmin(sampled, 0)
    print("parameter           balanced     front range         share of sampled range")
    for k in np.argsort(-np.ptp(values, axis=0) / np.where(full > 0, full, 1)):
        print(f"  {targets[1][k]:18s} {values[balanced, k]:10.4g}   {values[:, k].min():10.4g} {values[:, k].max():10.4g}"
              f"   {np.ptp(values[:, k]) / (full[k] or 1):.2f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
