# Sagehen calibration cases

Coupled SUMMA / MODFLOW 6 / mizuRoute calibrations on Sagehen Creek. They need the SUMMA checkout's
`bin/summa_modflow6_opt_sundials_mizuroute.exe`; set `SUMMA_ROOT` if it is not the default in the scripts.
How to configure a calibration is in SUMMA's `utils/test/test_calibration/README.md`.

## `test_calibration_sagehen.sh`

A coupled SUMMA / MODFLOW 6 / mizuRoute calibration on Sagehen Creek, scoring discharge (KGE) and
stream temperature (`T_reach` RMSE, degC) at USGS 10343500, the outlet of reach 9, and the water level
in a synthetic well (head RMSE, m, at its cell). NSGA-II scores water year 2018 after a water-year-2017
spin-up. It searches MODFLOW's K (with K33) and Sy alongside the SUMMA parameters, and
checks that each worker's model holds its last trial's K. The model is `ex-gwf-sagehen-gwe`: SFR carries
the streams and GWE the aquifer's heat, and with `deepTherml = aquiferTemp` GWE's water-table temperature
is the temperature of the baseflow the stream receives, so MODFLOW, not a SUMMA parameter, sets it. The
builder writes CND's KTS as `kts1.txt`, but it is not searched: conducting through 12 m of rock takes
years, so a year's score cannot see it (0.001 degC across x0.56 to x1.91). UZF's VKS is left at the
model's value too: the twin experiment below found these targets cannot see a VKS below the true one. The builder writes `kzone1.txt` (7 zones)
and `syzone1.txt` (3) from the model's own value classes, for zone multipliers.

Sagehen has no observation well, so the test makes one up. It sits on the lower valley side two cells
from the channel above the gauge (row 41, column 75, in GRU 7), where the water table is about 16 m
down, below 12 m of UZF. The 2017 snowmelt reaches it in February 2018, 0.30 m above its May low, and it
falls 0.12 m by the end of dry water year 2018; the valley floor
is held at land surface by the drains. Its record is the head there in a
reference coupled run at the default parameters, written by
`../test_mflow/tools/make_sagehen_synthetic_well.py` as a saved USGS daily-values response
(`sagehen/observations/SYNTHETIC-SAGEHEN-1_daily_values.json`). `acquire_groundwater_level.py --features`
turns that into the observation file exactly as it would a real well's. The target scores departures
from the water-year mean, since a well's datum is not the model's. The builder's `sagehen.obs` names the
cell as OBS6 observation `well1`, beside nine more cells the twin experiment below scores, and the target (`variable = "modflow_obs"`, `obs_name = "well1"`) reads
its head from the CSV MODFLOW writes in each trial's directory.

To remake the record, run the lumped domain once for water years 2017-2018 from the spun-up start below
(`strt1.txt` from `strt_spunup_lumped.txt`, `initConditionFile` the `state_spunup_lumped.nc`), with the
serial coupler and the MODFLOW OC saving `HEAD FREQUENCY 6`, then

```bash
make_sagehen_synthetic_well.py <run>/sagehen.hds "2016-10-01 00:00" sagehen/observations/SYNTHETIC-SAGEHEN-1_daily_values.json
../acquire_groundwater_level.py SYNTHETIC-SAGEHEN-1 2016-10-01 2018-09-30 \
    sagehen/observations/SYNTHETIC-SAGEHEN-1_daily_level.nc --model-utc-offset 0 --gauge-utc-offset -8 \
    --features sagehen/observations/SYNTHETIC-SAGEHEN-1_daily_values.json
```

```bash
./test_calibration_sagehen.sh [lumped|grid] [population] [generations] [n_ranks]   # defaults: lumped, 6, 3, 7
```

The domain is built at run time by `../test_mflow/tools/build_sagehen9_calibration.py` from the
9 D8 subcatchments of `ex-gwf-sagehen-gwe`, each with a stream HRU for its reach, and the bundled basin-mean
forcing in `sagehen/forcing/`. `lumped` gives each GRU one land HRU mapped onto all its cells (18 HRUs,
minutes a water year, UZF taking a week's drainage at a time); `grid` is `domain_sagehen9`, one land HRU per cell (3396 HRUs, hours a
trial). The MODFLOW model is copied with a TDIS as long as the forcing. The lumped layout starts SUMMA from
`sagehen/state_spunup_lumped.nc` and MODFLOW from `sagehen/strt_spunup_lumped.txt`, both made by
`../test_mflow/tools/spinup_sagehen9_cycles.py`: water years 2017-2018 cycled six times from the heads of
`spinup_sagehen9_heads.py`, carrying both models' states, UZF's included, with the heads extrapolated along
their decay after the fourth. Over the last cycle soil water changed 0.4 mm and the aquifer's volume 2.2 mm.
Every land HRU is evergreen needleleaf, and the builder sets UZF's `THTI` to carry that cover's mean
drainage, 85 mm a year, since a run started from `strt1.txt` builds UZF from it. Each
calibration's one-year spin-up then runs from that state rather than from `coldState.nc`. It needs
`bin/summa_modflow6_opt_sundials_mizuroute.exe` and a python3 with netCDF4 and pyproj.

At 4 x 2 on 5 ranks every trial runs from the spun-up start. The front runs from KGE 0.42 at 2.05 degC
and 0.026 m to KGE 0.41 at 2.12 degC and 0.019 m; eight trials of a random start are a smoke test, not a
calibration.

## `test_calibration_sagehen_twin.sh`

A twin experiment: whether discharge, stream temperature and heads can recover parameters the test
knows. On the lumped Sagehen domain of `test_calibration_sagehen.sh`, one trial of the calibration
driver with every searched bound pinned to a true value away from its default writes the reference
(the spinup still runs at the defaults, as every search trial's does):
the discharge and water temperature leaving reach 9 and the head at ten OBS6 cells, one per decile of
depth to water and spread over every GRU. `make_sagehen_twin_targets.py` takes them as daily means into
one observation file. The pinned trial, scored against them, must reproduce them. Then NSGA-II searches
from scratch against three objectives: discharge KGE, temperature RMSE, and the head RMSE pooled over
the ten cells.

| Parameter | Truth | Searched |
| --- | --- | --- |
| `K_mult` (kh1, kv1) | 2.0 | 0.2-5, log |
| `Sy_mult` | 0.7 | 0.5-1.5 |
| `rhk_mult` (SFR streambed K) | 3.0 | 0.1-10, log |
| `k_soil` | 3e-6 m s-1 | 1e-7 to 1e-5, log |
| `k_macropore` | 1e-2 m s-1 | 1e-7 to 1e-1, log |

It passes when the balanced front member, whose worst objective is nearest that objective's best on the
front, is within 10% of each parameter's searched range of the truth, on its search scale. The member
nearest the ideal point, each objective scaled by its spread over the front, is the wrong pick: heads spread
from 0.04 to 0.85 m, so it takes a member that fits discharge and temperature and misses heads by 0.55 m,
with Sy and streambed K wrong.

```bash
./test_calibration_sagehen_twin.sh [population] [generations] [n_ranks]   # defaults: 20, 15, 10
```

From the cold start it used before the spun-up state (13 October 2026), with the five parameters above and
15 generations, it passes: 300 trials in about 21 hours of computing,
none failed. The balanced member (KGE 0.964, 0.23 degC, 0.039 m) has K within 4.0% of its range, Sy 9.6%,
streambed K 3.1%, `k_soil` 0.2% and `k_macropore` 0.0%. Discharge and temperature alone are fooled by Sy
and streambed K: the best-temperature and best-KGE trials have them 20-64% off, and the heads pin them.

With VKS also searched (truth 0.4, 0.1-10) and 10 generations, it took about 13 hours and did not pass.
The true parameters reproduced their targets exactly, and the front member nearest the ideal point
recovered K (6.4% of its range), `k_soil` (2.3%) and `k_macropore` (0.7%) but missed Sy (12.8%), streambed
K (18.4%) and VKS (27.5%). The trial with the best heads (0.021 m) had Sy within 1% and streambed K within
3%, so heads do recover them; the front traded them against temperature and discharge. VKS was not
recoverable below its truth: every good trial sat at 0.11, since the targets see VKS only when it is too
high.

