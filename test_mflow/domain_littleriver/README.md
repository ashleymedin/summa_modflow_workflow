# Little River, Georgia

The Little River Experimental Watershed above USGS 02317797 (Little River near Tifton), 334 km²,
the basin Bailey et al. (2020, *Hydrology* 7:75) ran SWAT+ gwflow on. It is low relief (82–146 m),
with a shallow surficial aquifer on the Hawthorn and streams tied to it: the seepage-dominated regime
Sagehen lacks, and real USGS wells to score heads against.

## Layouts

| Directory | HRUs | Map | Period | Run |
| --- | --- | --- | --- | --- |
| `domain_littleriver_lumped` | 21, one per GRU | area intersection of the GRU polygons (`build_map_file.py`) | WY2011–2020, WY2010 spin-up | `run_littleriver_lumped.sh` |
| `domain_littleriver` | 8349, one per cell | identity | February 2013, January spin-up | `run_littleriver.sh` |

The 21 GRUs are the D8 subcatchments of a 190-cell channel network, about 400 cells each, as
Sagehen's 9 hold about 380. The cell-per-HRU layout is for structure and short runs; a decade of it
would take days. The lumped layout is for the decade and calibration. Each has a `_noLatflow`
variant (`groundwatr = modflow` in place of `modLatflow`).

## MODFLOW model, `ex-gwf-littleriver`

- 200 m cells in EPSG:32617, 168 × 127, 8349 active. `TOP` is the base of a 2 m soil column; `BOTM`
  is land surface less SoilGrids' depth to bedrock, at least 5 m below `TOP` (6–95 m, mean 21.6).
- K 5 m/d, SY 0.12 and SS 2 × 10⁻⁶ everywhere: the middle of Bailey's calibrated zones, and GSFLOW's SS.
- UZF below the column, a land-surface DRN (`surface_discharge`), SFR along the 373 channel cells
  (`baseflow`), and OBS head at the cell of each surficial-aquifer USGS well, every step, in
  `littleriver.head.csv`.
- A steady first period under 184 mm/yr of recharge, Bailey's, then hourly steps.

Uncalibrated, the steady water table is a median 4.8 m below land, 0.4–6.4 m deeper than the wells'
median readings. The calibration is what closes that.

## Data

`tools/build_littleriver.py` reads `~/Research/MODFLOW/littleriver_data`: the NLDI basin of 02317797, a
30 m 3DEP DEM, SoilGrids depth to bedrock, and USGS well metadata. `tools/make_littleriver_forcing.py`
writes the forcing from AORC v1.1 (NOAA, hourly, 1 km), clipped by `littleriver_data/fetch_aorc.py`;
each GRU takes the mean of its pixels. Forcing is not committed.

Observations for calibration, in `utils/test/test_calibration/littleriver/observations`: daily
discharge at 02317797, whose USGS record starts in 2010, and field readings at 8 surficial-aquifer
wells, 6–9 each over WY2011–2020.

Soil is ROSETTA loamy sand and vegetation IGBP cropland/natural mosaic everywhere; land cover and
soil maps are not yet used.
