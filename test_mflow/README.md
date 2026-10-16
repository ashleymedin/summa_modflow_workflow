# SUMMA / MODFLOW 6 coupled test cases

A thin BMI coupler runs SUMMA as the land model and MODFLOW 6 as the groundwater model.
This folder holds the coupler driver script, two MODFLOW 6 models, and eleven cases.

## Cases

| script | case | HRUs | `groundwatr` | what it is for |
|---|---|---|---|---|
| `./run_sagehen1.sh` | `domain_sagehen1` | 1 | `modflow` | the basic one-way/two-way coupling |
| `./run_sagehen4.sh` | `domain_sagehen4` | 4 | `modLatflow` | lateral flow between HRUs; guards the cascade ordering in `run_oneGRU` |
| `./run_sagehen4_noLatflow.sh` | `domain_sagehen4` | 4 | `modflow` | the same four HRUs without lateral flow |
| `./run_sagehen1_steady.sh` | `domain_sagehen1` | 1 | `modflow` | steady-state first stress period, and a DRN at land surface returning groundwater discharge as surface runoff |
| `./run_sagehen1_deeproot.sh` | `domain_sagehen1` | 1 | `modflow` | groundwater evapotranspiration through an EVT package |
| `./run_sagehen9.sh` | `domain_sagehen9` | 3396 | `modLatflow` | one HRU per MODFLOW cell in 9 GRUs; HRU area equals cell area, so there is no aggregation or area error left to correct |
| `./run_sagehen9_noLatflow.sh` | `domain_sagehen9` | 3396 | `modflow` | the same 3396 HRUs without lateral flow |
| `./run_sagehen9_mizuroute.sh` | `domain_sagehen9` | 3396 | `modLatflow` | `run_sagehen9.sh` with the 9 reaches routed by mizuRoute; needs a build with both couplers |
| `./run_sagehen9_wel.sh` | `domain_sagehen9` | 3396 | `modflow` | `run_sagehen9_noLatflow.sh` with one pumping well; `tools/check_sagehen_wel.py` shows its drawdown reach `lowerBoundHead` |
| `./run_sagehen9_deeproot.sh` | `domain_sagehen9` | 3396 | `modflow` | `run_sagehen1_deeproot.sh`'s EVT mechanism over one HRU per cell |
| `./run_sagehen1_wet.sh` | `domain_sagehen1` | 1 | `modflow` | the lumped half of `run_sagehen9_noLatflow.sh`: same event, decisions and MODFLOW model |
| `./run_sagehen1_wet_deeproot.sh` | `domain_sagehen1` | 1 | `modflow` | the lumped half of `run_sagehen9_deeproot.sh` |

`domain_sagehen4` carries two file managers, `fileManager_latflow.txt` and
`fileManager_noLatflow.txt`, whose decision files differ on the `groundwatr` line alone, so
running the pair isolates what lateral flow contributes. They write `run1_latflow*` and
`run1_noLatflow*`, so neither overwrites the other.

`run_sagehen1_steady.sh` and the three deeproot cases use `ex-gwf-sagehen-ss`, which adds a
steady-state first stress period and an EVT package; the rest share the MODFLOW 6 model
in `ex-gwf-sagehen`. Both carry a DRN at land surface, returned to SUMMA as surface runoff, so
the water table cannot rise above the ground over a spin-up. All of them run the same 72 hourly steps. Each case carries its own
`summa_modflow6.config`, so they can differ in package names, roles, HRU→cell map and feedback.

Every model carries a UZF package below the soil column, one cell per active cell bar the CHD
outlet, written by `tools/build_sagehen_uzf.py` with the upstream example's parameters and 3
trailing waves. Its starting water content carries the long-term mean recharge, so the
unsaturated zone starts in step with the settled heads rather than empty.

Every model also carries the base exchange's GHB (`GHBB`) and DRN (`DRNB`), one row per UZF cell,
written by `tools/build_sagehen_base.py` at `DIS/TOP` with zero conductance; the coupler sets them each
step. Only the `domain_sagehen9` configs name them. A lumped HRU's flux answers the mean head of its
cells, so a per-cell correction is not its linearisation, and the `domain_sagehen1` and
`domain_sagehen4` configs leave the exchange explicit, as lumped layouts never bounced. Every config
holds drainage a week before UZF (`uzf_hold_days = 7`).

Every MODFLOW model starts from the same `strt1.txt`: the heads `tools/spinup_sagehen9_heads.py`
settles on the lumped calibration domain, a steady state then 20 coupled water years 2017. The
native `ex-gwf-sagehen` heads sit near land surface, and the uplands drain from them for decades.

The February 2017 cases — `run_sagehen1_wet*.sh` and every `run_sagehen9*.sh` — start the event
from a spin-up, 2017-01-15 to 2017-02-06, 552 hourly steps with the same decisions, parameters and
MODFLOW model. The first run of a case makes its spin-up, which takes about twenty minutes for
`domain_sagehen9`; later runs restart from it. The spin-up lands in the domain's `simulations/spinup/`
— the aquifer heads as `heads_<case>.bin` and SUMMA's state as `spinup_<case>_restart_2017020623_*.nc`
— and is only remade once `heads_<case>.bin` is gone, so delete that folder after changing anything
the spin-up depends on. `run_sagehen9_mizuroute.sh` restarts from `run_sagehen9.sh`'s spin-up.

`run_sagehen9_wel.sh` runs `ex-gwf-sagehen-wel`, which is `ex-gwf-sagehen` by symbolic link bar its name
file, its output control and a WEL package pumping 500 m³ a day from one valley cell. The coupler does not
know the well is there: the drawdown reaches SUMMA through the head alone. Over the event the cell's water
table falls 0.86 m below the soil-column base, `lowerBoundHead` follows it a step behind, and the HRU's
drainage turns from rise to free drainage into UZF. Demand-driven pumping through AG has no SUMMA
counterpart, since SUMMA has no irrigation demand.

`run_sagehen9_mizuroute.sh` is the only case that also needs a TOML configuration file, which
`coupler_commands.sh` takes with `-t` and passes to the coupler as a third argument. It
configures the river network; everything SUMMA needs is still in the file manager. That case
needs `bin/summa_modflow6_mizuroute.exe`, built with `-DUSE_MODFLOW6=ON -DUSE_MIZUROUTE=ON`.

`run_sagehen1_deeproot.sh` is a mechanism test rather than a realistic Sagehen setup: it sets
`rootingDepth` to 6 m against a 4 m soil column with `rootProfil = doubleExp` and both root
scale factors at their lower bound, because the bundled parameters otherwise place 99.98% of
roots inside the soil column and groundwater ET never switches on. Note that `powerLaw` clamps
the rooting depth to the soil depth, so it can never place roots below the column.

Build the coupler first, with `-DUSE_MODFLOW6=ON` (see `build/cmake/build_mflow.mac.bash`);
it lands in `srcextern/summa/bin/summa_modflow6.exe`.

## Lumped against distributed

`run_sagehen1_wet*.sh` put the one-HRU domain on `domain_sagehen9`'s event, 2017-02-07 to
2017-02-09, and its decisions (`hc_profile = exp_prof`, `infRateMax = topmodel_GA`), so the pairs
differ in structure alone. The forcing is identical, bit for bit. `tools/compare_lumped_distributed.py`
sums each flux over the run and weights it by area, in mm, because `domain_sagehen1` is 5.47%
larger than the grid.

Every case starts the event from its own spin-up, so the soil column is in equilibrium with
its aquifer and there is about 440 mm of snow on the ground.

| mm over 72 h | lumped | distributed | + latflow | lumped, deeproot | distributed, deeproot |
|---|---:|---:|---:|---:|---:|
| rain + melt | 142.62 | 141.91 | 141.91 | 142.62 | 141.90 |
| infiltration | 117.27 | 134.58 | 128.32 | 114.00 | 139.43 |
| surface runoff | 25.34 | 7.33 | 13.59 | 28.62 | 2.47 |
| soil drainage (negative: up from the aquifer) | 0.35 | −1.96 | 0.84 | 0.24 | −2.48 |
| aquifer seepage (DRN) | 1.34 | 0.01 | 0.83 | 1.71 | 0.05 |
| lateral export to the reaches | – | – | 3.12 | – | – |
| `basin__TotalRunoff` | 26.72 | 7.35 | 17.54 | 30.36 | 2.54 |
| `averageRoutedRunoff` | 21.03 | 4.32 | 13.03 | 25.12 | 1.89 |
| `scalarTranspireLimAqfr` (mean) | – | – | – | 0.541 | 0.522 |

Drainage is small, so the runoff differences are structure. The lumped HRU turns 18% of the rain
and melt into surface runoff against 5% distributed, on either MODFLOW model. The distributed seepage
is nearly gone, 0.01 mm against 0.55 before the base exchange: the hourly bounce had been pushing
water out through the drain. The infiltration
closure is already evaluated per cell, so what remains is `surfRun_SE = homegrown_SE` taking the
infiltrating area from how full the root zone is, which the lumped column holds as one state. ET is
1.87 mm in every run, and aquifer transpiration is below 10⁻⁴ mm, since February is energy-limited.

## How the coupling works

Each SUMMA data step:

1. **(feedback)** the MODFLOW 6 water-table head from the previous step is written into SUMMA
   as the prescribed-head lower boundary of the soil column
   (BMI input `soil_water_sat-zone_top__head`, parameter `lowerBoundHead`).
2. SUMMA advances one step.
3. the drainage out the base of the SUMMA soil column
   (BMI output `soil_water__drainage_volume_flux`, flux `scalarSoilDrainage`)
   is regridded onto the MODFLOW 6 grid. Where a cell's water table is below its top, the soil-column
   base, a positive flux enters the UZF package's infiltration, as GSFLOW's gravity drainage enters
   UZF1, held `uzf_hold_days` and then taken at a constant rate over the next so hourly pulses do
   not each start a UZF wave; elsewhere it enters the RCH package `RECHARGE` array. The aquifer transpiration demand goes
   into the EVT package `RATE` array.
4. with the base exchange named, the coupler sets its GHB stage to the last head where that head
   was at or above the top, its DRN at the top where it was below, and both conductances to SUMMA's
   base conductance (BMI output `soil_bottom_surface__hydraulic_conductance`), so MODFLOW solves
   its correction to SUMMA's flux against its own new head.
5. MODFLOW 6 advances one step (prepare / do / finalize_time_step). Leading steady-state
   stress periods are solved out first, using the RCH package's own recharge rather than
   SUMMA's.
6. the new head field, each named boundary package's flow, UZF's rejected infiltration and what the
   base exchange moved are read back and aggregated per HRU, ready for step 1 of the next iteration —
   so SUMMA sees MODFLOW **with a one-step lag**. Rejected infiltration and the base exchange return to the base of the soil column
   (BMI input `soil_water__rejected_recharge_volume_flux`, `scalarAquiferReject`); SUMMA reports its
   drainage net of it, and the coupler adds it back before regridding.

Below the soil-column base SUMMA's lower boundary is free drainage into UZF; at or above it, the
head-dependent `presHead` flux carries water both ways through RCH. Over the first 10 cm below the
base one gives way linearly to the other, so the flux is continuous in the water table. Without the
base exchange, or with a jump there, cell-matched HRUs whose head sits at `DIS/TOP` reverse their
drainage every step: a trace of water past the top lands in confined storage and reads as metres of
head.

The scatter happens after `prepare_time_step`, not before: `prepare_time_step` reloads a
package's `PERIOD` block at the start of each stress period, which would overwrite it.

With `feedback = .true.` the coupler also writes back, each step, so that SUMMA's water
balance and routed streamflow include the aquifer:

    scalarAquiferStorage    Sy * (MODFLOW water table - soil-column base), diagnostic only
    scalarAquiferBaseflow   role=baseflow package outflow over the HRU footprint
    mfSurfaceDischarge      role=surface_discharge outflow, added to SUMMA's surface runoff
    scalarAquiferTranspire  role=gw_et extraction, against the demand SUMMA sent
    scalarTranspireLimAqfr  aquifer transpiration limiting factor, per MODFLOW cell

`scalarAquiferRecharge` is not exchanged — SUMMA sets it to its own soil drainage.

Boundary packages are named with a **role** rather than individually, so what SUMMA does with
a returned flux is separate from which MODFLOW package supplied it:

    bnd_package_names  = 'CHD', 'DRN', 'EVTA'
    bnd_package_roles  = 'baseflow', 'surface_discharge', 'gw_et'

Packages sharing a role are summed. `bflow_package_name` still works as a one-entry
`baseflow` table. `evt_package_name` names the EVT package whose `RATE` array SUMMA's demand
is written into; its `SURFACE` should sit at the base of the soil column and its `DEPTH` be
large, because the demand has already been reduced by `scalarTranspireLimAqfr` and a second
extinction ramp in MODFLOW would apply the same limit twice.

`head_restart_read` / `head_restart_write` save and reload the MODFLOW head field, so an
aquifer spun up once can start every later run. The calibration driver wires this up
automatically, per rank. The serial coupler also takes them on its command line, and writing
heads there makes SUMMA write its restart on the last step as well, so a spin-up leaves a
matched pair. A restarted run hands SUMMA the saved water table on step 1 rather than
`lowerBoundHead`, and restart heads replace a leading steady-state solve rather than being
overwritten by it. Boundary fluxes are not saved, so step 1 of a restarted run carries no
baseflow, seepage or groundwater ET back to SUMMA.

Start-up reports each HRU's area against its effective mapped cell area, and the run ends with
a coupled water budget: what SUMMA sent, what MODFLOW's UZF and RCH received, what UZF rejected
back to the soil column, what came back by role, and the residual.

## SUMMA model decisions

| decision | `modflow` | `modLatflow` |
|---|---|---|
| `groundwatr` | `modflow` | `modLatflow` |
| `bcLowrSoiH` | `presHead` | `presHead` |
| `hc_profile` | any | **`exp_prof`** |
| `infRateMax` | any valid for `hc_profile` | `topmodel_GA` or `noInfExc` |

`modLatflow` does everything `modflow` does and additionally runs TOPMODEL-style lateral
flow through the soil column above the water table, for hillslopes where water moves
downslope through the soil as well as recharging the aquifer. It requires `exp_prof`
because the lateral transmissivity is the vertical integral of conductivity over the soil
column alone — MODFLOW carries everything below — and `exp_prof` is the profile that
integrates to a finite base rather than assuming a shallow aquifer of its own. The lateral
flow is reported as `basin__ColumnOutflow` and added to total runoff alongside the MODFLOW
baseflow.

Both options need SUMMA built with `-DUSE_MODFLOW6=ON` (which sets `MODFLOW_ACTIVE`) and
must run through `summa_modflow6.exe`; a plain `summa` run is rejected at start-up.

## Requirements on the MODFLOW 6 model

Read from `mfsim.nam` in the working directory, it must:

* use length unit metres and `TDIS TIME_UNITS SECONDS`;
* have exactly **one MODFLOW time step per SUMMA forcing data step** — MODFLOW `delt` must
  equal SUMMA's `data_step`, taken from the forcing file attribute and constant for the run;
* contain an RCH package with `READASARRAYS`;
* contain a UZF package with one landflag-1 cell per column it covers, and none of `SIMULATE_ET`,
  `LINEAR_GWET`, `SQUARE_GWET`, `SIMULATE_GWSEEP` or `MOVER` (EVT, SUMMA and the coupler supply those);
  with a GWE model, a UZE package for it;
* be a single GWF model discretised with `DIS`.

All of these are checked by the coupler and reported by name if violated — the units, grid and
RCH and UZF checks at start-up, the time-step check on the first step — so you do not have to discover
them from a run that merely looks wrong.  Units left undefined in `mfsim.nam` draw a warning
rather than a stop, since a missing declaration is not proof of the wrong unit.

## Running

    ./coupler_commands.sh -c CONFIG [-t TOML] [-r HEADS [-s SPINUP_FILEMANAGER]] [-w HEADS]
                          MODFLOW_CASE SUMMA_FILEMANAGER [summa_modflow6.exe]

`-c/--config` is **required**: each case keeps its own config beside its settings. The
executable defaults to `../bin/summa_modflow6.exe`.

`-r HEADS` starts MODFLOW from a saved head field, and SUMMA from whatever restart its file manager
names as `initConditionFile`. `-w HEADS` writes the head field at the end, and SUMMA's restart with
it. `-s SPINUP_FILEMANAGER`, given with `-r`, runs that file manager first as the spin-up writing
`HEADS` whenever `HEADS` does not exist yet; the spin-up never routes, so `-t` is left out of it.

### MPI (`summa_modflow6_mpi.exe`)

Built alongside the plain coupler when both `-DUSE_MODFLOW6=ON` and `-DUSE_MPI=ON` are set.
It splits SUMMA's GRUs across MPI ranks the same way `summa_mpi.exe` does for a plain
(non-coupled) run; MODFLOW 6 itself stays a serial singleton on rank 0, since this libmf6
build has no PETSc/MPI support and the coupled model is one shared aquifer grid under
potentially many GRUs, not one MODFLOW grid per GRU (see the header of
`summa_modflow6_mpi.f90` for the full rationale). In `build/CMakeLists.txt` it links both
halves of that: `parallel_utils` for the GRU split (as `summa_mpi.exe` does) and
`bmif`/`MF6_LIB` for driving MODFLOW 6 (as the plain `summa_modflow6.exe` does). Run it
directly (it takes the same file manager/config arguments as the plain coupler, not through
`coupler_commands.sh`):

    cd MODFLOW_CASE && mpirun -n <nranks> /path/to/summa_modflow6_mpi.exe SUMMA_FILEMANAGER CONFIG

`<nranks>` must not exceed the number of GRUs in the run domain - a rank left with zero
GRUs hits a pre-existing SUMMA-core limitation in forcing-file reading, independent of
MODFLOW coupling (the same constraint applies to the plain `summa_mpi.exe`).

### Parameter calibration (`summa_modflow6_opt.exe`)

The calibration driver evaluates each parameter sample as a full coupled SUMMA/MODFLOW 6
simulation. It is the same program as the uncoupled `summa_opt.exe` - with `-DUSE_MODFLOW6=ON`
the calibration executable is simply renamed - and it couples whenever `simulation.use_modflow`
is set in the TOML configuration, so one binary covers both. See
`utils/test/test_calibration/README.md` for how to configure and run it.

## Where the coupling lives

The MODFLOW side is one module, `build/source/driver/mf6_coupling.f90`: the libmf6 bindings,
the HRU-to-cell map, the start-up checks, and the per-step exchange. Its interface is per-HRU
arrays in and per-HRU arrays out, so it knows nothing about SUMMA data structures, and all
three drivers share it - the serial coupler, the MPI coupler, and the calibration driver.
(The two couplers previously carried the same ~600 lines twice, verbatim.)

The SUMMA side is `build/source/driver/summa_mf6_exchange.f90`: the per-HRU getters and
setters for the exchanged quantities, plus the HRU geometry the cell map is built from. It
too has one implementation with two users - `summa_bmi.f90`'s `get_value`/`set_value` for the
coupler variable names delegate to it, and the calibration driver, which runs SUMMA without
the BMI, calls it directly.

## Coupler configuration

A Fortran namelist, conventionally `summa_modflow6.config` inside the case directory:

    &coupler
      mf6_model_name     = 'MYMODEL' ! GWF model name, as in mfsim.nam (upper case)
      rch_package_name   = 'RCHA'    ! RCH package name, as in the GWF name file (upper case)
      uzf_package_name   = 'UZF'     ! UZF package below the soil column (upper case)
      uze_package_name   = ''        ! GWE UZE package, required with gwe_model_name
      uzf_hold_days      = 1.0       ! days of drainage UZF takes at a time, at a constant rate
      base_ghb_package_name = ''     ! GHB and DRN, one row per top cell, that MODFLOW solves the
      base_drn_package_name = ''     !   base exchange through; set both or neither
      bflow_package_name = 'CHD'     ! head-dependent boundary package (CHD/DRN/RIV/GHB) whose
                                     !   simulated flow feeds back per HRU as scalarAquiferBaseflow
                                     !   ('' => skip the baseflow feedback)
      map_file           = ''        ! optional HRU->cell weight file; if blank a nearest-cell
                                     !   map is built from the MODFLOW 6 DIS grid geometry
      mf6_epsg           = 0         ! EPSG code of the MODFLOW grid's projected CRS, used only to
                                     !   reproject SUMMA's lon/lat HRU centres before the built-in
                                     !   nearest-cell map (nHRU>1).  Only WGS84 UTM is supported
                                     !   (32601-32660 N, 32701-32760 S); 0 = no reprojection, which
                                     !   is fine when nHRU==1 or with an explicit map_file
      feedback           = .true.    ! .false. => one-way (SUMMA drainage -> MODFLOW only)
    /

### `map_file` format

One `iHRU  cell  weight` triple per line, whitespace separated; blank lines and `#` comments
ignored. `cell` is the row-major horizontal MODFLOW index `(irow-1)*ncol + icol`. Weights are
normalised per HRU, so `1.0` on every line spreads an HRU evenly over its cells. An HRU may
span any number of lines.

`iHRU` is SUMMA's internal HRU index, which runs GRU by GRU. With one GRU that is the row
order of `attributes.nc`; with several it is that order only if the file is itself written
GRU by GRU, which `domain_sagehen9` is.

## Mapping HRUs to cells

HRUs may each own a disjoint set of cells — down to one cell per HRU. They need not overlap
and need not cover the whole grid; `domain_sagehen4` uses four disjoint regions.

The one requirement is that an HRU's `elevation` in `attributes.nc` be the mean land surface
of the cells it maps to, because the coupler forms

    lowerBoundHead = h_mf6 - (z_surface_HRU - soil_thickness)

so a mismatch hands the soil column a water table metres above or below it. The coupler
checks this against `DIS/TOP` plus the column depth at start-up and stops with the offending HRU named, rather
than letting SUMMA fail to converge for no visible reason.

`DIS/TOP` must be the base of SUMMA's soil column, as GSFLOW puts its MODFLOW top at the
soil-zone base, so the two models share no pore space. Every bundled grid has it 4 m below
the native land surface, which is kept as `land1.txt` for the build tools.
