# SUMMA / MODFLOW 6 workflow

The coupled SUMMA / MODFLOW 6 cases, their build tools and the Sagehen calibration experiments, kept out
of the SUMMA repository. SUMMA itself carries one small test of the coupling and calibration,
`utils/test/test_mflow/sagehen_calibration/`, whose domain `test_mflow/tools/build_sagehen_calibration_test.py`
builds from the files here.

- [`test_mflow/`](test_mflow/README.md) -- the coupled cases: Sagehen from 1 to 3396 HRUs (lateral flow,
  a steady-state stress period, groundwater ET, a pumping well, SFR, GWE, mizuRoute routing), Little River,
  the MODFLOW models, and the tools that build them.
- [`hpc/fir/`](hpc/fir/README.md) -- the Sagehen calibration on the Alliance Fir cluster: build, timing and
  SLURM scripts, with run directories prepared on your own machine.
- [`test_calibration/`](test_calibration/README.md) -- the Sagehen calibration (discharge, stream
  temperature, a synthetic well) and the twin experiment, with their forcing and observations, and the
  Little River well and gauge records.

The scripts find the executables in the SUMMA checkout's `bin/`, `SUMMA_ROOT`, which defaults to
`/Users/amedin/Research/MODFLOW/modflow6/srcextern/summa`:

```bash
SUMMA_ROOT=/path/to/summa ./test_mflow/run_sagehen1.sh
```

Run output (`simulations/`, `*_run/`, MODFLOW listing and binary files) is not meant to be kept; see
`.gitignore`.

Scripts that build inputs for SUMMA–MODFLOW runs, beside the `summa` checkout.

- `acquire_streamflow.py`: daily streamflow from the Water Survey of Canada or the USGS, in m3 s-1
- `acquire_stream_temperature.py`: daily water temperature at a USGS gauge, in degC
- `acquire_grace_tws.py`: the JPL GRACE/GRACE-FO mascon anomaly for a basin or a GRU, in mm (downloading needs a NASA Earthdata login)
- `acquire_groundwater_level.py`: daily or field-measured water level in a USGS well, relative to land surface in m, positive up
- `observation_netcdf.py`: the writer they share, in the form a SUMMA calibration target reads
- `add_stream_hru.py`: give a single-GRU domain a stream HRU and a one-reach mizuRoute topology
- `add_stream_network.py`: split tributary GRUs off that domain, each with its own reach, from a JSON spec

Each script's `--help` gives its arguments and an example.
