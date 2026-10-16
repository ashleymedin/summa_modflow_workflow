# SUMMA / MODFLOW 6 workflow

The coupled SUMMA / MODFLOW 6 cases, their build tools and the Sagehen calibration experiments, kept out
of the SUMMA repository. SUMMA itself carries one small test of the coupling and calibration,
`utils/test/test_mflow/sagehen_calibration/`, whose domain `test_mflow/tools/build_sagehen_calibration_test.py`
builds from the files here.

- [`test_mflow/`](test_mflow/README.md) -- the coupled cases: Sagehen from 1 to 3396 HRUs (lateral flow,
  a steady-state stress period, groundwater ET, a pumping well, SFR, GWE, mizuRoute routing), Little River,
  the MODFLOW models, and the tools that build them.
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
