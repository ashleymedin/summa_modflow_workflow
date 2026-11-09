# Sagehen calibration on Fir

Scripts for running the coupled Sagehen calibration (`test_calibration/test_calibration_sagehen.sh`) on the
Alliance Fir cluster (192 cores a node). A run directory is built on your own machine, where the Python
builders already work, and copied to Fir's scratch. Fir needs only the executables.

| File | Use |
| --- | --- |
| `env_sagehen.fir.bash` | modules, paths (`MF6_ROOT`, `SUMMA_ROOT`, `WORKFLOW`, `RUNS`) and one GRU thread per core a rank is given |
| `build_coupled.fir.bash` | build libmf6 and SUMMA's coupled executables (MPI, SUNDIALS, mizuRoute, MODFLOW 6, OpenMP) |
| `prepare_sagehen_calib.sh` | on your machine: build a lumped or grid run directory with Fir's paths in it |
| `sagehen9_timing.fir.slurm` | time the 3396-HRU grid's 72 h event at 1, 2, 4, 5 and 8 threads |
| `sagehen_calib.fir.slurm` | run one calibration from a run directory |

## Once

Clone onto Fir, matching the defaults in `env_sagehen.fir.bash` (or export your own paths):

```bash
git clone --recurse-submodules https://github.com/ashleymedin/modflow6 $HOME/models/modflow6
git -C $HOME/models/modflow6/srcextern/summa checkout multi-calib   # the submodule pins an older commit
git clone <this repository> $HOME/models/summa_modflow_workflow
```

SUNDIALS is expected at `$HOME/local/sundials`, as for SUMMA's own `utils/hpc/fir`. Then, on a login node:

```bash
$HOME/models/summa_modflow_workflow/hpc/fir/build_coupled.fir.bash
```

It needs `summa_modflow6_opt_sundials_mizuroute.exe` (the calibration) and
`summa_modflow6_serial_sundials_mizuroute.exe` (the timing) in `$SUMMA_ROOT/bin`.

## 1. Time a trial

```bash
sbatch $WORKFLOW/hpc/fir/sagehen9_timing.fir.slurm
```

Prints seconds per 72 h at each thread count. A trial scores one water year (WY2018; WY2017 is the
calibration's shared spin-up, run once), so a trial takes about 122 times the 72 h figure. The thread count
where the time stops falling is the one to give each trial; with 9 GRUs, the largest holding 721 of the
3387 cells, expect it near 5.

## 2. Prepare a run, on your machine

```bash
PYTHON=~/miniconda3/envs/modflow6/bin/python \
  ./prepare_sagehen_calib.sh grid 120 8 grid_p120 /scratch/$USER/sagehen
rsync -a runs/grid_p120/ fir:/scratch/$USER/sagehen/grid_p120/
```

The arguments are layout (`lumped` or `grid`), population, generations, a run name and Fir's `$RUNS`. The
script builds the domain with `test_mflow/tools/build_sagehen9_calibration.py`, writes the same
`config.toml` as the test (`test_calibration/sagehen_config.bash`), and puts every path the run reads
under `$RUNS/<run>`.

## 3. Run it

Ranks are population + 1: rank 0 coordinates, every other rank runs one trial at a time. NSGA-II waits
for a whole generation, so the run takes about generations x one trial when every trial of a generation
has a rank.

```bash
# the 9-HRU pre-calibration: minutes a trial
sbatch --ntasks=121 --time=1-00:00:00 $WORKFLOW/hpc/fir/sagehen_calib.fir.slurm lumped_p120

# the grid, one core a trial: fits one node at pop 120 (~generations x 5-10 h)
sbatch --ntasks=121 $WORKFLOW/hpc/fir/sagehen_calib.fir.slurm grid_p120

# the grid, 5 threads a trial, a small population on one node: checks the threading speedup in a real search
sbatch --ntasks=13 --cpus-per-task=5 $WORKFLOW/hpc/fir/sagehen_calib.fir.slurm grid_p12

# the grid, 5 threads a trial over 4 nodes (pop 120 x 5 = 600 cores)
sbatch --nodes=4 --ntasks=121 --cpus-per-task=5 $WORKFLOW/hpc/fir/sagehen_calib.fir.slurm grid_p120
```

Output, including `run.log` and the calibration file `output/sagehen_run0001_calibration.nc`, stays in
`$RUNS/<run>`.

## Seeding the grid from the pre-calibration

To start the grid search from the lumped search's final population, add to `[calibration.nsga2]` in the
grid run's `config.toml`, with `$RUNS` written out:

```toml
seed_file = "$RUNS/lumped_p120/output/sagehen_run0001_calibration.nc"
seed_fraction = 0.5
```

The best half of the grid's first generation comes from that file and the rest is random. The run stops if
a calibrated parameter is missing from the file or lies outside the grid run's bounds.

## Not in SUMMA

- **No restart.** A calibration cannot resume, so give `--time` room: Fir allows 7 days.
