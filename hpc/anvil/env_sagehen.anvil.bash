#!/bin/bash
# Anvil environment for the Sagehen coupled runs; sourced by the scripts beside it.
# Override any path by exporting it first.

module --force purge
module load gcc/14.2.0 openmpi/4.1.6
module load netcdf-fortran/4.5.3
module load openblas/0.3.27

export SUNDIALS_DIR=${SUNDIALS_DIR:-$HOME/Summa-Actors/utils/dependencies/install/sundials}
export MF6_ROOT=${MF6_ROOT:-$HOME/models/modflow6}             # ashleymedin/modflow6, libmf6 installed to its bin/
export SUMMA_ROOT=${SUMMA_ROOT:-$MF6_ROOT/srcextern/summa}     # multi-calib
export WORKFLOW=${WORKFLOW:-$HOME/models/summa_modflow_workflow}
export RUNS=${RUNS:-/anvil/scratch/$USER/sagehen}              # run directories, made by prepare_sagehen_calib.sh
export LD_LIBRARY_PATH=$SUNDIALS_DIR/lib64:${LD_LIBRARY_PATH:-}

# one core a trial; start MPI with mpirun, as this OpenMPI has no PMI for srun
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# sbatch, refused while this account has 8 jobs queued or running, the limit asked of us on Anvil
export MAX_JOBS=${MAX_JOBS:-8}
submit() {
  local n; n=$(squeue -h -u "$USER" | wc -l)
  if [ "$n" -ge "$MAX_JOBS" ]; then echo "submit: $n jobs queued or running, the limit is $MAX_JOBS" >&2; return 1; fi
  sbatch "$@"
}
