#!/bin/bash
# Fir environment for the Sagehen coupled runs; sourced by the scripts beside it.
# Override any path by exporting it first.

module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module load openblas/0.3.24
module load netcdf-mpi/4.9.2
module load netcdf-fortran-mpi/4.6.1

export SUNDIALS_DIR=${SUNDIALS_DIR:-$HOME/local/sundials}
export MF6_ROOT=${MF6_ROOT:-$HOME/models/modflow6}             # ashleymedin/modflow6, libmf6 installed to its bin/
export SUMMA_ROOT=${SUMMA_ROOT:-$MF6_ROOT/srcextern/summa}     # multi-calib
export WORKFLOW=${WORKFLOW:-$HOME/models/summa_modflow_workflow}
export RUNS=${RUNS:-$SCRATCH/sagehen}                          # run directories, made by prepare_sagehen_calib.sh

# one GRU thread per core the rank was given; BLAS stays on the calling thread
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_STACKSIZE=256M
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
