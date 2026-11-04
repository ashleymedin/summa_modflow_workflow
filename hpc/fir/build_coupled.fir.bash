#!/bin/bash
# Builds libmf6 and SUMMA's coupled executables on Fir (MPI, SUNDIALS, mizuRoute, MODFLOW 6, OpenMP).
# Run on a login node:  ./build_coupled.fir.bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/env_sagehen.fir.bash"

# keep build temporaries off /tmp
mkdir -p "$SCRATCH/tmp"
export TMPDIR="$SCRATCH/tmp"

# meson and ninja for MODFLOW 6, from a python venv
module load python/3.11
VENV=${VENV:-$HOME/venvs/mf6build}
[ -d "$VENV" ] || python -m venv "$VENV"
source "$VENV/bin/activate"
pip install --quiet meson ninja

# libmf6 into $MF6_ROOT/bin
cd "$MF6_ROOT"
if [ -d builddir ]; then
  meson setup --reconfigure --prefix="$MF6_ROOT" --libdir=bin --buildtype=release builddir
else
  meson setup --prefix="$MF6_ROOT" --libdir=bin --buildtype=release builddir
fi
meson install -C builddir

# SUMMA; USE_OPENMP is ignored, with a warning, by a SUMMA checkout older than its OpenMP option
cmake -B "$SUMMA_ROOT/build/cmake_build_fir_mf6" -S "$SUMMA_ROOT/build" \
  -DUSE_MPI=ON \
  -DUSE_SUNDIALS=ON \
  -DUSE_MIZUROUTE=ON \
  -DUSE_MODFLOW6=ON \
  -DMODFLOW6_LIB_DIR="$MF6_ROOT/bin" \
  -DUSE_OPENMP=ON \
  -DSPECIFY_LAPACK_LINKS=OFF \
  -DCMAKE_BUILD_TYPE=Release
cmake --build "$SUMMA_ROOT/build/cmake_build_fir_mf6" --target all -j 8

ls -l "$SUMMA_ROOT"/bin/summa_modflow6_*sundials_mizuroute.exe
