#!/bin/bash
# Builds libmf6 and SUMMA's coupled executables on Anvil (MPI, SUNDIALS, mizuRoute, MODFLOW 6).
# Run on a login node:  ./build_coupled.anvil.bash
set -eo pipefail

source "$(cd "$(dirname "$0")" && pwd)/env_sagehen.anvil.bash"

# meson and ninja for MODFLOW 6, from a python venv; the system python is too old for meson
VENV=${VENV:-$HOME/models/venv_mf6build}
PYTHON_BASE=${PYTHON_BASE:-$(bash -lc 'module load anaconda/2024.02-py311 >/dev/null 2>&1; which python')}
[ -d "$VENV" ] || "$PYTHON_BASE" -m venv "$VENV"
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

cmake -B "$SUMMA_ROOT/build/cmake_build_anvil_mf6" -S "$SUMMA_ROOT/build" \
  -DUSE_MPI=ON \
  -DUSE_SUNDIALS=ON \
  -DUSE_MIZUROUTE=ON \
  -DUSE_MODFLOW6=ON \
  -DMODFLOW6_LIB_DIR="$MF6_ROOT/bin" \
  -DSPECIFY_LAPACK_LINKS=OFF \
  -DCMAKE_BUILD_TYPE=Release
cmake --build "$SUMMA_ROOT/build/cmake_build_anvil_mf6" --target all -j 8

ls -l "$SUMMA_ROOT"/bin/summa_modflow6_*sundials_mizuroute.exe
