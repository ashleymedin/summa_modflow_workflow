#!/bin/bash
# domain_sagehen9 again, with the 9 reaches routed by mizuRoute; see domain_sagehen9/README.md.
# Needs bin/summa_modflow6_mizuroute.exe, built with -DUSE_MODFLOW6=ON -DUSE_MIZUROUTE=ON.
# Starts from run_sagehen9.sh's spin-up, making it if missing; the reaches start empty.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen9/summa_modflow6.config \
                      -t domain_sagehen9/settings/mizuroute/mizu_control_sagehen9.toml \
                      -r domain_sagehen9/simulations/spinup/heads_latflow.bin \
                      -s domain_sagehen9/settings/SUMMA/fileManager_latflow_spinup.txt \
                      ex-gwf-sagehen \
                      domain_sagehen9/settings/SUMMA/fileManager_mizuroute.txt \
                      "${SUMMA_ROOT:-/Users/amedin/Research/MODFLOW/modflow6/srcextern/summa}/bin/summa_modflow6_mizuroute.exe"
