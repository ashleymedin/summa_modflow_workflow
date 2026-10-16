#!/bin/bash
# run_sagehen9_sfr.sh with a GWE model (ex-gwf-sagehen-gwe) setting SUMMA's aquifer temperature;
# see domain_sagehen9/README.md. Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen9/summa_modflow6_gwe.config \
                      -r domain_sagehen9/simulations/spinup/heads_gwe.bin \
                      -s domain_sagehen9/settings/SUMMA/fileManager_gwe_spinup.txt \
                      ex-gwf-sagehen-gwe \
                      domain_sagehen9/settings/SUMMA/fileManager_gwe.txt
