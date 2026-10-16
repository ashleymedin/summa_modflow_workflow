#!/bin/bash
# run_sagehen9.sh with SFR in place of CHD (ex-gwf-sagehen-sfr); see domain_sagehen9/README.md.
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen9/summa_modflow6_sfr.config \
                      -r domain_sagehen9/simulations/spinup/heads_sfr.bin \
                      -s domain_sagehen9/settings/SUMMA/fileManager_sfr_spinup.txt \
                      ex-gwf-sagehen-sfr \
                      domain_sagehen9/settings/SUMMA/fileManager_sfr.txt
