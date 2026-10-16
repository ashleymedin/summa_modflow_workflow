#!/bin/bash
# One HRU with 6 m roots drawing on MODFLOW's EVT over the domain_sagehen9 event and decisions:
# the lumped half of run_sagehen9_deeproot.sh, compared by tools/compare_lumped_distributed.py.
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen1/summa_modflow6_deeproot.config \
                      -r domain_sagehen1/simulations/spinup/heads_wet_deeproot.bin \
                      -s domain_sagehen1/settings/SUMMA/fileManager_wet_deeproot_spinup.txt \
                      ex-gwf-sagehen-ss \
                      domain_sagehen1/settings/SUMMA/fileManager_wet_deeproot.txt
