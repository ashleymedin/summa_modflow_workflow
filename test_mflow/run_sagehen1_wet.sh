#!/bin/bash
# One HRU coupled to MODFLOW 6 over the domain_sagehen9 event and decisions: the lumped half of
# run_sagehen9_noLatflow.sh, compared by tools/compare_lumped_distributed.py.
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen1/summa_modflow6.config \
                      -r domain_sagehen1/simulations/spinup/heads_wet.bin \
                      -s domain_sagehen1/settings/SUMMA/fileManager_wet_spinup.txt \
                      ex-gwf-sagehen \
                      domain_sagehen1/settings/SUMMA/fileManager_wet.txt
