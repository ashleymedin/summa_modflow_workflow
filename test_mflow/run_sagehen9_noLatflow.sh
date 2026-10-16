#!/bin/bash
# One HRU per MODFLOW cell, 9 GRUs, without lateral flow (groundwatr = modflow); see domain_sagehen9/README.md.
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen9/summa_modflow6.config \
                      -r domain_sagehen9/simulations/spinup/heads_noLatflow.bin \
                      -s domain_sagehen9/settings/SUMMA/fileManager_noLatflow_spinup.txt \
                      ex-gwf-sagehen \
                      domain_sagehen9/settings/SUMMA/fileManager_noLatflow.txt
