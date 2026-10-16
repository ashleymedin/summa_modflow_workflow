#!/bin/bash
# One HRU per MODFLOW cell, 9 GRUs, lateral flow routed by SUMMA (groundwatr = modLatflow); see domain_sagehen9/README.md.
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen9/summa_modflow6.config \
                      -r domain_sagehen9/simulations/spinup/heads_latflow.bin \
                      -s domain_sagehen9/settings/SUMMA/fileManager_latflow_spinup.txt \
                      ex-gwf-sagehen \
                      domain_sagehen9/settings/SUMMA/fileManager_latflow.txt
