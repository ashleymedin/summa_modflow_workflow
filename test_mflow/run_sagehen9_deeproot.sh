#!/bin/bash
# One HRU per MODFLOW cell, 9 GRUs, 6 m roots drawing on MODFLOW's EVT, without lateral flow;
# the distributed half of run_sagehen1_wet_deeproot.sh. See domain_sagehen9/README.md.
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen9/summa_modflow6_deeproot.config \
                      -r domain_sagehen9/simulations/spinup/heads_deeproot.bin \
                      -s domain_sagehen9/settings/SUMMA/fileManager_deeproot_spinup.txt \
                      ex-gwf-sagehen-ss \
                      domain_sagehen9/settings/SUMMA/fileManager_deeproot.txt
