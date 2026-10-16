#!/bin/bash
# run_sagehen9_noLatflow.sh with one pumping well (ex-gwf-sagehen-wel), from the same spin-up, making it if missing.
# tools/check_sagehen_wel.py then compares the two at the well.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen9/summa_modflow6.config \
                      -r domain_sagehen9/simulations/spinup/heads_noLatflow.bin \
                      -s domain_sagehen9/settings/SUMMA/fileManager_noLatflow_spinup.txt \
                      ex-gwf-sagehen-wel \
                      domain_sagehen9/settings/SUMMA/fileManager_wel.txt
