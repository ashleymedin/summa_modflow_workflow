#!/bin/bash
# Little River GA, domain_littleriver, February 2013, without lateral flow (groundwatr = modflow).
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_littleriver/summa_modflow6.config \
                      -r domain_littleriver/simulations/spinup/heads_noLatflow.bin \
                      -s domain_littleriver/settings/SUMMA/fileManager_noLatflow_spinup.txt \
                      ex-gwf-littleriver \
                      domain_littleriver/settings/SUMMA/fileManager_noLatflow.txt
