#!/bin/bash
# Little River GA, domain_littleriver_lumped, WY2011-2020, without lateral flow (groundwatr = modflow).
# Starts from its spin-up, which the first run makes.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_littleriver_lumped/summa_modflow6.config \
                      -r domain_littleriver_lumped/simulations/spinup/heads_noLatflow.bin \
                      -s domain_littleriver_lumped/settings/SUMMA/fileManager_noLatflow_spinup.txt \
                      ex-gwf-littleriver \
                      domain_littleriver_lumped/settings/SUMMA/fileManager_noLatflow.txt
