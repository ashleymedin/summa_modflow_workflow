#!/bin/bash
# Four HRUs over the same MODFLOW 6 model, WITH lateral flow (groundwatr = modLatflow).
# Guards the HRU cascade ordering in run_oneGRU; see domain_sagehen4/README.md.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen4/summa_modflow6.config \
                      ex-gwf-sagehen \
                      domain_sagehen4/settings/SUMMA/fileManager_latflow.txt
