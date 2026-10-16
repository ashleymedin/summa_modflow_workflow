#!/bin/bash
# One HRU coupled to MODFLOW 6, groundwatr = modflow.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen1/summa_modflow6.config \
                      ex-gwf-sagehen \
                      domain_sagehen1/settings/SUMMA/fileManager.txt
