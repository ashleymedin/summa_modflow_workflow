#!/bin/bash
# run_sagehen1.sh on GSFLOW's two-layer aquifer (ex-gwf-sagehen-2layer, tools/build_sagehen_2layer.py).
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen1/summa_modflow6.config \
                      ex-gwf-sagehen-2layer \
                      domain_sagehen1/settings/SUMMA/fileManager.txt
