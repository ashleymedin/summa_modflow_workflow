#!/bin/bash
# One HRU coupled to MODFLOW 6, groundwatr = modflow, with a STEADY-STATE first stress period.
#
# Same SUMMA side as run_sagehen1.sh; the MODFLOW model (ex-gwf-sagehen-ss) adds a steady-state
# period 1 that equilibrates the water table before the coupled transient period begins.
# The coupler should report "ran 1 steady-state MODFLOW step(s)" and then couple 72 steps.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen1/summa_modflow6_steady.config \
                      ex-gwf-sagehen-ss \
                      domain_sagehen1/settings/SUMMA/fileManager.txt
