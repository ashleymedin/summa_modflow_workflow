#!/bin/bash
# One HRU coupled to MODFLOW 6 with GROUNDWATER EVAPOTRANSPIRATION (section 8.1).
#
# rootingDepth is 6 m against a 4 m soil column, so scalarAquiferRootFrac > 0 and the deep roots
# want water from below the soil column.  soilResist ramps that demand off as the MODFLOW water
# table falls past the deepest root, the demand is imposed on MODFLOW's EVT package, and what
# MODFLOW can actually supply comes back as scalarAquiferTranspire.
#
# The MODFLOW model is ex-gwf-sagehen-ss: steady-state first period, DRN at land surface, EVT.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen1/summa_modflow6_deeproot.config \
                      ex-gwf-sagehen-ss \
                      domain_sagehen1/settings/SUMMA/fileManager_deeproot.txt
