#!/bin/bash
# The same four HRUs WITHOUT lateral flow (groundwatr = modflow); the decisions differ from
# run_sagehen4.sh on that one line only, so the pair isolates what lateral flow does.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen4/summa_modflow6.config \
                      ex-gwf-sagehen \
                      domain_sagehen4/settings/SUMMA/fileManager_noLatflow.txt
