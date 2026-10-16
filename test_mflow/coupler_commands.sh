#!/bin/bash
# ---------------------------------------------------------------------------------------
# Expected in  MODFLOW_CASE  (the coupler runs here; MODFLOW reads mfsim.nam from cwd):
#   mfsim.nam                     MODFLOW 6 simulation name file (single GWF model, DIS,
#                                 RCH with READASARRAYS, TDIS TIME_UNITS SECONDS with one
#                                 time step per SUMMA data step, length unit metres)
#   (optional) the map_file named in summa_modflow6.config
#
# Expected one directory ABOVE MODFLOW_CASE:
#   summa_modflow6.config        &coupler namelist, e.g.
#                                   &coupler
#                                     mf6_model_name   = 'SAGEHEN'
#                                     rch_package_name = 'RCHA'
#                                     soil_thickness   = 2.0
#                                     map_file         = ''        ! blank => nearest-cell map
#                                     feedback         = .true.
#                                   /
#
# Expected in  SUMMA_DIR :
#   fileManager.txt              SUMMA file manager (its modelDecisions.txt must set
#                                 groundwatr = modflow  and  bcLowrSoiH = presHead)
#
# map_file format: one "iHRU  cell  weight" triple per line (whitespace separated;
# blank lines and '#' comments ignored).  An HRU may span any number of lines.
# cell is the row-major horizontal MODFLOW index (irow-1)*ncol + icol; weights are
# normalised per HRU, so put e.g. 1.0 on every line to spread an HRU over its cells.
#
# Usage:  ./coupler_commands.sh -c CONFIG [-t TOML] [-r HEADS [-s SPINUP_FILEMANAGER]] [-w HEADS]
#                               MODFLOW_CASE SUMMA_FILEMANAGER [summa_modflow6.exe]
#
#   -c, --config CONFIG   path to this case's summa_modflow6.config (required).  Each case keeps its
#                         own config beside its settings, so several cases can share one MODFLOW model
#                         directory while differing in model/package names, HRU->cell map_file and
#                         feedback.
#   -t, --toml TOML       path to a SUMMA TOML configuration file, needed only to route with mizuRoute
#                         as well (an executable built -DUSE_MIZUROUTE=ON, e.g. run_sagehen9_mizuroute.sh).
#                         Paths inside it are relative to MODFLOW_CASE, like those in CONFIG.
#   -r, --head-read HEADS start MODFLOW from the head field in HEADS; SUMMA starts from the restart
#                         its file manager names as initConditionFile.
#   -w, --head-write HEADS  write the final head field to HEADS, and SUMMA's restart on the last step.
#   -s, --spinup SPINUP_FILEMANAGER  with -r: if HEADS does not exist yet, first run this file manager
#                         as the spin-up that writes it (never routed, so -t is left out of it).
# ---------------------------------------------------------------------------------------
set -euo pipefail

CONFIG_ARG=""
TOML_ARG=""
HEAD_READ_ARG=""
HEAD_WRITE_ARG=""
SPINUP_ARG=""
while [ $# -gt 0 ]; do
  case $1 in
    -c|--config) CONFIG_ARG=${2:?"$0: -c/--config needs a path"}; shift 2 ;;
    -t|--toml)   TOML_ARG=${2:?"$0: -t/--toml needs a path"}; shift 2 ;;
    -r|--head-read)  HEAD_READ_ARG=${2:?"$0: -r/--head-read needs a path"}; shift 2 ;;
    -w|--head-write) HEAD_WRITE_ARG=${2:?"$0: -w/--head-write needs a path"}; shift 2 ;;
    -s|--spinup)     SPINUP_ARG=${2:?"$0: -s/--spinup needs a file manager"}; shift 2 ;;
    --)          shift; break ;;
    -*)          echo "$0: unknown option $1"; exit 1 ;;
    *)           break ;;
  esac
done

MODFLOW_CASE=${1:?"usage: $0 -c CONFIG MODFLOW_CASE SUMMA_FILEMANAGER [summa_modflow6_exe]"}
SUMMA_FILEMANAGER=${2:?"usage: $0 -c CONFIG MODFLOW_CASE SUMMA_FILEMANAGER [summa_modflow6_exe]"}
[ -n "$CONFIG_ARG" ] || { echo "$0: -c/--config is required (each case has its own summa_modflow6.config)"; exit 1; }
if [ $# -ge 3 ]; then
  EXE=$3
else
  # default: the coupler built into the SUMMA checkout's bin (SUMMA_ROOT)
  EXE="${SUMMA_ROOT:-/Users/amedin/Research/MODFLOW/modflow6/srcextern/summa}/bin/summa_modflow6.exe"
fi

# Resolve to absolute paths so they survive the cd into MODFLOW_CASE.
EXE=$(cd "$(dirname "$EXE")" 2>/dev/null && pwd)/$(basename "$EXE") || true
FILE_MANAGER=$(cd "$(dirname "$SUMMA_FILEMANAGER")" 2>/dev/null && pwd)/$(basename "$SUMMA_FILEMANAGER") || true
CONFIG=$(cd "$(dirname "$CONFIG_ARG")" 2>/dev/null && pwd)/$(basename "$CONFIG_ARG") || true
[ -f "$CONFIG" ] || CONFIG=$CONFIG_ARG   # unresolvable: keep what was typed, for the error below
TOML=""
if [ -n "$TOML_ARG" ]; then
  TOML=$(cd "$(dirname "$TOML_ARG")" 2>/dev/null && pwd)/$(basename "$TOML_ARG") || true
  [ -f "$TOML" ] || TOML=$TOML_ARG
fi
# restart files need only their directory to exist, since a spin-up is what writes them
HEAD_READ=""; HEAD_WRITE=""; SPINUP=""
[ -z "$HEAD_READ_ARG" ]  || HEAD_READ=$(cd "$(dirname "$HEAD_READ_ARG")" && pwd)/$(basename "$HEAD_READ_ARG")
[ -z "$HEAD_WRITE_ARG" ] || HEAD_WRITE=$(cd "$(dirname "$HEAD_WRITE_ARG")" && pwd)/$(basename "$HEAD_WRITE_ARG")
[ -z "$SPINUP_ARG" ]     || SPINUP=$(cd "$(dirname "$SPINUP_ARG")" 2>/dev/null && pwd)/$(basename "$SPINUP_ARG") || true

[ -x "$EXE" ]                    || { echo "coupler executable not found/executable: $EXE"; exit 1; }
[ -f "$MODFLOW_CASE/mfsim.nam" ] || { echo "missing $MODFLOW_CASE/mfsim.nam"; exit 1; }
[ -f "$CONFIG" ]                 || { echo "missing coupler config: $CONFIG"; exit 1; }
[ -f "$FILE_MANAGER" ]           || { echo "missing $SUMMA_FILEMANAGER"; exit 1; }
[ -z "$TOML" ] || [ -f "$TOML" ] || { echo "missing TOML configuration: $TOML"; exit 1; }
[ -z "$SPINUP" ] || [ -n "$HEAD_READ" ] || { echo "$0: -s/--spinup needs -r/--head-read"; exit 1; }
[ -z "$SPINUP" ] || [ -f "$SPINUP" ] || { echo "missing spin-up file manager: $SPINUP_ARG"; exit 1; }

# MODFLOW 6 is initialized from mfsim.nam in the working directory
cd "$MODFLOW_CASE"
if [ -n "$SPINUP" ] && [ ! -f "$HEAD_READ" ]; then
  echo "no spin-up at $HEAD_READ yet: running $SPINUP_ARG first"
  "$EXE" "$SPINUP" "$CONFIG" --head-restart-write "$HEAD_READ"
fi
"$EXE" "$FILE_MANAGER" "$CONFIG" ${TOML:+"$TOML"} \
       ${HEAD_READ:+--head-restart-read "$HEAD_READ"} ${HEAD_WRITE:+--head-restart-write "$HEAD_WRITE"}

echo "done. check (SUMMA output NetCDF vs the MODFLOW 6 listing budget):"
echo "  - scalarSoilDrainage    <-> RCH (RCHA) inflow          [coupler imposes this on MODFLOW]"
echo "  - scalarAquiferBaseflow <-> bflow package (CHD/DRN/..) discharge over the same HRUs"
echo "  - scalarAquiferStorage   = Sy * (MODFLOW water table - soil-column base); a copy for"
echo "                             inspecting the water-table position, not an independent check"
